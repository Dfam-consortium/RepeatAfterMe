//! The banded DP row kernel and the fit-preferred extension driver
//! (ports of `compute_nw_row` in `bnw_extend.c` and `extend_alignment` in
//! `ram_extend.c`).
//!
//! # Model
//!
//! Extension proceeds one consensus position ("row") at a time, in one
//! direction. For each row, all four candidate consensus bases are tried: for
//! every extendable core a banded Needleman-Wunsch row is computed under the
//! candidate, and the per-core best row scores — each floored at 0 and capped
//! from below at `high_score[core] + cap_penalty` — are summed. The candidate
//! with the highest sum wins the row (ties: earlier base in A,C,G,T order),
//! the row is recomputed under the winner, and per-core high scores advance.
//! The extension is finally trimmed back to the last row where the total
//! score kept up with `min_improvement` per column, and abandoned after
//! `when_to_stop` rows without such progress.
//!
//! # DP row layout
//!
//! Each core has a band of `2*bandwidth+1` cells per row; each cell holds a
//! substitution-state and a gap-state score (ins/del are merged: an insertion
//! directly following a deletion is never counted). Only the current and
//! previous rows are kept, in a flat arena indexed
//! `[row_parity][core][band_col][state]`.
//!
//! The band's matrix-edge condition supports "pre-alignment": when a cell
//! falls off the start of the loaded sequence and the row is still inside the
//! band (`offset < 0 && row_idx < bandwidth`), it is scored as
//! `gapopen + (row_idx+1)*gapextn` rather than the out-of-bounds sentinel,
//! letting alignments begin slightly before the nominal core edge.
//!
//! PORT NOTE: sequence indices use `u64` wrapping arithmetic exactly like the
//! C `uint64_t` — a below-zero index wraps huge and is caught by the
//! upper-bound check; the one case the C guards explicitly (small
//! `start_seq_pos`) is guarded identically here.

use rayon::prelude::*;

use crate::alphabet::compl;
use crate::library::{CoreAlignment, CoreBoundFlag, SequenceLibrary};
use crate::matrix::ScoringSystem;

/// Below this many cores the per-row rayon dispatch overhead outweighs the
/// parallel kernel work and the driver stays sequential.
const PAR_CORE_THRESHOLD: usize = 64;

/// Set the number of worker threads for parallel extension (call once,
/// before the first `extend_alignment`). 0 or an Err leaves rayon's default
/// (all logical CPUs). Returns false if the global pool was already built.
pub fn set_threads(n: usize) -> bool {
    if n == 0 {
        return true;
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .is_ok()
}

/// Out-of-bounds cells outside the pre-alignment edge get this score so they
/// are never chosen (C `OOBSENTINEL`).
pub const OOB_SENTINEL: i32 = -987_654_321;
const NEG_INF: i32 = -1_000_000_000;

/// Extension direction (C `direction`: 0 = left, 1 = right).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Left,
    Right,
}

impl Direction {
    fn is_right(self) -> bool {
        matches!(self, Direction::Right)
    }
}

/// Flat two-row score arena, core-major: each core's two row-parity bands
/// are contiguous (`[core][row_parity][band_col][sub|gap]`) so parallel
/// passes can split the arena into disjoint per-core chunks.
#[derive(Clone)]
pub struct ScoreArena {
    band_cols: usize,
    data: Vec<[i32; 2]>,
}

impl ScoreArena {
    pub fn new(ncores: usize, bandwidth: i32) -> ScoreArena {
        let band_cols = (2 * bandwidth + 1) as usize;
        ScoreArena {
            band_cols,
            data: vec![[0; 2]; 2 * ncores * band_cols],
        }
    }

    #[inline]
    fn idx(&self, row_ff: usize, n: usize, col: usize) -> usize {
        (n * 2 + row_ff) * self.band_cols + col
    }

    #[inline]
    pub fn get(&self, row_ff: usize, n: usize, col: usize) -> [i32; 2] {
        self.data[self.idx(row_ff, n, col)]
    }

    #[inline]
    pub fn set(&mut self, row_ff: usize, n: usize, col: usize, v: [i32; 2]) {
        let i = self.idx(row_ff, n, col);
        self.data[i] = v;
    }

    /// Full band row as a slice (for tight kernel loops).
    #[inline]
    pub fn row(&self, row_ff: usize, n: usize) -> &[[i32; 2]] {
        let start = self.idx(row_ff, n, 0);
        &self.data[start..start + self.band_cols]
    }

    /// C boundary conditions, written into the "previous" row for row 0
    /// (parity 1): centered diagonal 0, off-diagonal
    /// `gapopen + |offset|*gapextn` in both states.
    pub fn init_boundary_row(&mut self, n: usize, bandwidth: i32, scoring: &ScoringSystem) {
        for offset in -bandwidth..=bandwidth {
            let col = (offset + bandwidth) as usize;
            let v = if offset == 0 {
                [0, 0]
            } else {
                let s = offset.abs() * scoring.gapextn + scoring.gapopen;
                [s, s]
            };
            self.set(1, n, col, v);
        }
    }
}

/// Compute one banded NW row for core `n` under consensus hypothesis
/// `cons_base`. Writes the row at parity `row_idx % 2`; reads the previous
/// parity (and the current row to its left, for insertions).
///
/// Returns `(best_row_score, max_score_sequence_idx)` where the index is
/// `row_idx + offset` of the best cell — i.e. the 0-based number of sequence
/// bases consumed at the best cell, minus one.
#[allow(clippy::too_many_arguments)]
pub fn compute_nw_row(
    direction: Direction,
    row_idx: i32,
    n: usize,
    cons_base: u8,
    core: &CoreAlignment,
    arena: &mut ScoreArena,
    lower_seq_bound: u64,
    upper_seq_bound: u64,
    sequence: &[u8],
    scoring: &ScoringSystem,
    bandwidth: i32,
) -> (i32, i32) {
    let curr_ff = (row_idx % 2) as usize;
    let prev_ff = 1 - curr_ff;

    // First base outside the core, per direction and orientation.
    let start_seq_pos: u64 = if direction.is_right() {
        if core.orient {
            core.right_seq_pos.wrapping_sub(1)
        } else {
            core.right_seq_pos.wrapping_add(1)
        }
    } else if core.orient {
        core.left_seq_pos.wrapping_add(1)
    } else {
        core.left_seq_pos.wrapping_sub(1)
    };

    let mut best_row_score = NEG_INF;
    let mut max_score_sequence_idx = 0;

    for offset in -bandwidth..=bandwidth {
        let mut ins_score = NEG_INF;
        let mut del_score = NEG_INF;
        let mut sub_score = NEG_INF;

        // Sequence progression: decreasing for Right/Rev and Left/Fwd,
        // increasing for Right/Fwd and Left/Rev.
        let seq_offset: i32 = if direction.is_right() == core.orient {
            -offset - row_idx
        } else {
            offset + row_idx
        };
        let seq_idx = start_seq_pos.wrapping_add_signed(seq_offset as i64);

        let out_of_bounds = (start_seq_pos < bandwidth as u64
            && seq_offset < 0
            && seq_offset.unsigned_abs() as u64 > start_seq_pos)
            || seq_idx > upper_seq_bound
            || seq_idx < lower_seq_bound;

        if !out_of_bounds {
            // Deletion: consensus base aligns to '-'; from previous row,
            // band column offset+1.
            if offset < bandwidth {
                let old = arena.get(prev_ff, n, (offset + 1 + bandwidth) as usize);
                del_score = del_score
                    .max(old[0] + scoring.gapopen + scoring.gapextn)
                    .max(old[1] + scoring.gapextn);
            }

            // Substitution: diagonal move, same band column in previous row.
            let seq_base = sequence[seq_idx as usize];
            let subvalue = if core.orient {
                scoring.score(cons_base, compl(seq_base))
            } else {
                scoring.score(cons_base, seq_base)
            };
            let old = arena.get(prev_ff, n, (offset + bandwidth) as usize);
            sub_score = sub_score.max(old[0] + subvalue).max(old[1] + subvalue);

            // Insertion: within-row move from band column offset-1.
            if offset > -bandwidth {
                let old = arena.get(curr_ff, n, (offset - 1 + bandwidth) as usize);
                ins_score = ins_score
                    .max(old[0] + scoring.gapopen + scoring.gapextn)
                    .max(old[1] + scoring.gapextn);
            }
        } else if offset < 0 && row_idx < bandwidth {
            // Matrix edge: pre-alignment region.
            let s = scoring.gapopen + (row_idx + 1) * scoring.gapextn;
            sub_score = s;
            ins_score = s;
            del_score = s;
        } else {
            sub_score = OOB_SENTINEL;
            ins_score = OOB_SENTINEL;
            del_score = OOB_SENTINEL;
        }

        let gap_score = ins_score.max(del_score);
        arena.set(curr_ff, n, (offset + bandwidth) as usize, [sub_score, gap_score]);

        let cell_score = sub_score.max(gap_score);
        if cell_score > best_row_score {
            max_score_sequence_idx = row_idx + offset;
            best_row_score = cell_score;
        }
    }

    (best_row_score, max_score_sequence_idx)
}

/// Per-column lane cell for the fused kernel: `[sub|gap][consensus lane]`.
pub type LaneCell = [[i32; 4]; 2];

/// Compute one banded NW row for core `n` under **all four** consensus
/// hypotheses at once, writing lane cells into `scratch` (one `LaneCell` per
/// band column). Reads the previous row from the canonical arena; never
/// writes the arena — the caller adopts the winning lane afterwards.
///
/// Equivalence with four scalar `compute_nw_row` calls holds because the
/// hypotheses share every input except the substitution value: the deletion
/// score reads only the (shared) previous row, the substitution score is
/// `max(prev_sub, prev_gap) + matrix[a][base]`, and the insertion score reads
/// the same lane's own cells to the left — exactly the scalar in-row
/// dependency. `sv_table[base]` holds `matrix[a][base]` for the four lanes.
///
/// Returns per-lane `(best_row_score, max_score_sequence_idx)`.
#[allow(clippy::too_many_arguments)]
pub fn compute_nw_row_fused(
    direction: Direction,
    row_idx: i32,
    n: usize,
    core: &CoreAlignment,
    arena: &ScoreArena,
    scratch: &mut [LaneCell],
    lower_seq_bound: u64,
    upper_seq_bound: u64,
    sequence: &[u8],
    sv_table: &[[i32; 4]; 100],
    scoring: &ScoringSystem,
    bandwidth: i32,
) -> ([i32; 4], [i32; 4]) {
    let prev_ff = 1 - (row_idx % 2) as usize;
    let go_ge = scoring.gapopen + scoring.gapextn;
    let ge = scoring.gapextn;

    let start_seq_pos: u64 = if direction.is_right() {
        if core.orient {
            core.right_seq_pos.wrapping_sub(1)
        } else {
            core.right_seq_pos.wrapping_add(1)
        }
    } else if core.orient {
        core.left_seq_pos.wrapping_add(1)
    } else {
        core.left_seq_pos.wrapping_sub(1)
    };

    let mut best = [NEG_INF; 4];
    let mut pos = [0i32; 4];
    let band_cols = (2 * bandwidth + 1) as usize;
    let decreasing = direction.is_right() == core.orient;
    let prev_row = arena.row(prev_ff, n);

    // Fast path: when every band column lands inside the sequence bounds
    // (the overwhelmingly common case), run a branchless lane loop with the
    // left-neighbor cell carried in registers.
    let s = start_seq_pos as i64;
    let (lo_idx, hi_idx) = if decreasing {
        (s - row_idx as i64 - bandwidth as i64, s - row_idx as i64 + bandwidth as i64)
    } else {
        (s + row_idx as i64 - bandwidth as i64, s + row_idx as i64 + bandwidth as i64)
    };
    if lo_idx >= lower_seq_bound as i64 && hi_idx <= upper_seq_bound as i64 {
        let mut left_sub = [NEG_INF; 4];
        let mut left_gap = [NEG_INF; 4];
        for col in 0..band_cols {
            let offset = col as i32 - bandwidth;
            let seq_off: i64 = if decreasing {
                -(offset + row_idx) as i64
            } else {
                (offset + row_idx) as i64
            };
            let seq_base = sequence[(s + seq_off) as usize];
            let sv = &sv_table[if core.orient {
                compl(seq_base) as usize
            } else {
                seq_base as usize
            }];

            let del = if col + 1 < band_cols {
                let up = prev_row[col + 1];
                (up[0] + go_ge).max(up[1] + ge)
            } else {
                NEG_INF
            };
            let diag = prev_row[col];
            let prev_max = diag[0].max(diag[1]);

            let mut sub = [0i32; 4];
            let mut gap = [0i32; 4];
            for a in 0..4 {
                sub[a] = prev_max + sv[a];
                let ins = (left_sub[a] + go_ge).max(left_gap[a] + ge);
                gap[a] = ins.max(del);
            }
            scratch[col] = [sub, gap];
            let base_pos = row_idx + offset;
            for a in 0..4 {
                let c = sub[a].max(gap[a]);
                if c > best[a] {
                    best[a] = c;
                    pos[a] = base_pos;
                }
            }
            left_sub = sub;
            left_gap = gap;
        }
        return (best, pos);
    }

    // Careful path: per-column bounds checks, matrix edge case, sentinels.
    for offset in -bandwidth..=bandwidth {
        let col = (offset + bandwidth) as usize;

        let seq_offset: i32 = if decreasing {
            -offset - row_idx
        } else {
            offset + row_idx
        };
        let seq_idx = start_seq_pos.wrapping_add_signed(seq_offset as i64);

        let out_of_bounds = (start_seq_pos < bandwidth as u64
            && seq_offset < 0
            && seq_offset.unsigned_abs() as u64 > start_seq_pos)
            || seq_idx > upper_seq_bound
            || seq_idx < lower_seq_bound;

        let cell: LaneCell;
        if !out_of_bounds {
            // Deletion (shared across lanes): previous row, column offset+1.
            let del = if offset < bandwidth {
                let up = prev_row[col + 1];
                (up[0] + go_ge).max(up[1] + ge)
            } else {
                NEG_INF
            };

            // Substitution: shared prev max plus the per-lane matrix value.
            let seq_base = sequence[seq_idx as usize];
            let sv = &sv_table[if core.orient {
                compl(seq_base) as usize
            } else {
                seq_base as usize
            }];
            let diag = prev_row[col];
            let prev_max = diag[0].max(diag[1]);

            // Insertion: this lane's own cells one column to the left.
            let left = if offset > -bandwidth {
                Some(scratch[col - 1])
            } else {
                None
            };

            let mut sub = [0i32; 4];
            let mut gap = [0i32; 4];
            for a in 0..4 {
                sub[a] = prev_max + sv[a];
                let ins = match &left {
                    Some(l) => (l[0][a] + go_ge).max(l[1][a] + ge),
                    None => NEG_INF,
                };
                gap[a] = ins.max(del);
            }
            cell = [sub, gap];
        } else if offset < 0 && row_idx < bandwidth {
            let e = scoring.gapopen + (row_idx + 1) * ge;
            cell = [[e; 4], [e; 4]];
        } else {
            cell = [[OOB_SENTINEL; 4], [OOB_SENTINEL; 4]];
        }
        scratch[col] = cell;

        for a in 0..4 {
            let c = cell[0][a].max(cell[1][a]);
            if c > best[a] {
                best[a] = c;
                pos[a] = row_idx + offset;
            }
        }
    }

    (best, pos)
}

/// Tunables for one `extend_alignment` call (all mirror C globals).
pub struct ExtendParams {
    pub bandwidth: i32,
    /// Cap on the penalty a non-extending core contributes (negative).
    pub cap_penalty: i32,
    /// Required total-score improvement per consensus column.
    pub min_improvement: i32,
    /// Maximum extension length (the C `L`).
    pub l_max: i32,
    /// Give up after this many columns without a new trimmed maximum.
    pub when_to_stop: i32,
}

/// Result of one directional extension.
pub struct ExtensionResult {
    /// Consensus bases added (the C return value, `max_score_row_idx + 1`).
    pub bp: i32,
    /// True when the C tool would print its "Extended sequence ... to the
    /// limit" warning (break on the final row; see PORT NOTE in the code).
    ///
    /// **This is not the flag you want for "did the extension run out of
    /// room".** It is reproduced for C fidelity and is wrong for that purpose:
    /// the C tests `row_idx == L - 1` *after* the loop, so an extension that
    /// runs the full `l_max` rows to natural completion — the ordinary
    /// runaway — leaves it `false`. Use [`Self::capped`].
    pub hit_limit: bool,
    /// True when the loop ended because it exhausted `l_max` rows rather than
    /// because the score stopped improving.
    ///
    /// This is the honest "ran out of room" signal: the model never found an
    /// edge, it was simply not allowed to look further. Extending to the cap
    /// in *both* directions means the copies go on agreeing indefinitely,
    /// which is a satellite array or a segmental duplication rather than an
    /// element with boundaries.
    ///
    /// Note this cannot be inferred from `bp`: that reports where the best
    /// score sat (`max_score_row_idx + 1`), which may be well short of where
    /// the loop stopped.
    pub capped: bool,
}

fn extendable(direction: Direction, core: &CoreAlignment) -> bool {
    if direction.is_right() {
        core.right_extendable
    } else {
        core.left_extendable
    }
}

/// Extend the alignment in one direction, writing chosen consensus bases into
/// `master` and updating each core's extension length and score
/// (port of `extend_alignment`; `l` is fixed at 1 as in the C tool).
pub fn extend_alignment(
    direction: Direction,
    cores: &mut [CoreAlignment],
    arena: &mut ScoreArena,
    lib: &SequenceLibrary,
    master: &mut [u8],
    params: &ExtendParams,
    scoring: &ScoringSystem,
) -> ExtensionResult {
    let n_cores = cores.len();
    let l_max = params.l_max;
    let l_spacer = 1i32; // the C `l`
    let band_cols = (2 * params.bandwidth + 1) as usize;

    let mut overall_high = vec![0i32; n_cores];
    let mut overall_pos = vec![0i32; n_cores];
    let mut trimmed_high = vec![0i32; n_cores];
    let mut trimmed_pos = vec![0i32; n_cores];

    // Fused-lane workspace: one row of lane cells per core, plus the
    // per-lane best score/position the row computation returns.
    let mut scratch = vec![[[0i32; 4]; 2]; n_cores * band_cols];
    let mut best4 = vec![[0i32; 4]; n_cores];
    let mut pos4 = vec![[0i32; 4]; n_cores];
    // Transposed matrix slice: sv_table[base][a] = matrix[a][base].
    let mut sv_table = [[0i32; 4]; 100];
    for (b, entry) in sv_table.iter_mut().enumerate() {
        for (a, sv) in entry.iter_mut().enumerate() {
            *sv = scoring.matrix[a][b];
        }
    }

    for n in 0..n_cores {
        arena.init_boundary_row(n, params.bandwidth, scoring);
    }

    let mut max_score = 0i32;
    let mut max_score_row_idx = -1i32;
    let mut last_row_idx = 0i32;

    // Parallelize across cores only when there is enough work per row to
    // amortize the rayon dispatch. Results are identical either way: every
    // work item writes a disjoint scratch/output slot from shared inputs.
    let parallel = n_cores >= PAR_CORE_THRESHOLD && rayon::current_num_threads() > 1;

    let mut row_idx = 0i32;
    while row_idx < l_max {
        last_row_idx = row_idx;

        // One fused pass per core evaluates all four consensus hypotheses.
        let run_kernel = |n: usize,
                          core: &CoreAlignment,
                          scr: &mut [LaneCell],
                          best: &mut [i32; 4],
                          pos: &mut [i32; 4]| {
            if !extendable(direction, core) {
                return;
            }
            let (b, p) = compute_nw_row_fused(
                direction,
                row_idx,
                n,
                core,
                arena,
                scr,
                core.lower_seq_bound,
                core.upper_seq_bound,
                &lib.sequence,
                &sv_table,
                scoring,
                params.bandwidth,
            );
            *best = b;
            *pos = p;
        };
        if parallel {
            // with_min_len keeps work-stealing granularity coarse enough
            // that large pools don't thrash on ~µs-sized per-core items.
            cores
                .par_iter()
                .zip(scratch.par_chunks_mut(band_cols))
                .zip(best4.par_iter_mut())
                .zip(pos4.par_iter_mut())
                .enumerate()
                .with_min_len(8)
                .for_each(|(n, (((core, scr), best), pos))| {
                    run_kernel(n, core, scr, best, pos)
                });
        } else {
            for (n, core) in cores.iter().enumerate() {
                let scr = &mut scratch[n * band_cols..(n + 1) * band_cols];
                let (best, pos) = (&mut best4[n], &mut pos4[n]);
                run_kernel(n, core, scr, best, pos);
            }
        }

        // Choose the best consensus base for this row (same fold order and
        // tie-breaking as the sequential C loops).
        let mut curr_extension_score = 0i32;
        let mut besta = 0u8;
        for a in 0..4usize {
            let mut score_given_cons = 0i32;
            for (n, core) in cores.iter().enumerate() {
                if !extendable(direction, core) {
                    continue;
                }
                let best = best4[n][a].max(0);
                if best >= overall_high[n] + params.cap_penalty {
                    score_given_cons += best;
                } else {
                    score_given_cons += overall_high[n] + params.cap_penalty;
                }
            }
            if score_given_cons > curr_extension_score {
                curr_extension_score = score_given_cons;
                besta = a as u8;
            }
        }

        let midx = if direction.is_right() {
            (l_max + l_spacer + row_idx) as usize
        } else {
            (l_max - row_idx - 1) as usize
        };
        master[midx] = besta;

        // Adopt the winning lane as the canonical row (replaces the C
        // recompute pass) and advance per-core records. The core-major arena
        // layout makes each core's rows one disjoint chunk.
        let curr_ff = (row_idx % 2) as usize;
        let lane = besta as usize;
        let adopt = |core: &CoreAlignment, core_rows: &mut [[i32; 2]], scr: &[LaneCell]| {
            if !extendable(direction, core) {
                return;
            }
            let row = &mut core_rows[curr_ff * band_cols..(curr_ff + 1) * band_cols];
            for (dst, cell) in row.iter_mut().zip(scr.iter()) {
                *dst = [cell[0][lane], cell[1][lane]];
            }
        };
        if parallel {
            cores
                .par_iter()
                .zip(arena.data.par_chunks_mut(2 * band_cols))
                .zip(scratch.par_chunks(band_cols))
                .with_min_len(8)
                .for_each(|((core, core_rows), scr)| adopt(core, core_rows, scr));
        } else {
            for (n, core) in cores.iter().enumerate() {
                let core_rows = &mut arena.data[n * 2 * band_cols..(n + 1) * 2 * band_cols];
                adopt(core, core_rows, &scratch[n * band_cols..(n + 1) * band_cols]);
            }
        }
        for (n, core) in cores.iter().enumerate() {
            if !extendable(direction, core) {
                continue;
            }
            let best = best4[n][lane];
            if best > overall_high[n] {
                overall_high[n] = best;
                overall_pos[n] = pos4[n][lane];
            }
        }

        // Accept as the new trimmed maximum only with min_improvement per
        // column since the previous maximum.
        if curr_extension_score
            >= max_score + (max_score_row_idx - row_idx).abs() * params.min_improvement
        {
            max_score_row_idx = row_idx;
            max_score = curr_extension_score;
            trimmed_high.copy_from_slice(&overall_high);
            trimmed_pos.copy_from_slice(&overall_pos);
        }

        if (row_idx - max_score_row_idx).abs() >= params.when_to_stop {
            break;
        }
        row_idx += 1;
    }

    // PORT NOTE: C tests `row_idx == L - 1` after the loop, so a loop that
    // runs to natural completion (row_idx == L) does NOT warn; only a break
    // on exactly the final row does. Replicated via last_row_idx.
    let hit_limit = last_row_idx == l_max - 1 && row_idx != l_max;
    // The loop is `while row_idx < l_max` with an early `break` on
    // `when_to_stop`; reaching `l_max` therefore means it was cut off.
    let capped = row_idx == l_max;

    for (n, core) in cores.iter_mut().enumerate() {
        if trimmed_high[n] > 0 && trimmed_pos[n] >= 0 {
            if direction.is_right() {
                core.right_extension_len = trimmed_pos[n] + 1;
            } else {
                core.left_extension_len = trimmed_pos[n] + 1;
            }
            core.score += trimmed_high[n];
        }
    }

    ExtensionResult {
        bp: max_score_row_idx + 1,
        hit_limit,
        capped,
    }
}

/// One overlap-avoidance adjustment (for caller-side reporting).
pub struct OverlapEvent {
    pub s_seq_idx: usize,
    pub extended_to: u64,
    pub r_seq_idx: usize,
    pub old_bound: u64,
    pub pos_in_r: u64,
    /// True when the tightened bound was the upper bound (reverse-strand r).
    pub upper: bool,
}

/// After the right extension, shrink the extension bounds of any core whose
/// left flank would reuse sequence claimed by another core's right extension
/// (port of the block between the two `extend_alignment` calls in `main`).
///
/// Bounds are tightened in place, sequentially, exactly like the C loops —
/// a bound tightened by an earlier event is the bound a later containment
/// test sees.
pub fn apply_overlap_avoidance(
    cores: &mut [CoreAlignment],
    lib: &SequenceLibrary,
) -> Vec<OverlapEvent> {
    let mut events = Vec::new();
    for s in 0..cores.len() {
        let s_core = &cores[s];
        let s_seq_idx = s_core.seq_idx;
        let s_lower = lib.lower_bound(s_seq_idx);
        let extended_pos = if s_core.orient {
            s_core
                .right_seq_pos
                .wrapping_sub(s_core.right_extension_len as u64)
        } else {
            s_core
                .right_seq_pos
                .wrapping_add(s_core.right_extension_len as u64)
        };
        let seqid_extended_pos =
            lib.offsets[s_seq_idx].wrapping_add(extended_pos.wrapping_sub(s_lower).wrapping_add(1));

        for r in 0..cores.len() {
            let r_core = &cores[r];
            let r_seq_idx = r_core.seq_idx;
            if lib.identifiers[s_seq_idx] != lib.identifiers[r_seq_idx]
                || seqid_extended_pos <= lib.offsets[r_seq_idx]
            {
                continue;
            }
            let r_lower = lib.lower_bound(r_seq_idx);
            let pos_in_r = r_lower + (seqid_extended_pos - lib.offsets[r_seq_idx]);

            let r_core = &mut cores[r];
            if r_core.orient {
                if pos_in_r >= r_core.left_seq_pos && pos_in_r <= r_core.upper_seq_bound {
                    events.push(OverlapEvent {
                        s_seq_idx,
                        extended_to: seqid_extended_pos,
                        r_seq_idx,
                        old_bound: r_core.upper_seq_bound,
                        pos_in_r,
                        upper: true,
                    });
                    r_core.upper_seq_bound = pos_in_r;
                    r_core.upper_seq_bound_flag = CoreBoundFlag::ExtBoundary;
                }
            } else if pos_in_r >= r_core.lower_seq_bound && pos_in_r <= r_core.left_seq_pos {
                events.push(OverlapEvent {
                    s_seq_idx,
                    extended_to: seqid_extended_pos,
                    r_seq_idx,
                    old_bound: r_core.lower_seq_bound,
                    pos_in_r,
                    upper: false,
                });
                r_core.lower_seq_bound = pos_in_r;
                r_core.lower_seq_bound_flag = CoreBoundFlag::ExtBoundary;
            }
        }
    }
    events
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{char_to_num, num_to_char};
    use crate::library::CoreBoundFlag;

    fn encode(s: &str) -> Vec<u8> {
        s.bytes().map(|b| char_to_num(b).unwrap()).collect()
    }

    fn simple_lib(seq: &str) -> SequenceLibrary {
        SequenceLibrary {
            sequence: encode(seq),
            identifiers: vec!["seq1".into()],
            boundaries: vec![seq.len() as u64],
            offsets: vec![0],
        }
    }

    fn core_at(left: u64, right: u64, upper: u64) -> CoreAlignment {
        CoreAlignment {
            seq_idx: 0,
            left_seq_pos: left,
            right_seq_pos: right,
            left_extendable: true,
            right_extendable: true,
            lower_seq_bound: 0,
            upper_seq_bound: upper,
            lower_seq_bound_flag: CoreBoundFlag::SeqBoundary,
            upper_seq_bound_flag: CoreBoundFlag::SeqBoundary,
            left_extension_len: 0,
            right_extension_len: 0,
            score: 0,
            orient: false,
        }
    }

    /// Two identical copies right of their cores must extend along the exact
    /// shared sequence, reproducing it in the consensus.
    #[test]
    fn identical_copies_extend_to_the_shared_sequence() {
        //            core      ext
        let seq = "ACGTACGTGATTACAGATTACAXACGTACGTGATTACAGATTACA"
            .replace('X', "T");
        let lib = simple_lib(&seq);
        let n = lib.sequence.len() as u64;
        let mut cores = vec![core_at(0, 7, n - 1), core_at(23, 30, n - 1)];
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let params = ExtendParams {
            bandwidth: 5,
            cap_penalty: -90,
            min_improvement: 9, // one positively-scoring sequence suffices here
            l_max: 50,
            when_to_stop: 100,
        };
        let mut arena = ScoreArena::new(cores.len(), params.bandwidth);
        let mut master = vec![99u8; (2 * params.l_max + 1 + 1) as usize];
        let res = extend_alignment(
            Direction::Right,
            &mut cores,
            &mut arena,
            &lib,
            &mut master,
            &params,
            &scoring,
        );
        assert!(res.bp > 0, "no extension happened");
        let got: String = master[(params.l_max + 1) as usize..]
            [..res.bp as usize]
            .iter()
            .map(|&b| num_to_char(b) as char)
            .collect();
        // Both flanks read GATTACAGATTACAT... after the core.
        assert!(
            "GATTACAGATTACAT".starts_with(&got[..got.len().min(15)]),
            "consensus was {got}"
        );
        assert_eq!(cores[0].right_extension_len, res.bp);
        assert!(cores[0].score > 0);
    }

    /// A reverse-orient core must extend using complemented bases.
    #[test]
    fn reverse_orient_core_reads_complemented_flank() {
        // Forward copy: core ACGTACG then flank GATTACA.
        // Reverse copy embeds revcomp(ACGTACG GATTACA) = TGTAATC CGTACGT.
        let seq = "ACGTACGGATTACANNNNTGTAATCCGTACGT";
        let lib = simple_lib(seq);
        let n = lib.sequence.len() as u64;
        let mut fwd = core_at(0, 6, n - 1);
        fwd.left_extendable = false;
        // Reverse core occupies positions 25..=31 (CGTACGT = revcomp of core);
        // logical left = 31, right = 25; its right extension reads 24,23,...
        let mut rev = core_at(31, 25, n - 1);
        rev.orient = true;
        rev.left_extendable = false;
        let mut cores = vec![fwd, rev];
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let params = ExtendParams {
            bandwidth: 3,
            cap_penalty: -90,
            min_improvement: 9,
            l_max: 20,
            when_to_stop: 100,
        };
        let mut arena = ScoreArena::new(cores.len(), params.bandwidth);
        let mut master = vec![99u8; (2 * params.l_max + 1 + 1) as usize];
        let res = extend_alignment(
            Direction::Right,
            &mut cores,
            &mut arena,
            &lib,
            &mut master,
            &params,
            &scoring,
        );
        assert!(res.bp >= 7, "extended only {} bp", res.bp);
        let got: String = master[(params.l_max + 1) as usize..][..7]
            .iter()
            .map(|&b| num_to_char(b) as char)
            .collect();
        assert_eq!(got, "GATTACA");
        assert_eq!(cores[1].right_extension_len, cores[0].right_extension_len);
    }

    /// The fused 4-lane kernel must agree cell-for-cell with four scalar
    /// `compute_nw_row` calls, across rows, directions, and orientations.
    #[test]
    fn fused_kernel_matches_scalar_kernel_exactly() {
        // Deterministic pseudo-random sequence (LCG), including some Ns.
        let mut state = 0x2545F491u64;
        let seq: Vec<u8> = (0..120)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                match (state >> 33) % 9 {
                    0..=1 => 0u8,
                    2..=3 => 1,
                    4..=5 => 2,
                    6..=7 => 3,
                    _ => crate::alphabet::SYM_N,
                }
            })
            .collect();
        let lib = SequenceLibrary {
            sequence: seq,
            identifiers: vec!["s".into()],
            boundaries: vec![120],
            offsets: vec![0],
        };
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let bandwidth = 6;
        let band_cols = (2 * bandwidth + 1) as usize;
        let mut sv_table = [[0i32; 4]; 100];
        for (b, entry) in sv_table.iter_mut().enumerate() {
            for (a, sv) in entry.iter_mut().enumerate() {
                *sv = scoring.matrix[a][b];
            }
        }

        // Wide bounds exercise the branchless fast path; tight bounds force
        // out-of-bounds columns through the careful path (including the
        // matrix-edge case in early rows).
        for (lo, hi) in [(0u64, 119u64), (46, 63)] {
        for direction in [Direction::Left, Direction::Right] {
            for orient in [false, true] {
                let mut core = core_at(50, 60, hi);
                core.lower_seq_bound = lo;
                core.orient = orient;
                if orient {
                    // Logical left/right swap for a reverse-strand core.
                    core.left_seq_pos = 60;
                    core.right_seq_pos = 50;
                }

                // Canonical arena evolves by adopting lane (row % 4) each
                // row, exercising all lanes as the in-row dependency source.
                let mut canonical = ScoreArena::new(1, bandwidth);
                canonical.init_boundary_row(0, bandwidth, &scoring);
                let mut scratch = vec![[[0i32; 4]; 2]; band_cols];

                for row_idx in 0..12i32 {
                    let (best4, pos4) = compute_nw_row_fused(
                        direction, row_idx, 0, &core, &canonical, &mut scratch,
                        lo, hi, &lib.sequence, &sv_table, &scoring, bandwidth,
                    );
                    for a in 0u8..4 {
                        let mut trial = canonical.clone();
                        let (best, pos) = compute_nw_row(
                            direction, row_idx, 0, a, &core, &mut trial,
                            lo, hi, &lib.sequence, &scoring, bandwidth,
                        );
                        assert_eq!(best, best4[a as usize],
                            "best mismatch: dir {direction:?} orient {orient} row {row_idx} lane {a}");
                        assert_eq!(pos, pos4[a as usize],
                            "pos mismatch: dir {direction:?} orient {orient} row {row_idx} lane {a}");
                        let curr_ff = (row_idx % 2) as usize;
                        for col in 0..band_cols {
                            let s = trial.get(curr_ff, 0, col);
                            let f = scratch[col];
                            assert_eq!(
                                s, [f[0][a as usize], f[1][a as usize]],
                                "cell mismatch at col {col}: dir {direction:?} orient {orient} row {row_idx} lane {a}"
                            );
                        }
                    }
                    // Adopt one lane into the canonical arena for the next row.
                    let lane = (row_idx % 4) as usize;
                    let curr_ff = (row_idx % 2) as usize;
                    for col in 0..band_cols {
                        let f = scratch[col];
                        canonical.set(curr_ff, 0, col, [f[0][lane], f[1][lane]]);
                    }
                }
            }
        }
        }
    }

    /// The when_to_stop guard must abandon extension into unrelated sequence.
    #[test]
    fn extension_stops_in_unrelated_sequence() {
        let seq = "ACGTACGTAAAACCCCGGGGTTTTACGTACGTCTCTCTAGAGAGGCGCGC";
        let lib = simple_lib(seq);
        let n = lib.sequence.len() as u64;
        let mut cores = vec![core_at(0, 7, n - 1), core_at(24, 31, n - 1)];
        cores[0].left_extendable = false;
        cores[1].left_extendable = false;
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let params = ExtendParams {
            bandwidth: 3,
            cap_penalty: -90,
            min_improvement: 18, // require both copies to keep scoring
            l_max: 40,
            when_to_stop: 5,
        };
        let mut arena = ScoreArena::new(cores.len(), params.bandwidth);
        let mut master = vec![99u8; (2 * params.l_max + 1 + 1) as usize];
        let res = extend_alignment(
            Direction::Right,
            &mut cores,
            &mut arena,
            &lib,
            &mut master,
            &params,
            &scoring,
        );
        // The two flanks diverge immediately; nothing should be kept.
        assert_eq!(res.bp, 0, "kept {} bp of unrelated extension", res.bp);
    }
}
