//! Profile-scored consensus extension: an alternate to `engine::extend_alignment`.
//!
//! Imports go through `super` so the file also compiles as a `#[path]`
//! submodule of a shim that re-exports `alphabet` and `library` from ram-core
//! (how ram-denovo uses it ahead of a ram-core release).
//!
//! Same outer loop as the fit-preferred extension (grow the consensus one
//! column at a time away from the seed word, every core held in a banded
//! alignment to the growing consensus, stop after `when_to_stop` columns
//! without a new maximum, trim back to the maximum), with a different notion
//! of what a column is worth:
//!
//! * A core's base `x` at a column scores `log2(p(x) / bg(x))` bits, where `p`
//!   is the column's base distribution over the cores still active, computed
//!   with the core itself left out (jackknife, so a random column has
//!   expectation <= 0) and smoothed by a Dirichlet prior of weight `alpha`
//!   centred on the plurality base at `prior_div` divergence (transition
//!   aware). Small families therefore score like a fixed-divergence matrix,
//!   large ones are data-driven.
//! * `bg` is the genome composition, except inside a low-complexity run:
//!   when the mean smoothed column distribution of the last `bg_window`
//!   columns has entropy below `BG_ENTROPY_FULL` bits it replaces the genome
//!   composition (fading in up to `BG_ENTROPY_NONE`). A homopolymer therefore
//!   carries no information: its columns agree, but no better than the run
//!   predicts (without this, a poly-A tail shared by a fifth of a large Alu
//!   family extends the consensus a hundred bases into A-rich flank), while
//!   an element's own skewed composition is left alone (raising the
//!   background to it halves the information of every column of an old
//!   family and trims it to its most conserved stretch). Inside such a run
//!   the jackknife estimate sits strictly below the smoothed background, so
//!   the run counts against the extension without any per-core floor (a
//!   floor scales with the family size and ends a weakly informative
//!   stretch of a large old family early).
//! * Gaps are affine in bits. Masked and N bases neither vote nor match.
//! * The stop statistic is the cumulative column information (sum over
//!   active cores of their own jackknife score), not the alignment total.
//! * A core is frozen when its running score falls `cap_bits` below its peak,
//!   or when it has gained nothing over the last `window` columns (it no
//!   longer shares sequence here). Frozen cores stop voting and stop counting
//!   as support; the extension stops when fewer than
//!   `max(min_active, min_active_frac * n)` cores remain active.
//! * A synchronized freeze is an element boundary, once it has been probed.
//!   Among the cores that had gained at least `burst_min_bits` (so a seed
//!   family's hangers-on, which never share the extension, do not count),
//!   when at least `burst_frac` of those active `burst_window` columns ago
//!   froze within those columns, a tentative boundary is set where that
//!   cohort's scores peaked (median over the cohort, never past the
//!   cumulative information peak). The extension goes on with the survivors
//!   for `burst_probe` more columns while the cohort, still frozen, keeps
//!   aligning to the consensus the survivors build. If half the cohort
//!   regains `REVIVE_BITS` over its freeze score, the element had only
//!   dropped in conservation (a fade): the cohort rejoins and the boundary
//!   is dropped. Otherwise the cohort really ended there, and the extension
//!   is trimmed back to the boundary whatever the survivors went on to say.
//!   This is what keeps a minority with a shared neighbour (the copies of a
//!   segmental duplication, an Alu in another Alu's tail, the internal
//!   region behind a solo-LTR majority) from carrying the consensus past
//!   the element on their own say-so, without cutting an old family at its
//!   least conserved stretch. Fragments of an old family end at scattered
//!   positions and never trip it.
//!
//! Nothing in this module is used by the fit-preferred path; callers opt in
//! by calling [`profile_extend`] instead of `extend_alignment`.

use super::alphabet::SYM_N;
use super::library::{CoreAlignment, SequenceLibrary};
use rayon::prelude::*;
use std::collections::VecDeque;

/// Tunables for one [`profile_extend`] call.
#[derive(Debug, Clone)]
pub struct ProfileParams {
    /// Half-width of the alignment band, in bp.
    pub bandwidth: i32,
    /// Maximum extension length per direction.
    pub l_max: i32,
    /// Give up after this many columns without a new information maximum.
    pub when_to_stop: i32,
    /// Freeze a core whose running score is this many bits below its peak.
    pub cap_bits: f32,
    /// Freeze a core that gained nothing over this many columns (0 = off).
    pub window: usize,
    /// Stop when fewer than this fraction of the cores are still active.
    pub min_active_frac: f32,
    /// ... and never continue with fewer than this many.
    pub min_active: usize,
    /// Dirichlet prior weight on the column distribution.
    pub alpha: f32,
    /// Divergence the prior is centred at (fraction of non-plurality bases).
    pub prior_div: f32,
    /// Affine gap costs, bits.
    pub gap_open: f32,
    pub gap_ext: f32,
    /// Bits are multiplied by this and rounded into `CoreAlignment::score`.
    pub score_scale: f32,
    /// Raise the background to the composition of this many recent columns
    /// when that is higher than the genome's (0 = genome composition only).
    pub bg_window: usize,
    /// Synchronized-freeze stop: look back this many columns (0 = off) ...
    pub burst_window: usize,
    /// ... and stop when this fraction of the cores active at the start of
    /// the look-back froze within it ...
    pub burst_frac: f32,
    /// ... counting only cores that had gained at least this many bits ...
    pub burst_min_bits: f32,
    /// ... and confirmed unless the cohort revives within this many columns.
    pub burst_probe: usize,
}

impl Default for ProfileParams {
    fn default() -> Self {
        ProfileParams {
            bandwidth: 14,
            l_max: 7000,
            when_to_stop: 100,
            cap_bits: 15.0,
            window: 60,
            min_active_frac: 0.2,
            min_active: 3,
            alpha: 10.0,
            prior_div: 0.25,
            gap_open: 4.0,
            gap_ext: 1.5,
            score_scale: 10.0,
            bg_window: 20,
            burst_window: 30,
            burst_frac: 0.5,
            burst_min_bits: 8.0,
            burst_probe: 120,
        }
    }
}

/// Result of a two-sided [`profile_extend`].
pub struct ProfileExtension {
    /// Consensus bases added on the left, in consensus order (leftmost first);
    /// numeric alphabet, `SYM_N` where no core voted.
    pub left: Vec<u8>,
    /// Consensus bases added on the right, in consensus order.
    pub right: Vec<u8>,
    /// Cores still active when each direction stopped.
    pub left_active: usize,
    pub right_active: usize,
    /// Whether each direction was ended by the synchronized-freeze rule.
    pub left_burst: bool,
    pub right_burst: bool,
    /// Peak information per core (bits), left + right.
    pub core_bits: Vec<f32>,
}

const BG: [f32; 4] = [0.285, 0.215, 0.215, 0.285];
/// Local-composition background applies fully below this entropy (bits of
/// the window-mean composition) and not at all above the second.
const BG_ENTROPY_FULL: f32 = 1.0;
const BG_ENTROPY_NONE: f32 = 1.5;
/// A burst needs at least this many frozen cores to count.
const BURST_MIN_CORES: usize = 5;
/// A probed cohort core has revived when it is this many bits above the
/// score it froze at.
const REVIVE_BITS: f32 = 8.0;
const NEG: f32 = -1.0e9;
const PAR_THRESHOLD: usize = 256;

/// No base: masked or N (scores like a mismatch, the column is still
/// consumed).
const BASE_NONE: u8 = 4;
/// No position: outside the core's bounds (a sequence end or a neighbour
/// bound). Not alignable: the band may not consume it. Scoring it like a
/// masked base (-2) was cheaper than a gap, so the band walked past the
/// sequence end and the reported extension followed (hs1 chr1 telomere).
const BASE_OUT: u8 = 5;

/// Base of `core` at consensus offset `t` from its logical left core base
/// (negative = left of the core), in consensus orientation. `BASE_NONE` for
/// a masked base or N, `BASE_OUT` outside the core's bounds.
#[inline]
fn base_at(lib: &SequenceLibrary, core: &CoreAlignment, t: i64) -> u8 {
    let pos = if core.orient {
        core.left_seq_pos as i64 - t
    } else {
        core.left_seq_pos as i64 + t
    };
    if pos < core.lower_seq_bound as i64 || pos > core.upper_seq_bound as i64 {
        return BASE_OUT;
    }
    let c = lib.sequence[pos as usize];
    if c > 3 {
        return 4;
    }
    if core.orient {
        3 - c
    } else {
        c
    }
}

/// Jackknife score table for one column: `table[own][x]` is the score of base
/// `x` for a core whose own frontier base is `own` (4 = none, no leave-out).
fn column_prior(counts: &[f32; 4], p: &ProfileParams) -> [f32; 4] {
    let plur = (0..4).max_by(|&a, &b| counts[a].partial_cmp(&counts[b]).unwrap()).unwrap();
    let tr = [2usize, 3, 0, 1][plur];
    let mut prior = [p.prior_div / 4.0; 4];
    prior[plur] = 1.0 - p.prior_div;
    prior[tr] = p.prior_div / 2.0;
    prior
}

/// The column's smoothed distribution (no leave-out): what it contributes
/// to the local background of the columns after it.
fn column_smoothed(counts: &[f32; 4], n_act: f32, p: &ProfileParams) -> [f32; 4] {
    let prior = column_prior(counts, p);
    let denom = (n_act + p.alpha).max(1.0);
    let mut q = [0f32; 4];
    for x in 0..4 {
        q[x] = (counts[x] + p.alpha * prior[x]) / denom;
    }
    q
}

fn column_table(counts: &[f32; 4], n_act: f32, bg: &[f32; 4], p: &ProfileParams) -> [[f32; 4]; 5] {
    let denom = (n_act - 1.0 + p.alpha).max(1.0);
    let prior = column_prior(counts, p);
    let mut t = [[0f32; 4]; 5];
    for own in 0..5 {
        for x in 0..4 {
            let cnt = counts[x] - if own == x { 1.0 } else { 0.0 };
            t[own][x] = ((cnt.max(0.0) + p.alpha * prior[x]) / denom / bg[x]).log2();
        }
    }
    t
}

/// Outcome of one direction of [`extend_dir`].
struct DirResult {
    /// Consensus bases in growth order.
    cons: Vec<u8>,
    /// Bases consumed per core at the trimmed maximum.
    consumed: Vec<i32>,
    /// Cores active at the end.
    active: usize,
    /// Peak bits per core.
    peak: Vec<f32>,
    /// Ended by the synchronized-freeze rule.
    burst: bool,
}

/// One direction. `start_off` is the consensus offset of the first new
/// column (`k` for right, `-1` for left) and `dir` the offset step.
fn extend_dir(
    cores: &[CoreAlignment],
    lib: &SequenceLibrary,
    start_off: i64,
    dir: i64,
    p: &ProfileParams,
) -> DirResult {
    let n = cores.len();
    let w = p.bandwidth.max(1) as usize;
    let b = 2 * w + 1;
    let mut m = vec![NEG; n * b];
    let mut d = vec![NEG; n * b];
    let mut ins = vec![NEG; n * b];
    for i in 0..n {
        m[i * b + w] = 0.0;
    }
    let mut mn = vec![NEG; n * b];
    let mut dn = vec![NEG; n * b];
    let mut inn = vec![NEG; n * b];
    let mut peak = vec![0f32; n];
    let mut peak_c = vec![0usize; n];
    let mut cur = vec![0f32; n];
    let mut cur_d = vec![0i32; n];
    // consumed bases a core reports: live while active, held at the freeze
    // point while frozen (its alignment keeps running for revival checks only)
    let mut rep_d = vec![0i32; n];
    let mut nxt = vec![4u8; n];
    let mut frozen = vec![false; n];
    let mut frozen_at = vec![usize::MAX; n];
    let mut freeze_score = vec![0f32; n];
    let mut hist: VecDeque<Vec<f32>> = VecDeque::new();
    let mut cons: Vec<u8> = Vec::new();
    let mut total = 0f32;
    let mut best_total = 0f32;
    let mut best_c = 0usize;
    let mut best_d = vec![0i32; n];
    let min_active = p.min_active.max((p.min_active_frac * n as f32).ceil() as usize);
    let parallel = n >= PAR_THRESHOLD && rayon::current_num_threads() > 1;
    let mut active_now = n;
    // local composition: smoothed column distributions over the last
    // bg_window columns; columns not yet seen count as genome composition
    let mut bg_hist: VecDeque<[f32; 4]> = VecDeque::new();
    let mut bg_sum = [0f32; 4];
    // synchronized-freeze bookkeeping: supporting cores (peak >= burst_min_bits)
    // frozen per column, supporters active at each column's start, and
    // reported-consumed snapshots deep enough to trim back to where a burst
    // cohort peaked
    let mut newly_frozen: VecDeque<usize> = VecDeque::new();
    let mut supporters_at_start: VecDeque<usize> = VecDeque::new();
    let mut d_hist: VecDeque<Vec<i32>> = VecDeque::new();
    let d_depth = p.burst_window + p.window + 2;
    // a tentative boundary under probe: (column, consumed snapshot, probe end)
    let mut tentative: Option<(usize, Vec<i32>, usize)> = None;
    let mut probing = vec![false; n];
    let mut burst = false;

    let mut c = 0usize;
    while (c as i32) < p.l_max {
        active_now = frozen.iter().filter(|f| !**f).count();
        if active_now < min_active {
            break;
        }
        let supporters_now = (0..n).filter(|&i| !frozen[i] && peak[i] >= p.burst_min_bits).count();
        // frontier base per core: continuation of its best cell
        for i in 0..n {
            nxt[i] = if frozen[i] { 4 } else { base_at(lib, &cores[i], start_off + dir * cur_d[i] as i64).min(BASE_NONE) };
        }
        let mut counts = [0f32; 4];
        for i in 0..n {
            if !frozen[i] && nxt[i] < 4 {
                counts[nxt[i] as usize] += 1.0;
            }
        }
        let n_act: f32 = counts.iter().sum();
        let mut bg = BG;
        if p.bg_window > 0 && !bg_hist.is_empty() {
            let win = p.bg_window as f32;
            let missing = win - bg_hist.len() as f32;
            let mut local = [0f32; 4];
            let mut h = 0f32;
            for x in 0..4 {
                local[x] = (bg_sum[x] + missing * BG[x]) / win;
                if local[x] > 0.0 {
                    h -= local[x] * local[x].log2();
                }
            }
            let wgt = ((BG_ENTROPY_NONE - h) / (BG_ENTROPY_NONE - BG_ENTROPY_FULL)).clamp(0.0, 1.0);
            if wgt > 0.0 {
                for x in 0..4 {
                    let raised = BG[x] + wgt * (local[x] - BG[x]);
                    if raised > bg[x] {
                        bg[x] = raised;
                    }
                }
            }
        }
        let table = column_table(&counts, n_act, &bg, p);
        let colinfo: f32 = (0..n)
            .filter(|&i| !frozen[i] && nxt[i] < 4)
            .map(|i| table[nxt[i] as usize][nxt[i] as usize])
            .sum();

        // DP step for every core; frozen cores no longer vote (own = 4, no
        // leave-out) but keep aligning so a burst cohort can show revival
        let cc = c as i64;
        let step = |i: usize,
                    core: &CoreAlignment,
                    mo: &[f32],
                    dd: &[f32],
                    io: &[f32],
                    mn: &mut [f32],
                    dn: &mut [f32],
                    inn: &mut [f32]|
         -> (f32, i32) {
            if frozen[i] && !probing[i] {
                mn.copy_from_slice(mo);
                dn.copy_from_slice(dd);
                inn.copy_from_slice(io);
                return (cur[i], cur_d[i]);
            }
            let own = nxt[i] as usize;
            // band column j consumes sequence position `consumed`; past the
            // core's bounds there is no position, so neither a match nor an
            // insertion may take it (a deletion consumes nothing)
            let mut out = [false; 64];
            debug_assert!(b <= out.len());
            for j in 0..b {
                let consumed = cc + j as i64 - w as i64;
                if consumed < 0 {
                    mn[j] = NEG;
                    continue;
                }
                let base = base_at(lib, core, start_off + dir * consumed);
                if base == BASE_OUT {
                    out[j] = true;
                    mn[j] = NEG;
                    continue;
                }
                let sub = if base < 4 { table[own][base as usize] } else { -2.0 };
                let prev = mo[j].max(dd[j]).max(io[j]);
                mn[j] = prev + sub;
            }
            for j in 0..b - 1 {
                dn[j] = (mo[j + 1].max(io[j + 1]) - p.gap_open).max(dd[j + 1] - p.gap_ext);
            }
            dn[b - 1] = NEG;
            inn[0] = NEG;
            for j in 1..b {
                inn[j] = if out[j] { NEG } else { (mn[j - 1].max(dn[j - 1]) - p.gap_open).max(inn[j - 1] - p.gap_ext) };
            }
            let mut best = NEG;
            let mut bj = w;
            for j in 0..b {
                let v = mn[j].max(dn[j]).max(inn[j]);
                if v > best {
                    best = v;
                    bj = j;
                }
            }
            if best <= NEG {
                // no reachable cell (the core has run off its bounds): the
                // consumption stays where it was instead of following the
                // diagonal one base a column (chr1's end: 592 reported for
                // 577 bases of room)
                return (NEG, cur_d[i]);
            }
            (best, (cc + 1 + bj as i64 - w as i64) as i32)
        };
        let mut out = vec![(0f32, 0i32); n];
        if parallel {
            out.par_iter_mut()
                .enumerate()
                .zip(mn.par_chunks_mut(b))
                .zip(dn.par_chunks_mut(b))
                .zip(inn.par_chunks_mut(b))
                .for_each(|((((i, o), mnc), dnc), innc)| {
                    *o = step(i, &cores[i], &m[i * b..(i + 1) * b], &d[i * b..(i + 1) * b], &ins[i * b..(i + 1) * b], mnc, dnc, innc);
                });
        } else {
            for i in 0..n {
                let (mnc, dnc, innc) = (&mut mn[i * b..(i + 1) * b], &mut dn[i * b..(i + 1) * b], &mut inn[i * b..(i + 1) * b]);
                out[i] = step(i, &cores[i], &m[i * b..(i + 1) * b], &d[i * b..(i + 1) * b], &ins[i * b..(i + 1) * b], mnc, dnc, innc);
            }
        }
        std::mem::swap(&mut m, &mut mn);
        std::mem::swap(&mut d, &mut dn);
        std::mem::swap(&mut ins, &mut inn);
        for i in 0..n {
            cur[i] = out[i].0;
            cur_d[i] = out[i].1;
            if !frozen[i] {
                rep_d[i] = cur_d[i];
            }
        }
        c += 1;
        total += colinfo;
        for i in 0..n {
            if frozen[i] {
                continue;
            }
            if cur[i] > peak[i] {
                peak[i] = cur[i];
                peak_c[i] = c;
            }
            if cur[i] < peak[i] - p.cap_bits {
                frozen[i] = true;
                frozen_at[i] = c;
                freeze_score[i] = cur[i];
            }
        }
        if p.window > 0 {
            hist.push_back(cur.clone());
            if hist.len() > p.window {
                let old = hist.pop_front().unwrap();
                for i in 0..n {
                    if !frozen[i] && cur[i] - old[i] <= 0.0 {
                        frozen[i] = true;
                        frozen_at[i] = c;
                        freeze_score[i] = cur[i];
                    }
                }
            }
        }
        if p.bg_window > 0 && n_act > 0.0 {
            let f = column_smoothed(&counts, n_act, p);
            for x in 0..4 {
                bg_sum[x] += f[x];
            }
            bg_hist.push_back(f);
            if bg_hist.len() > p.bg_window {
                let old = bg_hist.pop_front().unwrap();
                for x in 0..4 {
                    bg_sum[x] = (bg_sum[x] - old[x]).max(0.0);
                }
            }
        }
        let plur = (0..4).max_by(|&a, &x| counts[a].partial_cmp(&counts[x]).unwrap()).unwrap();
        cons.push(if n_act > 0.0 { plur as u8 } else { SYM_N });
        if total > best_total {
            best_total = total;
            best_c = c;
            best_d.copy_from_slice(&rep_d);
        } else if (c - best_c) as i32 > p.when_to_stop {
            break;
        }
        if p.burst_window > 0 {
            d_hist.push_back(rep_d.clone());
            if d_hist.len() > d_depth {
                d_hist.pop_front();
            }
            if let Some((_, _, probe_end)) = tentative {
                // the cohort under probe: revived if it regained score against
                // the consensus the survivors went on to build
                let cohort = probing.iter().filter(|x| **x).count();
                let revived = (0..n).filter(|&i| probing[i] && cur[i] >= freeze_score[i] + REVIVE_BITS).count();
                if revived * 2 >= cohort {
                    // a fade, not an end: the cohort rejoins
                    for i in 0..n {
                        if probing[i] {
                            probing[i] = false;
                            frozen[i] = false;
                            frozen_at[i] = usize::MAX;
                            rep_d[i] = cur_d[i];
                            if cur[i] > peak[i] {
                                peak[i] = cur[i];
                                peak_c[i] = c;
                            }
                        }
                    }
                    tentative = None;
                    newly_frozen.clear();
                    supporters_at_start.clear();
                } else if c >= probe_end {
                    burst = true;
                    break;
                }
            } else {
                let newly = (0..n).filter(|&i| frozen_at[i] == c && peak[i] >= p.burst_min_bits).count();
                newly_frozen.push_back(newly);
                supporters_at_start.push_back(supporters_now);
                if newly_frozen.len() > p.burst_window {
                    newly_frozen.pop_front();
                    supporters_at_start.pop_front();
                }
                if newly_frozen.len() == p.burst_window {
                    let cohort: usize = newly_frozen.iter().sum();
                    let base = *supporters_at_start.front().unwrap();
                    if base > 0 && cohort >= BURST_MIN_CORES && cohort as f32 >= p.burst_frac * base as f32 {
                        // tentative boundary: where the cohort's scores peaked,
                        // never past the cumulative information peak
                        let mut pcs = Vec::new();
                        for i in 0..n {
                            if frozen_at[i] != usize::MAX && frozen_at[i] + p.burst_window > c && peak[i] >= p.burst_min_bits {
                                pcs.push(peak_c[i]);
                                probing[i] = true;
                            }
                        }
                        pcs.sort_unstable();
                        let bc = pcs[pcs.len() / 2].min(best_c);
                        let oldest = c + 1 - d_hist.len();
                        let snap = if bc == 0 {
                            vec![0i32; n]
                        } else {
                            let bc = bc.max(oldest);
                            d_hist[bc - oldest].clone()
                        };
                        let bc = if bc == 0 { 0 } else { bc.max(oldest) };
                        tentative = Some((bc, snap, c + p.burst_probe));
                    }
                }
            }
        }
    }
    if let Some((bc, snap, _)) = tentative {
        // the extension ended with a boundary still under probe and the
        // cohort never revived: it stands
        best_c = bc;
        best_d = snap;
        burst = true;
    }
    cons.truncate(best_c);
    DirResult { cons, consumed: best_d, active: active_now, peak, burst }
}

/// Extend `cores` (all sharing a core word of `k` bases, oriented) in both
/// directions. Sets each core's `left_extension_len`, `right_extension_len`
/// and `score` (bits x `score_scale`, rounded).
pub fn profile_extend(
    cores: &mut [CoreAlignment],
    lib: &SequenceLibrary,
    k: usize,
    params: &ProfileParams,
) -> ProfileExtension {
    let r = extend_dir(cores, lib, k as i64, 1, params);
    let mut l = extend_dir(cores, lib, -1, -1, params);
    l.cons.reverse();
    let mut core_bits = Vec::with_capacity(cores.len());
    for (i, core) in cores.iter_mut().enumerate() {
        core.left_extension_len = l.consumed[i];
        core.right_extension_len = r.consumed[i];
        let bits = l.peak[i] + r.peak[i];
        core.score = (bits * params.score_scale).round() as i32;
        core_bits.push(bits);
    }
    ProfileExtension {
        left: l.cons,
        right: r.cons,
        left_active: l.active,
        right_active: r.active,
        left_burst: l.burst,
        right_burst: r.burst,
        core_bits,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::library::CoreBoundFlag;

    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    /// What follows each planted copy, in element orientation.
    #[derive(Clone, Copy)]
    enum Flank {
        Random,
        /// A second 200 bp element right after copies `i % 10 < 3`.
        NeighbourIn3of10,
        /// 150 A's after every copy.
        PolyA,
        /// The element continues for 150 bp at 45% divergence, then 100 bp
        /// at the copy's own divergence again (a poorly conserved stretch).
        Fade,
    }

    /// Genome with `n` copies of a 300 bp element at `div` divergence, each
    /// flanked by random sequence; cores on the element's bases 100..114.
    fn planted(n: usize, div: f64, seed: u64) -> (SequenceLibrary, Vec<CoreAlignment>) {
        planted_with(n, div, seed, Flank::Random)
    }

    fn planted_with(n: usize, div: f64, seed: u64, flank: Flank) -> (SequenceLibrary, Vec<CoreAlignment>) {
        let mut s = seed;
        let elem: Vec<u8> = (0..300).map(|_| (rng(&mut s) % 4) as u8).collect();
        let elem2: Vec<u8> = (0..200).map(|_| (rng(&mut s) % 4) as u8).collect();
        let elem3: Vec<u8> = (0..250).map(|_| (rng(&mut s) % 4) as u8).collect();
        let mut seq: Vec<u8> = Vec::new();
        let mut cores = Vec::new();
        for i in 0..n {
            for _ in 0..500 {
                seq.push((rng(&mut s) % 4) as u8);
            }
            let start = seq.len();
            let flip = i % 3 == 2;
            let mut block: Vec<u8> = elem
                .iter()
                .map(|&b| if (rng(&mut s) % 1000) as f64 / 1000.0 < div { (b + 1 + (rng(&mut s) % 3) as u8) % 4 } else { b })
                .collect();
            // keep the core word exact so all copies share it
            block[100..114].copy_from_slice(&elem[100..114]);
            match flank {
                Flank::Random => {}
                Flank::NeighbourIn3of10 => {
                    if i % 10 < 3 {
                        block.extend_from_slice(&elem2);
                    }
                }
                Flank::PolyA => block.extend(std::iter::repeat(0u8).take(150)),
                Flank::Fade => {
                    for (k, &b) in elem3.iter().enumerate() {
                        let dv = if k < 150 { 0.45 } else { div };
                        block.push(if (rng(&mut s) % 1000) as f64 / 1000.0 < dv { (b + 1 + (rng(&mut s) % 3) as u8) % 4 } else { b });
                    }
                }
            }
            let block_len = block.len();
            if flip {
                block = block.iter().rev().map(|&b| 3 - b).collect();
            }
            seq.extend_from_slice(&block);
            let (lpos, rpos) = if flip {
                // element base 100 sits at block index block_len-1-100 once flipped
                let base100 = start + (block_len - 1 - 100);
                (base100 as u64, (base100 - 13) as u64)
            } else {
                ((start + 100) as u64, (start + 113) as u64)
            };
            cores.push(CoreAlignment {
                seq_idx: 0,
                left_seq_pos: lpos,
                right_seq_pos: rpos,
                left_extendable: true,
                right_extendable: true,
                lower_seq_bound: 0,
                upper_seq_bound: 0, // set below
                lower_seq_bound_flag: CoreBoundFlag::SeqBoundary,
                upper_seq_bound_flag: CoreBoundFlag::SeqBoundary,
                left_extension_len: 0,
                right_extension_len: 0,
                score: 0,
                orient: flip,
            });
        }
        for _ in 0..500 {
            seq.push((rng(&mut s) % 4) as u8);
        }
        let len = seq.len() as u64;
        for c in cores.iter_mut() {
            c.upper_seq_bound = len - 1;
        }
        let lib = SequenceLibrary { sequence: seq, identifiers: vec!["t".into()], boundaries: vec![len], offsets: vec![0] };
        (lib, cores)
    }

    /// A copy whose sequence ends a few bases past its core must not report
    /// more consumed bases than it has. Out-of-bounds positions used to
    /// score like masked bases (-2 a column), cheaper than a gap, so the
    /// band walked past the sequence end and the reported extension
    /// followed it (hs1 chr1: telomeric cores 18 bp from the chromosome
    /// start reported left extensions of 21; the pipeline's extent
    /// arithmetic wrapped).
    #[test]
    fn a_core_at_the_sequence_end_reports_no_more_than_its_room() {
        for &(seed, room) in &[(7u64, 40u64), (11, 3), (13, 0), (17, 12), (19, 120)] {
            let (mut lib, mut cores) = planted(61, 0.10, seed);
            // the last copy (index 60, forward) now runs off the end of
            // the sequence `room` bases past its core
            let last = cores[60].clone();
            assert!(!last.orient);
            let len = last.right_seq_pos + 1 + room;
            lib.sequence.truncate(len as usize);
            lib.boundaries = vec![len];
            for c in cores.iter_mut() {
                c.upper_seq_bound = len - 1;
            }
            let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
            assert!(ext.right.len() > 150, "seed {seed}: right consensus only {} bp", ext.right.len());
            let others = cores[..60].iter().filter(|c| c.right_extension_len > room as i32 + 50).count();
            assert!(others >= 30, "seed {seed}: only {others} other copies extended past the room");
            let ext = cores[60].right_extension_len as u64;
            assert!(ext <= room, "seed {seed}: the end copy reports {ext} consumed with {room} bases of room");
        }
    }

    #[test]
    fn recovers_a_planted_element_at_20_percent_divergence() {
        let (lib, mut cores) = planted(40, 0.20, 99);
        let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
        // element = 100 bp left of the core, 186 bp right of it
        assert!((85..=110).contains(&ext.left.len()), "left {}", ext.left.len());
        assert!((170..=200).contains(&ext.right.len()), "right {}", ext.right.len());
        assert!(cores.iter().filter(|c| c.right_extension_len > 150).count() >= 30);
    }

    #[test]
    fn random_cores_do_not_extend() {
        let (lib, mut cores) = planted(40, 0.0, 7);
        // move every core to a random locus outside its element
        let mut s = 12345u64;
        for c in cores.iter_mut() {
            let p = 20 + (rng(&mut s) % (lib.sequence.len() as u64 - 400)) as u64;
            c.left_seq_pos = p;
            c.right_seq_pos = p + 13;
            c.orient = false;
        }
        let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
        assert!(ext.left.len() <= 3 && ext.right.len() <= 3, "{} / {}", ext.left.len(), ext.right.len());
    }

    #[test]
    fn a_mixed_seed_family_is_not_a_boundary() {
        // 25 of 40 cores sit at random loci (a cloud seed's hangers-on) and
        // freeze at once; the 15 that share the element must still get it
        let (lib, mut cores) = planted(40, 0.05, 23);
        let mut s = 777u64;
        for (i, c) in cores.iter_mut().enumerate() {
            if i % 8 < 5 {
                let p = 20 + (rng(&mut s) % (lib.sequence.len() as u64 - 400)) as u64;
                c.left_seq_pos = p;
                c.right_seq_pos = p + 13;
                c.orient = false;
            }
        }
        let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
        assert!((85..=110).contains(&ext.left.len()), "left {}", ext.left.len());
        assert!((170..=200).contains(&ext.right.len()), "right {}", ext.right.len());
    }

    #[test]
    fn a_poorly_conserved_stretch_is_a_fade_not_a_boundary() {
        // element = 300 bp + 150 bp at 45% divergence + 100 bp conserved:
        // the cohort that freezes in the middle stretch revives after it
        let (lib, mut cores) = planted_with(40, 0.05, 41, Flank::Fade);
        let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
        // 186 bp of element right of the core, then 150 + 100
        assert!(ext.right.len() >= 400, "right {} (burst {})", ext.right.len(), ext.right_burst);
    }

    #[test]
    fn a_shared_neighbour_in_a_minority_does_not_carry_the_extension() {
        // 12 of 40 copies have the same 200 bp element right after them:
        // above the 20% floor, so only the synchronized-freeze rule stops it
        let (lib, mut cores) = planted_with(40, 0.05, 31, Flank::NeighbourIn3of10);
        let off = ProfileParams { burst_window: 0, ..ProfileParams::default() };
        let ext_off = profile_extend(&mut cores, &lib, 14, &off);
        assert!(ext_off.right.len() >= 300, "without the burst rule the neighbour is voted in: {}", ext_off.right.len());
        let ext = profile_extend(&mut cores, &lib, 14, &ProfileParams::default());
        assert!(ext.right_burst, "burst rule should have fired");
        assert!((170..=215).contains(&ext.right.len()), "right {}", ext.right.len());
        assert!((85..=110).contains(&ext.left.len()), "left {}", ext.left.len());
    }

    #[test]
    fn a_homopolymer_tail_is_not_element() {
        // every copy is followed by 150 A's: the columns agree, but no better
        // than the local composition predicts
        let (lib, mut cores) = planted_with(40, 0.05, 17, Flank::PolyA);
        let off = ProfileParams { bg_window: 0, burst_window: 0, ..ProfileParams::default() };
        let ext_off = profile_extend(&mut cores, &lib, 14, &off);
        assert!(ext_off.right.len() >= 300, "with a fixed background the tail is element: {}", ext_off.right.len());
        let on = ProfileParams { burst_window: 0, ..ProfileParams::default() };
        let ext = profile_extend(&mut cores, &lib, 14, &on);
        assert!((170..=225).contains(&ext.right.len()), "right {}", ext.right.len());
    }
}
