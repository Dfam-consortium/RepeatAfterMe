//! Glocal (seed-anchored) consensus-vs-genome alignment: a port of
//! `cons_seed_extend` and `complexityAdjust` from `bnw_extend.c`.
//!
//! Given a consensus (query) with a seed at a known position and the seed's
//! occurrence in the sequence library, extend the alignment left and right of
//! the seed with the same banded 2-state DP row kernel used by the extension
//! engine. Both extensions start from score 0 (the seed itself is not
//! scored), so the result is the sum of the two directional maxima, and the
//! subject start/end fall out of the best cells' band offsets without any
//! traceback. A caller locates a shared word between the consensus and the
//! library, extends from it, and keeps the hit if the combined score clears
//! a threshold.
//!
//! PORT NOTE (bug fixed): the C right-hand extension terminates on a
//! hard-coded score drop of 180 (`bnw_extend.c:430`) even though `y_drop`
//! was parameterized — only the left-hand check used the parameter. Here
//! `y_drop` governs both directions; a negative value disables the check and
//! extends until sequence or query exhaustion (this matches the C golden
//! tests, which run with `y_drop = -1`).

use crate::engine::{compute_nw_row, Direction, ScoreArena};
use crate::library::{CoreAlignment, CoreBoundFlag, SequenceLibrary};
use crate::matrix::ScoringSystem;
use aln_coord::Span;

/// Result of one glocal search (C `struct glocalSearchResult`).
///
/// `query` is on the query slice; `subj` is on the concatenated library
/// sequence. Both always contain the seed, so neither is empty.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GlocalResult {
    pub score: i32,
    pub query: Span,
    pub subj: Span,
}

/// Locate the subsequence containing `pos` (first boundary beyond it).
fn seq_idx_for_pos(lib: &SequenceLibrary, pos: u64) -> usize {
    lib.boundaries
        .iter()
        .position(|&b| b > pos)
        .unwrap_or(lib.boundaries.len().saturating_sub(1))
}

/// Extend a seed match between `query_seq` (encoded, plus-strand) and its
/// occurrence at `subj_seed_start..=subj_seed_end` in the library, in both
/// directions. `orient` refers to the subject strand. The `arena` must have
/// been created for at least one core at this `bandwidth`.
///
/// `y_drop`: stop a directional extension once the row score falls more than
/// this far below that direction's running maximum; negative disables.
/// `use_complexity_adjust`: apply Phil Green's composition adjustment to the
/// combined score.
#[allow(clippy::too_many_arguments)]
pub fn cons_seed_extend(
    query_seq: &[u8],
    query_seed_start: usize,
    query_seed_end: usize,
    lib: &SequenceLibrary,
    subj_seed_start: u64,
    subj_seed_end: u64,
    orient: bool,
    arena: &mut ScoreArena,
    scoring: &ScoringSystem,
    bandwidth: i32,
    use_complexity_adjust: bool,
    y_drop: i32,
) -> GlocalResult {
    let seq_idx = seq_idx_for_pos(lib, subj_seed_start);
    let lower_seq_bound = lib.lower_bound(seq_idx);
    let upper_seq_bound = lib.boundaries[seq_idx] - 1;

    // Synthetic single core anchored on the seed occurrence. Left/right are
    // logical (repeat-relative): swapped for a reverse-strand subject.
    let core = CoreAlignment {
        seq_idx,
        left_seq_pos: if orient { subj_seed_end } else { subj_seed_start },
        right_seq_pos: if orient { subj_seed_start } else { subj_seed_end },
        left_extendable: true,
        right_extendable: true,
        lower_seq_bound,
        upper_seq_bound,
        lower_seq_bound_flag: CoreBoundFlag::SeqBoundary,
        upper_seq_bound_flag: CoreBoundFlag::SeqBoundary,
        left_extension_len: 0,
        right_extension_len: 0,
        score: 0,
        orient,
    };

    // The DP bounds exclude the seed itself so each extension leg sees
    // standard NW edge conditions on its own side of the seed. A bound that
    // would go below zero wraps like the C uint64_t and is caught by the
    // kernel's range checks.
    let run_leg = |arena: &mut ScoreArena, direction: Direction| -> (i32, i32, i32) {
        let (lw_bnd, up_bnd) = match (direction, orient) {
            (Direction::Left, false) => (lower_seq_bound, subj_seed_start.wrapping_sub(1)),
            (Direction::Left, true) => (subj_seed_end.wrapping_add(1), upper_seq_bound),
            (Direction::Right, false) => (subj_seed_end.wrapping_add(1), upper_seq_bound),
            (Direction::Right, true) => (lower_seq_bound, subj_seed_start.wrapping_sub(1)),
        };
        let rows: usize = match direction {
            Direction::Left => query_seed_start,
            Direction::Right => query_seq.len() - (query_seed_end + 1),
        };

        arena.init_boundary_row(0, bandwidth, scoring);

        let mut max_score = 0i32;
        let mut seq_max_idx = -1i32;
        let mut cons_max_idx = -1i32;
        for row_idx in 0..rows as i32 {
            let query_base = match direction {
                Direction::Left => query_seq[query_seed_start - row_idx as usize - 1],
                Direction::Right => query_seq[query_seed_end + row_idx as usize + 1],
            };
            let (row_score, row_seq_idx) = compute_nw_row(
                direction,
                row_idx,
                0,
                query_base,
                &core,
                arena,
                lw_bnd,
                up_bnd,
                &lib.sequence,
                scoring,
                bandwidth,
            );
            // C order: the drop check precedes the max update.
            if y_drop >= 0 && max_score - row_score > y_drop {
                break;
            }
            if row_score > max_score {
                max_score = row_score;
                seq_max_idx = row_seq_idx;
                cons_max_idx = row_idx;
            }
        }
        (max_score, seq_max_idx, cons_max_idx)
    };

    let (left_score, seq_left_idx, cons_left_idx) = run_leg(arena, Direction::Left);
    let (right_score, seq_right_idx, cons_right_idx) = run_leg(arena, Direction::Right);

    // Convert the best-cell indices (bases consumed minus one; -1 = none)
    // into spans. The seed is inclusive on both ends, so the subject span
    // ends one past `subj_seed_end`.
    let left_bp = (seq_left_idx + 1) as u64;
    let right_bp = (seq_right_idx + 1) as u64;
    let (lo_bp, hi_bp) = if orient { (right_bp, left_bp) } else { (left_bp, right_bp) };
    let subj = Span::new(subj_seed_start - lo_bp, subj_seed_end + 1 + hi_bp)
        .expect("the seed keeps start below end");
    let query = Span::new(
        (query_seed_start - (cons_left_idx + 1) as usize) as u64,
        (query_seed_end + 1 + (cons_right_idx + 1) as usize) as u64,
    )
    .expect("the seed keeps start below end");

    let mut combined_score = left_score + right_score;
    if use_complexity_adjust {
        let mut comp = [0i64; 4];
        for &base in &lib.sequence[subj.range_usize()] {
            // Masked (>3) and N bases are excluded from the composition.
            if base <= 3 {
                comp[base as usize] += 1;
            }
        }
        combined_score = complexity_adjust(
            combined_score,
            &comp,
            scoring.lambda,
            &scoring.bg_freqs,
        );
    }

    GlocalResult {
        score: combined_score,
        query,
        subj,
    }
}

/// Phil Green's cross_match complexity score adjustment: rescale a raw score
/// by the entropy of the aligned subject composition relative to the matrix
/// background frequencies (port of C `complexityAdjust`, including the
/// cross_match `+ .999` rounding and the floor at 0).
pub fn complexity_adjust(score: i32, comp: &[i64; 4], lambda: f64, bg_freqs: &[f64; 4]) -> i32 {
    let mut t_factor = 0.0f64;
    let mut t_sum = 0.0f64;
    let mut t_counts = 0.0f64;
    for i in 0..4 {
        if comp[i] != 0 {
            let c = comp[i] as f64;
            t_factor += c * c.ln();
            t_sum += c * bg_freqs[i].ln();
            t_counts += c;
        }
    }
    if t_counts != 0.0 {
        t_factor -= t_counts * t_counts.ln();
    }
    t_sum -= t_factor;

    let adj_score = (score as f64 + t_sum / lambda + 0.999) as i32;
    adj_score.max(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::char_to_num;
    use crate::matrix::ScoringSystem;

    fn encode(s: &str) -> Vec<u8> {
        s.bytes().map(|c| char_to_num(c).unwrap()).collect()
    }

    /// Fixture from the C `MU_TEST(test_cons_seed_extend)` in bnw_extend.c:
    /// an Alu-derived consensus seeded against three embedded Alu copies
    /// (two forward, one reverse), golden scores/coordinates from the C
    /// suite. All three run with y_drop disabled, so the values are valid
    /// for both the C behavior and the fixed y-drop handling here.
    #[test]
    fn cons_seed_extend_matches_c_golden_values() {
        let cons_seq = encode(concat!(
            "TAAAAATAAAAATAGGCTTGGGCGCGGTGGCTCACGCCTGTAATCCCAGC",
            "ACTTAGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGAC",
            "GGCAGGAGATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGAT",
            "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA",
            "AAAAAAAAAAAAAAATAA",
        ));
        assert_eq!(cons_seq.len(), 218);

        let alu_seqs = encode(concat!(
            "ATATGGATAGCTAGCGTGCGACGTACGGCGATTGGTATATGAGCGATATC",
            "GGAGGAAATATTTATAGTTGGGCGCGGTGGCTCACGCCTGTAATCCCAGC",
            "ACTTAGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGAC",
            "GGCAGGAGATCCCTTGAACCCGGCGGAGGTTGCAGTGCATTAGCCGAGAT",
            "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA",
            "AAAAATTTAAAAAGAAAA",
            "GTAGTAGCATGCGAGCGTTAGCGATATATTAA",
            "ATATGAGATGCGGGCTAATGCATATATTTTATGCGCGAGCACAAACGATT",
            "GATTAGCGCGCATCACGATTAGGAGGATATGAATTATGCGGCGAGAAGAT",
            "TAAAAATAAAAATAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA",
            "CTTTGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGACC",
            "AGCCTGGCCAACATGTGAAACCCCGTGCTCTACTAAAAATACAAAAATTA",
            "GCCGGGCATGGTGGCACGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAG",
            "GCAGGAGTAATCGCTTGAACCCGGGAGGCGAGGTTGCAGTGAGCCGAGAT",
            "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA",
            "AAAAAAATAAAAAAAAA",
            "GATATCTACTATCATCTATAAATCTATTATTTA",
            "ACGATTGGGCGGCGCTTTATTGCAGTAGTATCTAGCTCTCTCAGGCAAAC",
            "TTATTTTATTTTTTTTTTTTGAGACGCAGTCTCGCTCTGTTGCCCAGGCT",
            "GGAGCAGTGGCGCGATCTCGGCACACTGCAACCTAGCCGCCTCCCGGGTT",
            "CAAGCGATCTCCTGCCGTCTCGAACTCCTGACCTCTAGTGATCCGCGCGC",
            "AATCGGCCTCCCTAAGTGCTGGGATTTCAGGCGTGAGCCACCGCGCCCAA",
            "GTATTATTATTTTTA",
            "AGTAGAGGTGCGGATATGCGGAGTAGCGGATATAG",
        ));
        assert_eq!(alu_seqs.len(), 1050);

        let lib = SequenceLibrary {
            sequence: alu_seqs,
            identifiers: vec!["seq1".to_string()],
            boundaries: vec![1050],
            offsets: vec![0],
        };

        let bandwidth = 5;
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let mut arena = ScoreArena::new(1, bandwidth);

        // Forward copy 1, no complexity adjustment.
        let r = cons_seed_extend(
            &cons_seq, 50, 60, &lib, 100, 110, false, &mut arena, &scoring, bandwidth, false, -1,
        );
        assert_eq!(
            r,
            GlocalResult {
                score: 1573,
                // C golden values are 0-based inclusive: query 17..=214,
                // subject 67..=267.
                query: Span::new(17, 215).unwrap(),
                subj: Span::new(67, 268).unwrap(),
            }
        );

        // Forward copy 2, complexity-adjusted.
        let r = cons_seed_extend(
            &cons_seq, 50, 60, &lib, 449, 459, false, &mut arena, &scoring, bandwidth, true, -1,
        );
        assert_eq!(
            r,
            GlocalResult {
                score: 743,
                query: Span::new(0, 100).unwrap(),
                subj: Span::new(400, 499).unwrap(),
            }
        );

        // Reverse-strand copy, complexity-adjusted.
        let r = cons_seed_extend(
            &cons_seq, 50, 60, &lib, 956, 966, true, &mut arena, &scoring, bandwidth, true, -1,
        );
        assert_eq!(
            r,
            GlocalResult {
                score: 1462,
                query: Span::new(0, 213).unwrap(),
                subj: Span::new(803, 1015).unwrap(),
            }
        );
    }

    /// A y-drop must terminate BOTH directions (the C code only honored it on
    /// the left; the right side used a literal 180). With a tiny y_drop and a
    /// long non-matching right flank, the right extension must stop early
    /// instead of scanning the full query.
    #[test]
    fn y_drop_applies_to_right_extension() {
        // Query: 10bp seed then 40bp that mismatch the subject completely.
        let query = encode(&("ACGTACGTAC".to_string() + &"A".repeat(40)));
        let subject = encode(&("ACGTACGTAC".to_string() + &"C".repeat(40)));
        let lib = SequenceLibrary {
            sequence: subject,
            identifiers: vec!["s".to_string()],
            boundaries: vec![50],
            offsets: vec![0],
        };
        let scoring = ScoringSystem::by_name("20p43g").unwrap();
        let bandwidth = 5;
        let mut arena = ScoreArena::new(1, bandwidth);

        let strict = cons_seed_extend(
            &query, 0, 9, &lib, 0, 9, false, &mut arena, &scoring, bandwidth, false, 30,
        );
        // Nothing scores above 0 to the right, and the drop check must fire
        // long before the query is exhausted, leaving the zero-length result.
        assert_eq!(strict.query.end(), 10);
        assert_eq!(strict.subj.end(), 10);
        assert_eq!(strict.score, 0);
    }
}
