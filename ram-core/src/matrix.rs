//! Scoring systems: the four internally-coded RepeatMasker-derived DNA
//! matrices (divergence-tuned, 43% GC) and the original RepeatScout
//! match/mismatch/linear-gap scheme (`score_system.c`).
//!
//! | name    | tuned div | gap open | gap ext |
//! |---------|-----------|----------|---------|
//! | 14p43g  | 14%       | -33      | -7      |
//! | 18p43g  | 18%       | -30      | -6      |
//! | 20p43g  | 20%       | -28      | -5      |
//! | 25p43g  | 25%       | -25      | -5      |
//!
//! Matrices are asymmetric and indexed `matrix[consensus_base][sequence_base]`
//! over the numeric alphabet (0–7 plus 99 = N), exactly as in C: a 100x100
//! table so code 99 indexes directly. Masked codes (4–7) score -1 against
//! everything (the "soft-masked behaves like N" rule); in the RepeatScout
//! scheme masked-vs-unmasked scores the mismatch penalty instead.
//!
//! Lambda is derived from the matrix + background frequencies by the same
//! bisection as the C `calculateLambda` and is used only for Phil Green
//! complexity adjustment (not by the extension path itself).

use crate::alphabet::{SYM_A, SYM_MASKED_A, SYM_MASKED_T, SYM_N, SYM_T};

const MSIZE: usize = 100;

/// A substitution matrix plus affine gap penalties and derived statistics.
pub struct ScoringSystem {
    pub name: &'static str,
    /// `matrix[consensus_base][sequence_base]`, numeric alphabet.
    pub matrix: Box<[[i32; MSIZE]; MSIZE]>,
    pub gapopen: i32,
    pub gapextn: i32,
    pub lambda: f64,
    /// Background frequencies `[A, C, G, T]`.
    pub bg_freqs: [f64; 4],
}

impl ScoringSystem {
    /// Look up one of the internally-coded matrices by name.
    pub fn by_name(name: &str) -> crate::Result<ScoringSystem> {
        // (gapopen, gapextn, bg = [A,C,G,T],
        //  rows in C source order: cons A, G, C, T vs seq [A, C, G, T, N])
        let (sname, gapopen, gapextn, a_row, g_row, c_row, t_row): (
            &'static str,
            i32,
            i32,
            [i32; 5],
            [i32; 5],
            [i32; 5],
            [i32; 5],
        ) = match name {
            "14p43g" => (
                "14p43g",
                -33,
                -7,
                [9, -18, -10, -21, -1],
                [-7, -18, 11, -18, -1],
                [-18, 11, -18, -7, -1],
                [-21, -10, -18, 9, -1],
            ),
            "18p43g" => (
                "18p43g",
                -30,
                -6,
                [9, -15, -8, -18, -1],
                [-5, -16, 10, -16, -1],
                [-16, 10, -16, -5, -1],
                [-18, -8, -15, 9, -1],
            ),
            "20p43g" => (
                "20p43g",
                -28,
                -5,
                [9, -15, -8, -17, -1],
                [-4, -15, 10, -15, -1],
                [-15, 10, -15, -4, -1],
                [-17, -8, -15, 9, -1],
            ),
            "25p43g" => (
                "25p43g",
                -25,
                -5,
                [8, -13, -6, -15, -1],
                [-2, -13, 9, -13, -1],
                [-13, 9, -13, -2, -1],
                [-15, -6, -13, 8, -1],
            ),
            other => {
                return Err(crate::Error::Format(format!(
                    "Custom matrices not supported ( yet ).  {other} is not an internally coded matrix!"
                )))
            }
        };

        let mut m = Box::new([[0i32; MSIZE]; MSIZE]);
        let seq_syms = [0usize, 1, 2, 3, SYM_N as usize];
        // C fills rows in the order A, G, C, T; cons index order is A=0,C=1,G=2,T=3.
        for (cons, row) in [(0usize, a_row), (2, g_row), (1, c_row), (3, t_row)] {
            for (j, &sym) in seq_syms.iter().enumerate() {
                m[cons][sym] = row[j];
            }
        }
        for &sym in &seq_syms {
            m[SYM_N as usize][sym] = -1;
        }
        // Masked codes score -1 against everything.
        for i in SYM_MASKED_A..=SYM_MASKED_T {
            let i = i as usize;
            for j in SYM_MASKED_A as usize..=SYM_MASKED_T as usize {
                m[i][j] = -1;
                m[j][i] = -1;
            }
            for j in SYM_A as usize..=SYM_T as usize {
                m[i][j] = -1;
                m[j][i] = -1;
            }
            m[i][SYM_N as usize] = -1;
            m[SYM_N as usize][i] = -1;
        }

        let bg = [0.285, 0.215, 0.215, 0.285];
        let mut sys = ScoringSystem {
            name: sname,
            matrix: m,
            gapopen,
            gapextn,
            lambda: 0.0,
            bg_freqs: bg,
        };
        sys.lambda = calculate_lambda(&sys);
        Ok(sys)
    }

    /// The original RepeatScout scheme: uniform match/mismatch, linear gaps
    /// (`gapopen` = 0, `gapextn` = gap).
    pub fn repeat_scout(match_score: i32, mismatch: i32, gap: i32) -> ScoringSystem {
        let mut m = Box::new([[0i32; MSIZE]; MSIZE]);
        for cons in 0..4usize {
            for seq in 0..4usize {
                m[cons][seq] = if cons == seq { match_score } else { mismatch };
            }
            m[cons][SYM_N as usize] = mismatch;
            m[SYM_N as usize][cons] = mismatch;
        }
        m[SYM_N as usize][SYM_N as usize] = mismatch;
        for i in SYM_MASKED_A..=SYM_MASKED_T {
            let i = i as usize;
            for j in SYM_MASKED_A as usize..=SYM_MASKED_T as usize {
                m[i][j] = -1;
                m[j][i] = -1;
            }
            for j in SYM_A as usize..=SYM_T as usize {
                m[i][j] = mismatch;
                m[j][i] = mismatch;
            }
            m[i][SYM_N as usize] = mismatch;
            m[SYM_N as usize][i] = mismatch;
        }
        let mut sys = ScoringSystem {
            name: "repeatscout",
            matrix: m,
            gapopen: 0,
            gapextn: gap,
            lambda: 0.0,
            bg_freqs: [0.25; 4],
        };
        sys.lambda = calculate_lambda(&sys);
        sys
    }

    /// Override the per-matrix default gap penalties (C `getMatrixUsingGapPenalties`).
    pub fn with_gap_penalties(mut self, gapopen: i32, gapextn: i32) -> ScoringSystem {
        self.gapopen = gapopen;
        self.gapextn = gapextn;
        self
    }

    #[inline]
    pub fn score(&self, cons_base: u8, seq_base: u8) -> i32 {
        self.matrix[cons_base as usize][seq_base as usize]
    }
}

/// Karlin-Altschul lambda by bisection, port of the C `calculateLambda`.
/// Returns -1.0 if the background frequencies do not sum to ~1.
pub fn calculate_lambda(sys: &ScoringSystem) -> f64 {
    let sum_for = |lambda: f64| -> Option<f64> {
        let mut sum = 0.0;
        let mut check = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                sum += sys.bg_freqs[i] * sys.bg_freqs[j] * (lambda * sys.matrix[i][j] as f64).exp();
                check += sys.bg_freqs[i] * sys.bg_freqs[j];
            }
        }
        if !(0.999..=1.001).contains(&check) {
            None
        } else {
            Some(sum)
        }
    };

    let mut lambda = 0.5;
    let mut lambda_lower = 0.0;
    loop {
        match sum_for(lambda) {
            None => return -1.0,
            Some(s) if s >= 1.0 => break,
            Some(_) => {
                lambda_lower = lambda;
                lambda *= 2.0;
            }
        }
    }
    let mut lambda_upper = lambda;

    while lambda_upper - lambda_lower > 0.00001 {
        lambda = (lambda_lower + lambda_upper) / 2.0;
        match sum_for(lambda) {
            None => return -1.0,
            Some(s) if s >= 1.0 => lambda_upper = lambda,
            Some(_) => lambda_lower = lambda,
        }
    }
    lambda
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{SYM_C, SYM_G};

    #[test]
    fn the_20p43g_matrix_matches_the_c_source() {
        let s = ScoringSystem::by_name("20p43g").unwrap();
        assert_eq!(s.gapopen, -28);
        assert_eq!(s.gapextn, -5);
        assert_eq!(s.score(SYM_A, SYM_A), 9);
        assert_eq!(s.score(SYM_A, SYM_G), -8); // transition, asymmetric
        assert_eq!(s.score(SYM_G, SYM_A), -4);
        assert_eq!(s.score(SYM_C, SYM_T), -4);
        assert_eq!(s.score(SYM_T, SYM_A), -17);
        assert_eq!(s.score(SYM_A, SYM_N), -1);
        assert_eq!(s.score(SYM_N, SYM_T), -1);
        assert_eq!(s.score(SYM_MASKED_A, SYM_A), -1);
    }

    #[test]
    fn repeatscout_matrix_scores_masked_as_mismatch() {
        let s = ScoringSystem::repeat_scout(1, -1, -5);
        assert_eq!(s.score(SYM_A, SYM_A), 1);
        assert_eq!(s.score(SYM_A, SYM_C), -1);
        assert_eq!(s.score(SYM_MASKED_A, SYM_A), -1);
        assert_eq!(s.gapopen, 0);
        assert_eq!(s.gapextn, -5);
    }

    #[test]
    fn lambda_for_plus_one_minus_one_is_ln3() {
        // The C unit test pins 1/-1 at lambda ~= 1.09861 (= ln 3).
        let mut s = ScoringSystem::repeat_scout(1, -1, -5);
        // repeat_scout uses mismatch for N; the lambda calc only reads 0..4.
        s.bg_freqs = [0.25; 4];
        let l = calculate_lambda(&s);
        assert!((l - 1.09861).abs() < 0.0001, "lambda = {l}");
    }

    #[test]
    fn unknown_matrix_names_are_rejected() {
        assert!(ScoringSystem::by_name("42p99g").is_err());
    }
}
