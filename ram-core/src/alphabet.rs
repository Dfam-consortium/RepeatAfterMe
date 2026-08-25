//! Numeric DNA encoding shared with the C tool (`sequence.h`).
//!
//! | code | meaning              |
//! |------|----------------------|
//! | 0–3  | A, C, G, T           |
//! | 4–7  | a, c, g, t (masked)  |
//! | 99   | N / any IUPAC code   |
//!
//! The masked codes 4–7 let a scoring matrix treat soft-masked sequence like
//! `N` without losing the underlying base. The 2bit loading path in this crate
//! uppercases everything (as the C `loadSequenceSubsetMinimal` does via
//! `toUpperN`), so masked codes only appear if a caller constructs a library
//! by hand.

pub const SYM_A: u8 = 0;
pub const SYM_C: u8 = 1;
pub const SYM_G: u8 = 2;
pub const SYM_T: u8 = 3;
pub const SYM_MASKED_A: u8 = 4;
pub const SYM_MASKED_C: u8 = 5;
pub const SYM_MASKED_G: u8 = 6;
pub const SYM_MASKED_T: u8 = 7;
pub const SYM_N: u8 = 99;

/// Encode one ASCII base. Any IUPAC ambiguity code (upper or lower case),
/// `N`, `n`, or `x` becomes [`SYM_N`]; anything else is an error (the C tool
/// exits with the same message).
pub fn char_to_num(c: u8) -> crate::Result<u8> {
    Ok(match c {
        b'A' => SYM_A,
        b'C' => SYM_C,
        b'G' => SYM_G,
        b'T' => SYM_T,
        b'a' => SYM_MASKED_A,
        b'c' => SYM_MASKED_C,
        b'g' => SYM_MASKED_G,
        b't' => SYM_MASKED_T,
        b'N' | b'n' | b'x' => SYM_N,
        c if is_iupac(c) => SYM_N,
        _ => {
            return Err(crate::Error::Format(format!(
                "ERROR: Cannot interpret input symbol '{}' [{}] as a DNA base.",
                c as char, c
            )))
        }
    })
}

fn is_iupac(c: u8) -> bool {
    matches!(
        c.to_ascii_uppercase(),
        b'R' | b'Y' | b'M' | b'K' | b'W' | b'S' | b'B' | b'D' | b'H' | b'V'
    )
}

/// Decode to ASCII; anything not 0–7 renders as `N`.
pub fn num_to_char(z: u8) -> u8 {
    match z {
        SYM_A => b'A',
        SYM_C => b'C',
        SYM_G => b'G',
        SYM_T => b'T',
        SYM_MASKED_A => b'a',
        SYM_MASKED_C => b'c',
        SYM_MASKED_G => b'g',
        SYM_MASKED_T => b't',
        _ => b'N',
    }
}

/// Complement an encoded base, preserving masking; non-bases stay [`SYM_N`].
pub fn compl(c: u8) -> u8 {
    match c {
        SYM_A => SYM_T,
        SYM_C => SYM_G,
        SYM_G => SYM_C,
        SYM_T => SYM_A,
        SYM_MASKED_A => SYM_MASKED_T,
        SYM_MASKED_C => SYM_MASKED_G,
        SYM_MASKED_G => SYM_MASKED_C,
        SYM_MASKED_T => SYM_MASKED_A,
        _ => SYM_N,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_the_eight_real_codes() {
        for c in [b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't'] {
            assert_eq!(num_to_char(char_to_num(c).unwrap()), c);
        }
    }

    #[test]
    fn iupac_codes_become_n() {
        for c in *b"RYMKWSBDHVrymkwsbdhvNnx" {
            assert_eq!(char_to_num(c).unwrap(), SYM_N);
        }
    }

    #[test]
    fn complement_preserves_masking_and_n() {
        assert_eq!(compl(SYM_A), SYM_T);
        assert_eq!(compl(SYM_MASKED_C), SYM_MASKED_G);
        assert_eq!(compl(SYM_N), SYM_N);
    }

    #[test]
    fn garbage_symbols_are_rejected() {
        assert!(char_to_num(b'-').is_err());
        assert!(char_to_num(b'*').is_err());
    }
}
