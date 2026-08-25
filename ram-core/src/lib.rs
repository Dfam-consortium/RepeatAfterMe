//! ram-core — fit-preferred banded MSA extension (Rust port of RAMExtend).
//!
//! This crate is a semantics-preserving port of the RepeatAfterMe `RAMExtend`
//! C tool (v0.0.7): simultaneous extension of a multiple alignment anchored by
//! a set of genomic core ranges, using the RepeatScout "fit-preferred" scoring
//! model (Price, Jones & Pevzner 2005) generalized to substitution matrices
//! and affine gaps.
//!
//! Port ground rules, in priority order:
//!  1. Bit-identical results to the C binary on identical inputs — including
//!     tie-breaking order, sentinels, and integer-division quirks. Documented
//!     deviations are marked `PORT NOTE` at the site.
//!  2. Depend only on aln-core (2bit reading) and rayon (row-kernel
//!     parallelism), and nothing else outside std.
//!  3. Structure-of-code may differ (flat arenas instead of pointer forests,
//!     `Vec` instead of linked lists), behavior may not.
//!
//! Module map:
//!  - [`alphabet`]  — the numeric base encoding shared with the C tool
//!  - [`matrix`]    — scoring systems (14/18/20/25p43g + RepeatScout simple)
//!  - [`twobit`]    — minimal UCSC 2bit random-access reader
//!  - [`library`]   — sequence library + core ranges loader (ranges TSV + 2bit)
//!  - [`engine`]    — the banded DP row kernel and the extension driver
//!  - [`glocal`]    — seed-anchored glocal alignment (`cons_seed_extend`)

pub mod alphabet;
pub mod engine;
pub mod glocal;
pub mod library;
pub mod matrix;
pub mod twobit;

use std::fmt;

/// Errors surfaced by ram-core.
///
/// PORT NOTE: the C tool calls `exit(1)` with a printf at every failure site;
/// the messages here preserve the C wording where a caller might reasonably
/// grep for it.
#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    /// Input violated a format expectation (ranges TSV, 2bit, matrix name...).
    Format(String),
    /// A documented C-tool limit was exceeded (e.g. more ranges than -maxoccurrences).
    Limit(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Format(m) => write!(f, "{m}"),
            Error::Limit(m) => write!(f, "{m}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;
