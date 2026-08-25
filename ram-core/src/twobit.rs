//! UCSC 2bit reading, re-exported from `aln_core::twobit`.
//!
//! This crate's implementation became the canonical one in dfam-lib — it was
//! the only one of the three that rejected inverted ranges instead of
//! underflowing — so what remains here is a re-export. Existing
//! `ram_core::twobit::TwoBitReader` paths keep working.

pub use aln_core::twobit::*;
