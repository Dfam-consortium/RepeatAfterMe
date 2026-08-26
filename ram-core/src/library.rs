//! Sequence library and core-range loading (port of `sequence.c`'s
//! `readBEDRanges` + `loadSequenceSubsetMinimal`).
//!
//! The library holds one encoded subsequence per core range: the core itself
//! plus up to `max_flanking_bp` (= L + bandwidth) of flanking sequence on each
//! extendable side. Flanks are trimmed when a neighboring core on the same
//! input sequence is closer than the full flank:
//!
//! | boundary flag  | meaning                                             |
//! |----------------|-----------------------------------------------------|
//! | `LBoundary`    | limited by the maximum extension length (L)         |
//! | `SeqBoundary`  | limited by the input sequence start/end             |
//! | `CoreBoundary` | limited by a neighboring core                       |
//! | `ExtBoundary`  | limited by a phase-1 extension (set by the engine)  |
//!
//! Ranges file: 6 tab-separated columns, `seqid  start  end  leftExt
//! rightExt  strand`, 0-based half-open, strand `+`/`-`. The left/right
//! extendable flags refer to the edges of the core MSA (strand-relative).
//!
//! # Neighbor flank splitting and C compatibility
//!
//! When two cores extend toward each other in the same directional pass
//! (which can only happen for an opposite-strand pair), the gap between them
//! is split at the midpoint; a same-strand neighbor, or one whose facing
//! extendable flag is off, cedes the full inter-core distance.
//!
//! PORT NOTE: the C code expresses the same-strand test as
//! `s->strand == neighbor->strand`, comparing `char*` pointers from separate
//! `cloneString` allocations — always false — so the C tool midpoints
//! same-strand pairs too whenever the neighbor's facing flag is set. Passing
//! `c_compat = true` replicates that bug for differential testing; the
//! default applies the intended semantics.
//! PORT NOTE: `ceil(prev_core_dist/2)` in C is `ceil()` of an already-floored
//! integer division — i.e. plain `dist/2`, same as the `floor` branch.

use std::path::Path;

use aln_coord::Span;

use crate::alphabet::SYM_N;
use crate::twobit::TwoBitReader;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoreBoundFlag {
    /// Limited by the maximum extension length permitted (L parameter).
    LBoundary,
    /// Limited by the sequence length.
    SeqBoundary,
    /// Limited by the presence of a neighboring core.
    CoreBoundary,
    /// Limited by the first extension phase.
    ExtBoundary,
}

/// A core aligned region from which extension is performed.
///
/// The engine walks single positions, not ranges, so this struct holds
/// positions: `left_seq_pos`/`right_seq_pos` are the indices of the core's
/// first and last bases in the in-memory `SequenceLibrary::sequence` array,
/// in the *logical* left/right order of the core MSA. For a reverse-strand
/// core the left position is numerically larger than the right. Read the
/// core as a range through [`CoreAlignment::core_span`] and
/// [`CoreAlignment::extended_span`] rather than adding 1 by hand.
#[derive(Debug, Clone)]
pub struct CoreAlignment {
    pub seq_idx: usize,
    pub left_seq_pos: u64,
    pub right_seq_pos: u64,
    pub left_extendable: bool,
    pub right_extendable: bool,
    /// Lowest in-memory position this core may extend into.
    pub lower_seq_bound: u64,
    /// Highest in-memory position this core may extend into. A position,
    /// not a range end, so it is the last usable index.
    pub upper_seq_bound: u64,
    pub lower_seq_bound_flag: CoreBoundFlag,
    pub upper_seq_bound_flag: CoreBoundFlag,
    /// After extension, the length of the left extension in bp.
    pub left_extension_len: i32,
    /// After extension, the length of the right extension in bp.
    pub right_extension_len: i32,
    /// After extension, the summed left+right extension score.
    pub score: i32,
    /// True = reverse strand.
    pub orient: bool,
}

impl CoreAlignment {
    /// The core's bases in the library buffer, ascending regardless of
    /// strand.
    pub fn core_span(&self) -> Span {
        let lo = self.left_seq_pos.min(self.right_seq_pos);
        let hi = self.left_seq_pos.max(self.right_seq_pos);
        Span::new(lo, hi + 1).expect("lo <= hi")
    }

    /// The core plus whatever the engine extended on each side, in the
    /// library buffer. The extension lengths are strand-relative; a
    /// reverse-strand core's left extension grows the upper end.
    ///
    /// The engine keeps the extension inside `lower_seq_bound`, so the start
    /// cannot go below 0; the subtraction saturates rather than wrapping so
    /// that a negative extension length (which the engine never produces)
    /// cannot silently become a huge coordinate.
    pub fn extended_span(&self) -> Span {
        let core = self.core_span();
        let (lo_ext, hi_ext) = if self.orient {
            (self.right_extension_len, self.left_extension_len)
        } else {
            (self.left_extension_len, self.right_extension_len)
        };
        Span::new(
            core.start().saturating_sub(lo_ext.max(0) as u64),
            core.end() + hi_ext.max(0) as u64,
        )
        .expect("core.start <= core.end")
    }
}

/// Concatenated encoded subsequences plus per-subsequence metadata.
///
/// `boundaries[i]` is the cumulative end (exclusive) of subsequence `i` in
/// `sequence`; subsequence `i` starts at `boundaries[i-1]` (0 for the first).
/// `offsets[i]` is the position in the original input sequence where
/// subsequence `i` begins (i.e. the flanking start).
pub struct SequenceLibrary {
    pub sequence: Vec<u8>,
    pub identifiers: Vec<String>,
    pub boundaries: Vec<u64>,
    pub offsets: Vec<u64>,
}

impl SequenceLibrary {
    pub fn length(&self) -> u64 {
        self.sequence.len() as u64
    }

    /// Start of subsequence `i` in `sequence` (the C `seqLowerBound`).
    pub fn lower_bound(&self, seq_idx: usize) -> u64 {
        if seq_idx > 0 {
            self.boundaries[seq_idx - 1]
        } else {
            0
        }
    }

    /// Translate a span in the library buffer, lying within subsequence
    /// `seq_idx`, to the source sequence that subsequence was cut from.
    pub fn to_source(&self, seq_idx: usize, span: Span) -> Span {
        let shift = |p: u64| p - self.lower_bound(seq_idx) + self.offsets[seq_idx];
        Span::new(shift(span.start()), shift(span.end())).expect("shifting both ends keeps order")
    }
}

/// One parsed line of the ranges TSV. Flag fields keep their raw `atoi`
/// values because the loader distinguishes `== 1` from `== 0`.
#[derive(Debug, Clone)]
pub struct RangeRecord {
    pub name: String,
    /// The core on the source sequence. The file is 0-based half-open, so
    /// this is the columns as written.
    pub span: Span,
    pub left_flag: i32,
    pub right_flag: i32,
    /// True = reverse strand.
    pub minus: bool,
}

/// C `atoi`: optional leading whitespace and sign, digits until the first
/// non-digit, 0 on no digits.
fn atoi_like(s: &str) -> i32 {
    let s = s.trim_start();
    let mut chars = s.chars().peekable();
    let mut neg = false;
    if let Some(&c) = chars.peek() {
        if c == '+' || c == '-' {
            neg = c == '-';
            chars.next();
        }
    }
    let mut val: i64 = 0;
    for c in chars {
        match c.to_digit(10) {
            Some(d) => val = (val * 10 + d as i64).min(i32::MAX as i64),
            None => break,
        }
    }
    let val = if neg { -val } else { val };
    val.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// C `strtol(_, _, 0)`: decimal, or hex/octal by prefix; 0 on no digits.
fn strtol_like(s: &str) -> i64 {
    let t = s.trim_start();
    let (neg, t) = match t.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    let (radix, digits) = if let Some(hex) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        (16, hex)
    } else if t.starts_with('0') && t.len() > 1 {
        (8, &t[1..])
    } else {
        (10, t)
    };
    let mut val: i64 = 0;
    for c in digits.chars() {
        match c.to_digit(radix) {
            Some(d) => val = val.saturating_mul(radix as i64).saturating_add(d as i64),
            None => break,
        }
    }
    if neg {
        -val
    } else {
        val
    }
}

/// Parse the 6-column ranges TSV (C `readBEDRanges`). Blank lines are
/// skipped; a row whose 6th field does not start with `+`/`-` is an error
/// with the C tool's message.
pub fn read_ranges(path: &Path) -> crate::Result<Vec<RangeRecord>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.splitn(7, '\t').collect();
        if fields.len() < 6 || !(fields[5].starts_with('+') || fields[5].starts_with('-')) {
            return Err(crate::Error::Format(
                "Error: ranges file does not appear to be in the correct format!".to_string(),
            ));
        }
        let (start, end) = (strtol_like(fields[1]), strtol_like(fields[2]));
        let span = u64::try_from(start)
            .ok()
            .zip(u64::try_from(end).ok())
            .and_then(|(s, e)| Span::new(s, e).ok())
            .ok_or_else(|| {
                crate::Error::Format(format!(
                    "Error: range {start}-{end} for {} is not a 0-based half-open interval",
                    fields[0]
                ))
            })?;
        out.push(RangeRecord {
            name: fields[0].to_string(),
            span,
            left_flag: atoi_like(fields[3]),
            right_flag: atoi_like(fields[4]),
            minus: fields[5].starts_with('-'),
        });
    }
    Ok(out)
}

/// Port of `loadSequenceSubsetMinimal`: load one subsequence per core with
/// enough flanking sequence to support extension by `max_flanking_bp`
/// (typically L + bandwidth), trimming flanks against neighboring cores.
///
/// `c_compat` selects the C tool's buggy neighbor-splitting (see the module
/// docs); warnings about overlapping cores go to stdout, matching the C tool.
pub fn load_sequence_subset_minimal(
    twobit: &TwoBitReader,
    ranges: &[RangeRecord],
    max_flanking_bp: i64,
    c_compat: bool,
) -> crate::Result<(SequenceLibrary, Vec<CoreAlignment>)> {
    // C sorts with bedNameStartRevEndCmp: name asc, start asc, end desc.
    let mut sorted: Vec<&RangeRecord> = ranges.iter().collect();
    sorted.sort_by(|a, b| {
        a.name
            .cmp(&b.name)
            .then(a.span.start().cmp(&b.span.start()))
            .then(b.span.end().cmp(&a.span.end()))
    });

    let mut lib = SequenceLibrary {
        sequence: Vec::new(),
        identifiers: Vec::new(),
        boundaries: Vec::new(),
        offsets: Vec::new(),
    };
    let mut cores: Vec<CoreAlignment> = Vec::new();

    for (i, s) in sorted.iter().enumerate() {
        let seq_size = twobit.seq_len(&s.name).ok_or_else(|| {
            crate::Error::Format(format!("chromosome {:?} not found in 2bit file", s.name))
        })? as i64;

        // The flank arithmetic below is signed, as in the C.
        let (s_start, s_end) = (s.span.start() as i64, s.span.end() as i64);
        let mut lower_flank_len: i64 = 0;
        let mut flanking_start: i64 = s_start;
        let mut flanking_end: i64 = s_end;
        let mut lower_bound_flag = CoreBoundFlag::LBoundary;
        let mut upper_bound_flag = CoreBoundFlag::LBoundary;

        // Neighboring cores on the same input sequence (in sorted order).
        let prev = if i > 0 && sorted[i - 1].name == s.name {
            let p = sorted[i - 1];
            let (p_start, p_end) = (p.span.start() as i64, p.span.end() as i64);
            let dist = if s_start < p_end {
                eprintln!(
                    "WARNING: core sequences overlap  {}:{}-{} and previous {}:{}-{}",
                    s.name, s_start, s_end, p.name, p_start, p_end
                );
                0
            } else {
                s_start - p_end
            };
            Some((p, dist))
        } else {
            None
        };
        let next = if i + 1 < sorted.len() && sorted[i + 1].name == s.name {
            let n = sorted[i + 1];
            let (n_start, n_end) = (n.span.start() as i64, n.span.end() as i64);
            let dist = if n_start < s_end {
                eprintln!(
                    "WARNING: core sequences overlap  {}:{}-{} and next {}:{}-{}",
                    s.name, s_start, s_end, n.name, n_start, n_end
                );
                0
            } else {
                n_start - s_end
            };
            Some((n, dist))
        } else {
            None
        };

        // Full distance when the neighbor cannot extend toward us in the
        // same pass (facing flag off, or — unless in C-compat mode, see the
        // module PORT NOTE — same strand); midpoint otherwise.
        let grants_full =
            |neighbor: &RangeRecord, facing_flag: i32| -> bool {
                facing_flag == 0 || (!c_compat && neighbor.minus == s.minus)
            };
        if s.minus {
            // Reverse strand: left = end side, right = start side.
            if s.left_flag == 1 {
                match next {
                    Some((n, dist)) if dist <= max_flanking_bp => {
                        let flank = if grants_full(n, n.left_flag) { dist } else { dist / 2 };
                        flanking_end = s_end + flank;
                        upper_bound_flag = CoreBoundFlag::CoreBoundary;
                    }
                    _ if s_end + max_flanking_bp < seq_size => {
                        flanking_end = s_end + max_flanking_bp;
                        upper_bound_flag = CoreBoundFlag::LBoundary;
                    }
                    _ => {
                        flanking_end = seq_size;
                        upper_bound_flag = CoreBoundFlag::SeqBoundary;
                    }
                }
            }
            if s.right_flag == 1 {
                match prev {
                    Some((p, dist)) if dist <= max_flanking_bp => {
                        let flank = if grants_full(p, p.right_flag) { dist } else { dist / 2 };
                        flanking_start = s_start - flank;
                        lower_flank_len = flank;
                        lower_bound_flag = CoreBoundFlag::CoreBoundary;
                    }
                    _ if max_flanking_bp < s_start => {
                        flanking_start = s_start - max_flanking_bp;
                        lower_flank_len = max_flanking_bp;
                        lower_bound_flag = CoreBoundFlag::LBoundary;
                    }
                    _ => {
                        flanking_start = 0;
                        lower_flank_len = s_start;
                        lower_bound_flag = CoreBoundFlag::SeqBoundary;
                    }
                }
            }
        } else {
            // Forward strand: left = start side, right = end side.
            if s.left_flag == 1 {
                match prev {
                    Some((p, dist)) if dist <= max_flanking_bp => {
                        let flank = if grants_full(p, p.left_flag) { dist } else { dist / 2 };
                        flanking_start = s_start - flank;
                        lower_flank_len = flank;
                        lower_bound_flag = CoreBoundFlag::CoreBoundary;
                    }
                    _ if max_flanking_bp < s_start => {
                        flanking_start = s_start - max_flanking_bp;
                        lower_flank_len = max_flanking_bp;
                        lower_bound_flag = CoreBoundFlag::LBoundary;
                    }
                    _ => {
                        flanking_start = 0;
                        lower_flank_len = s_start;
                        lower_bound_flag = CoreBoundFlag::SeqBoundary;
                    }
                }
            }
            if s.right_flag == 1 {
                match next {
                    Some((n, dist)) if dist <= max_flanking_bp => {
                        let flank = if grants_full(n, n.right_flag) { dist } else { dist / 2 };
                        flanking_end = s_end + flank;
                        upper_bound_flag = CoreBoundFlag::CoreBoundary;
                    }
                    _ if s_end + max_flanking_bp < seq_size => {
                        flanking_end = s_end + max_flanking_bp;
                        upper_bound_flag = CoreBoundFlag::LBoundary;
                    }
                    _ => {
                        flanking_end = seq_size;
                        upper_bound_flag = CoreBoundFlag::SeqBoundary;
                    }
                }
            }
        }

        let ascii = twobit
            .fetch(&s.name, flanking_start as u64, flanking_end as u64)
            .map_err(crate::Error::Io)?;
        let sub_start = lib.sequence.len() as u64; // start of this subsequence
        lib.sequence.extend(ascii.iter().map(|&c| match c {
            b'A' => 0u8,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => SYM_N,
        }));
        let sub_end_excl = lib.sequence.len() as u64;

        lib.identifiers.push(s.name.clone());
        lib.boundaries.push(sub_end_excl);
        lib.offsets.push(flanking_start as u64);

        let core_len = s.span.len();
        let (orient, left_seq_pos, right_seq_pos) = if s.minus {
            let right = sub_start + lower_flank_len as u64;
            (true, right + core_len - 1, right)
        } else {
            let left = sub_start + lower_flank_len as u64;
            (false, left, left + core_len - 1)
        };

        cores.push(CoreAlignment {
            seq_idx: lib.identifiers.len() - 1,
            left_seq_pos,
            right_seq_pos,
            left_extendable: s.left_flag == 1,
            right_extendable: s.right_flag == 1,
            lower_seq_bound: sub_start,
            upper_seq_bound: sub_end_excl - 1,
            lower_seq_bound_flag: lower_bound_flag,
            upper_seq_bound_flag: upper_bound_flag,
            left_extension_len: 0,
            right_extension_len: 0,
            score: 0,
            orient,
        });
    }

    Ok((lib, cores))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoi_matches_c_semantics() {
        assert_eq!(atoi_like("1"), 1);
        assert_eq!(atoi_like("0"), 0);
        assert_eq!(atoi_like("  42abc"), 42);
        assert_eq!(atoi_like("abc"), 0);
        assert_eq!(atoi_like("-3"), -3);
    }

    #[test]
    fn strtol_matches_c_base_zero_semantics() {
        assert_eq!(strtol_like("1234"), 1234);
        assert_eq!(strtol_like("0x10"), 16);
        assert_eq!(strtol_like("010"), 8);
        assert_eq!(strtol_like(""), 0);
        assert_eq!(strtol_like("12x"), 12);
    }

    /// Write a minimal little-endian, version-0, single-sequence 2bit file.
    fn write_test_2bit(path: &std::path::Path, name: &str, seq: &str) {
        let mut buf: Vec<u8> = Vec::new();
        buf.extend(0x1A412743u32.to_le_bytes()); // signature (LE layout)
        buf.extend(0u32.to_le_bytes()); // version
        buf.extend(1u32.to_le_bytes()); // sequence count
        buf.extend(0u32.to_le_bytes()); // reserved
        buf.push(name.len() as u8);
        buf.extend(name.as_bytes());
        let record_offset = (buf.len() + 4) as u32;
        buf.extend(record_offset.to_le_bytes());
        buf.extend((seq.len() as u32).to_le_bytes()); // dnaSize
        buf.extend(0u32.to_le_bytes()); // nBlockCount
        buf.extend(0u32.to_le_bytes()); // maskBlockCount
        buf.extend(0u32.to_le_bytes()); // reserved
        let code = |c: u8| match c {
            b'T' => 0u8,
            b'C' => 1,
            b'A' => 2,
            b'G' => 3,
            _ => 0,
        };
        for chunk in seq.as_bytes().chunks(4) {
            let mut byte = 0u8;
            for (i, &c) in chunk.iter().enumerate() {
                byte |= code(c) << ((3 - i) * 2);
            }
            buf.push(byte);
        }
        std::fs::write(path, buf).unwrap();
    }

    fn range(name: &str, start: u64, end: u64, l: i32, r: i32, minus: bool) -> RangeRecord {
        RangeRecord {
            name: name.to_string(),
            span: Span::new(start, end).unwrap(),
            left_flag: l,
            right_flag: r,
            minus,
        }
    }

    /// Same-strand neighbors must cede the full inter-core distance in fixed
    /// mode but only the midpoint in C-compat mode (the pointer-compare bug).
    #[test]
    fn same_strand_neighbors_split_full_vs_midpoint() {
        let dir = std::env::temp_dir().join(format!("ramcore-2bit-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.2bit");
        let seq: String = "ACGT".repeat(50); // 200 bp
        write_test_2bit(&path, "seq1", &seq);
        let tb = TwoBitReader::open(&path).unwrap();

        // Two forward cores 20 bp apart, both fully extendable.
        let ranges = vec![
            range("seq1", 10, 20, 1, 1, false),
            range("seq1", 40, 50, 1, 1, false),
        ];

        let (_, cores_fixed) = load_sequence_subset_minimal(&tb, &ranges, 100, false).unwrap();
        let (_, cores_compat) = load_sequence_subset_minimal(&tb, &ranges, 100, true).unwrap();

        // Core 0's right flank runs toward core 1 (gap = 20 bp): the loaded
        // subsequence is [flanking_start, flanking_end) and core 0 starts at
        // input position 0 (10 bp left flank capped by sequence start).
        // Fixed: right flank = 20 → subsequence length 10 + 10 + 20 = 40.
        // Compat: right flank = 10 → 10 + 10 + 10 = 30.
        let len = |c: &CoreAlignment| c.upper_seq_bound - c.lower_seq_bound + 1;
        assert_eq!(len(&cores_fixed[0]), 40, "fixed mode should cede the full gap");
        assert_eq!(len(&cores_compat[0]), 30, "compat mode should midpoint the gap");
        assert_eq!(
            cores_fixed[0].upper_seq_bound_flag,
            CoreBoundFlag::CoreBoundary
        );

        // An opposite-strand neighbor with its facing flag set is midpointed
        // in both modes.
        let ranges_op = vec![
            range("seq1", 10, 20, 1, 1, false),
            range("seq1", 40, 50, 1, 1, true),
        ];
        let (_, fixed_op) = load_sequence_subset_minimal(&tb, &ranges_op, 100, false).unwrap();
        let (_, compat_op) = load_sequence_subset_minimal(&tb, &ranges_op, 100, true).unwrap();
        assert_eq!(len(&fixed_op[0]), 30);
        assert_eq!(len(&compat_op[0]), 30);

        std::fs::remove_dir_all(&dir).ok();
    }
}

// ── Core edge report ──────────────────────────────────────────────────────────

/// One row of the core-edges view: what a copy looks like at the boundary
/// between its aligned core and the flank the extension may consume.
pub struct CoreEdge {
    pub index: usize,
    pub identifier: String,
    /// The core on the source sequence.
    pub span: Span,
    pub minus: bool,
    pub left_extendable: bool,
    pub right_extendable: bool,
    /// Up to 10 bases of flank, in consensus orientation. A trailing `*` marks
    /// a flank cut short by a sequence or neighbour boundary rather than by
    /// the 10-base window.
    pub left_flank: String,
    pub right_flank: String,
    /// 24-column preview of the core: centred when it is short, otherwise the
    /// first ten bases, `....`, and the last ten.
    pub core: String,
}

/// Build the core-edges view for a loaded set of cores.
///
/// A port of `printCoreEdges` in the C tool's `report.c`, which the Rust port
/// had not carried over. The value of the view is that it shows, per copy, the
/// exact bases the extension is about to reason over — a family whose copies
/// all share flanking sequence is a segmental duplication rather than a
/// transposable element, and that is visible here and almost nowhere else.
///
/// Everything is rendered in *consensus* orientation: a minus-strand core
/// reads complemented and walks backwards, so its "left" flank is at a higher
/// coordinate in the source sequence.
pub fn core_edges(lib: &SequenceLibrary, cores: &[CoreAlignment]) -> Vec<CoreEdge> {
    let sym = |z: u8| crate::alphabet::num_to_char(z) as char;
    let comp = |z: u8| crate::alphabet::num_to_char(crate::alphabet::compl(z)) as char;

    cores
        .iter()
        .enumerate()
        .map(|(n, c)| {
            let lower = lib.lower_bound(c.seq_idx);
            let upper = lib.boundaries[c.seq_idx].saturating_sub(1);
            let s = &lib.sequence;

            // Flank window: ten bases, or as many as the boundary allows.
            // Walk outward from the core edge. The in-memory bounds are the
            // only limit needed and they are direction-agnostic; an earlier
            // version also compared against `lower`, which is the wrong end
            // for a minus core and truncated every one of them to nothing.
            let take = |from: i64, step: i64| -> (String, bool) {
                let mut out = String::new();
                let mut got = 0;
                for k in 0..10i64 {
                    let p = from + step * k;
                    if p < lower as i64 || p > upper as i64 {
                        break;
                    }
                    out.push(if c.orient { comp(s[p as usize]) } else { sym(s[p as usize]) });
                    got += 1;
                }
                (out, got < 10)
            };

            let core_len = if c.orient {
                c.left_seq_pos.saturating_sub(c.right_seq_pos) + 1
            } else {
                c.right_seq_pos.saturating_sub(c.left_seq_pos) + 1
            };

            // Walk the core in consensus order: outward from left_seq_pos.
            let step: i64 = if c.orient { -1 } else { 1 };
            let at = |i: i64| -> char {
                let p = c.left_seq_pos as i64 + step * i;
                if p < 0 || p as usize >= s.len() {
                    return ' ';
                }
                if c.orient { comp(s[p as usize]) } else { sym(s[p as usize]) }
            };
            let core = if core_len < 20 {
                let pad = (24 - core_len as usize) / 2;
                let body: String = (0..core_len as i64).map(at).collect();
                format!("{}{}", " ".repeat(pad), body)
            } else {
                let head: String = (0..10).map(at).collect();
                let tail: String = ((core_len as i64 - 10)..core_len as i64).map(at).collect();
                format!("{head}....{tail}")
            };

            // Left flank lies before the core start in consensus order.
            let (mut lf, lcut) = take(c.left_seq_pos as i64 - step, -step);
            let lf: String = { lf = lf.chars().rev().collect(); lf };
            let right_start = c.left_seq_pos as i64 + step * core_len as i64;
            let (rf, rcut) = take(right_start, step);

            CoreEdge {
                index: n,
                identifier: lib.identifiers[c.seq_idx].clone(),
                span: lib.to_source(c.seq_idx, c.core_span()),
                minus: c.orient,
                left_extendable: c.left_extendable,
                right_extendable: c.right_extendable,
                left_flank: if lcut { format!("*{lf}") } else { lf },
                right_flank: if rcut { format!("{rf}*") } else { rf },
                core,
            }
        })
        .collect()
}
