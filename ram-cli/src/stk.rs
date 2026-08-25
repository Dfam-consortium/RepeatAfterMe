//! Native Stockholm mode — the parts of `extend-stk.pl` that belong with
//! the binary (items 1-3 of its workflow):
//!
//! 1. **Stockholm plumbing**: parse a (multi-family) RepeatModeler seed
//!    alignment via `dfam-stk-io`, taking member genomic coordinates from
//!    the Smitten identifiers.
//! 2. **Extendable-edge detection**: a member whose alignment starts
//!    within 10 columns of the MSA edge is extendable on that side
//!    (the Perl's leading/trailing `.`-padding rule); members extendable
//!    on neither side are dropped, and a family needs more than
//!    `min_aligning_seqs` extendable members to be attempted.
//! 3. **Adaptive parameterization**: average Kimura divergence of the
//!    members against the reference (the `#=GC RF` line, or a recalled
//!    consensus when RF is absent or contains `x`), overridable by an
//!    `mDiv=` note in the `#=GF DE` description, buckets into a matrix —
//!    `<16 -> 14p43g, <19 -> 18, <22.5 -> 20, else 25` — and
//!    `minimprovement = min_aligning_seqs x diagonal_average`
//!    (10, or 9 for 25p43g), exactly as `extend-stk.pl` computes them.
//!
//! Steps 5-6 of the script (Refiner-based MSA refinement and Stockholm
//! reassembly) are RepeatModeler-workflow specifics and stay in the
//! script.

use std::path::Path;

use aln_core::consensus::{build_consensus_from_sequences, ConsensusParams};
use aln_core::stats::{kimura_divergence, Masking};
use ram_core::library::RangeRecord;

/// Everything needed to run one family from a Stockholm record.
pub struct FamilyPlan {
    pub label: String,
    /// 1-based record number in the file (stable output suffix `rec<N>`).
    pub record_num: usize,
    /// The ungapped reference used for divergence — also the piece the
    /// caller sandwiches between the left/right extension consensi to
    /// form the combined extended consensus.
    pub reference_ungapped: String,
    pub ranges: Vec<RangeRecord>,
    pub matrix_name: &'static str,
    pub min_improvement: i32,
    pub divergence: f64,
    pub divergence_overridden: bool,
    pub n_members: usize,
    pub n_extendable: usize,
    /// Set when the family is skipped (too few extendable members).
    pub skip_reason: Option<String>,
}

/// The Perl's divergence buckets: matrix name and its diagonal average.
fn matrix_for_divergence(div: f64) -> (&'static str, i32) {
    if div >= 22.5 {
        ("25p43g", 9)
    } else if div >= 19.0 {
        ("20p43g", 10)
    } else if div >= 16.0 {
        ("18p43g", 10)
    } else {
        ("14p43g", 10)
    }
}

/// Extendable-edge rule on the aligned row: at most 10 columns of `.`
/// padding between the MSA edge and the first/last aligned base.
fn edge_flags(aligned: &str) -> (bool, bool) {
    let bytes = aligned.as_bytes();
    let lead = bytes.iter().take_while(|&&b| b == b'.').count();
    let trail = bytes.iter().rev().take_while(|&&b| b == b'.').count();
    let has_base = bytes.iter().any(|&b| b != b'.');
    (has_base && lead <= 10, has_base && trail <= 10)
}

/// Parse every record in a Stockholm file into a per-family plan with
/// the adaptive matrix/minimprovement choices. The caller substitutes
/// explicit `-matrix`/`-minimprovement` flags when the user gave them.
pub fn plan_families(
    path: &Path,
    min_aligning_seqs: i32,
) -> ram_core::Result<Vec<FamilyPlan>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut plans = Vec::new();

    for record in dfam_stk_io::iter_records(reader) {
        let record = record.map_err(|e| {
            ram_core::Error::Format(format!("Stockholm parse error in {}: {e}", path.display()))
        })?;
        let label = record.label();

        // Member rows with genomic coordinates; skip consensus-like rows.
        let mut ranges: Vec<RangeRecord> = Vec::new();
        let mut aligned_rows: Vec<&str> = Vec::new();
        let mut n_extendable = 0usize;
        for row in &record.sequences {
            let (Some(seq_id), Some(s0), Some(e0), Some(orient)) =
                (&row.sequence_id, row.seq_start, row.seq_end, row.orient)
            else {
                continue;
            };
            aligned_rows.push(&row.aligned_seq);
            let (left, right) = edge_flags(&row.aligned_seq);
            if !(left || right) {
                continue; // extendable on neither side: dropped, as in the Perl
            }
            n_extendable += 1;
            ranges.push(RangeRecord {
                // RepeatModeler seed names are 1-based inclusive and
                // dfam-stk-io returns them verbatim; RangeRecord wants
                // 0-based half-open (enforced by the calibration test).
                name: seq_id.clone(),
                start: s0 as i64 - 1,
                end: e0 as i64,
                left_flag: i32::from(left),
                right_flag: i32::from(right),
                minus: orient == '-',
            });
        }
        let n_members = record.sequences.len();

        // Seed alignments pad with '.'; aln-core's conventions are ' '
        // for not-reached edge columns and '-' for interior gaps. Without
        // this conversion every column looks covered and the consensus
        // keeps all 330 columns of a 163-match-column family (MultAln
        // yields 163 — verified against the original script on
        // ce10-fam1).
        let normalized: Vec<Vec<u8>> = aligned_rows
            .iter()
            .map(|s| {
                let b = s.as_bytes();
                let lead = b.iter().take_while(|&&c| c == b'.').count();
                let trail = b.iter().rev().take_while(|&&c| c == b'.').count();
                b.iter()
                    .enumerate()
                    .map(|(i, &c)| {
                        if c == b'.' {
                            if i < lead || i >= b.len() - trail {
                                b' '
                            } else {
                                b'-'
                            }
                        } else {
                            c
                        }
                    })
                    .collect()
            })
            .collect();

        // Reference: #=GC RF unless absent or containing 'x'.
        let rf = record.gc.get("RF").cloned();
        let reference: String = match rf {
            Some(r) if !r.to_ascii_lowercase().contains('x') => r,
            _ => {
                let rows: Vec<&[u8]> = normalized.iter().map(|v| v.as_slice()).collect();
                String::from_utf8_lossy(&build_consensus_from_sequences(
                    &rows,
                    &ConsensusParams::default(),
                ))
                .into_owned()
            }
        };

        // Average unadjusted Kimura divergence of members vs the reference.
        let mut divs: Vec<f64> = Vec::new();
        for row in &normalized {
            if row.len() != reference.len() {
                continue;
            }
            if let Ok(d) =
                kimura_divergence(row, reference.as_bytes(), false, Masking::Ignore)
            {
                if let Some(v) = d.value {
                    divs.push(v);
                }
            }
        }
        let mut divergence = if divs.is_empty() {
            0.0
        } else {
            divs.iter().sum::<f64>() / divs.len() as f64
        };

        // mDiv= override in the description, as in the Perl.
        let mut overridden = false;
        if let Some(de) = record.gf_first("DE") {
            if let Some(idx) = de.find("mDiv=") {
                let tail = &de[idx + 5..];
                let num: String = tail
                    .chars()
                    .take_while(|c| c.is_ascii_digit() || *c == '.')
                    .collect();
                if let Ok(v) = num.parse::<f64>() {
                    divergence = v;
                    overridden = true;
                }
            }
        }

        let (adaptive_matrix, diag_avg) = matrix_for_divergence(divergence);
        let min_improvement = min_aligning_seqs * diag_avg;

        let skip_reason = if n_extendable <= min_aligning_seqs as usize {
            Some(format!(
                "too few extendable sequences ({n_extendable}) for extension"
            ))
        } else {
            None
        };

        let reference_ungapped: String = reference
            .chars()
            .filter(|c| *c != '-' && *c != '.')
            .collect();
        plans.push(FamilyPlan {
            label,
            record_num: record.record_num,
            reference_ungapped,
            ranges,
            matrix_name: adaptive_matrix,
            min_improvement,
            divergence,
            divergence_overridden: overridden,
            n_members,
            n_extendable,
            skip_reason,
        });
    }
    Ok(plans)
}

/// Append a `>combined` record — left extension + ungapped reference +
/// right extension — to a `-cons` file written by the extension run
/// (which contains `>left-extension` / `>right-extension` sections).
pub fn append_combined_consensus(cons_path: &str, reference: &str) -> std::io::Result<()> {
    let text = std::fs::read_to_string(cons_path)?;
    let mut left = String::new();
    let mut right = String::new();
    let mut cur: Option<&mut String> = None;
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            cur = match name.trim() {
                n if n.starts_with("left-extension") => Some(&mut left),
                n if n.starts_with("right-extension") => Some(&mut right),
                _ => None,
            };
        } else if let Some(buf) = cur.as_deref_mut() {
            buf.push_str(line.trim());
        }
    }
    let mut out = std::fs::OpenOptions::new().append(true).open(cons_path)?;
    use std::io::Write;
    writeln!(out, ">combined")?;
    let combined = format!("{left}{reference}{right}");
    for chunk in combined.as_bytes().chunks(80) {
        out.write_all(chunk)?;
        writeln!(out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_stk(tag: &str, content: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!(
            "ram_stk_test_{}_{tag}.stk",
            std::process::id()
        ));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        p
    }

    /// Coordinate calibration: a Smitten name `chr1:101-200_+` (1-based
    /// inclusive in RepeatModeler files) must produce the 0-based
    /// half-open RangeRecord {start:100, end:200} that the loader and
    /// `stk2ranges.py` agree on.
    #[test]
    fn coordinates_match_stk2ranges_convention() {
        let stk = "# STOCKHOLM 1.0\n\
                   #=GF ID testfam\n\
                   chr1:101-200_+  ACGTACGTAC\n\
                   chr1:301-400_-  ACGTACGTAC\n\
                   chr1:501-600_+  ACGTACGTAC\n\
                   chr1:701-800_+  ACGTACGTAC\n\
                   chr1:901-999_+  ACGTACGTAC\n\
                   //\n";
        let p = write_stk("coords", stk);
        let plans = plan_families(&p, 3).unwrap();
        std::fs::remove_file(&p).unwrap();
        assert_eq!(plans.len(), 1);
        let r = &plans[0].ranges[0];
        assert_eq!((r.start, r.end), (100, 200), "calibration: got {r:?}");
        assert!(plans[0].ranges[1].minus);
        assert!(plans[0].skip_reason.is_none());
    }

    #[test]
    fn edge_rule_and_skip_logic() {
        // Second row has 11 leading dots and 11 trailing: not extendable
        // either side -> dropped; family then has 3 extendable members,
        // which is NOT more than min_aligning_seqs=3 -> skipped.
        let stk = "# STOCKHOLM 1.0\n\
                   #=GF ID testfam\n\
                   chr1:101-130_+  ACGTACGTACGTACGTACGTACGTACGTAC\n\
                   chr1:301-308_+  ...........ACGTACGT...........\n\
                   chr1:501-530_+  ACGTACGTACGTACGTACGTACGTACGTAC\n\
                   chr1:701-730_+  ACGTACGTACGTACGTACGTACGTACGTAC\n\
                   //\n";
        let p = write_stk("edge", stk);
        let plans = plan_families(&p, 3).unwrap();
        std::fs::remove_file(&p).unwrap();
        assert_eq!(plans[0].n_extendable, 3);
        assert!(plans[0].skip_reason.is_some());
    }

    #[test]
    fn divergence_buckets_and_mdiv_override() {
        assert_eq!(matrix_for_divergence(10.0), ("14p43g", 10));
        assert_eq!(matrix_for_divergence(16.0), ("18p43g", 10));
        assert_eq!(matrix_for_divergence(19.0), ("20p43g", 10));
        assert_eq!(matrix_for_divergence(22.5), ("25p43g", 9));
        // Identical rows -> divergence 0 -> 14p43g; mDiv override flips it.
        let stk = "# STOCKHOLM 1.0\n\
                   #=GF ID testfam\n\
                   #=GF DE Source:gsa, mDiv=23.22, something\n\
                   chr1:101-140_+  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                   chr1:301-340_+  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                   chr1:501-540_+  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                   chr1:701-740_+  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                   //\n";
        let p = write_stk("mdiv", stk);
        let plans = plan_families(&p, 3).unwrap();
        std::fs::remove_file(&p).unwrap();
        let f = &plans[0];
        assert!(f.divergence_overridden);
        assert_eq!(f.matrix_name, "25p43g");
        assert_eq!(f.min_improvement, 3 * 9);
    }
}
