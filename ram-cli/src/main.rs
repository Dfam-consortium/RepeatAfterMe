//! ram-extend — flag-compatible Rust port of the RAMExtend CLI.
//!
//! Accepts the same single-dash flags as the C tool (`-twobit`, `-ranges`,
//! `-L`, `-bandwidth`, `-matrix`, ...) and produces byte-identical `-outtsv`,
//! `-cons` and `-outfa` files, plus the stdout lines downstream tools scrape
//! (`Extended right: N bp` / `Extended left : N bp`).
//!
//! Three C bugs are fixed by default; `-ccompat` restores exact C behavior
//! for differential testing:
//!  - the loader midpoints same-strand neighbor gaps (a `char*` pointer
//!    comparison that can never be true) — fixed to cede the full distance;
//!  - the right-extension consensus block wraps at a phase offset by the
//!    left extension length — fixed to clean 80-column wrapping;
//!  - the `-outtsv` anchor_range omits the subsequence offset and the
//!    orientation swap that the stdout report applies — fixed to match.
//!
//! PORT NOTES (deviations from C, all deliberate):
//!  - `-outmat` (DP path debugging) is not implemented; passing it is an
//!    error rather than silently diverging.
//!  - More ranges than `-maxoccurrences` is a clean error; the C tool indexes
//!    out of bounds.
//!  - The cosmetic core-edge pretty-print at startup is omitted.
//!  - Like the C tool, `-gapopen` only takes effect together with a gap
//!    extension flag; the usage text's `-gapextn` spelling (which the C tool
//!    documents but ignores) is accepted as an alias unless `-ccompat`.

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use ram_core::alphabet::{compl, num_to_char, SYM_N};
use ram_core::engine::{
    apply_overlap_avoidance, extend_alignment, Direction, ExtendParams, ExtensionResult,
    ScoreArena,
};
use ram_core::library::{load_sequence_subset_minimal, read_ranges, CoreBoundFlag, RangeRecord};

mod stk;
use ram_core::matrix::ScoringSystem;
use ram_core::twobit::TwoBitReader;

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "-rust");

/// C `cmd_line_opts` semantics: scan argv for an exact flag match; the value
/// is the following argument. Unknown arguments are ignored.
struct Args(Vec<String>);

impl Args {
    fn get_string(&self, flag: &str) -> Option<&str> {
        self.0
            .iter()
            .position(|a| a == flag)
            .and_then(|i| self.0.get(i + 1))
            .map(|s| s.as_str())
    }
    fn get_int(&self, flag: &str) -> Option<i64> {
        self.get_string(flag).map(|v| v.parse().unwrap_or(0))
    }
    fn get_bool(&self, flag: &str) -> bool {
        self.0.iter().any(|a| a == flag)
    }
}

fn usage() -> ExitCode {
    println!(
        "RAMExtend Version {VERSION}\n\n\
         Usage:\n  ram-extend -twobit <seq.2bit> -ranges <ranges.tsv> [opts]\n\
  ram-extend -twobit <seq.2bit> -stk <seed.stk> [opts]\n\n\
         -stk mode (native extend-stk.pl items 1-3): parses RepeatModeler\n\
         seed alignments (multi-family supported; outputs suffixed by\n\
         record label), derives extendable flags from the alignment edges,\n\
         and picks matrix/minimprovement adaptively from the MSA's Kimura\n\
         divergence (-min_aligning_seqs, default 3; explicit -matrix or\n\
         -minimprovement override; defaults -L 20000 -bandwidth 40).\n\n\
         See the C RAMExtend usage for full option documentation; this port\n\
         accepts: -cons -outtsv -outfa -addflanking -L -minimprovement\n\
         -maxoccurrences -cappenalty -stopafter -matrix -gapopen -gapext\n\
         -match -mismatch -gap -bandwidth -version -v[v[v[v]]]\n\n\
         Port-specific:\n\
         -ccompat   # reproduce three C bugs byte-for-byte (neighbor-gap\n\
                    # midpointing, right-extension wrap phase, TSV anchor_range)\n"
    );
    ExitCode::from(1)
}

fn main() -> ExitCode {
    let args = Args(std::env::args().skip(1).collect());

    if args.get_bool("-version") {
        println!("RAMExtend Version {VERSION}");
        return ExitCode::SUCCESS;
    }
    if args.get_bool("-outmat") {
        eprintln!("ram-extend: -outmat is not supported by this port");
        return ExitCode::from(1);
    }

    let ranges_file = args.get_string("-ranges");
    let stk_file = args.get_string("-stk");
    if ranges_file.is_none() && stk_file.is_none() {
        return usage();
    }
    if ranges_file.is_some() && stk_file.is_some() {
        eprintln!("ram-extend: -ranges and -stk are mutually exclusive");
        return ExitCode::from(1);
    }
    let Some(twobit_file) = args.get_string("-twobit") else {
        if args.get_string("-sequence").is_some() {
            println!("-sequence is deprecated!....may return someday");
            return ExitCode::from(1);
        }
        return usage();
    };

    let flanking = args.get_int("-addflanking").unwrap_or(0);
    let outtsv = args.get_string("-outtsv");
    let outfa = args.get_string("-outfa");
    let cons_file = args.get_string("-cons");
    let l_max = args.get_int("-L").unwrap_or(10000) as i32;
    let bandwidth = args.get_int("-bandwidth").unwrap_or(14) as i32;
    let maxn = args.get_int("-maxoccurrences").unwrap_or(10000) as usize;
    let when_to_stop = args.get_int("-stopafter").unwrap_or(100) as i32;

    // Unlike the C tool (which parses -threads but never uses it), this port
    // parallelizes the per-row DP across cores. 0 = all logical CPUs.
    let threads = args.get_int("-threads").unwrap_or(0).max(0) as usize;
    ram_core::engine::set_threads(threads);

    let c_compat = args.get_bool("-ccompat");
    let matrix_name = args.get_string("-matrix").unwrap_or("20p43g");
    let scoring = if matrix_name == "repeatscout" {
        match (
            args.get_int("-match"),
            args.get_int("-mismatch"),
            args.get_int("-gap"),
        ) {
            (Some(m), Some(mm), Some(g)) => {
                ScoringSystem::repeat_scout(m as i32, mm as i32, g as i32)
            }
            _ => ScoringSystem::repeat_scout(1, -1, -5),
        }
    } else {
        let base = match ScoringSystem::by_name(matrix_name) {
            Ok(s) => s,
            Err(e) => {
                println!("{e}");
                return ExitCode::from(1);
            }
        };
        let gapext = args.get_int("-gapext").or(if c_compat {
            None
        } else {
            args.get_int("-gapextn")
        });
        match (args.get_int("-gapopen"), gapext) {
            (Some(go), Some(ge)) => base.with_gap_penalties(go as i32, ge as i32),
            _ => base,
        }
    };

    // Matrix-dependent defaults (C main).
    let (def_minimp, def_cap) = match matrix_name {
        "25p43g" => (24, -90),
        "repeatscout" => (3, -20),
        _ => (27, -90),
    };
    let min_improvement = args.get_int("-minimprovement").unwrap_or(def_minimp as i64) as i32;
    let cap_penalty = args.get_int("-cappenalty").unwrap_or(def_cap as i64) as i32;

    let verbose = if args.get_bool("-vvvvvv") {
        20
    } else if args.get_bool("-vvvvv") {
        12
    } else if args.get_bool("-vvvv") {
        10
    } else if args.get_bool("-vvv") {
        3
    } else if args.get_bool("-vv") {
        2
    } else {
        i32::from(args.get_bool("-v"))
    };

    if let Some(stk_path) = stk_file {
        // Stockholm mode: extend-stk.pl items 1-3 natively. Script-mode
        // defaults for L and bandwidth unless given explicitly; adaptive
        // per-family matrix/minimprovement unless given explicitly.
        let l_max = args.get_int("-L").unwrap_or(20000) as i32;
        let bandwidth = args.get_int("-bandwidth").unwrap_or(40) as i32;
        let min_aligning_seqs = args.get_int("-min_aligning_seqs").unwrap_or(3) as i32;
        let explicit_matrix = args.get_string("-matrix");
        let explicit_minimp = args.get_int("-minimprovement").map(|v| v as i32);

        let plans = match stk::plan_families(Path::new(stk_path), min_aligning_seqs) {
            Ok(p) => p,
            Err(e) => {
                println!("{e}");
                return ExitCode::from(1);
            }
        };
        let multi = plans.len() > 1;
        // Deterministic per-record suffix so wrapper scripts can pair
        // outputs with records without re-deriving label sanitization.
        let suffixed = |base: Option<&str>, rec: usize| -> Option<String> {
            base.map(|b| {
                if multi {
                    format!("{b}.rec{rec}")
                } else {
                    b.to_string()
                }
            })
        };
        let mut failures = 0usize;
        for plan in &plans {
            println!("\n== Family {} ==", plan.label);
            println!(
                "  members: {} ({} extendable); Kimura divergence: {:.2}{}",
                plan.n_members,
                plan.n_extendable,
                plan.divergence,
                if plan.divergence_overridden { " (mDiv override)" } else { "" }
            );
            if let Some(reason) = &plan.skip_reason {
                println!("  ** skipped: {reason} **");
                continue;
            }
            let fam_matrix = explicit_matrix.unwrap_or(plan.matrix_name);
            let fam_minimp = explicit_minimp.unwrap_or(plan.min_improvement);
            let fam_scoring = match ScoringSystem::by_name(fam_matrix) {
                Ok(s) => match (args.get_int("-gapopen"), args.get_int("-gapext")) {
                    (Some(go), Some(ge)) => s.with_gap_penalties(go as i32, ge as i32),
                    _ => s,
                },
                Err(_) => {
                    println!("Matrix name not found!");
                    return ExitCode::from(1);
                }
            };
            println!("  matrix: {fam_matrix}; minimprovement: {fam_minimp}");
            let outtsv_s = suffixed(outtsv, plan.record_num);
            let outfa_s = suffixed(outfa, plan.record_num);
            let cons_s = suffixed(cons_file, plan.record_num);
            let result = run(
                RunConfig {
                    twobit_file,
                    ranges_desc: format!("{stk_path} [{}]", plan.label),
                    flanking,
                    outtsv: outtsv_s.as_deref(),
                    outfa: outfa_s.as_deref(),
                    cons_file: cons_s.as_deref(),
                    l_max,
                    bandwidth,
                    maxn,
                    when_to_stop,
                    matrix_name: fam_matrix,
                    scoring: fam_scoring,
                    min_improvement: fam_minimp,
                    cap_penalty,
                    verbose,
                    c_compat,
                },
                plan.ranges.clone(),
            );
            match result {
                Err(e) => {
                    println!("{e}");
                    failures += 1;
                }
                Ok(()) => {
                    // Append the combined extended consensus (left
                    // extension + ungapped reference + right extension)
                    // so wrappers need no sequence reassembly of their
                    // own.
                    if let Some(cons_path) = &cons_s {
                        if let Err(e) =
                            stk::append_combined_consensus(cons_path, &plan.reference_ungapped)
                        {
                            println!("WARNING: could not append combined consensus: {e}");
                        }
                    }
                }
            }
        }
        return if failures == 0 {
            ExitCode::SUCCESS
        } else {
            ExitCode::from(1)
        };
    }

    let ranges_file = ranges_file.unwrap();
    let ranges = match read_ranges(Path::new(ranges_file)) {
        Ok(r) => r,
        Err(e) => {
            println!("{e}");
            return ExitCode::from(1);
        }
    };
    match run(
        RunConfig {
            twobit_file,
            ranges_desc: ranges_file.to_string(),
            flanking,
            outtsv,
            outfa,
            cons_file,
            l_max,
            bandwidth,
            maxn,
            when_to_stop,
            matrix_name,
            scoring,
            min_improvement,
            cap_penalty,
            verbose,
            c_compat,
        },
        ranges,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            println!("{e}");
            ExitCode::from(1)
        }
    }
}

struct RunConfig<'a> {
    twobit_file: &'a str,
    ranges_desc: String,
    flanking: i64,
    outtsv: Option<&'a str>,
    outfa: Option<&'a str>,
    cons_file: Option<&'a str>,
    l_max: i32,
    bandwidth: i32,
    maxn: usize,
    when_to_stop: i32,
    matrix_name: &'a str,
    scoring: ScoringSystem,
    min_improvement: i32,
    cap_penalty: i32,
    verbose: i32,
    c_compat: bool,
}

fn run(cfg: RunConfig, ranges: Vec<RangeRecord>) -> ram_core::Result<()> {
    let twobit = TwoBitReader::open(Path::new(cfg.twobit_file))?;
    let (lib, mut cores) = load_sequence_subset_minimal(
        &twobit,
        &ranges,
        (cfg.l_max + cfg.bandwidth) as i64,
        cfg.c_compat,
    )?;
    let n = cores.len();
    if n > cfg.maxn {
        return Err(ram_core::Error::Limit(format!(
            "Number of ranges ({n}) exceeds -maxoccurrences ({}); the C tool \
             would corrupt memory here — raise -maxoccurrences instead.",
            cfg.maxn
        )));
    }

    println!("\nRAMExtend Version {VERSION}");
    print_parameters(&cfg);
    println!("Read in {n} ranges, and {} bp of sequence\n", lib.length());

    let l = cfg.l_max;
    let l_spacer = 1usize;
    let mut master = vec![0u8; 2 * l as usize + l_spacer + 1];
    for x in 0..l_spacer {
        master[l as usize + x] = SYM_N;
    }

    let params = ExtendParams {
        bandwidth: cfg.bandwidth,
        cap_penalty: cfg.cap_penalty,
        min_improvement: cfg.min_improvement,
        l_max: l,
        when_to_stop: cfg.when_to_stop,
    };
    let mut arena = ScoreArena::new(n, cfg.bandwidth);

    // Extend RIGHT first.
    let right: ExtensionResult = extend_alignment(
        Direction::Right,
        &mut cores,
        &mut arena,
        &lib,
        &mut master,
        &params,
        &cfg.scoring,
    );
    // PORT NOTE: deliberate deviation. The C gates this warning on
    // `row_idx == L - 1` tested *after* the loop, so an extension that runs
    // the full L rows to natural completion — the ordinary runaway — prints
    // nothing, and only a break on exactly the final row warns. That is close
    // to the complement of what the message claims. We warn on `capped`, which
    // is true when the loop exhausted L rather than converging.
    //
    // The `Extended left :` / `Extended right:` lines below are unchanged:
    // Refiner greps them to recover the lengths, precisely because it could
    // not rely on this warning.
    if right.capped {
        println!("WARNING: Extended sequence right to the limit ( L={l} ).");
    }
    println!("Extended right: {} bp", right.bp);

    for ev in apply_overlap_avoidance(&mut cores, &lib) {
        println!(
            "OVERLAP AVOIDANCE: seqid {} extended to {}, limits seqid {} with existing {}_bound = {} because it's pos_in_r={}",
            ev.s_seq_idx,
            ev.extended_to,
            ev.r_seq_idx,
            if ev.upper { "upper" } else { "lower" },
            ev.old_bound,
            ev.pos_in_r
        );
    }

    let masterend = (l + 1 + right.bp) as usize;

    // Extend LEFT second.
    let left = extend_alignment(
        Direction::Left,
        &mut cores,
        &mut arena,
        &lib,
        &mut master,
        &params,
        &cfg.scoring,
    );
    if left.capped {
        println!("WARNING: Extended sequence left to the limit ( L={l} ).");
    }
    println!("Extended left : {} bp", left.bp);

    // Both directions capped means the model never found an edge in either:
    // the copies go on agreeing for as far as they are allowed to look. That
    // is a satellite array or a segmental duplication, not an element with
    // boundaries, and the extension is arbitrary. Reported here so a caller
    // does not have to infer it from the two lengths.
    if left.capped && right.capped {
        println!(
            "WARNING: Extension reached the limit in BOTH directions ( L={l} ) - \
probably a satellite or segmental duplication; treat the extension as unreliable."
        );
    }

    let masterstart = (l - left.bp) as usize;

    if right.bp > 0 || left.bp > 0 {
        // PORT NOTE: the C right-extension block wraps at a phase offset by
        // the left extension length; fixed mode wraps each block cleanly.
        let right_start = l as usize + l_spacer;
        let right_phase = if cfg.c_compat { masterstart } else { right_start };

        let mut stdout = std::io::stdout().lock();
        if left.bp > 0 {
            writeln!(stdout, ">left-extension").ok();
            write_wrapped(&mut stdout, &master, masterstart, l as usize, masterstart);
        }
        if right.bp > 0 {
            writeln!(stdout, ">right-extension").ok();
            write_wrapped(&mut stdout, &master, right_start, masterend, right_phase);
        }
        drop(stdout);

        if let Some(path) = cfg.cons_file {
            let mut fp = std::fs::File::create(path).map_err(|_| {
                ram_core::Error::Format(format!("Could not open input file {path}"))
            })?;
            if left.bp > 0 {
                writeln!(fp, ">left-extension {} bp", left.bp)?;
                write_wrapped(&mut fp, &master, masterstart, l as usize, masterstart);
            }
            if right.bp > 0 {
                writeln!(fp, ">right-extension {} bp", right.bp)?;
                write_wrapped(&mut fp, &master, right_start, masterend, right_phase);
            }
        }

        write_report(&cfg, &lib, &cores)?;
    }

    let _ = cfg.verbose;
    Ok(())
}

/// The C tool's 80-column wrapping: newline whenever
/// `(x - masterstart) % 80 == 79`, plus a trailing newline unless the block
/// ended exactly on a wrap. For the right-extension block the phase is offset
/// by the left extension length — replicated faithfully.
fn write_wrapped<W: Write>(w: &mut W, master: &[u8], start: usize, end: usize, phase_origin: usize) {
    for x in start..end {
        w.write_all(&[num_to_char(master[x])]).ok();
        if (x - phase_origin) % 80 == 79 {
            w.write_all(b"\n").ok();
        }
    }
    if (end - phase_origin) % 80 > 0 {
        w.write_all(b"\n").ok();
    }
}

fn print_parameters(cfg: &RunConfig) {
    println!("--------------------------------------------------------------");
    println!("Parameters:");
    println!("  VERBOSE {}", cfg.verbose);
    println!("  SEQUENCE_FILE {}", cfg.twobit_file);
    println!("  RANGES_FILE {}", cfg.ranges_desc);
    println!("  L {}", cfg.l_max);
    println!("  BANDWIDTH (bandwidth) {}", cfg.bandwidth);
    println!("  MAXN {}", cfg.maxn);
    if cfg.matrix_name == "repeatscout" {
        println!("  SCORING SYSTEM: Original RepeatScout method");
        println!("     - GAP = {}", cfg.scoring.gapextn);
        println!("     - CAPPENALTY {}", cfg.cap_penalty);
        println!("     - MINIMPROVEMENT {}", cfg.min_improvement);
    } else {
        println!("  SCORING SYSTEM: Internally coded matrix '{}'", cfg.matrix_name);
        println!("     - GAP_OPEN = {}", cfg.scoring.gapopen);
        println!("     - GAP_EXT = {}", cfg.scoring.gapextn);
        println!("     - CAPPENALTY {}", cfg.cap_penalty);
        println!("     - MINIMPROVEMENT {}", cfg.min_improvement);
    }
    println!("  WHEN_TO_STOP {}", cfg.when_to_stop);
    println!("--------------------------------------------------------------");
}

/// The final report: stdout table plus optional -outtsv / -outfa files,
/// byte-compatible with the C `main` (including the TSV's un-offset
/// anchor_range and the `-addflanking` start clamp to 1).
fn write_report(
    cfg: &RunConfig,
    lib: &ram_core::library::SequenceLibrary,
    cores: &[ram_core::library::CoreAlignment],
) -> ram_core::Result<()> {
    let mut fp_tsv = match cfg.outtsv {
        Some(p) => Some(std::fs::File::create(p).map_err(|_| {
            ram_core::Error::Format(format!("Could not create the TSV output file {p}"))
        })?),
        None => None,
    };
    let mut fp_fa = match cfg.outfa {
        Some(p) => Some(std::fs::File::create(p).map_err(|_| {
            ram_core::Error::Format(format!("Could not create the FASTA output file {p}"))
        })?),
        None => None,
    };

    println!("\n\nExtended Sequences Report:");
    println!("  *** Extended sequences are in 1-based, fully closed coordinates ***");
    println!("SEQID  SEQSTART SEQEND ORIENT  EXTENSION_DETAILS");

    for (x, s) in cores.iter().enumerate() {
        let seq_idx = s.seq_idx;
        let ident = &lib.identifiers[seq_idx];
        let seq_lower = lib.lower_bound(seq_idx);
        let seq_upper = lib.boundaries[seq_idx];

        // The report is 1-based fully closed on the source sequence.
        let orient_ch = if s.orient { '-' } else { '+' };
        let extended = s.extended_span();
        let (extended_start, extended_end) = lib
            .to_source(seq_idx, extended)
            .as_1b_closed()
            .expect("a core covers at least one base");
        let (core_start, core_end) = lib
            .to_source(seq_idx, s.core_span())
            .as_1b_closed()
            .expect("a core covers at least one base");
        let extended_length = extended.len() as i64;

        let mut line = format!(
            "{ident}\t{extended_start}\t{extended_end}\t{orient_ch}\tn={x},anchor_range={core_start}-{core_end}"
        );
        line.push_str(&if s.left_extendable {
            format!(",extended_left={}", s.left_extension_len)
        } else {
            ",extended_left=*".to_string()
        });
        line.push_str(&if s.right_extendable {
            format!(",extended_right={}", s.right_extension_len)
        } else {
            ",extended_right=*".to_string()
        });
        line.push_str(&format!(",len={extended_length},score={}", s.score));

        // Boundary-limit annotations (stdout only in C).
        if orient_ch == '+' {
            if s.left_extendable
                && s.left_seq_pos
                    .wrapping_sub(s.left_extension_len as u64)
                    .wrapping_sub(s.lower_seq_bound)
                    < 20
            {
                match s.lower_seq_bound_flag {
                    CoreBoundFlag::SeqBoundary => line.push_str(",leftSeqLimit"),
                    CoreBoundFlag::LBoundary => line.push_str(",leftExtLimit"),
                    CoreBoundFlag::CoreBoundary => line.push_str(",leftCoreLimit"),
                    CoreBoundFlag::ExtBoundary => {}
                }
            }
            if s.right_extendable
                && s.upper_seq_bound
                    .wrapping_sub(s.right_seq_pos.wrapping_add(s.right_extension_len as u64))
                    < 20
            {
                match s.upper_seq_bound_flag {
                    CoreBoundFlag::SeqBoundary => line.push_str(",rightSeqLimit"),
                    CoreBoundFlag::LBoundary => line.push_str(",rightExtLimit"),
                    CoreBoundFlag::CoreBoundary => line.push_str(",rightCoreLimit"),
                    CoreBoundFlag::ExtBoundary => {}
                }
            }
        } else {
            if s.left_extendable
                && s.upper_seq_bound
                    .wrapping_sub(s.left_seq_pos.wrapping_add(s.left_extension_len as u64))
                    < 20
            {
                line.push_str(",leftCoreLimit");
            }
            if s.right_extendable
                && s.right_seq_pos
                    .wrapping_sub(s.right_extension_len as u64)
                    .wrapping_sub(s.lower_seq_bound)
                    < 20
            {
                line.push_str(",rightCoreLimit");
            }
        }
        println!("{line}");

        if let Some(fp) = fp_tsv.as_mut() {
            // PORT NOTE: the C TSV writes the anchor range without the
            // subsequence offset and without orient swapping (so it is
            // descending for a minus core), disagreeing with its own stdout
            // report; fixed mode matches stdout.
            let (anchor_a, anchor_b) = if cfg.c_compat {
                (
                    s.left_seq_pos - seq_lower + 1,
                    s.right_seq_pos - seq_lower + 1,
                )
            } else {
                (core_start, core_end)
            };
            let mut t = format!(
                "{ident}\t{extended_start}\t{extended_end}\t{orient_ch}\tn={x},anchor_range={anchor_a}-{anchor_b}"
            );
            t.push_str(&if s.left_extendable {
                format!(",extended_left={}", s.left_extension_len)
            } else {
                ",extended_left=*".to_string()
            });
            t.push_str(&if s.right_extendable {
                format!(",extended_right={}", s.right_extension_len)
            } else {
                ",extended_right=*".to_string()
            });
            t.push_str(&format!(",len={extended_length},score={}\n", s.score));
            fp.write_all(t.as_bytes())?;
        }

        if let Some(fp) = fp_fa.as_mut() {
            if s.left_extension_len < 0 || s.right_extension_len < 0 {
                eprintln!("Error: Negative extension length detected");
                std::process::exit(1);
            }
            // Positions of the first and last base to emit, in the library
            // buffer. The C walks these inclusively, and its flank clamps
            // below are stated on positions, so they stay positions here.
            let (mut int_start, mut int_end) = (extended.start(), extended.end() - 1);
            let flanking = cfg.flanking;
            // Computed after the clamps below; `int_end` may equal the
            // subsequence boundary, which the C also reads.
            let emitted_1b = |a: u64, b: u64| {
                lib.to_source(seq_idx, ram_core::Span::new(a, b + 1).expect("a <= b"))
                    .as_1b_closed()
                    .expect("non-empty")
            };
            if flanking > 0 {
                // PORT NOTE: the C clamps to 1 (not 0) at a sequence start.
                if seq_lower == 0 && flanking as u64 > int_start {
                    int_start = 1;
                } else {
                    int_start -= flanking as u64;
                }
                if int_end + flanking as u64 > seq_upper {
                    int_end = seq_upper;
                } else {
                    int_end += flanking as u64;
                }
                writeln!(fp,
                    ">{ident}:{}-{}_{orient_ch}  n={x},anchor_range={}-{},extended_left={},extended_right={},len={extended_length},flanking={flanking},score={}",
                    emitted_1b(int_start, int_end).0,
                    emitted_1b(int_start, int_end).1,
                    s.left_seq_pos - seq_lower + 1,
                    s.right_seq_pos - seq_lower + 1,
                    s.left_extension_len,
                    s.right_extension_len,
                    s.score
                )?;
            } else {
                writeln!(fp,
                    ">{ident}:{}-{}_{orient_ch}  n={x},anchor_range={}-{},extended_left={},extended_right={},len={extended_length},score={}",
                    emitted_1b(int_start, int_end).0,
                    emitted_1b(int_start, int_end).1,
                    s.left_seq_pos - seq_lower + 1,
                    s.right_seq_pos - seq_lower + 1,
                    s.left_extension_len,
                    s.right_extension_len,
                    s.score
                )?;
            }
            let mut out = Vec::with_capacity((int_end - int_start + 2) as usize);
            if orient_ch == '-' {
                let mut j = int_end;
                loop {
                    out.push(num_to_char(compl(lib.sequence[j as usize])));
                    if j == int_start {
                        break;
                    }
                    // C decrements with a do-while on unsigned; int_start > 0
                    // always holds here in practice.
                    j -= 1;
                }
            } else {
                for j in int_start..=int_end {
                    out.push(num_to_char(lib.sequence[j as usize]));
                }
            }
            if out.is_empty() {
                eprintln!(
                    "Error: No sequence emitted for {ident}:{int_start}-{int_end}_{orient_ch}"
                );
                std::process::exit(1);
            }
            out.push(b'\n');
            fp.write_all(&out)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_wrapped;

    fn wrapped(master: &[u8], start: usize, end: usize, phase: usize) -> Vec<u8> {
        let mut out = Vec::new();
        write_wrapped(&mut out, master, start, end, phase);
        out
    }

    /// Fixed mode wraps a block at clean 80-column boundaries; C-compat mode
    /// (phase_origin = masterstart) shifts the first wrap by the phase.
    #[test]
    fn right_block_wraps_cleanly_in_fixed_mode_and_phased_in_compat() {
        // 100 A's (encoded 0) as a right-extension block starting at 30,
        // with a compat phase origin of 7 (i.e. leftbp+l = 23 earlier).
        let master = vec![0u8; 200];
        let fixed = wrapped(&master, 30, 130, 30);
        let lines: Vec<usize> = fixed.split(|&b| b == b'\n').map(|l| l.len()).collect();
        assert_eq!(lines, vec![80, 20, 0], "fixed: 80 then 20 then trailing");

        let compat = wrapped(&master, 30, 130, 7);
        let lines: Vec<usize> = compat.split(|&b| b == b'\n').map(|l| l.len()).collect();
        // First wrap fires when (x - 7) % 80 == 79, i.e. after 57 chars.
        assert_eq!(lines, vec![57, 43, 0], "compat: phase-shifted wrap");
    }

    /// A block ending exactly on a wrap must not emit a double newline.
    #[test]
    fn exact_multiple_of_eighty_has_single_trailing_newline() {
        let master = vec![1u8; 200];
        let out = wrapped(&master, 0, 160, 0);
        assert!(out.ends_with(b"C\n"), "no blank line at the end");
        assert_eq!(out.iter().filter(|&&b| b == b'\n').count(), 2);
    }
}
