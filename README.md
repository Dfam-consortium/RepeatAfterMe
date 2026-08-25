# RepeatAfterMe

A package for the extension of repetitive DNA cores.

`RAMExtend` takes a multiple sequence alignment that may cover only a fragment
of a longer repetitive family and extends it outward using flanking genomic
sequence. Given per-sequence coordinates, it performs a local alignment
extension, generates a consensus for the left and right extensions, and can
write the full set of extended sequences as FASTA.

The extension algorithm builds on the RepeatScout approach of Alkes Price,
Neil Jones, and Pavel Pevzner (see History). It adds multiple scoring schemes
and affine gap penalties for sensitivity, and it detects satellites and tries
not to extend past a single unit.

Robert Hubley, 2022-2026, Institute for Systems Biology.

## Implementations

Releases through `RepeatAfterMe_V0.0.7` were C. Those sources are frozen under
[`c/`](c/) and still build; the tags that produced them are unchanged, so
anything pinned to a release tarball keeps resolving to the same bytes.

This tree is a Rust reimplementation. It is flag-compatible with C v0.0.7,
produces byte-identical `-outtsv`, `-cons`, and `-outfa` output, and prints the
same stdout lines that RepeatModeler's `Refiner` scrapes. `harness/diff-c-rust.sh` checks that claim against the
C binary built from `c/`.

## Building

```bash
make            # builds the Rust RAMExtend at the top of the tree
make install    # installs to /usr/local/RepeatAfterMe-<version>
make c          # builds the frozen C tool under c/
make test       # unit tests plus the C-vs-Rust differential harness
```

Building the Rust tool needs a Rust toolchain (built and tested against 1.86;
no minimum is declared) and network access on the first build, since two
dependencies resolve from git tags. Building the C tool needs a GNU11 C
compiler and `make`.

Cargo names the binary `ram-extend`. RepeatModeler probes for `RAMExtend` and
dfam-tetools runs `RAMExtend -version`, so `make` leaves that name at the top
of the tree and `make install` installs both. They are the same binary.

## RAMExtend

Two inputs: a core alignment range file in modified BED-6, and a genome in
UCSC 2bit format.

| BED field  | Use by RAMExtend                            |
|------------|---------------------------------------------|
| chrom      | Sequence identifier                         |
| chromStart | Lower aligned position (0-based)            |
| chromEnd   | Upper aligned position (half-open)          |
| name       | Left extendable flag (`0` = no, `1` = yes)  |
| score      | Right extendable flag (`0` = no, `1` = yes) |
| strand     | `'+'` forward, `'-'` reverse                |

Fields are tab-separated, coordinates 0-based half-open. The extendable flags
say whether a given sequence takes part in left and right extension, which
matters for core-alignment fragments that do not reach the core edges.

```bash
./RAMExtend -ranges c/test/extension-test2.tsv \
            -twobit c/test/extension-test2.2bit \
            -cons consensuses.fasta \
            -outtsv ext_ranges.tsv \
            -outfa ext_sequences.fasta
```

Run it with no arguments for the option list. [`c/README.md`](c/README.md)
documents every option in full.

### Stockholm input

The Rust tool reads RepeatModeler seed alignments directly:

```bash
./RAMExtend -twobit c/test/ce10.2bit -stk c/test/ce10-fam1.stk
```

This covers the first three jobs of `extend-stk.pl`: parsing the alignment,
deriving extendable flags from its edges, and choosing a matrix and
`-minimprovement` from the alignment's Kimura divergence. Files holding several
families process in one invocation, with output suffixed by record label. The script at [`util/extend-stk.pl`](util/extend-stk.pl)
still rebuilds the MSA.

## Layout

- `ram-core/`: the library. Numeric alphabet, scoring systems
  (14/18/20/25p43g and RepeatScout), ranges and flank loading, the banded DP
  extension engine, and seed-anchored glocal alignment.
- `ram-cli/`: the `RAMExtend` binary and its Stockholm front end.
- `harness/`: `diff-c-rust.sh` runs C and Rust over the `c/test/` families
  under four parameter sets and compares byte-for-byte. `stk2ranges.py`
  converts Stockholm to ranges TSV. `truncation-bench.py` trims curated
  alignments and scores how much of each end the tool recovers.
- `util/extend-stk.pl`: the reduced driver that calls `RAMExtend -stk`. The
  original full script is at `c/util/extend-stk.pl`.
- `c/`: the frozen C implementation through V0.0.7, with its own README,
  Makefile, and test suite.

## History

RepeatAfterMe traces its roots to the work of Alkes Price, Neil Jones, and
Pavel Pevzner, who developed an automated method for detecting repetitive DNA
by building and extending multiple sequence alignments anchored on abundant
exact words:

> Price A.L., Jones N.C., Pevzner P.A. (2005).
> _De novo identification of repeat families in large genomes_.
> *Proceedings of the 13th Annual International Conference on Intelligent
> Systems for Molecular Biology (ISMB-05)*, Detroit, Michigan.

RepeatScout scored with a simple match/mismatch/gap model. This project began
as an attempt to give it affine gap penalties and custom scoring matrices, and
became a separate tool as the changes grew. `RAMExtend` and `extend-stk.pl`
were the first tools built on the result.

## License

CC0 1.0 Universal. See [LICENSE](LICENSE). The bundled UCSC kent sources under
`c/kentsrc/` and the minunit header under `c/minunit/` carry their own licenses.
