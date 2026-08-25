# RepeatAfterMe

A package for the extension of repetitive DNA cores

The **RepeatAfterMe** RAMExtend tool automatically extends a multiple
sequence alignment (MSA) that may represent only a fragment of a longer
repetitive sequence family. Provided that the MSA includes detailed
information on the coordinates for each sequence, RAMExtend will 
perform a local alignment extension of the MSA using flanking sequence.
The tool generates consensus sequences for both the left and right extensions 
and can optionally output the full set of extended sequences in FASTA format.

The extension algorithm is an enhanced version of the RepeatScout 
approach developed by Alkes Price, Neil Jones and Pavel Pevzner 
(See _History_ below).  Enhancements include:

- Support for multiple scoring schemes
- Affine gap penalties for improved sensitivity
- Satellite detection to prevent extending repeats beyond a single unit

**Author**: Robert Hubley (2022–2025)  
**Institution**: Institute for Systems Biology

## RAMExtend

### Overview

`RAMExtend` requires two input files:

1. **Core alignment ranges** (in a modified BED-6 format)
2. **Genome sequence** (in UCSC `.2bit` format)

The range file defines the alignment region for each sequence and which ends 
are eligible for extension.

#### Modified BED-6 Format

| BED Field | Use by RAMExtend                             |
|-----------|----------------------------------------------|
| chrom     | Sequence identifier                          |
| chromStart| Lower aligned position (0-based)             |
| chromEnd  | Upper aligned position (half-open)           |
| name      | Left extendable flag (`0` = no, `1` = yes)   |
| score     | Right extendable flag (`0` = no, `1` = yes)  |
| strand    | `'+'` for forward, `'-'` for reverse         |

All fields are **tab-separated** and use **0-based, half-open** coordinate 
conventions. The extendable flags control whether individual sequences 
participate in left/right extension—this is especially useful for sequence
fragments in the core alignment which are not proximal to the core
alignment edges.

#### Genome File

The genome must be in [UCSC 2bit format](https://genome.ucsc.edu/goldenPath/help/twoBit.html), 
and should contain the sequences referenced by the BED file.

---

### Example Usage

```bash
./RAMExtend -ranges test/extension-test2.tsv -twobit test/extension-test2.2bit
```

To output results to files:

```bash
./RAMExtend -ranges test/extension-test2.tsv -twobit test/extension-test2.2bit \
            -cons consensuses.fasta \
            -outtsv ext_ranges.tsv \
            -outfa ext_sequences.fasta
```

If no options are provided, usage help will be displayed.

---

## extend-stk.pl

A wrapper script for `RAMExtend` that automatically extends and refines 
RepeatModeler seed alignments in **Stockholm** format.  This helper
script translates the Stockholm file into a TSV file for RAMExtend,
invokes the tool, combines the left/right consensuses with the core 
alignment consensus, and rebuilds a new multiple sequence alignment 
using the core and extended sequences.

### Example Usage

```bash
./util/extend-stk.pl -assembly test/ce10.2bit \
                     -input test/ce10-fam1.stk \
                     -output ce10-fam1-extended.stk
```

Example log:

```
##
## extend-stk.pl
##
##   Program Version      : 0.2
##   RAMExtend Version    : 0.0.6
##   RepeatModeler Version: 2.0.6
##   Genome               : test/ce10.2bit
##   Input                : test/ce10-fam1.stk
##   Output               : ce10-fam1-extended.stk
##   Min Aligning Seqs    : 3
##
Working on rnd-1_family-45...
  - Temporary directory: /u3/home/rhubley/projects/RepeatAfterMe-public/FDLO9azqC5
  - Consensus length [recalculated]: 163
  - Kimura divergence: 15.16% (no CpG adjustment)
  - Instances: 100
  - Running RAMExtend [bandwidth=40, matrix=14p43g, minimprovement=30]...
    - Estimated extensions: left 354 bp, right 86 bp, total 440
  - Rebuilding MSA with extensions...
    - Final consensus length = 603 [440 bp change]
```

---

## Installation

### Requirements

To build **RepeatAfterMe**, you will need:

- A C compiler that supports GNU11 (e.g. `gcc`)
- `make`
- POSIX-compatible environment (e.g., Linux/macOS)

### Building

To compile the `RAMExtend` binary and the test suite:

```bash
make
```

This will:

- Compile the `RAMExtend` executable
- Build `kentsrc/libTwoBit.a` from source (included)
- Build a `test_suite` for validation and regression testing

### Running Tests

Run the internal test suite:

```bash
./test_suite
```

### Installing

Install the compiled binary and README to a system-wide or user-defined directory:

```bash
make install
```

By default, this installs into:

```
/usr/local/RepeatAfterMe-0.0.7/
```

You can add the binary to your `PATH` by:

```bash
export PATH=/usr/local/RepeatAfterMe-0.0.7/bin:$PATH
```

> To change the install location, edit the `INSTDIR` variable near the top of the `Makefile`.

### Cleaning the Build

To remove compiled artifacts:

```bash
make clean
```

---

## History

RepeatAfterMe traces its roots to the pioneering work of Alkes Price, 
Neil Jones, and Pavel Pevzner, who developed an automated method for 
detecting repetitive DNA by building and extending multiple sequence 
alignments based on abundant exact words (k-mers). Their method was 
introduced in:

> Price A.L., Jones N.C., Pevzner P.A. (2005).  
> _De novo identification of repeat families in large genomes_.  
> *Proceedings of the 13th Annual International Conference on Intelligent Systems for Molecular Biology (ISMB-05)*, Detroit, Michigan.

One limitation of the original **RepeatScout** was its use of a simple 
match/mismatch/gap scoring system. This project began as an effort to 
add affine gap penalties and custom scoring matrices to RepeatScout, 
but evolved into a standalone tool when the scope of improvements grew.

**RepeatAfterMe** is intended as an experimental platform for applying 
this enhanced extension algorithm to various core types (e.g., k-mers, 
partial alignments). `RAMExtend` and `extend-stk.pl` are the first 
tools to leverage these algorithms.

