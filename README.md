# RepeatAfterMe

A package for the extension of repetitive DNA cores

The RepeatAfterMe RAMExtend tool automatically extends a multiple
sequence alignment (MSA) that may represent only a fragment of a longer
repetitive sequence family. Provided that the MSA includes detailed
information on the coordinates for each sequence, RAMExtend will 
perform a local alignment extension of the MSA using flanking sequence.
The consensus sequence of both the left and right extension is
generated and optionally the full set of extended sequences is
output in FASTA format.

The extension algorithm is an enhanced version of the RepeatScout 
approach developed by Alkes Price, Neil Jones and Pavel Pevzner 
(See history below).  The new algorithm supports multiple scoring
schemes, and affine gap penalties for greater sensitivity.  In addition, 
the tool attempts to detect satellites and attempts to avoid extending 
these sequences beyond one unit.  

Robert Hubley 2022-2025
Institute for Systems Biology


RAMExtend
=========

The tool minimally requires two input files.  The first is a file
containing the core alignment to extend from.  Ranges are supplied in 
the form of a modified BED-6 format:

      BED-6 field name  : RAMExtend use
      field-1:chrom     : sequence identifier
      field-2:chromStart: lower aligned position ( 0 based )
      field-3:chromEnd  : upper aligned position ( 0 based, half open )
      field-4:name      : left extendable flag ( 0 = no, 1 = yes )
      field-5:score     : right extendable flag 
      field-6:strand    : strand ( '+' = forward, '-' = reverse )

The fields are tab separated. Coordinates are zero-based, half-open. 
The 'extendable?' flags are used to limit the use of individual sequences
in the left or right extension phases.  This is useful when sequences
in the core MSA are not of uniform length.

The second input file is the genome itself in 2bit format.  This is where
the tool extracts the flanking regions using the identifiers provided
in the BED file.

For example:

  ./RAMExtend -ranges test/extension-test2.tsv -twobit test/extension-test2.2bit

The output will go to the screen by default, however it may be redirected to
files using various options:

  ./RAMExtend -ranges test/extension-test2.tsv -twobit test/extension-test2.2bit \
              -cons consensuses.fasta -outtsv ext_ranges.tsv -outfa ext_sequences.fasta

Additional help is displayed if no options are specified.


extend-stk.pl
=============

A wrapper script for RAMExtend that will automatically extend and refine
RepeatModeler seed alignments in Stockholm format.

For example:

  ./util/extend-stk.pl -assembly test/ce10.2bit \
                       -input test/ce10-fam1.stk \
                       -output ce10-fam1-extended.stk


    ##
    ## extend-stk.pl
    ##
    ##   Program Version      : 0.2
    ##   RAMExtend Version    : 0.0.6-dev
    ##   RepeatModeler Version: 2.0.6
    ##   Genome               : ce10.2bit
    ##   Input                : ce10-fam1.stk
    ##   Output               : /dev/null
    ##   Min Aligning Seqs    : 5
    ##
    Working on rnd-1_family-45..
      - Temporary directory: /u3/home/rhubley/notebooks/2024/1004-repeat_after_me_satellites/RepeatAfterMe-dev/test/CYmEAyYOQt
      - Consensus length [recalculated]: 163
      - Kimura divergence: 15.16 % (no CpG adjustment)
      - Instances: 100
      - Running RAMExtend [bandwidth=40, matrix=14p43g, minimprovement=50]..
        - Estimated extensions: left 354 bp, right 86 bp, total 440
      - Rebuilding MSA with extensions...
        - Final consensus length = 603 [ 440 bp change ]



History
=======
RepeatAfterMe traces it roots to the pioneering work of Alkes Price, 
Neil Jones and Pavel Pevzner who introduced an automated way to detect
repetitive DNA sequences by building and extending a multiple sequence
alignment centered on abundant exact words (k-mers). The details of their
method are described in the following paper:

Price A.L., Jones N.C. and Pevzner P.A. 2005.  De novo identification of 
repeat families in large genomes.  To appear in Proceedings of the
13 Annual International conference on Intelligent Systems for
Molecular Biology (ISMB-05).  Detroit, Michigan.

One of the drawbacks for RepeatScout is the simple scoring system used
by the program ( Match/Mismatch/Gap ).  It was our initial goal to 
simply augment the RepeatScout package with custom matrix support and
full affine gap penalties.  As work progressed and the possibilities
of further enhancing the code became clear, we decided that a new
project should be created for this effort. RepeatAfterMe is designed
as an experimental workbench for the application of this powerful 
extension algorithm to various types of aligned cores (kmers, alignment
fragments etc).  RAMExtend and extend-stk.pl are the first tools to be
based on this updated set of algorithms.


