# RepeatAfterMe

A package for the extension of repetitive DNA sequences

RepeatAfterMe is an experimental package containing (currently) 
a single tool called RAMExtend.  RAMExtend employs a method for
extending a multiple sequence alignment (MSA) that represents a 
fragment of a larger repetitive sequence found in a given genome.
The input MSA is treated as an immutable core while the flanking 
sequences are used to extend the MSA and generate a left/right
consensus extension for the original core.

The extension algorithm is an enhanced version of the RepeatScout 
approach developed by Alkes Price, Neil Jones and Pavel Pevzner 
(See history below).  Unlike the original RepeatScout program this 
tool may be used on complete genome sequences.  In addition it has been
re-written to use custom scoring matrices, and affine gap penalties.

An improved de novo tool is envisioned and will be made available
in the new future.

Robert Hubley 2022
Institute for Systems Biology


RAMExtend
=========

The tool minimally requires two input files.  The first is a file
the core alignment to extend from.  Ranges are supplied in the form 
of a modified BED-6 format:

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

There are additional options that are displayed if no options are specified.


History
=======
The genesis for this project was with the great work by Alkes Price,
Neil Jones and Pavel Pevzner.  They developed a method for automatically
extending aligned sequences ( abundant exact words ) allowing for 
the development of a dynamic multiple alignment.  The details of their
method are described in the following paper:

Price A.L., Jones N.C. and Pevzner P.A. 2005.  De novo identification of 
repeat families in large genomes.  To appear in Proceedings of the
13 Annual International conference on Intelligent Systems for
Molecular Biology (ISMB-05).  Detroit, Michigan.

One of the drawbacks for RepeatScout is the simple scoring system used
by the program ( Match/Mismatch/Gap ).  It was our initial goal to 
simply augment the RepeatScout package with custom matrix support and
full affine gap penalties.  As work progressed and the possibilities
of further enhancing the code became clear it was decided that a new
project should be created for this effort. RepeatAfterMe is designed
as an experimental workbench for the application of this powerful 
extension algorithm to various types of aligned cores (kmers, alignment
fragments etc).

