#!/usr/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) extend-stk.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************

=head1 NAME

extend-stk - Use RAMExtend to automatically extend RepeatModeler seed alignments

NOTE: This tool is designed to be installed in the RepeatModeler/util subdirectory
and the path "RepeatAfterMeDir" defined in the script to point to the location
of the RAMExtend binary.

=head1 SYNOPSIS

extend-stk.pl  [-debug|d]
               [-msaout|m <*.out>]
               [-min_aligning_seqs|m #]
               [-save_tsv_files|s]
               -assembly|a <*.2bit>
               -input|i <*.stk>
               -output|o <*.stk>

=head1 DESCRIPTION

Use the RepeatAfterMe ExtendAlign tool to try and extend seed alignment
that likely fragements of longer TE.  This is often the case with the
output of de novo TE discovery programs such as RepeatScout and RECON
the primary algorithms used in RepeatModeler.  

The options are:

=over 4

=item -assembly <*.2bit>

The genome assembly in UCSC 2bit format which contains the sequences
referenced in the input seed alignment.

=item -input <*.stk>

A stockholm file containing one or more seed alignments.

=item -output <*.stk>

The file to hold the results of the seed alignments after extension.

=item -min_aligning_seqs|m #

The approximate number of sequences that are continuing to align
as extension progresses [default = 3].  This is converted into the 
-minimprovement parameter of RAMExtend using the formula:

  minimprovement = min_aligning_seqs * average_matrix_diagonal_score.

For genomes containing a large fraction of segmental duplications
it may be necessary to increase this value.

=item -msaout <*.out>

Optional output file for the refined *.out file.

=item -save_tsv_files|s

Save the TSV files generated from the seed alignment to the current working
directory.  This is useful for debugging and testing RAMExtend parameters.

=back

=head1 ALSO

RepeatMasker, RepeatModeler, Dfam.org

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Bootstrap Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin;
use File::Path 'rmtree';
use File::Temp qw/ tempfile tempdir /;
use Cwd;

my $Version = 0.3;
my $RepeatAfterMeDir = "$FindBin::RealBin/..";

# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-assembly|a=s',
                    '-input|i=s',
                    '-output|o=s',
                    '-min_aligning_seqs|m=i',
                    '-msaout|m=s',
                    '-onlyextend|e',
                    '-outprefix=s',
                    '-repeatmodeler_dir=s',
                    '-bandwidth|b=i',
                    '-save_tsv_files|s',
                    '-debug|d'
);
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit( 1 );
}

# Resolve dependencies
my $rmodDir = "";
if ( exists $ENV{'REPEATMODELER_DIR'} && -d $ENV{'REPEATMODELER_DIR'} ) {
  $rmodDir = $ENV{'REPEATMODELER_DIR'};
}
if ( $options{'repeatmodeler_dir'} && -d $options{'repeatmodeler_dir'} ) {
  $rmodDir = $options{'repeatmodeler_dir'};
}

if ( $rmodDir eq "" ) {
  die "\n\nERROR: Could not find RepeatModeler installation.\n" .
      "      Please either set the REPEATMODELER_DIR environment variable\n" .
      "      or use provide the path to RepeatModeler with the\n" .
      "      -repeatmodeler_dir command line option.\n\n";
}else { print "rmodDir = $rmodDir\n"; }

require lib;
lib->import($rmodDir);
require RepModelConfig;
RepModelCOnfig->import();
require MultAln;
require EMBL;
require SeedAlignment;
require SeedAlignmentCollection;

# RepeatModeler version
my $RMOD_VERSION = $RepModelConfig::VERSION;

my $DEBUG = 0;
$DEBUG = 1 if ( $options{'debug'} );

unless ( $options{'assembly'} ) {
  print "\n\n-assembly is a required field!\n\n";
  usage();
}

unless ( $options{'output'} ) {
  print "\n\n-output is a required field!\n\n";
  usage();
}

unless ( $options{'input'} ) {
  print "\n\n-input is a required field!\n\n";
  usage();
}

my $genomeFile = $options{'assembly'};
my $inputFile = $options{'input'};
my $outputFile = $options{'output'};
my $min_aligning_seqs = 3;
$min_aligning_seqs = $options{'min_aligning_seqs'} if ( $options{'min_aligning_seqs'} );


# Was 500 for Dfam/Server/extend-stk.pl
my $re_bandwidth = 40;
if ( $options{'bandwidth'} ) {
  $re_bandwidth = $options{'bandwidth'};
}
my $re_version = "Unknown";
open IN,"$RepeatAfterMeDir/RAMExtend -h|" or die;
while (<IN>) {
  # RAMExtend Version 0.0.4 - build 436 date 20241015
  if ( /RAMExtend Version (\S+)/ ) {
    $re_version = $1;
    last;
  }
}
close IN;

print "##\n";
print "## extend-stk.pl\n";
print "##\n";
print "##   Program Version      : $Version\n";
print "##   RAMExtend Version    : $re_version\n";
print "##   RepeatModeler Version: $RMOD_VERSION\n";
print "##   Genome               : $genomeFile\n";
print "##   Input                : $inputFile\n";
print "##   Output               : $outputFile\n";
print "##   Min Aligning Seqs    : $min_aligning_seqs\n";
print "##\n";

# Open up the stockholm file and read in data
my $stockholmColl = SeedAlignmentCollection->new();
open my $IN, "<$inputFile"
     or die "Could not open up stockholm file $inputFile for reading!\n";
$stockholmColl->read_stockholm( $IN );
close $IN;
if ( $stockholmColl->size() == 0 ) {
  die "No seed alignments found in $inputFile!\n";
}

my $CLEANUP = 1;
$CLEANUP = 0 if ( $DEBUG );

open OUTPUT, ">$outputFile" or die;

my $failed = 0;
my $segmental = 0;
my $triggerSave = 0;
for ( my $i = 0 ; $i < $stockholmColl->size() ; $i++ ) {
  my $seedAlign   = $stockholmColl->get( $i );
  my $id          = $seedAlign->getId();
  my $accession   = $seedAlign->getAccession();

  if ( $id eq "" ) {
    if ( $accession ne "" ) {
      $id = $accession;
    }else {
      $id = "Unnamed_Family";
    }
  }

  print "Working on $id..\n";
 
  # NOTE: This currently ignores the stockholm reference line and generates
  #       one by calling a consensus.
  my $tdiv = 14;
  my $mAlign = MultAln->new( seedAlignment => $seedAlign );

  # NOTE: Older stockholm files simply contained reference positions marked with "x" or "X". 
  #       Newer ones have the consensus pre-caculated and used instead of this sentinel. If
  #       we find "x"'s we recalculate the consensus.  
  my $reference = $mAlign->getReferenceSeq();
  if ( $reference =~ /x/i ) {
    $reference = $mAlign->consensus();
  }

  #my ( $null, $totDiv, $avgDiv ) = $mAlign->kimuraDivergence();
  # New API
  my ( $avgKimura, $avgKimuraMod, $totTransitions, $totTransitionsMod, $totTransversions, $wellCharacterized, $CpGSites, $numHigh ) =
     $mAlign->calculateKimuraDivergence( consensus => $reference );

  # Remove gaps from the reference sequence
  $reference =~ s/\-//g;
  
  # TODO: make command line override
  $tdiv = $avgKimura;
  
  # This overrides the multaln method
  # Divergence is coming from generateSeedAlignment.pl:
  # Source:gsa, mDiv=23.22, GCA_000699945.1_ASM69994v1_genomic.fna.2bit:108
  if ( $seedAlign->getDescription() =~ /mDiv=(\d+\.\d+)/ )
  {
    $tdiv = $1;
  }

  my $div = 14;
  my $diagAvg = 10;
  if ($tdiv >= 16) {
    $div = 18;
    if ($tdiv >= 19) {
      $div = 20;
      if ($tdiv >= 22.5) {
        $div = 25;
        $diagAvg = 9;
      }
    }
  }

  my $minimprovement = $min_aligning_seqs * $diagAvg;
 
  my $tdir = tempdir( DIR => cwd, CLEANUP => $CLEANUP );

  print "  - Temporary directory: $tdir\n";
  my $recalcConsLen = length($reference);
  print "  - Consensus length [recalculated]: $recalcConsLen\n";
  print "  - Kimura divergence: $avgKimura % (no CpG adjustment)\n";
  print "  - Instances: " . $seedAlign->alignmentCount() . "\n";

  open TSV,">$tdir/linup.tsv" or die "Could not open $tdir/linup.tsv for writing!\n";

  my $extendable_count = 0;
  for ( my $j = 0; $j < $seedAlign->alignmentCount(); $j++ ){
    my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
      $seedAlign->getAlignment( $j );

    # TODO: This disqualifies user-generated identifiers like gi|1234, perhaps this check should be an option we use internally
    if ( $sequenceName =~ /^gi\|\d+$/ ) {
      warn "File contains unresolved sequenceNames ($sequenceName) from RepeatModeler -- this is a RepeatModeler bug.  Deleting sequence.\n";
      next;   
    }
    $start--; # convert to zero-based half open

    # insertions...
    if ( 0 ) {
      if ( $assemblyName ne "" ){
        $sequenceName = $assemblyName.":".$sequenceName;
      }
    }
    # If insertions against assembly...insertions against insertions is already fine
    if ( 0 ) {
      print "Warning accounting for 50bp padding for 200Mammals insertions\n";
      $start -= 50;
      $end -= 50;
    }

    # TODO: Add sanity check for sequence ranges through comparison with the 2bit file to the seed alignment data

    # Tolerate up to 10bp from edge
    if ( $sequence =~ /^[\.]{0,10}[^\.]/ && $sequence =~ /[^\.][\.]{0,10}$/ ) {
      print TSV "$sequenceName\t".$start."\t$end\t1\t1\t$orient\n";
      $extendable_count++;
    }elsif ( $sequence =~ /^[\.]{0,10}[^\.]/ ) 
    {
      # can extend left
      print TSV "$sequenceName\t".$start."\t$end\t1\t0\t$orient\n";
      $extendable_count++;
    }elsif ( $sequence =~ /[^\.][\.]{0,10}$/ ) 
    {
      # can extend right
      print TSV "$sequenceName\t".$start."\t$end\t0\t1\t$orient\n";
      $extendable_count++;
    }else {
      print TSV "$sequenceName\t".$start."\t$end\t0\t0\t$orient\n";
    }
  } 
  close TSV;

  my $t1 = "Too few sequences";
  my $t2 = "Too few sequences";
  my $cmd = "";
  if ( $extendable_count > 3 ) {
    $cmd = "$RepeatAfterMeDir/RAMExtend -twobit $genomeFile -L 20000 -bandwidth $re_bandwidth -matrix $div" . "p43g -ranges $tdir/linup.tsv -outtsv $tdir/repam-ranges.tsv -outfa $tdir/repam-repseq.fa -cons $tdir/repam-cons.fa -vvv -minimprovement $minimprovement > $tdir/repam.log 2>&1";
    print "  - Running RAMExtend [bandwidth=$re_bandwidth, matrix=$div" . "p43g, minimprovement=$minimprovement]..\n";
    print "     - Command: $cmd\n" if ( $DEBUG );
    if ( $options{'save_tsv_files'} ) {
      system("cp $tdir/linup.tsv $id.tsv");
    }
    system($cmd);
    my $retCode = $? >> 8;
    if ( $retCode != 0 ) {
      print "  RAMExtend failed! [$retCode]:\n";
      my $log = `cat $tdir/repam.log`;
      print "  repam.log: $log\n";
      die;
    }
  }else {
    print "  **Too few extendable sequences ($extendable_count) for RAMExtend**\n";
    my $cnotes = $seedAlign->getCuratorComments();
    $cnotes .= "Too few extendable sequences for auto-extension\n";
    $seedAlign->setCuratorComments( $cnotes );
    print OUTPUT "" . $seedAlign->toString();
    unless($DEBUG) {
      rmtree([ $tdir ]);
    }
    next;
  }

  if ( $options{'onlyextend'} ) {
    $t1 = `fgrep "Extended right:" $tdir/repam.log`;
    $t2 = `fgrep "Extended left :" $tdir/repam.log`;
    if ( $options{'outprefix'} ) {
      $id = $options{'outprefix'};
    }
    open IN,"<$tdir/repam-cons.fa" or die;
    open OUT,">$id-combined-cons.fa" or die;
    my $found_right = 0;
    print OUT ">combined\n";
    while (<IN>) {
      if ( /^>left-extension/ ) {
        next;
      }
      if ( /^>right-extension/ ) {
        $found_right = 1;
        print OUT "$reference\n";
        next;
      }
      print OUT $_;
    }
    if ( ! $found_right ) {
      print OUT "$reference\n";
    }
    close IN;
    close OUT;
    system("mv $tdir/repam.log $id-repam.log");
    system("mv $tdir/linup.tsv $id-linup.tsv");
    system("mv $tdir/repam-ranges.tsv $id-repam-ranges.tsv") if ( -e "$tdir/repam-ranges.tsv" );
    system("mv $tdir/repam-cons.fa $id-ext-cons.fa") if ( -e "$tdir/repam-cons.fa" );
    print "RAMExtend Results:\n  $t2  $t1\n";
    unless($DEBUG) {
      rmtree([ $tdir ]);
    }
    print "\nExtension Only: files saved to $id-repam-ranges.tsv, $id-linup.tsv, $id-ext-cons.fa, $id-repam.log $id-combined-cons.fa\n\n";
    next;
  }
    
  #Program duration is 4.0 sec = 0.1 min = 0.0 hr
  my $newCons = "";
  
  # Test for no extension
  $t1 = `fgrep "Extended right: 0 bp" $tdir/repam.log`;
  $t2 = `fgrep "Extended left : 0 bp" $tdir/repam.log`;
  if ( $t1 && $t2 ) {
    #Extended right: 0 bp
    #Extended left : 0 bp
    # TODO: Use this example for previous skip case
    print "  **Could not extend**\n";
    my $cnotes = $seedAlign->getCuratorComments();
    $cnotes .= "RAMExtend[$re_version]: left: 0 bp, right: 0bp\n";
    $seedAlign->setCuratorComments( $cnotes );
    print OUTPUT "" . $seedAlign->toString();
    unless($DEBUG) {
      rmtree([ $tdir ]);
    }
    next;
  }elsif ( -e "$tdir/repam-cons.fa") {
    open IN,"<$tdir/repam-cons.fa" or die;
    my $ttid;
    my %seqs = ();
    while (<IN>){
      if ( />(\S+)/ )
      {
        $ttid = $1;
        next;
      }
      s/[\n\r\s]+//g;
      $seqs{$ttid} .= $_;
    }
    close IN;
    $newCons = $seqs{'left-extension'} . $reference . $seqs{'right-extension'};
    #print "newCons = $newCons\n";
  }else {
    print "Something went wrong with extend align:\n";
    print "  $cmd\n";
    $failed++;
    $triggerSave = 1;
    next;
  } 

  ##
  ## 
  ##
  $t1 = `fgrep "Extended right:" $tdir/repam.log`;
  $t2 = `fgrep "Extended left :" $tdir/repam.log`;
  my $rightExt = 0;
  if ( $t1 =~ /Extended\s+(left|right)\s*:\s+(\d+)/ ){
    $rightExt = $2;
  }
  my $leftExt = 0;
  if ( $t2 =~ /Extended\s+(left|right)\s*:\s+(\d+)/ ){
    $leftExt = $2;
  }
  print "    - Estimated extensions: left $leftExt bp, right $rightExt bp, total " . ($leftExt + $rightExt) . "\n";
  if ( $leftExt > 9990 && $rightExt > 9990 )
  {
    print "    - Extension hit limits in both directions...probably a segmental duplication...keeping unextended\n";
    $segmental++;
    my $desc = $seedAlign->getDescription();
    $desc =~ s/[\n\r]+//g;
    $desc .= " [possibly part of segmental duplication]";
    $seedAlign->setDescription( $desc );
    print OUTPUT "" . $seedAlign->toString();
    # Cleanup
    unless($DEBUG) {
      rmtree([ $tdir ]);
    }                
    next;
  }

  open OUT,">$tdir/repam-newrep.fa" or die;
  print OUT ">repam-newrep\n$newCons\n";
  close OUT;
  #print "  extended by " . (length($newCons)-length($reference)) . " bp\n";
  my $prevCons;

  # RAMExtend sometimes produces duplicate identical sequences (same position, same nucleotides,
  # different source locations). Report and remove these cases.
  #
  # TODO: This could be further generalized, e.g. deduplicating sequence ranges that overlap by >80%.
  #
  open my $repseq_in, "<$tdir/repam-repseq.fa" or die;
  open my $nodups_out, ">$tdir/repam-repseq-nodups.fa" or die;
  my %seen_seqs = ();
  my $seqhead = undef;
  my $seqname = undef;
  my $seq = '';
  while (<$repseq_in>) {
    chomp;
    if (/^>(\S+)/) {
      if (defined $seqhead) {
        my $key = "$seqname=$seq";
        if (exists $seen_seqs{$key}) {
          print "Warning: ignored duplicate (identical sequence) for $seqname\n";
        } else {
          print $nodups_out "$seqhead\n$seq\n";
          $seen_seqs{$key} = 1;
        }
        $seq = '';
      }

      $seqhead = $_;
      $seqname = $1;
    } else {
      $seq .= $_;
    }
  }
  # NB: duplicated code to handle last sequence in file
  if (defined $seqhead) {
    my $key = "$seqname=$seq";
    if (exists $seen_seqs{$key}) {
      print "Warning: ignored duplicate (identical sequence) for $seqname\n";
    } else {
      print $nodups_out "$seqhead\n$seq\n";
      $seen_seqs{$key} = 1;
    }
    $seq = '';
  }
  close $repseq_in;
  close $nodups_out;

  print "  - Rebuilding MSA with extensions...\n";
  my $maxIter = 10;
  my $cmd = "$rmodDir/util/alignAndCallConsensus.pl -c $tdir/repam-newrep.fa -e $tdir/repam-repseq.fa -refine $maxIter -st -q";
  print "     - Running $cmd\n" if ( $DEBUG );
  system($cmd);
  # User wants the out file saved
  if ( $options{'msaout'} && -s "$tdir/repam-newrep.out" ) {
    system("cp $tdir/repam-newrep.out " . $options{'msaout'});
  }
  # TODO: Add exit case (similar to above patterns) whereby, if we don't generate a stockholm file, we should report produce an error ( perhaps mention in curation notes ) and move on to next family.
  #
  # TODO: Catch errors from alignAndCallConsensus.pl -- namely the one where it can't find the reference sequence in the alignment output.
  #

  # Now rebuild the stockholm entry and save to ...
  my $singleSeedColl = SeedAlignmentCollection->new();
  open my $INN, "<$tdir/repam-newrep.stk"
    or die "Could not open up stockholm file $tdir/repam-newrep.stk for reading!\n";
  $singleSeedColl->read_stockholm( $INN );
  close $INN;
  my $finalAlign   = $singleSeedColl->get( 0 );
  $finalAlign->setId( $seedAlign->getId() );
  $finalAlign->setDescription( $seedAlign->getDescription() );
  $finalAlign->clearCitations();
  # Correct a strange mapping artefact in RepeatModeler data
  my $tClass = $seedAlign->getClassification();
  if ( $tClass eq "Unknown" ) {
    $tClass = "Interspersed_Repeat;Unknown";
  }
  $finalAlign->setClassification( $tClass );
  $finalAlign->clearClades();
  for ( my $k = 0; $k < $seedAlign->cladeCount(); $k++ ){
    $finalAlign->addClade( $seedAlign->getClade($k) );
  }
  my $cnotes = $seedAlign->getCuratorComments();
  $cnotes .= "RAMExtend[$re_version]: left: $leftExt bp, right: $rightExt bp\n";
  $finalAlign->setCuratorComments( $cnotes );
  my $finalCons = $finalAlign->getRfLine();
  $finalCons =~ s/[\n\r\s\.\-]+//g;
  print "    - Final consensus length = " . length($finalCons) . " [ " . (length($finalCons) - $recalcConsLen) . " bp change ]\n";

  print OUTPUT "" . $finalAlign->toString();

  # Cleanup
  unless($DEBUG) {
    rmtree([ $tdir ]);
  }
}
close OUTPUT;


sub readConsensusFromFile {
   my $file = shift;
   open IN,"<$file" or die;
   my $seq;
   while (<IN>){
     if ( />\S+/ )
     {
       next;
     }
     s/[\n\r\s]+//g;
     $seq .= $_;
   }
   close IN;
   return $seq;
}


1;
