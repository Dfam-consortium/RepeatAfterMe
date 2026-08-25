#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) extend-stk.pl  (reduced, drives ram-extend -stk)
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Reduction: Claude (2026-08) — items 1-4 of the original workflow
##      (Stockholm parsing, extendable-edge detection, adaptive matrix /
##      minimprovement selection, RAMExtend invocation) now live in the
##      `ram-extend -stk` binary; this script retains only the
##      RepeatModeler-specific steps: Refiner-based MSA refinement and
##      Stockholm reassembly with curation notes.
##---------------------------------------------------------------------------##

=head1 NAME

extend-stk - Use ram-extend to automatically extend RepeatModeler seed alignments

=head1 SYNOPSIS

extend-stk.pl  [-debug|d]
               [-msaout <*.out>]
               [-min_aligning_seqs|m #]
               [-onlyextend|e]
               [-outprefix <prefix>]
               [-repeatmodeler_dir <dir>]
               [-bandwidth|b #]
               -assembly|a <*.2bit>
               -input|i <*.stk>
               -output|o <*.stk>

=head1 DESCRIPTION

Reduced version of the original extend-stk.pl.  The heavy lifting —
Stockholm parsing, per-member extendable-flag derivation, Kimura-
divergence-based matrix and minimprovement selection, and the extension
itself — is one invocation of `ram-extend -stk`.  This script then, per
family: builds the refined MSA via RepeatModeler's
alignAndCallConsensus.pl and reassembles the output Stockholm record
with the original metadata plus curation notes.

Differences from the original, both deliberate:
  * duplicate extended sequences are actually removed before refinement
    (the original built repam-repseq-nodups.fa but then refined the
    un-deduplicated file);
  * all families run in a single ram-extend process (per-record outputs
    are suffixed `.recN`).

=cut

use strict;
use Getopt::Long;
use FindBin;
use File::Path 'rmtree';
use File::Temp qw/ tempdir /;
use Cwd;

my $Version = "0.4-reduced";

# ram-extend is expected beside this script's parent (workspace target),
# overridable with RAM_EXTEND.
my $ramExtend = $ENV{'RAM_EXTEND'} || "$FindBin::RealBin/../target/release/ram-extend";

my @getopt_args = (
    '-assembly|a=s', '-input|i=s', '-output|o=s',
    '-min_aligning_seqs|m=i', '-msaout=s', '-onlyextend|e',
    '-outprefix=s', '-repeatmodeler_dir=s', '-bandwidth|b=i', '-debug|d',
);
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
GetOptions( \%options, @getopt_args ) or usage();

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit(1);
}

my $DEBUG = $options{'debug'} ? 1 : 0;
usage() unless ( $options{'assembly'} && $options{'input'} && $options{'output'} );
die "Could not find ram-extend at $ramExtend (set RAM_EXTEND)\n" unless -x $ramExtend;

my $rmodDir = $ENV{'REPEATMODELER_DIR'} || "";
$rmodDir = $options{'repeatmodeler_dir'} if ( $options{'repeatmodeler_dir'} );
if ( !$options{'onlyextend'} ) {
  die "\nERROR: Could not find RepeatModeler (set REPEATMODELER_DIR or " .
      "-repeatmodeler_dir)\n" unless ( $rmodDir && -d $rmodDir );
  require lib;
  lib->import($rmodDir);
  require SeedAlignment;
  require SeedAlignmentCollection;
}

my $genomeFile = $options{'assembly'};
my $inputFile  = $options{'input'};
my $outputFile = $options{'output'};
my $mas        = $options{'min_aligning_seqs'} || 3;

my $re_version = "Unknown";
open my $VIN, "$ramExtend -version|" or die;
while (<$VIN>) { $re_version = $1 if (/RAMExtend Version (\S+)/); }
close $VIN;

##
## Step 1-4: one ram-extend -stk invocation for every family.
##
my $tdir = tempdir( DIR => cwd, CLEANUP => !$DEBUG );
my $bw = $options{'bandwidth'} ? "-bandwidth $options{'bandwidth'}" : "";
my $cmd = "$ramExtend -twobit $genomeFile -stk $inputFile $bw " .
          "-min_aligning_seqs $mas -outtsv $tdir/ranges.tsv " .
          "-outfa $tdir/repseq.fa -cons $tdir/cons.fa -vvv " .
          "> $tdir/repam.log 2>&1";
print "Running: $cmd\n" if $DEBUG;
system($cmd);

##
## Parse the per-family log sections: record number, label, skip
## messages, and the Extended left/right lines.
##
my %fam = ();    # recN -> { label, skipped, leftExt, rightExt }
{
  my $rec = 0;
  open my $LOG, "<$tdir/repam.log" or die "ram-extend produced no log\n";
  while (<$LOG>) {
    if (/^== Family record:(\d+)\/(.*) ==/) { $rec = $1; $fam{$rec}{'label'} = $2; }
    elsif (/\*\* skipped: (.*) \*\*/)       { $fam{$rec}{'skipped'} = $1; }
    elsif (/Extended right:\s+(\d+)/)       { $fam{$rec}{'rightExt'} = $1; }
    elsif (/Extended left\s*:\s+(\d+)/)     { $fam{$rec}{'leftExt'} = $1; }
  }
  close $LOG;
}
my $multi = ( scalar( keys %fam ) > 1 );
sub famFile { my ($base, $rec) = @_; return $multi ? "$base.rec$rec" : $base; }

##
## -onlyextend: save the raw artifacts (the combined consensus is
## already in the cons file's ">combined" record).
##
if ( $options{'onlyextend'} ) {
  my $prefix = $options{'outprefix'} || "extended";
  for my $rec ( sort { $a <=> $b } keys %fam ) {
    next if $fam{$rec}{'skipped'};
    my $id = $multi ? "$prefix-rec$rec" : $prefix;
    system("cp " . famFile("$tdir/ranges.tsv", $rec) . " $id-repam-ranges.tsv")
      if -e famFile("$tdir/ranges.tsv", $rec);
    system("cp " . famFile("$tdir/cons.fa", $rec) . " $id-ext-cons.fa")
      if -e famFile("$tdir/cons.fa", $rec);
    print "rec$rec ($fam{$rec}{'label'}): left $fam{$rec}{'leftExt'} bp, " .
          "right $fam{$rec}{'rightExt'} bp -> $id-*\n";
  }
  system("cp $tdir/repam.log " . ($options{'outprefix'} || "extended") . "-repam.log");
  exit(0);
}

##
## Steps 5-6 per family: refine the extended MSA and reassemble the
## Stockholm record.
##
my $seedColl = SeedAlignmentCollection->new();
open my $IN, "<$inputFile" or die "Could not open $inputFile\n";
$seedColl->read_stockholm($IN);
close $IN;

open my $OUTPUT, ">$outputFile" or die "Could not open $outputFile for writing\n";
my ( $extended, $unchanged, $segmental, $failed ) = ( 0, 0, 0, 0 );

for ( my $i = 0; $i < $seedColl->size(); $i++ ) {
  my $seedAlign = $seedColl->get($i);
  my $rec = $i + 1;
  my $info = $fam{$rec} || {};
  my $label = $info->{'label'} || $seedAlign->getId();
  print "\n== Family $label (record $rec) ==\n";

  my $leftExt  = $info->{'leftExt'}  || 0;
  my $rightExt = $info->{'rightExt'} || 0;

  if ( $info->{'skipped'} ) {
    print "  ** $info->{'skipped'} **\n";
    my $cnotes = $seedAlign->getCuratorComments();
    $cnotes .= "Too few extendable sequences for auto-extension\n";
    $seedAlign->setCuratorComments($cnotes);
    print $OUTPUT "" . $seedAlign->toString();
    $unchanged++;
    next;
  }
  if ( $leftExt == 0 && $rightExt == 0 ) {
    print "  ** Could not extend **\n";
    my $cnotes = $seedAlign->getCuratorComments();
    $cnotes .= "RAMExtend[$re_version]: left: 0 bp, right: 0bp\n";
    $seedAlign->setCuratorComments($cnotes);
    print $OUTPUT "" . $seedAlign->toString();
    $unchanged++;
    next;
  }
  # The original's segmental-duplication guard: limits hit both ways.
  if ( $leftExt > 9990 && $rightExt > 9990 ) {
    print "  ** Extension hit limits in both directions...probably a " .
          "segmental duplication...keeping unextended **\n";
    my $desc = $seedAlign->getDescription();
    $desc =~ s/[\n\r]+//g;
    $desc .= " [possibly part of segmental duplication]";
    $seedAlign->setDescription($desc);
    print $OUTPUT "" . $seedAlign->toString();
    $segmental++;
    next;
  }

  my $consFile = famFile( "$tdir/cons.fa", $rec );
  my $repFile  = famFile( "$tdir/repseq.fa", $rec );
  my $newCons  = "";
  if ( -e $consFile ) {
    # ram-extend -stk appends the assembled ">combined" record.
    open my $CIN, "<$consFile" or die;
    my $inCombined = 0;
    while (<$CIN>) {
      if (/^>combined/)   { $inCombined = 1; next; }
      if (/^>/)           { $inCombined = 0; next; }
      if ($inCombined)    { s/[\n\r\s]+//g; $newCons .= $_; }
    }
    close $CIN;
  }
  if ( $newCons eq "" || !-e $repFile ) {
    print "  ** extension outputs missing; keeping unextended **\n";
    print $OUTPUT "" . $seedAlign->toString();
    $failed++;
    next;
  }
  print "  - Estimated extensions: left $leftExt bp, right $rightExt bp\n";

  # Deduplicate identical extended sequences before refinement (fixed:
  # the original wrote the dedup file but refined the original).
  my $nodups = "$tdir/repseq-nodups-rec$rec.fa";
  my %seen = ();
  open my $RIN, "<$repFile" or die;
  open my $ROUT, ">$nodups" or die;
  my ( $head, $name, $seq ) = ( undef, undef, "" );
  my $flush = sub {
    return unless defined $head;
    my $key = "$name=$seq";
    if ( $seen{$key}++ ) {
      print "  Warning: ignored duplicate (identical sequence) for $name\n";
    } else {
      print $ROUT "$head\n$seq\n";
    }
    $seq = "";
  };
  while (<$RIN>) {
    chomp;
    if (/^>(\S+)/) { $flush->(); $head = $_; $name = $1; }
    else           { $seq .= $_; }
  }
  $flush->();
  close $RIN; close $ROUT;

  print "  - Rebuilding MSA with extensions...\n";
  # alignAndCallConsensus.pl names its outputs after the consensus
  # record's ID, not the -c filename — use a per-record ID so outputs
  # are unique in the shared temp dir and findable afterwards.
  my $newrep = "$tdir/newrep-rec$rec.fa";
  open my $NOUT, ">$newrep" or die;
  print $NOUT ">newrep-rec$rec\n$newCons\n";
  close $NOUT;
  my $rcmd = "$rmodDir/util/alignAndCallConsensus.pl -c $newrep -e $nodups " .
             "-refine 10 -st -q";
  print "     - Running $rcmd\n" if $DEBUG;
  system($rcmd);
  my $newrepBase = "$tdir/newrep-rec$rec";
  if ( $options{'msaout'} && -s "$newrepBase.out" ) {
    system( "cp $newrepBase.out " . $options{'msaout'} );
  }
  unless ( -s "$newrepBase.stk" ) {
    print "  ** refinement produced no Stockholm; keeping unextended **\n";
    my $cnotes = $seedAlign->getCuratorComments();
    $cnotes .= "RAMExtend[$re_version]: refinement failed after " .
               "left: $leftExt bp, right: $rightExt bp\n";
    $seedAlign->setCuratorComments($cnotes);
    print $OUTPUT "" . $seedAlign->toString();
    $failed++;
    next;
  }

  my $singleColl = SeedAlignmentCollection->new();
  open my $SIN, "<$newrepBase.stk" or die;
  $singleColl->read_stockholm($SIN);
  close $SIN;
  my $finalAlign = $singleColl->get(0);
  $finalAlign->setId( $seedAlign->getId() );
  $finalAlign->setDescription( $seedAlign->getDescription() );
  $finalAlign->clearCitations();
  my $tClass = $seedAlign->getClassification();
  $tClass = "Interspersed_Repeat;Unknown" if ( $tClass eq "Unknown" );
  $finalAlign->setClassification($tClass);
  $finalAlign->clearClades();
  for ( my $k = 0; $k < $seedAlign->cladeCount(); $k++ ) {
    $finalAlign->addClade( $seedAlign->getClade($k) );
  }
  my $cnotes = $seedAlign->getCuratorComments();
  $cnotes .= "RAMExtend[$re_version]: left: $leftExt bp, right: $rightExt bp\n";
  $finalAlign->setCuratorComments($cnotes);
  my $finalCons = $finalAlign->getRfLine();
  $finalCons =~ s/[\n\r\s\.\-]+//g;
  print "    - Final consensus length = " . length($finalCons) . "\n";
  print $OUTPUT "" . $finalAlign->toString();
  $extended++;
}
close $OUTPUT;

print "\nDone: $extended extended, $unchanged unchanged, " .
      "$segmental segmental-guard, $failed failed -> $outputFile\n";

1;
