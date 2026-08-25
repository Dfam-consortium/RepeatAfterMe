#!/usr/local/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempfile tempdir /;

# NOTE: This tool has the following dependencies set in the code below:
#     RepeatAfter/
#     Phrap/
#     RepeatModeler

##
## Localization
##
my $ucscToolsDir = "/usr/local/ucscTools";
# 2.0.2a or higher
my $repeatModelerDir = "/usr/local/RepeatModeler-2.0.2a";
my $repeatAfterMeDir = "/u1/home/rhubley/projects/RepeatAfterMe/forDavid/RepeatAfterMe-0fef7fa7d69b960f9dc53cd1b8f8e8a50ce506eb";
my $phrapDir = "/usr/local/phrap";


my @getopt_args = (
                    '-genome=s',
                    '-family=s',
                    '-outdir=s',
                    '-div=s',
                    '-h'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  die; 
}

if ( $options{'h'} || ! -s $options{'family'} || ! -s $options{'genome'} || ! $options{'outdir'} ) {
  print "./davidExtendConsRAM.pl -genome <*.2bit> -family <catTEFile> -outdir <dir> [-div 14|18|20|25 (default:18)]\n";
  print "   This version of the script is customized for the Ray lab where they start with a 'catTEFile'\n";
  exit(0);
}

my $matrixDir = "$repeatModelerDir/Matrices";
my $catTEfile = $options{'family'};
my $outdir = $options{'outdir'};
my $genome_file = File::Spec->rel2abs($options{'genome'});

if ( ! -d $outdir ){
  mkdir $outdir;
}

print "davidExtendConsRAM: Starting. Working on $catTEfile\n";

# Make a 2bit version of the genome if we don't already have one
# This will be used to extend flanking sequences
if ( $genome_file =~ /^.*\.(fa|fasta|FASTA)$/  )
{
  print "Generating a 2bit version of $genome_file...\n";
  system("$ucscToolsDir/faToTwoBit $genome_file $genome_file.2bit");
  $genome_file .= ".2bit";
}

##
## Process a catTEfile
##
die "Could not find catTEfile $catTEfile\n" if ( ! -s $catTEfile);
open IN,"<$catTEfile" or die "Could not open $catTEfile for reading!\n";
my $repseq_data;
my $rep_data;
my $family_name;
my $raw_id;
my $seq;
my $hStr = "";
my $copies = 0;
while ( <IN> ){
  if ( /^>(\S+)/ ) {
    my $tmp_id = $1;
    if ( $seq ) {
      # >CONSENSUS-rnd-1_family-130__LINE___L1
      if ( $raw_id =~ /^CONSENSUS-(\S+).*$/ ) {
        $family_name = $1;
        $rep_data = ">$1\n$hStr$seq$hStr\n";
      # zero-based, half-open
      # >AllBul_scaffold_4546:33483-35601(+)
      }elsif ( $raw_id =~ /^(\S+):(\d+)-(\d+)\(([+-])\).*$/ ) {
        # Convert back to 1-based, fully-closed for Arian's software
        $repseq_data .= ">$1:". ($2 + 1 ) . "-$3";
        $repseq_data .= "_R" if ( $4 eq "-" );
        $repseq_data .= "\n$seq\n";
        $copies++;
      }
    }
    $raw_id = $tmp_id;
    $seq = "";
    next;
  }
  s/[\n\r\s]//g;
  $seq .= $_;
}
my $coreConsensus = "";
if ( $seq ) {
  # >CONSENSUS-rnd-1_family-130__LINE___L1
  if ( $raw_id =~ /^CONSENSUS-(\S+).*$/ ) {
    $family_name = $1;
    $rep_data = ">$1\n$hStr$seq$hStr\n";
    $coreConsensus = $seq;
  # >AllBul_scaffold_4546:33483-35601(+)
  }elsif ( $raw_id =~ /^(\S+):(\d+)-(\d+)\(([+-])\).*$/ ) {
    # Convert back to 1-based, fully-closed for Arian's software
    $repseq_data .= ">$1:". ($2 + 1 ) . "-$3";
    $repseq_data .= "_R" if ( $4 eq "-" );
    $repseq_data .= "\n$seq\n";
    $copies++;
  }
}
close IN;

# Generate initial set of alignments for RepeatAfterMe to work on
open  OUT,">$outdir/rep" or die "Could not open $outdir/rep for writing!\n";
print OUT $rep_data;
close OUT;
open  OUT,">$outdir/repseq" or die "Could not open $outdir/repseq for writing!\n";
print OUT $repseq_data;
close OUT;

die "Something went wrong creating rep file!\n" if ( ! -s "$outdir/rep" );

## TODO: Why don't we just call alignAndCallConsensus HERE.

#
# Generate out
#
my $div = 18;
$div = $options{'div'} if ( $options{'div'} );
my $minmatch = 7;
my $minscore = 200;
my $parameters;
if ($div == 25 ) {
  $parameters = "-M $matrixDir/crossmatch/25p41g.matrix -gap_init -25 -gap_ext -5 -minscore $minscore -minmatch $minmatch";
} elsif ($div == 20 ) {
  $parameters = "-M $matrixDir/crossmatch/20p41g.matrix -gap_init -28 -gap_ext -6 -minscore $minscore -minmatch $minmatch";
} elsif($div == 18 ) {
  $parameters = "-M $matrixDir/crossmatch/18p41g.matrix -gap_init -30 -gap_ext -6 -minscore $minscore -minmatch $minmatch";
} elsif($div == 14 ) {
  $parameters = "-M $matrixDir/crossmatch/14p41g.matrix -gap_init -33 -gap_ext -6 -minscore $minscore -minmatch $minmatch";
} else {
  die "Wrong divergence level. $0 only accepts 14, 18, 20 or 25\n";
}

print "Aligning elements to consensus using crossmatch with $parameters ...\n";
system "$phrapDir/cross_match $outdir/repseq $outdir/rep $parameters -alignments 2>/dev/null | egrep -v -a '[[:cntrl:]]' > $outdir/out";

if ( -s "$outdir/out" ){
  my $result = `fgrep "matching entries" $outdir/out`;
  if ( $result =~ /^0 matching/ )
  {
    die "Could not find any matching entries for rep in repseq!\n";
  }
}else {
    die "Crossmatch failed to generate an out file!\n";
}

my $curDir = `pwd`;
chdir($outdir);

my $padLength = 0;
my $genomeFile = $genome_file;

open IN,"<rep" or die;
my $coreCons = "";
while (<IN>){
  next if ( /^>/ );
  s/[\n\r\s]+//g;
  $coreCons .= $_;
}
close IN;
  
open IN,"$repeatModelerDir/util/Linup -stockholm out|" or die;
open OUT,">linup-ranges.tsv" or die;
my $id;
my $genomeStart;
my $genomeEnd;
my $aseq;
while ( <IN> ) {
  # 1-based, fully-closed
  if ( /^(\S+):(\d+)-(\d+)(_R)?:(\d+)-(\d+)\s+(\S+)/ ) { 
    $id = $1;
    my $outerStart = $2;
    my $outerEnd = $3;
    # These could also be reversed ( not sure why that happens in this dataset though )
    my $insideStart = $5;
    my $insideEnd = $6;
    my $insideOrient = "+";
    if ( $5 > $6 ) {
      $insideOrient = "-";
      $insideStart = $6;
      $insideEnd = $5;
    }
    my $or = $4;
    $aseq = $7;
    my $orient = "+";
    if ( $or ne "" ){ 
      if ( $insideOrient eq "+" ) {
      $orient = "-";
      }
      $genomeEnd = $outerEnd - $insideStart + 1; # zero-based, half-open
      $genomeStart = $outerEnd - $insideEnd; # zero-based
    }else {
      if ( $insideOrient eq "-" ) {
        $orient = "-"; 
      }
      $genomeStart = $outerStart + $insideStart - 1; # zero-based
      $genomeEnd = $outerStart + $insideEnd; # zero-based, half-open
    }

    my $soffset = 0;
    if ( $aseq =~ /^(\.+)/ ) {
      $soffset = length($1);
    }
    my $eoffset = 0;
    if ( $aseq =~ /[^\.](\.+)$/ ) {
      $eoffset = length($1);
    }
    if ( $soffset < 5 && $eoffset < 5 ) {
      print OUT "$id\t$genomeStart\t$genomeEnd\t1\t1\t$orient\n";
    }elsif ( $soffset < 5 ) 
    {
      # can extend left
      print OUT "$id\t$genomeStart\t$genomeEnd\t1\t0\t$orient\n";
    }elsif ( $eoffset < 5 ) 
    {
      # can extend right
      print OUT "$id\t$genomeStart\t$genomeEnd\t0\t1\t$orient\n";
    }else {
      # can't extend either
      print OUT "$id\t$genomeStart\t$genomeEnd\t0\t0\t$orient\n";
    }
  } 
}
close IN;
system ("rm out.malign") if ( -e "out.malign" );
system ("rm repseq.log");

system("mv rep rep.unextended");
system("mv repseq repseq.unextended");
system("mv out out.unextended");

my $cmd = "$repeatAfterMeDir/RAMExtend -twobit $genomeFile -ranges linup-ranges.tsv -bandwidth 500 -matrix $div" . "p43g -outtsv repam-ranges.tsv -outfa repam-repseq.fa -cons repam-cons.fa >& repam.log";
print "Running: $cmd\n";
system($cmd);
if ( -e "repam-cons.fa") {
 open IN,"<repam-cons.fa" or die;
 my $id;
 my %seqs = ();
 while (<IN>){
   if ( />(\S+)/ )
   {
     $id = $1;
     next;
   }
   s/[\n\r\s]+//g;
   $seqs{$id} .= $_;
 }
 close IN;
 open OUT,">repam-newrep.fa" or die;
 print OUT ">repam-newrep\n" . $seqs{'left-extension'} . $coreCons . $seqs{'right-extension'} . "\n";
 close OUT;
 system("cp repam-newrep.fa rep");
 system("cp repam-repseq.fa repseq");
 print "Calling alignAndCallConsensus after extension...\n";
 system("$repeatModelerDir/util/alignAndCallConsensus.pl -defaults ");
 print "Producing extended-cons.fa\n";
 system("$repeatModelerDir/util/Linup -consensus repam-newrep.out > extended-cons.fa");
 print "Producing MSA-extended.fa\n";
 system("$repeatModelerDir/util/Linup -malignOut out.malign.extended -genome $genomeFile -includeFlanking 50 -msa repam-newrep.out > MSA-extended.fa");
 # Save unaltered repseq cons to ...
 #   makes out, out.malign.extended, ali
 print "Backing up files out/ali/rep/repseq\n";
 system("mv repam-newrep.out out.extended");
 system("mv repam-newrep.ali ali.extended");
 system("cp rep rep.extended");
 system("cp repseq repseq.extended");

 # Produce a special alignment with the CORE consensus represented as one of the members.
 #   NOTE: This is useful for visualization but not for consensus calling.
 open OUT,">>repseq" or die;
 print OUT ">CORECONS\n$coreCons\n";
 close OUT;
 print "Calling alignAndCallConsensus with CORE consensus represented for visualization\n";
 system("$repeatModelerDir/util/alignAndCallConsensus.pl --defaults");
 system("$repeatModelerDir/util/Linup -malignOut out.malign.w_rmod_cons -i -genome $genomeFile -includeFlanking 50 -msa repam-newrep.out > MSA-extended_with_rmod_cons.fa");
 system("$repeatModelerDir/util/visualizeAlignPNG.pl -genome $genome_file -outfile repam-newrep.out");
 system("mv repam-newrep.ali ali.w_rmod_cons");
 system("mv repam-newrep.out out.w_rmod_cons");
 system("cp rep.extended rep");
 system("cp repseq.extended repseq");
 system("cp out.extended out");
 system("cp out.malign.extended out.malign");
}else {
 print "Could not extend $catTEfile.\n";
 system("cp rep.unextended rep");
 system("cp repseq.unextended repseq");
 system("cp out.unextended out");
 system("$repeatModelerDir/util/visualizeAlignPNG.pl -genome $genome_file -outfile out");
}

chdir($curDir);

1;
