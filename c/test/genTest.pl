#!/usr/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) genTest.pl
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
#
#

=head1 NAME

genTest.pl - Generate a test file for a given pattern

=head1 SYNOPSIS

genTest.pl -pattern <pattern> [-seed <seed> | -random] -out <file>

=head1 DESCRIPTION

The options are:

=over 4

=item -pattern <ltr_aligned|ltr_misaligned|ltr_random|satellite>

The pattern to simulate

=item -seed #

The random number generator seed to use.

=item -random

Use time() to set the random number generator seed. This will produce a different result each time the program is run.

=back

=head1 ALSO

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-pattern=s',   
                    '-seed=s',
                    '-random',
                    '-out=s',
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) )
{
  usage();
}

sub usage
{
  exec "pod2text $0";
  exit;
}

if ( -s $options{'out'} || -d $options{'out'} )
{
  print "Directory or file $options{'out'} exists.  Please remove and rerun.\n";
  exit 1;
}

if ( !defined $options{'pattern'} )
{
  print "Pattern is required.\n";
  usage();
}

my $seed = 123456789;
if ( exists $options{'seed'} )
{
  $seed = $options{'seed'}; 
}elsif ( exists $options{'random'} )
{
  $seed = time();
}
srand( $seed );

print "#\n";
print "# genTest.pl\n";
print "#\n";
print "#   pattern: $options{'pattern'}\n";
print "#   seed: $seed\n";
print "#   out: $options{'out'}\n";

my @bases = ('A', 'C', 'G', 'T');
my $outDir = $options{'out'};
mkdir( $outDir, 0777 ) || die "Can't create directory $outDir: $!\n";
if ( $options{'pattern'} eq "satellite" )
{
  my $flanking_len = 1000;
  my $satellite_len = 1000;
  my $satellite_num = 5;
  my $core_pos = 250;
  my $core_len = 10;

  my $canonical_satellite;
  $canonical_satellite .= $bases[int rand scalar @bases] for 1..$satellite_len;
  my $final_seq;
  ## flanking + satellite + flanking
  # left flank
  my $random_seq;
  $random_seq .= $bases[int rand scalar @bases] for 1..$flanking_len;
  $final_seq = $random_seq;
  # satellite x $satellite_num
  $final_seq .= ($canonical_satellite)x$satellite_num;
  # right flank
  $random_seq = "";
  $random_seq .= $bases[int rand scalar @bases] for 1..$flanking_len;
  $final_seq .= $random_seq;
  open OUT,">$outDir/test.fasta" or die "Can't open $outDir/test.fasta: $!\n";
  print OUT ">seq1\n$final_seq\n";
  close OUT;
  system("/usr/local/ucscTools/faToTwoBit $outDir/test.fasta $outDir/test.2bit");

  my $pos = $flanking_len + $core_pos;
  open OUT,">$outDir/test.tsv" or die "Can't open $outDir/test.tsv: $!\n";
  for ( my $i = 0; $i < $satellite_num; $i++ ) {
    print OUT "seq1\t".($pos-1)."\t".($pos+$core_len)."\t1\t1\t+\n";
    $pos += $satellite_len;
  }
  close OUT;

  print "Now run:\n";
  print "../RAMExtend -twobit $outDir/test.2bit -ranges $outDir/test.tsv\n\n";

}

