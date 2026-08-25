#!/usr/local/bin/perl
use strict;
use lib '/home/rhubley/projects/RepeatMasker';
use SearchResult;

my $direction = $ARGV[0];
my $file = $ARGV[1];

if ( $direction ne "right" && $direction ne "left" ) {
  print "Usage: calcPaths.pl <direction>  <path file>\n";
  print "  direction = right or left\n";
  exit;
}
if ( $direction eq "left" ) {
  print "Only right direction supported currently\n";
  exit;
}

my $verbose = 0;

my %seqs;
my $inSeq = 0;
my $overall_cons = "";
my %seq_cons = ();
open DAT,"<$file" or die "Can't open $file\n";
while (<DAT>) {
  #dir=right:n=2:row=0: -C--C---------Tgtgtgcggctgct best:score=-7:offset=13:cons=C
  if ( /dir=(right|left):n=(\d+):row=(\d+):\s+(\S+) best:score=([\d-]+):offset=([-\d]+):cons=(\S)/ ) {
    my ($dir, $n, $row, $path, $score, $offset, $cons) = ($1, $2, $3, $4, $5, $6, $7);
    my @calls = split(//,$path);
    next if ( $dir eq $direction );
    # Only doing right currently...
    push @{$seqs{$n}}, [ $score, $offset, $cons, @calls ];
    if ( exists $seq_cons{$n} ) {
      $seq_cons{$n} .= $cons;
    }else {
      $seq_cons{$n} = $cons;
    }
  }
  if ( /^>right-extension/ ) {
    $inSeq = 1;
    next;
  }
  if ( $inSeq ) {
    if ( /^([ACGTNX]+)\s*$/ ) {
      s/[\n\r\s]+//g;
      $overall_cons .= $_;
    }else {
      $inSeq = 0;
    }
  }
}
close DAT;



my @seq_keys = keys(%seqs);
# override
#@seq_keys = ( 37 );

foreach my $id ( @seq_keys ) {
  my $seq = $seqs{$id};
  my $seqC = $seq_cons{$id};

  my $bestRow = -1;
  my $bestScore = -10000000000;
  for (my $i = $#{$seq}; $i >= 0; $i--) {
    my $row = $seq->[$i];
    my $score = $row->[0];
    if ( $score > $bestScore ) {
      $bestScore = $score;
      $bestRow = $i;
    }
  }
  
  if ( $verbose ) { 
    print "Overall Consensus = $overall_cons\n";
    print "Sequence # $id\n";
    print "Seq Consensus = $seqC\n";
    print "bestScore = $bestScore, bestRow = $bestRow\n";
  }

## TODO: This doesn't work for left currently
##       Need to check orientation...only tested on fwd
  my $currOffset = $seq->[$bestRow]->[1]+3;
  my $align = "";
  for (my $i = $bestRow; $i >= 0; $i--) {
    my $row = $seq->[$i];
    my $call = $row->[$currOffset];
    $align = $call . $align;
    if ( $call eq "-" ) {
      $currOffset++;
    }elsif ( $call =~ /[ACGTNX]/ ) {
    }else {
      $currOffset--;
      $i++;
    }
  }
  print "----------\n" if ( $verbose );
  
  my $cons_segment;
  my $cons_pos = 0;
  for ( my $i = 0; $i < length($align); $i++ ) {
    my $align_base = substr($align, $i, 1);
    if ( $align_base =~/[acgtn]/ ) {
      $cons_segment .= "-";
    }else {
      $cons_segment .= substr($seqC, $cons_pos, 1);
      $cons_pos++;
    }
  }
  
  if ( $verbose ) {
    print "Alignment = $align\n";
    print "ConsSegmt = $cons_segment\n";
    print "----------\n";
  }
  
  $align = uc($align);
  my $consSegBaseLen = ($cons_segment =~ tr/ACGTNX//);
  my $consBaseLen = ($seqC =~ tr/ACGTNX//);
  my $alignBaseLen = ($align =~ tr/ACGTNX//);
  
  my $result = SearchResult->new(
           queryName      => "cons",     
           queryStart     => 1,
           queryEnd       => $consSegBaseLen,
           queryRemaining => ( $consBaseLen - $consSegBaseLen ),
           queryString    => $cons_segment,             
           subjString     => $align,              
           orientation    => "+",                  
           subjName  => "n=$id",
           subjStart => 1,
           subjEnd   => $alignBaseLen,
           subjRemaining => 0,
           pctDiverge    => 0,
           pctInsert     => 0, 
           pctDelete     => 0,
           score         => 0
        );                                         
  
  print "" . $result->toStringFormatted( SearchResult::AlignWithQuerySeq) . "\n";
}
  
  

