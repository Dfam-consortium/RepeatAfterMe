/*
 * bnw_extend.c
 *
 * Banded Needleman-Wunsch Glocal Alignment
 *
 * This is essentially a banded implementation of the
 * Needleman Wunsch global alignment algorithm with support
 * for custom matrix scoring and affine gap penalties.
 *
 *   Needleman, Saul B. & Wunsch, Christian D. (1970). "A general
 *   method applicable to the search for similarities in the amino
 *   acid sequence of two proteins". Journal of Molecular
 *   Biology. 48 (3): 443–53
 *
 * There is a subtle difference in this implementatation.
 * The edge case for the dynamic programming matrix supports
 * "pre-alignment" up to the 1/2 the width-1 of the bandwidth.
 * Pre-alignment allows the algorithm to consider starting
 * flanking sequence in the overall alignment. If there isn't
 * an flanking sequence or less than the bandwidth would support,
 * the algorithm falls back on the standard NW edge case based
 * on the gap penalties and the row of the matrix.
 *
 * Authors:
 *    Robert Hubley
 *    Based on concepts from RepeatScout code
 *     by Alkes Price 2005
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "minunit.h"
#include "common.h"
#include "sequence.h"
#include "report.h"
#include "score_system.h"
#include "kentsrc/common.h"
#include "kentsrc/dnautil.h"
#include "bnw_extend.h"

// TODO: Formalize max negative sentinel for row scores

// TODO: develop a routine to initialize the first row of the matrix

// TODO: return the score of the extension

// TODO: make the edge case treatment ( allowing for outside sequence influences ) optional

// TODO: code a search strategy that extends based on either two positions within the genome
//       or a position within the genome and a position within a consensus ( align/mask )
//       The search code could also calculate a complexity adjusted score based on the
//       genomic composition of matched sequence.


//
//
//
//void initialize_score(int ****score, int num_align, int bandwidth, struct *scoreParams)
//{
//}

//
// Score is a 4 dimensional array:
//
//    score[current/previous row][kmer idx][offset][sub/gap]
//
//        current/previous row: The two rows of the dynamic
//                                programming matrix that we preserve
//                                as we progress. [2]
//        kmer idx:             The index for the matrix specific
//                                to the kmer copy we are extending.
//                                This should be as big as the max
//                                number of alignments allowed. [MAXN]
//        offset:               The amount of unbalanced insertions or
//                                deletions that are tolerated (i.e. 1/2
//                                the bandwidth of the matrix calculation).
//                                This sets the number of columns calculated
//                                in each matrix row. [(2*MAXOFFSET)+1]
//        sub/gap:              The two values that make up each cell
//                                of the matrix. [2]
//
//  ****score = 2 * MAXN * (2*MAXOFFSET+1) * 2 integers
//
int ****
allocate_score(int num_align, int bandwidth)
{
  int r, x, n;
  int ****score;

  if ( num_align < 1 )
    return ((int ****)NULL);

  // Start by allocating 2 int pointers for current and previous rows
  if ((score = (int ****) malloc(2 * sizeof(int***))) == NULL)
  {
    printf("allocate_score: Error could not allocate memory for curr/prev score rows\n");
    exit(1);
  }
  for (x = 0; x < 2; x++)
  {
    // Allocate pointers for the maximum number of alignments allowed
    if ((score[x] = (int ***) malloc(num_align * sizeof(int**))) == NULL)
    {
      printf("allocate_score: Error could not allocate score rows for"
             " %d alignments.\n", num_align);
      exit(1);
    }
    for (n = 0; n < num_align; n++)
    {
      // Allocate pointers for the columns ( up to 2*bandwidth+1 ).
      if ((score[x][n] =
           (int **) malloc((2 * bandwidth + 1) * sizeof(int*))) ==
          NULL)
      {
        printf("allocate_score: Error could not allocate %d columns ( bandwidth = %d )"
               " in each score row.\n", (2 * bandwidth + 1), bandwidth);
        exit(1);
      }
      // Finally allocate 2 values per matrix cell to hold sub/gap scores
      for (r = 0; r < (2 * bandwidth + 1); r++)
      {
        if ((score[x][n][r] =
             (int *) malloc(2 * sizeof(int))) == NULL)
        {
          printf("allocate_score: Error could not allocate memory for gap/match score rows\n");
          exit(1);
        }
      }
    }
  }
  return( score );
}

// TODO: score should be refactored to be a structure containing info about it's size
void
free_score(int num_align, int bandwidth, int ****score)
{
  int r, x, n;
  for (x = 0; x < 2; x++)
  {
    for (n = 0; n < num_align; n++)
    {
      for (r = 0; r < (2 * bandwidth + 1); r++)
      {
        free(score[x][n][r]);
      }
      free(score[x][n]);
    }
    free(score[x]);
  }
  free(score);
}

//  pair_seed_extend()
//
//  Perform a BLAST-like seed/core extension on two sequences
//  found in the sequenceLibrary.  In addition to providing the
//  location of the seed, the alignable ranges of each sequence
//  are given.  This is a convenience method that calls cons_seed_extend()
//
//   TODO: Shouldn't this be expressed in terms of cores e.g
//         core start/end
//
int
pair_seed_extend( int orient,
            uint64_t seq1_start, uint64_t seq1_end, uint64_t seq1_seed,
            uint64_t seq2_start, uint64_t seq2_end, uint64_t seq2_seed,
            int ****score, int l, struct sequenceLibrary *seqLib,
            struct scoringSystem *scoreParams, int bandwidth,
            struct glocalSearchResult *result, int use_complexity_adjust, int y_drop)
{
   uint64_t subj_seed_start = seq2_seed;
   uint64_t subj_seed_end = seq2_seed + l - 1;

   return cons_seed_extend( &(seqLib->sequence[seq1_start]), seq1_end - seq1_start + 1,
                 seq1_seed - seq1_start, seq1_seed - seq1_start + l - 1,
                 seqLib, subj_seed_start, subj_seed_end, orient, score, scoreParams, bandwidth,
                 result, use_complexity_adjust, y_drop );
}


//
//  cons_seed_extend()
//
//  Perform a BLAST-like seed/core extension on a given query
//  sequence and a core found in the sequenceLibrary.
//
//     - Computes only the local alignment that must
//       contain the pre-aligned core as a subalignment.
//
//     - The core alignment forms the starting point
//       of each extension making it possible to know
//       start/end positions of the final alignment
//       without traceback.
//
//     - Unlike BLAST and like local alignment each
//       extension left/right starts from a 0 score.
//       Therefore the total score is only a sum of
//       the left/right extensions and does not
//       account for the seed itself.  This makes it
//       possible to have different seed approaches
//       which may be scored independently of the
//       alignment.
//
//     - This uses a y-drop (See lastz) feature to
//       terminate extension if the score drops by
//       more than this threshold. If y_drop is negative
//       extend until either one or both sequences
//       are exhausted and return the global high score.
//
int
cons_seed_extend(char *query_seq, int query_len, int query_seed_start, int query_seed_end,
             struct sequenceLibrary *subj_lib, uint64_t subj_seed_start,
             uint64_t subj_seed_end, int orient, int ****score,
             struct scoringSystem *scoreParams, int bandwidth,
             struct glocalSearchResult *result, int use_complexity_adjust, int y_drop )
{
  int x;
  int offset;
  int row_idx;
  uint64_t lower_seq_bound;
  uint64_t upper_seq_bound;

  if (0)
  {
    // TODO: Ensure this doesn't go outside sequence length and 0 boundaries.
    printf("queryPtr = %p subjectPtr = %p\n", query_seq, subj_lib->sequence);
    printf("cons_seed_extend orient = %d ( query len: %d  seed: %d - %d ), "
           "( subject seed %ld - %ld ) : Called\n", orient,
           query_len, query_seed_start, query_seed_end, subj_seed_start, subj_seed_end );
   if ( 0 ) {
    printf("  query    : ");
    for (x = 0; x < query_seed_start; x++)
      printf("%c", num_to_char(query_seq[x]));
    printf("[");
    for (x = query_seed_start; x <= query_seed_end; x++)
      printf("%c", num_to_char(query_seq[x]));
    printf("]");
    for (x = query_seed_end+1; x < query_len; x++)
      printf("%c", num_to_char(query_seq[x]));
    printf("\n");
    if (orient)
    {
      printf("  sbjct (-): ...");
      for (x = subj_seed_end+10; x > subj_seed_end; x--)
        printf("%c", num_to_char(compl(subj_lib->sequence[x])));
      printf("[");
      for (x = subj_seed_end; x >= subj_seed_start; x--)
        printf("%c", num_to_char(compl(subj_lib->sequence[x])));
      printf("]");
      for (x = subj_seed_start - 1; x >= 0 && x >= subj_seed_start-10; x--)
        printf("%c", num_to_char(compl(subj_lib->sequence[x])));
    }
    else
    {
      printf("  sbjct    : ...");
      for (x = subj_seed_start-10; x < subj_seed_start; x++)
        printf("%c", num_to_char(subj_lib->sequence[x]));
      printf("[");
      for (x = subj_seed_start; x <= subj_seed_end; x++)
        printf("%c", num_to_char(subj_lib->sequence[x]));
      printf("]");
      for (x = subj_seed_end + 1; x < subj_seed_end+10; x++)
        printf("%c", num_to_char(subj_lib->sequence[x]));
    }
    printf("...\n");
   }
  }

  for (offset = -bandwidth; offset <= bandwidth; offset++)
  {
    // Affine GAP penalties
    if (offset == 0)
    {
      score[1][0][offset + bandwidth][0] = 0;
      score[1][0][offset + bandwidth][1] = 0;
    }else{
      score[1][0][offset + bandwidth][1] =
        (abs(offset) * scoreParams->gapextn) + scoreParams->gapopen;
      score[1][0][offset + bandwidth][0] =
        (abs(offset) * scoreParams->gapextn) + scoreParams->gapopen;
    }
  }

  int seed_seqIdx = 0;
  while (subj_lib->boundaries[seed_seqIdx] != 0)
  {
    if (subj_lib->boundaries[seed_seqIdx] > subj_seed_start)
    {
      break;
    }
    seed_seqIdx++;
  }
  lower_seq_bound = 0;
  if ( seed_seqIdx > 0 )
    lower_seq_bound = subj_lib->boundaries[seed_seqIdx - 1];
  upper_seq_bound = subj_lib->boundaries[seed_seqIdx] - 1;

  // Seq1 is assumed to be on the forward strand and orient refers
  // seq2.
  struct coreAlignment coreAlign;
  coreAlign.seqIdx = seed_seqIdx;
  // This structure holds the coordinates of the first sequence
  // positions for extension.  Add/Subtract from the seed position
  // accordingly.
  if ( orient ) {
    coreAlign.leftSeqPos = subj_seed_end;
    coreAlign.rightSeqPos = subj_seed_start;
  }else {
    coreAlign.leftSeqPos = subj_seed_start;
    coreAlign.rightSeqPos = subj_seed_end;
  }
  coreAlign.orient = orient;
  coreAlign.leftExtendable = 1;
  coreAlign.rightExtendable = 1;

  // Left extension
  int seq_left_maxscore_matrix_idx = -1;
  int cons_left_maxscore_matrix_idx = -1;
  int left_max_score = 0;

  // Setup the bounds so that the glocal alignment
  // has standard NW edge conditions
  uint64_t lwBnd = lower_seq_bound;
  uint64_t upBnd = subj_seed_start-1;
  if ( orient ){
    lwBnd = subj_seed_end+1;
    upBnd = upper_seq_bound;
  }

  if (0)
    printf("LEFT:\n");

  for ( row_idx = 0; row_idx < query_seed_start; row_idx++ ) {
    char query_base = query_seq[query_seed_start - row_idx - 1];
    int lscore_idx = 0;
    int lscore = compute_nw_row(0, row_idx, 0, query_base,
                   &coreAlign, score,
                   lwBnd, upBnd,
                   subj_lib->sequence, &lscore_idx,
                   scoreParams, bandwidth, 10000, 0, NULL);
    if (0) {
      printf("%d         Gap: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][1]);
      printf("\n");
      printf("%d         Sub: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][0]);
      printf("\n");
    }
    if (0)
      printf("  row %d: score = %d\n", row_idx, lscore);

    // Terminate extension if score drop is > y_drop
    // (disable threshold check if y_drop is negative)
    if ( y_drop >= 0 && left_max_score - lscore > y_drop )
      break;

    if ( lscore > left_max_score ){
      left_max_score = lscore;
      seq_left_maxscore_matrix_idx = lscore_idx;
      cons_left_maxscore_matrix_idx = row_idx;
    }

  }

  //printf("left max score = %d @ considx %d\n", left_max_score, cons_left_maxscore_matrix_idx);

  //
  // initialize boundary conditions in score[1][][]
  //
  for (offset = -bandwidth; offset <= bandwidth; offset++)
  {
    // Affine GAP penalties
    if (offset == 0)
    {
      score[1][0][offset + bandwidth][0] = 0;
      score[1][0][offset + bandwidth][1] = 0;
    }else{
      score[1][0][offset + bandwidth][1] =
        (abs(offset) * scoreParams->gapextn) + scoreParams->gapopen;
      score[1][0][offset + bandwidth][0] =
        (abs(offset) * scoreParams->gapextn) + scoreParams->gapopen;
    }
  }

  // Again, do not consider the core alignment in
  // this glocal extension.
  lwBnd = subj_seed_end+1;
  upBnd = upper_seq_bound;
  if ( orient ){
    lwBnd = lower_seq_bound;
    upBnd = subj_seed_start-1;
  }

  // Right extension
  if (0)
    printf("RIGHT:\n");
  int seq_right_maxscore_matrix_idx = -1;
  int cons_right_maxscore_matrix_idx = -1;
  int right_max_score = 0;
  for ( row_idx = 0; row_idx < query_len-(query_seed_end+1); row_idx++ ) {
    char query_base = query_seq[query_seed_end + row_idx + 1 ];
    int rscore_pos = 0;
    int rscore = compute_nw_row(1, row_idx, 0, query_base,
                   &coreAlign, score,
                   lwBnd, upBnd,
                   subj_lib->sequence, &rscore_pos,
                   scoreParams, bandwidth, 10000, 0, NULL);

    if ( 0 ) {
      printf("%d         Gap: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][1]);
      printf("\n");
      printf("%d         Sub: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][0]);
      printf("\n");
    }
    if (0)
      printf("  row %d: score = %d\n", row_idx, rscore);

    // Terminate extension if score drop is > y_drop
    // (disable threshold check if y_drop is negative)
    if ( y_drop >= 0 && right_max_score - rscore > 180 )
      break;

    if ( rscore > right_max_score ){
      right_max_score = rscore;
      seq_right_maxscore_matrix_idx = rscore_pos;
      cons_right_maxscore_matrix_idx = row_idx;
    }

  }
  //printf("right max score = %d @ considx %d\n", right_max_score, cons_right_maxscore_matrix_idx);

  uint64_t subj_start = subj_seed_start - (seq_left_maxscore_matrix_idx+1);
  uint64_t subj_end = subj_seed_end + (seq_right_maxscore_matrix_idx+1);
  if ( orient ) {
    subj_end = subj_seed_end + (seq_left_maxscore_matrix_idx+1);
    subj_start = subj_seed_start - (seq_right_maxscore_matrix_idx+1);
  }

  int combined_score = left_max_score + right_max_score;
  if ( use_complexity_adjust ) {
    // TODO: is left/right always in ascending order?
    // Get sequence composition of alignment from the genome
    int comp[4] = { 0, 0, 0, 0 };
    //for (x = seq_left_maxscore_matrix_idx; x <= seq_right_maxscore_matrix_idx; x++)
    //printf("composition adjusting over subject range: %ld - %ld\n", subj_start, subj_end);
    uint64_t xx = 0;
    int base;
    for (xx = subj_start; xx <= subj_end; xx++ )
    {
      base = (int)subj_lib->sequence[xx];
      // Screen out N's (99) and translate acgt to ACGT
      if ( base <= 3 )
        comp[base&3]++;
    }
    combined_score = complexityAdjust(combined_score,
                                     comp,
                                     scoreParams->m_lambda,
                                     scoreParams->m_bg_freqs);
  }

  /*
  printf("cons_seed_extend() scores = left:%d/right:%d, cons_ext:%d/%d, seq_ext:%d/%d\n"
         "                   cons_range: %d - %d, seq_range: %ld - %ld\n",
         left_max_score, right_max_score,
         cons_left_maxscore_matrix_idx + 1, cons_right_maxscore_matrix_idx + 1,
         seq_left_maxscore_matrix_idx + 1, seq_right_maxscore_matrix_idx + 1,
         ( query_seed_start - (cons_left_maxscore_matrix_idx+1) ), (query_seed_end + (cons_right_maxscore_matrix_idx+1)),
         subj_start, subj_end
         );
  */

  if ( result != NULL )
  {
    result->score = combined_score;
    result->query_start = query_seed_start - (cons_left_maxscore_matrix_idx+1);
    result->query_end = query_seed_end + (cons_right_maxscore_matrix_idx+1);
    result->subj_start = subj_start;
    result->subj_end = subj_end;
  }

  return ( combined_score );
}



//
//         (low index)                                                   (high index)
//       seqLib->sequence[seq1_start] -----> [seed]--->[seed+l] ----> seqLib->sequence[seq1_end]
//       seqLib->sequence[seq2_start] -----> [seed]--->[seed+l] ----> seqLib->sequence[seq2_end]
//
//  Returns: result
int
global_align( int orient,
            uint64_t seq1_start, uint64_t seq1_end, uint64_t seq1_seed,
            uint64_t seq2_start, uint64_t seq2_end, uint64_t seq2_seed,
            int ****score, int l, struct sequenceLibrary *seqLib,
            struct scoringSystem *scoreParams, int bandwidth)
            //struct glocalSearchResult *result)
{
  int x;
  int offset;
  int row_idx;

  if (0)
  {
    printf("global_align orient = %d ( start %ld, end %ld, seed %ld ), "
           "( start %ld, end %ld, seed %ld ) : Called\n", orient,
           seq1_start, seq1_end, seq1_seed, seq2_start, seq2_end, seq2_seed);
    printf("seq1    :");
    for (x = seq1_start; x < seq1_seed; x++)
      printf("%c", num_to_char(seqLib->sequence[x]));
    printf("[");
    for (x = seq1_seed; x < seq1_seed+l; x++)
      printf("%c", num_to_char(seqLib->sequence[x]));
    printf("]");
    for (x = seq1_seed + l; x < seq1_end; x++)
      printf("%c", num_to_char(seqLib->sequence[x]));
    printf("\n");
    if (orient)
    {
      printf("seq2(-) :");
      for (x = seq2_end; x > seq2_seed; x--)
        printf("%c", num_to_char(compl(seqLib->sequence[x])));
      printf("[");
      for (x = seq2_seed; x >= seq2_seed-l+1; x--)
        printf("%c", num_to_char(compl(seqLib->sequence[x])));
      printf("]");
      for (x = seq2_seed-l; x >= seq2_start; x--)
        printf("%c", num_to_char(compl(seqLib->sequence[x])));
      printf("seed range: %ld to %ld\n", seq2_seed, seq2_seed-l+1);
    }
    else
    {
      printf("seq2    :");
      for (x = seq2_start; x < seq2_seed; x++)
        printf("%c", num_to_char(seqLib->sequence[x]));
      printf("[");
      for (x = seq2_seed; x < seq2_seed+l; x++)
        printf("%c", num_to_char(seqLib->sequence[x]));
      printf("]");
      for (x = seq2_seed + l; x < seq2_end; x++)
        printf("%c", num_to_char(seqLib->sequence[x]));
    }
    printf("\n");
  }

  //
  // initialize boundary conditions in score[1][][]
  //
  for (offset = -bandwidth; offset <= bandwidth; offset++)
  {
    // Affine GAP penalties
    if (offset < 0)
    {
      // orig
      score[1][0][offset + bandwidth][1] = 0;
      score[1][0][offset + bandwidth][1] +=
        ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
      score[1][0][offset + bandwidth][0] = 0;
      score[1][0][offset + bandwidth][0] +=
        ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
    }
    else
    {
      if (offset > 0)
      {
        // orig
        score[1][0][offset + bandwidth][1] = 0;
        score[1][0][offset + bandwidth][1] +=
          ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
        score[1][0][offset + bandwidth][0] = 0;
        score[1][0][offset + bandwidth][0] +=
          ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
      }
      else
      {
        score[1][0][offset + bandwidth][0] = 0;
        score[1][0][offset + bandwidth][1] = 0;
      }
    }
  }

  // Seq1 is assumed to be on the forward strand and orient refers
  // seq2.
  struct coreAlignment coreAlign;
  // TODO: All this needs to be validated
  coreAlign.seqIdx = 0;
  coreAlign.leftSeqPos = seq2_start-1;
  coreAlign.rightSeqPos = seq2_start-1;
  coreAlign.leftExtendable = 0;
  coreAlign.rightExtendable = 1;
  coreAlign.orient = orient;
  uint64_t lower_seq_bound = seq2_start;
  uint64_t upper_seq_bound = seq2_end + 1;

  int max_local_score_pos = -1;
  int max_local_score = 0;
  int max_global_score = 0;
  for ( row_idx = 0; row_idx < seq1_end-seq1_start+1; row_idx++ ) {
    char seq1_base = seqLib->sequence[seq1_start+row_idx];
    // Right extension from start to end
    int row_local_score_pos = 0;
    // word size 0
    int row_local_score = compute_nw_row(1, row_idx, 0, seq1_base,
                   &coreAlign, score,
                   lower_seq_bound, upper_seq_bound,
                   seqLib->sequence, &row_local_score_pos,
                   scoreParams, bandwidth, 10000, 0, NULL);
    if ( row_local_score > max_local_score ){
      max_local_score = row_local_score;
      max_local_score_pos = row_local_score_pos;
    }
    if (0) {
      printf("%d         Gap: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][1]);
      printf("\n");
      printf("%d         Sub: ", row_idx);
      for (offset = -bandwidth; offset <= bandwidth; offset++)
        printf(" %d", score[row_idx % 2][0][offset + bandwidth][0]);
      printf("\n");
    }
  }

  // Obtain global score from max of last row
  for (offset = -bandwidth; offset <= bandwidth; offset++)
  {
    if ( max_global_score < score[(row_idx-1)%2][0][offset + bandwidth][0] )
      max_global_score = score[(row_idx-1)%2][0][offset + bandwidth][0];
    if ( max_global_score < score[(row_idx-1)%2][0][offset + bandwidth][1] )
      max_global_score = score[(row_idx-1)%2][0][offset + bandwidth][1];
  }

  printf("global score = %d local score = %d @ %d\n", max_global_score, max_local_score, max_local_score_pos);

  // TODO: Complexity adjust score if requested

  return ( max_global_score );
}


//
// compute_nw_row
// --------------
//
// Calculate the scores for a row in a given Needleman
// Wunsch matrix from the top left to the bottom right.
//
// NW Matrix Layout
//
// l
//  m
//   e
//    r Sequence (both directions)-------->
//     +----------------------------------+
//  0  |[][][][][][]                      |
//     |[][][][][][][]                    |
//  C  |[][][][][][][][]                  |
//  o  |[][][][][][][][][]                |
//  n  |[][][][][][][][][][]              |
//  s  |[][][][][][][][][][][]            |
//  e  |  [][][][][][][][][][][]          |
//  n  |    [][][][][][][][][][][]        |
//  s  |      [][][][][][][][][][][]      |
//  u  |                                  |
//  s  |       Ins----->                  |
//  |  |      D\                          |
//  |  |      e \                         |
//  |  |      l  \ Sub                    |
//  v  |      |   v                       |
//     |      v                       [][]|[][][][][][][][][]
//     |                                []|[][][][][][][][][][]
// L-1 +----------------------------------+
//
// Edge Case:
//   Given: match/mismatch = 1/-1 and gapopen/gapext = -2/-1
//
//     *  C  G  A  C  T  T  G  A  T  C
//    -5 -4 -3  0 -3 -4 -5 -6 -7 -8 -9
//               +--------------------
//  C -3  1 -1   |
//
//   Edge Row: Diagonal center = 0, off diagonal scored as gapopen + distance
//             off-diagonal (col) * gapext.
//   Edge Col: If flanking sequence exists scored as if it's inside the matrix.
//             Otherwise scored as a gapopen + distance off diagonal (row) *
//             gapext.
//
// Truely, rows in the picture above are stored as two rows each ( gap values,
// substituion values ):
//
//      score[current/previous row][kmer idx][offset][sub/gap]
//
//   Where:
//      current/previous row : A flip flop based (0 or 1) representing the
//                               current/previous rows.  The inital value of
//                               this is row_idx=0%2 = 0.  Therefore the boundary
//                               values are stored in -1%2 = 1.
//      kmer idx             : The index to lookup the starting position in
//                               pos[]/rev[].
//      offset               : The bandwidth offset from -bandwidth to bandwidth
//      sub/gap              : A boolean value, selecting which vector/array is
//                               being used.  0 = substition scores, and 1 = gap
//                               scores.
// PARAMETERS:
//   direction               : Direction for the extension 0=left, 1=right
//   row_idx                 : The matrix row index
//   n                       : The index of the kmer/left_edge being extended
//   cons_base               : Consensus base
//   coreAlign               : A coreAlignment structure containing the particular
//                               aligned region to extend.
//   score                   : The SW matrix (current and previous banded rows
//                               only) for all kmers/left_edges.
//   lower_seq_bound         : The lower bound on the sequence containing kmer
//                               /left_edge "n".
//   upper_seq_bound         : The upper bound on the sequence containing kmer
//                               /left_edge "n".
//   sequence                : The sequence array.
//   max_score_sequence_idx  : The sequence index for the max scoring cell in the
//                               calculated row. The value stored at this pointer
//                               is updated in this routine.
//   scoreParams             : The scoring system to use ( Matrix + Gap Penalties )
//   L                       : The max extension allowed in either direction
//   bandwidth               : Half the bandwidth. The number of cells per row is
//                               2*bandwidth+1.
//   VERBOSE                 : Verbosity level
//   pathString              : A string to store the alignment path for this row
//                               This should be capable of storing up to bandwidth*2
//                               symbols or a NULL pointer if this data is not
//                               required.
//
// RETURNS: Best row score for a given sequence (n), consensus base (cons_base),
//          and row index (row_index).
// MODIFIES: score[row_idx%2][][][], and max_score_sequence_idx
//
// TODO: This routine doesn't have any way of knowing if 'n' is too large
//       for the score datastructure.  Unfortunately it will cause a segfault
//       if it goes out-of-bounds on the score DS.  Perhaps the score DS
//       should encapsulate metadata on it's bounds.
int
compute_nw_row(int direction, int row_idx, int n, char cons_base,
               struct coreAlignment *coreAlign, int ****score,
               uint64_t lower_seq_bound, uint64_t upper_seq_bound,
               char *sequence, int *max_score_sequence_idx,
               struct scoringSystem *scoreParams, int bandwidth, int L, int VERBOSE,
               char *pathString)
{
  int oldoffset, tempscore, ins_score, del_score, gap_score, sub_score, cell_score, subvalue;
  int offset;
  int best_row_score = -1000000000;
  int seq_offset;
  int out_of_bounds;
  uint64_t seq_idx;
  uint64_t start_seq_pos;

  // Out of bounds score sentinel
  int OOBSENTINEL = -987654321;

  // Resolve this here for efficiency
  int curr_row_ff = row_idx % 2;
  int prev_row_ff = !curr_row_ff;

  // Identify the starting sequence position
  // as the first base outside of the core. This
  // depends on the direction we are going and the
  // orientation of the core.  See table below for
  // a description of this.
  if ( direction ) {
    if ( coreAlign->orient )
      start_seq_pos = coreAlign->rightSeqPos - 1;
    else
      start_seq_pos = coreAlign->rightSeqPos + 1;
  }else {
    if ( coreAlign->orient )
      start_seq_pos = coreAlign->leftSeqPos + 1;
    else
      start_seq_pos = coreAlign->leftSeqPos - 1;
  }

  if (0) {
    if ( row_idx == 0 && n == 0 ) {
      if ( coreAlign->orient )
        printf("n==%d/row==%d considering sequence base(compl): %c\n", n, row_idx, num_to_char(compl((int) sequence[start_seq_pos])));
      else
        printf("n==%d/row==%d considering sequence base: %c\n", n, row_idx, num_to_char((int) sequence[start_seq_pos]));
    }
  }

  // Must go from left to right in this SW matrix
  for (offset = -bandwidth; offset <= bandwidth; offset++)
  {
    ins_score = -1000000000;
    del_score = -1000000000;
    sub_score = -1000000000;

    // row_idx is the row (consensus) index from 0-L
    // This is the index from the top-left corner to the bottom-right of the
    // matrix.
    // It is also the offset from the starting sequence position in the
    // top-left
    // corner increasing or decreasing based on the direction and the strand
    // of the sequence.
    //
    //   ie.:
    //        Direction  Orientation  Seq_progression  Complement?
    //         1=Right      0=Fwd          Increasing       No
    //         1=Right      1=Rev          Decreasing       Yes
    //         0=Left       0=Fwd          Decreasing       No
    //         0=Left       1=Rev          Increasing       Yes
    //
    seq_idx = 0;
    out_of_bounds = 0;

    if ( direction == coreAlign->orient )
    {
      // Right/Rev or Left/Fwd
      // Decreasing sequence direction
      seq_offset = -offset-row_idx;
    }
    else
    {
      // Right/Fwd or Left/Rev
      // Increasing sequence direction
      seq_offset = offset+row_idx;
    }

    seq_idx = start_seq_pos + seq_offset;


    // Overflow (should be called underflow but the world is a strange
    // place ) of an unsigned int turns very very large.  Must be careful
    // with the indices.
    if ( start_seq_pos < bandwidth &&
         seq_offset < 0 && abs(seq_offset) > start_seq_pos )
    {
      out_of_bounds = 1;
      if (VERBOSE > 12 && VERBOSE < 100)
        printf
          ("Short circuiting right extension n=%d ( start_seq_pos = %ld, "
           "seq_offset = %d ) because of lower sequence bound %ld [overflow avoided]\n", n,
            start_seq_pos, seq_offset, lower_seq_bound);
    }else if ( seq_idx > upper_seq_bound )
    {
      out_of_bounds = 1;
      if (VERBOSE > 12 && VERBOSE < 100)
        printf
          ("Short circuiting right extension n=%d ( start_seq_pos = %ld, "
           "seq_offset = %d ) because of upper sequence bound %ld\n", n,
            start_seq_pos, seq_offset, upper_seq_bound);
    }else if ( seq_idx < lower_seq_bound )
    {
      out_of_bounds = 1;
      if (VERBOSE > 12 && VERBOSE < 100)
        printf
          ("Short circuiting right extension n=%d ( start_seq_pos = %ld, "
           "seq_offset = %d ) because of lower sequence bound %ld\n", n,
           start_seq_pos, seq_offset, lower_seq_bound);
    }

    if (!out_of_bounds)
    {
      /*
       * Deletion Case:
       * This means that cons_base aligns to "-":
       *
       * [0][1][2]
       * | Deletion moves down in SW matrix.
       * | In this representation...the
       * | previous row is in score[(row_idx-1)%2]...
       * | so the index would be offset by +1.
       * | Thus ( oldoffset + 1 is directly above
       * | our cell ).
       * v
       * [0][1][2]
       *
       * In this new recursion we account for the score coming
       * from either a substition state or from a gap state.
       * The difference being the tax of a gap_open or a gap_ext
       * penalty.
       */
      // Can a deletion occur?
      if (offset < bandwidth)
      {
        oldoffset = offset + 1;
        // Which is better a gap open? Substitution state = score[][][][0].
        tempscore =
          score[prev_row_ff][n][oldoffset + bandwidth][0] + scoreParams->gapopen +
          scoreParams->gapextn;
        if (tempscore > del_score)
          del_score = tempscore;
        // Or a gap extension? Gap state = score[][][][1].
        tempscore = score[prev_row_ff][n][oldoffset + bandwidth][1] + scoreParams->gapextn;
        if (tempscore > del_score)
          del_score = tempscore;
      }

      /*
       * Substitution Case:
       *
       * This means that cons_base aligns to sequence[seq_idx].
       *
       * [0][1][2]
       *  | Substitution moves diagonally in SW matrix.
       *  | In this representation...the
       *  | previous row is in score[(row_idx-1)%2]...
       *   \ so the index would be the same offset.
       *    \ Thus ( oldoffset = offset is the top
       *     | left diagonal to our cell ).
       *     v
       *    [0][1][2]
       *
       * Again we need to determine which does better Gap state
       * to Substitution state or Substitution state to Substitution
       * to Substitution state. Blah blah.
       */
      oldoffset = offset;
      // Negative offsets extend into the unaligned sequence
      if (coreAlign->orient)
      {
        // Reverse Complement Alignment
        // Moving right in the sequence from the starting position
        // and complementing the base as we go.
        subvalue = scoreParams->matrix[(int) cons_base][(int) compl(sequence[seq_idx])];
      }
      else
      {
        // Forward Alignment
        // Moving left in the sequence from the starting position.
//
//if ( cons_base < 0 || cons_base > 3 ){
//  printf("ERROR: cons_base = %d\n", cons_base); fflush(stdout);
//}
//if ( sequence[seq_idx] < 0  || (sequence[seq_idx] > 7 && sequence[seq_idx] != 99 ) ) { 
//  printf("ERROR: sequence_base = %d @ %lu\n", sequence[seq_idx], seq_idx); fflush(stdout);
//}
        subvalue = scoreParams->matrix[(int) cons_base][(int) sequence[seq_idx]];
      }

      // Which is better a Substitution state to Substitution state?
      tempscore = score[prev_row_ff][n][oldoffset + bandwidth][0] + subvalue;
      if (tempscore > sub_score)
        sub_score = tempscore;
      // or a Gap state to Substitution state?
      tempscore = score[prev_row_ff][n][oldoffset + bandwidth][1] + subvalue;
      if (tempscore > sub_score)
        sub_score = tempscore;

      /*
       * Insertion Case
       *
       * Insertion moves right in the SW matrix.
       * In this representation...the
       * current row is in score[row_idx%2] and the previous
       * row in score[(row_idx-1)%2]:
       * so the index would be offset + 1.
       * Thus ( oldoffset = offset + 1 is the right
       * cell ).
       *
       * [0][1][2]
       * [0]-->[1][2]
       */
      if (offset > -bandwidth)
      {
        oldoffset = offset - 1;
        // Which is better a Substitution state to Gap state?
        tempscore =
          score[curr_row_ff][n][oldoffset + bandwidth][0] + scoreParams->gapopen +
          scoreParams->gapextn;
        if (tempscore > ins_score)
          ins_score = tempscore;
        // or a Gap state to Gap state?
        tempscore = score[curr_row_ff][n][oldoffset + bandwidth][1] + scoreParams->gapextn;
        if (tempscore > ins_score)
          ins_score = tempscore;
      }
    }
    else
    {
      // out_of_bounds == 1 ---> OUT OF SEQUENCE BOUNDS
      if ( offset < 0 && row_idx < bandwidth ) {
        // Matrix edge case
        sub_score = scoreParams->gapopen + (( row_idx + 1 )*scoreParams->gapextn);
        ins_score = sub_score;
        del_score = sub_score;
      }else {
        // Outer edge doesn't impact matrix calc but should be
        // filled in with something really low to avoid picking
        // it in extreme cases.
        sub_score = OOBSENTINEL;
        ins_score = OOBSENTINEL;
        del_score = OOBSENTINEL;
      }
    }
    // New Substitution state score
    score[curr_row_ff][n][offset + bandwidth][0] = sub_score;

    if ( ins_score > del_score )
      gap_score = ins_score;
    else
      gap_score = del_score;

    // New Gap state score
    score[curr_row_ff][n][offset + bandwidth][1] = gap_score;

    if ( gap_score > sub_score )
      cell_score = gap_score;
    else
      cell_score = sub_score;

    if (cell_score > best_row_score)
    {
      *max_score_sequence_idx = row_idx + offset;
      best_row_score = cell_score;
    }

    // if ( direction && VERBOSE == 100 )
    if ( pathString != NULL )
    {
      char base = 'X';
      if ( seq_idx >= lower_seq_bound && seq_idx <= upper_seq_bound )
        base = num_to_char(sequence[seq_idx]);
      // Print the highest scoring cell state
      //   E.g "A","C","G","T","-" or "a","c","g","t".
      //   Where uppercase bases are substitutions, "-" are deletions
      //   and lowercase bases are insertions.
// IMPORTANT!
//   TODO: Take into account orientation if we are going to keep this debuging code
      if ( sub_score == cell_score )
          pathString[offset + bandwidth] = base;
      else if ( del_score == cell_score )
          pathString[offset + bandwidth] = '-';
      else
          pathString[offset + bandwidth] = tolower(base);
    }
  } // for(offset=..

  return best_row_score;
}


//
// Fisher-Yates Shuffle of the coreAlignment
// structure.
// NOTE: This supports a linked list containing some
//       amount of blank records.  A blank record is
//       a record with seqIdx=0, leftSeqPos=0,
//       rightSeqPos=0, leftExtendable=0 and
//       rightExtendable=0
//
struct coreAlignment *
fisher_yates_shuffle_calign(struct coreAlignment *coreAlign)
{
  uint64_t i, j;

  // coreAlignment is a linked-list structure.
  // The first step is to determine the size of the structure.
  struct coreAlignment *tCAlignPtr = coreAlign;
  int coreAlignSize = 0;
  struct coreAlignment *blankEntries = NULL;
  while ( tCAlignPtr != NULL ){
    if ( tCAlignPtr->seqIdx == 0 && tCAlignPtr->leftSeqPos == 0 &&
         tCAlignPtr->rightSeqPos == 0 && tCAlignPtr->leftExtendable == 0 &&
         tCAlignPtr->rightExtendable == 0 ){
      blankEntries = tCAlignPtr;
      break;
    }
    coreAlignSize++;
    tCAlignPtr = tCAlignPtr->next;
  }

  // Nothing to do unless there are more than 1 entries
  if ( coreAlignSize < 2 )
    return coreAlign;

  // The second step is to generate a pointer array.
  struct coreAlignment **llistPtrs = (struct coreAlignment **)malloc(
       (coreAlignSize+1)*sizeof(struct coreAlignment *));
  tCAlignPtr = coreAlign;
  for ( i = 0; i < coreAlignSize; i++ ){
    llistPtrs[i] = tCAlignPtr;
    tCAlignPtr = tCAlignPtr->next;
  }
  llistPtrs[coreAlignSize] = NULL;

  // Finally use the Fisher Yates Shuffle on the pointer array
  struct coreAlignment *tmpPtr;
  for (i = coreAlignSize - 1; i > 0; i--)
  {
    j = rand_int(i + 1);
    // Now swap i/j where j < i
    tmpPtr = llistPtrs[j];
    llistPtrs[j] = llistPtrs[i];
    llistPtrs[i] = tmpPtr;
  }

  // Now rethread the linked list
  tCAlignPtr = llistPtrs[0];
  for ( i = 1; i < coreAlignSize; i++ ){
    tCAlignPtr->next = llistPtrs[i];
    tCAlignPtr = tCAlignPtr->next;
  }
  if ( blankEntries != NULL )
    tCAlignPtr->next = blankEntries;
  else
    tCAlignPtr->next = NULL;

  tCAlignPtr = llistPtrs[0];
  free(llistPtrs);

  return tCAlignPtr;

}

//
// Phil Green's Complexity Score Adjustment
//
// Given the composition of the target sequence, the
// raw score, and characteristics of the scoring matrix
// return the complexity adjusted score.
//
// E.g.:
// int comp[4] = { 0, 0, 0, 0 };
//  for (y = bestsequencew; y <= bestsequencey; y++)
//    comp[(int) sequence[y]]++;
// TODO: store the alphabet size somewhere
int
complexityAdjust(int score, int *comp, double lambda, double *matrix_freqs)
{
  double t_factor = 0;
  double t_sum = 0;
  double t_counts = 0;
  int adj_score = 0;
  int i;

  for (i = 0; i < 4; i++)
  {
    if (comp[i])
    {
      t_factor += comp[i] * log(comp[i]);
      t_sum += comp[i] * log(matrix_freqs[i]);
      t_counts += comp[i];
    }
  }

  if (t_counts)
    t_factor -= t_counts * log((double) t_counts);
  t_sum -= t_factor;

  // To more closely match the method employed in cross_match
  // we use their rounding scheme rather than sprintf()'s.
  adj_score = score + t_sum / lambda + .999;

  if (adj_score < 0)
    adj_score = 0;

  return (adj_score);
}


MU_TEST(test_global_align) {
   char alu_seqs[] ="ATATGGATAGCTAGCGTGCGACGTACGGCGATTGGTATATGAGCGATATC"
                    //ALU Portion
                   //xxxxxxxxxxx (50-60, 11bp)
                    "TAAAAATAAAAATAGGCTTGGGCGCGGTGGCTCACGCCTGTAATCCCAGC"
                    "ACTTAGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGAC"
                    "CAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATT"
                    "AGCCGGGCATGGTGGCACGCGCCTGTAATCCCAGCTACTCGGGAGGCTGA"
                    "GGCAGGAGATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGAT"
                    "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA"
                    "AAAAAAAAAAAAAAAAAA"
                    // End ALU Portion
                    "GTAGTAGCATGCGAGCGTTAGCGATATATTAA"
                    "ATATGAGATGCGGGCTAATGCATATATTTTATGCGCGAGCACAAACGATT"
                    "GATTAGCGCGCATCACGATTAGGAGGATATGAATTATGCGGCGAGAAGAT"
                    //ALU Portion
                   //xxxxxxxxxxx
                    "TAAAAATAAAAATAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA"
                    "CTTTGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGACC"
                    "AGCCTGGCCAACATGTGAAACCCCGTGCTCTACTAAAAATACAAAAATTA"
                    "GCCGGGCATGGTGGCACGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAG"
                    "GCAGGAGTAATCGCTTGAACCCGGGAGGCGAGGTTGCAGTGAGCCGAGAT"
                    "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA"
                    "AAAAAAATAAAAAAAAA"
                    // End ALU Portion
                    "GATATCTACTATCATCTATAAATCTATTATTTA"
                    "ACGATTGGGCGGCGCTTTATTGCAGTAGTATCTAGCTCTCTCAGGCAAAC";

  int idx = 0;
  while ( alu_seqs[idx] != '\0' ){
    alu_seqs[idx] = char_to_num(alu_seqs[idx]);
    idx++;
  }

  struct sequenceLibrary seqLib;
  seqLib.sequence = alu_seqs;
  seqLib.boundaries = (uint64_t *) malloc(2 * sizeof(uint64_t));
  seqLib.boundaries[0] = 900;
  seqLib.boundaries[1] = 0;
  char *predef_id[] = { "seq1" };
  seqLib.identifiers = predef_id;
  seqLib.count = 1;
  seqLib.length = 900;

  int bandwidth = 5;
  int ****scores = allocate_score(1, bandwidth);
  struct scoringSystem *scoreParams = getMatrix("20p43g");
  int orient = 0;
  uint64_t seq1_start = 50;
  uint64_t seq1_end = 367;
  uint64_t seq1_seed = 60;
  uint64_t seq2_start = 500;
  uint64_t seq2_end = 815;
  uint64_t seq2_seed = 510;
  int l = 8; // TODO: is this needed?

  // seq1_seed and seq2_seed are not really necessary
  int align_score =  global_align(orient, seq1_start, seq1_end, seq1_seed,
             seq2_start, seq2_end, seq2_seed,
             scores, l, &seqLib, scoreParams, bandwidth);

  free(seqLib.boundaries);
  free(scores);
  freeScoringSystem(scoreParams);
  //global score = 2746 local score = 2758 @ 316 --- See spreadsheet
  mu_check( align_score == 2732 );
}



MU_TEST(test_cons_seed_extend)
{
  //                           1         2         3         4        4
  //                 01234567890123456789012345678901234567890123456789
  char cons_seq[] = "TAAAAATAAAAATAGGCTTGGGCGCGGTGGCTCACGCCTGTAATCCCAGC"
  //                 xxxxxxxxxxx ( 50-60, 11bp )
                    "ACTTAGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGAC"
                    "GGCAGGAGATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGAT"
                    "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA"
                    "AAAAAAAAAAAAAAATAA";
  int cons_start = 50;
  int cons_end = 60;
  int cons_seq_len = 218;

  char alu_seqs[] ="ATATGGATAGCTAGCGTGCGACGTACGGCGATTGGTATATGAGCGATATC"
                 //ALU Portion
                   "GGAGGAAATATTTATAGTTGGGCGCGGTGGCTCACGCCTGTAATCCCAGC"
                 // xxxxxxxxxxx (100-110)
                   "ACTTAGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGAC"
                   "GGCAGGAGATCCCTTGAACCCGGCGGAGGTTGCAGTGCATTAGCCGAGAT"
                   "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA"
                   "AAAAATTTAAAAAGAAAA"
                 // End ALU Portion
                                     "GTAGTAGCATGCGAGCGTTAGCGATATATTAA"
                   "ATATGAGATGCGGGCTAATGCATATATTTTATGCGCGAGCACAAACGATT"
                   "GATTAGCGCGCATCACGATTAGGAGGATATGAATTATGCGGCGAGAAGAT"
                 //ALU Portion
                 //                                                  x
                   "TAAAAATAAAAATAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA"
                 // xxxxxxxxxx (449-459)
                   "CTTTGGGAGGCCGAGGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGACC"
                   "AGCCTGGCCAACATGTGAAACCCCGTGCTCTACTAAAAATACAAAAATTA"
                   "GCCGGGCATGGTGGCACGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAG"
                   "GCAGGAGTAATCGCTTGAACCCGGGAGGCGAGGTTGCAGTGAGCCGAGAT"
                   "CGCGCCACTGCACTCCAGCCTGGGCAACAGAGCGAGACTCCGTCTCAAAA"
                   "AAAAAAATAAAAAAAAA"
                 // End ALU Portion
                                    "GATATCTACTATCATCTATAAATCTATTATTTA"
                   "ACGATTGGGCGGCGCTTTATTGCAGTAGTATCTAGCTCTCTCAGGCAAAC"
                 // Reversed Alu Portion
                   "TTATTTTATTTTTTTTTTTTGAGACGCAGTCTCGCTCTGTTGCCCAGGCT"
                   "GGAGCAGTGGCGCGATCTCGGCACACTGCAACCTAGCCGCCTCCCGGGTT"
                   "CAAGCGATCTCCTGCCGTCTCGAACTCCTGACCTCTAGTGATCCGCGCGC"
                 //       xxxxxxxxxxx (956-966)
                   "AATCGGCCTCCCTAAGTGCTGGGATTTCAGGCGTGAGCCACCGCGCCCAA"
                   "GTATTATTATTTTTA"
                 // End ALU Portion
                                  "AGTAGAGGTGCGGATATGCGGAGTAGCGGATATAG";

  uint64_t seq_seed_start1 = 100;
  uint64_t seq_seed_end1 = 110;
  uint64_t seq_seed_start2 = 449;
  uint64_t seq_seed_end2 = 459;
  uint64_t seq_seed_start3 = 956;
  uint64_t seq_seed_end3 = 966;

  int idx = 0;
  while ( alu_seqs[idx] != '\0' ){
    alu_seqs[idx] = char_to_num(alu_seqs[idx]);
    idx++;
  }
  idx = 0;
  while ( cons_seq[idx] != '\0' ){
    cons_seq[idx] = char_to_num(cons_seq[idx]);
    idx++;
  }

  struct sequenceLibrary seqLib;
  seqLib.sequence = alu_seqs;
  seqLib.boundaries = (uint64_t *) malloc(2 * sizeof(uint64_t));
  seqLib.boundaries[0] = 1050;
  seqLib.boundaries[1] = 0;
  char *predef_id[] = { "seq1" };
  seqLib.identifiers = predef_id;
  seqLib.count = 1;
  seqLib.length = 1050;

  int bandwidth = 5;
  int ****scores = allocate_score(1, bandwidth);
  struct scoringSystem *scoreParams = getMatrix("20p43g");
  int orient = 0;

  struct glocalSearchResult gSearchResult;

  int sc = cons_seed_extend(cons_seq, cons_seq_len, cons_start, cons_end,
                        &seqLib, seq_seed_start1, seq_seed_end1,
                        orient, scores, scoreParams, bandwidth,
                        &gSearchResult, 0, -1 );
  // Validated DP matrices ( zero-based coordinates ):
  //   Consensus Start/End: 17/214,
  //   Subject Start/End:   67/267
  //   Score Left/Right/Total: 319/1254/1573
  //
  //   left max score = 319 @ considx 32
  //   right max score = 1254 @ considx 153
  //
  //
  //if ( sc == 1573 && gSearchResult.subj_start == 67 && gSearchResult.subj_end == 267 &&
  //                   gSearchResult.query_start == 17 && gSearchResult.query_end == 214 )
  //{
  //  printf("cons_seed_extend alignment passed!\n");
  //}else {
  //  printf("cons_seed_extend alignment failed! :  %d : %ld - %ld\n", sc, gSearchResult.subj_start, gSearchResult.subj_end );
  //}
  mu_check( sc == 1573 && gSearchResult.subj_start == 67 && gSearchResult.subj_end == 267 &&
                     gSearchResult.query_start == 17 && gSearchResult.query_end == 214 );

  sc = cons_seed_extend(cons_seq, cons_seq_len, cons_start, cons_end,
                        &seqLib, seq_seed_start2, seq_seed_end2,
                        orient, scores, scoreParams, bandwidth,
                        &gSearchResult, 1, -1 );
  // Validated matrices ( zero-based coordinates ):
  //   Consensus Start/End: 0/99,
  //   Subject Start/End:   400/498
  //   UNADJUSTED SCORES: Score Left/Right/Total: 416/376/792
  //   TODO: NOTE: 743 has not been validated!
  mu_assert_int_eq( 743, sc);
  mu_assert_int_eq( 400, gSearchResult.subj_start );
  mu_assert_int_eq( 498, gSearchResult.subj_end );
  mu_assert_int_eq( 0, gSearchResult.query_start );
  mu_assert_int_eq( 99, gSearchResult.query_end );

  orient = 1;
  sc = cons_seed_extend(cons_seq, cons_seq_len, cons_start, cons_end,
                        &seqLib, seq_seed_start3, seq_seed_end3,
                        orient, scores, scoreParams, bandwidth,
                        &gSearchResult, 1, -1 );
  // Validated matrices ( zero-based coordinates ):
  //   Consensus Start/End: 0/212,
  //   Subject Start/End:   803/1014
  //   UNADJUSTED SCORES: Score Left/Right/Total: 365/1180/1486
  //   TODO: NOTE: 1462 has not been validated!
  mu_assert_int_eq( 1462, sc);
  mu_assert_int_eq( 803, gSearchResult.subj_start );
  mu_assert_int_eq( 1014, gSearchResult.subj_end );
  mu_assert_int_eq( 0, gSearchResult.query_start );
  mu_assert_int_eq( 212, gSearchResult.query_end );

  free(seqLib.boundaries);
  free(scores);
  freeScoringSystem(scoreParams);
}


MU_TEST(test_compute_nw_row)
{
  int x, n, r;
  int offset;
  int conspos;
  int bestscore;
  struct sequenceLibrary l_seqLib;
  struct sequenceLibrary r_seqLib;
  int max_score_seq_idx = 0;
  int N = 0;
  int MAXOFFSET = 5;
  int VERBOSE = 10;
  int L = 10000;
  int ****l_score; // 2 * MAXN * (2*MAXOFFSET+1) * 2
  int ****r_score; // 2 * MAXN * (2*MAXOFFSET+1) * 2

  // Structure to hold the core alignment ranges for which we
  // are going to extend.
  struct coreAlignment *coreAlign = NULL;

  //
  // Structure to hold the sequence library
  struct sequenceLibrary *seqLib = NULL;
  dnaUtilOpen();
  seqLib = loadSequenceSubset("test/extension-test2.2bit", "test/extension-test2.tsv",
                              &coreAlign, &N);
  printf("Read in %d ranges, and %ld bp of sequence\n", N, seqLib->length);

  // Structure to hold the matrix
  struct scoringSystem *scoreParams = NULL;
  scoreParams = getMatrix("20p43g");

  // score[2] : two SW rows per kmer/anchor
  l_score = (int ****) malloc(2 * sizeof(*l_score));
  r_score = (int ****) malloc(2 * sizeof(*r_score));
  for (x = 0; x < 2; x++)
  {
    // score[2][5] : up to 5 kmers/anchors
    l_score[x] = (int ***) malloc(5 * sizeof(*l_score[x]));
    r_score[x] = (int ***) malloc(5 * sizeof(*r_score[x]));
    for (n = 0; n < 5; n++)
    {
      // score[2][5][2*MAXOFFSET+1] : Row length = 2*bandwidth+1
      l_score[x][n] =
        (int **) malloc((2 * MAXOFFSET + 1) * sizeof(*l_score[x][n]));
      r_score[x][n] =
        (int **) malloc((2 * MAXOFFSET + 1) * sizeof(*r_score[x][n]));
      for (r = 0; r < (2 * MAXOFFSET + 1); r++)
      {
        // score[2][5][2*MAXOFFSET+1][2] : Gap/Sub row per SW row
        l_score[x][n][r] = (int *) malloc(2 * sizeof(*l_score[x][n][r]));
        r_score[x][n][r] = (int *) malloc(2 * sizeof(*r_score[x][n][r]));
      }
    }
  }

  // Boundary conditions are stored in score[1][][][]
  // as first row ( w=0 ), and subsequent rows are stored
  // in: score[w%2][][][]
  for (n = 0; n < 3; n++)
  {
    for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    {
      if (offset < 0)
      {
        l_score[1][n][offset + MAXOFFSET][1] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        l_score[1][n][offset + MAXOFFSET][0] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        r_score[1][n][offset + MAXOFFSET][1] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        r_score[1][n][offset + MAXOFFSET][0] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
      }
      else
      {
        if (offset > 0)
        {
          l_score[1][n][offset + MAXOFFSET][1] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          l_score[1][n][offset + MAXOFFSET][0] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          r_score[1][n][offset + MAXOFFSET][1] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          r_score[1][n][offset + MAXOFFSET][0] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
        }
        else
        {
          l_score[1][n][offset + MAXOFFSET][0] = 0;
          l_score[1][n][offset + MAXOFFSET][1] = 0;
          r_score[1][n][offset + MAXOFFSET][0] = 0;
          r_score[1][n][offset + MAXOFFSET][1] = 0;
        }
      }
    }
  }

  l_seqLib.boundaries = (uint64_t *) malloc(4 * sizeof(uint64_t));
  l_seqLib.boundaries[0] = 16;
  l_seqLib.boundaries[1] = 32;
  l_seqLib.boundaries[2] = 48;
  l_seqLib.boundaries[3] = 0;
  uint64_t lower_bound[] = { 0, 16, 32 };
  uint64_t upper_bound[] = { 16, 32, 48 };
  char *predef_id[] = { "seq1", "seq2", "seq3" };
  l_seqLib.identifiers = predef_id;
  l_seqLib.count = 3;
  l_seqLib.length = 48;
  l_seqLib.offsets = NULL;
  /*
   * n=0  seq1:         9 :   11:+  NACCTGAATC [           T            ] AGGAC*
   * n=0  seq2:        25 :   27:+  GACGTGAATC [           T            ] AGGAC*
   * n=0  seq3:        41 :   43:+  TTGATGAATG [           T            ] AGGAC*
   *
   * 0-3 = ACGT and N = 99
   *
   * LEFT Extension Example1
   *
   *                      row=0  *  *  *  *  *  *  *  *  *  *  *
   *                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
   *                 G  C  A  A  T  G  C  T  G  A  A  T  C  t  a  g */
  char p_left[] = {  2, 1, 0, 0, 3, 2, 1, 3, 2, 0, 0, 3, 1, 3, 0, 2,
    /*              16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     *               c  t  a  G  A  T  T  C  A  G  C  A  T  T  G  N*/
                     1, 3, 0, 2, 0, 3, 3, 1, 0, 2, 1, 0, 3, 3, 2, 99,
    /*              32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
     *               G  C  T  G  A  C  C  T  C  A  A  T  G  t  a  g */
                     2, 1, 3, 2, 0, 1, 1, 3, 1, 0, 0, 3, 2, 3, 0, 2
  };
  //                   C  T  A  A  C  T  C  C  A  G  T  C  G  T  A  A  C  G  G
  char cons_base[] = { 1, 3, 0, 0, 1, 3, 1, 1, 0, 2, 3, 1, 2, 3, 0, 0, 1, 2, 2 };
  l_seqLib.sequence = p_left;

  // NextPtr, seqID, leftSeqPos, rightSeqPos,
  //    leftExtendable, rightExtendable, lowerSeqBound, upperSeqBound, lowerSeqBoundFlag,
  //    upperSeqBoundFlag, llen, rlen, score, orient
  struct coreAlignment l_coreAlign_0 = { NULL, 0, 12, 14, 1, 1, 0, 15, 0, 0, 0, 0, 0, 0 };
  struct coreAlignment l_coreAlign_1 = { NULL, 1, 16, 18, 1, 1, 16, 31, 0, 0, 0, 0, 0, 1 };
  struct coreAlignment l_coreAlign_2 = { NULL, 2, 44, 46, 1, 1, 32, 47, 0, 0, 0, 0, 0, 0 };
  struct coreAlignment *l_cores[] = { &l_coreAlign_0, &l_coreAlign_1, &l_coreAlign_2 };

  /* RIGHT Extension Example1-reciprocal
   *                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
   *                  G  a  t  c  T  A  A  G  T  C  G  T  A  A  C  G*/
  char p_right[] = {  2, 0, 3, 1, 3, 0, 0, 2, 3, 1, 2, 3, 0, 0, 1, 2,
    /*               16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     *                N  G  T  T  A  C  G  A  C  T  T  A  G  a  t  c */
                     99, 2, 3, 3, 0, 1, 2, 0, 1, 3, 3, 0, 2, 0, 3, 1,
    /*               32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
     *                G  a  t  g  T  A  A  C  T  C  C  A  G  T  C  G */
                      2, 0, 3, 2, 3, 0, 0, 1, 3, 1, 1, 0, 2, 3, 1, 2
  };
  struct coreAlignment r_coreAlign_0 = { NULL, 0, 1, 3, 1, 1, 0, 15, 0, 0, 0, 0, 0, 0 };
  struct coreAlignment r_coreAlign_1 = { NULL, 1, 29, 31, 1, 1, 16, 31, 0, 0, 0, 0, 0, 1 };
  struct coreAlignment r_coreAlign_2 = { NULL, 2, 33, 35, 1, 1, 32, 47, 0, 0, 0, 0, 0, 0 };
  struct coreAlignment *r_cores[] = { &r_coreAlign_0, &r_coreAlign_1, &r_coreAlign_2 };

  r_seqLib.boundaries = (uint64_t *) malloc(4 * sizeof(uint64_t));
  r_seqLib.boundaries[0] = 16;
  r_seqLib.boundaries[1] = 32;
  r_seqLib.boundaries[2] = 48;
  r_seqLib.boundaries[3] = 0;
  r_seqLib.identifiers = predef_id;
  r_seqLib.sequence = p_right;
  r_seqLib.count = 3;
  r_seqLib.length = 48;
  r_seqLib.offsets = NULL;

  if ( VERBOSE > 10 ) {
  printf("Boundary Conditions\n");
  printf("  Gap: ");
  for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    printf(" %d", l_score[1][0][offset + MAXOFFSET][1]);
  printf("\n");
  printf("  Sub: ");
  for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    printf(" %d", l_score[1][0][offset + MAXOFFSET][0]);
  printf("\n\n");
  }

  // Consensus row_idx=6
  int final_right_gap[] = { -35, -21, -26, -17, -22, -7, -12, -17, -22, -27, -32 };
  int final_right_sub[] = { -31, -45, -31, -25, 26, -32, -25, -37, -47, -27, -57 };
  //int final_right_gap[] = { -35, -29, -15, -1, -6, -11, -2, -7, -12, -17, -22 };
  //int final_right_sub[] = { -33, -1, -41, -27, -11, 31, -8, -6, -7, -32, -43 };
  for ( n = 0; n < 3; n++ )
  {
    printf("n = %d\n", n);
    printf("Cores Left/Right reciprocal:\n");
    printCoreEdges(l_cores[n], &l_seqLib, 0, 0);
    printCoreEdges(r_cores[n], &r_seqLib, 0, 0);
    for (conspos = 0; conspos < 19; conspos++)
    {
      if ( VERBOSE > 10 )
        printf("row: %d,  cons_base=%d:\n", conspos, cons_base[conspos]);

      // Calculate from right to left
      bestscore =
        compute_nw_row(1,conspos, n, cons_base[conspos], r_cores[n], r_score,
                           lower_bound[n], upper_bound[n]-1, p_right,
                           &max_score_seq_idx, scoreParams,
                           MAXOFFSET, L, VERBOSE, NULL);
      if ( VERBOSE > 10 ) {
      printf("  RIGHT: bestscore = %d\n", bestscore);
      printf("         Gap: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", r_score[conspos % 2][n][offset + MAXOFFSET][1]);
      printf("\n");
      printf("         Sub: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", r_score[(conspos) % 2][n][offset + MAXOFFSET][0]);
      printf("\n");
      }

      bestscore =
        compute_nw_row(0,conspos, n, cons_base[conspos], l_cores[n], l_score,
                           lower_bound[n], upper_bound[n]-1, p_left,
                           &max_score_seq_idx, scoreParams,
                           MAXOFFSET, L, VERBOSE, NULL);

      if ( VERBOSE > 10 ){
      printf("  LEFT: bestscore = %d\n", bestscore);
      printf("         Gap: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", l_score[conspos % 2][n][offset + MAXOFFSET][1]);
      printf("\n");
      printf("         Sub: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", l_score[(conspos) % 2][n][offset + MAXOFFSET][0]);
      printf("\n");
      }

      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
      {
        mu_assert(
                  (( l_score[conspos % 2][n][MAXOFFSET-offset][0] ==
                    r_score[conspos % 2][n][MAXOFFSET-offset][0] ) &&
                  ( l_score[conspos % 2][n][MAXOFFSET-offset][1] ==
                    r_score[conspos % 2][n][MAXOFFSET-offset][1] )),
                  "Failed reciprocal test" );
        //if ( ( l_score[conspos % 2][n][MAXOFFSET-offset][0] !=
        //       r_score[conspos % 2][n][MAXOFFSET-offset][0] ) ||
        //     ( l_score[conspos % 2][n][MAXOFFSET-offset][1] !=
        //       r_score[conspos % 2][n][MAXOFFSET-offset][1] ) ) {
        //   printf("Failed reciprocal test row=%d\n", n);
        //}
      }
      // Validated maximum 26 at consensus row 6
      if ( n == 0 && conspos == 6 ) {
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
          mu_assert(
           (( r_score[conspos%2][0][offset + MAXOFFSET][1] ==
             final_right_gap[offset+MAXOFFSET] ) &&
           ( r_score[conspos%2][0][offset + MAXOFFSET][0] ==
             final_right_sub[offset+MAXOFFSET] )), "Failed final comparison" );
          //if ( ( r_score[3%2][0][offset + MAXOFFSET][1] !=
          //   final_right_gap[offset+MAXOFFSET] ) ||
          // ( r_score[3%2][0][offset + MAXOFFSET][0] !=
          //   final_right_sub[offset+MAXOFFSET] ) )
          //{
          //  printf("Failed final comparison at %d\n", offset+MAXOFFSET);
          //}
        }
      }
    }
  }
}

MU_TEST_SUITE(test_suite) {
        MU_RUN_TEST(test_global_align);
        MU_RUN_TEST(test_cons_seed_extend);
        MU_RUN_TEST(test_compute_nw_row);
}

void bnw_extend_test()
{
  MU_RUN_SUITE(test_suite);
  MU_REPORT();
}
