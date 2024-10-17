/*
 * score_system.c
 *
 * Module to hold the scoring system for sequence alignment
 *
 * TODO: Require that the matrix be contiguous!  Right now we
 *       have N=99.  This should be set to the next column.
 *       Changes to the allocation code and seqeunce.c as well
 *       as the unit test need to be made.
 *
 * Authors:
 *    Robert Hubley
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "minunit.h"
#include "score_system.h"
#include "sequence.h"

struct scoringSystem *
getMatrixUsingGapPenalties(char *matrixName, int gapopen_over, int gapextn_over)
{
  struct scoringSystem *dMat = getMatrix(matrixName);
  dMat->gapopen = gapopen_over;
  dMat->gapextn = gapextn_over;
  return ( dMat );
}


//
// Calculate the substitition matrix lambda value
//
//
double
calculateLambda( struct scoringSystem *scoringSys )
{
  int i,j;
  double lambda = 0.5;
  double lambda_upper = 0;
  double lambda_lower = 0;
  double check = 0;
  double sum = 0;

  for (;;) {
    sum = 0;
    check = 0;
    // Assumes base encoding of ACGT to be 0-3
    for ( i = 0; i < 4; i++ ){
      for ( j = 0; j < 4; j++ ){
        sum += scoringSys->m_bg_freqs[i] * scoringSys->m_bg_freqs[j] *
               exp( lambda * scoringSys->matrix[i][j] );
        check += scoringSys->m_bg_freqs[i] * scoringSys->m_bg_freqs[j];
      }
    }
    if ( check > (double)1.001 || check < (double)0.999 ) {
      return ( (double)-1.0 );
    }
    if ( sum >= (double)1.0 ) break;
    lambda_lower = lambda;
    lambda *= (double)2.0;
  }
  lambda_upper = lambda;

  while( lambda_upper - lambda_lower > (double)0.00001 ) {
    lambda = ( lambda_lower + lambda_upper ) / (double) 2.0;

    sum = 0;
    check = 0;
    for ( i = 0; i < 4; i++ ) {
      for ( j = 0; j < 4; j++ ) {
        sum += scoringSys->m_bg_freqs[i] * scoringSys->m_bg_freqs[j] *
               exp( lambda * scoringSys->matrix[i][j] );
        check += scoringSys->m_bg_freqs[i] * scoringSys->m_bg_freqs[j];
      }
    }
    if ( check > (double)1.001 || check < (double)0.999 ) {
      return ( (double)-1.0 );
    }
    if ( sum >= 1.0 ) {
      lambda_upper = lambda;
    }else {
      lambda_lower = lambda;
    }
  }
  return (lambda);
}

struct scoringSystem *
getRepeatScoutMatrix(int match, int mismatch, int gap)
{
  int i;
  int j;

  // Currently set to support 0-99 because 'N' is hardcoded as 99.
  int **matrix = (int **)malloc(100*sizeof(int*));
  if ( matrix == NULL )
  {
    printf("Could not allocate space for dna scoring matrix!\n");
    exit(1);
  }

  for(i=0;i<100;i++){
    matrix[i] = (int *)malloc(100*sizeof(int));
    if ( matrix[i] == NULL ) {
      printf("Could not allocate space for dna scoring matrix row!\n");
      exit(1);
    }
  }

  struct scoringSystem *dMat = (struct scoringSystem *)malloc(sizeof(struct scoringSystem));
  dMat->matrix = matrix;
  dMat->name = "repeatscout";
  dMat->msize = 100;
  dMat->alphabet = "ACGTN";
  dMat->gapopen = 0;
  dMat->gapextn = gap;
  matrix[SYM_A][SYM_A] = match;
  matrix[SYM_A][SYM_C] = mismatch;
  matrix[SYM_A][SYM_G] = mismatch;
  matrix[SYM_A][SYM_T] = mismatch;
  matrix[SYM_A][SYM_N] = mismatch;

  matrix[SYM_G][SYM_A] = mismatch;
  matrix[SYM_G][SYM_C] = mismatch;
  matrix[SYM_G][SYM_G] = match;
  matrix[SYM_G][SYM_T] = mismatch;
  matrix[SYM_G][SYM_N] = mismatch;

  matrix[SYM_C][SYM_A] = mismatch;
  matrix[SYM_C][SYM_C] = match;
  matrix[SYM_C][SYM_G] = mismatch;
  matrix[SYM_C][SYM_T] = mismatch;
  matrix[SYM_C][SYM_N] = mismatch;

  matrix[SYM_T][SYM_A] = mismatch;
  matrix[SYM_T][SYM_C] = mismatch;
  matrix[SYM_T][SYM_G] = mismatch;
  matrix[SYM_T][SYM_T] = match;
  matrix[SYM_T][SYM_N] = mismatch;

  matrix[SYM_N][SYM_A] = mismatch;
  matrix[SYM_N][SYM_C] = mismatch;
  matrix[SYM_N][SYM_G] = mismatch;
  matrix[SYM_N][SYM_T] = mismatch;
  matrix[SYM_N][SYM_N] = mismatch;

  // Penalize matches from/to a->t (masked) the same as to/from an N
  //   - Assumes ordering from A->T and a->t
  for( i = SYM_a; i <= SYM_t; i++ ) {
    for( j = SYM_a; j <= SYM_t; j++ ) {
      matrix[i][j] = -1;
      matrix[j][i] = -1;
    }
    for( j = SYM_A; j <= SYM_T; j++ ) {
      matrix[i][j] = mismatch;
      matrix[j][i] = mismatch;
    }
    matrix[i][SYM_N] = mismatch;
    matrix[SYM_N][i] = mismatch;
  }

  dMat->m_bg_freqs[0] = 0.25;
  dMat->m_bg_freqs[1] = 0.25;
  dMat->m_bg_freqs[2] = 0.25;
  dMat->m_bg_freqs[3] = 0.25;
  dMat->m_lambda = calculateLambda(dMat);

  return ( dMat );
}

void freeScoringSystem(struct scoringSystem *score){
  int i;
  for(i=0;i<score->msize;i++)
    free(score->matrix[i]);
  free(score->matrix);
  free(score);
}

struct scoringSystem *
getMatrix(char *matrixName)
{
  int i,j;

  int **matrix = (int **)malloc(100*sizeof(int*));
  if ( matrix == NULL )
  {
    printf("Could not allocate space for dna scoring matrix!\n");
    exit(1);
  }

  for(i=0;i<100;i++){
    matrix[i] = (int *)malloc(100*sizeof(int));
    if ( matrix[i] == NULL ) {
      printf("Could not allocate space for dna scoring matrix row!\n");
      exit(1);
    }
  }

  struct scoringSystem *dMat = (struct scoringSystem *)malloc(sizeof(struct scoringSystem));
  dMat->matrix = matrix;
  dMat->msize = 100;
  dMat->alphabet = "ACGTN";

  if (strcmp(matrixName, "14p43g") == 0)
  {
    dMat->name = "14p43g";
    dMat->gapopen = -33;
    dMat->gapextn = -7;
    // 14p43g non-symmetric matrix
    // matrix[consensus_base][sequence_base]
    matrix[SYM_A][SYM_A] = 9;
    matrix[SYM_A][SYM_C] = -18;
    matrix[SYM_A][SYM_G] = -10;
    matrix[SYM_A][SYM_T] = -21;
    matrix[SYM_A][SYM_N] = -1;

    matrix[SYM_G][SYM_A] = -7;
    matrix[SYM_G][SYM_C] = -18;
    matrix[SYM_G][SYM_G] = 11;
    matrix[SYM_G][SYM_T] = -18;
    matrix[SYM_G][SYM_N] = -1;

    matrix[SYM_C][SYM_A] = -18;
    matrix[SYM_C][SYM_C] = 11;
    matrix[SYM_C][SYM_G] = -18;
    matrix[SYM_C][SYM_T] = -7;
    matrix[SYM_C][SYM_N] = -1;

    matrix[SYM_T][SYM_A] = -21;
    matrix[SYM_T][SYM_C] = -10;
    matrix[SYM_T][SYM_G] = -18;
    matrix[SYM_T][SYM_T] = 9;
    matrix[SYM_T][SYM_N] = -1;

    matrix[SYM_N][SYM_A] = -1;
    matrix[SYM_N][SYM_C] = -1;
    matrix[SYM_N][SYM_G] = -1;
    matrix[SYM_N][SYM_T] = -1;
    matrix[SYM_N][SYM_N] = -1;

    dMat->m_bg_freqs[0] = 0.285; // A
    dMat->m_bg_freqs[1] = 0.215; // C
    dMat->m_bg_freqs[2] = 0.215; // G
    dMat->m_bg_freqs[3] = 0.285; // T

  }
  else if (strcmp(matrixName, "18p43g") == 0)
  {
    dMat->gapopen = -30;
    dMat->gapextn = -6;
    // 18p43g non-symmetric matrix
    // matrix[consensus_base][sequence_base]
    matrix[SYM_A][SYM_A] = 9;
    matrix[SYM_A][SYM_C] = -15;
    matrix[SYM_A][SYM_G] = -8;
    matrix[SYM_A][SYM_T] = -18;
    matrix[SYM_A][SYM_N] = -1;

    matrix[SYM_G][SYM_A] = -5;
    matrix[SYM_G][SYM_C] = -16;
    matrix[SYM_G][SYM_G] = 10;
    matrix[SYM_G][SYM_T] = -16;
    matrix[SYM_G][SYM_N] = -1;

    matrix[SYM_C][SYM_A] = -16;
    matrix[SYM_C][SYM_C] = 10;
    matrix[SYM_C][SYM_G] = -16;
    matrix[SYM_C][SYM_T] = -5;
    matrix[SYM_C][SYM_N] = -1;

    matrix[SYM_T][SYM_A] = -18;
    matrix[SYM_T][SYM_C] = -8;
    matrix[SYM_T][SYM_G] = -15;
    matrix[SYM_T][SYM_T] = 9;
    matrix[SYM_T][SYM_N] = -1;

    matrix[SYM_N][SYM_A] = -1;
    matrix[SYM_N][SYM_C] = -1;
    matrix[SYM_N][SYM_G] = -1;
    matrix[SYM_N][SYM_T] = -1;
    matrix[SYM_N][SYM_N] = -1;

    dMat->m_bg_freqs[0] = 0.285; // A
    dMat->m_bg_freqs[1] = 0.215; // C
    dMat->m_bg_freqs[2] = 0.215; // G
    dMat->m_bg_freqs[3] = 0.285; // T

  }
  else if (strcmp(matrixName, "20p43g") == 0)
  {
    dMat->name = "20p43g";
    dMat->gapopen = -28;
    dMat->gapextn = -5;
    // 20p43g non-symmetric matrix
    // matrix[consensus_base][sequence_base]
    matrix[SYM_A][SYM_A] = 9;
    matrix[SYM_A][SYM_C] = -15;
    matrix[SYM_A][SYM_G] = -8;
    matrix[SYM_A][SYM_T] = -17;
    matrix[SYM_A][SYM_N] = -1;

    matrix[SYM_G][SYM_A] = -4;
    matrix[SYM_G][SYM_C] = -15;
    matrix[SYM_G][SYM_G] = 10;
    matrix[SYM_G][SYM_T] = -15;
    matrix[SYM_G][SYM_N] = -1;

    matrix[SYM_C][SYM_A] = -15;
    matrix[SYM_C][SYM_C] = 10;
    matrix[SYM_C][SYM_G] = -15;
    matrix[SYM_C][SYM_T] = -4;
    matrix[SYM_C][SYM_N] = -1;

    matrix[SYM_T][SYM_A] = -17;
    matrix[SYM_T][SYM_C] = -8;
    matrix[SYM_T][SYM_G] = -15;
    matrix[SYM_T][SYM_T] = 9;
    matrix[SYM_T][SYM_N] = -1;

    matrix[SYM_N][SYM_A] = -1;
    matrix[SYM_N][SYM_C] = -1;
    matrix[SYM_N][SYM_G] = -1;
    matrix[SYM_N][SYM_T] = -1;
    matrix[SYM_N][SYM_N] = -1;

    dMat->m_bg_freqs[0] = 0.285; // A
    dMat->m_bg_freqs[1] = 0.215; // C
    dMat->m_bg_freqs[2] = 0.215; // G
    dMat->m_bg_freqs[3] = 0.285; // T

  }
  else if (strcmp(matrixName, "25p43g") == 0)
  {
    dMat->name = "25p43g";
    dMat->gapopen = -25;
    dMat->gapextn = -5;
    // 25p43g non-symmetric matrix
    // matrix[consensus_base][sequence_base]
    matrix[SYM_A][SYM_A] = 8;
    matrix[SYM_A][SYM_C] = -13;
    matrix[SYM_A][SYM_G] = -6;
    matrix[SYM_A][SYM_T] = -15;
    matrix[SYM_A][SYM_N] = -1;

    matrix[SYM_G][SYM_A] = -2;
    matrix[SYM_G][SYM_C] = -13;
    matrix[SYM_G][SYM_G] = 9;
    matrix[SYM_G][SYM_T] = -13;
    matrix[SYM_G][SYM_N] = -1;

    matrix[SYM_C][SYM_A] = -13;
    matrix[SYM_C][SYM_C] = 9;
    matrix[SYM_C][SYM_G] = -13;
    matrix[SYM_C][SYM_T] = -2;
    matrix[SYM_C][SYM_N] = -1;

    matrix[SYM_T][SYM_A] = -15;
    matrix[SYM_T][SYM_C] = -6;
    matrix[SYM_T][SYM_G] = -13;
    matrix[SYM_T][SYM_T] = 8;
    matrix[SYM_T][SYM_N] = -1;

    matrix[SYM_N][SYM_A] = -1;
    matrix[SYM_N][SYM_C] = -1;
    matrix[SYM_N][SYM_G] = -1;
    matrix[SYM_N][SYM_T] = -1;
    matrix[SYM_N][SYM_N] = -1;

    dMat->m_bg_freqs[0] = 0.285; // A
    dMat->m_bg_freqs[1] = 0.215; // C
    dMat->m_bg_freqs[2] = 0.215; // G
    dMat->m_bg_freqs[3] = 0.285; // T

  }else {
    printf("Custom matrices not supported ( yet ).  %s is not an internally coded matrix!\n", matrixName);
    exit(1);
  }

  // Penalize matches from/to a->t (masked) the same as to/from an N
  //   - Assumes ordering from A->T and a->t
  for( i = SYM_a; i <= SYM_t; i++ ) {
    for( j = SYM_a; j <= SYM_t; j++ ) {
      matrix[i][j] = -1;
      matrix[j][i] = -1;
    }
    for( j = SYM_A; j <= SYM_T; j++ ) {
      matrix[i][j] = -1;
      matrix[j][i] = -1;
    }
    matrix[i][SYM_N] = -1;
    matrix[SYM_N][i] = -1;
  }

  dMat->m_lambda = calculateLambda(dMat);

  return ( dMat );
}


//
// UNIT TESTS
//

MU_TEST(test_calculate_lambda)
{
  int i;
  int j;
  //                   A    C    G    T    N
  int residues[5] = {  0,   1,   2,   3,   99 };
  int **matrix = (int **)malloc(100*sizeof(int*));
  for(i=0;i<100;i++)
    matrix[i] = (int *)malloc(100*sizeof(int));

  struct scoringSystem *dMat = (struct scoringSystem *)malloc(sizeof(struct scoringSystem));
  dMat->matrix = matrix;
  dMat->name = "plus_one_minus_one";
  dMat->msize = 100;
  dMat->alphabet = "ACGTN";
  dMat->gapopen = 0;
  dMat->gapextn = 0;
  dMat->m_bg_freqs[0] = 0.25;
  dMat->m_bg_freqs[1] = 0.25;
  dMat->m_bg_freqs[2] = 0.25;
  dMat->m_bg_freqs[3] = 0.25;

  int match = 1;
  int mismatch = -1;
  for ( i = 0; i < 5; i++ )
    for ( j = 0; i < 5; i++ )
      if ( i == j )
        matrix[residues[i]][residues[j]] = match;
      else
        matrix[residues[i]][residues[j]] = mismatch;

  dMat->m_lambda = calculateLambda(dMat);

  // 1.098609924316
  mu_check(abs(1.09860-dMat->m_lambda) <= 0.00001);

  match = 1;
  mismatch = -3;
  for ( i = 0; i < 5; i++ )
    for ( j = 0; i < 5; i++ )
      if ( i == j )
        matrix[residues[i]][residues[j]] = match;
      else
        matrix[residues[i]][residues[j]] = mismatch;

  dMat->m_lambda = calculateLambda(dMat);

  // 1.374061584473
  mu_check(abs(1.37406-dMat->m_lambda) <= 0.00001);

}

MU_TEST_SUITE(test_suite) {
        MU_RUN_TEST(test_calculate_lambda);
}

void score_system_test()
{
  MU_RUN_SUITE(test_suite);
  MU_REPORT();
}

