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
#include "score_system.h"

struct scoringSystem *
getMatrixUsingGapPenalties(char *matrixName, int gapopen_over, int gapextn_over)
{
  struct scoringSystem *dMat = getMatrix(matrixName);
  dMat->gapopen = gapopen_over;
  dMat->gapextn = gapextn_over;
  return ( dMat );
}


// TODO: This needs testing
double
calculateLambda( struct scoringSystem *scoringSys )
{
  int i,j;
  double lambda = 0.5;
  double lambda_upper = 0;
  double lambda_lower = 0;
  double check = 0;
  double sum = 0;

  do {
    check = 0;
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
    if ( sum < 1.0 ) {
      lambda_lower = lambda;
      lambda *= (double)2.0;
    }
  }while( sum < (double)1.0 );

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
  // Convenience table -- could be optimised
  static int alpha_convert[256];
  static char alphabet[] = "ACGT";
  static char lc_alphabet[] = "acgt";
  static int alphsize = 4;
  int i;
  for (i = 0; i < 256; i++)
    alpha_convert[i] = -1;

  for (i = 0; i < alphsize; i++)
  {
    alpha_convert[(int) alphabet[i]] = i;
    alpha_convert[(int) lc_alphabet[i]] = i;
  }

  alpha_convert['N'] = 99;
  alpha_convert['n'] = 99;

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
    matrix[alpha_convert['A']][alpha_convert['A']] = match;
    matrix[alpha_convert['A']][alpha_convert['C']] = mismatch;
    matrix[alpha_convert['A']][alpha_convert['G']] = mismatch;
    matrix[alpha_convert['A']][alpha_convert['T']] = mismatch;
    matrix[alpha_convert['A']][alpha_convert['N']] = mismatch;

    matrix[alpha_convert['G']][alpha_convert['A']] = mismatch;
    matrix[alpha_convert['G']][alpha_convert['C']] = mismatch;
    matrix[alpha_convert['G']][alpha_convert['G']] = match;
    matrix[alpha_convert['G']][alpha_convert['T']] = mismatch;
    matrix[alpha_convert['G']][alpha_convert['N']] = mismatch;

    matrix[alpha_convert['C']][alpha_convert['A']] = mismatch;
    matrix[alpha_convert['C']][alpha_convert['C']] = match;
    matrix[alpha_convert['C']][alpha_convert['G']] = mismatch;
    matrix[alpha_convert['C']][alpha_convert['T']] = mismatch;
    matrix[alpha_convert['C']][alpha_convert['N']] = mismatch;

    matrix[alpha_convert['T']][alpha_convert['A']] = mismatch;
    matrix[alpha_convert['T']][alpha_convert['C']] = mismatch;
    matrix[alpha_convert['T']][alpha_convert['G']] = mismatch;
    matrix[alpha_convert['T']][alpha_convert['T']] = match;
    matrix[alpha_convert['T']][alpha_convert['N']] = mismatch;

    matrix[alpha_convert['N']][alpha_convert['A']] = mismatch;
    matrix[alpha_convert['N']][alpha_convert['C']] = mismatch;
    matrix[alpha_convert['N']][alpha_convert['G']] = mismatch;
    matrix[alpha_convert['N']][alpha_convert['T']] = mismatch;
    matrix[alpha_convert['N']][alpha_convert['N']] = mismatch;

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
  // Convenience table -- could be optimised
  // making the following static causes segfault
  // in a multi-threaded context...why?
  //static int alpha_convert[256];
  int alpha_convert[256];
  static char alphabet[] = "ACGT";
  static char lc_alphabet[] = "acgt";
  static int alphsize = 4;
  int i;
  for (i = 0; i < 256; i++)
    alpha_convert[i] = -1;

  for (i = 0; i < alphsize; i++)
  {
    alpha_convert[(int) alphabet[i]] = i;
    alpha_convert[(int) lc_alphabet[i]] = i;
  }

  alpha_convert['N'] = 99;
  alpha_convert['n'] = 99;

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
    matrix[alpha_convert['A']][alpha_convert['A']] = 9;
    matrix[alpha_convert['A']][alpha_convert['C']] = -18;
    matrix[alpha_convert['A']][alpha_convert['G']] = -10;
    matrix[alpha_convert['A']][alpha_convert['T']] = -21;
    matrix[alpha_convert['A']][alpha_convert['N']] = -1;

    matrix[alpha_convert['G']][alpha_convert['A']] = -7;
    matrix[alpha_convert['G']][alpha_convert['C']] = -18;
    matrix[alpha_convert['G']][alpha_convert['G']] = 11;
    matrix[alpha_convert['G']][alpha_convert['T']] = -18;
    matrix[alpha_convert['G']][alpha_convert['N']] = -1;

    matrix[alpha_convert['C']][alpha_convert['A']] = -18;
    matrix[alpha_convert['C']][alpha_convert['C']] = 11;
    matrix[alpha_convert['C']][alpha_convert['G']] = -18;
    matrix[alpha_convert['C']][alpha_convert['T']] = -7;
    matrix[alpha_convert['C']][alpha_convert['N']] = -1;

    matrix[alpha_convert['T']][alpha_convert['A']] = -21;
    matrix[alpha_convert['T']][alpha_convert['C']] = -10;
    matrix[alpha_convert['T']][alpha_convert['G']] = -18;
    matrix[alpha_convert['T']][alpha_convert['T']] = 9;
    matrix[alpha_convert['T']][alpha_convert['N']] = -1;

    matrix[alpha_convert['N']][alpha_convert['A']] = -1;
    matrix[alpha_convert['N']][alpha_convert['C']] = -1;
    matrix[alpha_convert['N']][alpha_convert['G']] = -1;
    matrix[alpha_convert['N']][alpha_convert['T']] = -1;
    matrix[alpha_convert['N']][alpha_convert['N']] = -1;

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
    matrix[alpha_convert['A']][alpha_convert['A']] = 9;
    matrix[alpha_convert['A']][alpha_convert['C']] = -15;
    matrix[alpha_convert['A']][alpha_convert['G']] = -8;
    matrix[alpha_convert['A']][alpha_convert['T']] = -18;
    matrix[alpha_convert['A']][alpha_convert['N']] = -1;

    matrix[alpha_convert['G']][alpha_convert['A']] = -5;
    matrix[alpha_convert['G']][alpha_convert['C']] = -16;
    matrix[alpha_convert['G']][alpha_convert['G']] = 10;
    matrix[alpha_convert['G']][alpha_convert['T']] = -16;
    matrix[alpha_convert['G']][alpha_convert['N']] = -1;

    matrix[alpha_convert['C']][alpha_convert['A']] = -16;
    matrix[alpha_convert['C']][alpha_convert['C']] = 10;
    matrix[alpha_convert['C']][alpha_convert['G']] = -16;
    matrix[alpha_convert['C']][alpha_convert['T']] = -5;
    matrix[alpha_convert['C']][alpha_convert['N']] = -1;

    matrix[alpha_convert['T']][alpha_convert['A']] = -18;
    matrix[alpha_convert['T']][alpha_convert['C']] = -8;
    matrix[alpha_convert['T']][alpha_convert['G']] = -15;
    matrix[alpha_convert['T']][alpha_convert['T']] = 9;
    matrix[alpha_convert['T']][alpha_convert['N']] = -1;

    matrix[alpha_convert['N']][alpha_convert['A']] = -1;
    matrix[alpha_convert['N']][alpha_convert['C']] = -1;
    matrix[alpha_convert['N']][alpha_convert['G']] = -1;
    matrix[alpha_convert['N']][alpha_convert['T']] = -1;
    matrix[alpha_convert['N']][alpha_convert['N']] = -1;

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
    matrix[alpha_convert['A']][alpha_convert['A']] = 9;
    matrix[alpha_convert['A']][alpha_convert['C']] = -15;
    matrix[alpha_convert['A']][alpha_convert['G']] = -8;
    matrix[alpha_convert['A']][alpha_convert['T']] = -17;
    matrix[alpha_convert['A']][alpha_convert['N']] = -1;

    matrix[alpha_convert['G']][alpha_convert['A']] = -4;
    matrix[alpha_convert['G']][alpha_convert['C']] = -15;
    matrix[alpha_convert['G']][alpha_convert['G']] = 10;
    matrix[alpha_convert['G']][alpha_convert['T']] = -15;
    matrix[alpha_convert['G']][alpha_convert['N']] = -1;

    matrix[alpha_convert['C']][alpha_convert['A']] = -15;
    matrix[alpha_convert['C']][alpha_convert['C']] = 10;
    matrix[alpha_convert['C']][alpha_convert['G']] = -15;
    matrix[alpha_convert['C']][alpha_convert['T']] = -4;
    matrix[alpha_convert['C']][alpha_convert['N']] = -1;

    matrix[alpha_convert['T']][alpha_convert['A']] = -17;
    matrix[alpha_convert['T']][alpha_convert['C']] = -8;
    matrix[alpha_convert['T']][alpha_convert['G']] = -15;
    matrix[alpha_convert['T']][alpha_convert['T']] = 9;
    matrix[alpha_convert['T']][alpha_convert['N']] = -1;

    matrix[alpha_convert['N']][alpha_convert['A']] = -1;
    matrix[alpha_convert['N']][alpha_convert['C']] = -1;
    matrix[alpha_convert['N']][alpha_convert['G']] = -1;
    matrix[alpha_convert['N']][alpha_convert['T']] = -1;
    matrix[alpha_convert['N']][alpha_convert['N']] = -1;

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
    matrix[alpha_convert['A']][alpha_convert['A']] = 8;
    matrix[alpha_convert['A']][alpha_convert['C']] = -13;
    matrix[alpha_convert['A']][alpha_convert['G']] = -6;
    matrix[alpha_convert['A']][alpha_convert['T']] = -15;
    matrix[alpha_convert['A']][alpha_convert['N']] = -1;

    matrix[alpha_convert['G']][alpha_convert['A']] = -2;
    matrix[alpha_convert['G']][alpha_convert['C']] = -13;
    matrix[alpha_convert['G']][alpha_convert['G']] = 9;
    matrix[alpha_convert['G']][alpha_convert['T']] = -13;
    matrix[alpha_convert['G']][alpha_convert['N']] = -1;

    matrix[alpha_convert['C']][alpha_convert['A']] = -13;
    matrix[alpha_convert['C']][alpha_convert['C']] = 9;
    matrix[alpha_convert['C']][alpha_convert['G']] = -13;
    matrix[alpha_convert['C']][alpha_convert['T']] = -2;
    matrix[alpha_convert['C']][alpha_convert['N']] = -1;

    matrix[alpha_convert['T']][alpha_convert['A']] = -15;
    matrix[alpha_convert['T']][alpha_convert['C']] = -6;
    matrix[alpha_convert['T']][alpha_convert['G']] = -13;
    matrix[alpha_convert['T']][alpha_convert['T']] = 8;
    matrix[alpha_convert['T']][alpha_convert['N']] = -1;

    matrix[alpha_convert['N']][alpha_convert['A']] = -1;
    matrix[alpha_convert['N']][alpha_convert['C']] = -1;
    matrix[alpha_convert['N']][alpha_convert['G']] = -1;
    matrix[alpha_convert['N']][alpha_convert['T']] = -1;
    matrix[alpha_convert['N']][alpha_convert['N']] = -1;

    dMat->m_bg_freqs[0] = 0.285; // A
    dMat->m_bg_freqs[1] = 0.215; // C
    dMat->m_bg_freqs[2] = 0.215; // G
    dMat->m_bg_freqs[3] = 0.285; // T

  }else {
    printf("Custom matrices not supported ( yet ).  %s is not an internally coded matrix!\n", matrixName);
    exit(1);
  }

  dMat->m_lambda = calculateLambda(dMat);

  return ( dMat );
}

