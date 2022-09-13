#ifndef __SCORE_SYSTEM_H__
#define __SCORE_SYSTEM_H__

struct scoringSystem
{
    char *name;               // Matrix name
    int **matrix;             // 2D scoring matrix
    int msize;                // The size of the matrix/alphabet
    char *alphabet;           // Matrix alphabet
    double m_lambda;          // Matrix lambda
    double m_bg_freqs[4];     // Matrix background frequencies [ A, C, G, T ]
    int gapopen;              // Gap open penalty
    int gapextn;              // Gap extension penalty
};


double
calculateLambda( struct scoringSystem *scoringSys );

void freeScoringSystem(struct scoringSystem *score);

struct scoringSystem *
getMatrixUsingGapPenalties(char *matrixName, int gapopen_over, int gapextn_over);

struct scoringSystem *
getRepeatScoutMatrix(int match, int mismatch, int gap);

struct scoringSystem *
getMatrix(char *matrixName);

#endif
