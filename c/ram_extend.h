#ifndef __RAM_EXTEND_H__
#define __RAM_EXTEND_H__

//
// Signatures
//
void print_parameters();

int
extend_alignment(int direction, struct coreAlignment *coreAlign, int ****score,
                 struct sequenceLibrary *seqLib, char *master, int MAXOFFSET, int CAPPENALTY,
                 int MINIMPROVEMENT, int L, int N, 
                 struct scoringSystem *scoreParams, FILE *pathStringFile);

void printSeedRegion ( int seedLastOcc );

#endif
