#ifndef __EXTEND_ALIGN_H__
#define __EXTEND_ALIGN_H__

void print_parameters();

void test_compute_nw_row();

int
extend_alignment(int direction, struct coreAlignment *coreAlign, int ****score,
                 struct sequenceLibrary *seqLib, char *master, int MAXOFFSET, int CAPPENALTY,
                 int MINIMPROVEMENT, int L, int N, int *bestbestscore,
                 struct scoringSystem *scoreParams);

void printSeedRegion ( int seedLastOcc );

#endif
