#ifndef __BNW_EXTEND_H__
#define __BNW_EXTEND_H__

#include "common.h"
#include "sequence.h"
#include "score_system.h"

//
//
//
struct glocalSearchResult
{
    int score;              // Score (could be complex adjusted) of the result
    int query_start;        // Query start position relative to the query sequence provided
    int query_end;          // Query end position relative to the query sequence provided
    uint64_t subj_start;    // Subject start position relative to the sequenceLibrary used
    uint64_t subj_end;      // Subject end position relative to the sequenceLibrary used
};


int ****
allocate_score(int num_align, int bandwidth);

void
free_score(int num_align, int bandwidth, int ****score);

void
test_cons_seed_extend();

void
test_global_align();

int
pair_seed_extend( int orient,
            uint64_t seq1_start, uint64_t seq1_end, uint64_t seq1_seed,
            uint64_t seq2_start, uint64_t seq2_end, uint64_t seq2_seed,
            int ****score, int l, struct sequenceLibrary *seqLib,
            struct scoringSystem *scoreParams, int bandwidth,
            struct glocalSearchResult *result, int use_complexity_adjust);

int
cons_seed_extend(char *query_seq, int query_len, int query_seed_start, int query_seed_end,
             struct sequenceLibrary *subj_lib, uint64_t subj_seed_start,
             uint64_t subj_seed_end, int orient, int ****score,
             struct scoringSystem *scoreParams, int bandwidth,
             struct glocalSearchResult *result, int use_complexity_adjust );

int
global_align(int orient, uint64_t seq1_start, uint64_t seq1_end, uint64_t seq1_seed,
             uint64_t seq2_start, uint64_t seq2_end, uint64_t seq2_seed,
             int ****score, int l, struct sequenceLibrary *seqLib,
             struct scoringSystem *scoreParams, int bandwidth);


int
compute_nw_row(int direction, int row_idx, int n, char cons_base, struct coreAlignment *coreAlign,
                    int ****score, uint64_t lower_seq_bound,
                    uint64_t upper_seq_bound, char *sequence,
                    int *max_score_sequence_idx,
                    struct scoringSystem *dMatrix, int bandwidth, int L, int VERBOSE);

struct coreAlignment *
fisher_yates_shuffle_calign(struct coreAlignment *coreAlign);

int
complexityAdjust(int score, int *comp, double lambda, double *matrix_freqs);

#endif
