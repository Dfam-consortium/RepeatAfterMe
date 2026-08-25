#ifndef __BNW_EXTEND_H__
#define __BNW_EXTEND_H__

#include "sequence.h"
#include "score_system.h"


//
// Structures
//
struct glocalSearchResult
{
    int score;              // Score (could be complex adjusted) of the result
    int query_start;        // Query start position relative to the query sequence provided
    int query_end;          // Query end position relative to the query sequence provided
    uint64_t subj_start;    // Subject start position relative to the sequenceLibrary used
    uint64_t subj_end;      // Subject end position relative to the sequenceLibrary used
};

//
// Signatures
//

int ****
allocate_score(int num_align, int bandwidth);

void
free_score(int num_align, int bandwidth, int ****score);

int
pair_seed_extend( int orient,
            uint64_t seq1_start, uint64_t seq1_end, uint64_t seq1_seed,
            uint64_t seq2_start, uint64_t seq2_end, uint64_t seq2_seed,
            int ****score, int l, struct sequenceLibrary *seqLib,
            struct scoringSystem *scoreParams, int bandwidth,
            struct glocalSearchResult *result, int use_complexity_adjust, int y_drop );

int
cons_seed_extend(char *query_seq, int query_len, int query_seed_start, int query_seed_end,
             struct sequenceLibrary *subj_lib, uint64_t subj_seed_start,
             uint64_t subj_seed_end, int orient, int ****score,
             struct scoringSystem *scoreParams, int bandwidth,
             struct glocalSearchResult *result, int use_complexity_adjust, int y_drop );

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
                    struct scoringSystem *dMatrix, int bandwidth, int L, int VERBOSE,
                    char *pathString);

struct coreAlignment *
fisher_yates_shuffle_calign(struct coreAlignment *coreAlign);

int
complexityAdjust(int score, int *comp, double lambda, double *matrix_freqs);

//
// UNIT TESTS
//
/*
MU_TEST(test_global_align)
{
   char alu_seqs[] ="ATATGGATAGCTAGCGTGCGACGTACGGCGATTGGTATATGAGCGATATC"
                    //ALU Portion
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

  int bandwidth = 10;
  int ****scores = allocate_score(1, bandwidth);
  struct scoringSystem *scoreParams = getMatrix("20p43g");
  int orient = 0;
  uint64_t seq1_start = 50;
  uint64_t seq1_end = 368;
  uint64_t seq1_seed = 60;
  uint64_t seq2_start = 500;
  uint64_t seq2_end = 820;
  uint64_t seq2_seed = 510;
  int l = 8;

  int align_score =  global_align(orient, seq1_start, seq1_end, seq1_seed,
             seq2_start, seq2_end, seq2_seed,
             scores, l, &seqLib, scoreParams, bandwidth);

  free(seqLib.boundaries);
  free(scores);
  freeScoringSystem(scoreParams);
  mu_check( align_score == 2746 );
}
*/

//MU_TEST(test_global_align);
//void test_global_align(void);
void bnw_extend_test();

#endif
