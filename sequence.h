#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "common.h"

// Sequence encoding
#define SYM_A 0
#define SYM_C 1
#define SYM_G 2
#define SYM_T 3
#define SYM_a 4
#define SYM_c 5
#define SYM_g 6
#define SYM_t 7
#define SYM_N 99

//
// A structure to hold the library of sequences. Sequences
// are stored in a single array with each subsequent sequence
// concatenated to the last without the use of any sentinel
// or spacer positions.
//
// The boundaries array holds the cummulative sums of the
// sequence lengths in the seuence array.  This is useful
// for finding the start/end position of each concatenated
// sequence:
//     seq1: 0 to boundaries[0]-1
//     seq2: boundaries[0] to boundaries[1]-1
//     ..
//     seqN: boundaries[N-2] to boundaries[N-1]-1
// The boundaries array is terminated with a value of 0.
//
// The identifiers array holds the identifier strings
// for each sequence.
//
struct sequenceLibrary
{
    char * sequence;            // Encoded DNA sequences - concatenated
    char **identifiers;         // Sequence identifiers
    uint64_t *boundaries;       // The start index in the sequence array for sequence
    uint64_t *offsets;          // If library contains subsequences these are the offsets into the full sequences
    uint64_t length;            // The length of all sequences in bp
    int count;                  // The number of sequences loaded
    int markov_chain_order;              // The max order calculated in markov_chain_prob_tables (if > 0)
    uint32_t **markov_chain_prob_tables; // Frequencies of all kmers up to the Markov chain order
};

#define MARKOV_PROB_SCALE_FACTOR 10000

struct twoBitSpec *
readBEDRanges(char *twoBitFile, char *bedFile);

void seqlib_destroy(struct sequenceLibrary *seqLib);

struct sequenceLibrary *
loadSequences(char *twoBitName, char **softmasked, int markovChainOrder);

struct sequenceLibrary *
loadSequenceSubset(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *N);

struct sequenceLibrary *
loadSequenceSubsetMinimal(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *num_cores, int max_flanking_bp);

struct AwFmIndex *
createFMIndex(char *twoBitName);

#define IUPAC(c) c == 'R' || c=='r' || c=='Y' || c=='y' || c=='M' || c=='m' || c=='K' || c=='k' || c=='W' || c=='w' || c=='S' || c=='s' || c=='B' || c=='b' || c=='D' || c=='d' || c=='H' || c=='h' || c=='V' || c=='v'

char num_to_char(char z);
char num_to_char_compl(char z);
char char_to_num(char c);
char compl( char c );
char encoded_toupper(char c);
char mask(char c);
int
kmermatch(char *kmer1, char *kmer2, int l);     /* forward match */
int kmermatcheither(char *kmer1, char *kmer2, int mismatches, int l); /* forward or rc match */
int kmermatchrc(char *kmer1, char *kmer2, int mismatches, int l); /* rc match */


#endif
