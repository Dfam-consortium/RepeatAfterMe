#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "common.h"

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
    uint64_t length;            // The length of all sequences in bp
    int count;                  // The number of sequences loaded
};

struct twoBitSpec *
readBEDRanges(char *twoBitFile, char *bedFile);


struct sequenceLibrary *
loadSequences(char *twoBitName, char **softmasked);

struct sequenceLibrary *
loadSequenceSubset(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *N);

struct AwFmIndex *
createFMIndex(char *twoBitName);

#define IUPAC(c) c == 'R' || c=='r' || c=='Y' || c=='y' || c=='M' || c=='m' || c=='K' || c=='k' || c=='W' || c=='w' || c=='S' || c=='s' || c=='B' || c=='b' || c=='D' || c=='d' || c=='H' || c=='h' || c=='V' || c=='v'

char num_to_char(char z);
char char_to_num(char c);
char compl( char c );
#endif
