#ifndef __COMMON_H__
#define __COMMON_H__

// An enumeration of coreAlignment extension boundary types:
//   L_BOUNDARY: A boundary limited by the maximum extension length permited (L parameter).
//   SEQ_BOUNDARY: A boundary limited by the sequence length.
//   CORE_BOUNDARY: A boundary limited by the presense of a neighboring core.
//   EXT_BOUNDARY: A boundary limited by the first extension phase.
//
enum CoreBoundFlag { L_BOUNDARY = 0, SEQ_BOUNDARY = 1, CORE_BOUNDARY = 2, EXT_BOUNDARY = 3 };

//
// coreAlignment
//
// This structure holds details of a core aligned region from which
// extension will be performed.  This may be as simple as an kmer
// where all entries are perfect matches to a specific pattern of
// fixed length, or a complex set of sequences from a multiple
// alignment.  In the kmer case all cores could reasonably be
// extended in either direction, whereas in the later case sequences
// may be not aligned up to the left or right edges and may have
// their left/rightExtendable flags disabled.
//
// The left/right sequence positions are indices into the in-memory
// seqLib->sequence[] array (0-based,fully closed) and represent
// the logical left/right ordering of a single core region with
// respect to the larger set of core alignmehts.
//
// NOTE: To translate these back to seqID coordinates, the seqLib->offsets[]
// array should be used.
//
// E.g.
//
//  Kmer Case:
//                       LEFT          RIGHT
//             seqIdx=0  [10] ACCTTAAG [17]
//             seqIdx=1  [32] ACCTTAAG [25] Reverse Strand Core
//             seqIdx=2  [42] ACCTTAAG [49]
//
//                       Where the 7bp sequence at zero-based sequence
//                       postion 10 = "ACCTTAAG"
//
//          coreAlign #1: leftSeqPos = 10, rightSeqPos = 17, orient = 0
//                        leftExtendable = 1, rightExtendable = 1
//          coreAlign #2: leftSeqPos = 32, rightSeqPos = 25, orient = 1
//                        leftExtendable = 1, rightExtendable = 1
//          coreAlign #3: leftSeqPos = 42, rightSeqPos = 49, orient = 0
//                        leftExtendable = 1, rightExtendable = 1
//
//          SeqLib->boundaries[1] = 20
//          SeqLib->offsets[1] = 100
//
//          For the second core, the actual in-memory sequence would
//          be found at SeqLib->sequence[32] representing the position
//          (startPos - seqLib->boundaries[1]) + seqLib->offsets[1] =
//          (32-20)+100 = seqIdx[1]:112
//
//  Multiple Alignment Case:
//
//                  [10] GAGGATT............. [16]
//                  [22] GAGGATTAACACGTGGAATA [41]
//                  [51] ...............GAATA [55]
//
//          coreAlign #1: leftSeqPos = 10, rightSeqPos = 16, orient = 0
//                        leftExtendable = 1, rightExtendable = 0
//          coreAlign #2: leftSeqPos = 22, rightSeqPos = 41, orient = 0
//                        leftExtendable = 1, rightExtendable = 1
//          coreAlign #3: leftSeqPos = 51, rightSeqPos = 55, orient = 0
//                        leftExtendable = 0, rightExtendable = 1
//
//  The left/rightExtendable flags are set when extension does not make
//  sense based on the aligned position of the core.  A true value does
//  not indicate that an extension is feasible ( e.g. there is adequate
//  flanking sequence available ).
//
//  The orientation of the aligned core is apparent from the order of
//  the sequence positions, however the orient flag is currently
//  required.
//
struct coreAlignment
{
    struct coreAlignment *next; // For use as a linked-list ( NULL terminated )
    int seqIdx;                 // Index into the sequence_ident list
    uint64_t leftSeqPos;        // The left-most position of the core alignment (sequence[] coordinates)
    uint64_t rightSeqPos;       // The right-most position of the core alignment
    char leftExtendable;        // Is left side aligned to the end and thus extendable? (1 = true)
    char rightExtendable;       // Is right side aligned to the end and thus extendable? (1 = true)
    uint64_t lowerSeqBound;     // During extension, this stores the in-memory lower bound of the sequence.
    uint64_t upperSeqBound;     // During extension, this stores the in-memory upper bound of the sequence.
    enum CoreBoundFlag lowerSeqBoundFlag;     // A flag to indicate type of boundary represented by lowerSeqBound.
                                //    0 = L-limited boundary, 1 = sequence limited bound,
                                //    2 = Core limited boundary, or 3 = extension limited boundary.
    enum CoreBoundFlag upperSeqBoundFlag;     // A flag to indicate type of boundary represented by upperSeqBound.
    int  leftExtensionLen;      // After extension, the length of the left extension in bp.
    int  rightExtensionLen;     // After extension, the length of the right extension in bp.
    int  score;                 // After extension, the score of the right+left extension.
    char orient;                // Sequence orientation ( 1 = reverse )
};

//
// Prototypes
//
uint32_t small_kmer_to_int (char *kmer, int kmer_k, int rc);
uint32_t small_kmer_seq_to_int (char *kmer, int kmer_k, int rc);
uint64_t large_kmer_to_int (char *kmer, int kmer_k, int rc);
uint64_t large_kmer_seq_to_int (char *kmer, int kmer_k, int rc);
uint64_t large_kmer_pair_to_int (char *leftKmer, char *rightKmer, int kmerSize, int rc);
void print_small_kmer_int(uint32_t val, int kmer_k);




int rand_int(int n);

double compute_entropy(char *seq, int kmer_k);
double compute_seq_entropy(char *seq);

#endif
