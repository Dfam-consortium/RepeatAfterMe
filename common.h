#ifndef __COMMON_H__
#define __COMMON_H__

int rand_int(int n);
double compute_seq_entropy(char *seq);

//
// coreAlignment
//
// This structure holds details of a core aligned region from which
// extension will be performed.  This may be as simple as an lmer
// where all entries are perfect matches to a specific pattern of
// fixed length, or a complex set of sequences from a multiple
// alignment.  In the lmer case all cores could reasonably be
// extended in either direction, whereas in the later case sequences
// may be not aligned up to the left or right edges and may have
// their left/rightExtendable flags disabled.
//
// The left/right sequence positions are relative to a particular
// sequenceLibrary and represent the logical left/right sequence
// positions ( 0-based ) of a single core region with respect to the
// larger set of core alignments.  E.g.
//
//  Lmer Case:
//                       LEFT          RIGHT
//                       [10] ACCTTAAG [17]
//                       [32] ACCTTAAG [25] Reverse Strand Core
//                       [42] ACCTTAAG [49]
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
    uint64_t leftSeqPos;        // The logical left-most position of the core alignment
    uint64_t rightSeqPos;       // The logical right-most position of the core alignment
    char leftExtendable;        // Is left side aligned to the end and thus extendable? (1 = true)
    char rightExtendable;       // Is right side aligned to the end and thus extendable? (1 = true)
    int  leftExtensionLen;      // After extension, the length of the left extension in bp.
    int  rightExtensionLen;     // After extension, the length of the right extension in bp.
    char orient;                // Sequence orientation ( 1 = reverse )
};

#endif
