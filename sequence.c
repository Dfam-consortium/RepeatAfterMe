/*
 * sequence.c
 *
 * Authors:
 *    Robert Hubley
 *    Based on RepeatScout code by Alkes Price 2005
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//typedef unsigned int uint;
#include "bnw_extend.h"
#include "sequence.h"
#include "kentsrc/common.h"
#include "kentsrc/hash.h"
#include "kentsrc/linefile.h"
#include "kentsrc/dnaseq.h"
#include "kentsrc/twoBitNew.h"



/*
 * loadSequences()
 *
 *
 *
 */
struct sequenceLibrary *
loadSequences(char *twoBitName, char **softmasked)
{
  struct twoBitFile *tbf;
  uint64_t total_seq_size = 0;

  char *sequence = NULL;
  uint64_t *boundaries = NULL;
  char **sequence_idents = NULL;

  int i;
  // Must start from 1 otherwise it looks like the
  // NULL return value.
  int seq_idx = 1;

  tbf = twoBitOpen(twoBitName);
  struct twoBitIndex *index = NULL;
  for (index = tbf->indexList; index != NULL; index = index->next)
  {
    // Complete sequence
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, index->name, 0, 0);
    //printf("  seq size = %d\n", seq->size);

    sequence_idents =
      (char **) realloc(sequence_idents,
                        (seq_idx + 1) * sizeof(char *));
    if (NULL == sequence_idents)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "sequence_idents array\n");
      exit(1);
    }
    boundaries =
      (uint64_t *) realloc(boundaries,
                           (seq_idx + 1) * sizeof(uint64_t));
    if (NULL == boundaries)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "boundaries array\n");
      exit(1);
    }
    sequence =
      (char *) realloc(sequence,
                       (total_seq_size + seq->size + 1) * sizeof(char));
    if (NULL == sequence)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "sequence array\n");
      exit(1);
    }

    // The presence of **softmasked indicates that we desire to use
    // the presence of lowercase base pairs as a softmask for
    // downstream analysis.  The softmasked array ( currently
    // char *, but should be a bitarray ) holds a 1 if it's masked
    // or a 0 if not.
    if ( softmasked != NULL )
    {
      *softmasked =
        (char *) realloc(*softmasked,
                         (total_seq_size + seq->size + 1) * sizeof(char));
      if (NULL == *softmasked)
      {
        fprintf(stderr, "Could not allocated more space for the "
                "softmasked array\n");
        exit(1);
      }
      for (i = 0; i < seq->size; i++)
        if ((seq->dna)[i] == 'c' || (seq->dna)[i] == 'a' ||
             (seq->dna)[i] == 'g' || (seq->dna)[i] == 't' )
          (*softmasked)[total_seq_size + i] = 1;
        else
          (*softmasked)[total_seq_size + i] = 0;
    }

    // identifier
    sequence_idents[seq_idx - 1] = cloneString(seq->name);
    sequence_idents[seq_idx] = NULL;

    // start position of sequences ( first sequence is assumed to be 0 and
    // not recorded )
    boundaries[seq_idx - 1] = total_seq_size + seq->size;
    //printf("Storing in boundaries as: %ld\n", boundaries[hash->elCount - 1]);
    boundaries[seq_idx] = 0;

    // We want all uppercase
    toUpperN(seq->dna, seq->size);

    // sequence
    for (i = 0; i < seq->size; i++)
      if ((seq->dna)[i] == 'A')
        sequence[total_seq_size + i] = 0;
      else if ((seq->dna)[i] == 'C')
        sequence[total_seq_size + i] = 1;
      else if ((seq->dna)[i] == 'G')
        sequence[total_seq_size + i] = 2;
      else if ((seq->dna)[i] == 'T')
        sequence[total_seq_size + i] = 3;
      else
        sequence[total_seq_size + i] = 99;

    total_seq_size += seq->size;
    seq_idx++;
    dnaSeqFree(&seq);
  }
  twoBitClose(&tbf);

  printf("Read in %d seqs\n", seq_idx - 1);

  struct sequenceLibrary *seqLib = (struct sequenceLibrary *)malloc(sizeof(struct sequenceLibrary));
  seqLib->sequence = sequence;
  seqLib->identifiers = sequence_idents;
  seqLib->boundaries = boundaries;
  seqLib->length = total_seq_size;
  seqLib->count = seq_idx - 1;

  return seqLib;
}

/*
 * loadSequenceSubset()
 *
 *
 *
 *
 * TODO:
 *     The two bit spec doesn't need to contain the twoBitFilename over and over again
 *     ...just prepend before call to twoBitReadSeqFrag()
 *
 */
struct sequenceLibrary *
loadSequenceSubset(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *N)
{
  struct twoBitFile *tbf;
  struct twoBitSpec *tbs;
  struct hash *hash = newHash(6);
  uint64_t total_seq_size = 0;

  char *sequence = NULL;
  uint64_t *boundaries = NULL;
  char **sequence_idents = NULL;
  struct coreAlignment *curr_core_align = NULL;

  int i;
  // Must start from 1 otherwise it looks like the
  // NULL return value.
  int seq_idx = 1;
  int p_seq_idx = 0;
  int range_cnt = 0;

  tbs = readBEDRanges(twoBitName, rangeBEDName);
  if (tbs == NULL)
    errAbort("%s is not a twoBit file", twoBitName);

  tbf = twoBitOpen(tbs->fileName);
  struct twoBitSeqSpec *s;
  for (s = tbs->seqs; s != NULL; s = s->next)
  {
    // This is the full sequence per unique id approach
    if (hashLookup(hash, s->name) == NULL)
    {
      //printf("adding to hash %s as seq_idx = %d\n", s->name, seq_idx);
      hashAddInt(hash, s->name, seq_idx);

      // Complete sequence
      struct dnaSeq *seq = twoBitReadSeqFrag(tbf, s->name, 0, 0);
      //printf("  seq size = %d\n", seq->size);

      // We want all uppercase
      toUpperN(seq->dna, seq->size);

      sequence_idents =
        (char **) realloc(sequence_idents,
                          (hash->elCount + 1) * sizeof(char *));
      if (NULL == sequence_idents)
      {
        fprintf(stderr, "Could not allocated more space for the "
                "sequence_idents array\n");
        exit(1);
      }
      boundaries =
        (uint64_t *) realloc(boundaries,
                             (hash->elCount + 1) * sizeof(uint64_t));
      if (NULL == boundaries)
      {
        fprintf(stderr, "Could not allocated more space for the "
                "boundaries array\n");
        exit(1);
      }
      sequence =
        (char *) realloc(sequence,
                         (total_seq_size + seq->size + 1) * sizeof(char));
      if (NULL == sequence)
      {
        fprintf(stderr, "Could not allocated more space for the "
                "sequence array\n");
        exit(1);
      }

      // identifier
      sequence_idents[hash->elCount - 1] = cloneString(seq->name);
      sequence_idents[hash->elCount] = NULL;

      // start position of sequences ( first sequence is assumed to be 0 and
      // not recorded )
      boundaries[hash->elCount - 1] = total_seq_size + seq->size;
      //printf("Storing in boundaries as: %ld\n", boundaries[hash->elCount - 1]);
      boundaries[hash->elCount] = 0;

      // sequence
      for (i = 0; i < seq->size; i++)
        if ((seq->dna)[i] == 'A')
          sequence[total_seq_size + i] = 0;
        else if ((seq->dna)[i] == 'C')
          sequence[total_seq_size + i] = 1;
        else if ((seq->dna)[i] == 'G')
          sequence[total_seq_size + i] = 2;
        else if ((seq->dna)[i] == 'T')
          sequence[total_seq_size + i] = 3;
        else
          sequence[total_seq_size + i] = 99;

      total_seq_size += seq->size;
      seq_idx++;
      dnaSeqFree(&seq);
    }
    if ((p_seq_idx = hashIntVal(hash, s->name)) != 0)
    {
      range_cnt++;

      if (curr_core_align == NULL)
      {
        curr_core_align =
          (struct coreAlignment *) malloc(sizeof(struct coreAlignment));
        *core_align = curr_core_align;
      }
      else
      {
        curr_core_align->next =
          (struct coreAlignment *) malloc(sizeof(struct coreAlignment));
        curr_core_align = curr_core_align->next;
      }

      curr_core_align->next = NULL;
      curr_core_align->seqIdx = p_seq_idx - 1;
      // BED6:name = leftExtendable in our use
      if (atoi(s->bedName) == 1)
        curr_core_align->leftExtendable = 1;
      else
        curr_core_align->leftExtendable = 0;
      // BED6:score = rightExtendable in our use
      if (atoi(s->score) == 1)
        curr_core_align->rightExtendable = 1;
      else
        curr_core_align->rightExtendable = 0;
      curr_core_align->leftExtensionLen = 0;
      curr_core_align->rightExtensionLen = 0;

      uint64_t seq_offset = 0;
      if (p_seq_idx > 1)
        seq_offset = boundaries[p_seq_idx - 2];

      // Translate ascending ranges into logical
      // left/right positions for the coreAlignment
      // datastructure.
      if (strcmp(s->strand, "-") == 0)
      {
        curr_core_align->orient = 1;
        curr_core_align->leftSeqPos = s->end + seq_offset - 1;
        curr_core_align->rightSeqPos = s->start + seq_offset;
      }
      else
      {
        curr_core_align->orient = 0;
        curr_core_align->leftSeqPos = s->start + seq_offset;
        curr_core_align->rightSeqPos = s->end + seq_offset - 1;
      }

      //printCoreAlignmentStruct( curr_core_align );
      if (*core_align == NULL)
        *core_align = curr_core_align;
    }
  }
  twoBitSpecFree(&tbs);
  twoBitClose(&tbf);

  printf("Read in %d seqs\n", hash->elCount);
  *N = range_cnt;

  struct sequenceLibrary *seqLib = (struct sequenceLibrary *)malloc(sizeof(struct sequenceLibrary));
  seqLib->sequence = sequence;
  seqLib->identifiers = sequence_idents;
  seqLib->boundaries = boundaries;
  seqLib->length = total_seq_size;
  seqLib->count = hash->elCount;

  return seqLib;
}

/*
 * Parse a BED6 file containing a list of core alignment ranges.
 * Generate a list of twoBitSpecs using the BED ranges and the
 * name of twoBitFile they will be retrieved from.
 *
 *  BED6 Format
 *
 *   chrom     : sequence identifier
 *   chromStart: lower aligned position ( 0 based )
 *   chromEnd  : upper aligned position ( 0 based, half open )
 *   name      : left extendable flag ( 0 = no, 1 = yes )
 *   score     : right extendable flag
 *   strand    : strand ( '+' = forward, '-' = reverse )
 *
 * The fields are tab separated. Coordinates are zero-based half
 * open and the orientation is either "+" or "-".
 */
struct twoBitSpec *
readBEDRanges(char *twoBitFile, char *bedFile)
{
  struct lineFile *lf = lineFileOpen(bedFile, TRUE);
  char *fields[6];
  char *e;
  struct twoBitSpec *spec;
  struct twoBitSeqSpec *seq;
  AllocVar(spec);
  spec->fileName = cloneString(twoBitFile);
  while (lineFileChopCharNext(lf, '\t', fields, 6))
  {
    if ( !( fields[5][0] == '+' || fields[5][0] == '-' ) )
    {
      printf("Error: ranges file does not appear to be in the correct format!\n");
      exit(1);
    }

    AllocVar(seq);
    seq->name = cloneString(fields[0]);
    seq->start = strtol(fields[1], &e, 0);
    seq->end = strtol(fields[2], &e, 0);
    // Using BED6 to hold additional details
    //  BED6:name  = left extension flag ( 1 = extendable, 0 = not extendable )
    //  BED6:score = right extension flag ( 1 = extendable, 0 = not extendable )
    //  BED6:strand = orientation "+" or "-"
    seq->bedName = cloneString(fields[3]);
    seq->score = cloneString(fields[4]);
    seq->strand = cloneString(fields[5]);
    slSafeAddHead(&spec->seqs, seq);
  }
  slReverse(&spec->seqs);
  lineFileClose(&lf);
  return spec;
}


// Encoding from original RepeatScout codebase
char
char_to_num(char c)
{
  if (c == 'A')
    return 0;
  if (c == 'C')
    return 1;
  if (c == 'G')
    return 2;
  if (c == 'T')
    return 3;
  if (c == 'a')
    return 0;
  if (c == 'c')
    return 1;
  if (c == 'g')
    return 2;
  if (c == 't')
    return 3;
  if (c == 'N')
    return 99;
  if (IUPAC(c))
    return 99;
  if (c == 'n')
    return 99;
  if (c == 'x')
    return 99;
  printf("ERROR: Cannot interpret input symbol '%c' [%d] as a DNA base.\n", c, c);
  exit(1);
}

char
num_to_char(char z)
{
  if (z == 0)
    return (char) 'A';
  if (z == 1)
    return (char) 'C';
  if (z == 2)
    return (char) 'G';
  if (z == 3)
    return (char) 'T';
  return (char) 'N';
}

char
compl(char c)
{
  if (c == 0)
    return 3;
  if (c == 1)
    return 2;
  if (c == 2)
    return 1;
  if (c == 3)
    return 0;
  return 99;
}
