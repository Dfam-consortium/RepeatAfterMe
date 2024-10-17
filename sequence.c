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
#include <math.h>
#include "bnw_extend.h"
#include "sequence.h"
#include "kentsrc/common.h"
#include "kentsrc/hash.h"
#include "kentsrc/linefile.h"
#include "kentsrc/dnaseq.h"
#include "kentsrc/twoBitNew.h"

#ifdef AW_FM_INDEX
#include "AwFmIndex.h"
#endif


/*
 * POSSIBLY DEPRECATED: Use the Minimal loading routine now. This doesn't
 *                      update the new coreAlignment bounds/flags.
 * loadSequences()
 *
 * TODO: Add option to support softmasking rather than forcing everything to uppercase
 *
 *  markovChainOrder:   If > 0, calculate the counts of all words from 1-# and store
 *                      in sequenceLibrary->markovChainProbTables.
 */
struct sequenceLibrary *
loadSequences(char *twoBitName, char **softmasked, int markovChainOrder)
{
  struct twoBitFile *tbf;
  uint64_t total_seq_size = 0;

  char *sequence = NULL;
  uint64_t *boundaries = NULL;
  char **sequence_idents = NULL;

  int i, j, k, m;
  // Must start from 1 otherwise it looks like the
  // NULL return value.
  int seq_idx = 1;


  // Allocate space for Markov chain probabily tables
  uint32_t **markovChainProbTables = NULL;
  uint64_t markovSeqLen = 0;
  uint64_t markovSeqIndex = 0;
  if ( markovChainOrder > 0 ) {
    markovChainProbTables =
      (uint32_t **) malloc(markovChainOrder * sizeof(uint32_t *));
    if (!markovChainProbTables)
    {
      fprintf(stderr,
              "Unable to allocate memory for the markovChainProbTables!\n");
      exit(1);
    }

    for (i = 0; i < markovChainOrder; i++)
    {
      markovChainProbTables[i] =
        (uint32_t *) malloc(pow(4, (i + 1)) * sizeof(uint32_t));
      if (!markovChainProbTables[i])
      {
        fprintf(stderr,
                "Unable to allocate memory for the %d-order markovWordFreqTable !\n",
                i + 1);
        exit(1);
      }
      for (j = 0; j < pow(4, (i + 1)); j++)
        markovChainProbTables[i][j] = 0;
    }
  }

  tbf = twoBitOpen(twoBitName);
  struct twoBitIndex *index = NULL;
  for (index = tbf->indexList; index != NULL; index = index->next)
  {
    markovSeqLen = 0;
    markovSeqIndex = 0;

    // Complete sequence
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, index->name, 0, 0);

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

    // process the sequence
    for (i = 0; i < seq->size; i++) {
      if ((seq->dna)[i] == 'A')
        sequence[total_seq_size + i] = 0;
      else if ((seq->dna)[i] == 'C')
        sequence[total_seq_size + i] = 1;
      else if ((seq->dna)[i] == 'G')
        sequence[total_seq_size + i] = 2;
      else if ((seq->dna)[i] == 'T')
        sequence[total_seq_size + i] = 3;
      else {
        sequence[total_seq_size + i] = 99;
        // When calculating the Markov frequencies
        // an "N" is treated like the end of a sequence.
        // I.e don't include it in the markov chain and
        // instead treat it was a the start of a new one.
        markovSeqLen = 0;
        markovSeqIndex = 0;
      }
      if ( markovChainOrder > 0 ) {
        markovSeqIndex = (markovSeqIndex << 2) | seq->dna[i];

        // Update the Markov Chains
        for (k = 0; k < markovChainOrder; k++)
        {
          if (k <= markovSeqLen)
          {
            uint32_t wordIndex =
              markovSeqIndex & (uint64_t) (pow(4, k + 1) - 1);
            markovChainProbTables[k][wordIndex]++;
          }
        }
        markovSeqLen++;
      }
    }

    total_seq_size += seq->size;
    seq_idx++;
    dnaSeqFree(&seq);
  }
  twoBitClose(&tbf);

  struct sequenceLibrary *seqLib = (struct sequenceLibrary *)malloc(sizeof(struct sequenceLibrary));
  seqLib->sequence = sequence;
  seqLib->identifiers = sequence_idents;
  seqLib->boundaries = boundaries;
  seqLib->length = total_seq_size;
  seqLib->count = seq_idx - 1;

  // Store Markov chain analysis, if performed.
  if ( markovChainOrder ) {

    // Counts to probabilities
    for (k = 0; k < markovChainOrder; k++)
    {
      for (j = 0; j < pow(4, (k + 1)); j += 4)
      {
        // Calculate the sum
        double prefixSum = 0.0;
        for (m = j; m < j + 4; m++)
        {
          prefixSum += markovChainProbTables[k][m];
        }

        for (m = j; m < j + 4; m++)
        {
          markovChainProbTables[k][m] =
            MARKOV_PROB_SCALE_FACTOR * ((double) markovChainProbTables[k][m] / prefixSum);
        }
      }
    }

    seqLib->markov_chain_order = markovChainOrder;
    seqLib->markov_chain_prob_tables = markovChainProbTables;
  }else {
    seqLib->markov_chain_order = 0;
    seqLib->markov_chain_prob_tables = NULL;
  }

 return seqLib;
}

void seqlib_destroy(struct sequenceLibrary *seqLib)
{
  if ( seqLib == NULL )
    return;

  int i;
  if ( seqLib->markov_chain_prob_tables != NULL )
  {
    for (i = 0; i < seqLib->markov_chain_order; i++)
      free(seqLib->markov_chain_prob_tables[i]);
    free(seqLib->markov_chain_prob_tables);
  }
  if ( seqLib->sequence != NULL )
    free(seqLib->sequence);
  if ( seqLib->identifiers != NULL )
  {
    for (i = 0; i < seqLib->count; i++ )
      free(seqLib->identifiers[i]);
    free(seqLib->identifiers);
  }
  if ( seqLib->boundaries != NULL )
    free(seqLib->boundaries);
  free(seqLib);
}

/*
 * POSSIBLY DEPRECATED: Use the Minimal loading routine now. This doesn't
 *                      update the new coreAlignment bounds/flags.
 * loadSequenceSubset()
 *
 *   This routine parses the range BED file and loads in one copy of every
 *   reference sequence.  This is a simple approach but if there are large
 *   sequences in the assembly and a large number of cores it will be likely
 *   that this will have to read in the entire assembly into memory.
 *
 * TODO:
 *     The two bit spec doesn't need to contain the twoBitFilename over and over again
 *     ...just prepend before call to twoBitReadSeqFrag()
 *
 *     Add option to support softmasking rather than forcing everything to uppercase.
 */
struct sequenceLibrary *
loadSequenceSubset(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *num_cores)
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
      curr_core_align->score = 0;

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
  *num_cores = range_cnt;

  struct sequenceLibrary *seqLib = (struct sequenceLibrary *)malloc(sizeof(struct sequenceLibrary));
  seqLib->sequence = sequence;
  seqLib->identifiers = sequence_idents;
  seqLib->boundaries = boundaries;
  seqLib->offsets = NULL;
  seqLib->length = total_seq_size;
  seqLib->count = hash->elCount;

  return seqLib;
}

/*
 * loadSequenceSubsetMinimal()
 *
 *   This routine parses the range BED file and loads in subsequences, one per
 *   core with flanking sequences.  This provides the option to specify the amount
 *   of flanking sequence to load.  NOTE: It is possible that duplicate flanking
 *   sequences are loaded for multiple cores, although in practice this doesn't
 *   account for much overhead due to the TANDEMDIST requirement.
 *
 *   NOTE: This is designed for BED files with low-high sequence coordinates
 *         and relies on the orient field to define which is the left/right
 *         edge.
 *
 *   TODO:
 *     - Add support for extension overlap detection ( tandem repeats )
 *       In RAMDenovo short tandem repeats are avoided using the TANDEMDIST
 *       parameter.  This is a simple way to avoid picking words near each
 *       other.   In RAMExtend I would to avoid *and* detect these patterns
 *       at the same time using the coordinates of each core sequence *and*
 *       the extension limits from the first side extension to limit the
 *       second side extension.
 *
 *       So when loading sequences, we can either limit loading the left
 *       flank, if there is another core sequence <10,000bp upstream.
 *       Or we could record this position in a new left-limit field.
 *
 *       The right sequence needs to be loaded, regardless, as it will be
 *       limited by the left-extension limits.
 *       To determine the best approach, I need to see what the easiest way
 *       will be to catch this in the dynamic programming loops.
 *
 */

int bedNameStartRevEndCmp(const void *va, const void *vb)
{
  const struct twoBitSeqSpec *a = *((struct twoBitSeqSpec **)va);
  const struct twoBitSeqSpec *b = *((struct twoBitSeqSpec **)vb);
  int diff = strcmp(a->name, b->name);
  if (diff == 0)
    diff = a->start - b->start;
  if (diff == 0)
    diff = b->end - a->end;
  return diff;
}

struct sequenceLibrary *
loadSequenceSubsetMinimal(char *twoBitName, char *rangeBEDName,
                   struct coreAlignment **core_align,
                   int *num_cores, int max_flanking_bp)
{
  struct twoBitFile *tbf;
  struct twoBitSpec *tbs;
  uint64_t total_seq_size = 0;

  char *sequence = NULL;
  uint64_t *boundaries = NULL;
  uint64_t *offsets = NULL;
  char **sequence_idents = NULL;
  struct coreAlignment *curr_core_align = NULL;

  int i;
  // Must start from 1 otherwise it looks like the
  // NULL return value.
  int seq_idx = 1;
  int range_cnt = 0;

  tbs = readBEDRanges(twoBitName, rangeBEDName);
  if (tbs == NULL)
    errAbort("%s is not a twoBit file", twoBitName);

  tbf = twoBitOpen(tbs->fileName);
  struct twoBitSeqSpec *s;

  // Can't depend on the BED file being sorted...let's
  // sort it out ourselves.
  slSort(&(tbs->seqs), bedNameStartRevEndCmp );

  struct twoBitSeqSpec *prev_seq = NULL;
  for (s = tbs->seqs; s != NULL; s = s->next)
  {
    // Determine sequence size
    int seq_size = twoBitSeqSize(tbf, s->name);

    //printf("BED: Core sequence %s:%d-%d\n", s->name, s->start, s->end);
    //printf("  - twoBit size for %s is %d\n", s->name, seq_size);

    int lower_flank_len = 0;
    int upper_flank_len = 0;
    int flanking_start = s->start;
    int flanking_end = s->end;
    enum CoreBoundFlag lower_bound_flag = L_BOUNDARY;
    enum CoreBoundFlag upper_bound_flag = L_BOUNDARY;

    if (strcmp(s->strand, "-") == 0)
    {
      // Reverse strand
      // left = end, right = start

      // BED6:name = leftExtendable in our use
      if ( atoi(s->bedName) == 1 )
      {
        if ( s->end + max_flanking_bp < seq_size )
        {
          flanking_end = s->end + max_flanking_bp;
          upper_flank_len = max_flanking_bp;
          upper_bound_flag =  L_BOUNDARY;  // max_flanking (L) limited boundary
        }else {
          flanking_end = seq_size;
          upper_flank_len = seq_size - s->end;
          upper_bound_flag =  SEQ_BOUNDARY; // sequence limited boundary
        }
      }
      // BED6:score = rightExtendable in our use
      if ( atoi(s->score) == 1 )
      {
        if ( max_flanking_bp < s->start )
        {
          flanking_start = s->start - max_flanking_bp;
          lower_flank_len = max_flanking_bp;
          lower_bound_flag = L_BOUNDARY;  // max_flanking (L) limited boundary
        }else {
          flanking_start = 0;
          lower_flank_len = s->start;
          lower_bound_flag = SEQ_BOUNDARY; // sequence limited boundary
        }
      }
    }else {
      // Forward strand
      // left = start, right = end

      // BED6:name = leftExtendable in our use
      if ( atoi(s->bedName) == 1 )
      {
        if ( max_flanking_bp < s->start )
        {
          flanking_start = s->start - max_flanking_bp;
          lower_flank_len = max_flanking_bp;
          lower_bound_flag = L_BOUNDARY;  // max_flanking (L) limited boundary
        }else
        {
          flanking_start = 0;
          lower_flank_len = s->start;
          lower_bound_flag = SEQ_BOUNDARY;  // sequence limited boundary
        }
      }
      // BED6:score = rightExtendable in our use
      if ( atoi(s->score) == 1 )
      {
        if ( s->end+max_flanking_bp < seq_size )
        {
          flanking_end = s->end + max_flanking_bp;
          upper_flank_len = max_flanking_bp;
          upper_bound_flag = L_BOUNDARY;  // max_flanking (L) limited boundary
        }else
        {
          flanking_end = seq_size;
          upper_flank_len = seq_size - s->end;
          upper_bound_flag = SEQ_BOUNDARY; // sequence limited boundary
        }
      }
    }

    // Extract subsequence + flanking:
    //      |--- <= max_flanking_bp ----^^^^^^^core^^^^^^^^--- <= max_flanking_bp ----|
    //   Coordinates are zero-based, half-open
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, s->name, flanking_start, flanking_end);

    //printf("  - Extraction range %d-%d (max_flanking_bp = %d) size should be = %d\n",
    //       flanking_start, flanking_end, max_flanking_bp, flanking_end-flanking_start);
    //printf("  seq size = %d\n", seq->size);

    // We want all uppercase
    toUpperN(seq->dna, seq->size);

    sequence_idents =
      (char **) realloc(sequence_idents,
                        (seq_idx+1) * sizeof(char *));
    if (NULL == sequence_idents)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "sequence_idents array\n");
      exit(1);
    }
    boundaries =
      (uint64_t *) realloc(boundaries,
                           (seq_idx+1) * sizeof(uint64_t));
    if (NULL == boundaries)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "boundaries array\n");
      exit(1);
    }
    offsets =
      (uint64_t *) realloc(offsets,
                           (seq_idx+1) * sizeof(uint64_t));
    if (NULL == offsets)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "offsets array\n");
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

    // Store identifier of each sequence
    sequence_idents[seq_idx - 1] = cloneString(s->name);
    sequence_idents[seq_idx] = NULL;

    // start position of sequences ( first sequence is assumed to be 0 and
    // not recorded )
    boundaries[seq_idx-1] = total_seq_size + seq->size;
    //printf("Storing in boundaries as: %ld\n", boundaries[hash->elCount - 1]);
    boundaries[seq_idx] = 0;

    offsets[seq_idx-1] = flanking_start;
    offsets[seq_idx] = 0;

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
    curr_core_align->seqIdx = seq_idx-1;

    //
    // Set the core alignment bounds
    //
    curr_core_align->lowerSeqBound = total_seq_size - seq->size;
    curr_core_align->upperSeqBound = total_seq_size - 1;
    curr_core_align->lowerSeqBoundFlag = lower_bound_flag;
    curr_core_align->upperSeqBoundFlag = upper_bound_flag;
    //printf("Core %d: before: %d-%d ", range_cnt, curr_core_align->lowerSeqBound, curr_core_align->upperSeqBound);
    // Look for core conflicts
    if ( prev_seq != NULL && strcmp(s->name, prev_seq->name) == 0 )
    {
      if ( s->start - prev_seq->end <= lower_flank_len )
      {
        //printf("BED: Found potential lower bound conflict!\n");
        curr_core_align->lowerSeqBound = total_seq_size - seq->size + ( lower_flank_len - (s->start - prev_seq->end) );
        curr_core_align->lowerSeqBoundFlag = CORE_BOUNDARY; // core limited boundary
      }
    }
    if ( s->next != NULL && strcmp(s->name, s->next->name) == 0 )
    {
      if ( s->end + upper_flank_len >= s->next->start )
      {
        //printf("BED: Found potential upper bound conflict! %ld and next is at %ld diff %ld ... %ld upper_flank_len = %ld\n", s->end, s->next->start, s->next->start - s->end, (upper_flank_len - (s->next->start - s->end - 1)), upper_flank_len);
        //printf("BED:     total_seq_size = %ld, ( upper_flank_len - (s->next->start - s->end - 1) ) = %ld\n", total_seq_size, ( upper_flank_len - (s->next->start - s->end - 1) ));
        curr_core_align->upperSeqBound = total_seq_size - ( upper_flank_len - (s->next->start - s->end ) ) - 1;
        curr_core_align->upperSeqBoundFlag = CORE_BOUNDARY;
      }
    }
    //printf("CORE: %d %ld [%d] - %ld [%d]\n", range_cnt, curr_core_align->lowerSeqBound, curr_core_align->lowerSeqBoundFlag, curr_core_align->upperSeqBound, curr_core_align->upperSeqBoundFlag);
    prev_seq = s;

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
    curr_core_align->score = 0;

    uint64_t seq_offset = 0;
    if (seq_idx > 1)
      seq_offset = boundaries[seq_idx - 2];

    // Translate ascending ranges into logical
    // left/right order for the coreAlignment
    // datastructure.
    //  E.g. The TSV file was ordered like a BED file would be:
    //
    //         seq1  10  20  +
    //         seq2  30  40  -
    //         seq2  50  100 -
    //
    //  The coreAlignment struct should be ordered like:
    //
    //      core1:  leftSeqPos=10, rightSeqPos=20
    //      core2:  leftSeqPos=40, rightSeqPos=30
    //      core3:  leftSeqPos=100, rightSeqPos=50
    //
    if (strcmp(s->strand, "-") == 0)
    {
      curr_core_align->orient = 1;
      curr_core_align->rightSeqPos = seq_offset + lower_flank_len;
      curr_core_align->leftSeqPos = curr_core_align->rightSeqPos + (s->end - s->start) - 1;
    }
    else
    {
      curr_core_align->orient = 0;
      curr_core_align->leftSeqPos = seq_offset + lower_flank_len;
      curr_core_align->rightSeqPos = curr_core_align->leftSeqPos + (s->end - s->start) - 1;
    }

    //printCoreAlignmentStruct( curr_core_align );
    if (*core_align == NULL)
      *core_align = curr_core_align;

    seq_idx++;
    dnaSeqFree(&seq);
  }
  twoBitSpecFree(&tbs);
  twoBitClose(&tbf);

  //printf("Read in %d seqs\n", range_cnt);
  *num_cores = range_cnt;

  struct sequenceLibrary *seqLib = (struct sequenceLibrary *)malloc(sizeof(struct sequenceLibrary));
  seqLib->sequence = sequence;
  seqLib->identifiers = sequence_idents;
  seqLib->boundaries = boundaries;
  seqLib->offsets = offsets;
  seqLib->length = total_seq_size;
  seqLib->count = range_cnt;

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



#ifdef AW_FM_INDEX
/*
 * createFMIndex()
 *
 *   Build an FMIndex given a twoBitFile.
 *
 * TODO: Add option to support softmasking from twobit file.
 * TODO: Would this be faster if we set the kmerLengthInSeedTable to l?
 */
struct AwFmIndex *
createFMIndex(char *twoBitName)
{
  struct twoBitFile *tbf;
  char *sequence = NULL;
  struct AwFmIndex *fm_index = NULL;
  uint64_t total_seq_size = 0;

  tbf = twoBitOpen(twoBitName);
  struct twoBitIndex *index;
  for (index = tbf->indexList; index != NULL; index = index->next)
  {
    // Complete sequence
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, index->name, 0, 0);
    toUpperN(seq->dna, seq->size);

    sequence =
      (char *) realloc(sequence,
                       (total_seq_size + seq->size + 1) * sizeof(char));
    if (NULL == sequence)
    {
      fprintf(stderr, "Could not allocated more space for the "
              "sequence array\n");
      exit(1);
    }
    strncpy(&(sequence[total_seq_size]), seq->dna, seq->size);
    sequence[total_seq_size + seq->size] = '\0';

    total_seq_size += seq->size;
    dnaSeqFree(&seq);
  }
  twoBitClose(&tbf);

  printf("createFMIndex(): Read in %ld bp\n", total_seq_size);

  struct AwFmIndexConfiguration fm_config = { .suffixArrayCompressionRatio = 8,
                                              .kmerLengthInSeedTable = 12,
                                              // EXPERIMENTAL
                                              .keepSuffixArrayInMemory = TRUE,
                                              .alphabetType =AwFmAlphabetNucleotide };

  enum AwFmReturnCode creationReturnCode = awFmCreateIndex(&fm_index, &fm_config, (uint8_t *)sequence,
                                                           (size_t)total_seq_size, "testIndex.awfmi", 1);

  if ( creationReturnCode != AwFmFileWriteOkay ) {
    if ( creationReturnCode == AwFmNullPtrError )
      printf("createFMIndex(): awFmCreateIndex returned AwFmNullPtrError\n");
    else if ( creationReturnCode == AwFmAllocationFailure )
      printf("createFMIndex(): awFmCreateIndex returned AwFmAllocationFailure\n");
    else if ( creationReturnCode == AwFmFileAlreadyExists )
      printf("createFMIndex(): awFmCreateIndex returned AwFmFileAlreadyExists\n");
    else if ( creationReturnCode == AwFmSuffixArrayCreationFailure )
      printf("createFMIndex(): awFmCreateIndex returned AwFmSuffixArrayCreationFailure\n");
    else if ( creationReturnCode == AwFmFileWriteFail )
      printf("createFMIndex(): awFmCreateIndex returned AwFmFileWriteFail\n");
    else
      printf("createFMIndex(): awFmCreateIndex returned an error code: %d\n", creationReturnCode);
  }

  free(sequence);
  return fm_index;
}
#endif

// Encoding similar to the original RepeatScout codebase except that
// we no longer treat lowercase as equivalent to uppercase nucleotides.
// This allows us to use a softmasking approach by setting the matrix
// values for {'a','c','g','t'} as equivalent to "N".
char
char_to_num(char c)
{
  if (c == 'A')
    return SYM_A;
  if (c == 'C')
    return SYM_C;
  if (c == 'G')
    return SYM_G;
  if (c == 'T')
    return SYM_T;
  if (c == 'a')
    return SYM_a;
  if (c == 'c')
    return SYM_c;
  if (c == 'g')
    return SYM_g;
  if (c == 't')
    return SYM_t;
  if (c == 'N')
    return SYM_N;
  if (IUPAC(c))
    return SYM_N;
  if (c == 'n')
    return SYM_N;
  if (c == 'x')
    return SYM_N;
  printf("ERROR: Cannot interpret input symbol '%c' [%d] as a DNA base.\n", c, c);
  exit(1);
}


// Convert sequence encoding back to ASCII characters
//   - This now handles soft-mask encoding
char
num_to_char(char z)
{
  if (z == SYM_A)
    return (char) 'A';
  if (z == SYM_C)
    return (char) 'C';
  if (z == SYM_G)
    return (char) 'G';
  if (z == SYM_T)
    return (char) 'T';
  if (z == SYM_a)
    return (char) 'a';
  if (z == SYM_c)
    return (char) 'c';
  if (z == SYM_g)
    return (char) 'g';
  if (z == SYM_t)
    return (char) 't';
  return (char) 'N';
}

// Same as above but complemented
char
num_to_char_compl(char z)
{
  if (z == SYM_A)
    return (char) 'T';
  if (z == SYM_C)
    return (char) 'G';
  if (z == SYM_G)
    return (char) 'C';
  if (z == SYM_T)
    return (char) 'A';
  if (z == SYM_a)
    return (char) 't';
  if (z == SYM_c)
    return (char) 'g';
  if (z == SYM_g)
    return (char) 'c';
  if (z == SYM_t)
    return (char) 'a';
  return (char) 'N';
}



// Complement encoded nucleotide
//   - This now handles soft-mask encoding
char
compl(char c)
{
  if (c == SYM_A) // A->T
    return SYM_T;
  if (c == SYM_C) // C->G
    return SYM_G;
  if (c == SYM_G) // G->C
    return SYM_C;
  if (c == SYM_T) // T->A
    return SYM_A;
  if (c == SYM_a) // a->t
    return SYM_t;
  if (c == SYM_c) // c->g
    return SYM_g;
  if (c == SYM_g) // g->c
    return SYM_c;
  if (c == SYM_t) // t->a
    return SYM_a;
  return SYM_N;
}

// Return uppercase version (unmasked) of an encoded base
char
encoded_toupper(char c)
{
  if (c == SYM_a)
    return SYM_A;
  if (c == SYM_c)
    return SYM_C;
  if (c == SYM_g)
    return SYM_G;
  if (c == SYM_t)
    return SYM_T;
  return c;
}

// Mask encoded nucleotide
char
mask(char c)
{
  if (c == SYM_A)
    return SYM_a;
  if (c == SYM_C)
    return SYM_c;
  if (c == SYM_G)
    return SYM_g;
  if (c == SYM_T)
    return SYM_t;
  return c;
}

int
kmermatch(char *kmer1, char *kmer2, int l)     /* forward match */
{
  int x;

  for (x = 0; x < l; x++)
  {
    if ((kmer1[x]&3) != (kmer2[x]&3))
      return 0;
  }
  return 1;
}


int
kmermatchrc(char *kmer1, char *kmer2, int mismatches, int l)   /* rc match */
{
  int x;

  int rem_mismatches = mismatches;
  for (x = 0; x < l; x++)
  {
    // A=0 + T=3  == 3
    // C=1 + G=2  == 3
    //  ...
    if (((kmer1[x]&3) + (kmer2[l - 1 - x]&3) != 3) && rem_mismatches-- == 0)
      return 0;
  }
  return 1;
}

// Return 1 if forward match, 0 if no match and -1 if reverse match
// Operate on encoded strings, case agnostic
int
kmermatcheither(char *kmer1, char *kmer2, int mismatches, int l)       /* forward or rc match */
{
  int x;

  int rem_mismatches = mismatches;
  for (x = 0; x < l; x++)
  {
    if ((kmer1[x]&3) != (kmer2[x]&3) && rem_mismatches-- == 0)
      return kmermatchrc(kmer1, kmer2, mismatches, l);
  }
  return 1;
}

