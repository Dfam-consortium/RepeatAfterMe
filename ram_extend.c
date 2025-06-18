/*
 * RAMExtend
 *
 * Use a RepeatScout-like method to simultaneously extend
 * a multiple alignment anchored by a set of genomic ranges.
 *
 * Authors:
 *    Robert Hubley
 *    Based on RepeatScout code by Alkes Price 2005
 *
 */
#include <stdio.h>
#include <sys/types.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "kentsrc/common.h"
#include "kentsrc/hash.h"
#include "kentsrc/linefile.h"
#include "kentsrc/dnaseq.h"
#include "kentsrc/twoBitNew.h"
#include "cmd_line_opts.h"
#include "sequence.h"
#include "score_system.h"
#include "bnw_extend.h"
#include "report.h"
#include "ram_extend.h"


// From Version.c created by the Makefile
extern const char *Version;
extern const char *BuildNumber;
extern const char *BuildDate;

uint64_t seqLength;             // Length of input sequence.
int l = 1;                      // Historical...this is the length of
                                // the region between left/right that
                                // will be kept in master.
char *SEQUENCE_FILE = NULL;     // Where to get sequence data from.
char *MATRIX_FILE = NULL;       // Eventually a file to read for a custom
                                //   matrix.  For now this can be used to
                                //   specify an internally hardcoded matrix.
char *OUTTSV_FILE = NULL;       // Where to store the table of extended ranges
char *OUTFA_FILE = NULL;        // Where to store the final extended sequences
char *CONS_FILE = NULL;         // Where to put consensus output
char *RANGES_FILE = NULL;       // Where to store the extended ranges
char *OUTMAT_FILE = NULL;       // Where to dump the DP matrix for debugging
int VERBOSE;                    // How chatty? Real chatty? Really really
                                // chatty?
                                // Super extra really chatty?
int L;                          // The maximum distance to attempt extension of each
                                //   core range. Length of master is 2*L+l (L >>> l)
int flanking = 0;               // The number of bp to add to the left/right side
                                //   of each -outfa sequence.
int BANDWIDTH;                  // max offset (14)
int MAXN;                       // max #occ of kmer (10000)
int WHEN_TO_STOP;               // stop if no improvement after extending
                                // this far (100)
float MINENTROPY;               // ignore freq kmers with entropy less than
                                // (0.70)

// Parameters dependent on matrix used
int GAP_EXT;                    // gap extension score
int GAP_OPEN;                   // gap initiation score
int CAPPENALTY;                 // cap on penalty for exiting alignment
                                // Similar in concept to X-Drop
                                // Made this 20*avg match score?
int MINIMPROVEMENT;             // amount totalbestscore must improve for
                                // each extra letter
                                // In Alkes paper this is the "c"
                                // value used in the fit-preferred
                                // scoring system. This should be
                                // a score that represents the minimum
                                // number of positively scoring
                                // sequences against the extended consensus.
                                //
int MATCH;                      // Only relevant to simple scoring system
int MISMATCH;                   // Only relatvie to simple scoring system

char *master;                   // Fancy name for current consensus sequence
int masterstart;
int masterend;


// Structure to hold the core alignment ranges for which we
// are going to extend.
struct coreAlignment *coreAlign = NULL;

// Structure to hold the sequence library
struct sequenceLibrary *seqLib = NULL;

// Structure to hold the matrix
struct scoringSystem *scoreParams = NULL;

int N;

//
// Changes to support matrix scoring:
//      The score array stores two scores in every cell rather than one.
//      The first score (0) represents the Substitution stat
//      score.  The second (1) score is the Gap state score.  NOTE: It is not
//      worth breaking the score down into del/ins states.  The assumption is
//      made that we will never count an insertion followed by deletion as
//      a valid transiton.
//
int ****score;                  // 2 * MAXN * (2*BANDWIDTH+1) * 2

//
//     Reduce memory use by freeing once no longer needed
//
//TODO: Move this into the core range structure
//int *overall_sequence_high_score;             // best alignment score for this n for any y
                                // or w seen so far
// Unused
int num_threads = 0;            // Number of threads or 0 if not
                                // multithreaded


void
usage()
{
  printf(
          "RAMExtend Version %s - build %s date %s\n"
          "\n"
          "  Perform a multiple sequence alignment (MSA) extension given\n"
          "  an existing core MSA and the flanking sequences.  The core\n"
          "  may be a set of word matches, a previously development MSA,\n"
          "  or the the result of any process that defines a core set of\n"
          "  sequence relationships.  The only requirement is that the sequences are\n"
          "  aligned to nearly the same point along one or both of the core\n"
          "  edges.  The extension process is a form of anchored alignment\n"
          "  where the details of the core alignment are considered fixed.\n"
          "\n"
          "Usage: \n"
          "  RAMExtend -twobit <seq.2bit> -ranges <ranges.tsv> [opts]\n"
          "\n"
          "--------------------------------------------------------------------------------------------\n"
          "     -cons <seq.fa>        # Save the left/right consensus sequences to a FASTA file.\n"
          "     -outtsv <final.tsv>   # Save the final sequence ranges to a TSV file.\n"
          "     -outfa <seq.fa>       # Save all the final sequences ( original range + extension )\n"
          "                           #   to a file.\n"
          "     -addflanking <num>    # Include an additional<num> bp of sequence to each sequence\n"
          "                           #   when using the -outfa option.\n"
          "     -L <num>              # Size of region to extend left or right (10000). \n"
          "     -minimprovement <num> # Amount that a the alignment score needs to improve each step\n"
          "                           #   to be considered progress (original: 3, nucleotide matrix: 27).\n"
          "     -maxoccurrences <num> # Cap on the number of sequences to align (10,000). \n"
          "     -cappenalty <num>     # Cap on penalty for exiting alignment of a sequence \n"
          "                           #   (original: -20, nucleotide matrix:-90).\n"
          "     -stopafter <num>      # Stop the alignment after this number of no-progress columns (100).\n"
          "     -minlength <num>      # Minimum required length for a sequence to be reported (50).\n"
          "     -outmat <file>        # Dump the dp matrix paths to a file for debugging.\n"
          "     -v[v[v[v]]]           # How verbose do you want it to be?  -vvvv is super-verbose.\n"
          "\n"
          "   New Scoring System Options:\n"
          "     -matrix 14p43g |      # There are several internally-coded matrices available.\n"
          "             18p43g |      #   The DNA substitution matrices are borrowed from RepeatMasker\n"
          "             20p43g |      #   and are tuned for a 43%% GC background.  The four encoded\n"
          "             25p43g        #   matrices are are divergence tuned, ranging from 14%% to 25%%.\n"
          "                           #   For example the '14p43g' matrix is tuned for 14%% divergence\n"
          "                           #   43%% GC background.  Each matrix has it's own default gap_open,\n"
          "                           #   gap_extension parameters.  These may be overrided by using the\n"
          "                           #   -gapopen and -gapextn options.\n"
          "     -gapopen <num>        # Affine gap open penalty to override per-matrix default\n"
          "     -gapextn <num>        # Affine gap extension penalty to override per-matrix default\n"
          "\n"
          "     -bandwidth <num>      # The maximum number of unbalanced gaps allowed (14)\n"
          "                           #   Half the bandwidth of the banded Smith-Waterman\n"
          "                           #   algorithm.  The default allows for at most 14bp\n"
          "                           #   of unbalanced insertions/deletions in any aligned\n"
          "                           #   sequence.  The full bandwidth is 2*bandwidth+1.\n"
          " or\n"
          "   Original RepeatScout Scoring System:\n"
          "     -matrix repeatscout   # The original RepeatScout scoring system is enabled if this\n"
          "                           #   option is set.  This uses a simple match/mismatch penalty\n"
          "                           #   and a linear gap model.  The original RepeatScout values\n"
          "                           #   for these parameters may be overriden with the following\n"
          "                           #   additional options:\n"
          "     -match <num>          # If '-matrix repeatscout' used, apply this reward for match\n"
          "                           #   (default +1)  \n"
          "     -mismatch <num>       # If '-matrix repeatscout' used, apply this penalty for a mismatch\n"
          "                           #   (default: -1) \n"
          "     -gap <num>            # If '-matrix repeatscout' used, apply this penalty for a gap\n"
          "                           #   (default: -5)\n"
          "\n"
          "Ranges:\n\n"
          "   Ranges are supplied in the form of a modified BED-6 format:\n\n"
          "      field-1:chrom     : sequence identifier\n"
          "      field-2:chromStart: lower aligned position ( 0 based )\n"
          "      field-3:chromEnd  : upper aligned position ( 0 based, half open )\n"
          "      field-4:left-ext  : left extendable flag ( 0 = no, 1 = yes ) [BED-6 'name' field]\n"
          "      field-5:right-ext : right extendable flag ( 0 = no, 1 = yes ) [BED-6 'score' field]\n"
          "      field-6:strand    : strand ( '+' = forward, '-' = reverse )\n"
          "\n"
          "   The fields are tab separated. Coordinates are zero-based half\n"
          "   open.\n"
          "\n"
          "   Left/Right extension flags are used to turn off extension for sequences that\n"
          "   do not reach the edge of the core alignment (e.g. fragments aligning in the\n"
          "   center of the core MSA, or one or the other edge only).  These sequences are\n"
          "   less likely to include related flanking sequence and weigh down the extension\n"
          "   score uneccesarily.  For these sequences it is desirable to flag one edge or\n"
          "   the other as 'unextendable'.  The 'left'/'right' designations refer to the\n"
          "   edges of the core MSA alignment. The 'left' flag may refer to the start or\n"
          "   the end coordinate depending on the state of the orientation flag\n"
          "   ('+' left=start, '-' left=end etc.).\n"
          , Version, BuildNumber, BuildDate);
  exit(1);

}

int
main(int argc, char *argv[])
{
  time_t start, finish;
  double duration;
  uint64_t x;
  int opt;
  FILE *fp = NULL;
  FILE *fp_fa = NULL;
  FILE *fp_mat = NULL;
  start = time(0);

  if (co_get_string(argc, argv, "-ranges", &RANGES_FILE) == 0)
  {
    usage();
    exit(1);
  }

  co_get_int(argc, argv, "-addflanking", &flanking);
  co_get_string(argc, argv, "-outtsv", &OUTTSV_FILE);
  co_get_string(argc, argv, "-outfa", &OUTFA_FILE);
  co_get_string(argc, argv, "-outmat", &OUTMAT_FILE);
  co_get_string(argc, argv, "-cons", &CONS_FILE);
  co_get_int(argc, argv, "-L", &L) || (L = 10000);
  co_get_int(argc, argv, "-bandwidth", &BANDWIDTH) || (BANDWIDTH = 14);
  co_get_int(argc, argv, "-maxoccurrences", &MAXN) || (MAXN = 10000);
  co_get_int(argc, argv, "-stopafter", &WHEN_TO_STOP) || (WHEN_TO_STOP = 100);
  if (!co_get_int(argc, argv, "-threads", &num_threads))
    num_threads = 0;

  VERBOSE =
    co_get_bool(argc, argv, "-vvvvvv", &opt) ? 20 :
    co_get_bool(argc, argv, "-vvvvv", &opt) ? 12 :
    co_get_bool(argc, argv, "-vvvv", &opt) ? 10 :
    co_get_bool(argc, argv, "-vvv", &opt) ? 3 :
    co_get_bool(argc, argv, "-vv", &opt) ? 2 :
    co_get_bool(argc, argv, "-v", &opt) ? 1 : 0;

  //
  // RMH: Build simple or complex matrices
  //
  /*
  for (i = 0; i < 256; i++)
    alpha_convert[i] = -1;

  for (i = 0; i < alphsize; i++)
  {
    alpha_convert[(int) alphabet[i]] = i;
    alpha_convert[(int) lc_alphabet[i]] = i;
  }
  // Hack to be equivalent to RS numbering.  Revisit
  // this later.
  alpha_convert['N'] = 99;
  alpha_convert['n'] = 99;
  */

  // Default matrix
  if (!co_get_string(argc, argv, "-matrix", &MATRIX_FILE))
    MATRIX_FILE = "20p43g";

  if (strcmp(MATRIX_FILE, "repeatscout") == 0)
  {
    if ( co_get_int(argc, argv, "-match", &MATCH) &&
         co_get_int(argc, argv, "-mismatch", &MISMATCH) &&
         co_get_int(argc, argv, "-gap", &GAP_EXT) )
      scoreParams = getRepeatScoutMatrix(MATCH,MISMATCH,GAP_EXT);
    else
      scoreParams = getRepeatScoutMatrix(1,-1,-5);
  }else {
    if ( co_get_int(argc, argv, "-gapopen", &GAP_OPEN) &&
         co_get_int(argc, argv, "-gapext", &GAP_EXT) )
    {
      scoreParams = getMatrixUsingGapPenalties(MATRIX_FILE,GAP_OPEN,GAP_EXT);
    }else {
      scoreParams = getMatrix(MATRIX_FILE);
    }
  }

  if (strcmp(MATRIX_FILE, "14p43g") == 0)
  {
    //
    // The minimum improvement of the multiple alignment score per
    // position.
    // - If at a minimum we want 3 sequences positively extending
    //   then using the current matrix ( minimum positive value = 9 )
    //   the MINIMPROVEMENT should be 27.
    //
    // NOTES:
    //   The above may be too restrictive....
    //
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 27);
    //
    // The lowest score a poorly performing alignment may reach before
    // being capped.
    // - Max substitution penalty in our hardcoded matrix = -17
    //   Equivalent to Alkes original scoring schema ( -20 = 20 positions )
    //   would be 20*-17 = -340
    // - In practice this was too high.  Setting current default to
    //   -90.
    //
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -90);
  }else if (strcmp(MATRIX_FILE, "18p43g") == 0)
  {
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 27);
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -90);
  }else if (strcmp(MATRIX_FILE, "20p43g") == 0)
  {
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 27);
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -90);
  }else if (strcmp(MATRIX_FILE, "25p43g") == 0)
  {
    // TODO: here we diverge...why?
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 24);
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -90);
  }else if (strcmp(MATRIX_FILE, "repeatscout") == 0)
  {
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 3);
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -20);
  }else {
    printf("Matrix name not found!\n");
    exit(1);
  }

  // Consensus sequence which is kmer + left and right extension max (L)
  master = (char *) malloc((2 * L + l + 1) * sizeof(char));
  if (NULL == master)
  {
    fprintf(stderr, "Could not allocate space for master array\n");
    exit(1);
  }
  master[2 * L + l] = '\0';

  //overall_sequence_high_score = (int *) malloc(MAXN * sizeof(int));
  //if (NULL == overall_sequence_high_score)
  //{
  //  fprintf(stderr, "Could not allocated space for internal arrays\n");
  //  exit(1);
  //}

  if (co_get_string(argc, argv, "-twobit", &SEQUENCE_FILE))
  {
    dnaUtilOpen();
    // Load the cores and flanking sequence sufficient to
    // support extension by L+bandwidth ( assuming at most
    // bandwidth insertions that would be tolerated before the
    // alignment is terminated).
    seqLib = loadSequenceSubsetMinimal(SEQUENCE_FILE, RANGES_FILE,
                                        &coreAlign, &N, L+BANDWIDTH );

  }
  else if (co_get_string(argc, argv, "-sequence", &SEQUENCE_FILE))
  {
    printf("-sequence is deprecated!....may return someday\n");
    exit(1);
  }
  else
  {
    usage();
    exit(1);
  }

  if ( OUTMAT_FILE != NULL )
  {
    if ((fp_mat = fopen(OUTMAT_FILE, "w")) == NULL)
    {
      fprintf(stderr, "Could not create the output matrix file %s\n", OUTMAT_FILE);
      exit(1);
    }
  }


  //
  // Introduce ourselves
  //
  printf( "\nRAMExtend Version %s - build %s date %s\n"
          , Version, BuildNumber, BuildDate);
  print_parameters();
  printf("Read in %d ranges, and %ld bp of sequence\n\n", N, seqLib->length);

  // Initialize data structures
  score = allocate_score(MAXN,BANDWIDTH);

  finish = time(0);
  duration = difftime(finish, start);

  if ( VERBOSE ) {
    printCoreEdges(coreAlign, seqLib, 0, 1);
  }else {
    printCoreEdges(coreAlign, seqLib, 0, 0);
  }

  // Initialize master to single N ( since l=1 in this implementation )
  for (x = 0; x < l; x++)
    master[L + x] = 99;

  // NOTE: At this stage the coreAlign->score fields are initialized to
  // zero.

  //
  // Extend RIGHT first (direction = 1)
  //
  int rightbp = extend_alignment(1, coreAlign, score, seqLib, master,
                                 BANDWIDTH, CAPPENALTY, MINIMPROVEMENT, L, N,
                                 scoreParams, fp_mat);
  printf("Extended right: %d bp\n", rightbp);

  // Update the bounds of the core sequences should they overlap one
  // of these right extensions.
  //
  //  e.g
  //                        1          2           3
  //              123456789 012 345678901234 567 890123456
  //    seq1:     ATTAGCTGT[ata]ACGTATTTCGGT[ata]ACGTAGGTA
  //
  // right extension (with "|" denoting existing bounds):
  //                            >>>>>>>     |    >>>>>    |
  //    seq1:     ATTAGCTGT[ata]ACGTATTTCGGT[tat]ACGTAGGTA
  //
  // update bounds for left extension:
  //             |       <<           |   <<
  //    seq1:     ATTAGCTGT[ata]ACGTATTTCGGT[tat]ACGTAGGTA
  //
  struct coreAlignment *s;
  for (s = coreAlign; s != NULL; s = s->next)
  {
    int s_seqIdx = s->seqIdx;
    char *s_ident = seqLib->identifiers[s_seqIdx];
    uint64_t s_seqLowerBound = 0;
    if (s_seqIdx > 0)
      s_seqLowerBound = seqLib->boundaries[s_seqIdx - 1];

    // The position relative the core sequence
    uint64_t extended_pos = (s->rightSeqPos + s->rightExtensionLen);
    if (s->orient)
      extended_pos = (s->rightSeqPos - s->rightExtensionLen);
    // The aboslute position in the seq_id
    uint64_t seqid_extended_pos =
            seqLib->offsets[s_seqIdx] + (extended_pos - s_seqLowerBound + 1);

    //printf("EXT: s_idx = %d right_ext_len = %d to position %ld\n", s_seqIdx, s->rightExtensionLen, seqid_extended_pos);

    struct coreAlignment *r;
    for (r = coreAlign; r != NULL; r = r->next)
    {
      int r_seqIdx = r->seqIdx;
      char *r_ident = seqLib->identifiers[r_seqIdx];
      uint64_t r_seqLowerBound = 0;
      if (r_seqIdx > 0)
        r_seqLowerBound = seqLib->boundaries[r_seqIdx - 1];

      if ( strcmp(s_ident, r_ident) == 0 && seqid_extended_pos > seqLib->offsets[r_seqIdx] )
      {
        uint64_t pos_in_r = r_seqLowerBound + (seqid_extended_pos - seqLib->offsets[r_seqIdx]);

        if ( r->orient )
        {
          if ( pos_in_r >= r->leftSeqPos && pos_in_r <= r->upperSeqBound )
          {
            //printf("In seqid=%d, the extended pos maps to %ld and is in between left_core %ld and upperbound %ld\n", r_seqIdx, pos_in_r, r->leftSeqPos, r->upperSeqBound);
            //printf("Setting upper bound to %ld, was %ld\n", pos_in_r, r->upperSeqBound);
            printf("OVERLAP AVOIDANCE: seqid %d extended to %ld, limits seqid %d with existing upper_bound = %ld because it's pos_in_r=%ld\n", s_seqIdx, seqid_extended_pos, r_seqIdx, r->upperSeqBound, pos_in_r);
            r->upperSeqBound = pos_in_r;
            r->upperSeqBoundFlag = EXT_BOUNDARY;
          }
        }else {
          if ( pos_in_r >= r->lowerSeqBound && pos_in_r <= r->leftSeqPos )
          {
            //printf("In seqid=%d, the extended pos maps to %ld and is in between lowerbound %ld and right_core %ld\n", r_seqIdx, pos_in_r, r->lowerSeqBound, r->rightSeqPos);
            //printf("Setting lower bound to %ld, was %ld\n", pos_in_r, r->lowerSeqBound);
            printf("OVERLAP AVOIDANCE: seqid %d extended to %ld, limits seqid %d with existing lower_bound = %ld because it's pos_in_r=%ld\n", s_seqIdx, seqid_extended_pos, r_seqIdx, r->lowerSeqBound, pos_in_r);
            r->lowerSeqBound = pos_in_r;
            r->lowerSeqBoundFlag = EXT_BOUNDARY;
          }
        }
      }
    } // for r
  } // for s

  masterend = L + l + rightbp;

  //
  // Extend LEFT second (direction = 0)
  //
  int leftbp = extend_alignment(0, coreAlign, score, seqLib, master,
                                BANDWIDTH, CAPPENALTY, MINIMPROVEMENT, L, N,
                                scoreParams, fp_mat);
  printf("Extended left : %d bp\n", leftbp);

  // To get this back to the master index
  masterstart = L - leftbp;

  // No need to display
  if (rightbp > 0 || leftbp > 0)
  {
    // Some FASTA parsers are unhappy if an ID is present but the sequence is empty
    if ( leftbp > 0 ) {
      printf(">left-extension\n");
      for (x = masterstart; x < L; x++)
      {
        printf("%c", num_to_char(master[x]));
        if ((x - masterstart) % 80 == 79)
          printf("\n");
      }
      if ((x - masterstart) % 80 > 0)
        printf("\n");
    }

    if ( rightbp > 0 ) {
      printf(">right-extension\n");
      for (x = L + l; x < masterend; x++)
      {
        printf("%c", num_to_char(master[x]));
        if ((x - masterstart) % 80 == 79)
          printf("\n");
      }
      if ((x - masterstart) % 80 > 0)
        printf("\n");
    }

    // This was a confusing output type...removing
    //printf(">combined_w_N_spacer\n");
    //for (x = masterstart; x < masterend; x++)
    //{
    //  printf("%c", num_to_char(master[x]));
    //  if ((x - masterstart) % 80 == 79)
    //    printf("\n");
    //}
    //if ((x - masterstart) % 80 > 0)
    //  printf("\n");

    // Save consensus to a file, if requested
    if (CONS_FILE != NULL)
    {
      if ((fp = fopen(CONS_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not open input file %s\n", CONS_FILE);
        exit(1);
      }
      if ( leftbp > 0 ) {
        fprintf(fp, ">left-extension %d bp\n", leftbp);
        for (x = masterstart; x < L; x++)
        {
          fprintf(fp, "%c", num_to_char(master[x]));
          if ((x - masterstart) % 80 == 79)
            fprintf(fp, "\n");
        }
        if ((x - masterstart) % 80 > 0)
          fprintf(fp, "\n");
      }

      if ( rightbp > 0 ) {
        fprintf(fp, ">right-extension %d bp\n", rightbp);
        for (x = L + l; x < masterend; x++)
        {
          fprintf(fp, "%c", num_to_char(master[x]));
          if ((x - masterstart) % 80 == 79)
            fprintf(fp, "\n");
        }
        if ((x - masterstart) % 80 > 0)
          fprintf(fp, "\n");
        fclose(fp);
      }
    }

    // Write out TSV and FA files or screen output
    if (OUTTSV_FILE != NULL)
    {
      if ((fp = fopen(OUTTSV_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not create the TSV output file %s\n", OUTTSV_FILE);
        exit(1);
      }
    }

    if (OUTFA_FILE != NULL)
    {
      if ((fp_fa = fopen(OUTFA_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not create the FASTA output file %s\n", OUTFA_FILE);
        exit(1);
      }
    }

    printf("\n\nExtended Sequences Report:\n");
    printf("  *** Extended sequences are in 1-based, fully closed coordinates ***\n");
    printf("SEQID  SEQSTART SEQEND ORIENT  EXTENSION_DETAILS\n");
    struct coreAlignment *s;
    x = 0;
    for (s = coreAlign; s != NULL; s = s->next)
    {
      int seqIdx = s->seqIdx;
      char *ident;
      uint64_t seqLowerBound = 0;
      uint64_t seqUpperBound = 0;
      if (seqIdx >= 0)
      {
        ident = seqLib->identifiers[seqIdx];
        if (seqIdx > 0)
          seqLowerBound = seqLib->boundaries[seqIdx - 1];
      }
      else
      {
        ident = "unknown";
        seqLowerBound = 0;
      }
      seqUpperBound = seqLib->boundaries[seqIdx];

      // Check if this is a subsequence
      uint64_t subseq_offset = 0;
      if ( seqLib->offsets != NULL && seqLib->offsets[seqIdx] > 0 )
          subseq_offset = seqLib->offsets[seqIdx];

      uint64_t extended_start = s->leftSeqPos - s->leftExtensionLen - seqLowerBound + 1;
      uint64_t extended_end =
        (s->rightSeqPos + s->rightExtensionLen) - seqLowerBound + 1;
      char orient = '+';
      uint64_t core_start = s->leftSeqPos - seqLowerBound + 1 + subseq_offset;
      uint64_t core_end = s->rightSeqPos - seqLowerBound + 1 + subseq_offset;
      if (s->orient)
      {
        orient = '-';
        extended_end =
          (s->leftSeqPos + s->leftExtensionLen) - seqLowerBound + 1;
        extended_start =
          (s->rightSeqPos - s->rightExtensionLen) - seqLowerBound + 1;
        core_start = s->rightSeqPos - seqLowerBound + 1 + subseq_offset;
        core_end = s->leftSeqPos - seqLowerBound + 1 + subseq_offset;
      }
      int extended_length = extended_end - extended_start + 1;

      printf ("%s\t%ld\t%ld\t%c\tn=%ld,anchor_range=%ld-%ld",
              ident, subseq_offset+extended_start, subseq_offset+extended_end, orient,
              x, core_start, core_end);

      if ( s->leftExtendable )
        printf(",extended_left=%d", s->leftExtensionLen);
      else
        printf(",extended_left=*");

      if ( s->rightExtendable )
        printf(",extended_right=%d", s->rightExtensionLen);
      else
        printf(",extended_right=*");

      printf(",len=%d,score=%d", extended_length, s->score);

      if ( orient == '+' ) {
        if ( s->leftExtendable && ((s->leftSeqPos - s->leftExtensionLen) - s->lowerSeqBound) < 20 ) {
          if ( s->lowerSeqBoundFlag == SEQ_BOUNDARY )
            printf(",leftSeqLimit");
          else if ( s->lowerSeqBoundFlag == L_BOUNDARY )
            printf(",leftExtLimit");
          else if ( s->lowerSeqBoundFlag == CORE_BOUNDARY )
            printf(",leftCoreLimit");
        }
        if ( s->rightExtendable && (s->upperSeqBound - (s->rightSeqPos + s->rightExtensionLen)) < 20 ) {
          if ( s->upperSeqBoundFlag == SEQ_BOUNDARY )
            printf(",rightSeqLimit");
          else if ( s->upperSeqBoundFlag == L_BOUNDARY )
            printf(",rightExtLimit");
          else if ( s->upperSeqBoundFlag == CORE_BOUNDARY )
            printf(",rightCoreLimit");
        }
      }else {
        if ( s->leftExtendable && (s->upperSeqBound - (s->leftSeqPos + s->leftExtensionLen)) < 20 ) {
          printf(",leftCoreLimit");
        }
        if ( s->rightExtendable && ((s->rightSeqPos - s->rightExtensionLen) - s->lowerSeqBound) < 20 ) {
          printf(",rightCoreLimit");
        }
      }
      printf("\n");

      if (OUTTSV_FILE != NULL) {
        fprintf
          (fp,
           "%s\t%ld\t%ld\t%c\tn=%ld,anchor_range=%ld-%ld",
           ident, subseq_offset+extended_start, subseq_offset+extended_end, orient, x,
           s->leftSeqPos - seqLowerBound + 1,
           s->rightSeqPos - seqLowerBound + 1 );
        if ( s->leftExtendable )
          fprintf(fp,",extended_left=%d", s->leftExtensionLen);
        else
          fprintf(fp,",extended_left=*");
        if ( s->rightExtendable )
          fprintf(fp,",extended_right=%d", s->rightExtensionLen);
        else
          fprintf(fp,",extended_right=*");
        fprintf (fp, ",len=%d,score=%d\n", extended_length, s->score);
      }

      if (OUTFA_FILE != NULL) {
        // Sanity checks
        if ( s->leftExtensionLen < 0 || s->rightExtensionLen < 0 ) {
          fprintf(stderr, "Error: Negative extension length detected\n");
          exit(1);
        }
        uint64_t j;
        uint64_t intStart = s->leftSeqPos - s->leftExtensionLen;
        uint64_t intEnd = s->rightSeqPos + s->rightExtensionLen;
        if (s->orient)
        {
          intEnd = s->leftSeqPos + s->leftExtensionLen;
          intStart = s->rightSeqPos - s->rightExtensionLen;
        }
        if ( flanking > 0 ) {
          if ( seqLowerBound == 0 && flanking > intStart )
            intStart = 1;
          else
            intStart -= flanking;
          if ( intEnd + flanking > seqUpperBound )
            intEnd = seqUpperBound;
          else
            intEnd += flanking;
          fprintf(fp_fa,">%s:%ld-%ld_%c  n=%ld,anchor_range=%ld-%ld,extended_left=%d,extended_right=%d,len=%d,flanking=%d,score=%d\n",
                     ident, subseq_offset+intStart-seqLowerBound+1, subseq_offset+intEnd-seqLowerBound+1, orient,
                     x, s->leftSeqPos - seqLowerBound + 1,
                     s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
                     s->rightExtensionLen, extended_length, flanking, s->score);
        }else{
          fprintf(fp_fa,">%s:%ld-%ld_%c  n=%ld,anchor_range=%ld-%ld,extended_left=%d,extended_right=%d,len=%d,score=%d\n",
                     ident, subseq_offset+intStart-seqLowerBound+1, subseq_offset+intEnd-seqLowerBound+1, orient,
                     x, s->leftSeqPos - seqLowerBound + 1,
                     s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
                     s->rightExtensionLen, extended_length, s->score);
        }
        int seqEmitted = 0;
        if ( orient == '-' ){
          // j is a uint64_t so must be careful with decrements inclusive of zero:
          j = intEnd;
          do {
            fprintf(fp_fa,"%c",num_to_char(compl(seqLib->sequence[j])));
            seqEmitted++;
          }while ( j-- > intStart );
        }else{
          for ( j = intStart; j <= intEnd; j++ ) {
            fprintf(fp_fa,"%c",num_to_char(seqLib->sequence[j]));
            seqEmitted++;
          }
        }
        fprintf(fp_fa,"\n");
        // Sanity check to ensure we do not create any outputs with zero-length.  This
        // should never happen *and* is problematic for a known issue in rmblastn.
        if ( seqEmitted == 0 ) {
          fprintf(stderr, "Error: No sequence emitted for %s:%ld-%ld_%c\n", ident, intStart, intEnd, orient);
          exit(1);
        }
      }

      x++;
    }
    if (OUTTSV_FILE != NULL)
      fclose(fp);
    if (OUTFA_FILE != NULL)
      fclose(fp_fa);
    if (OUTMAT_FILE != NULL)
      fclose(fp_mat);

  }

  // Free up endpoint arrays
  finish = time(0);
  duration = difftime(finish, start);
  printf("Program duration is %.1f sec = %.1f min = %.1f hr\n", duration,
         duration / 60.0, duration / 3600.0);

  return 0;
}


void
print_parameters()
{
  printf("--------------------------------------------------------------\n");
  printf("Parameters:\n");
  printf("  VERBOSE %d\n", VERBOSE);
  printf("  SEQUENCE_FILE %s\n", SEQUENCE_FILE);
  printf("  RANGES_FILE %s\n", RANGES_FILE);
  if (num_threads)
    printf("  Multi-Threaded-Masking (EXPERIMENTAL): num_threads = %d\n",
           num_threads);
  printf("  L %d\n", L);
  printf("  BANDWIDTH (bandwidth) %d\n", BANDWIDTH);
  printf("  MAXN %d\n", MAXN);
  if (strcmp(MATRIX_FILE, "repeatscout") == 0)
  {
    printf("  SCORING SYSTEM: Original RepeatScout method\n");
    printf("     - GAP = %d\n", scoreParams->gapextn);
    printf("     - MATCH = %d\n", MATCH);
    printf("     - MISMATCH = %d\n", MISMATCH);
    printf("     - CAPPENALTY %d\n", CAPPENALTY);
    printf("     - MINIMPROVEMENT %d\n", MINIMPROVEMENT);
  }
  else
  {
    printf("  SCORING SYSTEM: Internally coded matrix '%s'\n", MATRIX_FILE);
    printf("     - GAP_OPEN = %d\n", scoreParams->gapopen);
    printf("     - GAP_EXT = %d\n", scoreParams->gapextn);
    printf("     - CAPPENALTY %d\n", CAPPENALTY);
    printf("     - MINIMPROVEMENT %d\n", MINIMPROVEMENT);
  }
  printf("  WHEN_TO_STOP %d\n", WHEN_TO_STOP);
  printf("--------------------------------------------------------------\n");
}


// TODO move this to common?
void
printCoreAlignmentStruct(struct coreAlignment *coreAlign)
{
  printf("coreAlignment:\n");
  printf("  seqIdx           : %d\n", coreAlign->seqIdx);
  printf("  leftSeqPos       : %ld (relative to sequence[])\n",
         coreAlign->leftSeqPos);
  printf("  rightSeqPos      : %ld (relative to sequence[]\n",
         coreAlign->rightSeqPos);
  if (coreAlign->leftExtendable == 1)
    printf("  leftExtendable   : %d = yes\n", coreAlign->leftExtendable);
  else
    printf("  leftExtendable   : %d = no\n", coreAlign->leftExtendable);
  if (coreAlign->rightExtendable)
    printf("  rightExtendable  : %d = yes\n", coreAlign->rightExtendable);
  else
    printf("  rightExtendable  : %d = no\n", coreAlign->rightExtendable);
  printf("  leftExtensionLen : %d bp\n", coreAlign->leftExtensionLen);
  printf("  rightExtensionLen: %d bp\n", coreAlign->rightExtensionLen);
  printf("  score            : %d\n", coreAlign->score);
  if (coreAlign->orient == 1)
    printf("  orient           : %d = reverse strand\n", coreAlign->orient);
  else
    printf("  orient           : %d = forward strand\n", coreAlign->orient);
}


//  Returns:
//    The index in master for the end of the consensus extension
int
extend_alignment(int direction, struct coreAlignment *coreAlign,
                 int ****score, struct sequenceLibrary *seqLib, char *master, int BANDWIDTH,
                 int CAPPENALTY, int MINIMPROVEMENT, int L, int N,
                 struct scoringSystem *scoreParams,
                 FILE *pathStringFile)
{
  int n, offset;
  int row_idx;                  /* The current row/col of the matrix on which
                                   the computation band is centered. */
  int curr_row_best_col_idx;    /* The sequence offset (from pos[n]) which
                                   has the highest score for a given row
                                   and matrix matrix. */
  int curr_row_best_score;      /* The highest score in a row of a particular
                                   matrix. */
  int score_given_cons;         /* The sum of curr_row_best_score for all
                                   matrices given a putative cosensus base. */
  int curr_extension_score = 0; /* The max sum of score_given_cons for all
                                   matrices for the current row. */
  int max_extension_score = 0;  /* The maximum curr_extension_score seen
                                   thus far ( given MINIMPROVEMENT requirement ). */
  int max_extension_score_row_idx;      /* The row idx for the best scoring.
                                           extension so far */
  struct coreAlignment *currCore;
  char a, besta;
  uint64_t lower_seq_bound, upper_seq_bound;

  if (VERBOSE >= 3)
  {
    if (direction)
      printf("extend_alignment(right): Called with %d edges\n", N);
    else
      printf("extend_alignment(left): Called with %d edges\n", N);
  }

  char *pathString = NULL;
  if ( pathStringFile != NULL ) {
    pathString = malloc(((BANDWIDTH*2) + 1) * sizeof(char));
  }

  // Allocate highscore position and score arrays
  int *overall_sequence_high_score = (int *)malloc(N*sizeof(int));
  int *overall_sequence_high_score_pos = (int *)malloc(N*sizeof(int));
  int *trimmed_sequence_high_score = (int *)malloc(N*sizeof(int));
  int *trimmed_sequence_high_score_pos = (int *)malloc(N*sizeof(int));

  //
  // initialize boundary conditions in score[1][][] and
  // sequence score/position arrays
  //
  for (n = 0; n < N; n++)
  {
    overall_sequence_high_score[n] = 0;
    overall_sequence_high_score_pos[n] = 0;
    trimmed_sequence_high_score[n] = 0;
    trimmed_sequence_high_score_pos[n] = 0;
    for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
    {
      // Affine GAP penalties
      if (offset < 0)
      {
        // orig
        score[1][n][offset + BANDWIDTH][1] = 0;
        score[1][n][offset + BANDWIDTH][1] +=
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        score[1][n][offset + BANDWIDTH][0] = 0;
        score[1][n][offset + BANDWIDTH][0] +=
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
      }
      else
      {
        if (offset > 0)
        {
          // orig
          score[1][n][offset + BANDWIDTH][1] = 0;
          score[1][n][offset + BANDWIDTH][1] +=
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          score[1][n][offset + BANDWIDTH][0] = 0;
          score[1][n][offset + BANDWIDTH][0] +=
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
        }
        else
        {
          score[1][n][offset + BANDWIDTH][0] = 0;
          score[1][n][offset + BANDWIDTH][1] = 0;
        }
      }
    }
    overall_sequence_high_score[n] = 0;

    if (VERBOSE >= 12)
    {
      printf("SW Matrix Boundary Conditions ( n = %d ):\n", n);
      printf("  GAP: ");
      for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
        printf(" %d", score[1][n][offset + BANDWIDTH][1]);
      printf("\n  SUB: ");
      for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
        printf(" %d", score[1][n][offset + BANDWIDTH][0]);
      printf("\n");
    }
  }

  //
  // Extend
  //

  // Initialize best extension index ( 0 to L-1 )
  max_extension_score_row_idx = -1;

  // Rows are the consensus positions
  for (row_idx = 0; row_idx < L; row_idx++)
  {
    // Find the best consensus base
    curr_extension_score = 0;
    besta = 0;
    for (a = 0; a < 4; a++)
    {
      score_given_cons = 0;

      n = 0;
      for (currCore = coreAlign; currCore != NULL; currCore = currCore->next)
      {
        // Only consider the sequences we are told are extendable in
        // this direction
        if ((direction && currCore->rightExtendable) ||
            (!direction && currCore->leftExtendable))
        {
          // Use boundaries defined in the core datastructure.  These may
          // now be hard or soft boundaries.
          lower_seq_bound = currCore->lowerSeqBound;
          upper_seq_bound = currCore->upperSeqBound;

          if (VERBOSE >= 10)
          {
            if ( direction )
              printf("RIGHT ROW %d with '%c': n = %d\n", row_idx, num_to_char(a), n);
            else
              printf("LEFT ROW %d with '%c': n = %d\n", row_idx, num_to_char(a), n);
          }

          // Compute a row of the banded DP matrix under a given consensus
          // hypothesis ("a").
//if ( (row_idx == 21) )  {
//    VERBOSE = 17;
//        }
          curr_row_best_score = compute_nw_row(direction, row_idx, n, a, currCore,
                                               score, lower_seq_bound,
                                               upper_seq_bound,
                                               seqLib->sequence,
                                               &curr_row_best_col_idx,
                                               scoreParams, BANDWIDTH, L, VERBOSE, NULL);
//VERBOSE = 0;

          if (VERBOSE >= 12)
//if ( (row_idx >= 20 && row_idx <= 21) ) 
          {
            printf("    Gap: ");
            for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
              printf(" %d", score[row_idx % 2][n][offset + BANDWIDTH][1]);
            printf("\n");
            printf("    Sub: ");
            for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
              printf(" %d", score[row_idx % 2][n][offset + BANDWIDTH][0]);
            printf("\n");
          }

          if (VERBOSE >= 10)
          {
            if ( curr_row_best_score < 0 )
              printf("  best score = %d @ column %d -- max(0,best_score) = 0!, prev best score = %d\n",
                   curr_row_best_score, curr_row_best_col_idx,
                   overall_sequence_high_score[n]);
            else
              printf("  best score = %d @ column %d, prev best score = %d\n",
                   curr_row_best_score, curr_row_best_col_idx,
                   overall_sequence_high_score[n]);
            // RMH: TODO...is this the best/only way to determine that a sequence has hit a limit?
            if (  ( direction && score[row_idx % 2][n][0][1] < -279000 ) ||
                  ( !direction && score[row_idx % 2][n][BANDWIDTH+BANDWIDTH][1] < -279000 ))
              printf(" **OUT_OF_SEQ**");
          }

          if ( curr_row_best_score < 0 )
             curr_row_best_score = 0;

//if (( row_idx >= 20 && row_idx <= 21) ) {
// printf("    c=%c  n=%d  => score = %d\n",num_to_char(a), n,
//            curr_row_best_score);
//}

          //       curr_row_best_score, n, overall_sequence_high_score[n],
          //       CAPPENALTY);
          if (curr_row_best_score >= (overall_sequence_high_score[n] + CAPPENALTY))
          {
            score_given_cons += curr_row_best_score;
          }
          else
          {
            if (VERBOSE >= 10)
              printf(" **CAPPED** contributing = %d",
                     (overall_sequence_high_score[n] + CAPPENALTY));
            score_given_cons += overall_sequence_high_score[n] + CAPPENALTY;
          }
          if (VERBOSE >= 10)
            printf("\n");

          if (VERBOSE >= 12)
            printExtensionRegion( currCore, direction, row_idx, 0, seqLib );

        }                       // if extendable
        n++;
      }                         // for ( currCore=...

//if (( row_idx >= 20 && row_idx <= 21) ) {
//  printf("  Total Score for '%c' = %d\n", num_to_char(a),
//               score_given_cons);
//}

      if (VERBOSE >= 10)
        printf("  Total Score for '%c' = %d\n", num_to_char(a),
               score_given_cons);
      if (score_given_cons > curr_extension_score)
      {
        curr_extension_score = score_given_cons;
        besta = a;
      }
    }   // for (a=0...

    if (VERBOSE >= 10)
      printf("ROW %d complete, '%c' chosen as the consensus. curr_ext_score = %d\n", row_idx,
             num_to_char(besta), curr_extension_score);

    if (direction)
      master[L + l + row_idx] = besta;
    else
      master[L - row_idx - 1] = besta;

    // Now recalculate scores for the chosen consensus base (besta).  NOTE
    // this is strictly not necessary for besta='T' as that was the last
    // one calculated above.  However to keep the code simple, we will simply
    // recalculate it here.
    int totalExtendable = 0;
    int numOutOfSeq = 0;
    int numHighScore = 0;
    n = 0;
    for (currCore = coreAlign; currCore != NULL; currCore = currCore->next)
    {
      if ((direction && currCore->rightExtendable) ||
          (!direction && currCore->leftExtendable))
      {
        // Use boundaries defined in the core datastructure.  These may
        // now be hard or soft boundaries.
        lower_seq_bound = currCore->lowerSeqBound;
        upper_seq_bound = currCore->upperSeqBound;

        // Rebuild the DP matrix row based on the new consensus call ('besta').
        curr_row_best_score =
            compute_nw_row(direction, row_idx, n, besta, currCore, score,
                               lower_seq_bound, upper_seq_bound, seqLib->sequence,
                               &curr_row_best_col_idx,
                               scoreParams, BANDWIDTH, L, 0, pathString);

        // Optional, save matrix path to file for debugging
        if ( pathString ) {
          pathString[BANDWIDTH*2] = '\0';
          if ( direction == 1 )
            fprintf(pathStringFile, "dir=right");
          else
            fprintf(pathStringFile, "dir=left");
          fprintf(pathStringFile, ":n=%d:row=%d: %s best:score=%d:offset=%d:cons=%c\n",
                  n, row_idx, pathString, curr_row_best_score,
                  curr_row_best_col_idx-row_idx+BANDWIDTH, num_to_char(besta));
        }

        // How many have run out of sequence?
        if ( (! direction && score[row_idx % 2][n][0][1] < -279000  ) ||
             (direction && score[row_idx % 2][n][BANDWIDTH+BANDWIDTH][1] < -279000 ))
          numOutOfSeq++;

        // Record setters
        if (curr_row_best_score > overall_sequence_high_score[n])
        {
          //  The current sequence has a higher score than its previous best.
          //  Save this value, even if the overall msa no longer meets the
          //  MINIMPROVEMENT criteria.  Should it start to meet the criteria
          //  again, we can simply grab these values and store them in the
          //  trimmed_sequence_high_score/pos arrays.
          overall_sequence_high_score[n] = curr_row_best_score;
          overall_sequence_high_score_pos[n] = curr_row_best_col_idx;
          numHighScore++;
        }
        totalExtendable++;
      } // if extendable

//if ( (row_idx >= 20 && row_idx <= 21) ) 
//{
//printf("    GapAFTER : ");
//for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
//printf(" %d", score[row_idx % 2][n][offset + BANDWIDTH][1]);
//printf("\n");
//printf("    SubAFTER: ");
//for (offset = -BANDWIDTH; offset <= BANDWIDTH; offset++)
//printf(" %d", score[row_idx % 2][n][offset + BANDWIDTH][0]);
//printf("\n");
//}


      n++;
    }
    float score_per_position_since_max = 0.0;
    if (abs(row_idx - max_extension_score_row_idx) > 0)
    {
      score_per_position_since_max =
        ((float) (curr_extension_score - max_extension_score) /
         abs(row_idx - max_extension_score_row_idx));
    }

    if (VERBOSE >= 10)
      printf("Alignment Extension: curr_extension_score = %d, "
             "max_extension_score = %d @ row %d\n                     "
             "score_per_position_overall = %0.1f, "
             "score_per_position_since_max = %0.1f\n"
             "                     total_edges = %d, num_extending = %d, num_out_of_seq = %d\n",
             curr_extension_score, max_extension_score,
             max_extension_score_row_idx,
             ((float) curr_extension_score / (row_idx + 1)),
             score_per_position_since_max, totalExtendable, numHighScore, numOutOfSeq );

//if ( row_idx >= 20 && row_idx <= 21 ) {
//printf("===> ext score %d ... row_idx = %d\n", curr_extension_score, row_idx);
//}

    // What is the best extension score we have seen with sequences meeting
    // our minmum improvement per position criteria?
    if ((curr_extension_score >=
         max_extension_score +
         (abs(max_extension_score_row_idx - row_idx) * MINIMPROVEMENT)))
    {
      if ( VERBOSE >= 10)
        printf("                     **This row is now the new max**\n");
      max_extension_score_row_idx = row_idx;
      max_extension_score = curr_extension_score;
      int nn;
      for (nn = 0; nn < N; nn++)
      {
        trimmed_sequence_high_score[nn] = overall_sequence_high_score[nn];
        trimmed_sequence_high_score_pos[nn] = overall_sequence_high_score_pos[nn];
      }
    }
    else
    {
      if (VERBOSE >= 10)
        printf("                     extensions since last max score %d\n",
               abs(row_idx - max_extension_score_row_idx));
    }

    if (abs(row_idx - max_extension_score_row_idx) >= WHEN_TO_STOP)
    {
      if (VERBOSE >= 3)
        printf
          ("Ending...due to row_idx=%d - max_extension_score_row_idx=%d <= -WHEN_TO_STOP=%d\n",
           row_idx, max_extension_score_row_idx, WHEN_TO_STOP);
      break;
    }
  } // for (row_idx
  if (row_idx == L - 1)
  {
     if ( direction )
       printf("WARNING: Extended sequence right to the limit ( L=%d ).\n", L);
     else
       printf("WARNING: Extended sequence left to the limit ( L=%d ).\n", L);
  }

  // Update the core extensionLen and score fields.
  n=0;
  for (currCore = coreAlign; currCore != NULL && n < N; currCore = currCore->next)
  {
    if ( trimmed_sequence_high_score[n] > 0  && trimmed_sequence_high_score_pos[n] >= 0)
    {
      if (direction) {
        currCore->rightExtensionLen = trimmed_sequence_high_score_pos[n] + 1;
      }else {
        currCore->leftExtensionLen = trimmed_sequence_high_score_pos[n] + 1;
      }
      currCore->score += trimmed_sequence_high_score[n];
    }
    n++;
  }

  // Release memory
  free(overall_sequence_high_score);
  free(overall_sequence_high_score_pos);
  free(trimmed_sequence_high_score);
  free(trimmed_sequence_high_score_pos);
  free( pathString );

  // Return the total number of extended consensus bases
  return max_extension_score_row_idx + 1;
}

