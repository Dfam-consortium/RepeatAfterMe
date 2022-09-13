/*
 * ExtendAlign
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
//typedef unsigned int uint;
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
#include "extend_align.h"
 

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
char *RANGES_FILE = NULL;       //
int VERBOSE;                    // How chatty? Real chatty? Really really
                                // chatty?
                                // Super extra really chatty?
int L;                          // The maximum distance to attempt extension of each
                                //   core range. Length of master is 2*L+l (L >>> l)
int flanking = 0;               // The number of bp to add to the left/right side
                                //   of each -outfa sequence.
int MAXOFFSET;                  // max offset (5)
int MAXN;                       // max #occ of lmer (10000)
int WHEN_TO_STOP;               // stop if no improvement after extending
                                // this far (100)
float MINENTROPY;               // ignore freq l-mers with entropy less than
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
int ****score;                  // 2 * MAXN * (2*MAXOFFSET+1) * 2

//
//     Reduce memory use by freeing once no longer needed
//
//TODO: Move this into the core range structure
int *bestbestscore;             // best alignment score for this n for any y
                                // or w seen so far
// Unused
int num_threads = 0;            // Number of threads or 0 if not
                                // multithreaded


void
usage()
{
  printf(
          "ExtendAlign Version %s - build %s date %s\n"
          "\n"
          "  Perform a multiple alignment extension given an existing\n"
          "  core multiple alignment.  The core may be a set of word\n"
          "  matches, a previously development multiple alignment, or the\n"
          "  the result of any process that defines a core set of sequence\n"
          "  relationships.  The only requirement is that the sequences are\n"
          "  aligned to nearly the same point along one or both of the core\n"
          "  edges.  The internal aspects of the core alignment itself are\n"
          "  not used taken into account during extension.\n"
          "\n"
          "Usage: \n"
          "  ExtendAlign -twobit <seq.2bit> -ranges <ranges.tsv> [opts]\n"
          "\n"
          "--------------------------------------------------------------------------------------------\n"
          "     -cons <seq.fa>        # Save the left/right consensus sequences to a FASTA file.\n"
          "     -outtsv <final.tsv>   # Save the final sequence ranges to a TSV file.\n"
          "     -outfa <seq.fa>       # Save all the final sequences ( original range + extension )\n"
          "                           #   to a file.\n"
          "     -addflanking <num>    # Include an additional<num> bp of sequence to each sequence\n"
          "                           #   when using the -outfa option.\n"
          "     -L <num>              # Size of region to extend left or right (10000). \n"
          "     -matrix 14p43g |      # There are several internally-coded matrices available.\n"
          "             18p43g |      #   The original RepeatScout scoring system is 'repeatscout'\n"
          "             20p43g |      #   and encodes a simple +1/-1 substitution matrix with a single\n"
          "             25p43g |      #   gap penalty.  The other choices are encodings of the\n"
          "             repeatscout   #   RepeatMasker matrices specified by divergence and GC\n"
          "                           #   background as '<divergence>p<GC>g' ( eg. 20p43g\n"
          "                           #   is tuned for 20%% divergence in a 43%% GC background ).\n"
          "                           #   Each matrix has it's own default gap_open, gap_extension\n"
          "                           #   parameters.  These may be overrided by using the -gapopen\n"
          "                           #   and -gapextn options.\n"
          "     -match <num>          # Reward for a match (+1: original scoring system only)  \n"
          "     -mismatch <num>       # Penalty for a mismatch (-1: original scoring system only) \n"
          "     -gap <num>            # Penalty for a gap (-5: original scoring system only)\n"
          " or\n"
          "     -gapopen <num>        # Gap open penalty (-28: nucleotide scoring matrix only)\n"
          "     -gapextn <num>        # Gap extension penalty (-6: nucleotide scoring matrix only)\n"
          "\n"
          "     -minimprovement <num> # Amount that a the alignment needs to improve each step\n"
          "                           #   to be considered progress (original: 3, nucleotide matrix: 27)\n"
          "     -bandwidth <num>      # The maximum number of unbalanced gaps allowed (14)\n"
          "                           #   Half the bandwidth of the banded Smith-Waterman\n"
          "                           #   algorithm.  The default allows for at most 14bp\n"
          "                           #   of unbalanced insertions/deletions in any aligned\n"
          "                           #   sequence.  The full bandwidth is 2*bandwidth+1.\n"
          "     -maxoccurrences <num> # Cap on the number of sequences to align (10,000) \n"
          "     -cappenalty <num>     # Cap on penalty for exiting alignment of a sequence \n"
          "                           #   (original: -20, nucleotide matrix:-90)\n"
          "     -stopafter <num>      # Stop the alignment after this number of no-progress columns (100)\n"
          "     -minlength <num>      # Minimum required length for a sequence to be reported (50)\n"
          "     -v[v[v[v]]]           # How verbose do you want it to be?  -vvvv is super-verbose\n"
          "     -unittests            # Run unit tests\n"
          "\n"
          "Ranges:\n"
          "   Ranges are supplied in the form of a modified BED-6 format:\n"
          "      field-1:chrom     : sequence identifier\n"
          "      field-2:chromStart: lower aligned position ( 0 based )\n"
          "      field-3:chromEnd  : upper aligned position ( 0 based, half open )\n"
          "      field-4:name      : left extendable flag ( 0 = no, 1 = yes )\n"
          "      field-5:score     : right extendable flag \n"
          "      field-6:strand    : strand ( '+' = forward, '-' = reverse )\n"
          "   The fields are tab separated. Coordinates are zero-based half\n"
          "   open.\n"
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
  start = time(0);

  if (co_get_string(argc, argv, "-ranges", &RANGES_FILE) == 0)
  {
    usage();
    exit(1);
  }

  co_get_int(argc, argv, "-addflanking", &flanking);
  co_get_string(argc, argv, "-outtsv", &OUTTSV_FILE);
  co_get_string(argc, argv, "-outfa", &OUTFA_FILE);
  co_get_string(argc, argv, "-cons", &CONS_FILE);
  co_get_int(argc, argv, "-L", &L) || (L = 10000);
  co_get_int(argc, argv, "-bandwidth", &MAXOFFSET) || (MAXOFFSET = 14);
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

  if (co_get_bool(argc, argv, "-unittests", &opt))
  {
    test_global_align();
    test_cons_seed_extend();
    test_compute_nw_row();
    exit(1);
  }

  // Consensus sequence which is lmer + left and right extension max (L)
  master = (char *) malloc((2 * L + l + 1) * sizeof(char));
  if (NULL == master)
  {
    fprintf(stderr, "Could not allocate space for master array\n");
    exit(1);
  }
  master[2 * L + l] = '\0';

  bestbestscore = (int *) malloc(MAXN * sizeof(int));
  if (NULL == bestbestscore)
  {
    fprintf(stderr, "Could not allocated space for internal arrays\n");
    exit(1);
  }

  if (co_get_string(argc, argv, "-twobit", &SEQUENCE_FILE))
  {
    dnaUtilOpen();
    seqLib = loadSequenceSubset(SEQUENCE_FILE, RANGES_FILE,
                                &coreAlign, &N);
    printf("Read in %d ranges, and %ld bp of sequence\n", N, seqLib->length);
  }
  else if (co_get_string(argc, argv, "-sequence", &SEQUENCE_FILE))
  {
    printf("-sequence is deprecated!....may return someday\n");
    exit(1);
    // Crudely estimate the sequence length.  I.e it can't be larger
    // than the length ( in bytes ) of the input file itself.  However,
    // this is a crude estimate because it is often much smaller than
    // the input file size after accounting for line terminations and
    // FASTA identifier strings.
    fp = fopen(SEQUENCE_FILE, "ro");
    if (NULL == fp)
    {
      fprintf(stderr, "Could not open sequence file %s\n", SEQUENCE_FILE);
      exit(1);
    }
    fseek(fp, 0, SEEK_END);
    off_t ft_size = ftello(fp);
    // Sanitize the cast to unsigned
    if (ft_size > 0)
      seqLength = (uint64_t) ft_size;
    else
    {
      printf("Error: input sequence file is empty!\n");
      exit(1);
    }

    //sequence = (char *) malloc((seqLength + 1) * sizeof(char));
    //if (NULL == sequence)
   // {
   //   fprintf(stderr, "Could not allocate space for sequence\n");
   //   exit(1);
   // }
    // build sequence : We store the forward sequence and calculate
    // the reverse on the fly. Calculate the true size of the input
    // sequence and update the global variable here.
    //seqLength = build_sequence(sequence, SEQUENCE_FILE);

    // DEBUG
    finish = time(0);
    duration = difftime(finish, start);
    printf
      ("Sequences read in (%ld bp): Program duration is %.1lf sec = %.1lf min = %.1lf hr\n",
       seqLength, duration, duration / 60.0, duration / 3600.0);
    // DEBUG

    // Setup the locations of the aligned sequences
    //read_ranges(MAXN, RANGES_FILE, sequence_idents, leftpos, rightpos, rev);
    printf("Read in %d anchor ranges\n", N);
  }
  else
  {
    usage();
    exit(1);
  }

  //
  // print parameters
  //
  print_parameters();

  // Initialize data structures
  score = allocate_score(MAXN,MAXOFFSET);

  finish = time(0);
  duration = difftime(finish, start);

// DEBUG
  printCoreEdges(coreAlign, seqLib, 0);

  // Initialize master to single N ( since l=1 in this implementation )
  for (x = 0; x < l; x++)
    master[L + x] = 99;

  // printf("Extending right: ");
  // Direction: 1=right, 0=left
  int rightbp = extend_alignment(1, coreAlign, score, seqLib, master,
                                 MAXOFFSET, CAPPENALTY, MINIMPROVEMENT, L, N,
                                 bestbestscore, scoreParams);
  printf("Extended right: %d bp\n", rightbp);
  // printf("Extending left: ");
  masterend = L + l + rightbp;
  // Returns number of bases extended
  int leftbp = extend_alignment(0, coreAlign, score, seqLib, master,
                                MAXOFFSET, CAPPENALTY, MINIMPROVEMENT, L, N,
                                bestbestscore, scoreParams);
  printf("Extended left : %d bp\n", leftbp);
  // To get this back to the master index
  masterstart = L - leftbp;

  // No need to display
  if (rightbp > 0 || leftbp > 0)
  {
    printf(">left-extension\n");
    for (x = masterstart; x < L; x++)
    {
      printf("%c", num_to_char(master[x]));
      if ((x - masterstart) % 80 == 79)
        printf("\n");
    }
    if ((x - masterstart) % 80 > 0)
      printf("\n");

    printf(">right-extension\n");
    for (x = L + l; x < masterend; x++)
    {
      printf("%c", num_to_char(master[x]));
      if ((x - masterstart) % 80 == 79)
        printf("\n");
    }
    if ((x - masterstart) % 80 > 0)
      printf("\n");

    printf(">combined_w_N_spacer\n");
    for (x = masterstart; x < masterend; x++)
    {
      printf("%c", num_to_char(master[x]));
      if ((x - masterstart) % 80 == 79)
        printf("\n");
    }
    if ((x - masterstart) % 80 > 0)
      printf("\n");

    // Now to a file, if requested
    if (CONS_FILE != NULL)
    {
      if ((fp = fopen(CONS_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not open input file %s\n", CONS_FILE);
        exit(1);
      }
      fprintf(fp, ">left-extension %d bp\n", leftbp);
      for (x = masterstart; x < L; x++)
      {
        fprintf(fp, "%c", num_to_char(master[x]));
        if ((x - masterstart) % 80 == 79)
          fprintf(fp, "\n");
      }
      if ((x - masterstart) % 80 > 0)
        fprintf(fp, "\n");

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


    // Write out TSV and FA files or screen output
    if (OUTTSV_FILE != NULL)
    {
      if ((fp = fopen(OUTTSV_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not open TSV output file %s\n", OUTTSV_FILE);
        exit(1);
      }
    }

    if (OUTFA_FILE != NULL)
    {
      if ((fp_fa = fopen(OUTFA_FILE, "w")) == NULL)
      {
        fprintf(stderr, "Could not open FASTA output file %s\n", OUTFA_FILE);
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

      //
      //if ( s->leftSeqPos - s->leftExtensionLen < seqLowerBound )
      //  extended_start = 1;
      //else

      uint64_t extended_start = s->leftSeqPos - s->leftExtensionLen - seqLowerBound + 1;
      uint64_t extended_end =
        (s->rightSeqPos + s->rightExtensionLen) - seqLowerBound + 1;
      char orient = '+';
      if (s->orient)
      {
        orient = '-';
        extended_end =
          (s->leftSeqPos + s->leftExtensionLen) - seqLowerBound + 1;
        extended_start =
          (s->rightSeqPos - s->rightExtensionLen) - seqLowerBound + 1;
      }
      int extended_length = extended_end - extended_start + 1;

      printf
        ("%s\t%ld\t%ld\t%c\tn=%ld,anchor_range=%ld-%ld,extended_left=%d,extended_right=%d,len=%d,score=%d\n",
         ident, extended_start, extended_end, orient,
         x, s->leftSeqPos - seqLowerBound + 1,
         s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
         s->rightExtensionLen, extended_length, bestbestscore[x]);

      if (OUTTSV_FILE != NULL)
        fprintf
          (fp,
           "%s\t%ld\t%ld\t%c\tn=%ld,anchor_range=%ld-%ld,extended_left=%d,extended_right=%d,len=%d,score=%d\n",
           ident, extended_start, extended_end, orient, x,
           s->leftSeqPos - seqLowerBound + 1,
           s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
           s->rightExtensionLen, extended_length, bestbestscore[x]);

      if (OUTFA_FILE != NULL) {
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
                     ident, intStart-seqLowerBound+1, intEnd-seqLowerBound+1, orient,
                     x, s->leftSeqPos - seqLowerBound + 1,
                     s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
                     s->rightExtensionLen, extended_length, flanking,bestbestscore[x]);
        }else{
          fprintf(fp_fa,">%s:%ld-%ld_%c  n=%ld,anchor_range=%ld-%ld,extended_left=%d,extended_right=%d,len=%d,score=%d\n",
                     ident, intStart-seqLowerBound+1, intEnd-seqLowerBound+1, orient,
                     x, s->leftSeqPos - seqLowerBound + 1,
                     s->rightSeqPos - seqLowerBound + 1, s->leftExtensionLen,
                     s->rightExtensionLen, extended_length, bestbestscore[x]);
        }
        if ( orient == '-' ){
          // j is a uint64_t so must be careful with decrements inclusive of zero:
          j = intEnd;
          do { 
            fprintf(fp_fa,"%c",num_to_char(compl(seqLib->sequence[j])));
          }while ( j-- > intStart );
        }else{
          for ( j = intStart; j <= intEnd; j++ )
            fprintf(fp_fa,"%c",num_to_char(seqLib->sequence[j]));
        }
        fprintf(fp_fa,"\n");
      }

      x++;
    }
    if (OUTTSV_FILE != NULL)
      fclose(fp);
    if (OUTFA_FILE != NULL)
      fclose(fp_fa);

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
  printf("  MAXOFFSET (bandwidth) %d\n", MAXOFFSET);
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
  if (coreAlign->orient == 1)
    printf("  orient           : %d = reverse strand\n", coreAlign->orient);
  else
    printf("  orient           : %d = forward strand\n", coreAlign->orient);
}


// bestbestscore[n] should be renamed
//  Returns:
//    The index in master for the end of the consensus extension
int
extend_alignment(int direction, struct coreAlignment *coreAlign,
                 int ****score, struct sequenceLibrary *seqLib, char *master, int MAXOFFSET,
                 int CAPPENALTY, int MINIMPROVEMENT, int L, int N,
                 int *bestbestscore, struct scoringSystem *scoreParams)
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

  //
  // initialize boundary conditions in score[1][][]
  //
  for (n = 0; n < N; n++)
  {
    for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    {
      // Affine GAP penalties
      if (offset < 0)
      {
        // orig
        score[1][n][offset + MAXOFFSET][1] = 0;
        score[1][n][offset + MAXOFFSET][1] +=
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        score[1][n][offset + MAXOFFSET][0] = 0;
        score[1][n][offset + MAXOFFSET][0] +=
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
      }
      else
      {
        if (offset > 0)
        {
          // orig
          score[1][n][offset + MAXOFFSET][1] = 0;
          score[1][n][offset + MAXOFFSET][1] +=
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          score[1][n][offset + MAXOFFSET][0] = 0;
          score[1][n][offset + MAXOFFSET][0] +=
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
        }
        else
        {
          score[1][n][offset + MAXOFFSET][0] = 0;
          score[1][n][offset + MAXOFFSET][1] = 0;
        }
      }
    }
    bestbestscore[n] = 0;

    if (VERBOSE >= 10)
    {
      printf("SW Matrix Boundary Conditions ( n = %d ):\n", n);
      printf("  GAP: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", score[1][n][offset + MAXOFFSET][1]);
      printf("\n  SUB: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", score[1][n][offset + MAXOFFSET][0]);
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
          // Sequence boundaries:
          //   boundaries[] contains the start position for each sequence. 
          //   Therefore the range of the sequence is
          //   boundaries[i] to boundaries[i+1]-1 inclusive.
          lower_seq_bound = 0;
          if ( currCore->seqIdx > 0 )
             lower_seq_bound = seqLib->boundaries[currCore->seqIdx - 1];
          upper_seq_bound = seqLib->boundaries[currCore->seqIdx] - 1;

          if (VERBOSE >= 10){
            if ( direction )
              printf("RIGHT ROW %d with '%c': n = %d\n", row_idx, num_to_char(a), n);
            else
              printf("LEFT ROW %d with '%c': n = %d\n", row_idx, num_to_char(a), n);
          }

          curr_row_best_score = compute_nw_row(direction, row_idx, n, a, currCore,
                                                      score, lower_seq_bound,
                                                      upper_seq_bound,
                                                      seqLib->sequence,
                                                      &curr_row_best_col_idx,
                                                      scoreParams, MAXOFFSET, L, VERBOSE);


          if (VERBOSE >= 12)
          {
            printf("    Gap: ");
            for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
              printf(" %d", score[row_idx % 2][n][offset + MAXOFFSET][1]);
            printf("\n");
            printf("    Sub: ");
            for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
              printf(" %d", score[row_idx % 2][n][offset + MAXOFFSET][0]);
            printf("\n");
          }

          if (VERBOSE >= 10) {
            printf("  best score = %d @ column %d, prev best score = %d",
                   curr_row_best_score, curr_row_best_col_idx,
                   bestbestscore[n]);
            if (  ( direction && score[row_idx % 2][n][0][1] < -279000 ) ||
                  ( !direction && score[row_idx % 2][n][MAXOFFSET+MAXOFFSET][1] < -279000 ))
              printf(" **OUT_OF_SEQ**");
          }

          if ((curr_row_best_score - bestbestscore[n]) > CAPPENALTY)
          {
            score_given_cons += curr_row_best_score;
          }
          else
          {
            if (VERBOSE >= 10)
              printf(" **CAPPED** contributing = %d",
                     (bestbestscore[n] + CAPPENALTY));
            score_given_cons += bestbestscore[n] + CAPPENALTY;
          }
          if (VERBOSE >= 10)
            printf("\n");

          if (VERBOSE >= 12)
            printExtensionRegion( currCore, direction, row_idx, 0, seqLib );

        }                       // if extendable
        n++;
      }                         // for ( currCore=...
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
      printf("ROW %d complete, '%c' chosen as the consensus\n", row_idx,
             num_to_char(besta));

    // TODO....drop the use of a contiguous master
    if (direction)
      master[L + l + row_idx] = besta;
    else
      master[L - row_idx - 1] = besta;

    // Now reset the row scores for besta unless besta == 3 (T) since this
    // is what the rows scores currently represent.
    int totalExtendable = 0;
    int numOutOfSeq = 0;
    int numHighScore = 0;
    n = 0;
    for (currCore = coreAlign; currCore != NULL; currCore = currCore->next)
    {
      if ((direction && currCore->rightExtendable) ||
          (!direction && currCore->leftExtendable))
      {

        // Sequence boundaries:
        lower_seq_bound = 0;
        if ( currCore->seqIdx > 0 )
           lower_seq_bound = seqLib->boundaries[currCore->seqIdx - 1];
        upper_seq_bound = seqLib->boundaries[currCore->seqIdx] - 1;
        // Rebuild the matrix row based on the new consensus call
        curr_row_best_score =
            compute_nw_row(direction, row_idx, n, besta, currCore, score,
                               lower_seq_bound, upper_seq_bound, seqLib->sequence,
                               &curr_row_best_col_idx,
                               scoreParams, MAXOFFSET, L, VERBOSE);

        // How many have run out of sequence?
        if ( (! direction && score[row_idx % 2][n][0][1] < -279000  ) ||
             (direction && score[row_idx % 2][n][MAXOFFSET+MAXOFFSET][1] < -279000 ))
          numOutOfSeq++;

        // Record setters
        if (curr_row_best_score > bestbestscore[n])
        {
          // TODO: have generic pre-allocated array passed in
          //Experimental - extra conditional to not allow the LEN to increase unless the whole thing does
          //  The purpose is to not allow increasing score of an individual alignment ( by chance ) to go
          //  beyond the collective extension.
          if ((curr_extension_score >= max_extension_score + (abs(max_extension_score_row_idx - row_idx) * MINIMPROVEMENT)))
          {
            if (direction)
              currCore->rightExtensionLen = curr_row_best_col_idx;
            else
              currCore->leftExtensionLen = curr_row_best_col_idx;
          }
          bestbestscore[n] = curr_row_best_score;
          numHighScore++;
        }
        totalExtendable++;
      } // if extendable
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

    // What is the best extension score we have seen with sequences meeting
    // our minmum improvement per position criteria?
    if ((curr_extension_score >=
         max_extension_score +
         (abs(max_extension_score_row_idx - row_idx) * MINIMPROVEMENT)))
    {
      if (VERBOSE >= 10)
        printf("                     **This row is now the new max**\n");
      max_extension_score_row_idx = row_idx;
      max_extension_score = curr_extension_score;
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

  // Return the total number of extended consensus bases
  return max_extension_score_row_idx + 1;
}


// Unit test for compute_nw_row() function
void
test_compute_nw_row()
{
  int x, n, r;
  int offset;
  int conspos;
  int bestscore;
  struct sequenceLibrary l_seqLib;
  struct sequenceLibrary r_seqLib;
  int max_score_seq_idx = 0;
  int ****l_score; // 2 * MAXN * (2*MAXOFFSET+1) * 2
  int ****r_score; // 2 * MAXN * (2*MAXOFFSET+1) * 2
  // score[2] : two SW rows per lmer/anchor
  l_score = (int ****) malloc(2 * sizeof(*l_score));
  r_score = (int ****) malloc(2 * sizeof(*r_score));
  for (x = 0; x < 2; x++)
  {
    // score[2][5] : up to 5 lmers/anchors
    l_score[x] = (int ***) malloc(5 * sizeof(*l_score[x]));
    r_score[x] = (int ***) malloc(5 * sizeof(*r_score[x]));
    for (n = 0; n < MAXN; n++)
    {
      // score[2][5][2*MAXOFFSET+1] : Row length = 2*bandwidth+1
      l_score[x][n] =
        (int **) malloc((2 * MAXOFFSET + 1) * sizeof(*l_score[x][n]));
      r_score[x][n] =
        (int **) malloc((2 * MAXOFFSET + 1) * sizeof(*r_score[x][n]));
      for (r = 0; r < (2 * MAXOFFSET + 1); r++)
      {
        // score[2][5][2*MAXOFFSET+1][2] : Gap/Sub row per SW row
        l_score[x][n][r] = (int *) malloc(2 * sizeof(*l_score[x][n][r]));
        r_score[x][n][r] = (int *) malloc(2 * sizeof(*r_score[x][n][r]));
      }
    }
  }

  // Boundary conditions are stored in score[1][][][]
  // as first row ( w=0 ), and subsequent rows are stored
  // in: score[w%2][][][]
  for (n = 0; n < 3; n++)
  {
    for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    {
      if (offset < 0)
      {
        l_score[1][n][offset + MAXOFFSET][1] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        l_score[1][n][offset + MAXOFFSET][0] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        r_score[1][n][offset + MAXOFFSET][1] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
        r_score[1][n][offset + MAXOFFSET][0] =
          ((-offset) * scoreParams->gapextn) + scoreParams->gapopen;
      }
      else
      {
        if (offset > 0)
        {
          l_score[1][n][offset + MAXOFFSET][1] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          l_score[1][n][offset + MAXOFFSET][0] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          r_score[1][n][offset + MAXOFFSET][1] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
          r_score[1][n][offset + MAXOFFSET][0] =
            ((offset) * scoreParams->gapextn) + scoreParams->gapopen;
        }
        else
        {
          l_score[1][n][offset + MAXOFFSET][0] = 0;
          l_score[1][n][offset + MAXOFFSET][1] = 0;
          r_score[1][n][offset + MAXOFFSET][0] = 0;
          r_score[1][n][offset + MAXOFFSET][1] = 0;
        }
      }
    }
  }

  l_seqLib.boundaries = (uint64_t *) malloc(4 * sizeof(uint64_t));
  l_seqLib.boundaries[0] = 16;
  l_seqLib.boundaries[1] = 32;
  l_seqLib.boundaries[2] = 48;
  l_seqLib.boundaries[3] = 0;
  uint64_t lower_bound[] = { 0, 16, 32 };
  uint64_t upper_bound[] = { 16, 32, 48 };
  char *predef_id[] = { "seq1", "seq2", "seq3" };
  l_seqLib.identifiers = predef_id;
  l_seqLib.count = 3;
  l_seqLib.length = 48;
  /*
   * n=0  seq1:         9 :   11:+  NACCTGAATC [           T            ] AGGAC*
   * n=0  seq2:        25 :   27:+  GACGTGAATC [           T            ] AGGAC*
   * n=0  seq3:        41 :   43:+  TTGATGAATG [           T            ] AGGAC*
   *
   * 0-3 = ACGT and N = 99
   *
   * LEFT Extension Example1
   *
   *                      row=0  *  *  *  *  *  *  *  *  *  *  *
   *                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
   *                 G  C  A  A  T  G  C  T  G  A  A  T  C  t  a  g */
  char p_left[] = {  2, 1, 0, 0, 3, 2, 1, 3, 2, 0, 0, 3, 1, 3, 0, 2,
    /*              16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     *               c  t  a  G  A  T  T  C  A  G  C  A  T  T  G  N*/
                     1, 3, 0, 2, 0, 3, 3, 1, 0, 2, 1, 0, 3, 3, 2, 99,
    /*              32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
     *               G  C  T  G  A  C  C  T  C  A  A  T  G  t  a  g */
                     2, 1, 3, 2, 0, 1, 1, 3, 1, 0, 0, 3, 2, 3, 0, 2
  };
  //                   C  T  A  A  C  T  C  C  A  G  T  C  G  T  A  A  C  G  G
  char cons_base[] = { 1, 3, 0, 0, 1, 3, 1, 1, 0, 2, 3, 1, 2, 3, 0, 0, 1, 2, 2 };
  l_seqLib.sequence = p_left;

  // NextPtr, seqID, leftSeqPos, rightSeqPos,
  //    leftExtendable, rightExtendable, llen, rlen, orient
  struct coreAlignment l_coreAlign_0 = { NULL, 0, 12, 14, 1, 1, 0, 0, 0 };
  struct coreAlignment l_coreAlign_1 = { NULL, 1, 16, 18, 1, 1, 0, 0, 1 };
  struct coreAlignment l_coreAlign_2 = { NULL, 2, 44, 46, 1, 1, 0, 0, 0 };
  struct coreAlignment *l_cores[] = { &l_coreAlign_0, &l_coreAlign_1, &l_coreAlign_2 };

  /* RIGHT Extension Example1-reciprocal
   *                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
   *                  g  a  t  C  T  A  A  G  T  C  G  T  A  A  C  G*/
  char p_right[] = {  2, 0, 3, 1, 3, 0, 0, 2, 3, 1, 2, 3, 0, 0, 1, 2,
    /*               16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
     *                N  G  T  T  A  C  G  A  C  T  T  A  G  a  t  c */
                     99, 2, 3, 3, 0, 1, 2, 0, 1, 3, 3, 0, 2, 0, 3, 1,
    /*               32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
     *                g  a  t  G  T  A  A  C  T  C  C  A  G  T  C  G */
                      2, 0, 3, 2, 3, 0, 0, 1, 3, 1, 1, 0, 2, 3, 1, 2
  };
  // NextPtr, seqID, leftSeqPos, rightSeqPos,
  //    leftExtendable, rightExtendable, llen, rlen, orient
  struct coreAlignment r_coreAlign_0 = { NULL, 0, 1, 3, 1, 1, 0, 0, 0 };
  struct coreAlignment r_coreAlign_1 = { NULL, 1, 29, 31, 1, 1, 0, 0, 1 };
  struct coreAlignment r_coreAlign_2 = { NULL, 2, 33, 35, 1, 1, 0, 0, 0 };
  struct coreAlignment *r_cores[] = { &r_coreAlign_0, &r_coreAlign_1, &r_coreAlign_2 };

  r_seqLib.boundaries = (uint64_t *) malloc(4 * sizeof(uint64_t));
  r_seqLib.boundaries[0] = 16;
  r_seqLib.boundaries[1] = 32;
  r_seqLib.boundaries[2] = 48;
  r_seqLib.boundaries[3] = 0;
  r_seqLib.identifiers = predef_id;
  r_seqLib.sequence = p_right;
  r_seqLib.count = 3;
  r_seqLib.length = 48;
 

  if ( VERBOSE > 10 ) {
  printf("Boundary Conditions\n");
  printf("  Gap: ");
  for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    printf(" %d", l_score[1][0][offset + MAXOFFSET][1]);
  printf("\n");
  printf("  Sub: ");
  for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
    printf(" %d", l_score[1][0][offset + MAXOFFSET][0]);
  printf("\n\n");
  }

  int bad = 0;
  int final_right_gap[] = { -35, -29, -15, -1, -6, -11, -2, -7, -12, -17, -22 };
  int final_right_sub[] = { -33, -1, -41, -27, -11, 31, -8, -6, -7, -32, -43 };
  for ( n = 0; n < 3; n++ )
  {
    printf("n = %d\n", n);
    printf("Cores Left/Right reciprocal:\n");
    printCoreEdges(l_cores[n], &l_seqLib, 0);
    printCoreEdges(r_cores[n], &r_seqLib, 0);
    for (conspos = 0; conspos < 19; conspos++)
    {
      if ( VERBOSE > 10 )
        printf("row: %d,  cons_base=%d:\n", conspos, cons_base[conspos]);

      // Calculate from right to left
      bestscore =
        compute_nw_row(1,conspos, n, cons_base[conspos], r_cores[n], r_score,
                           lower_bound[n], upper_bound[n]-1, p_right,
                           &max_score_seq_idx, scoreParams,
                           MAXOFFSET, L, VERBOSE);
      if ( VERBOSE > 10 ) {
      printf("  RIGHT: bestscore = %d\n", bestscore);
      printf("         Gap: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", r_score[conspos % 2][n][offset + MAXOFFSET][1]);
      printf("\n");
      printf("         Sub: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", r_score[(conspos) % 2][n][offset + MAXOFFSET][0]);
      printf("\n");
      }

      bestscore =
        compute_nw_row(0,conspos, n, cons_base[conspos], l_cores[n], l_score,
                           lower_bound[n], upper_bound[n]-1, p_left,
                           &max_score_seq_idx, scoreParams,
                           MAXOFFSET, L, VERBOSE);

      if ( VERBOSE > 10 ){
      printf("  LEFT: bestscore = %d\n", bestscore);
      printf("         Gap: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", l_score[conspos % 2][n][offset + MAXOFFSET][1]);
      printf("\n");
      printf("         Sub: ");
      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
        printf(" %d", l_score[(conspos) % 2][n][offset + MAXOFFSET][0]);
      printf("\n");
      }

      for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
//        if ( ( l_score[conspos % 2][n][offset + MAXOFFSET][0] !=
//               r_score[conspos % 2][n][MAXOFFSET-offset][0] ) ||
//             ( l_score[conspos % 2][n][offset + MAXOFFSET][1] !=
//               r_score[conspos % 2][n][MAXOFFSET-offset][1] ) )
        if ( ( l_score[conspos % 2][n][MAXOFFSET-offset][0] !=
               r_score[conspos % 2][n][MAXOFFSET-offset][0] ) ||
             ( l_score[conspos % 2][n][MAXOFFSET-offset][1] !=
               r_score[conspos % 2][n][MAXOFFSET-offset][1] ) )
      {
        bad = 1;
        printf("Failed reciprocal test row=%d\n", n);
      }
      if ( n == 0 && conspos == 5 ) {
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
          if ( ( r_score[3%2][0][offset + MAXOFFSET][1] !=
             final_right_gap[offset+MAXOFFSET] ) ||
           ( r_score[3%2][0][offset + MAXOFFSET][0] !=
             final_right_sub[offset+MAXOFFSET] ) )
          {
            bad = 1;
            printf("Failed final comparison at %d\n", offset+MAXOFFSET);
          }
      }
    }
  }

  if (!bad)
  {
    printf("All tests passed!\n");
    exit(0);
  }else {
    exit(1);
  }

}

