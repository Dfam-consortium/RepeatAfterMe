/*
 * report.c
 *
 * Routines for reporting on program status
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "bnw_extend.h"
#include "score_system.h"
#include "sequence.h"


// RMH: Corrected for core left/rightSeqPos change and TESTED
void
printExtensionRegion( struct coreAlignment *currCore, int direction, int row, int capped, struct sequenceLibrary *seqLib )
{
  int j = 0;
  char coreSeq[30];
  char extSeq[20];
  char base[2];
  uint64_t lastAligned = 0;
  base[1] = '\0';
  coreSeq[0] = '\0';
  extSeq[0] = '\0';

  uint64_t lowerBound = 0;
  if ( currCore->seqIdx > 0 )
     lowerBound = seqLib->boundaries[currCore->seqIdx - 1];
  uint64_t upperBound = seqLib->boundaries[currCore->seqIdx];


  printf("C_EDGE: sIdx=%d, orient=%d", currCore->seqIdx, currCore->orient );
  if ( direction ) {
    // Right Extension
    if ( currCore->orient ) {
      // Right + reverse
      lastAligned = currCore->rightSeqPos - row - 1;
      printf(", pos=%ld: ",lastAligned - lowerBound + 1);
      for (j = 10; j >= 1; j--) {
        if ( lastAligned + j == upperBound ) {
          strcat(coreSeq,"*");
        }else if ( lastAligned + j > upperBound ) {
          strcat(coreSeq," ");
        }else {
          base[0] = num_to_char(compl(seqLib->sequence[lastAligned+j]));
          strcat(coreSeq, base);
        }
      }
      for (j = -1; j >= -10; j--) {
        if ( lastAligned + j == lowerBound ) {
          strcat(extSeq,"*");
        }else if ( lastAligned + j < lowerBound ) {
          strcat(extSeq," ");
        }else {
          base[0] = num_to_char(compl(seqLib->sequence[lastAligned+j]));
          strcat(extSeq, base);
        }
      }
      base[0] = num_to_char(compl(seqLib->sequence[lastAligned]));
    }else {
      // Right + forward
      lastAligned = currCore->rightSeqPos + row + 1;
      printf(", pos=%ld: ",lastAligned - lowerBound + 1);
      for (j = -10; j < 0; j++) {
        if ( lastAligned + j == lowerBound ) {
          strcat(coreSeq,"*");
        }else if ( lastAligned + j < lowerBound ) {
          strcat(coreSeq," ");
        }else {
          base[0] = num_to_char(seqLib->sequence[lastAligned+j]);
          strcat(coreSeq, base);
        }
      }
      for (j = 1; j <= 10; j++) {
        if ( lastAligned + j == upperBound ) {
          strcat(extSeq,"*");
        }else if ( lastAligned + j > upperBound ) {
          strcat(extSeq," ");
        }else {
          base[0] = num_to_char(seqLib->sequence[lastAligned+j]);
          strcat(extSeq, base);
        }
      }
      base[0] = num_to_char(seqLib->sequence[lastAligned]);
    }
    printf("%s] {%c} %s\n", coreSeq, base[0], extSeq);
  }else {
    // Left Extension
    if ( currCore->orient ) {
      // Left + reverse
      lastAligned = currCore->leftSeqPos + row + 1;
      printf(", pos=%ld: ",lastAligned - lowerBound + 1);
      for (j = -1; j >= -10; j--) {
        if ( lastAligned + j == lowerBound ) {
          strcat(coreSeq,"*");
        }else if ( lastAligned + j < lowerBound ) {
          strcat(coreSeq," ");
        }else {
          base[0] = num_to_char(compl(seqLib->sequence[lastAligned+j]));
          strcat(coreSeq, base);
        }
      }
      for (j = 9; j > 0; j--) {
        if ( lastAligned + j == upperBound ) {
          strcat(extSeq,"*");
        }else if ( lastAligned + j > upperBound ) {
          strcat(extSeq," ");
        }else {
          base[0] = num_to_char(compl(seqLib->sequence[lastAligned+j]));
          strcat(extSeq, base);
        }
      }
      base[0] = num_to_char(compl(seqLib->sequence[lastAligned]));
    }else {
      // Left + forward
      lastAligned = currCore->leftSeqPos - row - 1;
      printf(", pos=%ld: ",lastAligned - lowerBound + 1);
      for (j = 1; j <= 10; j++) {
        if ( lastAligned + j == upperBound ) {
          strcat(coreSeq,"*");
        }else if ( lastAligned + j > upperBound ) {
          strcat(coreSeq," ");
        }else {
          base[0] = num_to_char(seqLib->sequence[lastAligned+j]);
          strcat(coreSeq, base);
        }
      }
      for (j = -10; j < 0; j++) {
        if ( lastAligned + j == lowerBound ) {
          strcat(extSeq,"*");
        }else if ( lastAligned + j < lowerBound ) {
          strcat(extSeq," ");
        }else {
          base[0] = num_to_char(seqLib->sequence[lastAligned+j]);
          strcat(extSeq, base);
        }
      }
      base[0] = num_to_char(seqLib->sequence[lastAligned]);
    }
    printf("%s {%c} [%s\n", extSeq, base[0], coreSeq);
 }
  //printf("%11s [%c] %-11s\n", leftExt, coreSeq[0], rightExt);
}

/*
 *  printCoreEdges
 *   Print the core datastructure to the screen
 *
 *   coreAlign:  The core datastructure
 *   seqLib:     The sequence library
 *   omitBlanks: Omit cores where the seqIdx, leftSeqPos, and rightSeqPos are all 0
 *                 (currently boolean 1 or 0)
 *   debug:      Print additional internal-sequence cordinates and boundaries
 *                 for debugging (currently 1 or 0 ).
 */
void
printCoreEdges(struct coreAlignment *coreAlign, struct sequenceLibrary *seqLib,
               char omitBlanks, char debug) {
  struct coreAlignment *currCore;
  int j = 0;
  int n = 0;
  int i = 0;
  int disp = 0;
  char leftExt[20];
  char coreSeq[30];
  char rightExt[20];
  char base[2];
  base[1] = '\0';
  char idBuff[51];
  idBuff[0] = '\0';
  char posBuffer[50];
  posBuffer[0] = '\0';
  char sArrBuffer[50];
  sArrBuffer[0] = '\0';

  // Get a sense for the magnitude of the individual sequence indexes
  // This will be used to format column widths
  uint64_t maxPos = 0;
  int maxIDLen = 0;

  int maxIdx = 1;
  for (currCore = coreAlign; currCore != NULL; currCore = currCore->next)
  {
     // Check if this is a subsequence
     uint64_t subseq_offset = 0;
     if ( seqLib->offsets != NULL && seqLib->offsets[currCore->seqIdx] > 0 )
       subseq_offset = seqLib->offsets[currCore->seqIdx];

     if ( omitBlanks == 1 && currCore->seqIdx == 0 && currCore->leftSeqPos == 0 &&
          currCore->rightSeqPos == 0 )
       break;

     uint64_t seqLowerBound = 0;
     if ( currCore->seqIdx > 0 )
       seqLowerBound = seqLib->boundaries[currCore->seqIdx - 1];

     if ( strlen(seqLib->identifiers[currCore->seqIdx]) > maxIDLen )
       maxIDLen = strlen(seqLib->identifiers[currCore->seqIdx]);

     if ( subseq_offset + currCore->leftSeqPos - seqLowerBound + 1 > maxPos )
       maxPos = subseq_offset + currCore->leftSeqPos - seqLowerBound + 1;
     if ( subseq_offset + currCore->rightSeqPos - seqLowerBound + 1 > maxPos )
       maxPos = subseq_offset + currCore->rightSeqPos - seqLowerBound + 1;
     maxIdx++;
  }
  // Maximum ID we will display is 50 characters.  After that we will truncate at 47 and add "..."
  if ( maxIDLen > 50 )
    maxIDLen = 50;
  if ( maxIDLen < 5 )
    maxIDLen = 5;

  // Sanity check the genomic coordinate lengths
  posBuffer[0] = '\0';
  if ( snprintf(posBuffer,50,"%ld",maxPos) > 50 ) {
    printf("printCoreEdges: warning this sequence position is absurdly high: %ld\n", maxPos);
  }
  int maxPosLen = strlen(posBuffer);
  posBuffer[0] = '\0';
  if ( snprintf(posBuffer,50,"%ld",seqLib->length) > 50 ) {
    printf("printCoreEdges: warning this sequence position is absurdly high: %ld\n", seqLib->length);
  }
  // There is a minimum so width so that the heading prints correctly
  if ( maxPosLen*2+1 < 9 )
    maxPosLen = 9;

  int maxSArrayLen = strlen(posBuffer);
  int maxIdxLen = floor(log10(maxIdx))+1;
  if ( maxIdxLen < 4 )
    maxIdxLen = 4;

  if ( debug == 1 ){
    printf("%-*s %-*s %-*s %-*s %-*s  Left-Flank           Core           Right-Flank  %-*s Hard-Bounds Soft-Bounds BoundFlags(Lower/Upper)\n",
         maxIdxLen, "Seq", maxIDLen, "Ident", maxPosLen*2+1,
         "BED-range", 6, "Orient", 4 , "L/R?", maxSArrayLen, "seq[]-core");
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }else {
    printf("%-*s %-*s %-*s %-*s %-*s  Left-Flank           Core           Right-Flank\n",
         maxIdxLen, "Seq", maxIDLen, "Ident", maxPosLen*2+1,
         "BED-range", 6, "Orient", 4 , "L/R?");
    printf("------------------------------------------------------------------------------------\n");
  }
  //
  // leftSeqPos/rightSeqPos are the core coordinates (0-based,fully closed ...ie. inclusive)
  // NOTE: leftSeqPos will be greater than rightSeqPos when the orientation of the word is reversed.
  //
  for (currCore = coreAlign; currCore != NULL; currCore = currCore->next)
  {
     if ( omitBlanks == 1 && currCore->seqIdx == 0 && currCore->leftSeqPos == 0 &&
          currCore->rightSeqPos == 0 )
       break;

     // Check if this is a subsequence
     uint64_t subseq_offset = 0;
     if ( seqLib->offsets != NULL && seqLib->offsets[currCore->seqIdx] > 0 )
       subseq_offset = seqLib->offsets[currCore->seqIdx];

    // Stats for sequence/core
    uint64_t lowerBound = 0;
    if (currCore->seqIdx > 0)
      lowerBound = seqLib->boundaries[currCore->seqIdx - 1];
    uint64_t upperBound = seqLib->boundaries[currCore->seqIdx] - 1;
    int coreWidth = abs(currCore->rightSeqPos - currCore->leftSeqPos) + 1;

    // Build out ID up to maxIDLen
    strncpy(idBuff, seqLib->identifiers[currCore->seqIdx], maxIDLen);
    idBuff[maxIDLen] = '\0';
    int idLen = strlen(idBuff);
    for ( i = 0; i < maxIDLen - idLen; i++ )
      strcat(idBuff," ");
    if ( strlen(seqLib->identifiers[currCore->seqIdx]) > 50 )
      for ( i = 0; i < 3; i++ )
        idBuff[maxIDLen-i-1] = '.';

    // Put core alignment coordinates back into BED format
    char tmpBuff[50];
    tmpBuff[0] = '\0';
    if (currCore->orient)
      sprintf(tmpBuff,"%ld-%ld",subseq_offset+currCore->rightSeqPos-lowerBound,
              subseq_offset+currCore->leftSeqPos-lowerBound+1);
    else
      sprintf(tmpBuff,"%ld-%ld",subseq_offset+currCore->leftSeqPos-lowerBound, 
              subseq_offset+currCore->rightSeqPos-lowerBound+1);
    int numPads = ((maxPosLen*2)+1) - strlen(tmpBuff);
    posBuffer[0] = '\0';
    for ( i = 0; i < numPads; i++ )
      strcat(posBuffer, " ");
    strcat(posBuffer, tmpBuff);

    // Format the internal sequence* coordinates
    tmpBuff[0] = '\0';
    sprintf(tmpBuff,"%ld-%ld",currCore->leftSeqPos, currCore->rightSeqPos);
    numPads = ((maxSArrayLen*2)+1) - strlen(tmpBuff);
    sArrBuffer[0] = '\0';
    for ( i = 0; i < numPads; i++ )
      strcat(sArrBuffer, " ");
    strcat(sArrBuffer, tmpBuff);

    // For cores up to 19bp display the full core sequence with up to 10bp
    // flanking:
    //
    //       AATACAAAAA [          ATT           ] AGCCGGGCG
    //
    // or if there is fewer than bp of flanking display the end with a '*':
    //
    //          *CAAAAA [          ATT           ] AGCCGGGCG
    //
    // If the core is 20bp or larger display only 10bp core edges as:
    //
    //       CTCTACTAAA [AATACAAAAA....AGCCGGGCGT] GGTGGCGGG
    //
    if (currCore->orient)
    {
      //
      // Reverse Strand Cores
      //    NOTE: in the reverse strand case the leftPos is the higher
      //          numerical coordinate:
      //               e.g: leftSeqPos = 32
      //                    rightSeqPos = 10


      // Left Flanking Region
      disp = 10;
      if (currCore->leftSeqPos + 1 + disp > upperBound)
        disp = upperBound - currCore->leftSeqPos - 1;

      leftExt[0] = '\0';
      if (disp < 10)
        strcat(leftExt, "*");
      for (j = disp; j > 0; j--)
      {
        // Possibly + 1
        base[0] = num_to_char(compl(seqLib->sequence[currCore->leftSeqPos + j]));
        strcat(leftExt, base);
      }

      // Core Region
      coreSeq[0] = '\0';
      if (coreWidth < 20)
      {
        // Short core case e.g. "[        CCATTCA        ]"
        for (j = 0; j < ((24 - coreWidth) / 2); j++)
          strcat(coreSeq, " ");
        for (j = 0; j > -coreWidth; j--)
        {
          base[0] = num_to_char(compl(seqLib->sequence[currCore->leftSeqPos + j]));
          strcat(coreSeq, base);
        }
      }
      else
      {
        // Long core case e.g. "[ACGGAATAGC....ACGTAGATTA]"
        for (j = 0; j > -10; j--)
        {
          base[0] = num_to_char(compl(seqLib->sequence[currCore->leftSeqPos + j]));
          strcat(coreSeq, base);
        }
        strcat(coreSeq, "....");
        for (j = 9; j >= 0; j--)
        {
          base[0] = num_to_char(compl(seqLib->sequence[currCore->rightSeqPos + j]));
          strcat(coreSeq, base);
        }
      }

      // Right Flanking Region
      disp = 10;
      if (currCore->rightSeqPos < disp)
        disp = currCore->rightSeqPos;
      else if (currCore->rightSeqPos - disp < lowerBound)
        disp = currCore->rightSeqPos - lowerBound;

      rightExt[0] = '\0';
      for (j = 0; j > -disp; j--)
      {
        base[0] = num_to_char(compl(seqLib->sequence[currCore->rightSeqPos + j - 1]));
        strcat(rightExt, base);
      }
      if (disp < 10)
        strcat(rightExt, "*");

      printf("%-*d %s %s -      %d/%d ", maxIdxLen, n,
             idBuff, posBuffer, currCore->leftExtendable, currCore->rightExtendable);
    }
    else
    {
      //
      // Forward Strand Cores
      //

      // Left Flanking Region
      disp = 10;
      if (currCore->leftSeqPos < disp)
        disp = currCore->leftSeqPos;
      else if (currCore->leftSeqPos - disp < lowerBound)
        disp = currCore->leftSeqPos - lowerBound;

      leftExt[0] = '\0';
      if (disp < 10)
        strcat(leftExt, "*");
      for (j = -disp; j <= -1; j++)
      {
        base[0] = num_to_char(seqLib->sequence[currCore->leftSeqPos + j]);
        strcat(leftExt, base);
      }

      // Core Region
      coreSeq[0] = '\0';
      if (coreWidth < 20)
      {
        // Short core case e.g. "[        CCATTCA        ]"
        for (j = 0; j < ((24 - coreWidth) / 2); j++)
          strcat(coreSeq, " ");
        for (j = 0; j < coreWidth; j++)
        {
          base[0] = num_to_char(seqLib->sequence[currCore->leftSeqPos + j]);
          strcat(coreSeq, base);
        }
      }
      else
      {
        // Long core case e.g. "[ACGGAATAGC....ACGTAGATTA]"
        for (j = 0; j < 10; j++)
        {
          base[0] = num_to_char(seqLib->sequence[currCore->leftSeqPos + j]);
          strcat(coreSeq, base);
        }
        strcat(coreSeq, "....");
        for (j = -9; j <= 0; j++)
        {
          base[0] = num_to_char(seqLib->sequence[currCore->rightSeqPos + j]);
          strcat(coreSeq, base);
        }
      }

      // Right Flanking Region
      disp = 10;
      if (currCore->rightSeqPos + disp > upperBound)
        disp = upperBound - currCore->rightSeqPos;

      rightExt[0] = '\0';
      for (j = 1; j <= disp; j++)
      {
        base[0] = num_to_char(seqLib->sequence[currCore->rightSeqPos + j]);
        strcat(rightExt, base);
      }
      if (disp < 10)
        strcat(rightExt, "*");

      printf("%-*d %s %s +      %d/%d ", maxIdxLen, n,
             idBuff, posBuffer, currCore->leftExtendable, currCore->rightExtendable);
    }
    if (currCore->leftExtendable)
      printf("%11s", leftExt);
    else
      printf("           ");
    printf(" [%-24s] ", coreSeq);
    if (currCore->rightExtendable)
      printf("%-11s", rightExt);
    else
      printf("           ");

    if ( debug == 1 )
    {
      printf(" %s %ld-%ld %ld-%ld", sArrBuffer, lowerBound, upperBound, currCore->lowerSeqBound, currCore->upperSeqBound);
      switch ( currCore->lowerSeqBoundFlag )  {
        case L_BOUNDARY:
          printf(" L_BOUNDARY");
          break;
        case SEQ_BOUNDARY:
          printf(" SEQ_BOUNDARY");
          break;
        case CORE_BOUNDARY:
          printf(" CORE_BOUNDARY");
          break;
        case EXT_BOUNDARY:
          printf(" EXT_BOUNDARY");
          break;
      }
      switch ( currCore->upperSeqBoundFlag )  {
        case L_BOUNDARY:
          printf("/L_BOUNDARY");
          break;
        case SEQ_BOUNDARY:
          printf("/SEQ_BOUNDARY");
          break;
        case CORE_BOUNDARY:
          printf("/CORE_BOUNDARY");
          break;
        case EXT_BOUNDARY:
          printf("/EXT_BOUNDARY");
          break;
      }
    }
    printf("\n");
    n++;
  }
}
