#ifndef __REPORT_H__
#define __REPORT_H__


void
printExtensionRegion( struct coreAlignment *currCore, int direction, int row, int capped, struct sequenceLibrary *seqLib );

void
printCoreEdges(struct coreAlignment *coreAlign, struct sequenceLibrary *seqLib, char omitBlanks);

#endif
