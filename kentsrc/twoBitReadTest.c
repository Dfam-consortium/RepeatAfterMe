#include "common.h"
#include "linefile.h"
#include "dnaseq.h"
#include "twoBitNew.h"

/* parse a file containing a list of ranges for sequences in the
 * specified twoBit file. Specifications are one per line in the
 * form:
 *     seqName start end orient
 *
 * The fields are tab separated. Coordinates are zero-based half
 * open and the orientation is either "+" or "-".
 *
 *  e.g A 1000bp from "chr1" on the forward strand would be:
 *     chr1	1000	1001	+
 */
struct twoBitSpec *readBEDRanges(char *twoBitFile, char *bedFile)
{
struct lineFile *lf = lineFileOpen(bedFile, TRUE);
char *fields[4];
char *e;
struct twoBitSpec *spec;
struct twoBitSeqSpec *seq;
AllocVar(spec);
spec->fileName = cloneString(twoBitFile);
while (lineFileChopCharNext(lf, '\t', fields, 4)){
    AllocVar(seq);
    seq->name = cloneString(fields[0]);
    seq->start = strtol(fields[1], &e, 0);
    seq->end = strtol(fields[2], &e, 0);
    if ( *fields[3] == '-' )
      seq->strand = "-";
    else
      seq->strand = "+";
    slSafeAddHead(&spec->seqs, seq);
}
slReverse(&spec->seqs);
lineFileClose(&lf);
return spec;
}

void loadSequenceSubset(char *inName)
{
struct twoBitFile *tbf;
struct twoBitSpec *tbs;

tbs = readBEDRanges(inName,"test.seqlist");
if (tbs == NULL)
    errAbort("%s is not a twoBit file", inName);

tbf = twoBitOpen(tbs->fileName);
struct twoBitSeqSpec *s;
for (s = tbs->seqs; s != NULL; s = s->next){
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, s->name, s->start, s->end);
    // We want all uppercase
    toUpperN(seq->dna, seq->size);
    if ( s->strand[0] == '-' ){
      printf(">%s -\n", seq->name);
      reverseComplement(seq->dna, seq->size);
    }else{
      printf(">%s +\n", seq->name);
    }
    printf("seq->dna[0] = %c\n", (seq->dna)[0]);
    dnaSeqFree(&seq);
}

twoBitSpecFree(&tbs);
twoBitClose(&tbf);
}

int main(int argc, char *argv[])
{
dnaUtilOpen();
loadSequenceSubset("/usr/local/genomes/hg38/hg38.unmasked.2bit");
return 0;
}
