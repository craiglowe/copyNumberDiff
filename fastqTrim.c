/*
 * fastqTrim.c
 */

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "fq.h"
#include <inttypes.h>

static struct optionSpec options[] = {
	{"minLength", OPTION_INT},
	{"maxLength", OPTION_INT},
	{"cutoff", OPTION_INT},
	{NULL, 0},
};

int optMaxLength = 0;
int optMinLength = 30;
int optCutoff = 3;

void usage()
{
	errAbort(
	"fastqTrim - Trims fastq entries at the 3-prime end based on quality and length\n"
	"usage:\n"
	"	fastqTrim input.fq output.fq\n"
	"options:\n"
	"	-maxLength (0;disabled) trim read to this length if longer\n"
	"	-minLength (36)         remove entry if shorter than this after trimming\n"
	"	-cutoff    (3)          the lowest passing score\n"
	);
}


unsigned int fqBestQualIntervalStart(struct fq *record, unsigned int window)
{
	unsigned int bestScore = 0, bestIdx = 0, x = 0, y = 0, qualitySum = 0, origLen = 0;

	origLen = strlen(record->dna);
	if(window > origLen){errAbort("Error: asking fqBestInterval to use a window longer than the starting length\n");}
	if(window == origLen){return(0);}

	for(x=0; x < origLen - window; x++)
	{
		for(y=0, qualitySum=0; y < window; y++)
		{
			qualitySum += record->quality[x+y];
		}
		if(qualitySum > bestScore)
		{
			bestScore = qualitySum;
			bestIdx = x;
		}
	}
	return(bestIdx);
}


void fqClipSeqAndQual(struct fq *record, unsigned int start, unsigned int length)
{
	unsigned int x = 0;

	for(x=0; x<length; x++)
	{
		record->dna[x] = record->dna[x+start];
		record->quality[x] = record->quality[x+start];
	}
	record->dna[length] = '\0';
	record->quality[length] = '\0';
}


boolean trim(struct fq *in, int minLength, int maxLength, int cutoff)
{
	int newLen = 0, base = 33, start = 0, end = 0, runningSum = 0, x = 0, len = 0, firstPositiveIdx = 0, maxSum = 0;

	len = strlen(in->dna);

	for(x=len-1; x>=0; x--)
	{
		runningSum += in->quality[x] - base - cutoff;
		if(runningSum < 0)
		{
			runningSum = 0;
			firstPositiveIdx = 0;
		}
		else if (runningSum > 0)
		{
			if(firstPositiveIdx == 0)
			{
				firstPositiveIdx = x;
			}
			if(runningSum > maxSum)
			{
				start = x;
				end = firstPositiveIdx;
			}
		}
	}
	verbose(3, "%s %d %d\n", in->quality, start, end);
	if(end - start + 1 > maxLength)
	{
		start = fqBestQualIntervalStart(in, maxLength);
		end = start + maxLength - 1;
	}
	if(start != 0){fqClipSeqAndQual(in, start, end - start + 1);}
	in->dna[end+1] = '\0';
	in->quality[end+1] = '\0';
	newLen = strlen(in->dna);
	verbose(3, "%s %d %d\n", in->quality, start, end);
	if(newLen < minLength){return(FALSE);}
	return(TRUE);
}

void trimFastqFile(char *inFilename, int minLength, int maxLength, int cutoff, char *outFilename)
{
	FILE *f = mustOpen(outFilename, "w");
	struct lineFile *lf = lineFileOpen(inFilename, TRUE);
	struct fq *record = NULL;
	unsigned int progress = 0;
	while ((record = fqReadNext(lf)) != NULL)
	{
		if(trim(record, minLength, maxLength, cutoff))
		{
			fqWriteNext(record, f);
		}
		fqFree(&record);
		progress++;
		if(progress % 1000000 == 0){verbose(2, "Processed %u records\n", progress);}
	}
	carefulClose(&f);
}

int main(int argc, char *argv[])
{
	optionInit(&argc, argv, options);
	optMinLength = optionInt("minLength", optMinLength);
	optMaxLength = optionInt("maxLength", optMaxLength);
	optCutoff = optionInt("cutoff", optCutoff);

	if(argc != 3){usage();}
	
	trimFastqFile(argv[1], optMinLength, optMaxLength, optCutoff, argv[2]);

	return(0);
}

