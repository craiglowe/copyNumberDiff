/*

faToGcStats.c
Written by: Craig Lowe

*/

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"
#include "fa.h"
#include "sqlNum.h"
#include "obscure.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"wiggle", OPTION_STRING},
	{NULL, 0}
};

char *optWiggle = NULL;

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	 "faToGcStats - count the gc content of windows\n"
	 "usage:\n"
	 "   faToGcStats in.fa noGap.bed windowLength out.stats\n"
	 "options:\n"
	 "     -wiggle  (filename)  Output number of G/Cs in window starting at each base\n"
	 "notes:\n"
	 );
}

/*---------------------------------------------------------------------------*/


struct hash *bedLoadNInHash(char *filename, int fields)
{
	struct bed *bed = NULL, *currList = NULL, *temp = NULL, *nextBed = NULL;
	struct hash *regionHash = newHash(6);
	struct bed *regions;

	regions = bedLoadNAll(filename, fields);
	slSort(&regions, bedCmp);
	currList = regions;
	for(bed = regions; bed != NULL; bed = nextBed)
	{
		nextBed = bed->next;
		if((bed->next == NULL) || (differentString(bed->chrom,bed->next->chrom)))
		{
			temp = bed->next;
			bed->next = NULL;
			hashAdd(regionHash, bed->chrom, currList);
			currList = temp;
		}
	}
	return(regionHash);
}


unsigned int reportGcCount(char *seq, unsigned int windowLength)
{
	unsigned int i = 0,  gcCount = 0;

	for(i=0; i<windowLength; i++)
	{
		if(seq[i] == 'c' || seq[i] == 'C' || seq[i] == 'g' || seq[i] == 'G'){gcCount++;}
		else if(seq[i] == 'a' || seq[i] == 'A' || seq[i] == 't' || seq[i] == 'T'){}
		else{errAbort("Error: ungapped regions should be only acgtACGT bases: %c\n", seq[i]);}
	}
        return(gcCount);
}


void outputGcStatsWiggle(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int windowLength, char *outFilename)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	unsigned int i = 0;
	char *currWindow = NULL;
	FILE *fout = mustOpen(outFilename, "w");

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			fprintf(fout, "fixedStep chrom=%s start=%u step=1\n", currRegion->chrom, currRegion->chromStart+1);
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - windowLength; i++)
			{
				currWindow = &(currSeq->dna[i]);
				fprintf(fout, "%u\n", reportGcCount(currWindow, windowLength));
			}
		}
	}
	carefulClose(&fout);
}


void incrementGcBin(char *seq, unsigned int windowLength, double *gcBins)
{
	unsigned int i = 0,  gcCount = 0;

	for(i=0; i<windowLength; i++)
	{
		if(seq[i] == 'c' || seq[i] == 'C' || seq[i] == 'g' || seq[i] == 'G'){gcCount++;}
		else if(seq[i] == 'a' || seq[i] == 'A' || seq[i] == 't' || seq[i] == 'T'){}
		else{errAbort("Error: ungapped regions should be only acgtACGT bases: %c\n", seq[i]);}
	}
	gcBins[gcCount]++;
}


double gatherGcStats(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int windowLength, double *gcBins)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	unsigned int i = 0;
	char *currWindow = NULL;
	double totalWindows = 0;

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - windowLength; i++)
			{
				currWindow = &(currSeq->dna[i]);
				incrementGcBin(currWindow, windowLength, gcBins);
				totalWindows++;
			}
		}
	}

	return(totalWindows);
}


void writeOutput(double *gcBins, unsigned int windowSize, double totalWindows, char *outFilename)
{
	unsigned int i = 0;
	FILE *fout = mustOpen(outFilename, "w");
	double percentGc = 0, percentTotal = 0;

	for(i=0; i<=windowSize; i++)
	{
		percentGc = (double)i/(double)windowSize;
		percentTotal = gcBins[i]/totalWindows;
		fprintf(fout, "%u %f %f\n", i, percentGc, percentTotal);
	}
	carefulClose(&fout);
}


/*---------------------------------------------------------------------------*/

void faToGcStats(char *inFaFile, char *noGapBedFile, unsigned int windowSize, char *outFile)
{
	struct dnaSeq *seqList = NULL;
	struct hash *noGapHash = NULL;
	double totalWindows = 0;
	double *gcBins = NULL;
	AllocArray(gcBins, windowSize+1);

	verbose(2, "reading in fa file\n");
	seqList = faReadAllDna(inFaFile);

	verbose(2, "reading in no gap file\n");
	noGapHash = bedLoadNInHash(noGapBedFile, 3);

	verbose(2, "finding the kmers that appear more than once\n");
	totalWindows = gatherGcStats(seqList, noGapHash, windowSize, gcBins);
	
	verbose(2, "writing uniqueness of each position\n");
	writeOutput(gcBins, windowSize, totalWindows, outFile);

	if(optWiggle != NULL)
	{
		verbose(2, "writing wiggle output\n");
		outputGcStatsWiggle(seqList, noGapHash, windowSize, optWiggle);
	}

	verbose(2, "done\n");
}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
/* Process command line. */
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 5){usage();}

	optWiggle = optionVal("wiggle", optWiggle);

	faToGcStats(argv[1],argv[2],sqlUnsigned(argv[3]),argv[4]);

	return(0);
}

