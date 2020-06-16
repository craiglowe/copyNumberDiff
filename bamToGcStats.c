/*

bamToGcStats.c
Written by Craig Lowe

*/

#include "common.h"
#include "options.h"
#include "memalloc.h"
#include "hash.h"
#include "sqlNum.h"
#include "linefile.h"
#include "htslib/sam.h"
#include "bamFile.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"minQual", OPTION_INT},
	{NULL, 0}
};

unsigned int optMinQual = 0;

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	"bamToGcStats - Computer the GC statitics for reads in bam file\n"
	"usage:\n"
	"   bamToGcStats in.bam readLength out.stats\n"
	"options:\n"
	);
}

/*---------------------------------------------------------------------------*/


double addReadCounts(char *filename, double *bins, uint32_t readLength)
{
	//samfile_t *samFp = NULL;
	samfile_t *samFp = bamMustOpenLocal(filename, "rb", NULL);
	bam_header_t *head = sam_hdr_read(samFp);
	if (head == NULL) {
		errAbort("Aborting ... bad BAM header in %s", filename);
	}
	//samFp = samopen(filename, "rb", 0);
	bam1_t* b = bam_init1();
	bam1_core_t *c = NULL;
	unsigned int i = 0;
	uint32_t base = 0;
	uint32_t at = 0, gc = 0, n = 0;
	double used = 0, tooShort = 0, badMapping = 0, hadN = 0;
	uint8_t *seq = NULL;

	while(sam_read1(samFp, head, b) >= 0)
	{
		c = &b->core;
		if((c->qual >= optMinQual) && (!(c->flag&(BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FUNMAP))))
		{
			if(readLength <= c->l_qseq)
			{
				seq = bam1_seq(b);
				for(i = 0, at = 0, gc = 0, n = 0;  i < readLength;  i++)
				{
					base = bam1_seqi(seq, i);
					if(base == 1 || base == 8){at++;}
					else if(base == 2 || base == 4){gc++;}
					else if(base == 15){n++;}
					else{errAbort("Base: %u not recognized\n", base);}
				}
				if(n == 0)
				{
					used++;
					bins[gc]++;
				}
				else{hadN++;}
			}
			else{tooShort++;}
		}
		else{badMapping++;}
	}
	samclose(samFp);
	bam_destroy1(b);
	verbose(3, "used:%e hadN:%e tooShort:%e badMapping:%e\n", used, hadN, tooShort, badMapping);
	return(used);
}


void printGcBins(double *gcBins, int32_t readLength, double totalReads, char *outFilename)
{
	uint32_t i = 0;
	double percentGc = 0, percentTotal = 0;
	FILE *fout = mustOpen(outFilename, "w");
	for(i=0; i<=readLength; i++)
	{
		percentGc = (double)i/(double)readLength;
		percentTotal = gcBins[i]/totalReads;
		fprintf(fout, "%u %f %e\n", i, percentGc, percentTotal);
	}
	carefulClose(&fout);
}


/*---------------------------------------------------------------------------*/


void bamToGcStats(char *bamFilename, uint32_t readLength, char *outFilename)
{
	double totalReads = 0;
	double *gcBins = NULL;
	AllocArray(gcBins, readLength+1);

	verbose(2, "Read bam file\n");
	totalReads = addReadCounts(bamFilename, gcBins, readLength);

	verbose(2, "Writing to the output file\n");
	printGcBins(gcBins, readLength, totalReads, outFilename);
}


/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 4){usage();}

	optMinQual = optionInt("minQual", optMinQual);

	bamToGcStats(argv[1], sqlUnsigned(argv[2]), argv[3]);

	return(0);
}

