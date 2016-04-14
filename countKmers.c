/*

countKmers.c
Written by Craig Lowe

*/

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"
#include "fa.h"
#include "sqlNum.h"
#include "obscure.h"
#include "murmur3.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
{
	{NULL, 0}
};

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	"countKmers - count the number of times each k-mer appears\n"
	"usage:\n"
	"   countKmers kmerLength in.fa noGap.bed kmers.wig\n"
	"options:\n"
	"notes:\n"
	);
}

/*---------------------------------------------------------------------------*/


struct kmerData
{
	unsigned int length;
	unsigned int  totalCount;
	struct hash *counts;
};


struct bloomFilter
{
	unsigned int *table;
	unsigned int tabLen;
	unsigned int numFunc;
};


struct kmerData *createKmerData(unsigned int kmerLength)
{
	struct kmerData *answer = NULL;
	AllocVar(answer);
	answer->length = kmerLength;
	answer->totalCount = 0;
	answer->counts = newHash(28);
	return(answer);
}


struct bloomFilter *createBloomFilter()
{
	struct bloomFilter *answer = NULL;
	AllocVar(answer);
	answer->tabLen = 1e9;
	answer->numFunc = 5;
	AllocArray(answer->table, answer->tabLen);
	return(answer);
}


void freeBloom(struct bloomFilter **pEl)
{
	struct bloomFilter *el = NULL;
	el = *pEl;
	if(el != NULL)
	{
		freeMem(el->table);
		freez(pEl);
	}
}


void putInBloom(char *seq, struct bloomFilter *bloom)
{
	uint64_t hash[2];
	unsigned int pos = 0, flag = 0, seed = 17, randNum = 0, i = 0, index = 0;
	
	for(i=0; i<bloom->numFunc; i++, seed++)
	{
		MurmurHash3_x64_128(seq, strlen(seq), seed, hash);
		randNum = hash[0] % (64 * bloom->tabLen);
		index = randNum / 64;
		pos = randNum % 64;
		flag = 1;
		flag = flag << pos;
		bloom->table[index] = bloom->table[index] | flag;
	}
}


int checkBloom(char *seq, struct bloomFilter *bloom)
{
	uint64_t hash[2];
	unsigned int pos = 0, flag = 0, seed = 17, randNum = 0, i = 0, index = 0;

	for(i=0; i<bloom->numFunc; i++, seed++)
	{
		MurmurHash3_x64_128(seq, strlen(seq), seed, hash);
		randNum = hash[0] % (64 * bloom->tabLen);
		index = randNum / 64;
		pos = randNum % 64;
		flag = 1;
		flag = flag << pos;
		if(!(bloom->table[index] & flag))
		{
			return(0);
		}
	}
	return(1);
}


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


struct hash *createBasewiseHash(struct hash *noGapHash)
{
	struct bed *headBed = NULL, *currBed = NULL;
	unsigned int max = 0;
	unsigned int *bases = NULL;
	struct hash *answer = newHash(8);

	struct hashEl *hel = NULL, *helList = hashElListHash(noGapHash);
	for (hel = helList; hel != NULL; hel = hel->next)
	{
		headBed = hel->val;
		max = 0;
		for(currBed = headBed, max=0; currBed != NULL; currBed=currBed->next)
		{
			if(currBed->chromEnd > max){max = currBed->chromEnd;}
		}
		AllocArray(bases, max);
		hashAdd(answer, headBed->chrom, bases);
	}
	hashElFreeList(&helList);
	return(answer);
}



void revComp(char *src, char *dest, unsigned int len)
{
	unsigned int i=0, j=0;

	dest[len] = '\0';
	for(i=0, j=len-1; i<len; i++, j--)
	{
		if(src[i] == 'A' || src[i] == 'a'){dest[j] = 't';}
		else if(src[i] == 'C' || src[i] == 'c'){dest[j] = 'g';}
		else if(src[i] == 'G' || src[i] == 'g'){dest[j] = 'c';}
		else if(src[i] == 'T' || src[i] == 't'){dest[j] = 'a';}
		else{errAbort("Error: unexpected character (%c)\n", src[i]);}
	}
}

void recordIndividualKmerSecondPass(char *currKmer, struct kmerData *answer)
{
	struct hashEl *hel = NULL;
	answer->totalCount++;
	hel = hashLookup(answer->counts, currKmer);
	if(hel != NULL)
	{
		unsigned int *counter = NULL;
		counter = hel->val;
		(*counter)++;
	}
}


void recordIndividualKmerFirstPass(char *currKmer, struct bloomFilter *bloom, struct kmerData *answer)
{
	struct hashEl *hel = NULL;

	if(checkBloom(currKmer, bloom))
	{
		hel = hashLookup(answer->counts, currKmer);
		if(hel == NULL)
		{
			unsigned int *temp = NULL;
			AllocVar(temp);
			*temp = 0;
			hashAdd(answer->counts, currKmer, temp);
		}
	}
	else
	{
		putInBloom(currKmer, bloom);
	}
}


unsigned int getKmerCount(char *currKmer, struct kmerData *kmers)
{
	struct hashEl *hel = hashLookup(kmers->counts, currKmer);
	if(hel == NULL){return(1);}
	else
	{
		unsigned int *count = hel->val;
		return(*count);
	}
}


struct kmerData *populateKmerHash(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	int i = 0;
	char *currKmer = NULL;
	struct bloomFilter *bloom = createBloomFilter();
	char revCompKmer[kmerLength+1];
	char temp;
	struct kmerData *answer = createKmerData(kmerLength);

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		verbose(3, " currently on %s\n", currSeq->name);
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - kmerLength; i++)
			{
				temp = currSeq->dna[i+kmerLength];
				currSeq->dna[i+kmerLength] = '\0';
				currKmer = &(currSeq->dna[i]);
				recordIndividualKmerFirstPass(currKmer, bloom, answer);
				revComp(currKmer, revCompKmer, kmerLength);
				recordIndividualKmerFirstPass(revCompKmer, bloom, answer);
				currSeq->dna[i+kmerLength] = temp;
			}
		}
	}
	freeBloom(&bloom);
	return(answer);
}


void recordKmers(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength, struct kmerData *answer)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	int i = 0;
	char *currKmer = NULL;
	char revCompKmer[kmerLength+1];
	char temp;

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		verbose(3, " currently on %s\n", currSeq->name);
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - kmerLength; i++)
			{
				temp = currSeq->dna[i+kmerLength];
				currSeq->dna[i+kmerLength] = '\0';
				currKmer = &(currSeq->dna[i]);
				recordIndividualKmerSecondPass(currKmer, answer);
				revComp(currKmer, revCompKmer, kmerLength);
				recordIndividualKmerSecondPass(revCompKmer, answer);
				currSeq->dna[i+kmerLength] = temp;
			}
		}
	}
}


void calcUniqueness(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength, struct kmerData *kmers, struct hash *mapsHash)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	char *currKmer = NULL;
	char temp;
	unsigned int *maps = NULL;
	unsigned int i = 0, counts = 0;

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		maps = hashFindVal(mapsHash, currSeq->name);
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - kmerLength; i++)
			{
				temp = currSeq->dna[i+kmerLength];
				currSeq->dna[i+kmerLength] = '\0';
				currKmer = &(currSeq->dna[i]);
				counts = getKmerCount(currKmer, kmers);
				maps[i]=counts;
				currSeq->dna[i+kmerLength] = temp;
			}
		}
	}
}


void writeOutput(struct hash *noGapHash, struct hash *mapsHash, struct kmerData *kd, unsigned int kmerLength, char *outFilename)
{
	struct bed *bunk = NULL;
	unsigned int *maps = NULL;
	unsigned int i = 0;
	struct hashEl *hel = NULL;
	struct hashCookie cookie = hashFirst(mapsHash);
	FILE *fout = mustOpen(outFilename, "w");

	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		maps = hel->val;
		for(bunk = hashMustFindVal(noGapHash, hel->name); bunk != NULL; bunk=bunk->next)
		{
			fprintf(fout, "fixedStep chrom=%s start=%u step=1 kmerLength %u totalKmers %u\n", bunk->chrom, bunk->chromStart+1, kd->length, kd->totalCount/2);
			for(i=bunk->chromStart; i <= bunk->chromEnd - kmerLength; i++)
			{
				if(maps[i] == 0){errAbort("How is this zero? chrom=%s base=%u\n", bunk->chrom, i);}
				fprintf(fout, "%u\n", maps[i] - 1);
			}
		}
	}
	carefulClose(&fout);
}


void printKmerCount(struct hash *kmers)
{
	struct hashEl *hel = NULL;
	unsigned int *count = NULL;

	struct hashCookie cookie = hashFirst(kmers);
	while ((hel = hashNext(&cookie)) != NULL)
	{
		count = hel->val;
		uglyf("%s %u\n", hel->name, *count);
	}
}


/*---------------------------------------------------------------------------*/

void countKmers(unsigned int kmerLength, char *inFaFile, char *noGapBedFile, char *outFile)
{
	struct dnaSeq *seqList = NULL;
	struct hash *mapsHash = NULL, *noGapHash = NULL;
	struct kmerData *kmers = NULL;

	verbose(2, "reading in fa file\n");
	seqList = faReadAllDna(inFaFile);

	verbose(2, "reading in no gap file\n");
	noGapHash = bedLoadNInHash(noGapBedFile, 3);

	verbose(2, "finding the kmers that appear more than once\n");
	kmers = populateKmerHash(seqList, noGapHash, kmerLength);

	verbose(2, "counting how many times each kmer appears\n");
	recordKmers(seqList, noGapHash, kmerLength, kmers);

	verbose(2, "creating struct for base-wise data\n");
	mapsHash = createBasewiseHash(noGapHash);

	verbose(2, "calculating uniqueness of each position\n");
	calcUniqueness(seqList, noGapHash, kmerLength, kmers, mapsHash);

	verbose(2, "writing uniqueness of each position\n");
	writeOutput(noGapHash, mapsHash, kmers, kmerLength, outFile);

	verbose(2, "done\n");
}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
/* Process command line. */
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 5){usage();}

	unsigned int kmerLength = sqlUnsigned(argv[1]);

	countKmers(kmerLength,argv[2],argv[3],argv[4]);

	return(0);
}

