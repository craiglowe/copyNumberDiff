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
	unsigned int totalCount;
	struct hash *counts;
};


struct bloomFilter
{
	uint64_t *table;
	uint64_t tabLen;
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
	answer->tabLen = 5e8;
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

// The following code I read on wikipedia to calculate
// the Hamming Weight / Pop Count of an int64
const uint64_t m1  = 0x5555555555555555;
const uint64_t m2  = 0x3333333333333333;
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
int popCount(uint64_t x) {
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    x += x >> 8;
    x += x >> 16;
    x += x >> 32;
    return(x & 0x7f);
}


double bloomPercentUsed(struct bloomFilter *bloom)
{
	unsigned int idx = 0;
	double count = 0;
	for(idx=0; idx<bloom->tabLen; idx++)
	{
		count += popCount(bloom->table[idx]);
	}
	return(count / (64.0 * (double)bloom->tabLen));
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


void packDna(char *src, char *dest, unsigned int len)
{
	unsigned int packedLen = (len + 2) / 3;
	dest[packedLen] = '\0';
	char twoBit = 0;
	unsigned x=0, y=0;

	for(x=0; x*3<len; x++)
	{
		dest[x] = 32;
		for(y=0; y<3 && x*3+y<len; y++)
		{
			if(src[x*3+y] == 'A' || src[x*3+y] == 'a'){twoBit=0;}
                	else if(src[x*3+y] == 'C' || src[x*3+y] == 'c'){twoBit=1;}
                	else if(src[x*3+y] == 'G' || src[x*3+y] == 'g'){twoBit=2;}
                	else if(src[x*3+y] == 'T' || src[x*3+y] == 't'){twoBit=3;}
                	else{errAbort("Error: unexpected character (%c)\n", src[x*3+y]);}
			dest[x] += twoBit << (4-2*y);
		}
	}
	verbose(5, "unpacked:%s, packed:%s\n", src, dest);
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
		/*
 		unsigned int *counter = NULL;
		counter = hel->val;
		(*counter)++;
		*/
		hashIncInt(answer->counts, currKmer);
	}
}


int recordIndividualKmerFirstPass(char *currKmer, struct bloomFilter *bloom, struct kmerData *answer)
{
	struct hashEl *hel = NULL;

	if(checkBloom(currKmer, bloom))
	{
		hel = hashLookup(answer->counts, currKmer);
		if(hel == NULL)
		{
			/*
			uint32_t *temp = NULL;
			AllocVar(temp);
			*temp = 0;
			hashAdd(answer->counts, currKmer, temp);
			*/
			hashAddInt(answer->counts, currKmer, 0);
			return(1);
		}
	}
	else
	{
		putInBloom(currKmer, bloom);
	}
	return(0);
}


unsigned int getKmerCount(char *currKmer, struct kmerData *kmers)
{
	/*
	struct hashEl *hel = hashLookup(kmers->counts, currKmer);
	if(hel == NULL){return(1);}
	else
	{
		unsigned int *count = hel->val;
		return(*count);
	}
	*/
	return(hashIntValDefault(kmers->counts, currKmer, 1));
}


struct kmerData *populateKmerHash(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	int i = 0;
	char *currKmer = NULL;
	struct bloomFilter *bloom = createBloomFilter();
	char revCompKmer[kmerLength+1];
	unsigned int packedLen = (kmerLength + 2) / 3;
	char packedKmer[packedLen+1];
	char temp;
	struct kmerData *answer = createKmerData(kmerLength);
	unsigned int itemsInHash = 0;

	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		verbose(3, " currently on %s\n", currSeq->name);
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - kmerLength; i++)
			{
				/*temp = currSeq->dna[i+kmerLength];
				currSeq->dna[i+kmerLength] = '\0';
				currKmer = &(currSeq->dna[i]);
				recordIndividualKmerFirstPass(currKmer, bloom, answer);
				revComp(currKmer, revCompKmer, kmerLength);
				recordIndividualKmerFirstPass(revCompKmer, bloom, answer);
				currSeq->dna[i+kmerLength] = temp;*/

				temp = currSeq->dna[i+kmerLength];
                                currSeq->dna[i+kmerLength] = '\0';
                                currKmer = &(currSeq->dna[i]);
				revComp(currKmer, revCompKmer, kmerLength);
				if(strcmp(currKmer, revCompKmer) <= 0)
				{
					packDna(currKmer, packedKmer, kmerLength);
	                                itemsInHash += recordIndividualKmerFirstPass(packedKmer, bloom, answer);
				}
				else
				{
					packDna(revCompKmer, packedKmer, kmerLength);
                                	itemsInHash += recordIndividualKmerFirstPass(packedKmer, bloom, answer);
				}
                                currSeq->dna[i+kmerLength] = temp;
			}
			//verbose(4, "  itemsInHash:%u percentUsed:%f\n", itemsInHash, bloomPercentUsed(bloom));
		}
		verbose(3, "  itemsInHash:%u percentUsed:%f\n", itemsInHash, bloomPercentUsed(bloom));
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
	unsigned int packedLen = (kmerLength + 2) / 3;
	char packedKmer[packedLen+1];
	char temp;
	
	for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
	{
		verbose(3, " currently on %s\n", currSeq->name);
		for(currRegion = hashFindVal(noGapHash, currSeq->name); currRegion != NULL; currRegion = currRegion->next)
		{
			for(i = currRegion->chromStart; i <= currRegion->chromEnd - kmerLength; i++)
			{
				/*temp = currSeq->dna[i+kmerLength];
				currSeq->dna[i+kmerLength] = '\0';
				currKmer = &(currSeq->dna[i]);
				recordIndividualKmerSecondPass(currKmer, answer);
				revComp(currKmer, revCompKmer, kmerLength);
				recordIndividualKmerSecondPass(revCompKmer, answer);
				currSeq->dna[i+kmerLength] = temp;*/

				temp = currSeq->dna[i+kmerLength];
                                currSeq->dna[i+kmerLength] = '\0';
                                currKmer = &(currSeq->dna[i]);
				revComp(currKmer, revCompKmer, kmerLength);
				if(strcmp(currKmer, revCompKmer) <= 0)
				{
					packDna(currKmer, packedKmer, kmerLength);
                                	recordIndividualKmerSecondPass(packedKmer, answer);
				}
				else
				{
					packDna(revCompKmer, packedKmer, kmerLength);
                                	recordIndividualKmerSecondPass(packedKmer, answer);
				}
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
	char revCompKmer[kmerLength+1];
	unsigned int packedLen = (kmerLength + 2) / 3;
	char packedKmer[packedLen+1];

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
				revComp(currKmer, revCompKmer, kmerLength);
				if(strcmp(currKmer, revCompKmer) <= 0)
				{
					packDna(currKmer, packedKmer, kmerLength);
					counts = getKmerCount(packedKmer, kmers);
				}
				else
				{
					packDna(revCompKmer, packedKmer, kmerLength);
					counts = getKmerCount(packedKmer, kmers);
				}
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
			fprintf(fout, "fixedStep chrom=%s start=%u step=1 kmerLength %u totalKmers %u\n", bunk->chrom, bunk->chromStart+1, kd->length, kd->totalCount);
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
		//uglyf("%s %u\n", hel->name, *count);
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

