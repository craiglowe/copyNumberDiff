/*

countKmersWithTree.c
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
	"countKmersWithTree - count the number of times each k-mer appears\n"
	"usage:\n"
	"   countKmersWithTree kmerLength in.fa noGap.bed kmers.wig\n"
	"options:\n"
	"notes:\n"
	);
}

/*---------------------------------------------------------------------------*/


struct seqTreeNode
{
	unsigned int count;           // keeps track of how many sequences end up at this leaf node
	struct seqTreeNode **children; // array of {0,1,2,3} for {A,C,G,T}
};


struct seqTreeNode *createSeqTreeNode()
{
	struct seqTreeNode *answer = NULL;
	AllocVar(answer);
	answer->count = 0;
	answer->children = NULL;
	return(answer);
}


void insertInTree(char *seq, struct seqTreeNode *node)
{
	unsigned int idx = 0;

	if(node == NULL){errAbort("Error: inserting tree into null node\n");}

	if(seq[0] == '\0'){node->count++;}
	else
	{
		if(seq[0] == 'A' || seq[0] == 'a'){idx=0;}
		else if(seq[0] == 'C' || seq[0] == 'c'){idx=1;}
		else if(seq[0] == 'G' || seq[0] == 'g'){idx=2;}
		else if(seq[0] == 'T' || seq[0] == 't'){idx=3;}
		else{errAbort("Error: unexpected character %c\n", seq[0]);}
		
		if(node->children == NULL){AllocArray(node->children, 4);}
		if(node->children[idx] == NULL){node->children[idx] = createSeqTreeNode();}
		insertInTree(&(seq[1]), node->children[idx]);
	}
}


unsigned int getKmerCount(char *seq, struct seqTreeNode *node)
{
	unsigned int idx = 0;

	if(node == NULL){errAbort("Error: given a null node when getting counts from tree\n");}

	if(seq[0] == '\0'){return(node->count);}
	else
	{
		if(seq[0] == 'A' || seq[0] == 'a'){idx=0;}
		else if(seq[0] == 'C' || seq[0] == 'c'){idx=1;}
		else if(seq[0] == 'G' || seq[0] == 'g'){idx=2;}
		else if(seq[0] == 'T' || seq[0] == 't'){idx=3;}
		else{errAbort("Error: unexpected character %c\n", seq[0]);}

		if(node->children == NULL || node->children[idx] == NULL){errAbort("Error: found a null node when getting counts from tree\n");}
		return(getKmerCount(&(seq[1]), node->children[idx]));
	}
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


struct seqTreeNode *populateTree(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	int i = 0;
	char *currKmer = NULL;
	char revCompKmer[kmerLength+1];
	char temp;
	struct seqTreeNode *root = NULL;
	root = createSeqTreeNode();

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
				revComp(currKmer, revCompKmer, kmerLength);
				if(strcmp(currKmer, revCompKmer) <= 0)
				{
	                                insertInTree(currKmer, root);
				}
				else
				{
                                	insertInTree(revCompKmer, root);
				}
                                currSeq->dna[i+kmerLength] = temp;
			}
		}
	}
	return(root);
}


void calcUniqueness(struct dnaSeq *seqList, struct hash *noGapHash, unsigned int kmerLength, struct seqTreeNode *root, struct hash *mapsHash)
{
	struct dnaSeq *currSeq = NULL;
	struct bed *currRegion = NULL;
	char *currKmer = NULL;
	char temp;
	unsigned int *maps = NULL;
	unsigned int i = 0, counts = 0;
	char revCompKmer[kmerLength+1];

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
					counts = getKmerCount(currKmer, root);
				}
				else
				{
					counts = getKmerCount(revCompKmer, root);
				}
				maps[i]=counts;
				currSeq->dna[i+kmerLength] = temp;
			}
		}
	}
}


void writeOutput(struct hash *noGapHash, struct hash *mapsHash, unsigned int kmerLength, char *outFilename)
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
			fprintf(fout, "fixedStep chrom=%s start=%u step=1 kmerLength %u totalKmers %u\n", bunk->chrom, bunk->chromStart+1, kmerLength, 0);
			for(i=bunk->chromStart; i <= bunk->chromEnd - kmerLength; i++)
			{
				if(maps[i] == 0){errAbort("How is this zero? chrom=%s base=%u\n", bunk->chrom, i);}
				fprintf(fout, "%u\n", maps[i] - 1);
			}
		}
	}
	carefulClose(&fout);
}


/*---------------------------------------------------------------------------*/

void countKmers(unsigned int kmerLength, char *inFaFile, char *noGapBedFile, char *outFile)
{
	struct dnaSeq *seqList = NULL;
	struct hash *mapsHash = NULL, *noGapHash = NULL;
	struct seqTreeNode *root = NULL;

	verbose(2, "reading in fa file\n");
	seqList = faReadAllDna(inFaFile);

	verbose(2, "reading in no gap file\n");
	noGapHash = bedLoadNInHash(noGapBedFile, 3);

	verbose(2, "finding the kmers that appear more than once\n");
	root = populateTree(seqList, noGapHash, kmerLength);

	verbose(2, "creating struct for base-wise data\n");
	mapsHash = createBasewiseHash(noGapHash);

	verbose(2, "calculating uniqueness of each position\n");
	calcUniqueness(seqList, noGapHash, kmerLength, root, mapsHash);

	verbose(2, "writing uniqueness of each position\n");
	writeOutput(noGapHash, mapsHash, kmerLength, outFile);

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

