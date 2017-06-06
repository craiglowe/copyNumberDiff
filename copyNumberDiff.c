/*

copyNumberDiff.c
Written by Craig Lowe

*/

#include "common.h"
#include "options.h"
#include "memalloc.h"
#include "hash.h"
#include "sqlNum.h"
#include "linefile.h"
#include "dystring.h"
#include "bed.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "sam.h"
#include "bamFile.h"
#define MATHLIB_STANDALONE
#include "Rmath.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"roi", OPTION_STRING},
	{"maxReadDepth", OPTION_INT},
	{"logProbs", OPTION_BOOLEAN},
	{"debugStats", OPTION_BOOLEAN},
	{"chrom", OPTION_STRING},
	{"minQual", OPTION_INT},
	{"gcPseudo", OPTION_DOUBLE},
	{"noHet", OPTION_BOOLEAN},
	{"epsilon", OPTION_DOUBLE},
	{NULL, 0}
};

char *optChrom = NULL;
char *optRoi = NULL;
int optMinQual = 0;
boolean optLogProbs = FALSE;
boolean optDebugStats = FALSE;
double optGcPseudo = 0;
boolean optNoHet = FALSE;
double optEpsilon = 0.01;

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	"copyNumberDiff - Compute probabilities for HMM states corresponding to relative reductions, consistent, and amplifications.\n"
	"usage:\n"
	"   copyNumberDiff noGap.bed kmerMismaps.wig gcContent.wig group1_1.bam,...,group1_N.bam group2_1.bam,...,group2_N.bam outputEmissionProbs.wig\n"
	"options:\n"
	"   -roi             (NULL)     A bed file listing regions of interest that will be analyzed instead of the entire ungapped region\n"
	"   -minQual         (0)        Mapping quality must be at least this\n"
	"   -chrom           (NULL)     Restrict both chrom info and bams to only this chromosome\n"
	"   -logProbs        (false)    Do some of the calculations and report the emission probs as log(probs)\n"
	"   -debugStats      (false)    Give output for debugging different distributions\n"
	"   -gcPseudo        (0)        Can get zero probabilities if GC-content of a genomic region was never seen in a read\n"
	"   -noHet           (false)    Eliminate the states for heterozygous deletions and duplications\n"
	"   -epsilon         (0.01)     Adjustment to expected mismapping\n"
	);
}

/*---------------------------------------------------------------------------*/

struct slChrom
{
	struct slChrom *next;
	char *name;
	unsigned int length;
};


struct stateProbs
{
	double twoLessProb;
	double twoSameProb;
	double twoMoreProb;
};


struct slChrom *createChrom(char *name, unsigned int length)
{
	struct slChrom *ret;
	AllocVar(ret);
	ret->name = cloneString(name);
	ret->length = length;
	return(ret);
}


double addLog(double x, double y)
{
	if(x == -INFINITY){return(y);}
	if(y == -INFINITY){return(x);}
	if(x >= y){return(x + log(1 + exp(y-x)));}
	else{return(y + log(1 + exp(x-y)));}
}


double multiplyLog(double x, double y)
{
	if(x == -INFINITY || y == -INFINITY){return(-INFINITY);}
	else{return(x + y);}
}


double binomPdfApproxPoisson(unsigned int k, double p, double n)
{
	return(gsl_ran_poisson_pdf(k, n * p));
}


double negBinomPdfMuSize(unsigned int k, double mu, double size)
{
	double p = size / (size + mu);
	double n = mu * p / (1.0 - p);
	return(gsl_ran_negative_binomial_pdf(k, p, n));
}


struct hash *bedLoadNInHash(char *filename, int fields, char *restrictToChrom)
{
	struct bed *regions = NULL, *bed = NULL, *currList = NULL, *temp = NULL, *nextBed = NULL;
	struct hash *regionHash = newHash(6);

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
			if(restrictToChrom == NULL || sameString(bed->chrom, restrictToChrom))
			{
				hashAdd(regionHash, bed->chrom, currList);
			}
			else
			{
				bedFreeList(&currList);
			}
			currList = temp;
		}
	}
	return(regionHash);
}


struct slChrom *chromListFromNoGapHash(struct hash *noGapHash)
{
	struct hashCookie cookie = hashFirst(noGapHash);
	unsigned int max = 0;
	struct bed *headBed = NULL, *bunk = NULL;
	struct slChrom *headChrom = NULL, *currChrom = NULL;
	struct hashEl *hel = NULL;

	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		max = 0;
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			if(bunk->chromEnd > max){max = bunk->chromEnd;}
		}
		currChrom = createChrom(headBed->chrom, max);
		currChrom->next = headChrom;
		headChrom = currChrom;
	}
	return(headChrom);
}


struct hash *chromListToUnsignedHash(struct slChrom *chromList)
{
	struct slChrom *currChrom = NULL;
	struct hash *answer = newHash(6);

	for(currChrom=chromList; currChrom != NULL; currChrom=currChrom->next)
	{
		unsigned int *baseData = NULL;
		AllocArray(baseData, currChrom->length);
		hashAddUnique(answer, currChrom->name, baseData);
	}

	return(answer);
}


void setMismapDataFromWig(char *wigFilename, struct hash *baseDataHash, unsigned int *retKmerLength, unsigned int *retTotalKmers)
{
	struct lineFile *lf = NULL;
	char *row[8];
	unsigned int base = 0;
	char *chrom = NULL;
	unsigned int *baseData = NULL;

	lf = lineFileOpen(wigFilename, TRUE);
	while(lineFileChop(lf, row))
	{
		if(sameString(row[0], "fixedStep"))
		{
			if(differentString(row[3], "step=1")){errAbort("expecting step=1 when fixedStep\n");}
			if(differentString(row[4], "kmerLength")){errAbort("expecting kmerLength to be after step=1\n");}
			if(differentString(row[6], "totalKmers")){errAbort("expecting totalKmers to be after the kmer length\n");}
			chrom = row[1] + 6;
			base = sqlUnsigned(row[2] + 6) - 1;
			*retKmerLength = sqlUnsigned(row[5]);
			*retTotalKmers = sqlUnsigned(row[7]);
			baseData = (unsigned int *)hashFindVal(baseDataHash, chrom);
		}
		else if(baseData != NULL)
		{
			if(baseData == NULL){errAbort("wig file not in right format, don't know the chrom\n");}
			baseData[base] = sqlUnsigned(row[0]);
			base++;
		}
	}
	lineFileClose(&lf);
}


void setGcDataFromWig(char *wigFilename, struct hash *baseDataHash, double *gcBins, unsigned int kmerLength)
{
	struct lineFile *lf = NULL;
	char *row[4];
	unsigned int i = 0, base = 0, data = 0;
	char *chrom = NULL;
	unsigned int *baseData = NULL;
	double totalDataPoints = 0;
	
	lf = lineFileOpen(wigFilename, TRUE);
	while(lineFileChop(lf, row))
	{
		if(sameString(row[0], "fixedStep"))
		{
			if(differentString(row[3], "step=1")){errAbort("expecting step=1 when fixedStep\n");}
			chrom = row[1] + 6;
			base = sqlUnsigned(row[2] + 6) - 1;
			baseData = (unsigned int *)hashFindVal(baseDataHash, chrom);
		}
		else
		{
			data = sqlUnsigned(row[0]);
			gcBins[data]++;
			totalDataPoints++;
			if(baseData != NULL)
			{
				baseData[base] = data;
				base++;
			}
		}
	}
	lineFileClose(&lf);

	for(i=0; i<=kmerLength; i++)
	{
		gcBins[i] /= totalDataPoints;
	}
}


void probOfDepthBasicBinom(unsigned int depth, unsigned int copyNumber, double totalReads, double *retProb, double *retDensity)
{
	double prob = 0;

	prob = copyNumber / 2.0 / 446627861.0 * 0.91 + 0.09 * 1.0 / 446627861.0;
	*retProb = prob;
	*retDensity = binomPdfApproxPoisson(depth, prob, totalReads*36.0);
}


void probOfDepthBasicNegBinom(unsigned int depth, unsigned int copyNumber, double totalReads, double *retMean, double *retSize, double *retDensity)
{
	double prob = 0, mean = 0, var = 0;

	prob = copyNumber / 2.0 / 446627861.0 * 0.91 + 0.09 * 1.0 / 446627861.0;
	mean = prob * totalReads * 36.0;
	if(copyNumber == 0){var = 0.788975469945228;}
	else if(copyNumber == 1){var = 1.95722512845021;}
	else if(copyNumber == 2){var = 3.64840107796587;}
	else if(copyNumber == 3){var = 4.84091327080995;}
	else if(copyNumber == 4){var = 7.06653026718155;}
	else{errAbort("Error: don't know size for negative binomial");}
	*retMean = mean;
	*retSize = mean * mean / (var - mean);
	*retDensity = negBinomPdfMuSize(depth, mean, mean * mean / (var - mean));
}


void probOfDepthBasicDynamic(unsigned int depth, unsigned int copyNumber, double totalReads, unsigned int totalKmers, unsigned int kmerLength, unsigned int *kmerMismaps, struct bed *noGapBed, unsigned int basePos, double *retProb, double *retDensity)
{
	unsigned int i = 0, start = 0;
	double prob = 0, numer = 0, denom = 0, correctedMisMaps = 0;

	if(noGapBed->chromStart + kmerLength > basePos){start = noGapBed->chromStart;}
	else{start = 1 + basePos - kmerLength;}

	for(i=start; i<=basePos; i++)
	{
		correctedMisMaps = (double)kmerMismaps[i] + 0.01 * ((double)kmerMismaps[i] + 1.0);
		numer = (double)copyNumber + correctedMisMaps * 2.0;
		denom = (double)totalKmers * 2.0 * (correctedMisMaps + 1.0);
		prob += numer / denom;
	}
	*retProb = prob;
	*retDensity = binomPdfApproxPoisson(depth, prob, totalReads);
}


double probOfDepth(unsigned int depth, unsigned int copyNumber, double totalReads, unsigned int totalKmers, unsigned int kmerLength, unsigned int *kmerMismaps, struct bed *noGapBed, unsigned int basePos, unsigned int *gcContent, double *gcCorrection)
{
	unsigned int i = 0, start = 0;
	double prob = 0, numer = 0, denom = 0, correctedMisMaps = 0;

	/*if(optBasicStats)
	{
		double newCN = (double)copyNumber + 0.0001;
		double temp = binomPdfApproxPoisson(depth, newCN/2.0/446627861.0, totalReads*36.0);
		verbose(4, "prob is %e : depth=%u copyNumber=%f totalReads=%f\n", temp, depth, newCN, totalReads);
		return(log(temp));
	}*/

	if(noGapBed->chromStart + kmerLength > basePos){start = noGapBed->chromStart;}
	else{start = 1 + basePos - kmerLength;}

	for(i=start; i<=basePos; i++)
	{
		correctedMisMaps = (double)kmerMismaps[i] + optEpsilon * ((double)kmerMismaps[i] + 1.0);
		numer = (double)copyNumber + correctedMisMaps * 2.0;
		denom = (double)totalKmers * 2.0 * (correctedMisMaps + 1.0);
		prob += numer / denom * gcCorrection[gcContent[basePos]];
	}
	verbose(4, "depth=%u prob=%e totalReads=%f\n", depth, prob, totalReads);
	//return(log(binomPdfApproxPoisson(depth, prob, totalReads)));
	return(dbinom(depth, totalReads, prob, 1));
}


int maxUnsignedInt(unsigned int a, unsigned int b)
{
	if(a>=b){return(a);}
	else{return(b);}
}


struct hash *initDataStructure(struct slChrom *chromList, unsigned int numberOfSamples)
{
	unsigned int x = 0;
	struct hash *coverageHash = newHash(6);
	struct slChrom *currChrom = NULL;
	unsigned int **baseCoverage = NULL;

	for(currChrom = chromList; currChrom != NULL; currChrom = currChrom->next)
	{
		AllocArray(baseCoverage, numberOfSamples);
		for(x=0; x<numberOfSamples; x++)
		{
			AllocArray(baseCoverage[x], currChrom->length);
		}
		hashAddUnique(coverageHash, currChrom->name, baseCoverage);
	}
	return(coverageHash);
}


double addReadCounts(struct hash *coverageHash, char *filename, int sampleNumber, double gcPseudoCounts, double **readGcBins, double *genomicGcBins, unsigned int kmerLength)
{
	samfile_t *samFp = NULL;
	samFp = samopen(filename, "rb", 0);
	bam1_t* b = bam_init1();
	bam1_core_t *c = NULL;
	char *chrom = NULL;
	int chromId = -1;
	unsigned int x=0, chromStart=0, queryCount=0, refCount=0, at=0, gc=0, ns=0;
	unsigned int **coverage = NULL;
	char op = '?';
	int n = 0, i = 0;
	unsigned int *cigarPacked = NULL;
	//double totalBaseCoverage = 0;
	double used = 0, notUsed = 0, usedGc = 0;
	uint8_t *seq = NULL;
	uint32_t base = 0;
	char secondTryChrom[1024];

	while(samread(samFp, b) >= 0)
	{
		c = &b->core;
		//added a check that read is at least kmerLength in length.  Longer reads will be trimmed
		//verbose(4, "query length: %d  kmerLength: %u\n", c->l_qseq, kmerLength);
		if((c->qual >= optMinQual) && (c->l_qseq >= kmerLength) && (!(c->flag&(BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FUNMAP))))
		{
			//verbose(4, "inside\n");
			used++;

			if(c->tid != chromId)
			{
				chromId = c->tid;
				chrom = samFp->header->target_name[chromId];
				sprintf(secondTryChrom, "chr%s", chrom);
				verbose(3, "Looking for coverage data structure for %s or %s\n", chrom, secondTryChrom);
				if(optChrom == NULL || sameString(optChrom, chrom) || sameString(optChrom, secondTryChrom))
				{
					coverage = hashFindVal(coverageHash, chrom);
					if(coverage == NULL){coverage = hashFindVal(coverageHash, secondTryChrom);}
				}
				else{coverage = NULL;}
				if(coverage == NULL){verbose(3, "Did not find coverage data structure for %s or %s\n", chrom, secondTryChrom);}
				else{verbose(3, "Found coverage data structure for %s or %s\n", chrom, secondTryChrom);}
			}

			seq = bam1_seq(b);
			for(i = 0, at = 0, gc = 0, ns = 0;  i < kmerLength;  i++)
			{
				base = bam1_seqi(seq, i);
				if(base == 1 || base == 8){at++;}
				else if(base == 2 || base == 4){gc++;}
				else if(base == 15){ns++;}
				else{errAbort("Base: %u not recognized\n", base);}
			}
			if(ns == 0)
			{
				usedGc++;
				readGcBins[sampleNumber][gc]++;
			}
			
			if(coverage != NULL)
			{
				chromStart = c->pos;
				cigarPacked = bam1_cigar(b);
				for (i = 0, refCount = 0, queryCount=0;  i < c->n_cigar;  i++)
				{
					n = bamUnpackCigarElement(cigarPacked[i], &op);
					if(op == 'M' || op == '=' || op == 'X')
					{
						for(x = 0; x < n; x++)
						{
							if(queryCount < kmerLength)
							{
								//coverage[sampleNumber][chromStart+refCount+x]++;
								coverage[sampleNumber][chromStart+refCount]++;
								queryCount++;
								refCount++;
							}
						}
						//totalBaseCoverage += x;
						//refCount += x;
					}
					else if(op == 'D' || op == 'N')
					{
						refCount += n;	
					}
					else if(op == 'I')
					{
						queryCount +=n;
					}
					else if(op == 'P')
					{
					}
					else if(op == 'S' || op == 'H')
					{
					}
					else
					{
						errAbort("bamShowCigarEnglish: unrecognized CIGAR op %c -- update me", op);
					}
				}
				//used++;
			}
			/*else
			{
				cigarPacked = bam1_cigar(b);
				//cigarCount += c->n_cigar;
				for (i = 0, refCount = 0;  i < c->n_cigar;  i++)
				{
					n = bamUnpackCigarElement(cigarPacked[i], &op);
					if(op == 'M' || op == '=' || op == 'X')
					{
						totalBaseCoverage += n;
					}
				}
			}*/
		}
		else{notUsed++;}
	}
	samclose(samFp);
	bam_destroy1(b);
	verbose(3, "used: %e notUsed: %e\n", used, notUsed);
	for(i=0; i<=kmerLength; i++)
	{
		readGcBins[sampleNumber][i] = ((gcPseudoCounts + readGcBins[sampleNumber][i]) / (gcPseudoCounts * (kmerLength+1) + usedGc)) / genomicGcBins[i];
		if(readGcBins[sampleNumber][i] == 0)
		{
			errAbort("Error: calculated a probability of zero for a read mapping based on GC-content.  You should increase pseudocounts.");
		}
	}
	return(used);
}


void calcEmissionProbs(unsigned int numSamplesOne, unsigned int numSamplesTwo, double *totalReads, unsigned int basePos, unsigned int **depth, unsigned int maxCopy, struct bed *noGapBed, unsigned int totalKmerCount, unsigned int kmerLength, unsigned int *kmerMisMaps, double *currProbsOne, double *currProbsTwo, double **currProbs, unsigned int *gcContent, double **gcCorrection)
{
	unsigned int copyNumber = 0, groupOneCopyNumber = 0, groupTwoCopyNumber = 0, sampleNumber = 0, numSamples = 0;
	double lastProb = 0;

	numSamples = numSamplesOne + numSamplesTwo;
	for(copyNumber=0; copyNumber<=maxCopy; copyNumber+=1)
	{
		if(optNoHet == FALSE || copyNumber % 2 == 0)
		{
			currProbsOne[copyNumber] = 0;
			currProbsTwo[copyNumber] = 0;
			for(sampleNumber=0; sampleNumber < numSamplesOne; sampleNumber++)
			{
				lastProb =  currProbsOne[copyNumber];
				currProbsOne[copyNumber] = multiplyLog(currProbsOne[copyNumber], probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber])); //pm->prob[sampleNumber][depth[sampleNumber][basePos]][copyNumber]);
				if(currProbsOne[copyNumber] == -INFINITY)
				{
					errAbort("Error: Probability has gone to zero: %e %e %u %u %u\n", lastProb, probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber]), depth[sampleNumber][basePos], copyNumber, basePos);
				}
				if(isnan(currProbsOne[copyNumber]))
				{
					errAbort("Error: Probability is not defined: %e %e %u %u %u\n", lastProb, probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber]), depth[sampleNumber][basePos], copyNumber, basePos);
				}
			}
			for(sampleNumber=numSamplesOne; sampleNumber < numSamples; sampleNumber++)
			{
				lastProb = currProbsTwo[copyNumber];
				currProbsTwo[copyNumber] = multiplyLog(currProbsTwo[copyNumber], probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber]));
				if(currProbsTwo[copyNumber] == -INFINITY)
				{
					errAbort("Error: Probability has gone to zero %e %e %u %u %u\n", lastProb, probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber]), depth[sampleNumber][basePos], copyNumber, basePos);
				}
				if(isnan(currProbsTwo[copyNumber]))
				{
					errAbort("Error: Probability is not defined: %e %e %u %u %u\n", lastProb, probOfDepth(depth[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, gcContent, gcCorrection[sampleNumber]), depth[sampleNumber][basePos], copyNumber, basePos);
				}
			}
		}
	}
	
	for(groupOneCopyNumber=0; groupOneCopyNumber <= maxCopy; groupOneCopyNumber+=1)
	{
		for(groupTwoCopyNumber=0; groupTwoCopyNumber <= maxCopy; groupTwoCopyNumber+=1)
		{
			if(optNoHet == FALSE || (groupOneCopyNumber % 2 == 0 && groupTwoCopyNumber % 2 == 0))
			{
				currProbs[groupOneCopyNumber][groupTwoCopyNumber] = multiplyLog(currProbsOne[groupOneCopyNumber], currProbsTwo[groupTwoCopyNumber]);
			}
		}
	}
}


struct bed *findNoGapRegion(struct hash *noGapHash, struct bed *roi)
{
	struct bed *curr = NULL;
	for(curr=hashMustFindVal(noGapHash, roi->chrom); curr != NULL; curr=curr->next)
	{
		if(curr->chromStart <= roi->chromStart && curr->chromEnd >= roi->chromEnd)
		{
			return(curr);
		}
	}
	errAbort("Error: no noGap region containing this roi was located\n");
	return(NULL);
}


void printIntermediateData(unsigned int numSamplesOne, unsigned int numSamplesTwo, unsigned int maxCopy, struct hash *noGapHash, struct hash *roiHash, struct hash *coverageHash, double *totalReads, unsigned int totalKmerCount, unsigned int kmerLength, struct hash *kmerMisMapsHash, struct hash *gcContentHash, double **gcCorrection, char *outFilename)
{
	unsigned int basePos = 0;
	unsigned int **coverage = NULL;
	double probBinom = 0, densityBinom = 0, meanNegBinom = 0, sizeNegBinom = 0, densityNegBinom = 0, probDynamicBinom = 0, densityDynamicBinom = 0;
	unsigned int numSamples = 0, sampleNumber = 0, copyNumber = 0;
	unsigned int *kmerMisMaps = NULL, *gcData = NULL;

	verbose(3, " Allocating memory\n");

	numSamples = numSamplesOne + numSamplesTwo;

	FILE *fout = mustOpen(outFilename, "w");

	struct hashCookie cookie = hashFirst(roiHash);
	struct bed *headBed = NULL, *bunk = NULL, *noGapBed = NULL;
	struct hashEl *hel = NULL;
	double *logLikeOne = NULL, *logLikeTwo = NULL, *logLikeThree = NULL, *mseOne = NULL, *mseTwo = NULL, *mseThree = NULL;
	AllocArray(logLikeOne, maxCopy + 1);
	AllocArray(logLikeTwo, maxCopy + 1);
	AllocArray(logLikeThree, maxCopy + 1);
	AllocArray(mseOne, maxCopy + 1);
	AllocArray(mseTwo, maxCopy + 1);
	AllocArray(mseThree, maxCopy + 1);
	for(copyNumber=0;copyNumber<=maxCopy;copyNumber++)
	{
		logLikeOne[copyNumber] = 0;
		logLikeTwo[copyNumber] = 0;
		logLikeThree[copyNumber] = 0;
		mseOne[copyNumber] = 0;
		mseTwo[copyNumber] = 0;
		mseThree[copyNumber] = 0;
	}

	verbose(3, " Writing output...\n");
	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		coverage = hashMustFindVal(coverageHash, headBed->chrom);
		kmerMisMaps = hashMustFindVal(kmerMisMapsHash, headBed->chrom);
		gcData = hashMustFindVal(gcContentHash, headBed->chrom);
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			verbose(3, " Working on %s:%u-%u\n", bunk->chrom, bunk->chromStart, bunk->chromEnd);
			fprintf(fout, "fixedStep chrom=%s start=%u step=1\n", bunk->chrom, bunk->chromStart + 1);
			noGapBed = findNoGapRegion(noGapHash, bunk);
			for(copyNumber=0;copyNumber<=maxCopy;copyNumber++)
			{
				logLikeOne[copyNumber] = 0;
				logLikeTwo[copyNumber] = 0;
				logLikeThree[copyNumber] = 0;
				mseOne[copyNumber] = 0;
				mseTwo[copyNumber] = 0;
				mseThree[copyNumber] = 0;
			}
			for(basePos=bunk->chromStart; basePos < bunk->chromEnd; basePos++)
			{
				sampleNumber = 0;
				fprintf(fout, "%s\t%u\t", bunk->chrom, basePos);
				for(copyNumber = 0; copyNumber <= maxCopy; copyNumber+=1)
				{
					probOfDepthBasicBinom(coverage[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], &probBinom, &densityBinom);
					probOfDepthBasicNegBinom(coverage[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], &meanNegBinom, &sizeNegBinom, &densityNegBinom);
					probOfDepthBasicDynamic(coverage[sampleNumber][basePos], copyNumber, totalReads[sampleNumber], totalKmerCount, kmerLength, kmerMisMaps, noGapBed, basePos, &probDynamicBinom, &densityDynamicBinom);
					fprintf(fout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t", probBinom, densityBinom, meanNegBinom, sizeNegBinom, densityNegBinom, probDynamicBinom, densityDynamicBinom);
					logLikeOne[copyNumber] = multiplyLog(logLikeOne[copyNumber], log(densityBinom));
					logLikeTwo[copyNumber] = multiplyLog(logLikeTwo[copyNumber], log(densityNegBinom));
					logLikeThree[copyNumber] = multiplyLog(logLikeThree[copyNumber], log(densityDynamicBinom));
					mseOne[copyNumber] += (1.0 - densityBinom) * (1.0 - densityBinom);
					mseTwo[copyNumber] += (1.0 - densityNegBinom) * (1.0 - densityNegBinom);
					mseThree[copyNumber] += (1.0 - densityDynamicBinom) * (1.0 - densityDynamicBinom);
				}
				fprintf(fout, "%u\n", coverage[sampleNumber][basePos]);
			}
			fprintf(fout, "summaryBinom\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mseOne[0], logLikeOne[0], mseOne[1], logLikeOne[1], mseOne[2], logLikeOne[2], mseOne[3], logLikeOne[3], mseOne[4], logLikeOne[4]);
			fprintf(fout, "summaryNegBinom\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mseTwo[0], logLikeTwo[0], mseTwo[1], logLikeTwo[1], mseTwo[2], logLikeTwo[2], mseTwo[3], logLikeTwo[3], mseTwo[4], logLikeTwo[4]);
			fprintf(fout, "summaryDynamicBinom\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", mseThree[0], logLikeThree[0], mseThree[1], logLikeThree[1], mseThree[2], logLikeThree[2], mseThree[3], logLikeThree[3], mseThree[4], logLikeThree[4]);
		}
	}
	carefulClose(&fout);
}


void printProbabilities(unsigned int numSamplesOne, unsigned int numSamplesTwo, unsigned int maxCopy, struct hash *noGapHash, struct hash *coverageHash, double *totalReads, unsigned int totalKmerCount, unsigned int kmerLength, struct hash *kmerMisMapsHash, struct hash *gcContentHash, double **gcCorrection, char *outFilename)
{
	unsigned int basePos = 0;
	unsigned int **coverage = NULL;
	double **currProbs = NULL;
	double *currProbsOne = NULL, *currProbsTwo = NULL;
	unsigned int numSamples = 0, i = 0, j = 0;
	unsigned int *kmerMisMaps = NULL, *gcData = NULL;

	verbose(3, " Allocating memory\n");
	AllocArray(currProbsOne, maxCopy+1);
	AllocArray(currProbsTwo, maxCopy+1);
	AllocArray(currProbs, maxCopy+1);
	for(i=0; i<=maxCopy; i++)
	{
		AllocArray(currProbs[i], maxCopy+1);
	}

	numSamples = numSamplesOne + numSamplesTwo;

	FILE *fout = mustOpen(outFilename, "w");

	struct hashCookie cookie = hashFirst(noGapHash);
	struct bed *headBed = NULL, *bunk = NULL;
	struct hashEl *hel = NULL;

	verbose(3, " Writing output...\n");
	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		coverage = hashMustFindVal(coverageHash, headBed->chrom);
		kmerMisMaps = hashMustFindVal(kmerMisMapsHash, headBed->chrom);
		gcData = hashMustFindVal(gcContentHash, headBed->chrom);
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			verbose(3, " Working on %s:%u-%u\n", bunk->chrom, bunk->chromStart, bunk->chromEnd);
			fprintf(fout, "fixedStep chrom=%s start=%u step=1\n", bunk->chrom, bunk->chromStart + 1);
			for(basePos=bunk->chromStart; basePos < bunk->chromEnd; basePos++)
			{
				calcEmissionProbs(numSamplesOne, numSamplesTwo, totalReads, basePos, coverage, maxCopy, bunk, totalKmerCount, kmerLength, kmerMisMaps, currProbsOne, currProbsTwo, currProbs, gcData, gcCorrection);
				for(i=0; i<=maxCopy; i+=1)
				{
					for(j=0; j<=maxCopy; j+=1)
					{
						if(optNoHet == FALSE || (i % 2 == 0 && j % 2 == 0))
						{
							if(optLogProbs){fprintf(fout, "%e\t", currProbs[i][j]);}
							else{fprintf(fout, "%e\t", exp(currProbs[i][j]));}
						}
					}
				}
				fprintf(fout, "\n");
			}
		}
	}
	carefulClose(&fout);
}


/*---------------------------------------------------------------------------*/


void copyNumberDiff(char *noGapFilename, char *misMapsWigFilename, char *gcContentWigFilename, char *inputBamFilenamesOne, char *inputBamFilenamesTwo, char *outFilename)
{
	struct hash *coverageHash = NULL, *roiHash = NULL, *noGapHash = NULL, *kmerMisMapsHash = NULL, *gcContentHash = NULL;
	struct slName *bamFilenamesListOne = NULL, *bamFilenamesListTwo = NULL, *currFilename = NULL;
	unsigned int numFilesOne=0, numFilesTwo=0, sampleNumber=0, numSamples=0;
	double *totalReads = NULL;
	unsigned int kmerLength = 0, totalKmerCount = 0;
	struct slChrom *chromList = NULL;
	unsigned int maxCopyNumber = 4;
	double *genomicGcBins = NULL;
	double **gcCorrection = NULL;

	verbose(2, "Reading in no gap file\n");
	noGapHash = bedLoadNInHash(noGapFilename, 3, optChrom);

	if(optRoi != NULL)
	{
		verbose(2, "Reading in ROI file\n");
		roiHash = bedLoadNInHash(optRoi, 3, optChrom);
	}

	verbose(2, "Creating a list of chroms from the no gap file\n");
	chromList = chromListFromNoGapHash(noGapHash);

	verbose(2, "Creating data structure for mismap data\n");
	kmerMisMapsHash = chromListToUnsignedHash(chromList);

	verbose(2, "Reading in mismap wig file\n");
	setMismapDataFromWig(misMapsWigFilename, kmerMisMapsHash, &kmerLength, &totalKmerCount);

	verbose(2, "Creating data structure for gcContent\n");
	gcContentHash = chromListToUnsignedHash(chromList);
	AllocArray(genomicGcBins, kmerLength+1);

	verbose(2, "Reading in gcContent wig file\n");
	setGcDataFromWig(gcContentWigFilename, gcContentHash, genomicGcBins, kmerLength);

	verbose(2, "Parsing bam filenames\n");
	bamFilenamesListOne = slNameListFromComma(inputBamFilenamesOne);
	bamFilenamesListTwo = slNameListFromComma(inputBamFilenamesTwo);
	numFilesOne = slCount(bamFilenamesListOne);
	numFilesTwo = slCount(bamFilenamesListTwo);
	numSamples = numFilesOne + numFilesTwo;
	AllocArray(totalReads, numSamples);
	AllocArray(gcCorrection, numSamples);

	verbose(2, "Setting up mapping data structure\n");
	coverageHash = initDataStructure(chromList, numSamples);

	verbose(2, "Reading bam files\n");
	sampleNumber = 0;
	for(currFilename=bamFilenamesListOne; currFilename!=NULL; currFilename=currFilename->next, sampleNumber++)
	{
		AllocArray(gcCorrection[sampleNumber], kmerLength+1);
		verbose(2, " Reading %s\n", currFilename->name);
		totalReads[sampleNumber] = addReadCounts(coverageHash, currFilename->name, sampleNumber, optGcPseudo, gcCorrection, genomicGcBins, kmerLength);
		verbose(3, "  %s has %f reads\n", currFilename->name, totalReads[sampleNumber]);
	}
	for(currFilename=bamFilenamesListTwo; currFilename!=NULL; currFilename=currFilename->next, sampleNumber++)
	{
		AllocArray(gcCorrection[sampleNumber], kmerLength+1);
		verbose(2, " Reading %s\n", currFilename->name);
		totalReads[sampleNumber] = addReadCounts(coverageHash, currFilename->name, sampleNumber, optGcPseudo, gcCorrection, genomicGcBins, kmerLength);
		verbose(3, "  %s has %f reads\n", currFilename->name, totalReads[sampleNumber]);
	}

	verbose(2, "Printing emission probabilities\n");
	if(optDebugStats && optRoi)
	{
		printIntermediateData(numFilesOne, numFilesTwo, maxCopyNumber, noGapHash, roiHash, coverageHash, totalReads, totalKmerCount, kmerLength, kmerMisMapsHash, gcContentHash, gcCorrection, outFilename);
	}
	else
	{
		if(optRoi == NULL){printProbabilities(numFilesOne, numFilesTwo, maxCopyNumber, noGapHash, coverageHash, totalReads, totalKmerCount, kmerLength, kmerMisMapsHash, gcContentHash, gcCorrection, outFilename);}
		else{printProbabilities(numFilesOne, numFilesTwo, maxCopyNumber, roiHash, coverageHash, totalReads, totalKmerCount, kmerLength, kmerMisMapsHash, gcContentHash, gcCorrection, outFilename);}
	}
}


/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 7){usage();}

	optChrom = optionVal("chrom", optChrom);
	optRoi = optionVal("roi", optRoi);
	optMinQual = optionInt("minQual", optMinQual);
	optLogProbs = optionExists("logProbs");
	optDebugStats = optionExists("debugStats");
	optGcPseudo = optionDouble("gcPseudo", optGcPseudo);
	optNoHet = optionExists("noHet");
	optEpsilon = optionDouble("epsilon", optEpsilon);

	copyNumberDiff(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	return(0);
}

