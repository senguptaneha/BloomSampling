#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "bloom.h"
#include "bloomTreeDynamic.c"

long daSample(long M, struct bloom *b){
	long i;
	long sample = -1;
	int numCandidates = 0;
	for (i=0;i<=M;i++){
		if (is_in(i,b)){
			numCandidates += 1;
			if (numCandidates > 0){
				int z = rand()%numCandidates;
				if (z==0) sample = i;
			}
		}
	}
	return sample;
}

struct inputSet{
	long *set;
	int size;
	int nextLoc;
};

struct inputSet* createSet(int initialSize){
	struct inputSet* iS = (struct inputSet*)malloc(sizeof(struct inputSet));
	iS->set = (long*)malloc(sizeof(long)*initialSize);
	int i;
	iS->size = initialSize;
	for (i=0;i< iS->size;i++) iS->set[i] = -1;
	iS->nextLoc = 0;
	return iS;
}

void reAllocate(struct inputSet* iS, int extraSize){
	long *newSet = (long*)malloc(sizeof(long)*(iS->size + extraSize));
	int i;
	for (i=0; i< iS->nextLoc; i++) newSet[i] = iS->set[i];
	for (i=iS->nextLoc; i< (iS->size + extraSize); i++) newSet[i] = -1;

	iS->size += extraSize;
	free(iS->set);
	iS->set = newSet;
}

void insertElemInSet(struct inputSet *iS, long elem){
	if (iS->nextLoc == iS->size) reAllocate(iS, 100);
	iS->set[iS->nextLoc] = elem;
	iS->nextLoc += 1;
}

int isInSet(struct inputSet *iS, long elem){
	int i;
	if (elem == -1) return 1;
	for (i=0; i< iS->nextLoc; i++)
		if (iS->set[i] == elem)
			return 1;
	return 0;
}

void resetSet(struct inputSet *iS){
	int i;
	for (i=0; i< iS->nextLoc; i++) iS->set[i] = -1;
	iS->nextLoc = 0;
}

int main(int argc, char* argv[]){
	if (argc<8){
		printf("Usage: ./dynamicBTtwitterTest M nSets inputFile k m T l\n");
		return 0;
	}
	struct timespec start,finish,mid;
	double elapsed=0.0;
	srand(time(NULL));
	setSeeds();	//these seeds are for murmur hash functions
	long nVertices = atol(argv[1]);
	int nSets = atoi(argv[2]); //number of hashTags
	K = atoi(argv[4]);
	VECT_SIZE = atoi(argv[5]);
	T = atoi(argv[6]);
	levelThreshold = atoi(argv[7]);

	seiveInitial();
	struct bloomTree *a = getBloomTree(0,nVertices);
	struct bloom ** bfList = (struct bloom**)malloc(nSets *sizeof(struct bloom*));
	int i;

	for (i=0;i<nSets;i++){
		bfList[i] = (struct bloom*)malloc(sizeof(struct bloom));
		bfList[i]->bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
		init(bfList[i]);
	}

	char inp[10000];
	int numSamplingRounds = 100, numInsertions = 0;
	double samplingTimeAvg, samplingTimeVar, insertionTime, daTime = 0, samplingAcc, numNodes;
	int lineNum = 0;
	double ratio = 0.2;
	int method = 0;
	int numConfigRounds = 3;
	struct inputSet **inputSets = (struct inputSet**)malloc(sizeof(struct inputSet*)*nSets);
	for (i=0; i<nSets; i++) inputSets[i] = createSet(1000);
	sampleAndPopulateLeaves(a,0,1);
	resetTreeNoFree(a);
	printf("Populated BST\n");
	while (ratio < 1){
		for (method = 1;method <=1; method++){
			int configRound = 0;
			samplingTimeAvg = 0;
			samplingTimeVar = 0;
			samplingAcc = 0;
			numNodes = 0;
			int nObs = 0;
			for (configRound = 0; configRound < numConfigRounds; configRound++){
				for (i=0;i<nSets; i++) resetSet(inputSets[i]);
				printf("Sampling for ratio = %lf\n",ratio);
				sampleLeavesOnly(a,method,ratio);
				printf("Done sampling for ratio = %lf\n",ratio);

				FILE *fi = fopen(argv[3],"r");
				int lineNum = 0;
				while (fgets(inp,10000,fi)){
					char *tok = strtok(inp,"\t");
					int hashId;
					long uId;
					if (tok){
						hashId = atoi(tok);
						tok = strtok(NULL,"\t");
						tok = strtok(NULL,"\t");
						if (tok)
							uId = atol(tok);
						else continue;
					}
					else continue;
					if (searchInTree(a,uId)){
						insert(uId,bfList[hashId]);
						insertElemInSet(inputSets[hashId],uId);
					}
					lineNum += 1;
				}
				int j;
				long sample;
				for (j=0;j<numSamplingRounds;j++){
					int h = rand()%nSets;
					int query_ones = num_ones(bfList[h]);
					clock_gettime(CLOCK_MONOTONIC, &start);
					sample = sampleFromTree(bfList[h],a,query_ones);
					clock_gettime(CLOCK_MONOTONIC, &finish);
					double sampleTime = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
					int isCorrect = isInSet(inputSets[h],sample);
					samplingAcc += isCorrect;
					nObs += 1;
					if (nObs == 1){
						samplingTimeAvg = sampleTime;
						samplingTimeVar = 0;
					}
					else{
						double prevAvg = samplingTimeAvg;
						samplingTimeAvg = samplingTimeAvg + (sampleTime - samplingTimeAvg)/nObs;
						samplingTimeVar = samplingTimeVar + (sampleTime - prevAvg)*(sampleTime - samplingTimeAvg);
					}
				}
				fclose(fi);

				resetTreeNoFree(a);
				for (i=0;i<nSets;i++) init(bfList[i]);
			}
			numNodes /= (numConfigRounds);
//			double memUsed = (numNodes * VECT_SIZE)/(8*1024*1024);	//memory used in MBs
//			printf("OUT: %d\t%lf\t%lf\n",method,ratio,memUsed);
			double d = (numConfigRounds*numSamplingRounds);
			printf("OUT: %d\t%lf\t%lf\t%lf\t%lf\n", method, ratio, samplingTimeAvg, sqrt(samplingTimeVar/(d-1)), (samplingAcc/d));
		}
		if (ratio < 0.1) ratio += 0.02;
		else ratio += 0.1;
	}


}

