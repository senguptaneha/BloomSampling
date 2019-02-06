#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../bloom.h"
#include "bloomTreeDynamic.c"
extern int DEBUG;

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
	double samplingTime = 0, insertionTime, daTime = 0;
	int lineNum = 0;
	double ratio = 0.6;
	int method = 0;
	int numTimeslots = 150;
	double *meanArr = (double*)malloc(sizeof(double)*numTimeslots);
	double *varArr = (double*)malloc(sizeof(double)*numTimeslots);
	int *numObs = (int*)malloc(sizeof(int)*numTimeslots);
	int numConfigRounds = 1;
	double threshold = 10;
	while (ratio < 1){
		for (method = 0;method <1; method++){
			int configRound = 0;
			for (i=0;i<numTimeslots;i++){
				meanArr[i] = 0;
				varArr[i] = 0;
				numObs[i] = 0;
			}
			for (configRound = 0; configRound < numConfigRounds; configRound++){
				sampleAndPopulateLeaves(a,method,ratio);
				DEBUG = 1;
				printf("HELLO CHECKPOINT 1\n");
				fflush(stdout);
				FILE *fi = fopen(argv[3],"r");
				int lineNum = 0;
				while (fgets(inp,10000,fi)){
					//printf("C1, %s\n",inp); fflush(stdout);
					char *tok = strtok(inp,"\t");
					int hashId;
					long uId;
					//printf("C2\n"); fflush(stdout);
					if (tok){
						hashId = atoi(tok);
						tok = strtok(NULL,"\t");
						tok = strtok(NULL,"\t");
						if (tok)
							uId = atol(tok);
						else continue;
					}
					else continue;
					//printf("C3\n"); fflush(stdout);
					if (searchInTree(a,uId)){
			//			printf("C4, %d, %d\n", hashId, nSets);fflush(stdout);
						insert(uId,bfList[hashId]);
			//			printf("C5\n");fflush(stdout);
					}
					//printf("C10\n"); fflush(stdout);
					lineNum += 1;
				}
				fclose(fi);
				printf("HELLO CHECKPOINT 2\n");
				fflush(stdout);
				for (i=0;i<nSets;i++){
					//check if there is a zero intersection with any of the levels in the BloomSampleTree
					struct bloomNode *bN = emptyIntersectionTree(a, bfList[i], threshold);
					printf("Done %d\n",i); fflush(stdout);
					if (bN){
						printf("%d \t %d \t %ld \t %ld \t %d\n",i,bN->level,bN->start, bN->end, numDescendantLeaves(bN)); fflush(stdout);
					}
				}
				resetTree(a);
				for (i=0;i<nSets;i++) init(bfList[i]);
			}
		}
		ratio *= 2;
	}


}
