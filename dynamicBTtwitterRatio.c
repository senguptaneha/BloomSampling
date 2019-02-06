#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../bloom.h"
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
	double ratio = 0.01;
	int method = 0;
	int numTimeslots = 150;
	double *meanArr = (double*)malloc(sizeof(double)*numTimeslots);
	double *varArr = (double*)malloc(sizeof(double)*numTimeslots);
	int *numObs = (int*)malloc(sizeof(int)*numTimeslots);
	int numConfigRounds = 5;

	while (ratio < 1){
		for (method = 0;method <=1; method++){
			int configRound = 0;
			for (i=0;i<numTimeslots;i++){
				meanArr[i] = 0;
				varArr[i] = 0;
				numObs[i] = 0;
			}
			for (configRound = 0; configRound < numConfigRounds; configRound++){
				sampleAndPopulateLeaves(a,method,ratio);

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
					if (searchInTree(a,uId))
						insert(uId,bfList[hashId]);
					lineNum += 1;
					if (lineNum %1000000 == 0){
						samplingTime = 0;
						int sample,j;
						for (j=0;j<numSamplingRounds;j++){
							int h = rand()%nSets;
							clock_gettime(CLOCK_MONOTONIC, &start);
							int query_ones = num_ones(bfList[h]);
							sample = sampleFromTree(bfList[h],a,query_ones);
							clock_gettime(CLOCK_MONOTONIC, &finish);
							samplingTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
						}
						samplingTime /= numSamplingRounds;
						i = lineNum/1000000;
						numObs[i] += 1;
						if (numObs[i] == 1){
							meanArr[i] = samplingTime;
						}
						else{
							double prevMean = meanArr[i];
							meanArr[i] = meanArr[i] + (samplingTime - meanArr[i])/numObs[i];
							varArr[i] = varArr[i] + (samplingTime - prevMean)*(samplingTime - meanArr[i]);
						}
					}
				}
				fclose(fi);
				resetTree(a);
				for (i=0;i<nSets;i++) init(bfList[i]);
			}
			printf("OUT: **************method = %d, ratio = %lf***************\n",method,ratio);
			for (i=0;i<numTimeslots;i++){
				printf("OUT: %d\t%lf\t%lf\n",i,meanArr[i],sqrt(varArr[i]));
			}
		}
		ratio *= 2;
	}


}
