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
        FILE *fi = fopen(argv[3],"r");
	int lineNum = 0;
	while (fgets(inp,10000,fi) && (lineNum < 1000000)){
		/*printf("Checkpoint 1, %s\n", inp);
		fflush(stdout);*/

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
		insert(uId,bfList[hashId]);
		clock_gettime(CLOCK_MONOTONIC,&start);
		insertInTree(uId,a);
		clock_gettime(CLOCK_MONOTONIC, &finish);
		insertionTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		numInsertions++;

		int r = rand() % 1000;
		if (r==0){
			samplingTime = 0;
			int sample,j;
			for (j=0;j<numSamplingRounds;j++){
				int h = rand()%nSets;
				//int h = hashId;
				clock_gettime(CLOCK_MONOTONIC, &start);
				int query_ones = num_ones(bfList[h]);
				sample = sampleFromTree(bfList[h],a,query_ones);
				clock_gettime(CLOCK_MONOTONIC, &finish);
				samplingTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
			}

                        daTime = 0;
                        int h = rand()%nSets;
			clock_gettime(CLOCK_MONOTONIC, &start);
                        sample = daSample(nVertices,bfList[h]);
                        clock_gettime(CLOCK_MONOTONIC, &finish);
                        daTime = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);

		//	samplingTime /= numSamplingRounds;
			double memUsed = (1.0*malloc_count_current())/(1024*1024);	//memory used in MBs
			insertionTime = insertionTime/numInsertions;
			printf("%d\t%lf\t%lf\t%.2lf\t%lf\t%d\t%d\n", lineNum, samplingTime/numSamplingRounds, daTime, memUsed, insertionTime, a->numNodesPopulated, a->numTotalNodes);
			numInsertions = 0;
			insertionTime = 0;
		}
		/*printf("Memory used = %d\n",malloc_count_current());*/
		lineNum += 1;
	}
	fclose(fi);
}
