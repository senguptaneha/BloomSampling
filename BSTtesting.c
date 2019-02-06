#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "bloomTree.c"
#include "../bloom.h"
extern int K;
extern int VECT_SIZE;
extern int T;
extern int levelThreshold;

int dictionarySample(struct bloom *a, int nVertices){
	int i = 0, currSample = -1, numCandidates = 0;
	for (i=0;i<nVertices;i++){
		if (is_in(i,a)){
			numCandidates++;
			int z = rand()%numCandidates;
			if (z==0) currSample = i;
		}
	}
	return currSample;
}


void updateTrueValue(struct bloomNode *a, int value){
	if (a==NULL) return;
	if ((a->start<=value)&&(a->end>=value)) a->trueValue++;
	updateTrueValue(a->lchild,value);
	updateTrueValue(a->rchild,value);
}

void updateFlags(struct bloomNode *a, struct bloom *q){
	if (a==NULL) return;
	a->flag = intersectNode(a,q);
	updateFlags(a->lchild,q);
	updateFlags(a->rchild,q);
}
int getRequiredm(int M, int n, int k, double a){
    double z = (n/a - n)/(M-n);
    z = 1 - pow(z,1.0/k);
    z = 1 - pow(z,1.0/(n*k));
    z = 1.0/z;
    return (int)(floor(z));
}

int bestLevel(int m, int n, int M){
    seiveInitial();

        struct bloom *k = (struct bloom *)malloc(sizeof(struct bloom));
        k->bloom_vector = (int*)malloc(sizeof(int)*(m/NUM_BITS + 1));
    init(k);
        struct bloom *l = (struct bloom *)malloc(sizeof(struct bloom));
        l->bloom_vector = (int*)malloc(sizeof(int)*(m/NUM_BITS + 1));
    init(l);

        int val = rand()%M;

        int i;
        for (i=0;i<n;i++){
                int v = rand()%M;
                insert(v,k);
                v = rand()%M;
                insert(v,l);
        }

        double m_cost = 0, i_cost=0;
        struct timespec start,finish;

        clock_gettime(CLOCK_MONOTONIC,&start);
        for (i=0;i<M;i++){
                is_in(i,k);
        }
        clock_gettime(CLOCK_MONOTONIC,&finish);
        m_cost += finish.tv_sec - start.tv_sec;
        m_cost += (finish.tv_nsec - start.tv_nsec)/pow(10,9);
//      m_cost /= 100000;

        int nk = num_ones(k), nl = num_ones(l);
        clock_gettime(CLOCK_MONOTONIC,&start);
        elemsIntersection(k,l,nk,nl);
        clock_gettime(CLOCK_MONOTONIC,&finish);
        i_cost += finish.tv_sec - start.tv_sec;
        i_cost += (finish.tv_nsec - start.tv_nsec)/pow(10,9);

        int totalLevels = 1;
        int nElems = M/2;
        while (nElems>0){
                totalLevels +=1;
                nElems /= 2;
        }
//      printf("M = %d, VECT_SIZE = %d, m_cost = %lf,i_cost = %lf, totalLevels = %d\n",M,VECT_SIZE,m_cost,i_cost,totalLevels);

        nElems = M;
        int d = 1;
        while (nElems > 0){
                nElems = M/pow(2,d);
                double m1 = nElems/log(nElems);
                double m2 = (M*i_cost/m_cost);
                if (m1 < m2) break;
                d = d+1;
        }
    return d;
}

int main1(int argc, char *argv[]){
	K = 3;
	VECT_SIZE = 10000;
	printf("bestLevel = %d\n", bestLevel(10000,1000,100000));
}
/*
*Usage: ./sampleNoThread<Hash> M nSets inputFile k m T l
*/
int main(int argc, char *argv[]){
	if (argc<7){
		printf("Usage: ./BSTtesting M nSets inputFile k T desAcc\n");
		return 0;
	}
	struct timespec start,finish,mid;
	double elapsed_bst = 0.0, elapsed_da = 0.0;
	srand(time(NULL));
	setSeeds();
	int nVertices = atoi(argv[1]);
	int nSets = atoi(argv[2]);
	K = atoi(argv[4]);
	T = atoi(argv[5]);
	double desAcc = atof(argv[6]);

	int setSize = 10000;
	VECT_SIZE = getRequiredm(nVertices, setSize, K, desAcc);
	levelThreshold = bestLevel(VECT_SIZE, setSize, nVertices);
	int maxSize = 15426, nHit = 0;
	seiveInitial();
	struct bloom *k = (struct bloom *)malloc(sizeof(struct bloom));
	k->bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
	init(k);
	int startSize = malloc_count_current();
	struct bloomTree *a = getBloomTree(0,nVertices);
	int endSize = malloc_count_current();
	//printf("created the Bloom tree %d\n",(endSize - startSize));fflush(stdout);
	int *realArr = (int*)malloc(sizeof(int)*maxSize);
	char inp[100000];

	FILE *fi = fopen(argv[3],"r");
	int i = 0;
	while (i<nSets){
		init(k);
		int val=0;
		fgets(inp,100000,fi);
		char *tok = strtok(inp,",");
		int setSize = 0;
		while (tok){
			val = atoi(tok);
			realArr[setSize++] = val;
			insert(val,k);
			tok = strtok(NULL,",");
		}
		int query_ones = num_ones(k);
		int sample;
		int j = 0;
		for (j=0;j < 100; j++){
			clock_gettime(CLOCK_MONOTONIC, &start);
			sample = sampleFromTree(k,a,query_ones);
			clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed_bst += 1000*((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9));
			int z = belongsInArr1(sample, realArr, setSize);
			if (z >= 1) nHit++;
			clock_gettime(CLOCK_MONOTONIC, &start);
			sample = dictionarySample(k,nVertices);
			clock_gettime(CLOCK_MONOTONIC,&finish);
			elapsed_da += 1000*((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9));

		}
/*		clock_gettime(CLOCK_MONOTONIC,&start);
		sample = sampleFromTree(k,a,query_ones);
		clock_gettime(CLOCK_MONOTONIC,&finish);
		elapsed_bst += 1000*((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9));

//		printf("sample = %d, elapsed = %lf, query_ones = %d\n", sample, elapsed_bst, query_ones);
		clock_gettime(CLOCK_MONOTONIC,&start);
		sample = dictionarySample(k,nVertices);
		clock_gettime(CLOCK_MONOTONIC,&finish);
		elapsed_da += 1000*((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9));
		int z = belongsInArr1(sample,realArr,setSize);
		if (z>=1) nHit++;*/
		i++;
//		printf("i = %d\n",i);fflush(stdout);
	}
	fclose(fi);
//	printf("%.1f %lf %lf\n", desAcc, elapsed_bst, elapsed_da);
//	printf("%.1f %lf \n", desAcc, (1.0*nHit)/(100*nSets));
	printf("(M=%d,n=%d,m=%d,thresh=%d) : TIME_BST (ms) = %lf, ACCURACY_BST = %lf, TIME_DA (ms) = %lf\n",
			nVertices, nSets, VECT_SIZE, levelThreshold, elapsed_bst/(100*nSets), (1.0*nHit)/nSets, elapsed_da/(100*nSets));
}
