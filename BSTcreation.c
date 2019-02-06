#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
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
//	m_cost /= 100000;

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
//	printf("M = %d, VECT_SIZE = %d, m_cost = %lf,i_cost = %lf, totalLevels = %d\n",M,VECT_SIZE,m_cost,i_cost,totalLevels);

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

struct tt{
	double elapsedS;
	double elapsedU;
};

void getTimeBST(int M, int n, int k, double a, struct tt *r){
    struct rusage usageS, usageE;
//    struct timeval start, finish;
    double elapsedS, elapsedU;
    K = k;
    VECT_SIZE = getRequiredm(M,n,k,a);
    T = 0;
//    printf("bestLevel: %d, %d, %d, %lf\n ", VECT_SIZE, n, M, a); fflush(stdout);
    levelThreshold = bestLevel(VECT_SIZE, n, M);
  //  printf("levelThreshold = %d\n",levelThreshold);
    seiveInitial();
    getrusage(RUSAGE_SELF, &usageS);
    //clock_gettime(CLOCK_MONOTONIC, &start);
    struct bloomTree *bT = getBloomTree(0,M);
    //clock_gettime(CLOCK_MONOTONIC, &finish);
    getrusage(RUSAGE_SELF, &usageE);

    //elapsed = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
    r->elapsedS = 1000*((usageE.ru_stime.tv_sec - usageS.ru_stime.tv_sec) +(usageE.ru_stime.tv_usec - usageS.ru_stime.tv_usec)/pow(10,6));
    r->elapsedU = 1000*((usageE.ru_utime.tv_sec - usageS.ru_utime.tv_sec) +(usageE.ru_utime.tv_usec - usageS.ru_utime.tv_usec)/pow(10,6));

//    return elapsed;
}
/*
*Usage: ./sampleNoThread<Hash> M nSets inputFile k m T l
*/
int main(int argc, char *argv[]){
	struct tt *r = (struct tt *)malloc(sizeof(struct tt));
	double acc;
	srand(time(NULL));
	setSeeds();
	int numRounds = 100;
	int M = 10000, n = 10000, k = 3, i, j;
	printf("M VECT_SIZE L acc meanS stdDevS meanU stdDevU\n");

	for (i = 5; i <= 7; i++){
	        M = 10 * M;
        	for (acc = 0.5; acc <= 0.9; acc += 0.1){
        		double meanS = 0, stdDevS = 0, meanU = 0, stdDevU = 0;
            		for (j = 0; j < numRounds; j++){
                		getTimeBST(M,n,k,acc,r);
	                	meanS += r->elapsedS; meanU += r->elapsedU;
        	        	stdDevS += r->elapsedS * r->elapsedS; stdDevU += r->elapsedU * r->elapsedU;
            		}
            		meanS /= numRounds; meanU /= numRounds;
	            	stdDevS /= numRounds; stdDevU /= numRounds;
        	    	stdDevS = pow(stdDevS - meanS*meanS, 0.5);
            		stdDevU = pow(stdDevU - meanU*meanU, 0.5);
            		printf("%d\t%d\t%d\t%.1f\t%lf\t%lf\t%lf\t%lf\n",M,VECT_SIZE,levelThreshold,acc,meanS,stdDevS,meanU,stdDevU);
            		fflush(stdout);
        	}
	}
}
