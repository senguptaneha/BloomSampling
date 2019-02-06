#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "bloomTree.c"
#include "../bloom.h"
extern int K;
extern int VECT_SIZE;
int T;
extern int levelThreshold;
int nMembership = 0;
int numIntersections = 0;

int sampleValue(struct bloom *k, struct bloomNode *r, int k1){
	if ((!r->lchild)&&(!r->rchild)){
	//if (r->level >= levelThreshold){
		int i=0,j;
		int numValues = r->end - r->start + 1;
		int *values = (int*)malloc(sizeof(int)*numValues);
		for (j=r->start;j<=r->end;j++){
			nMembership++;
			if (is_in(j,k)){
				values[i] = j;
				i++;
			}
		}
		if (i==0){
			free(values);
			return -1;
		}
		int index = rand()%i;

		int z = values[index];
		free(values);

		return z;
	}
	else{
/*		int flagL = intersectNode(r->lchild,k);
		int flagR = intersectNode(r->rchild,k);
		double lElems = -1.0*log(1 - ((double) flagL)/VECT_SIZE);
		double rElems = -1.0*log(1 - ((double) flagR)/VECT_SIZE);
*/
		double lElems = elemsIntersection(&(r->lchild->filter), k, r->lchild->nOnes ,k1);
		double rElems = elemsIntersection(&(r->rchild->filter), k, r->rchild->nOnes, k1);
		//printf("(%d,%d):(%lf,%lf), (%d,%d)\n",r->start,r->end,lElems,rElems,r->lchild->trueValue,r->rchild->trueValue);
		numIntersections += 2;
		//int flagL = intersectNode(r->lchild,k);
		//int flagR = intersectNode(r->rchild,k);

		//double lElems = r->lchild->trueValue;
		//double rElems = r->rchild->trueValue;
		int flagL = (lElems > T)?1:0;
		int flagR = (rElems > T)?1:0;

		if ((flagL)&&(!flagR)) return sampleValue(k,r->lchild,k1);
		if ((!flagL)&&(flagR)) return sampleValue(k,r->rchild,k1);
		if ((!flagL)&&(!flagR)) {
			return -1;
		}

		double z = ((double)rand())/((double)RAND_MAX)*(lElems+rElems);
		if (z<lElems){
			int lSample = sampleValue(k,r->lchild,k1);
			if (lSample==-1) return sampleValue(k,r->rchild,k1);
			return lSample;
		}
		int rSample = sampleValue(k,r->rchild,k1);
		if (rSample==-1) return sampleValue(k,r->lchild,k1);
		else return rSample;
	}
	return -1;
}

int sampleFromTree(struct bloom *k, struct bloomTree *a, int k1){
/*	int flagL = intersectNode(a->left,k);
	int flagR = intersectNode(a->right,k);*/
	//int flagL = a->left->flag;
	//int flagR = a->right->flag;
	double lElems = elemsIntersection(&(a->left->filter),k,a->left->nOnes,k1);
	double rElems = elemsIntersection(&(a->right->filter),k,a->right->nOnes,k1);

	int flagL = (lElems > T)?1:0;
	int flagR = (rElems > T)?1:0;


	if ((flagL)&&(!flagR)) return sampleValue(k,a->left,k1);
	if ((!flagL)&&(flagR)) return sampleValue(k,a->right,k1);
	if ((!flagL)&&(!flagR)){
		return -1;
	}
/*
	int z = rand()%2;
	if (z==0){
*/
/*
	int z = rand()%(flagL + flagR);
	if (z<flagL){
*/
	double z = ((double)rand())/((double)RAND_MAX)*(lElems+rElems);
	if (z<lElems){

		int lSample = sampleValue(k,a->left,k1);
		if (lSample==-1) return sampleValue(k,a->right,k1);
		return lSample;
	}
	else{
		int rSample = sampleValue(k,a->right,k1);
		if (rSample==-1) return sampleValue(k,a->left,k1);
		return rSample;
	}

	return -1;
}

int belongsInArr1(int val, int *arr, int size){
	int j;
	for (j=0;j<size;j++)
		if (arr[j]==val)
			return j;
	return -1;
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

/*
*Usage: ./sampleNoThread<Hash> M n inputFile k m T l
*/
int main(int argc, char *argv[]){
	struct timespec start,finish,mid;
	double elapsed=0.0;
	srand(time(NULL));
	setSeeds();
	int nVertices = atoi(argv[1]);
	int nNodes = atoi(argv[2]);
	K = atoi(argv[4]);
	VECT_SIZE = atoi(argv[5]);
	T = atoi(argv[6]);
	levelThreshold = atoi(argv[7]);

	seiveInitial();
	struct bloom *k = (struct bloom *)malloc(sizeof(struct bloom));
	k->bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
	init(k);
	struct bloomTree *a = getBloomTree(0,nVertices);
	int *realArr = (int*)malloc(sizeof(int)*nNodes);


	FILE *fi = fopen(argv[3],"r");
	int i = 0;
	while (i<nNodes){
		int val=0;
		fscanf(fi,"%d\n",&val);
		insert(val,k);
		realArr[i] = val;
		updateTrueValue(a->left,val);
		updateTrueValue(a->right,val);
		i++;
	}
	setSignature(k);
	//printSignature(k);
	fclose(fi);
	int query_ones = num_ones(k);
	/*clock_gettime(CLOCK_MONOTONIC,&start);
	int sample = sampleFromTree(k,a);
	clock_gettime(CLOCK_MONOTONIC,&finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec)/pow(10.0,9);
	printf("TIME TAKEN FOR SAMPLING(%d,%d) = %lf\n",belongsInArr1(sample,realArr,nNodes),sample,elapsed);*/
	int nSamples = 1000, nHit = 0,nInv = 0;
	int *obs = (int*)malloc(sizeof(int)*nNodes);
	for (i=0;i<nNodes;i++) obs[i] = 0;

	for (i=0;i<nSamples;i++){
		clock_gettime(CLOCK_MONOTONIC,&start);
		int sample = sampleFromTree(k,a,query_ones);
		if (sample==-1) nInv++;
		clock_gettime(CLOCK_MONOTONIC,&finish);
		elapsed += finish.tv_sec - start.tv_sec;
		elapsed += (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		int z = belongsInArr1(sample,realArr,nNodes);
		if (z>=0){
			nHit++;
			obs[z] += 1;
		}
		//printf("%d \t %lf\n",i,(1.0*nHit)/i);
		//fflush(stdout);
	}
	//double accuracy = (1.0*nHit)/nSamples;
	double expectedValue = (1.0*nHit)/nNodes, chi_square = 0.0;
	/*FILE *fout = fopen("Frequencies","w");*/
	for (i=0;i<nNodes;i++){
		double o = obs[i], e = expectedValue;
		double c = (o - e)*(o - e)/e;
		chi_square += c;
//		fprintf(fout,"%d\t%d\n",(i+1),obs[i]);
	}
//	fclose(fout);*/
	/*for (i=0;i<nVertices;i++){
		int j = belongsInArr1(i,realArr,nNodes);
		if (j>=0) printf("%d\t%d\n",i,obs[j]);
		else printf("%d\t%d\n",i,0);
	}*/
	//double chi_square = 0.0;
	printf("(M=%d,n=%d,m=%d,thresh=%d) : TIME TAKEN = %lf,CHI_SQUARE = %lf,ACCURACY = %lf,INTERSECTIONS = %lf,MEMBERSHIP = %lf\n",
			nVertices,nNodes,VECT_SIZE,levelThreshold,elapsed/nSamples,chi_square,(1.0*nHit)/nSamples,(1.0*numIntersections)/nSamples, (1.0*nMembership)/nSamples);

	//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",levelThreshold,elapsed/nSamples,chi_square, (1.0*nHit)/nSamples, (1.0*numIntersections)/nSamples + 2, (1.0*nMembership)/nSamples);
/*
	int nSamples = 130*nNodes, nHits = 0, *obs;
	obs = (int*)malloc(sizeof(int)*nNodes);
	for (i=0;i<nNodes;i++) obs[i] = 0;

	for (i=0;i<nSamples;i++){
		int sample = sampleFromTree(k,a);
		int z = belongsInArr1(sample,realArr,nNodes);
		if (z>=0){
			nHits++;
			obs[z] += 1;
		}
	}

	double expectedValue = (1.0*nHits)/nNodes, chi_square = 0.0;

	for (i=0;i<nNodes;i++){
		double o = obs[i], e = expectedValue;
		double c = (o - e)*(o - e)/e;
		chi_square += c;
	}
	printf("VALUE OF CHI SQUARE = %lf\n",chi_square);
*/
}
