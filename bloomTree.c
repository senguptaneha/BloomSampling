#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "../bloom.h"

extern int K;
extern int VECT_SIZE;
int levelThreshold;
int T;
int nMembership = 0, numIntersections = 0;

struct bloomNode{
	struct bloom filter;
	int start;		/*starting vertex for this node*/
	int end;		/*ending vertex for this node*/
	struct bloomNode *lchild;
	struct bloomNode *rchild;
	int trueValue;	//TO BE DELETED
	int flag;
	int leafID;
        int sample;
        double prevEstimate;
	int level;
	int nOnes;
};

struct bloomTree{		/*The root level*/
	struct bloomNode* left;
	struct bloomNode* right;
};

int numNodes = 0;
int leafID = 0;
struct bloomNode *createBloomTree(int startVertex, int endVertex, int level){
	numNodes++;

	struct bloomNode *r = (struct bloomNode*)malloc(sizeof(struct bloomNode));
	r->level = level;
	r->start = startVertex;
	r->end = endVertex;
	r->lchild = NULL;
	r->rchild = NULL;
	if (level <= levelThreshold){
		r->filter.bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
		init(&r->filter);
		r->trueValue = 0;	//TO BE DELETED
		int i;

		if (level < levelThreshold){
            int mid = ((endVertex + startVertex)/2);
            r->lchild = createBloomTree(startVertex,mid,level+1);
            r->rchild = createBloomTree(mid+1,endVertex,level+1);
            for(i=0;i< (VECT_SIZE/NUM_BITS) + 1;i++)
                (r->filter).bloom_vector[i] = (r->lchild->filter).bloom_vector[i] | (r->rchild->filter).bloom_vector[i];
        }
        else{
            r->leafID = leafID;
            leafID++;
            /*Inserting actual values in the bloom filter*/
            for (i=startVertex;i<=endVertex;i++){
                insert(i,&r->filter);
            }
        }
		r->nOnes = num_ones(&r->filter);
		setSignature(&(r->filter));
	}


	return r;
}

struct bloomTree *getBloomTree(int startVertex, int endVertex){
	int mid = (endVertex + startVertex)/2;
	struct bloomTree *bT = (struct bloomTree*)malloc(sizeof(struct bloomTree));

	bT->left = createBloomTree(startVertex,mid,1);
	bT->right = createBloomTree(mid+1,endVertex,1);
	return bT;
}

/*Count the number of ones in the intersection of r and k*/
int intersectNode(struct bloomNode *r,struct bloom *k){
	int j = 0;
	int flag = 0;
	flag = numOnes_intersection(&r->filter,k);
	return flag;
}

void traverseTree(int *realArr, int nNodes, struct bloomNode *r){
	if (r->lchild)
		traverseTree(realArr,nNodes,r->lchild);
	if ((!r->lchild)&&(!r->rchild)){
		int i,num=0;
		for (i=0;i<nNodes;i++){
			if ((realArr[i] <= r->end)&&(realArr[i] >= r->start)) num++;
		}
		if (num>0) printf("(%d,%d):%d\t",r->start,r->end,num);
	}
	if (r->rchild)
		traverseTree(realArr,nNodes,r->rchild);
}

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
/*              int flagL = intersectNode(r->lchild,k);
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
	nMembership = 0; numIntersections = 0;
	double lElems = elemsIntersection(&(a->left->filter),k,a->left->nOnes,k1);
        double rElems = elemsIntersection(&(a->right->filter),k,a->right->nOnes,k1);

        int flagL = (lElems > T)?1:0;
        int flagR = (rElems > T)?1:0;

//	printf("sampleFromTree: flagL = %d, flagR = %d, lElems = %lf, rElems = %lf, a->left->nOnes = %d, a->right->nOnes = %d\n",flagL, flagR, lElems, rElems, a->left->nOnes, a->right->nOnes);
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
