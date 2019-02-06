#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "../bloom.h"

extern int K;
extern int VECT_SIZE;
int levelThreshold;
int T;
int nMembership = 0, numIntersections = 0;

struct bloomNode{
	struct bloom filter;
	long start;		/*starting vertex for this node*/
	long end;		/*ending vertex for this node*/
	struct bloomNode *lchild;
	struct bloomNode *rchild;
	struct bloomNode *parent;
	int flag;
	int leafID;
        long sample;
        double prevEstimate;
	int level;
	int nOnes;
	int isCreated;
};

struct bloomTree{		/*The root level*/
	struct bloomNode* left;
	struct bloomNode* right;
	int numNodesPopulated;
	int numTotalNodes;
};

int numNodes = 0;
int numNodesPopulated = 0;
int leafID = 0;
int leafSize = -1;
struct bloomNode *createBloomTree(long startVertex, long endVertex, int level){
	numNodes++;
	struct bloomNode *r = (struct bloomNode*)malloc(sizeof(struct bloomNode));
	r->level = level;
	r->start = startVertex;
	r->end = endVertex;
	r->lchild = NULL;
	r->rchild = NULL;
	r->isCreated = 0;
	if (level <= levelThreshold){
	}
	if (level < levelThreshold){
		long mid = ((endVertex + startVertex)/2);
		r->lchild = createBloomTree(startVertex,mid,level+1);
		r->rchild = createBloomTree(mid+1,endVertex,level+1);
	}
	else{
                r->leafID = leafID;
                leafID++;
		if (leafSize == -1)
			leafSize = r->end - r->start + 1;
        }

	return r;
}

void setParentPointers(struct bloomNode *bT){
	if (bT == NULL) return;
	if (bT->lchild){
		(bT->lchild)->parent = bT;
		setParentPointers(bT->lchild);
	}
	if (bT->rchild){
		(bT->rchild)->parent = bT;
		setParentPointers(bT->rchild);
	}
}

struct bloomTree *getBloomTree(long startVertex, long endVertex){
	int mid = (endVertex + startVertex)/2;
	struct bloomTree *bT = (struct bloomTree*)malloc(sizeof(struct bloomTree));

	bT->left = createBloomTree(startVertex,mid,1);
	bT->right = createBloomTree(mid+1,endVertex,1);
	bT->numNodesPopulated = 0;
	bT->numTotalNodes = numNodes;
	setParentPointers(bT->left);
	setParentPointers(bT->right);
	(bT->left)->parent = NULL;
	(bT->right)->parent = NULL;
	return bT;
}

struct leafSelector{
	struct bloomNode *leaf;
	int isSelected;
	double prob;
};

struct leafSelector *leafList;

void populateLeafList(struct bloomNode *bT){
	if ((bT->lchild == NULL) && (bT->rchild == NULL)){
		int leafId = bT->leafID;
		leafList[leafId].leaf = bT;
		return;
	}
	populateLeafList(bT->lchild);
	populateLeafList(bT->rchild);
}

void selectUniform(double ratio, int numLeaves){
	int i;
	for (i = 0; i < numLeaves; i++){
		double rVal = ((double) rand())/((double) RAND_MAX);
		if (rVal <= ratio)
			leafList[i].isSelected = 1;
	}

}

int selectSkewed(double ratio, int numLeaves){
	int leavesToSelect = (int)(floor(ratio*numLeaves));
	int numLeavesSelected = 0, i;

	while (numLeavesSelected < leavesToSelect){
		double rVal = ((double) rand())/((double) RAND_MAX), cdf = 0;
		int index = -1;
		for (i=0; i<numLeaves; i++){
			cdf += leafList[i].prob;
			if ((cdf <= rVal)&&( leafList[i].isSelected == 0)){
				if (leafList[i].isSelected == 0) numLeavesSelected += 1;
				leafList[i].isSelected = 1;
//				numLeavesSelected += 1;
				index = i;
				break;
			}
		}
		if (index == -1) return 1; //error
		double prob = leafList[index].prob;
		int lIndex = index - 1;
		if (lIndex < 0) lIndex = numLeaves - 1;
		while (leafList[lIndex].prob == 0.0){
			lIndex -= 1;
			if (lIndex < 0) lIndex = numLeaves - 1;
		}
		int rIndex = index + 1;
		if (rIndex >= numLeaves) rIndex = 0;
		while (leafList[rIndex].prob == 0.0){
			rIndex += 1;
			if (rIndex >= numLeaves) rIndex = 0;
		}
		if ((lIndex >= 0) && (rIndex < numLeaves)){
			leafList[index].prob = 0.0;
			for (i = 0; i < numLeaves; i++){
				if ((i != index) && (i != lIndex) && (i != rIndex)){
					double p = 0.1* leafList[i].prob;
					leafList[i].prob -= p;
					prob += p;
				}
			}
			leafList[lIndex].prob += prob/2.0;
			leafList[rIndex].prob += prob/2.0;
		}
		else return 1;	//error
	}
//	printf("OUT: numLeavesSelected = %d\n",numLeavesSelected);
	return 0;
}

void populateBTByLeafList(int numLeaves){
	int i = 0;
	int numLeavesSelected = 0;
	for (i = 0; i < numLeaves; i++){
		if (leafList[i].isSelected == 0) continue;
		struct bloomNode *leaf = leafList[i].leaf;
		printf("--------------%d--------------\n",numLeavesSelected);
		fflush(stdout);
		numLeavesSelected += 1;
		while (leaf){
			if (!leaf->isCreated){
				printf("Populating %ld, %ld\n",leaf->start,leaf->end);
				fflush(stdout);
				int numStepsNoChange = 0;
				leaf->filter.bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
                		init(&leaf->filter);
                		long i;
                		/*Inserting actual values in the bloom filter*/
				int prevOnes = (leaf->filter).nOnes, newOnes = 0;
                		for (i=leaf->start;i<=leaf->end;i++){
                        		insert(i,&leaf->filter);
					newOnes = (leaf->filter).nOnes;
					if (newOnes == VECT_SIZE){
						printf("Breaking at %ld\n",i);
						break;
					}
					if (newOnes == prevOnes){
						numStepsNoChange += 1;
						if (numStepsNoChange >= 100000) break;
					}
					else numStepsNoChange = 0;
					prevOnes = newOnes;
                		}
                		leaf->nOnes = num_ones(&leaf->filter);
                		setSignature(&(leaf->filter));
                		leaf->isCreated = 1;
                		numNodesPopulated += 1;
			}
			leaf = leaf->parent;
		}
	}
}

void sampleAndPopulateLeaves(struct bloomTree *bT, int method, double ratio){
	printf("Sampling for method = %d, ratio = %lf\n",method,ratio);
	fflush(stdout);
	/*method = 0 for uniform, 1 for clustered*/
	int numLeaves = 1,i;
	for (i = 0; i < levelThreshold; i++) numLeaves *= 2;
	leafList = (struct leafSelector*)malloc(numLeaves *sizeof(struct leafSelector));
	double prob = 1.0/numLeaves;
	for (i = 0; i < numLeaves; i++){
		leafList[i].isSelected = 0;
		leafList[i].prob = prob;
	}

	populateLeafList(bT->left);
	populateLeafList(bT->right);


	if (method==0) selectUniform(ratio, numLeaves);
	else if (method == 1) {
		int z = 1;
		while (z!=0){
			z = selectSkewed(ratio, numLeaves);
			if (z!=0){
				for (i = 0; i < numLeaves; i++){
					leafList[i].isSelected = 0;
					leafList[i].prob = prob;
				}
			}
		}
	}

	populateBTByLeafList(numLeaves);
	/*remember to free leafList*/
	free(leafList);
	printf("Done Sampling for method = %d, ratio = %lf\n",method,ratio);
}

int sampleLeavesOnly(struct bloomTree *bT, int method, double ratio){
	printf("Sampling for method = %d, ratio = %lf\n",method,ratio);
	fflush(stdout);
	/*method = 0 for uniform, 1 for clustered*/
	int numLeaves = 1,i;
	for (i = 0; i < levelThreshold; i++) numLeaves *= 2;
	leafList = (struct leafSelector*)malloc(numLeaves *sizeof(struct leafSelector));
	double prob = 1.0/numLeaves;
	for (i = 0; i < numLeaves; i++){
		leafList[i].isSelected = 0;
		leafList[i].prob = prob;
	}

	populateLeafList(bT->left);
	populateLeafList(bT->right);


	if (method==0) selectUniform(ratio, numLeaves);
	else if (method == 1) {
		int z = 1;
		while (z!=0){
			z = selectSkewed(ratio, numLeaves);
			if (z!=0){
				for (i = 0; i < numLeaves; i++){
					leafList[i].isSelected = 0;
					leafList[i].prob = prob;
				}
			}
		}
	}
	int numNodes = 0;
	int numLeavesSelected = 0;
	for (i=0;i<numLeaves;i++)
		if (leafList[i].isSelected){
			numLeavesSelected += 1;
			struct bloomNode *leaf = leafList[i].leaf;
			while (leaf){
				if (leaf->isCreated == 0){
					leaf->isCreated = 1;
					numNodes += 1;
				}
				leaf = leaf->parent;
			}
		}
	/*remember to free leafList*/
	free(leafList);
	printf("Done Sampling for method = %d, ratio = %lf\n",method,ratio);
/*
	if (method == 1)
		printf("OUT:: numLeavesSelected = %d\n",numLeavesSelected);
*/
	return numNodes;
}

int searchInTree(struct bloomTree *bT, long value){
	if (value <  (bT->left)->start) return 0;
	if (value > (bT->right)->end) return 0;
	struct bloomNode *b = bT->left;
	if (value > (bT->left)->end) b = bT->right;

	while ((b->lchild) && (b->rchild)){
		if (value > (b->lchild)->end)
			b = b->rchild;
		else	b = b->lchild;
	}
	if (b->isCreated == 1) return 1;
	return 0;
}

void insertInNode(long value, struct bloomNode *a){
        if (a==NULL) return;    /*Only a protection mechanism*/
        if (a->isCreated == 0){
		a->filter.bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
                init(&a->filter);
                long i;
                /*Inserting actual values in the bloom filter*/
                for (i=a->start;i<=a->end;i++){
                        insert(i,&a->filter);
//			if (i%1000000 == 0) printf("done: %ld/%ld\n ",i,a->end);
                }
                a->nOnes = num_ones(&a->filter);
                setSignature(&(a->filter));
                a->isCreated = 1;
		numNodesPopulated += 1;
        }
//	printf("Done node (%ld, %ld)\n",a->start,a->end);
        if (a->lchild && a->rchild){    /*Not a leaf node*/
                if ((value >= (a->lchild)->start) && (value <= (a->lchild)->end)){
                        insertInNode(value, a->lchild);
                }
                else if ((value >= (a->rchild)->start) && (value <= (a->rchild)->end)){
                        insertInNode(value, a->rchild);
                }
        }
}

void insertInTree(long value, struct bloomTree *a){
	numNodesPopulated = 0;
//	printf("Yes, in here, value = %ld, (a->left)->start = %ld, (a->left)->end = %ld, (a->right)->start = %ld, (a->right)->end = %ld\n",
//	value,(a->left)->start,(a->left)->end,(a->right)->start, (a->right)->end);
        if ((value >= (a->left)->start) && (value <= (a->left)->end)){
                insertInNode(value,a->left);
        }
        else if ((value >= (a->right)->start) && (value <= (a->right)->end)){
                insertInNode(value,a->right);
        }
	a->numNodesPopulated += numNodesPopulated;
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
		if (num>0) printf("(%ld,%ld):%d\t",r->start,r->end,num);
	}
	if (r->rchild)
		traverseTree(realArr,nNodes,r->rchild);
}

long sampleValue(struct bloom *k, struct bloomNode *r, int k1){
	if (r->isCreated == 0) return -1;
        if ((!r->lchild)&&(!r->rchild)){
        //if (r->level >= levelThreshold){
                long j;
                long numValues = r->end - r->start + 1;
		int nCandidates = 0;
		long sample = -1;
                for (j=r->start;j<=r->end;j++){
                        nMembership++;
                        if (is_in(j,k)){
				nCandidates++;
				int rval = rand()%nCandidates;
				if (rval ==0) sample = j;
                        }
                }

                return sample;
        }
        else{
		double lElems = 0;
		double rElems = 0;
		if (r->lchild->isCreated) lElems = elemsIntersection(&(r->lchild->filter), k, r->lchild->nOnes ,k1);
                if (r->rchild->isCreated) rElems = elemsIntersection(&(r->rchild->filter), k, r->rchild->nOnes, k1);
                numIntersections += 2;
                int flagL = (lElems > T)?1:0;
                int flagR = (rElems > T)?1:0;

                if ((flagL)&&(!flagR)) return sampleValue(k,r->lchild,k1);
                if ((!flagL)&&(flagR)) return sampleValue(k,r->rchild,k1);
                if ((!flagL)&&(!flagR)) {
                        return -1;
                }


                double z = ((double)rand())/((double)RAND_MAX)*(lElems+rElems);
                if (z<lElems){
                        long lSample = sampleValue(k,r->lchild,k1);
                        if (lSample==-1) return sampleValue(k,r->rchild,k1);
                        return lSample;
                }
                long rSample = sampleValue(k,r->rchild,k1);
                if (rSample==-1) return sampleValue(k,r->lchild,k1);
                else return rSample;
        }
        return -1;
}

long sampleFromTree(struct bloom *k, struct bloomTree *a, int k1){
	double lElems = 0;
	double rElems = 0;
	if (a->left->isCreated == 1) lElems = elemsIntersection(&(a->left->filter),k,a->left->nOnes,k1);
        if (a->right->isCreated == 1) rElems = elemsIntersection(&(a->right->filter),k,a->right->nOnes,k1);

        int flagL = (lElems > T)?1:0;
        int flagR = (rElems > T)?1:0;


        if ((flagL)&&(!flagR)) return sampleValue(k,a->left,k1);
        if ((!flagL)&&(flagR)) return sampleValue(k,a->right,k1);
        if ((!flagL)&&(!flagR)){
                return -1;
        }
	double z = ((double)rand())/((double)RAND_MAX)*(lElems+rElems);
        if (z<lElems){

                long lSample = sampleValue(k,a->left,k1);
                if (lSample==-1) return sampleValue(k,a->right,k1);
                return lSample;
        }
        else{
                long rSample = sampleValue(k,a->right,k1);
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

void resetNode(struct bloomNode *bT){
	if (bT->lchild){
		resetNode(bT->lchild);
	}
	if (bT->rchild){
		resetNode(bT->rchild);
	}
	if(bT->isCreated){
		free(bT->filter.bloom_vector);
		bT->isCreated = 0;
		bT->nOnes = 0;
	}
}

void resetNodeNoFree(struct bloomNode *bT){
	if (bT->lchild){
		resetNodeNoFree(bT->lchild);
	}
	if (bT->rchild){
		resetNodeNoFree(bT->rchild);
	}
	if(bT->isCreated){
		bT->isCreated = 0;
		bT->nOnes = 0;
	}
}

void resetTree(struct bloomTree *bT){
	resetNode(bT->left);
	resetNode(bT->right);
}

void resetTreeNoFree(struct bloomTree *bT){
	resetNodeNoFree(bT->left);
	resetNodeNoFree(bT->right);
}

struct bloomNode *hasEmptyIntersection(struct bloomNode *bN, struct bloom *k, int k1, double threshold){
	if ((bN->lchild == NULL) && (bN->rchild == NULL)) return NULL;	//this is a leaf node, an empty intersection here doesn't help
	if (bN->isCreated == 0) return NULL;

	double estElems = elemsIntersection(&(bN->filter), k, bN->nOnes ,k1);
	if (estElems <= threshold) return bN;
	struct bloomNode *bLeft = hasEmptyIntersection(bN->lchild, k, k1, threshold);
	if (bLeft) return bLeft;
	struct bloomNode *bRight = hasEmptyIntersection(bN->rchild, k, k1, threshold);
	return bRight;
}

struct bloomNode *emptyIntersectionTree(struct bloomTree *bT, struct bloom *k, double threshold){
	int k1 = num_ones(k);

	struct bloomNode *bN = NULL;
	if ((bT->left)->isCreated) bN = hasEmptyIntersection(bT->left, k, k1, threshold);
	if ((bN == NULL) && ((bT->right)->isCreated)) bN = hasEmptyIntersection(bT->right, k, k1, threshold);
	return bN;
}

int numDescendantLeaves(struct bloomNode *bN){
	if (bN->isCreated == 0) return 0;
	if ((bN->lchild == NULL) && (bN->rchild == NULL)) return 1;
	int leftLeaves = numDescendantLeaves(bN->lchild);
	int rightLeaves = numDescendantLeaves(bN->rchild);
	return (leftLeaves + rightLeaves);
}
