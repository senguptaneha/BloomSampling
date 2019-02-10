#include <stdio.h>
#include "bloomTree.c"
#include "bloom.h"
#include <pthread.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define NUMTHREADS 48
#define QSIZE 100000

extern int K;
extern int VECT_SIZE;
/*Add conditions on queue(for correctness), and on mutexes (for performance)*/

int levelThreshold;
double overlapThreshold;

int numElements = 0;
pthread_mutex_t m_numElements;

int totalRange = 0;
pthread_mutex_t m_totalRange;
pthread_cond_t c_overflow = PTHREAD_COND_INITIALIZER;

int nIntersections = 0;
pthread_mutex_t m_numIntersections;

int numMembership = 0;
pthread_mutex_t m_numMembership;

int k1;	//Number of 1s in the query bloom filter

FILE *foutput;
pthread_mutex_t m_output;

int **output;

int nVertices;

struct bloom query;

struct task{
	struct bloomNode *root;
	struct task *next;
};

struct bloomNode *taskQueue[QSIZE];
pthread_mutex_t m_taskQueue;

int rear = -1;	//rear == front indicates empty queue
int front = -1;

pthread_t ids[NUMTHREADS];

struct bloomNode *extractTask(){
	pthread_mutex_lock(&m_taskQueue);
	if (rear==front) {
		pthread_mutex_unlock(&m_taskQueue);
		return NULL;
	}
	struct bloomNode *t = taskQueue[front];
	front = (front+1)%QSIZE;
	pthread_mutex_unlock(&m_taskQueue);
	pthread_cond_broadcast(&c_overflow);
	return t;
}

void addTask(struct bloomNode *t){
	pthread_mutex_lock(&m_taskQueue);
	if ((rear+1)%QSIZE == front){
		pthread_cond_wait(&c_overflow,&m_taskQueue);
		taskQueue[rear] = t;
		rear = (rear+1)%QSIZE;
		pthread_mutex_unlock(&m_taskQueue);
	}
	else if (rear==front){	//empty queue
		rear = 0;
		front = 0;
		taskQueue[rear] = t;
		rear = (rear+1)%QSIZE;
		pthread_mutex_unlock(&m_taskQueue);
	}
	else{
		taskQueue[rear] = t;
		rear = (rear+1)%QSIZE;
		pthread_mutex_unlock(&m_taskQueue);
	}
	// pthread_mutex_unlock(&m_taskQueue);
}

void searchLeaf(struct bloomNode *t){
	int i = 0;
	int localCount = 0;
	int leafID = t->leafID;

	int localMembership = 0;
	//pthread_mutex_lock(&m_output);
	for (i=t->start;i<=t->end;i++){
		localMembership++;
		if (is_in(i,&query)){
			// printf("OVER HERE %d. %d\n",leafID,localCount);
			// fflush(stdout);
			output[leafID][localCount] = i;
			// printf("DONE WRITING TO %d, %d\n",leafID,localCount);
			// fflush(stdout);
			localCount++;
			// printf("%d,",i);
		}
	}
	//fflush(foutput);
	//pthread_mutex_unlock(&m_output);
	pthread_mutex_lock(&m_numElements);
	numElements += localCount;
	pthread_mutex_unlock(&m_numElements);

	pthread_mutex_lock(&m_numMembership);
	numMembership += localMembership;
	pthread_mutex_unlock(&m_numMembership);
}

void descendLevels(struct bloomNode *root, int numLevels, int currLevel, double estimate){
	if (currLevel == numLevels){
		root->prevEstimate = estimate;
		addTask(root);
	}
	else if ((root->lchild==NULL) && (root->rchild==NULL)){
		root->prevEstimate = estimate;
		addTask(root);
	}
	else{
		descendLevels(root->lchild,numLevels,currLevel+1,estimate);
		descendLevels(root->rchild,numLevels,currLevel+1,estimate);
	}
}

void *executeTask(void *x){
	while (1){
		struct bloomNode *t = extractTask();
		if (t == NULL){
			pthread_mutex_lock(&m_totalRange);
			int tR = totalRange;
			pthread_mutex_unlock(&m_totalRange);
			if (tR>=nVertices)
				pthread_exit(NULL);
			//printf("I AM WAITING FOR A TASK %d, %d\n",((int)x),tR);
			//fflush(stdout);
		}
		else{
			//if ((!t->lchild)&&(!t->rchild)){
			if (t->level >= levelThreshold){
				searchLeaf(t);
				//printf("I SEARCHED A LEAF\n");
				pthread_mutex_lock(&m_totalRange);
				totalRange += (t->end - t->start + 1);
				//printf("totalRange = %d, (front = %d, rear = %d)\n",totalRange,front,rear);
				//fflush(stdout);
				pthread_mutex_unlock(&m_totalRange);

				continue;
			}

			pthread_mutex_lock(&m_numIntersections);
			nIntersections++;
			pthread_mutex_unlock(&m_numIntersections);

			double estimate = elemsIntersection(&(t->lchild->filter),&query,t->nOnes,k1);
			printf("estimate is %lf for (%d,%d)\n",estimate,t->start,t->end);
			int flag = (estimate>overlapThreshold)?1:0;
			if (flag>0){
				double ratio = estimate/t->prevEstimate;
				int numLevels;

				if (ratio <= 0.25) numLevels = 1;
				else if (ratio <= 0.5) numLevels = 2;
				else if (ratio <= 0.75) numLevels = 3;
				else numLevels = 4;

				descendLevels(t,numLevels,0,estimate);
			}
			else{
				pthread_mutex_lock(&m_totalRange);
                                totalRange += (t->end - t->start + 1);
				//printf("TotalRange = %d\n",totalRange);
				//fflush(stdout);
                                pthread_mutex_unlock(&m_totalRange);

			}
		}
	}
}


/*****************************config functions ********************************/
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
	// m_cost /= 100000;

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
	// printf("M = %d, VECT_SIZE = %d, m_cost = %lf,i_cost = %lf, totalLevels = %d\n",M,VECT_SIZE,m_cost,i_cost,totalLevels);

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
/********************************************************************************************************************/

int main(int argc, char *argv[]){
	if (argc<7){
		printf("Usage: ./BSTtesting M nSets inputFile k T desAcc\n");
		return 0;
	}
	
	double elapsed_bst = 0.0, elapsed_da = 0.0;
	srand(time(NULL));
	setSeeds(); //required for MD5 or murmur hash functions
	int nVertices = atoi(argv[1]);
	int nSets = atoi(argv[2]);
	K = atoi(argv[4]);
	overlapThreshold = atoi(argv[5]);
	double desPrec = atof(argv[6]);

	int setSize = 10000;

	VECT_SIZE = getRequiredm(nVertices, setSize, K, desPrec);
	levelThreshold = bestLevel(VECT_SIZE, setSize, nVertices);

	seiveInitial();

	struct bloomTree *a = getBloomTree(0,nVertices);
	query.bloom_vector = (int*)malloc(sizeof(int)*(VECT_SIZE/NUM_BITS + 1));
	init(&query);

	int *realArr = (int*)malloc(sizeof(int)*setSize);
	int j;
	for (j=0;j<setSize;j++) realArr[j] = 0;

	FILE *fi = fopen(argv[3],"r");
	int i = 0;
	while (i<setSize){
		int val=0;
		fscanf(fi,"%d\n",&val);
		printf("inserting %d\n", val);
		insert(val,&query);
		realArr[i] = val;
		i++;
	}
	fclose(fi);
	int rc;

	int actualElems = 0;
	for (j = 0; j < nVertices; j++)
		if(is_in(j, &query)) actualElems += 1;
	printf("acutal elems = %d\n", actualElems);
	
	(a->left)->prevEstimate = 1.0*setSize;
	(a->right)->prevEstimate = 1.0*setSize;
	/*Initialize all mutexes and condition variables*/
	pthread_mutex_init(&m_numElements, NULL);
	pthread_mutex_init(&m_taskQueue, NULL);
	pthread_mutex_init(&m_totalRange,NULL);

	pthread_mutex_init(&m_numIntersections,NULL);
	pthread_mutex_init(&m_numMembership,NULL);
//	pthread_mutex_init(&m_output,NULL);
	addTask(a->left);
	addTask(a->right);

	int M = nVertices;
	int nLeaves = 1;
	int currLevel = 0;
	while (currLevel<levelThreshold){
		nLeaves *= 2;
		M = M/2;
		currLevel += 1;
	}
	output = (int**)malloc(sizeof(int*)*nLeaves);
	for (i=0;i<nLeaves;i++){
		output[i] = (int*)malloc(sizeof(int)*M);
	}

	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC,&start);
	for (i=0;i<NUMTHREADS;i++){
		rc = pthread_create(&ids[i],NULL,executeTask,(void*)(i));
	}
	void *status;
	for (i=0;i<NUMTHREADS;i++){
		pthread_join(ids[i],&status);
	}
	clock_gettime(CLOCK_MONOTONIC,&finish);
//	fclose(foutput);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec)/pow(10,9);
	printf("(M=%d,n=%d,m=%d,o=%lf): RECONSTRUCTED SET OF SIZE %d  WITH PRECISION %lf IN TIME %lf sec, INTERSECTIONS=%d,MEMBERSHIP=%d\n",nVertices,setSize,VECT_SIZE,overlapThreshold,numElements,(1.0*setSize)/numElements,elapsed,nIntersections,numMembership-1);
	pthread_exit(NULL);
}
