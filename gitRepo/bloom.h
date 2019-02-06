#ifndef BLOOM_H
#define BLOOM_H
#define NUM_BITS 32

struct bloom{
	int *bloom_vector;
	long signature;	//view as a bit array, stores a 0 for a group of integers that are zero in bloom_vector.
	int nOnes;
};


#define a1 959
#define b1 925
#define a2 652
#define b2 653
#define a3 948
#define b3 490
#define ROT32(x,n)      ((x << n) | (x >> (32 - n)))

void setSeeds();
void seiveInitial();
void init(struct bloom *bl);
void setSignature(struct bloom *a);
void insert(long val, struct bloom* bl);
int is_in(long val, struct bloom* bl);
int is_inHash1(int val, struct bloom *bl);
int is_inHash2(int val, struct bloom *bl);
int is_inHash3(int val, struct bloom *bl);
int hash1(int a);
int hash2(int a);
int hash3(int a);
struct bloom * intersect_bloom(struct bloom *a, struct bloom *b);
struct bloom* union_bloom(struct bloom *a, struct bloom *b);
int dot_product(struct bloom *a, struct bloom *b);
int num_zero(struct bloom *a);
long num_ones(struct bloom *a);
int numOnes_intersection(struct bloom *a, struct bloom *b);
double elemsIntersection(struct bloom *a, struct bloom *b, int ta, int tb);
long num_elem(struct bloom *a);
int getInverse(int b, int m);
int chineseRemainder(int *a, int *m);
int solveVal(int v1, int v2, int v3);
int *reservoirSample(int k, struct bloom *a);
int sampleFromBloom(struct bloom *a, int nvertices);
int sampleFromBloom1(struct bloom *a, int nvertices);
int count_ones(int x);
#endif
