#include <stdio.h>
#include "bloom.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>	//for murmur hash function
#include "md5.c"

#define BLOOM_DEBUG 0
#define HASHTYPE 0	//0 - Simple, 1 - Murmur, 2 - MD5

#define D2 89
#define A2 567
#define D3 123
#define A3 239

extern int K=0;
extern int VECT_SIZE=40000;
extern int DEBUG = 0;

int seed[10];
//int aHash[10]={369, 448, 918, 211, 474, 944, 554, 271, 69, 747};
//int bHash[10]={261, 219, 865, 949, 624, 390, 64, 67, 945, 253};
int aHash[10]={373, 983, 383, 211, 474, 944, 554, 271, 69, 747};
int bHash[10]={269, 229, 929, 949, 624, 390, 64, 67, 945, 253};
int m[10];

void seiveInitial(){
	int *x = (int*)malloc(sizeof(int)*VECT_SIZE),i,j;
	for (i=0;i<VECT_SIZE;i++) x[i] = i+1;
	for (i=1;i<VECT_SIZE;i++){
		j = x[i];
		if (j==-1) continue;
		int q = 1;
		while ((i + q*j) < VECT_SIZE){
			x[i+q*j] = -1;
			q++;
		}
	}
	i = VECT_SIZE;
	for (j=0;j<10;j++){
		for (i=i-1;i>0;i--){
			if (x[i]!=-1){
				m[j] = x[i];
				break;
			}
		}
	}
}

void setSignature(struct bloom *a){
	int ngroups = sizeof(long)*8*4;
	int nVals = (VECT_SIZE + NUM_BITS - 1)/NUM_BITS;

	int groupSize = (nVals + ngroups-1)/ngroups,i,j;
	a->signature = 0;
	long bitSetter = 1;
	for (i=0;i<ngroups;i++){
		for (j=i*groupSize;j<(i+1)*groupSize;j++){
			if ((j<=VECT_SIZE/NUM_BITS + 1)&&(a->bloom_vector[j]!=0))
				a->signature = a->signature | (bitSetter);
		}
		bitSetter = bitSetter<<1;
	}

}

void printSignature(struct bloom *a){
	long x = 1;
	int i;
	for (i=0;i<sizeof(long)*8;i++){
		if (a->signature & x)	printf("%d",1);
		else			printf("%d",0);
		x= x<<1;
	}
	printf("\n");
}

void setSeeds(){	//Assuming that the random number generated has already been seeded
	int i;
	for (i=0;i<10;i++) seed[i] = rand();
}

void init(struct bloom *bl)//Initialise the bloom filter
{
	int i=0;
	for(i=0;i<(VECT_SIZE)/NUM_BITS + 1;i++)
	{
		bl->bloom_vector[i]=0;
	}
	bl->nOnes = 0;
}

long num_ones(struct bloom *a)
{
	int i = 0, n = 0;
	for (i=0;i<VECT_SIZE/NUM_BITS + 1;i++) n += count_ones(a->bloom_vector[i]);
	return n;
}

int count_ones(int x){
	int count = 0;
	while (x!=0){
		x = x & (x-1);
		count++;
	}
	return count;
}

void insert(long val, struct bloom* bl)// INsert the element
{
	int i;
	for (i=0;i<K;i++){
		int a = hash(val,i);
		int loc = a/NUM_BITS;
		int off = a%NUM_BITS;
		int new = 1<<off;

		int prevValue = bl->bloom_vector[loc];
		bl->bloom_vector[loc] = bl->bloom_vector[loc] | new;
		if (prevValue != bl->bloom_vector[loc])
			bl->nOnes += 1;
	}
}

int is_in(long val, struct bloom* bl)
{
	int i;
	for (i=0;i<K;i++){
		int a = hash(val,i);
		int loc = a/NUM_BITS;
		int off = a%NUM_BITS;
		int new = 1<<off;
		if (!(bl->bloom_vector[loc] & new)) return 0;
	}
	return 1;
}


int is_inHash1(int val, struct bloom *bl){
	int a = hash1(val);
	int loc = a/NUM_BITS;
	int off = a%NUM_BITS;
	if (bl->bloom_vector[loc] & (1<<off)) 	return 1;
	else					return 0;
}

int is_inHash2(int val, struct bloom *bl){
	int a = hash2(val);
	int loc = a/NUM_BITS;
	int off = a%NUM_BITS;
	if (bl->bloom_vector[loc] & (1<<off)) 	return 1;
	else					return 0;
}

int is_inHash3(int val, struct bloom *bl){
	int a = hash3(val);
	int loc = a/NUM_BITS;
	int off = a%NUM_BITS;
	if (bl->bloom_vector[loc] & (1<<off)) 	return 1;
	else					return 0;
}

int dot_product(struct bloom *a, struct bloom *b)
{
	int i,num=0;
	for(i=0;i< (VECT_SIZE/NUM_BITS) + 1;i++)
	{
		int j=0;
		int n1=a->bloom_vector[i];
		int n2 = b->bloom_vector[i];
		int iter=1;
		while(j<NUM_BITS)
		{
			if( (n1 & iter) & (n2 & iter) )
				num++;
			iter=iter<<1;
			j++;
		}
	}
	return num;
}

int num_zero(struct bloom *a)
{
	int i,num=0;
	for(i=0;i< (VECT_SIZE/NUM_BITS) + 1;i++)
	{
		int j=0;
		int n1=a->bloom_vector[i];
		int iter=1;
		while(j<NUM_BITS)
		{
			if( (n1 & iter) )
				num++;
			iter=iter<<1;
			j++;
		}
	}
	return VECT_SIZE - num;
}

struct bloom *create_bloom(){
	struct bloom *a = malloc(sizeof(struct bloom));
	a->bloom_vector = (int*)malloc(VECT_SIZE/NUM_BITS + 1);
}

//N: This should be a bitwise AND rather than just min?
struct bloom * intersect_bloom(struct bloom *a, struct bloom *b)
{
	struct bloom *c = (struct bloom*)malloc(sizeof(struct bloom));
	init(c);
	int i,num=0;
	for(i=0;i< (VECT_SIZE/NUM_BITS) + 1;i++)
	{
		/*if(a->bloom_vector[i]<b->bloom_vector[i])
			a->bloom_vector[i] = a->bloom_vector[i];
		else
			a->bloom_vector[i] = b->bloom_vector[i];*/
		c->bloom_vector[i] = a->bloom_vector[i] & b->bloom_vector[i];
	}
	return c;
}

int numOnes_union(struct bloom *a, struct bloom *b){
	int i, num = 0;
	for (i=0;i<(VECT_SIZE/NUM_BITS + 1);i++){
		num += count_ones(a->bloom_vector[i] | b->bloom_vector[i]);
	}
	return num;
}

double elemsEst(struct bloom *a){
	double e = (-1.0*VECT_SIZE)/K;
	int n = num_ones(a);
	//printf("n = %d, e = %lf",n,e);
	e = e*log(1 - ((double)n)/VECT_SIZE);
	return e;
}

int numOnes_intersection(struct bloom *a, struct bloom *b){
	int i, num = 0, numCounts = 0;
/*
	int ngroups = sizeof(long)*8;
        int nVals = (VECT_SIZE + NUM_BITS - 1)/NUM_BITS;

        int groupSize = (nVals + ngroups-1)/ngroups,j;

	long andSig = a->signature & b->signature;
	long bitChecker = 1;
	for (i=0;i<ngroups;i++){
		if (andSig & bitChecker){
			for (j=i*groupSize;j<(i+1)*groupSize;j++){
				if (j<(VECT_SIZE/NUM_BITS + 1))
					num += count_ones(a->bloom_vector[j] & b->bloom_vector[j]);
				numCounts++;
			}
		}
		bitChecker = bitChecker<<1;
	}

*/
	for (i=0;i<(VECT_SIZE/NUM_BITS + 1);i++){
		int temp = a->bloom_vector[i] & b->bloom_vector[i];
		num += count_ones(temp);
	}


	return num;
}


double elemsIntersection(struct bloom *a, struct bloom *b, int ta, int tb){

	double nA = elemsEst(a);
	double nB = elemsEst(b);

	double nAUB = ((-1.0*VECT_SIZE)/K)* log(1.0 - ((double) numOnes_union(a,b))/VECT_SIZE);
	//printf("elemsIntersection: nA = %lf, nB = %lf, nAUB = %lf\n",nA,nB,nAUB);
	double val = nA + nB - nAUB;
	if (val<0) return 0;
	return (nA + nB - nAUB);

/*	double t1 = ta;
	double t2 = tb;
	double th = numOnes_intersection(a,b);
	double m = VECT_SIZE;
	if ((m - t1 - t2 + th)==0) //return VECT_SIZE;
		return 0;
	double temp = (th*m - t1*t2)/(m - t1 -t2 + th);
	double numerator = log(m - temp) - log(m);
	double denominator = K*log(1 - 1.0/m);
	return numerator/denominator;
*/
}

double elemsIntersection1(struct bloom *a, struct bloom *b){
	double m = VECT_SIZE;
	double z1 = m - num_ones(a);
	double z2 = m - num_ones(b);
	double z12 = m - numOnes_intersection(a,b);
	double numerator = log(z1) + log(z2) - log(z1 + z2 - z12);
	double denominator = K * (log(1.0/m) - log(1.0 - 1.0/m));
	return numerator/denominator;
}

struct bloom* union_bloom(struct bloom *a, struct bloom *b)
{
	int i,num=0;
	for(i=0;i< (VECT_SIZE/NUM_BITS) + 1;i++)
	{
		a->bloom_vector[i] = a->bloom_vector[i] | b->bloom_vector[i];
	}
	return a;
}

long num_elem(struct bloom *bl)
{
	double a = ((double)(VECT_SIZE - num_zero(bl))*1.0)/VECT_SIZE;
	a = 1-a;
	a= log(a);
	a = -VECT_SIZE*a*1.0/K;
	return a;
}
//Different hash functions

uint32_t murmur3_32(const char *key, uint32_t len, uint32_t seed) {
	static const uint32_t c1 = 0xcc9e2d51;
	static const uint32_t c2 = 0x1b873593;
	static const uint32_t r1 = 15;
	static const uint32_t r2 = 13;
	static const uint32_t m = 5;
	static const uint32_t n = 0xe6546b64;

	uint32_t hash = seed;

	const int nblocks = len / 4;
	const uint32_t *blocks = (const uint32_t *) key;
	int i;
	uint32_t k;
	for (i = 0; i < nblocks; i++) {
		k = blocks[i];
		k *= c1;
		k = ROT32(k, r1);
		k *= c2;

		hash ^= k;
		hash = ROT32(hash, r2) * m + n;
	}

	const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
	uint32_t k1 = 0;

	switch (len & 3) {
	case 3:
		k1 ^= tail[2] << 16;
	case 2:
		k1 ^= tail[1] << 8;
	case 1:
		k1 ^= tail[0];

		k1 *= c1;
		k1 = ROT32(k1, r1);
		k1 *= c2;
		hash ^= k1;
	}

	hash ^= len;
	hash ^= (hash >> 16);
	hash *= 0x85ebca6b;
	hash ^= (hash >> 13);
	hash *= 0xc2b2ae35;
	hash ^= (hash >> 16);

	return hash;
}

int md5_getHash(int key, int salt, int mi){
//	printf("HELLO I AM IN HERE!\n");
//	fflush(stdout);
	BYTE *d = (BYTE *)(&key);
	MD5_CTX ctx;
	BYTE buf[16], res[16];

	//md5_init(&ctx);
	//md5_update(&ctx,d,strlen(d));
	//md5_final(&ctx,buf);

	BYTE *temp = strcat(d,(char*)(&salt));

	md5_init(&ctx);
	md5_update(&ctx,temp,strlen(temp));
	md5_final(&ctx,res);

	BYTE temp_res[4];
	int i;
	for (i=0;i<4;i++)
		temp_res[i] = res[12+i];

	unsigned int a = *((unsigned int*) temp_res);
	unsigned int z = a%mi;
//	printf("BYE, I AM LEAVING %d\n",z);
//	fflush(stdout);
	return z;
}

int hash(long a, int i){
	if (HASHTYPE == 0){
		long z = aHash[i];
		z = z*a + bHash[i];
		return (z% m[i]);
	}
	else if (HASHTYPE == 1)
		return (murmur3_32((const char *)(&a),sizeof(int),seed[i])%m[i]);
	else
		return md5_getHash(a,seed[i],m[i]);
}

int hash1(int a)
{
	long z = a1;
	z = z*a + b1;
	return (z)%VECT_SIZE;
	//return a%VECT_SIZE;
	//return (murmur3_32((const char *)(&a),sizeof(int),seed[0])%VECT_SIZE);
	//return md5_getHash(a,seed1);
}

int hash2(int a)
{
	long z = a2;
	z = z*a + b2;
	return (z)%VECT_SIZE;
	//return (a+a%D2+A2)%VECT_SIZE;
	//return (murmur3_32((const char *)(&a), sizeof(int), seed[1]) % VECT_SIZE);
	//return md5_getHash(a,seed2);
}

int hash3(int a)
{
	long z = a3;
	z = z*a + b3;
	return (z)%VECT_SIZE;
	//return (a+a%D3+A3)%VECT_SIZE;
	//return (murmur3_32((const char *)(&a), sizeof(int), seed[2]) % VECT_SIZE);
	//return md5_getHash(a,seed3);
}

int getInverse(int b, int m){
	int t = 0, newt = 1, r = m, newr = b;
	while (newr!=0){
		int quotient = r/newr;
		int temp = newt;
		newt = t - quotient*newt;
		t = temp;
		temp = newr;
		newr = r - quotient*newr;
		r = temp;
	}
	if (t<0) t = t+m;
	if (r>1){
		printf("ERROR IN MODULO INVERSE: %d is not invertible\n",b);
		exit(1);
	}
	return t;
}

int chineseRemainder(int *a, int *m){
	int n = 3;
	int b[3] = {0,0,0};
	int M = 1,i;
	for (i=0;i<n;i++) M = M*m[i];
	int x = 0;
	for (i=0;i<n;i++){
		b[i] = getInverse(M/m[i],m[i]);
		x = x + (a[i]*b[i]*M)/m[i];
	}
	return x;
}

int solveVal(int v1, int v2, int v3){
	int r1 = v1;
	int z = (v1 + A2)%VECT_SIZE;
	int p = (D2 - v2 + z)/VECT_SIZE;
	if ((D2-v2+z)%VECT_SIZE == 0) p--;
	int r2 = VECT_SIZE*p + v2 - z;

	z = (v1 + A3)%VECT_SIZE;
	p = (D3 - v3 + z)/VECT_SIZE;
	if ((D3 - v3 + z)%VECT_SIZE == 0) p--;
	int r3 = VECT_SIZE*p + v3 - z;
	int b[3] = {r1,r2,r3};
	int m[3] = {VECT_SIZE,D2,D3};
	return chineseRemainder(b,m);
}

int *reservoirSample(int k, struct bloom *a){
	//srand(time(NULL));
	long n = num_ones(a);
	if (n<k) return NULL;	//Nothing can be done

	int *indices = (int *)malloc(n*sizeof(int));
	int i,z=0;
	for (i=0;i<VECT_SIZE;i++){
		int v = a->bloom_vector[i/NUM_BITS] & (1 << (i%NUM_BITS));
		if (v>0){
			indices[z] = i;
			z++;
			if (z==n) break;
		}
	}
	//reservoir sampling
	int *s = (int *)malloc(k*sizeof(int));
	for (i=0;i<k;i++)
		s[i] = indices[i];
	for (i=k;i<n;i++){
		int j = rand()%(i+1);
		if (j<k) s[j] = indices[i];
	}
	free(indices);
	return s;
}

int sampleFromBloom(struct bloom *a, int nVertices){
/*
	int currSample = -1, nT;
	for (nT=1;nT<=100;nT++){
		int *s = reservoirSample(3,a);
		printf("\t%d,%d,%d\n",s[0],s[1],s[2]);
		int l[6];
		l[0] = solveVal(s[1],s[2],s[3]);
		l[1] = solveVal(s[1],s[3],s[2]);
		l[2] = solveVal(s[2],s[1],s[3]);
		l[3] = solveVal(s[2],s[3],s[1]);
		l[4] = solveVal(s[3],s[1],s[2]);
		l[5] = solveVal(s[3],s[2],s[1]);
		printf("\t%d,%d,%d,%d,%d,%d\n",l[0],l[1],l[2],l[3],l[4],l[5]);
		//srand(time(NULL));
		int i = rand()%6;
		free(s);
		if (rand()%nT == 0)
			currSample = l[i];
	}
	return currSample;
*/
	int *z = reservoirSample(1,a);
	if (z==NULL) return -1;	//empty bloom filter
	int j = z[0],numCandidates = 0,i,count = 0,currSample = -1,sample=0;

	while (sample<nVertices){
		sample = j - b1 + count*VECT_SIZE;
		if (sample%a1==0){
			int s = sample/a1;
			//if ((s>=nVertices)||(s<0)) break;
			//fflush(stdout);
			if (is_in(s,a)){
				numCandidates++;
				i = rand()%numCandidates;
				if (i==0) currSample = s;
			}
		}
		//printf("sample = %ld,b1=%d,j=%d,count=%d\n",sample,b1,j,count);
		count++;
		sample = sample/a1;
	}
	count=0;
	sample = 0;
	while (sample<nVertices){
                sample = j - b2 + count*VECT_SIZE;
                if (sample%a2==0){
                        int s = sample/a2;
			//if ((s>=nVertices)||(s<0)) break;
                        if (is_in(s,a)){
                                numCandidates++;
                                i = rand()%numCandidates;
                                if (i==0) currSample = s;
                        }
                }
                sample = sample/a2;
                count++;
        }

	count = 0;
	sample = 0;
	while (sample<nVertices){
                sample = j - b3 + count*VECT_SIZE;
                if (sample%a3==0){
                        int s = sample/a3;
			//if ((s>=nVertices)||(s<0)) break;
                        if (is_in(s,a)){
                                numCandidates++;
                                i = rand()%numCandidates;
                                if (i==0) currSample = s;
                        }
                }
                sample = sample/a3;
                count++;
        }
	//if (numCandidates==0) printf("j = %d\n",j);
	return currSample;
}

int sampleFromBloom1(struct bloom *a, int nVertices){
	int *z = reservoirSample(1,a);
	if (z==NULL) return -1;
	int j = z[0];
	int p = 0, nL = 0,currS = -1,l,y;
	//h1(x) set this bit

	while ((p*VECT_SIZE+j)<nVertices){
		int k = p*VECT_SIZE+j;
		if (is_in(k,a)&&(k<=nVertices)){
			nL++;
			int r = rand()%nL;
			if (r==0) currS = k;
		}
		p = p+1;
	}

	//h2(x) set this bit
	for (l=j;l<(nVertices+D2+A2);l=l+VECT_SIZE){
		int q = 0;
		int l1 = l-A2;
		while ((l1-D2*q)/2 > D2) q++;
		while ((l1-D2*q)>=0){
			if ((l1-D2*q)%2==0){
				int as = (l1-D2*q)/2;
				int y = D2*q + as;
				if (is_in(y,a)&&(y<=nVertices)){
					nL++;
					int r = rand()%nL;
					if (r==0) currS = y;
				}
			}
			q++;
		}
	}

	//h3(x) set this bit
	for (l=j;l<(nVertices+D3+A3);l=l+VECT_SIZE){
                /*for (y=(l-A3)/2;y<(l-A3);y=y+D3/2){
                        if (is_in(y,a)){
                                nL++;
                                int r = rand()%nL;
                                if (r==0) currS = y;
                        }
                }*/
		int q = 0;
                int l1 = l-A3;
                while ((l1-D3*q)/2 > D3) q++;
                while ((l1-D3*q)>=0){
                        if ((l1-D3*q)%2==0){
                                int as = (l1-D3*q)/2;
                                int y = D3*q + as;
                                if (is_in(y,a)&&(y<=nVertices)){
                                        nL++;
                                        int r = rand()%nL;
                                        if (r==0) currS = y;
                                }
                        }
                        q++;
                }

        }
	/*if (currS==-1){
		printf("j = %d, nL = %d\n",j,nL);
	}*/
	return currS;
}

int belongsInArr(int val, int *arr, int size){
        int j;
        for (j=0;j<size;j++)
                if (arr[j]==val)
                        return 1;
        return 0;
}

int mainS1(int argc, char *argv[]){
	srand(time(NULL));
	clock_t begin,end;
	double time_spent = 0;
        int nVertices = atoi(argv[1]), nNodes = atoi(argv[2]), nExp = 1000;
        struct bloom *k = (struct bloom *)malloc(sizeof(struct bloom));
        int i,j,nSuccess=0,nNull=0,nRange=0,nMystery=0;
        int *realArr = (int*)malloc(sizeof(int)*nNodes);
        for (i=0;i<nExp;i++){
		printf("*********************  i = %d ****************************\n",i);
                init(k);
                for (j=0;j<nNodes;j++){
                        int val = rand()%nVertices;
                        realArr[j] = val;
                        insert(val,k);
                }
		begin = clock();
                int sample = sampleFromBloom(k,nVertices);
		end = clock();
		time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
		/*int d = VECT_SIZE*D2*D3;
		printf("SAMPLED HERE = %d\n",sample);
		if (sample<0)
			while (sample<0) sample = sample + d;
		if (sample>nVertices)
			while (sample>nVertices) sample = sample - d;*/

		//printf("I SAMPLED %d\n",sample);
                if (belongsInArr(sample,realArr,nNodes)) nSuccess++;
                else if (sample==-1) nNull++;
		if (!is_in(sample,k)) nMystery++;
		if (0<=sample && sample<=nVertices) nRange++;
        }
        double acc = ((double)nSuccess)/nExp;
        free(realArr);
        printf("ACCURACY = %lf\n",acc);
        printf("SUCC = %d, NULL = %d, NRANGE = %d, NMYSTERY = %d\n",nSuccess,nNull,nRange,nMystery);
	printf("TOTAL TIME SPENT = %lf\n",time_spent);
}

int mainR1(int argc, char *argv[]){
        srand(time(NULL));
        int nVertices = atoi(argv[1]), nNodes = atoi(argv[2]);
        struct bloom *k = (struct bloom *)malloc(sizeof(struct bloom));
        int *realArr = (int*)malloc(sizeof(int)*nNodes);

        FILE *fi = fopen(argv[3],"r");
        int i = 0;
        while (i<nNodes){
                int val=0;
                fscanf(fi,"%d\n",&val);
                insert(val,k);
                realArr[i] = val;
                i++;
        }
        fclose(fi);

        FILE *fo = fopen("RepeatedSamples","w");
        for (i=0;i<1000000;i++){
                int sample = sampleFromBloom1(k,nVertices);
                if (belongsInArr(sample,realArr,nNodes))
                        fprintf(fo,"1,%d\n",sample);
		else	fprintf(fo,"2,%d\n",sample);
        }
        fclose(fo);

}
