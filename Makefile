BSTSample: malloc_count-0.7/malloc_count.o
	gcc BSTtesting.c bloom.c malloc_count-0.7/malloc_count.o -lm -lrt -ldl -o BSTSample
BSTReconstruct: malloc_count-0.7/malloc_count.o
	gcc bloomTreeReconstructE.c bloom.c malloc_count-0.7/malloc_count.o -lm -lrt -ldl -lpthread -o BSTReconstruct
dynamicBTtwitterRatio:
	gcc dynamicBTtwitterRatio1.c bloom.c -lm -lrt -o dynamicBTtwitterRatio
