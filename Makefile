BSTSample: malloc_count-0.7.1/malloc_count.o
	gcc BSTtesting.c bloom.c malloc_count-0.7.1/malloc_count.o -lm -lrt -ldl -o BSTSample
