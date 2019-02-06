./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 123174 0 7 0.5 >> hashFamilyMurmur
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 143343 0 7 0.6 >> hashFamilyMurmur
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 168630 0 7 0.7 >> hashFamilyMurmur
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 204937 0 7 0.8 >> hashFamilyMurmur
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 273404 0 6 0.9 >> hashFamilyMurmur
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 626895 0 5 1.0 >> hashFamilyMurmur


sed 's/HASHTYPE 1/HASHTYPE 2/' <../bloom.c > ../bloom.c.temp
mv ../bloom.c.temp ../bloom.c

gcc BSTtesting.c ../bloom.c ../malloc_count-0.7/malloc_count.o -lm -lrt -ldl -o BSTtesting

./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 123174 0 11 0.5 >> hashFamilyMD5
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 143343 0 11 0.6 >> hashFamilyMD5
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 168630 0 11 0.7 >> hashFamilyMD5
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 204937 0 11 0.8 >> hashFamilyMD5
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 273404 0 10 0.9 >> hashFamilyMD5
./BSTtesting 10000000 100 ../SimulatedDatasets/10M/Uniform10K_100 3 626895 0 9 1.0 >> hashFamilyMD5
