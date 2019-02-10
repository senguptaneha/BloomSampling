# BloomSampling
##This repo is under construction

## Introduction 
Bloom filters are summary data structures that store the elements of a set compactly in a bit array. 
Given a Bloom fitler B(S) that stores the set S, how can we sample elements from the original set S?
A related question is how to reconstruct the set S from the Bloom filter?

This repository implements algorithms to achieve both of the above tasks [1].

## Installation

To compile the code, download the malloc_count library from here: https://panthema.net/2013/malloc_count/
and place inside the repo directory.
Compile using 
```
make
```
## Sampling

To run the sampling tests, use
```
./BSTtesting namespaceSize nSets inputFile numHashFunctions threshold desiredAccuracy
```
### Parameters
* `namespaceSize` is the size of the universe that S is drawn from. If namespace size is `M`, the universe is assumed to be `[0 ... M-1]`.
* `nSets` is the number of such sets S drawn from the universe. The code samples 100 times from each set
* `inputFile` contains nSets lines. Each line is a comma separated list of values in the set.
* `numHashFunctions` is typically set to 3. It is the number of hash functions that a Bloom filter uses. For the same Bloom filter size, more hash functions imply greater accuracy.
* `threshold` is the intersection size threshold used by the algorithm. A smaller threshold will cause higher sampling time but fewer errors in sampling (i.e. it will avoid that valid elements in the Bloom filter are never sampled). (See our paper [1] for details.). For sampling, this threshold will typically be set to 0.
* `desiredAccuracy` is the desired accuracy of sampling. For example, if this parameter is 0.8, then the drawn sample belongs to the original set S with probability 80 %.

### Output
The code reports the average time taken to sample using the BloomSampleTree based algorithm ([1]), and the dictionary attack which is the baseline. The sampling time reported is averaged over 100 sampling rounds over all query sets. In addition, it outputs the measured memory footprint of the method.


## Datasets
Simulated datasets available at:  https://drive.google.com/file/d/1tOiwwnvl5uoFD0Zaaq-t8DlSJnycU8aE/view?usp=sharing
Twitter data available at https://drive.google.com/file/d/1YKng-ueZkj0jbWiNkjpWVwn72cGhDk9w/view?usp=sharing

## References
[1] Neha  Sengupta,  Amitabha  Bagchi,  Srikanta  Bedathur,  and  Maya  Ramanath. Sampling  and  Reconstruction  Using  Bloom  Filters. IEEE Transactions on Knowledge and Data Engineering, 2017.