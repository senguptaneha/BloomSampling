# BloomSampling
### This repository is under construction

## Introduction 
Bloom filters are summary data structures that store the elements of a set compactly in a bit array. 
Given a Bloom fitler BF(S) that stores the set S, how can we sample elements from the original set S?
A related question is how to reconstruct the set S from the Bloom filter?

This repository implements algorithms to achieve both of the above tasks [1].

## Installation

To compile the code, download the malloc_count library from here: https://panthema.net/2013/malloc_count/
and place inside the repo directory.

For sampling, compile using 
```
make BSTSample
```

For sampling with sparse namespaces, compile using
```
make dynamicBTtwitterRatio
```

For reconstruction, compile using
```
make BSTReconstruct
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

### Sampling with a sparsely occupied namespaces
In many real world applications, the namespace of elements from which the set S is drawn may be sparsely occupied. This means, that a majority of elements in `[0 ... M-1]` are not a part of any set. An example of this use-case is that of Twitter data, where a set corresponds to a hash tag and contains the user ids that tweeted with that hash tag. Since user ids are large alpha-numeric strings, the namespace of user ids is sparsely occupied.

In the scenario where the namespace is known to be sparsely occupied, the sampling process can be made faster and with smaller memory footprint. For such cases, use,
```
./dynamicBTtwiiterRatio namespaceSize nSets inputFile numHashFunctions BloomFilterSize threshold numLevels
```

The parameters are similar to that of sampling above. Additionally, `BloomFilterSize` should provide the size of the Bloom filters used in bits. For larger values of `namespaceSize` or larger sets to be sampled from, this size should be higher. `numLevels` gives the number of levels in `BloomSampleTree` [1].

## Reconstruction

To run the sampling tests, use
```
./BSTReconstruct namespaceSize inputFile numHashFunctions threshold desiredPrecision
```
Unlike `BSTSample`, the provided `BSTReconstruct` utility can be used to reconstruct only a given set at a time. However, it can also be easily extended to reconstruct multiple sets at once.

### Parameters
* `namespaceSize` is the size of the universe that S is drawn from. If namespace size is `M`, the universe is assumed to be `[0 ... M-1]`.
* `inputFile` contains elements of the set, one element per line.
* `numHashFunctions` is the number of hash functions that a Bloom filter uses.
* `threshold` is the intersection size threshold used by the algorithm. This is an important parameter for the reconstruction algorithm. A higher threshold will result in higher precision, smaller reconstruction time, but smaller recall. A small threshold will give lower precision, larger reconstruction time, and higher recall. Best values for this threshold for typical parameter settings can be found in the file `bestOverlapThresholds`.
* `desiredPrecision` is the desired precision of reconstructed set.

### Output
The code reports the size of the reconstructed set, the precision, reconstruction time in seconds, and the number of intersection and membership operations.

## Datasets
Simulated datasets available at:  https://drive.google.com/file/d/1tOiwwnvl5uoFD0Zaaq-t8DlSJnycU8aE/view?usp=sharing
Twitter data available at https://drive.google.com/open?id=1yIra0I7qQqNyvNdVshw9EugYM0WB0APW

## References
[1] Neha  Sengupta,  Amitabha  Bagchi,  Srikanta  Bedathur,  and  Maya  Ramanath. Sampling  and  Reconstruction  Using  Bloom  Filters. IEEE Transactions on Knowledge and Data Engineering, 2017.