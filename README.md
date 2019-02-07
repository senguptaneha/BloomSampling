# BloomSampling
##This repo is under construction

## Introduction 
Bloom filters are summary data structures that store the elements of a set compactly in a bit array. 
Given a Bloom fitler B(S) that stores the set S, how can we sample elements from the original set S?
A related question is how to reconstruct the set S from the Bloom filter?

This repository implements algorithms to achieve both of the above tasks.

## Installation

To compile the code, download the malloc_count library from here: https://panthema.net/2013/malloc_count/
and place inside the repo directory.
Compile using 
```
make
```
## Sampling

   
Simulated datasets available at:  https://drive.google.com/file/d/1tOiwwnvl5uoFD0Zaaq-t8DlSJnycU8aE/view?usp=sharing
Twitter data available at https://drive.google.com/file/d/1YKng-ueZkj0jbWiNkjpWVwn72cGhDk9w/view?usp=sharing
