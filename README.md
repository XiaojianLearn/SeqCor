# SeqCor

## Introduction

The sequence of gRNAs in the library has impacts on the investigation of CRISPR/Cas9-based screening, and building up a standard strategy for correcting results of all gRNA libraries is impracticable.

SeqCor is a program based on random forest algorithm, which enables researchers to address the result bias potentially caused by the composition of gRNA sequence by the organization of gRNA in the library
utilized in CRISPR/Cas9-based screening. In light of four public CRISPR/Cas9-based screening data, SeqCor
completely computerizes the extracting of sub-sequence that may influence sgRNA knockout efficacy in a machine learning manner.

## Compile

- Linux: g++

```
g++ SeqCor.cpp -o SeqCor -O3 -Wall -fopenmp -msse2
```

## Usage

### Optional parameters

```
chmod +x SeqCor
./SeqCor
```

Usage:	grna [options] grna input output

- -g:	covariant	default=NULL
- -d:	2^dims dimensions	default=9
- -t:	trees	default=896
- -m:	mtry	default=2
- -l:	leaf	default=cbrt(N)

### Example

```
SeqCor gRNA.txt input.txt output.txt

gRNA.txt	input.txt	73151
signal is NOT adjusted by gene
dims	512
trees	896
mtry	2
leaf	11
```

## Citation

JianXiao Liu, Qiurong Ding, Yi Wang et al. SeqCor: correct the effect of gRNA sequences on experiments by machine learning algorithm. 2020 (**Manuscript submitted**)

## Contacts

For questions, request and bug reports, please contact: [wangyi_fudan@fudan.edu.cn](mailto:wangyi_fudan@fudan.edu.cn) 

