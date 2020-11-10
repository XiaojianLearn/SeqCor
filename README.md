# SeqCor



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
#example data
```
$ unzip example\ data.zip

$ more example\ data/gRNA.txt
AAAGAGAGAAGAAATCTTAT
ACAACTGCCGATTATTTGTT
GCAGAAATGGGCGCCTCTTA
TGGTGGCCCTCCACCTGGTT
AAATATGGTGGCCCTCCACC

$ more example\ data/input.txt
5.3659760
5.6869754
6.0450053
5.6167711
0

$ ./SeqCor example\ data/gRNA.txt example\ data/input.txt output.txt

$ more output.txt
5.69994
5.66099
6.01038
5.68672
6.05828

```
## Contacts

For questions, request and bug reports, please contact: [wangyi_fudan@fudan.edu.cn](mailto:wangyi_fudan@fudan.edu.cn) 

