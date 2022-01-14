# Phenotype generation algorithm

# !There are some changes in parameters. I will update this readme on Monday, hopefully #

_"Working code" for phenotype generation, subject to further updates._

Dependents: gsl, blas

To compile, run 
```
gcc Main.c Support.c -o where/and/what/you/want/it/to/be -L. -lm -lgsl -fPIC -lcblas -lblas
```

To use, run
```
.\TheTool ParFile seed
```
where ```Parfile``` is a parameter file in plain text that looks like below: 
```
Heritability 0.5,0.6
Prevalence 0.1,0.2
a -0.4 
b -0.5
c 0.5
nComponent 3
CompWeight 1,2,4
PopCovar 1,1,1,1 
Polygenicity 0.05 
CausalList /path/to/the/causal/SNP/list
SampleList /path/to/the/sample/file
Reference /path/to/the/reference/file
GenoFile /path/to/the/input/genotype/file
Output /path/to/the/output/file
```
and ```seed``` is a numeric seed for randomness.

_This will be updated to connect with config_

Currently the parameters are as below

```Heritability``` is observed causal SNP heritability in each population, nPOP entries separated by comma.

```Prevalence``` is disease prevalence in each population, nPOP entries separated by comma.

```a```, ```b```, ```c``` have -0.4, -0.5 and 0.5 as their suggested values.

```nComponent``` and ```CompWeight``` are number of guassian mixture components and their weights.

```PopCovar``` is flattened correlation matrix for population genetic correlation, nPOP x nPOP entries separated by comma.

```Polygenicity``` set at 0.05 means around 5% SNPs will be causal.

```CausalList``` is a list of predefined SNPs to be used as causal, overrides ```Polygenicity``` parameter if specified.

```SampleList``` is a population code list, in the order of sample listed in the header of input .traw file.

```Reference``` **It also takes a reference file for LD, but it's a little too big for github (~36Mb in .gz). Need to put it somewhere**

```GenoFile``` Input genotype file, in .traw format. Can be generated using ```plink --recode A-transpose``` from other formats.

It output two files, ```.pheno``` includes the genetic effect, environmental effect and synthetic phenotype (in liability) for each individual.

```.causal``` includes the causal SNPs and their effect sizes in each population.



