# Phenotype generation algorithm

_"Working code" for phenotype generation, subject to further updates._

_2023.05.10 Implemented interaction effect; user specified causal effect size; user specified homo/heter effect; made annotation file optional; fixed a bug searching for SNPs in the reference file; and included some minor fix to (hopefully) make more robust. Tested and working. Please let us know if you notice anything abnormal with the new functions._

_2022.07.28 Now accepting plink binary files as input. Need to install plinkio library._

_2022.06.21 Now process genotype input file by chromosome, so that for larget data set generating traw file is more feasible._


_2022.01.18 Changed non positive definite matrices handling a bit, but the input (PopulationCorr, TraitCorr) are still tricky for users._
_Hmm.._

_2022.01.19 Made a little faster (hopefully?)._

Dependents: gsl, blas, lplinkio

## Basic usage

To compile, run 
```
gcc Main.c Support.c -o where/and/what/you/want/it/to/be -L. -lm -lgsl -fPIC -lcblas -lblas -lplinkio
```

To use, run
```
.\TheTool ParFile seed
```
where ```Parfile``` is a parameter file in plain text that looks like below: 
```
nPopulation 2
nTrait 2
a -0.4
b -1
c 0.5
nComponent 3
PropotionGeno 0.5,0.7,0.8,0.2
PropotionCovar 0.0,0.0,0.0,0.0
Polygenicity 0.1,0.08
Pleiotropy 1,0.5
TraitCorr 1,-0.7,-0.7,1
PopulationCorr 1,0.25,0.25,1
CompWeight 1,5,10
CausalList /path/to/the/causal/SNP/list
SampleList /path/to/the/sample/file
Reference /path/to/the/reference/file
GenoFile /path/to/the/input/genotype/file
Output /path/to/the/output/file
```
and ```seed``` is a numeric seed for randomness.

_This will be updated to connect with config_

Currently the parameters are as below

```nPopulation``` is the number of populations (nPop), integer.

```nTrait``` is the number of traits (nTrait), integer.

```PropotionGeno``` is observed causal SNP heritability in _each population, each trait_. Flatten nPop * nTrait matrix, entries separated by comma.

```PropotionCovar``` is observed proportion of variance contributed by the covariate (input in SampleList file) in _each population, each trait_. Flatten nPop * nTrait matrix, entries separated by comma.

```Prevalence``` is disease prevalence in _each population, each trait_. Flatten nPop * nTrait matrix, entries separated by comma. If prevalence is specified, output will include a column for binary case/control statues.

```a```, ```b```, ```c``` have -0.4, -1 and 0.5 as their suggested values.

```nComponent``` and ```CompWeight``` are number of guassian mixture components and their weights.

```PopulationCorr``` is a flattened correlation matrix for _population genetic correlation_ (symmetric positive definite). nPop x nPop entries separated by comma. Here _population genetic correlation_ is defined as correlation of SNP effects on the same trait across differnet populations. For the same trait, we assume the same set of causal variants shared across all populations, but each population can have their specific but overall correlated effect sizes.

```TraitCorr``` is a flattened correlation matrix for _traits correlation_ (symmetric positive definite). nTrait x nTrait entries separated by comma. Here _traits correlation_ parameter applies to 1. genetic effects of _shared_ causal variants between traits and 2. correlation of enviromental noise between traits.

```Pleiotropy``` is nTrait vector of trait's pleiotropy relationship **comparing to trait 1**. i.e. if trait 2 has Pleiotropy = 0.9, it means 90% of causal SNPs in trait 1 are also causal in trait 2. Therefore, first entry of ```Pleiotropy``` vector is always 1. Entries separated by comma.

```Polygenicity``` nTrait vector of trait polygenicity, measured by proportion of total SNPs being causal. e.g. Polygenicity = 0.05 meaning around 5% SNPs will be causal.

```CausalList``` is the prefix for lists of predefined SNPs to be used as causal, overrides ```Polygenicity``` parameter if specified. Each column contains causal SNPs for one trait, columns separated by comma. CausalList is specified per phenotype. **For each trait, the causal list file should be names as**  ```prefix```**n, where n is the trait index.**

```SampleList``` is a population and covariate (if any) list, in the order of sample listed in the header of input .traw file. Fist column should contain categorical population code for each sample (mandatory); rest of the columns should contains numerical covariates, columns separated by comma. Length of this file must match sample size in .traw file.

```Reference``` It also takes an **optional** reference file for LD. A example reference file with hapmap3 SNPs can be downloaded with HAPNEST. If this is not provided, LD and functional based effect size simulation will be disabled.

```GenoFile``` Input genotype file, in .traw format **by chromosome**, ie takes plink files by chromosome, named as _GenoFile_-chr.bim/bed/fam, where chr is the chromosome number (1-22). Can be generated using ```plink --make-bed``` from other formats

_For each trait_, it outputs two files, ```.pheno``` includes the genetic effect, environmental effect and synthetic phenotype for each individual.

```.causal``` includes the causal SNPs and their effect sizes in each population.

## User specified variant effects

As showed above, user can specified a list of causal variants for each trait, or trait polygenicity parameter to control the amount of variants being causal. In both cases, effect sizes of causal variants will be drawn randomly from a distribution. Alternatively, user can also specify the effect of each causal variant in the list. In this case, ```CausalList``` for each trait should contain 2 column: variant ID, and effect size in terms of beta for homo and homozygous genotypes, seperated by comma. Below is an example
```
rs56175625  0.5,0.25
rs57804877  0.2,0.15
rs201404203 0.3,0
```
The first row means causal variant ```rs56175625``` has effect 0.5 for genotype AA and effect 0.25 for Aa, were A is the risk allele. 
For user-sepcified effects to be applied, please make sure to add below in the parameter file:
```
CausalEffect  T
```
Otherwise the effect size column in the ```CausalList``` will be ignored.
**Note when user-sepcified effects are applied, trait and population genetic correlation will be disabled.** There will not be any population specific effect under user-sepcified effect option. 

## Adding interaction terms

You can now add interaction effects to the phenotype using the ```Interaction``` parameter. Similar to ```CausalList```, it takes prefix for lists of predefined interaction items to be applied on the phenotype. The list should contrain three columns: interaction term1, term2 and the effect size. An example looks like below:
```
rs202076079	covar1	0.1
covar1	covar2	0.3
rs16831418	rs115414386	0.2
```
Row 1 means that the genotype of variant rs202076079 x the first covariate in the ```SampleList``` (first column proceeding sample ID) has effect of 0.1 on the phenotype. Note the interaction can include SNP x SNP interaction (Epistasis effect), SNP x covariate effect or covariate x covariate effect. Covariate terms should be specified as ```covar```n, where n is the covariate index as ordered in the ```SampleList```.
**Same as ```CausalList```, for each trait, the interation item file should be names as**  ```prefix```**n, where n is the trait index.**



