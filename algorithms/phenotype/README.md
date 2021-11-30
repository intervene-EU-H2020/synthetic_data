# Phenotype generation algorithm

"Working code" for phenotype generation, subject to further updates.

Dependents: gsl, blas

To compile, run 
```
gcc Main.c Support.c -o where/and/what/you/want/it/to/be -L. -lm -lgsl -fPIC -lcblas -lblas
```

To use, run
```
.\TheTool ParFile seed
```
where ```Parfile``` is a parameter file in plain text, and ```seed``` is a numeric seed for randomness.

_This will be updated to connect with config_

**It also takes a reference file for LD, but it's a little too big for github (~36Mb in .gz)**



