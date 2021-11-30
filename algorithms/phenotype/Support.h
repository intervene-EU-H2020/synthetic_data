#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// DOES NOT deal with nan in Geno input! need to encode with something before used as input

// Defined parameters
#define nMaxCausal 500000
#define nMaxPop 10
#define nMaxInd 100000

gsl_rng * r;

// IO file names
char InCausal[1000];
char InRef[1000];
char InSample[1000];
char InGeno[1000];
char OutCausal[1000];
char OutPheno[1000];

// Mixture model parameters
double pCausal;
int PolyFlag, CausalFlag;
char tmpHerr[10000];
char tmpPrev[10000];
char tmpProb[10000];
char tmpWeight[10000];
char tmpCov[10000];
int nComp;
double * ProbComp;
double * wComp;

// GenoMat: nPop x nIndPerPop 
double GenoMat[nMaxPop][nMaxInd];
int PopIndicator[nMaxInd];
int nSamplePerPop[nMaxPop];
char PopList[nMaxPop][50];
char SampleList[nMaxInd][100];

char CausalList[nMaxCausal][50];

double CausalMAF[nMaxPop];
double CausalLDscore;
double CausalAnnot;
double CausalBeta[nMaxPop];

double * GenoEff[nMaxPop];
double * EnvEff[nMaxPop];
double * PhenoSim[nMaxPop];

double herr[nMaxPop];
double prev[nMaxPop];
double VarGeno[nMaxPop];
double VarEnv[nMaxPop];
long int nCausal, nSample, nPop;

gsl_matrix * Sigma;
gsl_matrix * L;
gsl_vector * mu;

double a, b, c, det, prob;

struct SNPinRef {
   char SNP[50];
   double AfricaMAF;
   double AfricaLDscore;
   int exon; // 1 -- yes, 0 -- no
   int DHS; // 1 -- yes, 0 -- no
};  
typedef struct SNPinRef SNPinfo;

SNPinfo RefSNP[26][200000];
long int SNPct[26];
char buffer[50000000];

void ReadParam(const char *ParIn);

void ExtractParam();

int PopIndex(char PopCode[50], int nPopCt);

void ReadPopulation();

void ReadCausal();

void ReadRef();

double IsCausal(char SNP[50]);

long int FindSNPinRef(char SNPname[50], int chr);

double GetMAF(int PopIndex);

double GetBeta(double BaseBeta, int PopIndex);

void AnalyzeSNP(int chr, char SNP[50]);
