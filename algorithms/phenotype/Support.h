#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// DOES NOT deal with nan in Geno input! need to encode with something before used as input

int status;
int statusPop;
int statusTrait;

// Defined parameters
#define nMaxCausal 500000
#define nMaxPop 10
#define nMaxTrait 10
#define nMaxInd 100000
#define nPCAsnp 10000
#define nMaxCovar 10
#define nMaxBetaGen 10000

gsl_rng * r;

// IO file names
char InCausal[1000];
char InRef[1000];
char InSample[1000];
char InGeno[1000];
char OutCausal[1000];
char OutPheno[1000];

int PolyFlag, BinaryFlag;
double pCausal[nMaxTrait];
int CausalFlag[nMaxTrait];
char tmpPCausal[10000];
char tmpPleio[10000];
char tmpGenoEffProp[10000];
char tmpCovarEffProp[10000];
char tmpPopCorr[10000];
char tmpPrev[10000];
char tmpProb[10000];
char tmpWeight[10000];
char tmpTraitCorr[10000];
int nComp;
double * ProbComp;
double * wComp;

// GenoMat: nPop x nIndPerPop 
double PopMatTmp[nMaxPop][3][nMaxInd];
double GenoMat[nMaxInd];
double CovarMat[nMaxCovar][nMaxInd];
int PopIndicator[nMaxInd];
int nSamplePerPop[nMaxPop];
char PopList[nMaxPop][50];
char SampleList[nMaxInd][100];
char CausalList[nMaxTrait][nMaxCausal][50];

double CausalMAF[nMaxPop];
double CausalLDscore;
double CausalAnnot;
double CausalBeta[nMaxPop][nMaxTrait];
// double CovarBeta[nMaxPop][nMaxTrait];

// Include effects from shared causal variants (effect size generated by genetic correlation) and trait specific effects
double BaseBeta[nMaxPop][nMaxTrait]; // Population * Trait
double * GenoEff[nMaxTrait]; 
double * CovarEff[nMaxTrait]; // Covar Eff includes effect from PCs and non genetic covariates (user input)
double * EnvEff[nMaxTrait]; // Include corrlated noise and trait specific noise 
double * PhenoSim[nMaxTrait][2]; // 1--Continious measurement/liability 2--binary diagnosis

double GenoEffProp[nMaxPop][nMaxTrait];
double CovarEffProp[nMaxPop][nMaxTrait];
double EnvEffProp[nMaxPop][nMaxTrait];
double Prev[nMaxPop][nMaxTrait];
// Pleiotropy model with first trait as reference (how many causal SNPs shared with the first trait), starting with 1
double Pleio[nMaxTrait];
double TraitCorr[nMaxTrait][nMaxTrait]; // Only support positive genetic correlation now
double PopCorr[nMaxPop][nMaxPop];

gsl_matrix * GenoBeta;
gsl_matrix * SigmaTrait;
gsl_matrix * SigmaPop;
gsl_matrix * L;
gsl_vector * mu;

double VarGeno[nMaxPop][nMaxTrait];
double VarCovar[nMaxPop][nMaxTrait];
double VarEnv[nMaxPop][nMaxTrait];
double GCEweight[nMaxPop][3];
long int nCausal[nMaxTrait], nSample, nBetaIndex;
int nPop, nTrait, nItem, nValidItem, nCovar;
double a, b, c, prob;

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

void MakeCovMat();

int PopIndex(char PopCode[50], int nPopCt);

void ReadPopulation();

void ReadCausal();

void ReadRef();

void BaseBetaGen();

void BaseBetaGet(double sigma);

void GetCovarEff();

void GetEnvEff();

double IsCausal(char SNP[50]);

long int FindSNPinRef(char SNPname[50], int chr);

double GetMAF(int PopIndex);

double GetBeta(double BaseBeta, int PopIndex);

void AnalyzeSNP(int chr, char SNP[50]);
