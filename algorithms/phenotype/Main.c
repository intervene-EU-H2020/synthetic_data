#include "Support.h"


int main(int argc, char const *argv[])
{
	const gsl_rng_type * T;
	gsl_rng_env_setup();
    gsl_rng_default_seed = atoi(argv[2]);
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    ReadParam(argv[1]);

    char *tok; char *p;
    long int i, j, k, n;
    long int PopSampleCt[nMaxPop], tmpPopCt;
    int SNPindex, PopIndex;
    ExtractParam();
    ReadPopulation();
    ReadRef();
    if (!PolyFlag)
		ReadCausal();

    for (n = 0; n < nMaxPop; n++) {
    	GenoEff[n] = malloc(sizeof(double) * nSamplePerPop[n]);
		memset(GenoEff[n], 0.0, sizeof(double) * nSamplePerPop[n]);
		EnvEff[n] = malloc(sizeof(double) * nSamplePerPop[n]);
		memset(EnvEff[n], 0.0, sizeof(double) * nSamplePerPop[n]);
		PhenoSim[n] = malloc(sizeof(double) * nSamplePerPop[n]);
		memset(PhenoSim[n], 0.0, sizeof(double) * nSamplePerPop[n]);
    }

	printf("memset, done.\n");

	FILE *OutFileCausal;
	OutFileCausal = fopen(OutCausal,"w");
	if (OutFileCausal == NULL) {
        printf("Cannot open output file.\n");
        exit(0);
    } //check first wether the output file can be opened
    fprintf(OutFileCausal, "SNP\ttBaseEff\tMAF\tLDscore\tBeta\n");
	fclose(OutFileCausal);

// takes traw as input
    FILE *InFileGeno;
    int CHR;
    char SNP[50];

    strcat(InGeno, ".traw");
    InFileGeno = fopen(InGeno, "r");
    if (InFileGeno == NULL) {
        printf("Cannot open the traw file %s.\n", InGeno);
        exit(0);
    }
    else {
		i = 0; // SNP counter
		j = 0; // Sample counter; PopSampleCt for population sample counter
		k = 0; // Causal SNP counter

		fgets(buffer, sizeof(buffer), InFileGeno); // read header
		p = buffer;
		tok = strtok_r(p, " ,\t", &p); // CHR
		tok = strtok_r(p, " ,\t", &p); // SNP
		tok = strtok_r(p, " ,\t", &p); // Morgan pos
	    tok = strtok_r(p, " ,\t", &p); // BP pos
	    tok = strtok_r(p, " ,\t", &p); // Ref Allele
	    tok = strtok_r(p, " ,\t", &p); // ALT Allele
	    while ((tok = strtok_r(p, " ,\t\n", &p))) {
	    	strcpy(SampleList[j++], tok);
	    }
	    if (j != nSample) {
	    	printf("Sample size in genotype file does not match population file, nSample = %ld, j = %ld.\n", nSample, j);
	    	exit(0);
	    }
		while (fgets(buffer, sizeof(buffer), InFileGeno) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t", &p); //CHR
			CHR = atoi(tok);
			if (CHR < 1 || CHR > 26) {
		    	printf("line %ld has invalid CHR code.\n", i+1);
		        exit(0);
		    }
			else {
				CHR = CHR-1;
				tok = strtok_r(p, " ,\t", &p); // SNP ID
				strcpy(SNP, tok);
				CausalFlag = IsCausal(SNP);
				if (CausalFlag) {
					k++;
					tok = strtok_r(p, " ,\t", &p); // Morgan pos
				    tok = strtok_r(p, " ,\t", &p); // BP pos
				    tok = strtok_r(p, " ,\t", &p); // Ref Allele
				    tok = strtok_r(p, " ,\t", &p); // ALT Allele
				    j = 0;
				    memset(PopSampleCt, 0, sizeof(long int)*nMaxPop);
				    while ((tok = strtok_r(p, " ,\t", &p))) {
				    	PopIndex = PopIndicator[j++];
				    	tmpPopCt = PopSampleCt[PopIndex];
				    	GenoMat[PopIndex][tmpPopCt] = atof(tok);
				    	PopSampleCt[PopIndex]++;
				    }
				    AnalyzeSNP(CHR, SNP);
				}
				i++;
			}
		}
		fclose(InFileGeno);
		printf("Input genomat: %ld SNPs; Using in total %ld causal SNPs.\n", i, k);
	}	

	for (n = 0; n < nPop; n++) {
		VarGeno[n] = gsl_stats_variance(GenoEff[n], 1, nSamplePerPop[n]);
		VarEnv[n] = VarGeno[n]/herr[n] - VarGeno[n];
		for (j = 0; j < nSamplePerPop[n]; j++) {
			EnvEff[n][j] = gsl_ran_gaussian(r, sqrt(VarEnv[n]));
			PhenoSim[n][j] = EnvEff[n][j] + GenoEff[n][j];
		}
		printf("Population %s: VarGeno = %lf, VarEnv = %lf, VarPheno = %lf\n", PopList[n], VarGeno[n], VarEnv[n], gsl_stats_variance(PhenoSim[n], 1, nSamplePerPop[n]));
	}

	memset(PopSampleCt, 0, sizeof(long int)*nMaxPop);
	FILE *OutFilePheno;
	OutFilePheno = fopen(OutPheno,"w");
	if (OutFilePheno == NULL) {
        printf("Cannot open output file %s.\n", OutCausal);
        exit(0);
    }
    fprintf(OutFilePheno, "Sample\tGenoEff\tEnvEff\tPhenotype\n");
    for (j = 0; j < nSample; j++) {
    	PopIndex = PopIndicator[j];
    	tmpPopCt = PopSampleCt[PopIndex];
		fprintf(OutFilePheno, "%s\t%lf\t%lf\t%lf\n", SampleList[j], GenoEff[PopIndex][tmpPopCt], EnvEff[PopIndex][tmpPopCt], PhenoSim[PopIndex][tmpPopCt]);
		PopSampleCt[PopIndex]++;
    }
	fclose(OutFilePheno);
	return 0;
}


