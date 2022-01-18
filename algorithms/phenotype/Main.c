#include "Support.h"


int main(int argc, char const *argv[])
{
	gsl_set_error_handler_off();
	const gsl_rng_type * T;
	gsl_rng_env_setup();
    gsl_rng_default_seed = atoi(argv[2]);
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    memset(GenoEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
    memset(CovarEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
    memset(EnvEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
    memset(VarGeno, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
    memset(VarCovar, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
    memset(VarEnv, 0.0, sizeof(double) * nMaxPop * nMaxTrait);

    ReadParam(argv[1]);

    char *tok; char *p;
    char tmpOutCausal[1000];
    char tmpOutPheno[1000];
    char tmpBuff[20];
    long int i, j, k, n;
    long int PopSampleCt[nMaxPop], tmpPopCt;
    double tmpPheMean, tmpPheVar;
    double tmpPheCutoff[nMaxPop];
    int SNPindex, PopIndex, flag;
    FILE *OutFileCausal;

    ExtractParam();
    ReadPopulation();
    ReadRef();

    if (!PolyFlag)
		ReadCausal();

    for (n = 0; n < nMaxTrait; n++) {
    	GenoEff[n] = malloc(sizeof(double) * nSample);
		memset(GenoEff[n], 0.0, sizeof(double) * nSample);
		CovarEff[n] = malloc(sizeof(double) * nSample);
		memset(CovarEff[n], 0.0, sizeof(double) * nSample);
		EnvEff[n] = malloc(sizeof(double) * nSample);
		memset(EnvEff[n], 0.0, sizeof(double) * nSample);
		for (i = 0; i < 2; i++) {
			PhenoSim[n][i] = malloc(sizeof(double) * nSample);
			memset(PhenoSim[n][i], 0.0, sizeof(double) * nSample);
		}
    }

	printf("memset, done.\n");

	for (k = 0; k < nTrait; k++) {
		strcpy(tmpOutCausal, OutCausal);
		sprintf(tmpBuff, "%d", (int)k+1);
		OutFileCausal = fopen(strcat(tmpOutCausal, tmpBuff),"w");
		if (OutFileCausal == NULL) {
	        printf("Cannot open output causal file.\n");
	        exit(0);
	    } //check first wether the output file can be opened
	    fprintf(OutFileCausal, "SNP\ttBaseEff\tMAF\tLDscore\tBeta\n");
		fclose(OutFileCausal);
	}

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
				flag = IsCausal(SNP);
				if (flag) {
					k++;
					tok = strtok_r(p, " ,\t", &p); // Morgan pos
				    tok = strtok_r(p, " ,\t", &p); // BP pos
				    tok = strtok_r(p, " ,\t", &p); // Ref Allele
				    tok = strtok_r(p, " ,\t", &p); // ALT Allele
				    j = 0;
				    memset(PopSampleCt, 0, sizeof(long int) * nMaxPop);
				    memset(PopMatTmp, 0, sizeof(double) * nMaxPop * nMaxInd * 3);
				    memset(GenoMat, 0, sizeof(double) * nMaxInd);
				    while ((tok = strtok_r(p, " ,\t", &p))) {
				    	GenoMat[j] = atof(tok);
				    	PopIndex = PopIndicator[j++];
				    	tmpPopCt = PopSampleCt[PopIndex];
				    	PopMatTmp[PopIndex][0][tmpPopCt] = atof(tok);
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

	if (nCovar) {
		printf("Getting Cov Effect...\n");
		GetCovarEff();
		printf("done.\n");
	}

	printf("Getting Env Effect...\n");
	GetEnvEff();
	printf("done.\n");

	for (k = 0; k < nTrait; k++) {
		memset(PopSampleCt, 0, sizeof(long int) * nMaxPop);
		memset(PopMatTmp, 0, sizeof(double) * nMaxPop * nMaxInd * 3);
		for (i = 0; i < nSample; i++) {
	    	PopIndex = PopIndicator[i];
	    	tmpPopCt = PopSampleCt[PopIndex];
	    	PopMatTmp[PopIndex][0][tmpPopCt] = GenoEff[k][i];
	    	PopMatTmp[PopIndex][1][tmpPopCt] = CovarEff[k][i];
	    	PopMatTmp[PopIndex][2][tmpPopCt] = EnvEff[k][i];
	    	PopSampleCt[PopIndex]++;
    	}
	    for (n = 0; n < nPop; n++) {
	    	VarGeno[n][k] = gsl_stats_variance(PopMatTmp[n][0], 1, PopSampleCt[n]);
	    	VarCovar[n][k] = gsl_stats_variance(PopMatTmp[n][1], 1, PopSampleCt[n]);
	    	VarEnv[n][k] = gsl_stats_variance(PopMatTmp[n][2], 1, PopSampleCt[n]);
	    }
	}

	FILE *OutFilePheno;
	for (k = 0; k < nTrait; k++) {
		strcpy(tmpOutPheno, OutPheno);
		sprintf(tmpBuff, "%d", (int)k+1);
		OutFilePheno = fopen(strcat(tmpOutPheno, tmpBuff),"w");
		if (OutFilePheno == NULL) {
	        printf("Cannot open output pheno file.\n");
	        exit(0);
	    }
	    if (BinaryFlag)
			fprintf(OutFilePheno, "Sample\tGenoEff\tCovarEff\tEnvEff\tPhenotype(liability)\tPhenotype(binary)\n");
		else 
			fprintf(OutFilePheno, "Sample\tGenoEff\tCovarEff\tEnvEff\tPhenotype\n");
		memset(GCEweight, 0.0, sizeof(double) * nMaxPop * 3);
		for (n = 0; n < nPop; n++) {
			GCEweight[n][0] = sqrt(GenoEffProp[n][k]/VarGeno[n][k]);
			GCEweight[n][1] = sqrt(CovarEffProp[n][k]/VarCovar[n][k]);
			GCEweight[n][2] = sqrt(EnvEffProp[n][k]/VarEnv[n][k]);
			// printf("GenoEffP = %lf, VarGeno = %lf, EnvEffP = %lf, VarEnv = %lf\n", GenoEffProp[n][k], VarGeno[n][k], EnvEffProp[n][k], VarEnv[n][k]);
			GCEweight[n][0] = (isnan(GCEweight[n][0]) ? 0.0 : GCEweight[n][0]);
			GCEweight[n][1] = (isnan(GCEweight[n][1]) ? 0.0 : GCEweight[n][1]);
			GCEweight[n][2] = (isnan(GCEweight[n][2]) ? 0.0 : GCEweight[n][2]);
		}
		memset(PopSampleCt, 0, sizeof(long int) * nMaxPop);
		memset(PopMatTmp, 0, sizeof(double) * nMaxPop * nMaxInd * 3);
		for (i = 0; i < nSample; i++) {
			PopIndex = PopIndicator[i];
			GenoEff[k][i] *= GCEweight[PopIndex][0];
			CovarEff[k][i] *= GCEweight[PopIndex][1];
			EnvEff[k][i] *= GCEweight[PopIndex][2];
			PhenoSim[k][0][i] = GenoEff[k][i] + CovarEff[k][i] + EnvEff[k][i];
	    	if (!BinaryFlag)
				fprintf(OutFilePheno, "%s\t%lf\t%lf\t%lf\t%lf\n", SampleList[i], GenoEff[k][i], CovarEff[k][i], EnvEff[k][i], PhenoSim[k][0][i]);
			else {
				tmpPopCt = PopSampleCt[PopIndex];
		    	PopMatTmp[PopIndex][0][tmpPopCt] = PhenoSim[k][0][i];
		    	PopSampleCt[PopIndex]++;
			}
		}

		if (BinaryFlag) {
			for (n = 0; n < nPop; n++) {
				tmpPheMean = gsl_stats_mean(PopMatTmp[n][0], 1, PopSampleCt[n]);
				tmpPheVar = gsl_stats_variance(PopMatTmp[n][0], 1, PopSampleCt[n]);
				tmpPheCutoff[n] = gsl_cdf_gaussian_Qinv(Prev[n][k], sqrt(tmpPheVar)) + tmpPheMean;
				printf("Trait %ld, Population %ld: Phe Mean = %lf, Var = %lf, Diagnosis Cutoff = %lf\n", k+1, n+1, tmpPheMean, tmpPheVar, tmpPheCutoff[n]);
			}
			for (i = 0; i < nSample; i++) {
				PopIndex = PopIndicator[i];
				PhenoSim[k][1][i] = ((PhenoSim[k][0][i] > tmpPheCutoff[PopIndex]) ? 1.0 : 0.0 );
				fprintf(OutFilePheno, "%s\t%lf\t%lf\t%lf\t%lf\t%d\n", SampleList[i], GenoEff[k][i], CovarEff[k][i], EnvEff[k][i], PhenoSim[k][0][i], (int) PhenoSim[k][1][i]); 
			}
		}
		fclose(OutFilePheno);
	}
	printf("Writing phenotypes, done.\n");
	return 0;
}


