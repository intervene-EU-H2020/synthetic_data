#include "Support.h"


int main(int argc, char const *argv[])
{
	if (argc < 2) {
		printf("No parameter file input. Terminate.\n");
		exit(0);
	}

	gsl_set_error_handler_off();
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	if (argc < 3) {
		printf("No random seed input. Using 142857.\n");
		gsl_rng_default_seed = 142857;
	}
	else
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
    if ( !NoRefFlag )
	    ReadRef();

	if (InterFlag)
		ReadInter();

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

	if (CausalEffFlag != 1) { // if user input causal effect sizes, do not output causal SNP file
		for (k = 0; k < nTrait; k++) {
			strcpy(tmpOutCausal, OutCausal);
			sprintf(tmpBuff, "%d", (int)k+1);
			OutFileCausal = fopen(strcat(tmpOutCausal, tmpBuff),"w");
			if (OutFileCausal == NULL) {
		        printf("Cannot open output causal file.\n");
		        exit(0);
		    } //check first wether the output file can be opened
		    fprintf(OutFileCausal, "SNP\tBaseEff\tMAF\tLDscore\tBeta\n");
			fclose(OutFileCausal);
		}
	}

	int CHR;
	char SNP[50];
	char InGenoCHR[1000][22];
	long int i_chr; // chr SNP counter

	i = 0; // SNP counter
    k = 0; // casual SNP counter

    for (CHR = 0; CHR < 22; CHR++) {
    	i_chr = 0;
    	sprintf(InGenoCHR[CHR],"%s-%d", InGeno, CHR+1);
	    if( pio_open( &InGenoPlink, InGenoCHR[CHR]) != PIO_OK ) {
			printf( "Cannot open plink file for chrmosome %d: %s, using other chromosomes.\n", CHR+1, InGenoCHR[CHR]);
		}
	    else if( !pio_one_locus_per_row(&InGenoPlink) ){
			printf( "This script requires that snps are rows and samples columns.\n" );
			exit(0);
		}
		else if (nSample != pio_num_samples(&InGenoPlink)) {
			printf("Plink file aample size does not match sample file.\n");
			exit(0);
		}
		else {
			j = 0; // sample counter
			BaseBetaGen();
			struct pio_sample_t *sample;
			struct pio_locus_t *locus;
			SNPbuffer = (snp_t *) malloc(pio_row_size(&InGenoPlink));
			for (j = 0; j < pio_num_samples(&InGenoPlink); j++) {
				sample = pio_get_sample( &InGenoPlink, j);
				strcpy(SampleList[j], sample->iid);
			}
			while (pio_next_row( &InGenoPlink, SNPbuffer) == PIO_OK) {
				locus = pio_get_locus( &InGenoPlink, i_chr);
				strcpy(SNP, locus->name);
				flag = IsCausal(SNP);
				if (flag) {
					memset(PopSampleCt, 0, sizeof(long int) * nMaxPop);
				    memset(PopMatTmp, 0, sizeof(double) * nMaxPop * nMaxInd * 3);
				    memset(GenoMat, 0, sizeof(double) * nMaxInd);
					for (j = 0; j < nSample; j++) {
						GenoMat[j] = ((SNPbuffer[j]==3.0) ? 0.0 : 2-SNPbuffer[j]); // counting the number of A1 allele in the bim file
						PopIndex = PopIndicator[j];
						tmpPopCt = PopSampleCt[PopIndex];
						PopMatTmp[PopIndex][0][tmpPopCt] = ((SNPbuffer[j]==3.0) ? 0.0 : 2-SNPbuffer[j]); // counting the number of A1 allele in the bim file
						PopSampleCt[PopIndex]++;
					}
				    AnalyzeSNP((int)locus->chromosome, SNP);
				    k++;
				}
				if (InterFlag) {
					flag = FindInterItem(SNP);
					if (flag != -1) {
						for (j = 0; j < nSample; j++) {
							InterItemMat[flag][j] = ((SNPbuffer[j]==3.0) ? 0.0 : 2-SNPbuffer[j]);
						}
					}
				}
				i++;
				i_chr++;
			}
			free(SNPbuffer);
			pio_close(&InGenoPlink);
		}
	}
	printf("Input genomat: %ld SNPs; Using in total %ld causal SNPs.\n", i, k);

	if (nCovar) {
		printf("Getting Cov Effect...\n");
		GetCovarEff();
		printf("done.\n");
	}


	if (InterFlag) {
		printf("Getting Interaction Effect...\n");
		GetInterEff();
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
			GCEweight[n][0] = (isnan(GCEweight[n][0]) || isinf(GCEweight[n][0]) ? 0.0 : GCEweight[n][0]);
			GCEweight[n][1] = (isnan(GCEweight[n][1]) || isinf(GCEweight[n][1]) ? 0.0 : GCEweight[n][1]);
			GCEweight[n][2] = (isnan(GCEweight[n][2]) || isinf(GCEweight[n][2]) ? 0.0 : GCEweight[n][2]);
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


