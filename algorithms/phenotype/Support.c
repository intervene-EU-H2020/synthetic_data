#include "Support.h"

void ReadParam(const char *ParIn) {
	FILE *ParFile;
	char *tok; char *p;
	ParFile = fopen(ParIn,"r");
	memset(herr, 0.0, sizeof(double) * nMaxPop);
	memset(prev, 0.0, sizeof(double) * nMaxPop);
	a = 0;
	b = 0; 
	c = 0;
	nComp = 1;
	pCausal = 0;

	if (ParFile == NULL) {
	    printf("Cannot open parameter file.\n");
	    exit(0);
	}
	else {
		while (fgets(buffer, sizeof(buffer), ParFile) != NULL) {
			p = buffer;
			tok = strtok_r(p, " \t", &p);
	    	if (tok != NULL) {
	    		if (strcmp(tok, "Heritability") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpHerr,tok);
				}
				else if (strcmp(tok, "Prevalence") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpPrev,tok);
	    		}
	    		else if (strcmp(tok, "Polygenicity") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			pCausal = atof(tok);
	    		}
	    		else if (strcmp(tok, "PopCovar") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpCov,tok);
	    		}
	    		else if (strcmp(tok, "a") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			a = atof(tok);
	    		}
	    		else if (strcmp(tok, "b") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			b = atof(tok);
	    		}
	    		else if (strcmp(tok, "c") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			c = atof(tok);
	    		}
	    		else if (strcmp(tok, "nComponent") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			nComp = atoi(tok);
	    		}
	    		else if (strcmp(tok, "CausalList") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InCausal,tok);
	    		}
	    		else if (strcmp(tok, "Reference") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InRef,tok);
	    		}
	    		else if (strcmp(tok, "SampleList") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InSample,tok);
	    		}
	    		else if (strcmp(tok, "GenoFile") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InGeno,tok);
	    		}
	    		else if (strcmp(tok, "Output") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(OutCausal, tok);
					strcat(OutCausal, ".causal");
					strcpy(OutPheno, tok);
					strcat(OutPheno, ".pheno");
	    		}
	    		else if (strcmp(tok, "CompWeight") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpWeight,tok);
	    		}
	    	}
		}
		if ( !strcmp(InRef, "") || !strcmp(InGeno, "") || !strcmp(InSample, "")) {
			printf("Missing input files!\n");
			exit(0);
		}
		if ( !strcmp(InCausal, "") && pCausal<=0) {
			printf("Needs positive polygenicity or causal SNP list!\n");
			exit(0);
		}
		if ( strcmp(InCausal, "") && pCausal>0) {
			printf("Both polygenicity and causal SNP list input, ignore polygenicity.\n");
		}
		if ( strcmp(InCausal, "") )
			PolyFlag = 0;
		else
			PolyFlag = 1;
		if (!strcmp(OutCausal, "")) {
			strcpy(OutCausal, "Output");
			strcat(OutCausal, ".causal");
			strcpy(OutPheno, "Output");
			strcat(OutPheno, ".pheno");
		}
		if (nComp == 1) {
			strcpy(tmpProb,"1");
			strcpy(tmpWeight,"1");
		}
		else if ( !strcmp(tmpWeight, "") ) {
			printf("Missing component weight!\n");
			exit(0);
		}
	}
}


void ExtractParam() {
	char *tok; char *p;
	long int i, m, n;
	double SumProb;

	p = tmpHerr;
    i = 0; //weights counter
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nMaxPop) {
    		if (atof(tok) > 0.0) {
	    		herr[i] = atof(tok);
    		}
    		else
    			printf("Heritability needs to be positive!");
    	}
    	else {
    		printf("Increase max Population limit!\n");
    		exit(0);
    	}
    	i++;
	}
	nPop = i;

	p = tmpPrev;
    i = 0; //weights counter
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nMaxPop) {
    		if (atof(tok) > 0.0) {
	    		prev[i] = atof(tok);
    		}
    		else
    			printf("Prevalence needs to be positive!");
    	}
    	else {
    		printf("Increase max Population limit!\n");
    		exit(0);
    	}
    	i++;
	}
	if (i != nPop) {
		printf("Mismatch population number.\n");
	}


	if (nPop > 1) {
		mu = gsl_vector_calloc(nPop);
		for (i = 0; i < nPop; i++)
			gsl_vector_set(mu, i, 0.0); 

		Sigma = gsl_matrix_calloc(nPop, nPop);
		p = tmpCov;
	    m = 0; n = 0;
	    while ((tok = strtok_r(p, ",", &p))) {
	    	if (atof(tok) < 0.0) {
	    		printf("Negative covariance!\n");
	    		exit(0);
	    	}
	    	else
	    		gsl_matrix_set(Sigma, m, n, atof(tok)); 
	    	n += 1;
	    	if (n == nPop) {
	    		m += 1;
	    		n = 0;
	    	}
		}
		if (m != nPop || n != 0) {
			printf("Mismatch population number.\n");
		}
		int signum; 
		gsl_matrix *tmpSig = gsl_matrix_alloc(nPop, nPop);
		gsl_permutation *perm = gsl_permutation_alloc(nPop);
		gsl_linalg_LU_decomp(tmpSig, perm, &signum);
		det = gsl_linalg_LU_det(tmpSig, signum);
		if (det > 0) {
			L = gsl_matrix_calloc(nPop, nPop);
			gsl_matrix_memcpy(L, Sigma);
			gsl_linalg_cholesky_decomp1(Sigma);
		}
		else {
			printf("Cov Matrix not positive definite. Assuming all corrPop = 1.0.\n");
		}
	}


	p = tmpWeight;
    i = 0; //weights counter
    SumProb = 0;
    wComp = malloc(sizeof(double) * nComp);
    memset(wComp, 0, sizeof(double) * nComp);
    ProbComp = malloc(sizeof(double) * nComp);
    memset(ProbComp, 0, sizeof(double) * nComp);
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nComp) {
    		if (atof(tok) > 0.0) {
	    		wComp[i] = atof(tok);
	    		ProbComp[i] = 1.0/wComp[i];
	    		SumProb += ProbComp[i];
    		}
    		else
    			printf("Component weight needs to be positive!");
    	}
    	else {
    		printf("Wrong number of component weights!\n");
    		exit(0);
    	}
    	i++;
	}
	for (i = 0; i < nComp; i++) {
    	ProbComp[i] = ProbComp[i]/SumProb;
    	if (i > 0)
    		ProbComp[i] += ProbComp[i-1];
    }
}


int PopIndex(char PopCode[50], int nPopCt) {
	int i;
	for (i = 0; i < nPopCt; i++) {
		if (strcmp(PopCode, PopList[i]) == 0) {
			return(i);
		}
	}
	return(nPopCt);
}


// Sample list: one column of population code, in the order of sample list in the header of .traw input
void ReadPopulation() {
	char *tok; char *p;
	char PopCode[50];
	int nPopCt;
	long int i;
	nPopCt = 0;
	memset(nSamplePerPop, 0, sizeof(int) * nMaxPop);
	
	FILE *InFileSample;
    InFileSample = fopen(InSample, "r");
    if (InFileSample == NULL) {
        printf("Cannot open the Sample list %s.\n", InSample);
        exit(0);
    }
    else {
		i = 0; // Line counter
		while (fgets(buffer, sizeof(buffer), InFileSample) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); // Population Code
			strcpy(PopCode, tok);
			if (nPopCt == 0) {
				strcpy(PopList[nPopCt], PopCode);
				PopIndicator[i] = nPopCt;
				nPopCt += 1;
			}
			else {
				PopIndicator[i] = PopIndex(PopCode, nPopCt);
				if (PopIndicator[i] == nPopCt) {
					strcpy(PopList[nPopCt], PopCode);
					printf("Pop = %s\n", PopCode);
					nPopCt += 1;
				}
			}
			if (nPopCt >= nMaxPop) {
				printf("Increase max Population limit!\n");
				exit(0);
			}
			nSamplePerPop[PopIndicator[i]] += 1;
			i++;
		}
		fclose(InFileSample);
		nSample = i;
		printf("Input %ld samples from population list\n", nSample);
	}
	if (nPopCt != nPop) {
		printf("Mismatch population number.\n");
	}
}


void ReadCausal() {
	long int i;
	char *tok; char *p;

	FILE *InFileCausal;
    InFileCausal = fopen(InCausal, "r");
    if (InFileCausal == NULL) {
        printf("Cannot open the Causal SNP list %s.\n", InCausal);
        exit(0);
    }
    else {
		i = 0; // Line counter
		while (fgets(buffer, sizeof(buffer), InFileCausal) != NULL) {
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); // Causal SNP
			if ( strlen(tok) != 0) {
			    strcpy(CausalList[i], tok);
	    	}
			i++;
		}
		fclose(InFileCausal);
		printf("Input %ld causal SNPs from list.\n", i);
	}
	nCausal = i;
}


void ReadRef() {
	long int i, SNPindex;
	int chr;
	char *tok; char *p;

	memset(SNPct, 0, sizeof(int) * 26);

	FILE *InFileRef;
    InFileRef = fopen(InRef, "r");
    if (InFileRef == NULL) {
        printf("Cannot open the Reference file %s.\n", InRef);
        exit(0);
    }
    else {
    	i = 0; // Line counter
    	fgets(buffer, sizeof(buffer), InFileRef); //header
		while (fgets(buffer, sizeof(buffer), InFileRef) != NULL) {
			SNPinfo snp;
			p = buffer;
			tok = strtok_r(p, " ,\t\n", &p); // Chr
			if ( strlen(tok) != 0) {
			    chr = atoi(tok)-1;
	    	}
	    	tok = strtok_r(p, " ,\t\n", &p); // Causal SNP
			if ( strlen(tok) != 0) {
			    strcpy(snp.SNP, tok);
	    	}
			tok = strtok_r(p, " ,\t\n", &p); // A1
			tok = strtok_r(p, " ,\t\n", &p); // A2
			tok = strtok_r(p, " ,\t\n", &p); // MAF
			if ( strlen(tok) != 0) {
			    snp.AfricaMAF = atof(tok);
	    	}
	    	tok = strtok_r(p, " ,\t\n", &p); // LDscore
			if ( strlen(tok) != 0) {
			    snp.AfricaLDscore = atof(tok);
	    	}
	    	tok = strtok_r(p, " ,\t\n", &p); // Annot
	    	tok = strtok_r(p, " ,\t\n", &p); // Gene
	    	tok = strtok_r(p, " ,\t\n", &p); // Exon
	    	if ( strlen(tok) != 0) {
			    if (atoi(tok) != 0) 
			    	snp.exon = 1;
			    else
			    	snp.exon = 0;
	    	}
	    	tok = strtok_r(p, " ,\t\n", &p); // DHS
	    	if ( strlen(tok) != 0) {
			    snp.DHS = atoi(tok);
	    	}
			SNPindex = SNPct[chr];
	    	RefSNP[chr][SNPindex] = snp;
	    	SNPct[chr]++;
			i++;
		}
		fclose(InFileRef);
		printf("Read reference file, done. %ld SNPs read from reference.\n", i);
		/*for (i = 0; i < 26; i++)
			printf("Chr %ld has %ld SNPs.\n", i, SNPct[i]);*/
	}
}


double IsCausal(char SNP[50]) {
	if (PolyFlag) {
		prob = gsl_rng_uniform(r);
		if (prob < pCausal) {
			return(1);
		}
		else
			return(0);
	}
	else {
		int i; 
		for (i = 0; i < nCausal; i++) {
			if (strcmp(SNP, CausalList[i]) == 0) {
				return(1);
			}
		}
		return(0);
	}
}


long int FindSNPinRef(char SNP[50], int chr) {
	int i; 
	for (i = 0; i < SNPct[chr]; i++) {
		// printf("i = %d, SNP = %s, total SNP on chr = %d\n", i, RefSNP[chr][i].SNP, SNPct[chr]);
		if (strcmp(SNP, RefSNP[chr][i].SNP) == 0) {
			return(i);
		}
	}
	return(-1);
}


double GetMAF(int PopIndex) {
	double freq;
	freq = gsl_stats_mean(GenoMat[PopIndex], 1, nSamplePerPop[PopIndex])/2.0;
	freq = (freq > 0.5) ? (1-freq) : freq;
	return(freq);
}


double GetBeta(double BaseBeta, int PopIndex) {
	if (CausalMAF[PopIndex] == 0.0) {
		return(0.0);
	}
	else if (CausalLDscore == 0.0) {
		CausalLDscore = 0.0001; // will run into problem if divide 0; therefore use 0.0001 as minimum
	}
	double fMAF = pow(CausalMAF[PopIndex]*(1-CausalMAF[PopIndex]), a);
	double fLD = pow(CausalLDscore, b);
	double fScore = pow(CausalAnnot, c);
	double sigma2 = fMAF*fLD*fScore;
	return(BaseBeta * sqrt(sigma2));
}


void AnalyzeSNP(int chr, char SNP[50]) {
	double prob;
	gsl_vector * BaseBeta = gsl_vector_calloc(nPop);
	long int i, j, k, n, SNPindex;
	char tmp[50];
	char tmpMAF[500];
	char tmpBeta[500];
	char tmpBaseBeta[500];
    prob = gsl_rng_uniform(r);
    if (prob < ProbComp[0]) {
    	k = 0;
    }
    else {
		for (j = 1; j < nComp; j++) {
			if ((prob > ProbComp[j-1]) && (prob < ProbComp[j])) {
				k = j;
				break;
			}
		}
    }

    if (nPop == 1)
    	gsl_vector_set(BaseBeta, 0, gsl_ran_gaussian(r, wComp[k]));
    else if (det > 0) {
	    gsl_ran_multivariate_gaussian(r, mu, L, BaseBeta);
	}
	else {
		double tmp = gsl_ran_gaussian(r, wComp[k]);
		for (i = 0; i < nPop; i++) {
			gsl_vector_set(BaseBeta, i, tmp);
		}
	}

    memset(tmpMAF, '\0', sizeof tmpMAF);
    memset(tmpBeta, '\0', sizeof tmpBeta);
    memset(tmpBaseBeta, '\0', sizeof tmpBaseBeta);

    SNPindex = FindSNPinRef(SNP, chr);
    if (SNPindex == -1) {
		CausalLDscore = 0.25;
		CausalAnnot = 1.0;
	}
	else {
    	CausalLDscore = RefSNP[chr][SNPindex].AfricaLDscore;
    	CausalAnnot = RefSNP[chr][SNPindex].DHS + 1;
    }
    for (n = 0; n < nPop; n++) {
    	CausalMAF[n] = GetMAF(n);
    	CausalBeta[n] = GetBeta(gsl_vector_get(BaseBeta, n), n);
    	for (i = 0; i < nSamplePerPop[n]; i++)
    		GenoEff[n][i] += CausalBeta[n]*GenoMat[n][i];

		sprintf(tmp, "%f,", CausalMAF[n]);
		strcat(tmpMAF, tmp);
		sprintf(tmp, "%f,", CausalBeta[n]);
		strcat(tmpBeta, tmp);
		sprintf(tmp, "%f,", gsl_vector_get(BaseBeta, n));
		strcat(tmpBaseBeta, tmp);
    }
    tmpMAF[strlen(tmpMAF)-1] = '\0';
    tmpBeta[strlen(tmpBeta)-1] = '\0';
    tmpBaseBeta[strlen(tmpBaseBeta)-1] = '\0';
    FILE *OutFileCausal;
	OutFileCausal = fopen(OutCausal,"a");
    fprintf(OutFileCausal, "%s\t%s\t%s\t%lf\t%s\n", SNP, tmpBaseBeta, tmpMAF, CausalLDscore, tmpBeta);
	fclose(OutFileCausal);
}



