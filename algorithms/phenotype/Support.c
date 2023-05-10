#include "Support.h"

void ReadParam(const char *ParIn) {
	PolyFlag = 0;
	BinaryFlag = 0;
	CausalEffFlag = 0;
	NoRefFlag = 0;
	InterFlag = 0;
	FILE *ParFile;
	char *tok; char *p;
	memset(GenoEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
	memset(CovarEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
	memset(pCausal, 0, sizeof(double) * nMaxTrait);
	a = 0;
	b = 0; 
	c = 0;
	nComp = 1;

	ParFile = fopen(ParIn,"r");
	if (ParFile == NULL) {
	    printf("Cannot open parameter file.\n");
	    exit(0);
	}
	else {
		while (fgets(buffer, sizeof(buffer), ParFile) != NULL) {
			p = buffer;
			tok = strtok_r(p, " \t", &p);
	    	if (tok != NULL) {
	    		if (strcmp(tok, "nPopulation") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			nPop = atoi(tok);
	    		}
	    		else if (strcmp(tok, "nTrait") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			nTrait = atoi(tok);
	    		}
	    		else if (strcmp(tok, "CausalEffect") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			if (strcmp(tok, "T") == 0) {
	    				CausalEffFlag = 1;
	    			}
	    		}
	    		// Heritability matrix: flatten nPop * nTrait matrix 
	    		else if (strcmp(tok, "PropotionGeno") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpGenoEffProp,tok);
				}
				// Covariate variance ratio matrix: flatten nPop * nTrait matrix 
				else if (strcmp(tok, "PropotionCovar") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpCovarEffProp,tok);
				}
				// Prevalence matrix: flatten nPop * nTrait matrix 
				else if (strcmp(tok, "Prevalence") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpPrev,tok);
	    		}
	    		// Polygenicity vector: nTrait entries 
	    		else if (strcmp(tok, "Polygenicity") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpPCausal,tok);
	    		}
	    		// Pleiotropy vector: nTrait entries 
	    		else if (strcmp(tok, "Pleiotropy") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpPleio,tok);
	    		}
	    		// Trait correlation matrix: flattern nTrait * nTrait matrix
	    		else if (strcmp(tok, "TraitCorr") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpTraitCorr,tok);
	    		}
	    		// Population correlation matrix: flattern nPop * nPop matrix
	    		else if (strcmp(tok, "PopulationCorr") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(tmpPopCorr,tok);
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
	    			strcpy(InCausal, tok);
	    		}
	    		// (hopefully) accounting for SNP x SNP interation, SNP x covariate interaction and covariate x covariate interatction
	    		// Input file looks like 3 col: term1, term2, non-optional effect size for term1 x term2. Term 1 and 2 can be SNP or covar
	    		// How to count the variance explained by SNP x covar interaction, more parameters..?
	    		else if (strcmp(tok, "Interaction") == 0) { 
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InInteract, tok);
	    			InterFlag = 1;
	    		}
	    		else if (strcmp(tok, "Reference") == 0) {
	    			tok = strtok_r(p, " \t\n", &p);
	    			strcpy(InRef,tok);
	    		}
	    		// Sample list: 1st column is string population code, folllowing columns are numeric covariates
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

		if (!strcmp(InGeno, "") || !strcmp(InSample, "")) {
			printf("Missing input files!\n");
			exit(0);
		}
		if ( !strcmp(InCausal, "") && CausalEffFlag) {
			printf("No causal SNP list input. CausalEff flag disabled.\n");
			CausalEffFlag = 0;
		}
		if ( !strcmp(InRef, "") ) {
			printf("No reference file specified. SNP LD and annotation effect disabled.\n");
			NoRefFlag = 1;
		}
		if ( !strcmp(InCausal, "") && !strcmp(tmpPCausal, "") ) {
			printf("Needs polygenicity parameter or causal SNP list input!\n");
			exit(0);
		}
		if ( strcmp(InCausal, "") && !strcmp(tmpPCausal, "") ) {
			printf("Both polygenicity and causal SNP list input, ignore polygenicity.\n");
		}
		if ( strcmp(InCausal, "") )
			PolyFlag = 0;
		else
			PolyFlag = 1;

		if ( !strcmp(tmpGenoEffProp, "") ) {
			printf("No heritability (PropotionGeno) input. Assuming 0.0. \n");
		}
		if ( !strcmp(tmpCovarEffProp, "") ) {
			printf("No covariate variance (PropotionCovar) input. Assuming 0.0. \n");
		}
		if (!strcmp(OutCausal, "")) {
			printf("No output prefix specified, using default prefix Output.\n");
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
			printf("Missing component weight! Assuming 1 component used.\n");
			nComp = 1;
			strcpy(tmpProb,"1");
			strcpy(tmpWeight,"1");
		}
	}
}


void ExtractParam() {
	char *tok; char *p;
	long int i, k, m, n;
	double SumProb, tmp;

	nItem = nTrait * nPop;
	statusPop = 1;
	statusTrait = 1;
	status = 1;

	// Pop(row) * Trait(col) genetic effect propotion matrix, flatten
	p = tmpGenoEffProp;
    i = 0; k = 0;
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nPop && k < nTrait) {
    		if (atof(tok) > 0.0 && atof(tok) < 1.0) {
	    		GenoEffProp[i][k] = atof(tok);
	    		k += 1;
	    		if (k == nTrait) {
	    			i += 1;
	    			k = 0;
	    		}
    		}
    		else {
    			printf("Genetic effect propotion needs to be between 0.0 and 1.0!");
    			exit(0);
    		}
    	}
    	else {
    		printf("Wrong dimension for input flatten heritability matrix!\n");
    		exit(0);
    	}
	}

	// Pop(row) * Trait(col) covar effect propotion matrix, flatten
	p = tmpCovarEffProp;
    i = 0; k = 0;
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nPop && k < nTrait) {
    		if (atof(tok) >= 0.0 && atof(tok) <= 1.0) {
	    		CovarEffProp[i][k] = atof(tok);
	    		k += 1;
	    		if (k == nTrait) {
	    			i += 1;
	    			k = 0;
	    		}
    		}
    		else {
    			printf("Covariate effect propotion needs to be between 0.0 and 1.0!");
    			exit(0);
    		}
    	}
    	else {
    		printf("Wrong dimension for input flatten covariate prop matrix!\n");
    		exit(0);
    	}
	}

	for (i = 0; i < nPop; i++) {
		for (k = 0; k < nTrait; k++) {
			tmp = 1 - GenoEffProp[i][k] - CovarEffProp[i][k];
			if (tmp < 0.0) {
				printf("Genetic effect  plus covariate effect propotion needs to be less than 1.0!");
				exit(0);
			}
			else
				EnvEffProp[i][k] = tmp;
		}
	}

	p = tmpPrev;
    i = 0; k = 0; tmp = 0;
    while ((tok = strtok_r(p, ",", &p))) {
    	if (i < nPop && k < nTrait) {
    		if (atof(tok) > 0.0 && atof(tok) < 1.0) {
	    		Prev[i][k] = atof(tok);
	    		tmp += Prev[i][k];
	    		k += 1;
	    		if (k == nTrait) {
	    			i += 1;
	    			k = 0;
	    		}
    		}
    		else {
    			printf("Prevalence needs to be between 0.0 and 1.0!");
    			exit(0);
    		}
    	}
    	else {
    		printf("Wrong dimension for input flatten prevalence matrix!\n");
    		exit(0);
    	}
	}
	BinaryFlag = (tmp ? 1 : 0);

	p = tmpPleio;
    k = 0;
    while ((tok = strtok_r(p, ",", &p))) {
    	if (k < nTrait) {
    		if (k == 0) {
    			if (atof(tok) != 1.0)
					printf("Trait 1 is reference, pleiotropy for trait 1 is set to be 1.0.\n");
    			Pleio[k] = 1.0;
				k += 1;
    		}
    		else if (atof(tok) >= 0.0 && atof(tok) <= 1.0) {
	    		Pleio[k] = atof(tok);
	    		k += 1;
    		}
    		else {
    			printf("Pleiotropy needs to be between 0.0 and 1.0!\n");
    			exit(0);
    		}
    	}
    	else {
    		printf("Wrong length for input pleiotropy vector!\n");
    		exit(0);
    	}
	}
	
	p = tmpPCausal;
    k = 0;
    while ((tok = strtok_r(p, ",", &p))) {
    	if (k < nTrait) {
    		pCausal[k] = atof(tok);
			pCausal[k] = ((k == 0) ? pCausal[k] : (pCausal[k] - pCausal[0]*Pleio[k]));
    		if (pCausal[k] < 0.0 || pCausal[k] > 1.0) {
    			printf("Given the pleiotropic model, polygenicity needs to be between 0.0 and 1.0!");
    			exit(0);
    		}
    		k += 1;
    	}
    	else {
    		printf("Wrong length for input polygenicity vector!\n");
    		exit(0);
    	}
	}

	p = tmpTraitCorr;
	m = 0; n = 0;
	while ((tok = strtok_r(p, ",", &p))) {
		if (m < nTrait && n < nTrait) {
			if (atof(tok) >= -1.0 && atof(tok) <= 1.0) {
				TraitCorr[m][n] = atof(tok);
				n += 1;
				if (n == nTrait) {
					m += 1;
					n = 0;
				}
	    	}
	    	else {
	    		printf("Trait correlation out of bound!\n");
	    		exit(0);
	    	}
	    }
	    else {
    		printf("Wrong dimension for input flatten trait correlation matrix!\n");
    		exit(0);
    	}
	}

	p = tmpPopCorr;
	m = 0; n = 0;
	while ((tok = strtok_r(p, ",", &p))) {
		if (m < nPop && n < nPop) {
			if (atof(tok) >= -1.0 && atof(tok) <= 1.0) {
				PopCorr[m][n] = atof(tok);
				n += 1;
				if (n == nPop) {
					m += 1;
					n = 0;
				}
	    	}
	    	else {
	    		printf("Population correlation out of bound!\n");
	    		exit(0);
	    	}
	    }
	    else {
    		printf("Wrong dimension for input flatten population correlation matrix!\n");
    		exit(0);
    	}
	}

	if (nItem > 1)
		MakeCovMat();
	else
		nValidItem = 1;

	GenoBeta = gsl_matrix_calloc(nValidItem, nMaxBetaGen);
	tmpCorrGenoBeta = gsl_matrix_calloc(nValidItem, nMaxBetaGen);
	tmpL = gsl_matrix_calloc(nValidItem, nValidItem);

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


// Order of variables: Trait 1 Pop1,2,3 ... Trait 2 Pop 1,2,3 ... Trait 3 ...
void MakeCovMat() {
	int i, j, m, n, nRow, nCol;
	double tmp; 

	gsl_matrix * Sigma = gsl_matrix_calloc(nItem, nItem);
	for (i = 0; i < nTrait; i++) {
		for (j = 0; j < nPop; j++){
			nRow = nPop * i + j;
			for (m = 0; m < nTrait; m++) {
				for (n = 0; n < nPop; n++){
					nCol = nPop * m + n;
					tmp = PopCorr[j][n] * TraitCorr[i][m];
					gsl_matrix_set(Sigma, nRow, nCol, tmp);
				}
			}
		}
	}
	status = gsl_linalg_cholesky_decomp1(Sigma);

	SigmaPop = gsl_matrix_calloc(nPop, nPop);
	for (i = 0; i < nPop; i++) {
		for (j = 0; j < nPop; j++)
			gsl_matrix_set(SigmaPop, i, j, PopCorr[i][j]);
	}
	gsl_permutation *permPop = gsl_permutation_alloc(nPop);
	statusPop = gsl_linalg_cholesky_decomp1(SigmaPop);

	SigmaTrait = gsl_matrix_calloc(nTrait, nTrait);
	for (m = 0; m < nTrait; m++) {
		for (n = 0; n < nTrait; n++)
			gsl_matrix_set(SigmaTrait, m, n, TraitCorr[m][n]);
	}
	statusTrait = gsl_linalg_cholesky_decomp1(SigmaTrait);

	if (!status) {
		printf("Overall correlation matrix is positive definite, ok!\n");
		nValidItem = nItem;
		L = gsl_matrix_calloc(nValidItem, nValidItem);
		gsl_matrix_memcpy(L, Sigma);
	}
	else if (statusPop && !statusTrait) {
		printf("Population correlation matrix is not positive definite. Assuming all population corr = 1.0.\n");
		nValidItem = nTrait;
		L = gsl_matrix_calloc(nValidItem, nValidItem);
		gsl_matrix_memcpy(L, SigmaTrait);
	}
	else if (statusTrait && !statusPop) {
		printf("Trait correlation matrix is not positive definite. Assuming all trait corr = 1.0.\n");
		nValidItem = nPop;
		L = gsl_matrix_calloc(nValidItem, nValidItem);
		gsl_matrix_memcpy(L, SigmaPop);
	}
	else if (statusPop && statusTrait) {
		printf("Neither correlation matrix is positive definite. Assuming all corr = 1.0.\n");
		nValidItem = 1;
	}
	else {
		// If not possible to satisfy both, always try to meet the trait correlation.
		printf("Both correlation matrices are positive definite, but overall correlation matrix is not. Assuming pop corr = 1.0.\n");
		nValidItem = nTrait;
		L = gsl_matrix_calloc(nTrait, nTrait);
		gsl_matrix_memcpy(L, SigmaTrait);
	}
	mu = gsl_vector_calloc(nValidItem);
	gsl_vector_set_zero(mu);
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
	long int i, j;
	nPopCt = 0;
	nCovar = 0;
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
					nPopCt += 1;
				}
			}
			if (nPopCt >= nMaxPop) {
				printf("Increase max Population limit!\n");
				exit(0);
			}
			nSamplePerPop[PopIndicator[i]] += 1;

			j = 0;
			while ((tok = strtok_r(p, " ,\t\n", &p))) {
				if (j < nMaxCovar)
					CovarMat[j++][i] = atof(tok);
				else
					printf("Increase max Covariate limit!\n");
		    }
		    if (nCovar == 0)
		    	nCovar = j;
		    else if (j != nCovar)
		    	printf("Line %ld does not have %d covariates.\n", i+1, nCovar);
		    i++;
		}
		fclose(InFileSample);
		nSample = i;
		printf("Input %ld samples from population list\n", nSample);
	}
	if (nPopCt != nPop) {
		printf("Mismatch population number.\n");
		exit(0);
	}
	if (nCovar == 0) {
		printf("No covariate input. Assigning covatiate effect variance to 0.0.\n");
		memset(CovarEffProp, 0.0, sizeof(double) * nMaxPop * nMaxTrait);
	}

}


void ReadCausal() {
	long int i, k;
	char *tok; char *p;
	char tmpInCausal[1000];
	FILE *InFileCausal;
	memset(nCausal, 0, sizeof(long int) * nMaxTrait);
	memset(CausalEff, 0.0, sizeof(double) * nMaxTrait * nMaxCausal * 2);
	for (k = 0; k < nTrait; k++) {
		strcpy(tmpInCausal, InCausal);
		sprintf(tmpBuff, "%d", (int)k+1);
	    InFileCausal = fopen(strcat(tmpInCausal, tmpBuff), "r");
	    if (InFileCausal == NULL) {
	        printf("Cannot open the Causal SNP list %s.\n", tmpInCausal);
	        exit(0);
	    }
	    else {
	    	while (fgets(buffer, sizeof(buffer), InFileCausal) != NULL) {
	    		p = buffer;
				tok = strtok_r(p, " ,\t\n", &p);
				if ( tok[0] != '\0' ) {
				    strcpy(CausalList[k][nCausal[k]], tok);
				    if (CausalEffFlag == 1) {
				    	tok = strtok_r(p, " ,\t\n", &p);
				    	if ( tok[0] != '\0' ) {
							CausalEff[k][nCausal[k]][0] = atof(tok);
						}
						else {
							printf("Trait %ld SNP %s missing homogenous causal effect size. Assigning 0.0.\n", k, CausalList[k][nCausal[k]]);
							CausalEff[k][nCausal[k]][0] = 0.0;
						}
						tok = strtok_r(p, " ,\t\n", &p);
				    	if ( tok[0] != '\0' ) {
							CausalEff[k][nCausal[k]][1] = atof(tok);
						}
						else {
							printf("Trait %ld SNP %s missing heterogenous causal effect size. Assigning 0.0.\n", k, CausalList[k][nCausal[k]]);
							CausalEff[k][nCausal[k]][1] = 0.0;
						}	
				    }
				    nCausal[k] += 1;
		    	}
	    	}
			printf("Trait %ld: input %ld causal SNPs.\n", k, nCausal[k]);
			fclose(InFileCausal);
		}
	}
	if (CausalEffFlag == 1)
		printf("Population and trait genetic correlation disabled.\nUsing user input trait causal effects.\n");
}


int FindInterItem(char item[50]) {
	int i;
	for (i = 0; i < nTotInterItem; i++) {
		if (strcmp(item, InterItemList[i]) == 0)
			return(i);
	}
	return(-1);
}


void ReadInter() {
	long int i, k;
	char *tok; char *p;
	char tmpInInteract[1000];
	FILE *InFileInteract;
	memset(nInter, 0, sizeof(int) * nMaxTrait);
	nTotInterItem = 0;
	for (k = 0; k < nTrait; k++) {
		strcpy(tmpInInteract, InInteract);
		sprintf(tmpBuff, "%d", (int)k+1);
	    InFileInteract = fopen(strcat(tmpInInteract, tmpBuff), "r");
	    if (InFileInteract == NULL) {
	        printf("Cannot open the interaction list %s.\n", tmpInInteract);
	    }
	    else {
	    	while (fgets(buffer, sizeof(buffer), InFileInteract) != NULL) {
	    		p = buffer;
				tok = strtok_r(p, " ,\t\n", &p);
				if ( tok[0] != '\0' ) {
				    strcpy(InterList[k][nInter[k]].Term1, tok);
				    if (FindInterItem(tok) == -1) {
				    	strcpy(InterItemList[nTotInterItem], tok);
				    	nTotInterItem += 1;
				    }
				}
				tok = strtok_r(p, " ,\t\n", &p);
				if ( tok[0] != '\0' ) {
				    strcpy(InterList[k][nInter[k]].Term2, tok);
				    if (FindInterItem(tok) == -1) {
				    	strcpy(InterItemList[nTotInterItem], tok);
				    	nTotInterItem += 1;
				    }
				}
				else {
					printf("Missing 2nd interaction term. Skip.\n");
					strcpy(InterList[k][nInter[k]].Term2, "NaN");
				}
				tok = strtok_r(p, " ,\t\n", &p);
				if ( tok[0] != '\0' )
				    InterList[k][nInter[k]].InterEff = atof(tok);
				else {
					printf("Missing interaction effect size. Assigning 0.0.\n");
					InterList[k][nInter[k]].InterEff = 0.0;
				}
				nInter[k] += 1;
	    	}
			printf("Trait %ld: input %d interaction terms.\n", k, nInter[k]);
			fclose(InFileInteract);
		}
	}
	if (nTotInterItem != 0) {
		printf("Total %d interaction terms.\n", nTotInterItem);
		InterItemMat = malloc(nSample * sizeof(double *));
		for (i = 0; i < nTotInterItem; i++) {
			InterItemMat[i] = malloc(nSample * sizeof(double));
			memset(InterItemMat[i], 0.0, sizeof(double) * nSample);
		}
	}
	else {
		printf("No interaction terms input. Interaction disabled.\n");
		InterFlag = 0;
	}
}



void ReadRef() {
	long int i, SNPindex;
	int chr;
	char *tok; char *p;

	memset(SNPct, 0, sizeof(int) * 26);

	FILE *InFileRef;
    InFileRef = fopen(InRef, "r");
    if (InFileRef == NULL) {
        printf("Cannot open the Reference file %s. Disable LD and annotation effect.\n", InRef);
        NoRefFlag = 1;
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
	}
}


long int findSNPinCausalList(char SNP[50], int k) {
	long int i; 
	for (i = 0; i < nCausal[k]; i++) {
		if (strcmp(SNP, CausalList[k][i]) == 0)
			return(i);
	}
	return(-1);
}


void BaseBetaGen() {
	int i, j;
	nBetaIndex = 0;
	gsl_matrix_set_zero(GenoBeta);
	for (i = 0; i < nValidItem; i++) {
		for (j = 0; j < nMaxBetaGen; j++)
			gsl_matrix_set(GenoBeta, i, j, gsl_ran_gaussian(r, 1.0));
	}

	if (nValidItem > 1) {
		gsl_matrix_set_zero(tmpCorrGenoBeta);
		gsl_matrix_set_zero(tmpL);
		for (i = 0; i < nValidItem; i++) {
			for (j = 0; j <= i; j++) {
				gsl_matrix_set(tmpL, i, j, gsl_matrix_get(L, i, j));
			}
		}
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmpL, GenoBeta, 0.0, tmpCorrGenoBeta);
		gsl_matrix_memcpy(GenoBeta, tmpCorrGenoBeta);
	}
}


void BaseBetaGet(double sigma) {
	int i, j;
	memset(BaseBeta, 0.0, sizeof(double)*nMaxPop*nMaxTrait);

	if (!status) {
		for (i = 0; i < nTrait; i++) {
			for (j = 0; j < nPop; j++)
				BaseBeta[j][i] = gsl_matrix_get(GenoBeta, nPop*i+j, nBetaIndex) * sigma;
		}
	}
	else if (!statusTrait) {
		for (i = 0; i < nTrait; i++) {
			for (j = 0; j < nPop; j++) {
				BaseBeta[j][i] = gsl_matrix_get(GenoBeta, i, nBetaIndex) * sigma;
			}
		}
	}
	else if (!statusPop) {
		for (i = 0; i < nPop; i++) {
			for (j = 0; j < nTrait; j++) {
				BaseBeta[i][j] = gsl_matrix_get(GenoBeta, i, nBetaIndex) * sigma;
			}
		}
	}
	else {
		for (i = 0; i < nPop; i++) {
			for (j = 0; j < nTrait; j++)
				BaseBeta[i][j] = gsl_matrix_get(GenoBeta, 0, nBetaIndex) * sigma;
		}
	}
	nBetaIndex++;
	if (nBetaIndex == nMaxBetaGen) 
		BaseBetaGen();
}



void GetCovarEff() {
	int i, j, k, popIndex;
	long int l;
	for (k = 0; k < nCovar; k++) {
		BaseBetaGet(1.0);
		for (i = 0; i < nTrait; i++) {
			for (l = 0; l < nSample; l++) {
				popIndex = PopIndicator[l];
				CovarEff[i][l] += BaseBeta[popIndex][i] * CovarMat[k][l];
			}
		}
	}
}


// SNP x SNP interaction is considered genetic effect; SNP x covariate, covariate x covariate interaction is considered covariate effect
// (Sensible..?)
void GetInterEff() {
	int i, j, k;
	int CovarIndex[2];
	int SNPindex[2];
	int Term1SNP; int Term2SNP;

	printf("Note: SNP x SNP interaction is considered genetic effect; SNP x covariate, covariate x covariate interactions are considered covariate effect.\n");
	for (k = 0; k < nTrait; k++) {
		for (i = 0; i < nInter[k]; i++) {
			memset(CovarIndex, 0, sizeof(int) * 2);
			memset(SNPindex, 0, sizeof(int) * 2);
			if (strncmp("covar", InterList[k][i].Term1, strlen("covar")) == 0) {
				CovarIndex[0] = atoi(InterList[k][i].Term1 + 5) - 1;
				Term1SNP = 0;
			}
			else {
				SNPindex[0] = FindInterItem(InterList[k][i].Term1);
				Term1SNP = 1;
			}

			if (strncmp("covar", InterList[k][i].Term2, strlen("covar")) == 0) {
				CovarIndex[1] = atoi(InterList[k][i].Term2 + 5) - 1;
				Term2SNP = 0;
			}
			else {
				SNPindex[1] = FindInterItem(InterList[k][i].Term2);
				Term2SNP = 1;
			}

			if (Term1SNP && Term2SNP) { // Epistasis, effect attributed to genetics
				for (j = 0; j < nSample; j++)
					GenoEff[k][j] += InterList[k][i].InterEff * InterItemMat[SNPindex[0]][j] * InterItemMat[SNPindex[1]][j];
			}
			if (Term1SNP && !Term2SNP) { // SNP x covariate interaction
				if (CovarIndex[1] < nCovar){
					for (j = 0; j < nSample; j++)
						CovarEff[k][j] += InterList[k][i].InterEff * InterItemMat[SNPindex[0]][j] * CovarMat[CovarIndex[1]][j];
				}
				else 
					printf("Covariate index our of bound! Input %d covatiates, index %d.\n", nCovar, CovarIndex[1]+1);
			}
			if (!Term1SNP && Term2SNP) { // SNP x covariate interaction
				if (CovarIndex[0] < nCovar)
					for (j = 0; j < nSample; j++)
						CovarEff[k][j] += InterList[k][i].InterEff * InterItemMat[SNPindex[1]][j] * CovarMat[CovarIndex[0]][j];
				else 
					printf("Covariate index our of bound! Input %d covatiates, index %d.\n", nCovar, CovarIndex[0]+1);
			}
			if (!Term1SNP && !Term2SNP) { // covariate x covariate interaction
				if ((CovarIndex[0] < nCovar) && (CovarIndex[1] < nCovar)) {
					for (j = 0; j < nSample; j++)
						CovarEff[k][j] += InterList[k][i].InterEff * CovarMat[CovarIndex[0]][j] * CovarMat[CovarIndex[1]][j];
				}
				else 
					printf("Covariate index our of bound! Input %d covatiates, index %d and %d.\n", nCovar, CovarIndex[0]+1, CovarIndex[1]+1);
			}
		}
	}
}


void GetEnvEff() {
	int i, j, k;
	long int l;
	gsl_matrix * tmpEnvEff = gsl_matrix_calloc(nTrait, nSample);
	for (l = 0; l < nSample; l++) {
		for (k = 0; k < nTrait; k++)
			gsl_matrix_set(tmpEnvEff, k, l, gsl_ran_gaussian(r, 1.0));
	}

	if (!statusTrait) {
		gsl_matrix * tmpCorrEnvEff = gsl_matrix_calloc(nTrait, nSample);
		gsl_matrix * tmpLenv = gsl_matrix_calloc(nTrait, nTrait);
		gsl_matrix_set_zero(tmpL);
		for (i = 0; i < nTrait; i++) {
			for (j = 0; j <= i; j++)
				gsl_matrix_set(tmpLenv, i, j, gsl_matrix_get(SigmaTrait, i, j));
		}
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmpLenv, tmpEnvEff, 0.0, tmpCorrEnvEff);
		gsl_matrix_memcpy(tmpEnvEff, tmpCorrEnvEff);
		gsl_matrix_free(tmpCorrEnvEff);
		gsl_matrix_free(tmpLenv);
	}
	else if (statusTrait && nTrait > 1) {
		printf("Trait correlation matrix is not positive definite, assuming enveromental effect independent.\n");
	}

	for (l = 0; l < nSample; l++) {
		for (k = 0; k < nTrait; k++)
			EnvEff[k][l] = gsl_matrix_get(tmpEnvEff, k, l);
	}
	gsl_matrix_free(tmpEnvEff);
}


double IsCausal(char SNP[50]) {
	memset(CausalFlag, 0, sizeof(int)*nMaxTrait);
	int i, k, flag;
	flag = 0;
	if (PolyFlag) {
		prob = gsl_rng_uniform(r);
		if (prob < pCausal[0]) {
			CausalFlag[0] = 1;
			for (k = 1; k < nTrait; k++) {
				prob = gsl_rng_uniform(r);
				CausalFlag[k] = ((prob < Pleio[k]) ? 1 : 0);
			}
		}
		else {
			CausalFlag[0] = 0;
			for (k = 1; k < nTrait; k++) {
				prob = gsl_rng_uniform(r);
				CausalFlag[k] = ((prob < pCausal[k]) ? 1 : 0);
			}
		}
	}
	else {
		for (k = 0; k < nTrait; k++) {		
			for (i = 0; i < nCausal[k]; i++) {
				if (!strcmp(SNP, CausalList[k][i])) {
					CausalFlag[k] = 1;
					break;
				}
			}
		}
	}
	for (k = 0; k < nTrait; k++)
		flag += CausalFlag[k];
	flag = ((flag > 0) ? 1.0 : 0.0);
	return(flag);
}


long int FindSNPinRef(char SNP[50], int chr) {
	int i; 
	for (i = 0; i < SNPct[chr]; i++) {
		// printf("i = %d, SNP = %s, total SNP on chr = %d\n", i, RefSNP[chr][i].SNP, SNPct[chr]);
		if (strcmp(SNP, RefSNP[chr][i].SNP) == 0)
			return(i);
	}
	return(-1);
}


double GetMAF(int popIndex) {
	double freq;
	long int i, j;
	double * tmpGeno;
	tmpGeno = malloc(sizeof(double) * nSamplePerPop[popIndex]);
	j = 0;
	for (i = 0; i < nSamplePerPop[popIndex]; i++) {
		if (PopMatTmp[popIndex][0][i] != 3.0) {
			tmpGeno[j] = PopMatTmp[popIndex][0][i];
			j++;
		}
	}
	freq = gsl_stats_mean(tmpGeno, 1, j)/2.0;
	freq = (freq > 0.5) ? (1-freq) : freq;
	free(tmpGeno);
	return(freq);
}


double GetBeta(double BaseBeta, int popIndex) {
	if (CausalEffFlag == 1)
		return(BaseBeta);
	else {
		if (CausalMAF[popIndex] == 0.0)
			return(0.0);
		else if (CausalLDscore == 0.0)
			CausalLDscore = 0.0001; // will run into problem if divide 0; therefore use 0.0001 as minimum
		double fMAF = pow(CausalMAF[popIndex]*(1-CausalMAF[popIndex]), a);
		double fLD = pow(CausalLDscore, b);
		double fScore = pow(CausalAnnot, c);
		double sigma2 = fMAF*fLD*fScore;
		return(BaseBeta * sqrt(sigma2));
	}
}


void AnalyzeSNP(int chr, char SNP[50]) {
	double prob;
	long int SNPindex;
	int i, j, n, k, popIndex;
	char tmp[50];

	char tmpMAF[500];
	char tmpBeta[500];
	char tmpBaseBeta[500];
	char tmpOutCausal[1000];
	
	FILE *OutFileCausal;
	memset(CausalBeta, 0.0, sizeof(double)*nMaxPop*nMaxTrait);
	memset(CausalMAF, 0.0, sizeof(double)*nMaxPop);

	if (CausalEffFlag != 1) {
		if ( !NoRefFlag )
			SNPindex = FindSNPinRef(SNP, chr);
		else 
			SNPindex = -1;

	    if (SNPindex == -1) {
			CausalLDscore = 0.25;
			CausalAnnot = 1.0;
		}
		else {
	    	CausalLDscore = RefSNP[chr][SNPindex].AfricaLDscore;
	    	CausalAnnot = RefSNP[chr][SNPindex].DHS + 1;
	    }
	    memset(tmpMAF, '\0', sizeof tmpMAF);
	    for (n = 0; n < nPop; n++) {
	    	CausalMAF[n] = GetMAF(n);
	    	sprintf(tmp, "%f,", CausalMAF[n]);
			strcat(tmpMAF, tmp);
	    }
	    tmpMAF[strlen(tmpMAF)-1] = '\0';

	    prob = gsl_rng_uniform(r);
	    if (prob < ProbComp[0])
	    	k = 0;
	    else {
			for (j = 1; j < nComp; j++) {
				if ((prob > ProbComp[j-1]) && (prob < ProbComp[j])) {
					k = j;
					break;
				}
			}
	    }
		BaseBetaGet(wComp[k]);

		for (k = 0; k < nTrait; k++) {
			if (CausalFlag[k]) {
				memset(tmpBeta, '\0', sizeof tmpBeta);
			    memset(tmpBaseBeta, '\0', sizeof tmpBaseBeta);
				for (i = 0; i < nPop; i++) {
					CausalBeta[k][i] = GetBeta(BaseBeta[i][k], i);
					sprintf(tmp, "%f,", CausalBeta[k][i]);
					strcat(tmpBeta, tmp);
					sprintf(tmp, "%f,", BaseBeta[i][k]);
					strcat(tmpBaseBeta, tmp);
				}
				tmpBeta[strlen(tmpBeta)-1] = '\0';
		    	tmpBaseBeta[strlen(tmpBaseBeta)-1] = '\0';
		    	strcpy(tmpOutCausal, OutCausal);
		    	sprintf(tmpBuff, "%d", (int)k+1);
				OutFileCausal = fopen(strcat(tmpOutCausal, tmpBuff),"a");
			    fprintf(OutFileCausal, "%s\t%s\t%s\t%lf\t%s\n", SNP, tmpBaseBeta, tmpMAF, CausalLDscore, tmpBeta);
				fclose(OutFileCausal);
			}
		    else {
		    	for (i = 0; i < nPop; i++)
					CausalBeta[k][i] = 0.0;
		    }
		}
	    for (i = 0; i < nSample; i++) {
	    	popIndex = PopIndicator[i];
	    	for (k = 0; k < nTrait; k++)
	    		GenoEff[k][i] += CausalBeta[k][popIndex]*GenoMat[i];
	    }
	}

	else {
		for (k = 0; k < nTrait; k++) {
			SNPindex = findSNPinCausalList(SNP, k);
			if (SNPindex != -1) {
				for (i = 0; i < nSample; i++) {
		    		if (GenoMat[i] == 2)
		    			GenoEff[k][i] += CausalEff[k][SNPindex][0]*GenoMat[i];
		    		if (GenoMat[i] == 1)
		    			GenoEff[k][i] += CausalEff[k][SNPindex][1]*GenoMat[i];
				}
			}
		}
	}
}

