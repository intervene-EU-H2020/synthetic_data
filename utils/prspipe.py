"""Converts geno+pheno+gwas files to format expected by prspipe

(note that this code is not formally part of the pipeline, 
but is for illustrative purposes to show how you can run this conversion)

Inputs: 
- genotype files for 22 chromosomes
- phenotype output 
- gwas evaluation output

Outputs:
- genotype output same as above, doesn't change
- phenotype and gwas output written to the specified output directory
"""

PHENO_PREFIX = "data/outputs/test/test.pheno1"
GWAS_PREFIX = "data/outputs/test/evaluation/test1.trait1.PHENO1.glm.linear"
OUTPUT_DIR = "data/outputs/test/prspipe"

import pandas as pd

pheno_df = pd.read_csv(PHENO_PREFIX, sep="\t")
pheno_df[0] = pheno_df['Sample'].apply(lambda x : x.split("_")[0])
pheno_df[1] = pheno_df['Sample'].apply(lambda x : x.split("_")[1])
pheno_df[2] = pheno_df['Phenotype']
pheno_df = pheno_df[[0,1,2]]
pheno_df.to_csv("{}/pheno.tsv".format(OUTPUT_DIR), header=None, sep='\t', index=None)

gwas_df = pd.read_csv(GWAS_PREFIX, sep="\t")
gwas_df = gwas_df[gwas_df.keys()[[0, 1, 2, 3, 4, 7, 8, 9, 13]].tolist()]
gwas_df = gwas_df.set_axis(["CHR", "POS", "SNP", "A2", "A1", "N", "BETA", "SE", "P"], axis=1, inplace=False)
gwas_df.to_csv("{}/test1.trait1.PHENO1.glm.linear".format(OUTPUT_DIR), sep='\t', index=None)
