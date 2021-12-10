# Pre-processing

The pre-processing scripts prepare the raw data files to be used as input by the pipeline to generate synthetic datasets.

## Standard (default) pre-processing

The data pre-processing procedure is rather time-consuming, so for users who want to quickly generate synthetic datasets, we provide an already pre-processed dataset:

- Based on 1000 Genomes dataset downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- Filtered to retain the variants listed in the files located at `data/external/variant_lists`. Please note that the synthetic dataset will contain the same set of variants as the pre-processed input dataset

### Usage

Download the raw datasets by running
```
preprocessing/data_download.sh data/inputs/raw
``` 
This downloads the data to the `data/inputs/raw` directory.

Configure the `config.yml` file and run the preprocessing pipeline

```
julia run_program.jl --config config.yml --preprocessing
```

## Custom pre-processing

There may be cases where you want to use your own input dataset. In this case, you can either pre-process your input dataset using the pipeline we provide or do the pre-processing yourself.

### If using the provided pipeline:
- See the `config.yml` file for inputs that need to be specified 
    - The program expects the genotype data input to be phased genotype data in VCF format, following the convention of encoding the phased genotype as allele values separated by the "|" character
    - Provide a variant list specifying which variants to retain in the data
    - Genetic distance maps can be downloaded using the `preprocessing/data_download.sh` script
- Run the preprocessing pipeline using

```
julia run_program.jl --config config.yml --preprocessing
```

### If pre-processing the data yourself

You will need to create the following files:
- Haplotype matrix files: the `convert_vcf_to_hap()` function from `preprocessing/utils.jl` shows how to convert a (phased) VCF file to the haplotype matrix input required by the synthetic data algorithm
- Genetic distance files: specifying conversion from basepair to centimorgan distances
- Metadata files: information about the genetic variants, to be stored with the synthetic dataset

There should be one set of inputs for each chromosome. Please see examples of the expected input data formats in the `data/inputs/processed` directory
