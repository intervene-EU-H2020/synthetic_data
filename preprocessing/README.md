# Pre-processing

The pre-processing scripts prepare the raw data files to be used as input by the pipeline to generate synthetic datasets. The data pre-processing procedure is rather time-consuming, so for users who want to quickly generate synthetic datasets, we provide an already pre-processed dataset (based on the 1000 Genomes dataset) in the `data/inputs/processed` directory. However, there may be cases where you want to pre-process your own input dataset - for example, if you would like to generate synthetic data for a different set of variants, or if you would like to use a larger input dataset. We provide specifications for how to prepare custom input datasets below.

## Standard (default) pre-processing

The data that we provide already pre-processed as a default option is prepared according to the following specifications:

- 1000 Genomes dataset for 22 chromosomes downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- Filtered to retain only the (HapMap) variants listed in the file `data/inputs/raw/variant_list.txt`. Please note that the synthetic dataset will contain the same set of variants as the input dataset
- All 3202 samples are retained, for 27 populations, corresponding to 5 superpopulations (AFR, AMR, EAS, EUR, SAS). For the synthetic dataset you can specify the number of samples, and custom population structure/admixture according to the given list of populations/superpopulations

### Usage

- The raw datasets can be downloaded by running `preprocessing/data_download.sh data/inputs/raw`, which will download the raw data to the `data/inputs/raw` directory
- The raw datasets can be preprocessed by running `preprocessing/main.jl --config config/preprocessing.yml`, where input parameters are configured in the `config/preprocessing.yml` file

## Custom pre-processing

If using your own dataset, you can either use the `preprocessing/main.jl` script to pre-process your dataset or pre-process the dataset yourself.

If using the `preprocessing/main.jl` script:
- The program expects the input to be phased genotype data in VCF format, following the VCF convention of encoding the phased genotype as allele values separated by the "|" character
- See the `config/preprocessing.yml` file for other inputs that need to be specified 

If not using the `preprocessing/main.jl` script, you will need to create the following inputs:
- Haplotype matrix files: the `convert_vcf_to_hap()` function in `preprocessing/utils.jl` shows how to convert a (phased) VCF file to the haplotype matrix input required by the synthetic data algorithm
- Genetic distance files: specifying conversion from basepair to centimorgan distances
- Metadata files: information about the genetic variants, to be stored with the synthetic dataset

There should be one set of inputs for each chromosome. Please see examples of the expected input data formats in the `data/inputs/processed` directory
