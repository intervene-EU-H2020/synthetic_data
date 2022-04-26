# Data pre-processing

The pre-processing scripts prepare the training data for generating synthetic datasets. There are options of using the provided 1000 Genomes dataset or pre-processing your own training dataset.

## Standard pre-processing (default) 

Use this option if you want to jump straight into generating synthetic datasets without having to run pre-processing yourself. Using the `fetch` command, the container downloads data that we've already pre-processed for use in the synthetic data algorithm. Based on the 1000 Genomes Project and Human Genomic Diversity Project datasets, the data has been filtered to retain only the HapMap3 variants. Please note that synthetic datasets contain the same set of variants as the pre-processed training dataset.

### How this data was created

Raw datasets were downloaded by running 

```
preprocessing/download.sh data/inputs/raw/1KG+HGDP
``` 

This downloads the data to the `data/inputs/raw/1KG+HGDP`.

The preprocessing pipeline is run using

```
julia run_program.jl --config config.yml --preprocessing
```

Note that setting `chromosome` to `all` in the `config.yml` will run the code for each chromosome sequentially. Alternatively you can run the preprocessing pipeline separately for each chromosome on different compute nodes.

## Custom pre-processing

Use this option if you want to use your own training dataset. You will need to edit the `filepaths` section of `config.yml`. In particular, the pre-processing pipeline requires the following files as input:

- `vcf_input_raw`: Genotype data given in phased VCF format (following the convention of encoding the phased genotype as allele values separated by the "|" character)
- `variant_list`: `.txt` file giving a list of variants (one variant per line, no header) to retain during pre-processing. This corresponds to the variants that will be included in the synthetic dataset
- `genetic_mapfile`: Genetic distance maps of basepair to centimorgan distances. Columns are variant id, basepair distance, centimorgan distance (no headers)
- `mutation_mapfile`: Age of mutation for each variant
- `popfile_raw`: A population file similar to https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

Note that for all of the above inputs, except `popfile_raw`, you will need to provide *one file for each chromosome*. In the `config.yml` file this is specified by using the `{chromosome}` wildcard in the filepath name.

Please see examples of how these inputs should be formatted by viewing the default dataset provided.

After updating the `config.yml` you can run the preprocessing pipeline using
```
julia run_program.jl --config config.yml --preprocessing
```

