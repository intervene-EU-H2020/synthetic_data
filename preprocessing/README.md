# Data pre-processing

The pre-processing scripts prepare the reference data for generating synthetic datasets. There are options of using the provided dataset or pre-processing your own reference dataset.

## Standard pre-processing (default) 

Use this option if you want to jump straight into generating synthetic datasets without having to run pre-processing yourself. Using the `fetch` command, the container downloads data that we've already pre-processed for use in the synthetic data algorithm. Based on the 1000 Genomes Project and Human Genomic Diversity Project datasets, the data has been filtered to retain the HapMap3 variants. Please note that synthetic datasets contain the same set of variants as the pre-processed training dataset.

### How this data was created

Raw datasets for the 1KG+HGDP reference were downloaded from https://gnomad.broadinstitute.org/downloads. Genetic maps were obtained from https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_from_hapmap, and age of mutation from https://human.genome.dating/download/index.

The preprocessing pipeline is run using the `preprocess` command:

```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif preprocess data/config.yaml
```

## Custom pre-processing

Use this option if you want to use your own training dataset. You will need to edit the `filepaths` section of `config.yaml`. In particular, the pre-processing pipeline requires the following files as input:

- **`vcf_input_raw`**: Genotype data given in phased VCF format (following the convention of encoding the phased genotype as allele values separated by the "|" character)
- **`variant_list`**: File giving a list of variants to retain during pre-processing. This corresponds to the variants that will be included in the synthetic dataset
- `remove_list`: File giving a list of samples to remove during pre-preprocessing 
- **`rsid_list`**: File giving a mapping of variant ids (that appear in raw VCF input) to RSids
- **`genetic_mapfile`**: Genetic distance maps of basepair to centimorgan distances
- **`mutation_mapfile`**: Age of mutation for each variant
- `popfile_raw`: A population file

Note that for inputs marked in bold, you will need to provide *one file for each chromosome*. In the `config.yaml` file this is specified by using the `{chromosome}` wildcard in the filepath name.

Please see examples of how these inputs should be formatted by viewing the default dataset provided. Note that if you use the default `genetic_mapfile` and `mutation_mapfile`, the `preprocessing` pipeline will automatically interpolate values for the remaining variants.

After updating the `config.yaml` you can run the preprocessing pipeline using the `preprocess` command:
```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif preprocess data/config.yaml
```

Note that setting `chromosome` to `all` in the `config.yaml` will run the code for each chromosome sequentially. Alternatively you can run the preprocessing pipeline separately for each chromosome on different compute nodes.

