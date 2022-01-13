# Synthetic data generator for genotypes and phenotypes

Software program for generating large and realistic genotype+phenotype datasets. 

Insert image giving overview of pipeline here

## Quickstart

1. Setup Docker container ...
2. Create a synthetic dataset by running

```
julia run_program.jl --config config.yml --genotype --phenotype
```

## Detailed usage

### Training data

By default, the code uses the publicly available 1000 Genomes dataset as the training dataset. This data is downloaded with the container so you can proceed directly to the next steps for creating synthetic datasets. However, if you would like to use your own training datasets, you will need to follow the instructions for [preprocessing training data](preprocessing/README.md).

### Creating synthetic datasets 

The code implements several different pipelines, which can be run separately or together:
- Preprocessing of training data
- Creation of synthetic genotype and/or phenotype datasets
- Evaluation of synthetic data quality
- Optimisation of synthetic data algorithm parameters

There is one main command for running the program, which has a number of optional flags for indicating which pipeline/s you would like to run:
```
julia run_program.jl --config config.yml [--preprocessing] [--genotype] [--phenotype] [--evaluation] [--optimisation]
```

- `--config config.yml` gives the path to the configuration file
- `--preprocessing` is used for pre-processing training data
- `--genotype` is used for creating synthetic genotype data
- `--phenotype` is used for creating synthetic phenotype data
- `--evaluation` is used for evaluating synthetic data quality
- `--optimisation` is used for optimising synthetic data algorithm parameters

For example, if you wanted to create a synthetic genotype dataset and evaluate the data quality, you would run the command

```
julia run_program.jl --config config.yml --genotype --evaluation
```


## Code contributors

* Add contributors here

This work is part of the [INTERVENE project](https://www.interveneproject.eu/)

## Cite as

Link to publication here

## Acknowledgments

* Acknowledgments for data, code, etc.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
