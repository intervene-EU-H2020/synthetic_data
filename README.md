# Synthetic data generation and evaluation

- Efficiently generate large-scale, diverse and realistic datasets for genotypes and phenotypes
- Easily analyse data quality with an extensive workflow for evaluating synthetic data reliability and generalisability

TODO architecture diagram

TODO add citation

## Quickstart

This quickstart tutorial will show you the simplest approach for generating and evaluating synthetic datasets.

1. **Download one of the software containers by choosing one of the [download links](#preface-on-containerisation). The rest of these instructions will be based on the Singularity container, but the same idea also applies for Docker, or running the software without a container.**
2. **Choose a location where you want to generate synthetic data, and create a `data` directory containing a copy of the `config.yaml` file downloaded from this repository and a `containers` directory containing the downloaded container file (your container version may be different, so you should update this in any subsequent commands):**

```
.
├── data
    └── config.yaml
└── containers
    └── synthetic-data-v1.0.0.sif
```

The `config.yaml` file contains the parameter values used for synthetic data generation. In this tutorial we will use the default configuration.

3. **From the directory you selected in the previous step, run the `init` command to complete the setup of software dependencies:**

```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif init
```

By default this will store Julia package files in the `~/.julia` directory. If you don't want the files stored at this location, add `--env JULIA_DEPOT_PATH={path-to-your-preffered-location}` to the above command (and all subsequent singularity commands). Also notice that the above command binds the `data` directory you created in the previous step. This directory is used by the synthetic data generation program to access input reference files and store output data files.

4. **The first time using the container, you need to fetch the reference dataset using the `fetch` command:**

```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif fetch
```

This will download the reference files used as input to the synthetic data generation program to the `data/inputs/processed` directory (see the [preprocessing/README.md](preprocessing/README.md) file for details about how this data was created).

5. **Now generate a synthetic dataset using the `generate_geno` and `generate_pheno` commands:**

```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif generate_geno 8 --config data/config.yaml
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif generate_pheno --config data/config.yaml
```

The number given before the `--config` flag in the `generate_geno` command in the number of computing threads you want the software to use. Running these commands should generate a small synthetic dataset in the `data/outputs/test` directory. 

Now it would be useful to evaluate the quality of this data.

6. **Evaluate synthetic data quality using the `validate` command:**

```
singularity exec --bind data/:/data/ containers/synthetic-data-v1.0.0.sif validate --config data/config.yaml
```

This will store visualisations in the `data/outputs/test/evaluation` directory and print results for quantitative metrics to stdout. See our manuscript for details about the evaluation workflow.

Now that you understand the basics, you can read about [how to customise your synthetic datasets](#customising-your-synthetic-datasets). You may be also be interested in [how to generate very large datasets](#large-scale-synthetic-data-generation).


## Preface on containerisation

For ease of portability and reproducibility, we've made this software available as both Docker and Singularity containers. These containerisation systems manage all the software dependencies, to make it easier for you to get started with generating synthetic datasets. 

You can download your choice of container from one of the following links:
- Docker: TODO add url
- Singularity: TODO add url

Alternatively, you can run this software without a container by manually installing the software dependencies. If you prefer to use this approach, you can view the required dependencies in the `Dockerfile` of this repository. 

## Customising your synthetic datasets

You can customize the specifications for synthetic data generation by modifying the `config.yml` file.

TODO

## Large scale synthetic data generation

In this section we provide instructions for utilising a SLURM-based HPC cluster to efficiently generate very large datasets.

TODO