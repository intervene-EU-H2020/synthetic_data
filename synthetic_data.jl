using YAML, Printf


"""Parses the .yaml config files and executes the code for generating synthetic datasets
"""

general_options = YAML.load_file("config/general.yml")
data_options = YAML.load_file("config/data.yml")

if data_options["pipelines"]["generate_genotype_data"]
    @info "generating genotype data"
    # TODO generate genotype data - assume that datasets are stored according to the directory structure specified in the config file
end

if data_options["pipelines"]["generate_phenotype_data"]
    @info "generating phenotype data"
    # TODO generate phenotype data - assume that datasets are stored according to the directory structure specified in the config file
end

if data_options["pipelines"]["evaluation"]
    data_filepath = join([general_options["filepaths"]["output_directory"],general_options["filepaths"]["synthetic_data_prefix"]], "/")
    @info @sprintf("evaluating data quality for datasets available at %s", data_filepath)
    # TODO run evaluation - assume that datasets are stored according to the directory structure specified in the config file, and only run evaluation metrics requested in config
end

