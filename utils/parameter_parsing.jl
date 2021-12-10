"""Utility functions for parsing the options from the config file
"""


"""Struct specifying the filepaths for inputs and outputs
"""
struct Filepaths
    input_dir::String
    output_dir::String
    vcf_input::String
    vcf_output::String
    variant_list::String
    vcftools::String
    mapfile::String
    distfile::String
    hap1_output::String
    hap2_output::String
    meta_output::String
end


function parse_chromosome(options)
    chromosome = options["global_parameters"]["chromosome"]

    if chromosome == "all"
        @info "Processing all 22 chromosomes"
    elseif isa(chromosome, Int) && chromosome >= 1 && chromosome <= 22
        @info @sprintf("Processing chromosome %i", chromosome)
    else
        throw(error("Config error: chromosome value must be either `all` or an integer between 1 and 22"))
    end

    return chromosome
end


function parse_superpopulation(options)
    superpopulation = options["global_parameters"]["superpopulation"]

    if superpopulation == "none"
        @info "No superpopulation specified"
    elseif superpopulation ∈ ["AFR", "AMR", "EAS", "EUR", "SAS"]
        @info @sprintf("Using the superpopulation %s", superpopulation)
    else
        throw(error("Config error: superpopulation value must be either `none` or one of AFR, AMR, EAS, EUR or SAS"))
    end
end


format_filepath(filepath, chromosome, superpopulation) = replace(replace(filepath, "{superpopulation}" => superpopulation), "{chromosome}" => chromosome)


function parse_filepaths(options, chromosome, superpopulation)
    input_dir = format_filepath(options["filepaths"]["input_dir"], chromosome, superpopulation)
    output_dir = format_filepath(options["filepaths"]["output_dir"], chromosome, superpopulation)
    vcf_input = format_filepath(options["filepaths"]["vcf_input"], chromosome, superpopulation)
    vcf_output = format_filepath(options["filepaths"]["vcf_output"], chromosome, superpopulation)  
    variant_list = format_filepath(options["filepaths"]["variant_list"], chromosome, superpopulation)
    vcftools = format_filepath(options["software_paths"]["vcftools"], chromosome, superpopulation)
    mapfile = format_filepath(options["filepaths"]["mapfile"], chromosome, superpopulation)
    distfile = format_filepath(options["filepaths"]["distfile"], chromosome, superpopulation)
    hap1_output = format_filepath(options["filepaths"]["hap1_output"], chromosome, superpopulation)
    hap2_output = format_filepath(options["filepaths"]["hap2_output"], chromosome, superpopulation)
    meta_output = format_filepath(options["filepaths"]["meta_output"], chromosome, superpopulation)

    return Filepaths(input_dir, output_dir, vcf_input, vcf_output, variant_list, vcftools, mapfile, distfile, hap1_output, hap2_output, meta_output)
end