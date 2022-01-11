using Base: String

using Printf

include("utils.jl")


"""Struct specifying the filepaths for pre-processing inputs and outputs
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


"""Run the full pre-processing pipeline
"""
function run_preprocessing(filepaths)
    isdir(filepaths.output_dir) || mkdir(filepaths.output_dir)

    @info "Filtering SNPs"
    extract_variants(filepaths.vcftools, filepaths.vcf_input, filepaths.vcf_output, filepaths.variant_list)
    
    @info "Creating genetic distance files"
    get_genetic_distances(filepaths.vcf_output, filepaths.mapfile, filepaths.distfile)

    @info "Storing haplotype matrices"
    convert_vcf_to_hap(filepaths.vcf_output, filepaths.hap1_output, filepaths.hap2_output)

    @info "Storing metadata files"
    convert_vcf_to_meta(filepaths.vcf_output, filepaths.meta_output)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "config"
            help = "path to the YAML configuration file"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    fp = Filepaths() # TODO read from configuration file
    run_preprocessing(fp)
end


main()
