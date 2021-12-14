using Printf

include("../utils/parameter_parsing.jl")
include("utils.jl")


"""Implements the full sequence of pre-processing steps
"""
function preprocessing_pipeline(filepaths)
    isdir(filepaths.output_dir) || mkdir(filepaths.output_dir)

    @info "Filtering SNPs"
    extract_variants(filepaths.vcftools, filepaths.vcf_input_raw, filepaths.vcf_input_processed_prefix, filepaths.vcf_input_processed, filepaths.variant_list)
    
    @info "Creating genetic distance files"
    get_genetic_distances(filepaths.vcf_input_processed, filepaths.genetic_mapfile, filepaths.genetic_distfile)

    @info "Storing haplotype matrices"
    convert_vcf_to_hap(filepaths.vcf_input_processed, filepaths.hap1_matrix_output, filepaths.hap2_matrix_output)

    @info "Storing metadata files"
    convert_vcf_to_meta(filepaths.vcf_input_processed, filepaths.metadata_output)
    cp(filepaths.popfile_raw, filepaths.popfile_processed)
end


"""Entry point to running the pre-processing pipeline
"""
function run_preprocessing(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            preprocessing_pipeline(fp)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        preprocessing_pipeline(fp)
    end
end
