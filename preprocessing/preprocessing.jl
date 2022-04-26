using Printf

include("../utils/parameter_parsing.jl")
include("utils.jl")


"""Implements the full sequence of pre-processing steps
"""
function preprocessing_pipeline(filepaths)

    # TODO what to do about this... already converted 1KG+HGDP to hapmap format (though some snps appear to be missing)
    # TODO the ID field in the VCFs is missing - need to fill this (add a code check for this and throw error if not present)
    @info "Filtering SNPs"
    # extract_variants(filepaths.vcftools, filepaths.vcf_input_raw, filepaths.vcf_input_processed_prefix, filepaths.vcf_input_processed, filepaths.variant_list)
    cp(filepaths.vcf_input_raw, filepaths.vcf_input_processed, force=true)
    
    @info "Creating genetic distance files"
    get_genetic_distances(filepaths.vcf_input_processed, filepaths.genetic_mapfile, filepaths.rsid_list, filepaths.genetic_distfile)

    @info "Creating mutation age files"
    get_mutation_ages(filepaths.vcf_input_processed, filepaths.mutation_mapfile, filepaths.rsid_list, filepaths.mutation_agefile)
    
    @info "Storing haplotype matrices"
    convert_vcf_to_hap(filepaths.vcf_input_processed, filepaths.hap1_matrix_output, filepaths.hap2_matrix_output)

    @info "Storing metadata files"
    convert_vcf_to_meta(filepaths.vcf_input_processed, filepaths.metadata_output)
    # TODO the popfile was not given for 1KG+HGDP
    cp(filepaths.popfile_raw, filepaths.popfile_processed, force=true)
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
