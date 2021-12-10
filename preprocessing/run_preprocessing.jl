using Printf

include("../utils/parameter_parsing.jl")
include("preprocessing_utils.jl")


"""Run the full pre-processing pipeline
"""
function preprocessing_pipeline(filepaths)
    isdir(filepaths.output_dir) || mkdir(filepaths.output_dir)

    @info "Filtering SNPs"
    extract_variants(filepaths.vcftools, filepaths.vcf_input, filepaths.vcf_output, filepaths.variant_list)
    
    @info "Creating genetic distance files"
    get_genetic_distances(filepaths.vcf_output, filepaths.mapfile, filepaths.distfile)

    @info "Storing haplotype matrices"
    convert_vcf_to_hap(filepaths.vcf_output, filepaths.hap1_output, filepaths.hap2_output)

    @info "Storing metadata files"
    convert_vcf_to_meta(filepaths.vcf_output, filepaths.meta_output)

    # TODO delete files not needed
end


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
