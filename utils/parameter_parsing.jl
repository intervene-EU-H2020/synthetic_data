"""Utility functions for parsing the options from the config file
"""

using Mmap, Printf

"""Struct specifying the filepaths for inputs and outputs
"""
struct Filepaths
    vcf_input_raw::String
    vcf_input_processed_prefix::String
    vcf_input_processed::String
    variant_list::String
    genetic_mapfile::String
    genetic_distfile::String
    hap1_matrix_output::String
    hap2_matrix_output::String
    metadata_output::String
    popfile_raw::String
    popfile_processed::String
    synthetic_data_prefix::String
    evaluation_output::String
    optimisation_output::String
    reference_dir::String
    vcftools::String
    plink::String
    plink2::String
    king::String
    mapthin::String
end


"""Struct specifying the metadata for constructing synthetic genotype data
"""
mutable struct GenomicMetadata
    nsamples::Integer
    nvariants::Integer
    H1::Matrix
    H2::Matrix
    fixed_fields::Vector # [list of metadata strings for each variant]
    haplotypes::Dict # {pop : [list of ids for haplotypes in the reference set]}
    index_map::Dict # {sample_id : sample_index}
    population_groups::Vector # [list of population groups for haplotypes in the synthetic set]
    population_weights::Dict # {popgroup: {pop : frac}}
    population_Ns::Dict # {pop : N}
    population_Nes::Dict # {pop : Ne}
    population_rhos::Dict # {pop : rho}
    genetic_distances::Vector # [list of genetic distances (in centimorgans) at each variant position]
    outfile_type::String
    outfile_prefix::String
    batchsize::Integer
    plink::String
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
    return superpopulation
end


function format_filepath(filepath, chromosome, superpopulation, chr_required)
    if chr_required && !occursin("{chromosome}", filepath)
        throw(error(@sprintf("Config error: no chromosome wildcard {chromosome} was specified for the filepath %s", filepath)))
    end

    formatted_filepath = replace(replace(filepath, "{superpopulation}" => superpopulation), "{chromosome}" => chromosome)
    return formatted_filepath
end


function parse_filepaths(options, chromosome, superpopulation)
    vcf_input_raw = format_filepath(options["filepaths"]["vcf_input_raw"], chromosome, superpopulation, true)
    vcf_input_processed_prefix = format_filepath(options["filepaths"]["vcf_input_processed_prefix"], chromosome, superpopulation, true)  
    vcf_input_processed = format_filepath(options["filepaths"]["vcf_input_processed"], chromosome, superpopulation, true)  
    variant_list = format_filepath(options["filepaths"]["variant_list"], chromosome, superpopulation, true)
    genetic_mapfile = format_filepath(options["filepaths"]["genetic_mapfile"], chromosome, superpopulation, true)
    genetic_distfile = format_filepath(options["filepaths"]["genetic_distfile"], chromosome, superpopulation, true)
    hap1_matrix_output = format_filepath(options["filepaths"]["hap1_matrix"], chromosome, superpopulation, true)
    hap2_matrix_output = format_filepath(options["filepaths"]["hap2_matrix"], chromosome, superpopulation, true)
    metadata_output = format_filepath(options["filepaths"]["vcf_metadata"], chromosome, superpopulation, true)
    popfile_raw = format_filepath(options["filepaths"]["popfile_raw"], chromosome, superpopulation, false)
    popfile_processed = format_filepath(options["filepaths"]["popfile_processed"], chromosome, superpopulation, false)
    synthetic_data_prefix = format_filepath(options["filepaths"]["synthetic_data_prefix"], chromosome, superpopulation, true)
    evaluation_output = format_filepath(options["filepaths"]["evaluation_output"], chromosome, superpopulation, true)
    optimisation_output = format_filepath(options["filepaths"]["optimisation_output"], chromosome, superpopulation, true)
    reference_dir = format_filepath(options["filepaths"]["reference_dir"], chromosome, superpopulation, false)

    vcftools = format_filepath(options["software_paths"]["vcftools"], chromosome, superpopulation, false)
    plink = format_filepath(options["software_paths"]["plink"], chromosome, superpopulation, false)
    plink2 = format_filepath(options["software_paths"]["plink2"], chromosome, superpopulation, false)
    king = format_filepath(options["software_paths"]["king"], chromosome, superpopulation, false)
    mapthin = format_filepath(options["software_paths"]["mapthin"], chromosome, superpopulation, false)

    return Filepaths(vcf_input_raw, vcf_input_processed_prefix, vcf_input_processed, variant_list, genetic_mapfile, genetic_distfile, hap1_matrix_output, hap2_matrix_output, metadata_output, popfile_raw, popfile_processed, synthetic_data_prefix, evaluation_output, optimisation_output, reference_dir, vcftools, plink, plink2, king, mapthin)
end


function open_hapfile(filepath)
    s = open(filepath)
    m = read(s, Int)
    n = read(s, Int)
    H = Mmap.mmap(s, Matrix{Int8}, (n, m))
    return H
end


function get_haplotype_list(popfile, poplist)
    pop_df = DataFrame(CSV.File(popfile))
    haplist = Dict{String,Vector{String}}()
    for pop in Set(pop_df.Superpopulation)
        if pop ∉ poplist
            throw(error(@sprintf("Config error: Population %s in %s is not supported", pop, popfile)))
        end
        haplist[pop] = pop_df[pop_df.Superpopulation .== pop, :].SampleID
    end
    return haplist
end


function get_index_map(popfile)
    pop_df = DataFrame(CSV.File(popfile))
    index_map = Dict(zip(pop_df.SampleID, 1:length(pop_df.SampleID)))
    return index_map
end


function get_population_sizes(haplotypes, poplist)
    population_sizes = Dict{String, Integer}()
    for pop ∈ poplist
        population_sizes[pop] = haskey(haplotypes, pop) ? length(haplotypes) : 0
    end
    return population_sizes
end


function get_genetic_distances(distfile)
    dist_df = DataFrame(CSV.File(distfile))
    sort!(dist_df, [:Distance])
    variant_dist = dist_df.Distance
    return variant_dist
end


function get_fixed_fields(metadata)
    fixed_fields = readlines(metadata)
    return fixed_fields
end


function get_population_structure(superpopulation, options, poplist)
    population_groups = []
    population_weights = Dict()
    
    use_default = options["genotype_data"]["samples"]["use_default"]
    
    if use_default
        @info "Using default population structure"
        nsamples = options["genotype_data"]["samples"]["default"]["nsamples"]
        if superpopulation == "none"
            population_groups = repeat(poplist, Int(ceil(nsamples/length(poplist))))[1:nsamples]
            for pop in poplist
                population_weights[pop] = Dict(pop=>100/length(poplist))
            end
        else
            population_groups = vcat(population_groups, repeat([superpopulation],nsamples*2))
            population_weights[superpopulation] = Dict(superpopulation=>100)
        end
    else
        custom_populations = options["genotype_data"]["samples"]["custom"]
        for pop in custom_populations
            population_groups = vcat(population_groups, repeat([pop["id"]],pop["nsamples"]*2))
            population_weights[pop["id"]] = merge(pop["populations"]...)
            if sum(values(population_weights[pop["id"]])) != 100
                throw(error("Config error: population weights do not sum to 100"))
            end
        end
    end
    
    nsamples = Int(length(population_groups)/2)
    
    return nsamples, population_groups, population_weights    
end


function parse_genomic_metadata(options, superpopulation, filepaths)

    poplist = ["AFR", "AMR", "EAS", "EUR", "SAS"]

    nsamples, population_groups, population_weights =  get_population_structure(superpopulation, options, poplist)
    H1 = open_hapfile(filepaths.hap1_matrix_output)
    H2 = open_hapfile(filepaths.hap2_matrix_output)
    fixed_fields = get_fixed_fields(filepaths.metadata_output)
    haplotypes = get_haplotype_list(filepaths.popfile_processed, poplist)
    index_map = get_index_map(filepaths.popfile_processed)
    population_N = get_population_sizes(haplotypes, poplist) 
    population_Nes = Dict(pop=>options["genotype_data"]["Ne"][pop] for pop in poplist)
    population_rhos = Dict(pop=>options["genotype_data"]["rho"][pop] for pop in poplist)
    genetic_distances = get_genetic_distances(filepaths.genetic_distfile)
    outfile_type = options["genotype_data"]["filetype"]
    outfile_prefix = filepaths.synthetic_data_prefix
    batchsize = outfile_type=="plink" ? options["genotype_data"]["batchsize"] : -1
    plink = filepaths.plink

    nvariants = length(genetic_distances)

    return GenomicMetadata(nsamples, nvariants, H1, H2, fixed_fields, haplotypes, index_map, population_groups, population_weights, population_N, population_Nes, population_rhos, genetic_distances, outfile_type, outfile_prefix, batchsize, plink)
end
