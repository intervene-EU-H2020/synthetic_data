"""Utility functions for parsing the options from the config file
"""

using Mmap, Printf

"""Struct specifying the filepaths for inputs and outputs
"""
mutable struct Filepaths
    vcf_input_raw::String
    vcf_input_processed_prefix::String
    vcf_input_processed::String
    variant_list::String
    remove_list::String
    rsid_list::String
    genetic_mapfile::String
    genetic_distfile::String
    mutation_mapfile::String
    mutation_agefile::String
    hap1_matrix_output::String
    hap2_matrix_output::String
    metadata_output::String
    popfile_raw::String
    popfile_processed::String
    synthetic_data_prefix::String
    evaluation_output::String
    optimisation_output::String
    reference_dir::String
    prspipe_dir::String
    phenotype_causal_list::String
    phenotype_sample_list::String
    phenotype_reference::String
    vcftools::String
    plink::String
    plink2::String
    king::String
    mapthin::String
    phenoalg::String
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
    mutation_ages::Vector # [list of mutation ages (in years) at each variant position]
    outfile_prefix::String
    batchsize::Integer
    plink::String
    memory::Integer
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
    elseif superpopulation ∈ ["AFR", "AMR", "EAS", "EUR", "CSA", "MID"]
        @info @sprintf("Using the superpopulation %s", superpopulation)
    else
        throw(error("Config error: superpopulation value must be either `none` or one of AFR, AMR, EAS, EUR, CSA or MID"))
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
    vcf_input_raw = format_filepath(options["filepaths"]["genotype"]["vcf_input_raw"], chromosome, superpopulation, true)
    vcf_input_processed = format_filepath(options["filepaths"]["genotype"]["vcf_input_processed"], chromosome, superpopulation, true) 
    vcf_input_processed_prefix = endswith(vcf_input_processed,".recode.vcf") ? chop(vcf_input_processed, tail=11) : chop(vcf_input_processed, tail=4)
    variant_list = format_filepath(options["filepaths"]["genotype"]["variant_list"], chromosome, superpopulation, true)
    remove_list = format_filepath(options["filepaths"]["genotype"]["remove_list"], chromosome, superpopulation, false)
    rsid_list = format_filepath(options["filepaths"]["genotype"]["rsid_list"], chromosome, superpopulation, true)
    genetic_mapfile = format_filepath(options["filepaths"]["genotype"]["genetic_mapfile"], chromosome, superpopulation, true)
    genetic_distfile = format_filepath(options["filepaths"]["genotype"]["genetic_distfile"], chromosome, superpopulation, true)
    mutation_mapfile= format_filepath(options["filepaths"]["genotype"]["mutation_mapfile"], chromosome, superpopulation, true)
    mutation_agefile= format_filepath(options["filepaths"]["genotype"]["mutation_agefile"], chromosome, superpopulation, true)
    hap1_matrix_output = format_filepath(options["filepaths"]["genotype"]["hap1_matrix"], chromosome, superpopulation, true)
    hap2_matrix_output = format_filepath(options["filepaths"]["genotype"]["hap2_matrix"], chromosome, superpopulation, true)
    metadata_output = format_filepath(options["filepaths"]["genotype"]["vcf_metadata"], chromosome, superpopulation, true)
    popfile_raw = format_filepath(options["filepaths"]["genotype"]["popfile_raw"], chromosome, superpopulation, false)
    popfile_processed = format_filepath(options["filepaths"]["genotype"]["popfile_processed"], chromosome, superpopulation, false)

    synthetic_data_prefix = format_filepath(string(options["filepaths"]["general"]["output_dir"],"/",options["filepaths"]["general"]["output_prefix"]), chromosome, superpopulation, true)
    evaluation_output = format_filepath(string(options["filepaths"]["general"]["output_dir"],"/evaluation/",options["filepaths"]["general"]["output_prefix"]), chromosome, superpopulation, true)
    optimisation_output = format_filepath(string(options["filepaths"]["general"]["output_dir"],"/optimisation/",options["filepaths"]["general"]["output_prefix"]), chromosome, superpopulation, true)
    reference_dir = format_filepath(string(options["filepaths"]["general"]["output_dir"],"/reference"), chromosome, superpopulation, false)
    prspipe_dir = format_filepath(string(options["filepaths"]["general"]["output_dir"],"/prspipe"), chromosome, superpopulation, false)
    
    # check format of output prefix
    if !endswith(synthetic_data_prefix, @sprintf("-%s",chromosome))
        throw(error("Config error: output_prefix must use the naming convention {prefix}-{chromosome} where {prefix} is replaced by your choice and {chromosome} is left as a wildcard"))
    end

    phenotype_causal_list = format_filepath(options["filepaths"]["phenotype"]["causal_list"], chromosome, superpopulation, false)
    phenotype_sample_list = @sprintf("%s.sample", synthetic_data_prefix)
    phenotype_reference = format_filepath(options["filepaths"]["phenotype"]["reference"], chromosome, superpopulation, false)
    
    vcftools = format_filepath(options["filepaths"]["software"]["vcftools"], chromosome, superpopulation, false)
    plink = format_filepath(options["filepaths"]["software"]["plink"], chromosome, superpopulation, false)
    plink2 = format_filepath(options["filepaths"]["software"]["plink2"], chromosome, superpopulation, false)
    king = format_filepath(options["filepaths"]["software"]["king"], chromosome, superpopulation, false)
    mapthin = format_filepath(options["filepaths"]["software"]["mapthin"], chromosome, superpopulation, false)
    phenoalg = format_filepath(options["filepaths"]["software"]["phenoalg"], chromosome, superpopulation, false)

    return Filepaths(vcf_input_raw, vcf_input_processed_prefix, vcf_input_processed, variant_list, remove_list, rsid_list, genetic_mapfile, genetic_distfile, mutation_mapfile, mutation_agefile, hap1_matrix_output, hap2_matrix_output, metadata_output, popfile_raw, popfile_processed, synthetic_data_prefix, evaluation_output, optimisation_output, reference_dir, prspipe_dir, phenotype_causal_list, phenotype_sample_list, phenotype_reference, vcftools, plink, plink2, king, mapthin, phenoalg)
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
        population_sizes[pop] = haskey(haplotypes, pop) ? length(haplotypes[pop]) : 0
    end
    return population_sizes
end


function get_genetic_distances(distfile)
    dist_df = DataFrame(CSV.File(distfile))
    sort!(dist_df, [:Distance])
    variant_dist = dist_df.Distance
    return variant_dist
end


function get_mutation_ages(agefile)
    age_df = DataFrame(CSV.File(agefile))
    variant_age = age_df.Age
    return variant_age
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
            # use all superpopulations in equal ratio
            population_groups = repeat(repeat(poplist,inner=2), Int(ceil(nsamples/length(poplist))))[1:nsamples*2]
            for pop in poplist
                population_weights[pop] = Dict(pop=>100)
            end
        else
            # use a single population group, specified in the config
            population_groups = vcat(population_groups, repeat([superpopulation],nsamples*2))
            population_weights[superpopulation] = Dict(superpopulation=>100)
        end
    else
        # use custom population structure, specified in the config
        custom_populations = options["genotype_data"]["samples"]["custom"]
        for pop in custom_populations
            population_groups = vcat(population_groups, repeat([pop["id"]],pop["nsamples"]*2))
            population_weights[pop["id"]] = merge([Dict(k => i[k] for k in keys(i)) for i in pop["populations"]]...)
            if sum(values(population_weights[pop["id"]])) != 100
                throw(error("Config error: population weights do not sum to 100"))
            end
        end
    end
    
    nsamples = Int(length(population_groups)/2)
    
    return nsamples, population_groups, population_weights    
end


function get_batchsize(nsamples, batchsize)
    if nsamples > batchsize
        return batchsize
    else
        return nsamples
    end
end


function parse_genomic_metadata(options, superpopulation, filepaths)

    poplist = ["AFR", "AMR", "EAS", "EUR", "CSA", "MID"]
    
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
    mutation_ages = get_mutation_ages(filepaths.mutation_agefile)
    outfile_prefix = filepaths.synthetic_data_prefix
    batchsize = get_batchsize(nsamples, options["global_parameters"]["batchsize"])
    plink = filepaths.plink
    memory = options["global_parameters"]["memory"]

    nvariants = length(genetic_distances)

    return GenomicMetadata(nsamples, nvariants, H1, H2, fixed_fields, haplotypes, index_map, population_groups, population_weights, population_N, population_Nes, population_rhos, genetic_distances, mutation_ages, outfile_prefix, batchsize, plink, memory)
end
