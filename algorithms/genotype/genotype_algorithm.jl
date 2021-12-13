using Random, Distributions, DataFrames, CSV, StatsBase, ProgressMeter

include("../../utils/parameter_parsing.jl")
include("write_output.jl")


"""Struct specifying the metadata for constructing synthetic data
"""
struct GenomicMetadata
    nsamples::Integer
    nvariants::Integer
    haplotypes::Dict # {pop : [list of ids for haplotypes in the reference set]}
    population_groups::Vector # [list of population groups for haplotypes in the reference set]
    population_weights::Dict # {popgroup: {pop : frac}}
    population_Ns::Dict # {pop : N}
    population_Nes::Dict # {pop : Ne}
    population_rhos::Dict # {pop : rho}
    genetic_distances::Dict # {chr : [list of genetic distances (in centimorgans) at each variant position]}
    filetype::String
    filename::String
end


"""Sample effective population size
"""
function sample_T(N, Ne)
    T_dist = Exponential(Ne/N)
    T = rand(T_dist)
    return T
end


"""Sample segment length
"""
function sample_L(T, ρ)
    L_dist = Exponential(1/(2*T*ρ))
    L = rand(L_dist)
    return L
end


"""Updates the variant position in the chromosome by adding the genetic distance L to the current position
"""
function update_variant_position(cur_pos, L, genetic_distances)
    var_pos = cur_pos
    dist = genetic_distances[var_pos]
    while dist <= dist+L
        var_pos += 1
        dist = genetic_distances[var_pos]
    end
    return var_pos
end


"""Creates reference table for synthetic genotype construction

Each row of the reference table represents a segment to be copied
from the reference set to the synthetic genotype dataset.

The reference table has the following columns:
- H: the identifier for the haplotype (numbered from 1...2*nsamples)
- I: the identifier for the individual to be copied from the reference haplotype set
- S: the index of the starting variant position for the segment (numbered from 1...nvariants)
- E: the index of the ending variant position for the segment (numbered from 1...nvariants)
- P: the population group for the synthetic haplotype 
- Q: the population group for the segment

Note that the start and end variant positions are included in the segment
"""
function create_reference_table(chromosome, metadata)
    ref_df = DataFrame(H=Int[], I=String[], S=Int[], E=Int[], P=String[], Q=String[])
    @showprogress for hap in 1:(2*nsamples)
        happop = metadata.population_groups[haplotype]
        pos = 1
        # create segments until all variant positions are filled
        while pos <= nvariants
            # sample a population group for the segment
            segpop = sample(keys(metadata.population_weights[happop]), Weights(values(metadata.population_weights[happop])))
            # sample the segment distance
            T = sample_T(metadata.population_Ns[segpop], metadata.population_Nes[segpop])
            L = sample_L(T, metadata.population_rhos[segpop])
            # sample the haplotype to be copied
            seghap = rand(metadata.haplotypes[segpop])
            # update the start and end (variant) positions of the segment
            start_pos = pos
            end_pos = min(update_variant_position(pos, L, genetic_distances[chromosome]), nvariants)
            ref_df = push!(ref_df, [hap, seghap, start_pos, end_pos, happop, segpop])
            pos = end_pos+1
        end
    end
    return ref_df
end


"""Create the metadata struct containing the information for constructing synthetic genotypes
"""
function create_metadata(options, filepaths)
    # TODO - create metadata struct and also do error checking on inputs
    nsamples = 
    nvariants = 
    return GenomicMetadata(nsamples, ...)
end


"""Create the synthetic data for the specified chromosome
"""
function create_synthetic_genotype_for_chromosome(chromosome, superpopulation, options)
    # create the reference dataframe then store the data using the specified file format
    fp = parse_filepaths(options, chromosome, superpopulation)
    metadata = create_metadata(options, fp)
    ref_df = create_reference_table(chromosome_i, metadata)

    if metadata.filetype == "plink"
        write_to_plink(ref_df, metadata) # TODO
    elseif metadata.filetype == "vcf"
        write_to_vcf(ref_df, metadata) # TODO
    else
        throw(error("Config error: filetype not supported"))
    end
end


"""Entry point to running the pre-processing pipeline
"""
function create_synthetic_genotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)

    # synthetic genotype files are generated chromosome-by-chromosome
    if chromosome == "all"
        for chromosome_i in 1:22
            create_synthetic_genotype_for_chromosome(chromosome_i, superpopulation, options)
        end
    else
        create_synthetic_genotype_for_chromosome(chromosome, superpopulation, options)
    end
end
