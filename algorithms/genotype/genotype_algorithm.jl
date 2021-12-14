using Random, Distributions, DataFrames, CSV, StatsBase, ProgressMeter

include("../../utils/parameter_parsing.jl")
include("write_output.jl")


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
function create_reference_table(metadata)
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
            end_pos = min(update_variant_position(pos, L, genetic_distances), nvariants)
            ref_df = push!(ref_df, [hap, seghap, start_pos, end_pos, happop, segpop])
            pos = end_pos+1
        end
    end
    return ref_df
end


"""Create the synthetic data for the specified chromosome
"""
function create_synthetic_genotype_for_chromosome(chromosome, superpopulation, options)
    # create the reference dataframe then store the data using the specified file format
    fp = parse_filepaths(options, chromosome, superpopulation)
    metadata = parse_genomic_metadata(options, superpopulation, fp)
    ref_df = create_reference_table(metadata)

    if metadata.outfile_type == "plink"
        write_to_plink(ref_df, metadata.batchsize, metadata)
    elseif metadata.outfile_type == "vcf"
        write_to_vcf(ref_df, metadata)
    else
        throw(error("Config error: outfile_type not supported"))
    end
end


"""Entry point to running the algorithm for generating a synthetic genotype dataset
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
