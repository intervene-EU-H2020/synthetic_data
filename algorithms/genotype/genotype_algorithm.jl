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
function update_variant_position(cur_pos, L, genetic_distances, nvariants)
    var_pos = cur_pos
    dist = genetic_distances[var_pos]
    objective = dist+L
    while dist <= objective && var_pos<nvariants
        var_pos += 1
        dist = genetic_distances[var_pos]
    end
    return var_pos
end


"""Generates mutations in the synthetic haplotype by selecting SNP positions
"""
function get_mutations(N, Ne, μ, hap, nvariants)
    # sample number of mutations, based on coalescence time
    coal_T = sample_T(N, Ne)
    M_dist = Poisson(coal_T * μ)
    M = rand(M_dist)
    while M >= nvariants
        M = rand(M_dist)
    end
    # choose positions at random 
    # TODO could do this non-randomly using allele age
    variants = sample(1:nvariants, M, replace=false)
    mut_df = DataFrame(H=repeat([hap], length(variants)), P=variants)
    return mut_df
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
    mut_df = DataFrame(H=Int[], P=Int[])
    @showprogress for hap in 1:(2*metadata.nsamples)
        happop = metadata.population_groups[hap]
        pos = 1
        # create segments until all variant positions are filled
        while pos <= metadata.nvariants
            # sample a population group for the segment
            segpop = sample(collect(keys(metadata.population_weights[happop])), Weights(collect(values(metadata.population_weights[happop]))))
            # sample the segment distance
            T = sample_T(metadata.population_Ns[segpop], metadata.population_Nes[segpop])
            L = sample_L(T, metadata.population_rhos[segpop])
            # sample the haplotype to be copied
            seghap = rand(metadata.haplotypes[segpop])
            # update the start and end (variant) positions of the segment
            start_pos = pos
            end_pos = min(update_variant_position(pos, L, metadata.genetic_distances, metadata.nvariants), metadata.nvariants)
            push!(ref_df, [hap, seghap, start_pos, end_pos, happop, segpop])
            pos = end_pos+1
        end
        # sample mutations
        mutpop = sample(collect(keys(metadata.population_weights[happop])), Weights(collect(values(metadata.population_weights[happop]))))
        mut_df = vcat(mut_df, get_mutations(metadata.population_Ns[mutpop], metadata.population_Nes[mutpop], metadata.population_mus[mutpop], hap, metadata.nvariants))
    end
    return ref_df, mut_df
end


"""Create the synthetic data for the specified chromosome
"""
function create_synthetic_genotype_for_chromosome(metadata)
    @info "Creating the algorithm mapping table"
    ref_df, mut_df = create_reference_table(metadata)
    
    @info "Writing the population file"
    create_sample_list_file(metadata)
    
    if metadata.outfile_type == "plink"
        @info "Writing genotype output to PLINK"
        write_to_plink(ref_df, mut_df, metadata.batchsize, metadata)
    elseif metadata.outfile_type == "vcf"
        @info "Writing genotype output to VCF"
        write_to_vcf(ref_df, mut_df, metadata)
    else
        throw(error("Config error: outfile_type not supported"))
    end
end


"""Creates the population list file for synthetic samples
"""
function create_sample_list_file(metadata)
    outfile = @sprintf("%s.sample", metadata.outfile_prefix)
    open(outfile, "w") do io
        writedlm(io, metadata.population_groups[1:2:end])
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
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            metadata = parse_genomic_metadata(options, superpopulation, fp)
            create_synthetic_genotype_for_chromosome(metadata)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        metadata = parse_genomic_metadata(options, superpopulation, fp)
        create_synthetic_genotype_for_chromosome(metadata)
    end
end
