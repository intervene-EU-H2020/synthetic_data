using Random, Distributions, DataFrames, CSV, StatsBase, ProgressMeter

include("../../utils/parameter_parsing.jl")
include("write_output.jl")


"""Sample effective population size
"""
function sample_T(N, Ne)
    T_dist = Gamma(2, Ne/N)
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
- T: the coalescent age sampled for the segment

Note that the start and end variant positions are included in the segment
"""
function create_reference_table(metadata, start_haplotype, end_haplotype, batchsize)
    total_num_hap = 2*batchsize
    # to multi-thread the operation we create a dataframe for each haplotype, then concat these together
    ref_df_samples = Vector{DataFrame}(undef, total_num_hap)
    p = Progress(total_num_hap)
    Threads.@threads for hap in start_haplotype:end_haplotype
        happop = metadata.population_groups[hap]
        pos = 1
        # create segments until all variant positions are filled
        ref_df_hap = DataFrame(H=Int[], I=String[], S=Int[], E=Int[], P=String[], Q=String[], T=Float64[])
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
            push!(ref_df_hap, [hap, seghap, start_pos, end_pos, happop, segpop, T])
            pos = end_pos+1
        end
        ref_df_samples[hap-start_haplotype+1] = ref_df_hap
        next!(p)
    end
    # concat the arrays of dataframes
    ref_df = vcat(ref_df_samples...)
    return ref_df
end


"""Adjusts the batch size if not a divisor of the number of samples
"""
function get_adjusted_batchsize(batchsize, number_of_batches, batch_number, metadata)
    adjusted_batchsize = batchsize
    if batch_number == number_of_batches
        adjusted_batchsize = metadata.nsamples - batchsize*(number_of_batches-1)
    end
    return adjusted_batchsize
end


"""Returns the subset of ref_df containing only rows relevant to the current batch
"""
function get_start_end_haplotypes(batch_number, prev_batchsize, cur_batchsize)
    start_haplotype = ((batch_number-1)*prev_batchsize)*2+1
    end_haplotype = start_haplotype+cur_batchsize*2-1
    return start_haplotype, end_haplotype
end


"""Create the synthetic data for the specified chromosome
"""
function create_synthetic_genotype_for_chromosome(metadata)
    @info "Writing the population file"
    create_sample_list_file(metadata)

    number_of_batches = Int(ceil(metadata.nsamples/metadata.batchsize))
    batch_files = Vector{String}(undef, number_of_batches)

    total_samples_written = 0
    for batch_number in 1:number_of_batches
        adjusted_batchsize = get_adjusted_batchsize(metadata.batchsize, number_of_batches, batch_number, metadata)
        @info @sprintf("Creating the reference table for batch %i", batch_number)
        start_haplotype, end_haplotype = get_start_end_haplotypes(batch_number, metadata.batchsize, adjusted_batchsize)
        batch_ref_df = create_reference_table(metadata, start_haplotype, end_haplotype, adjusted_batchsize)
        @info @sprintf("Writing the plink file for batch %i", batch_number)
        batch_file = write_to_plink_batch(batch_ref_df, metadata.batchsize, adjusted_batchsize, batch_number, metadata)
        batch_files[batch_number] = batch_file[1:end-4]
        total_samples_written += adjusted_batchsize
    end

    @assert total_samples_written == metadata.nsamples # check all synthetic samples were written to the output

    @info "Merging all batch files"
    merge_batch_files(batch_files, metadata.outfile_prefix, metadata.plink, metadata.memory)
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
