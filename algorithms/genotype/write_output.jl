using PyCall

@pyimport bed_reader


"""Writes genotype data to PLINK BED/BIM/FAM output

Output is written in batches then merged back together
"""
function write_to_plink(ref_df, mut_df, batchsize, metadata)

    number_of_batches = Int(ceil(metadata.nsamples/batchsize))

    batch_files = Vector{String}(undef, number_of_batches)
    
    for batch_number in 1:number_of_batches
        batch_ref_df, batch_mut_df = get_batch_dfs(ref_df, mut_df, batch_number, batchsize)
        batch_file = write_to_plink_batch(batch_ref_df, batch_mut_df, batchsize, batch_number, metadata)
        batch_files[batch_number] = batch_file[1:end-4]
    end

    merge_batch_files(batch_files, metadata.outfile_prefix, metadata.plink)
end


"""Returns the subset of ref_df containing only rows relevant to the current batch
"""
function get_batch_dfs(ref_df, mut_df, batch_number, batchsize)
    start_haplotype = ((batch_number-1)*batchsize)*2+1
    end_haplotype = start_haplotype+batchsize*2-1
    batch_ref_df = ref_df[(ref_df.H .<= end_haplotype).&(ref_df.H .>= start_haplotype), :]
    batch_mut_df = mut_df[(mut_df.H .<= end_haplotype).&(mut_df.H .>= start_haplotype), :]
    return batch_ref_df, batch_mut_df
end


"""Merge all the batch files together
"""
function merge_batch_files(batch_files, outfile, plink)
    if length(batch_files) > 1
        # merge together the batch files
        mergelist = @sprintf("%s_mergefile.txt", outfile)
        basefile = batch_files[1]
        io = open(mergelist, "w") do io
        for x in batch_files[2:end]
            println(io, x)
        end
        end
        
        run(`$plink --bfile $basefile --merge-list $mergelist --make-bed --out $outfile`)

        # remove the batch files
        for file_prefix in batch_files
            rm(string(file_prefix,".bed"))
            rm(string(file_prefix,".bim"))
            rm(string(file_prefix,".fam"))
        end
        rm(mergelist)
    else
        # no files to merge - reformat and remove single batch file
        file_prefix = batch_files[1]
        run(`$plink --bfile $file_prefix --make-bed --out $outfile`)
        rm(string(file_prefix,".bed"))
        rm(string(file_prefix,".bim"))
        rm(string(file_prefix,".fam"))
    end
end


"""Add mutations by only copying if T < age of mutation
"""
function add_mutations(segment, T, mutation_ages, start_idx)
    new_segment = Vector{Int8}(undef, length(segment))
    for pos in 1:length(segment)
        if T < mutation_ages[start_idx+pos-1]
            new_segment[pos] = segment[pos]
        else
            new_segment[pos] = 0
        end
    end
    return new_segment
end


"""Construct the data into the format for writing to plink output
"""
function get_genostr(batch_ref_df, batch_mut_df, batchsize, start_haplotype, metadata)
    I_hap = Dict(1=>Array{Int8}(undef, batchsize, metadata.nvariants), 2=>Array{Int8}(undef, batchsize, metadata.nvariants))

    p = Progress(batchsize)
    Threads.@threads for genotype in 1:batchsize
        for hap in [1,2]
            # construct the synthetic haplotypes, using the coordinates from the reference dataframe
            I = Vector{Int8}(undef, metadata.nvariants)

            if hap == 1
                true_hap = (start_haplotype + (genotype-1)*2)
            else
                true_hap = (start_haplotype + (genotype-1)*2 + 1)
            end
            
            hap_df = batch_ref_df[batch_ref_df.H .== true_hap, :]
            
            I_pos = 1
            segment_sums = 0
            for row in eachrow(hap_df)
                if hap == 1
                    segment = metadata.H1[metadata.index_map[row.I], row.S:row.E]
                    segment = add_mutations(segment, row.T, metadata.mutation_ages, I_pos)
                else
                    segment = metadata.H2[metadata.index_map[row.I], row.S:row.E]
                    segment = add_mutations(segment, row.T, metadata.mutation_ages, I_pos)
                end
                segment_sums += sum(segment)
                length_of_segment = length(segment)
                I[I_pos:I_pos+length_of_segment-1] = segment
                I_pos += length_of_segment
            end
            
            I = vcat(I...)
            @assert segment_sums == sum(I)
            
            # add mutations
            for mut_pos in batch_mut_df[batch_mut_df.H .== true_hap, :].P 
                I[mut_pos] = 1 - I[mut_pos]
            end

            I_hap[hap][genotype,:] = I
        end
        next!(p)
    end

    genostr = I_hap[1] + I_hap[2]
    @assert length(genostr) == metadata.nvariants*batchsize
    return genostr
end


"""Writes the plink output for a single batch, using the Python package bed_reader
"""
function write_to_plink_batch(batch_ref_df, batch_mut_df, batchsize, batch_number, metadata)
    properties = Dict("fid"=>[string("syn",x) for x in ((batch_number-1)*batchsize+1):((batch_number-1)*batchsize+batchsize)],
            "iid"=>[string("syn",x) for x in ((batch_number-1)*batchsize+1):((batch_number-1)*batchsize+batchsize)],
            "chromosome"=>[split(f, "\t")[1][5:end] for f in metadata.fixed_fields], 
            "sid"=>[split(f, "\t")[3] for f in metadata.fixed_fields],
            "bp_position"=>[split(f, "\t")[2] for f in metadata.fixed_fields], 
            "allele_1"=>[split(f, "\t")[5] for f in metadata.fixed_fields], 
            "allele_2"=>[split(f, "\t")[4] for f in metadata.fixed_fields])

    start_haplotype = ((batch_number-1)*batchsize)*2+1
    genostr = get_genostr(batch_ref_df, batch_mut_df, batchsize, start_haplotype, metadata)
    batch_file = @sprintf("%s_%i.bed", metadata.outfile_prefix, (batch_number-1))

    bed_reader.to_bed(batch_file, genostr, properties=properties)

    return batch_file
end


"""Writes genotype data to VCF output

Warning: only use VCF format for smaller datasets, because it is much slower than PLINK
"""
function write_to_vcf(ref_df, mut_df, metadata)

    outfile = @sprintf("%s.vcf", metadata.outfile_prefix)

    # the file format line is always required
    file_format = "##fileformat=VCFv4.1\n"

    # the first 9 columns are fixed
    fixed_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"

    # one column for each synthetic sample
    data_fields = string(join([string("syn",x) for x in 1:metadata.nsamples], "\t"), "\n") 

    # note that we don't include any metadata lines
    
    out_string = string(file_format, fixed_columns, data_fields)
    
    open(outfile, "w") do io
        
        # write the initial lines to the file
        write(io, out_string)
        
        # write the remaining lines, variant-by-variant
        @showprogress for snp_number in 1:metadata.nvariants-1
            out_string = replace(metadata.fixed_fields[snp_number], "\"" => "") # append fixed fields at start of line
            snp_df = sort!(ref_df[(ref_df.S .<= snp_number).&(ref_df.E .>= snp_number), :], :H)
            I1_map = [string(metadata.H1[metadata.index_map[i],snp_number]) for i in snp_df.I[1:2:end]]
            I2_map = [string(metadata.H2[metadata.index_map[i],snp_number]) for i in snp_df.I[2:2:end]]
            # TODO add mutations for VCF
            out_string = string(out_string, join(broadcast.(*, I1_map, "|", I2_map, "\t"))[1:end-1], "\n")
            
            write(io, out_string)
            
        end
    
    end
end
