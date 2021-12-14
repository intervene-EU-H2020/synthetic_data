"""Writes genotype data to PLINK BED/BIM/FAM output

Output is written in batches then merged back together
"""
function write_to_plink(ref_df, batchsize, metadata)

    number_of_batches = Int(ceil(metadata.nsamplesmples/batchsize))

    batch_files = Vector{String}(undef, number_of_batches)
    
    for batch_number in 1:number_of_batches
        batch_ref_df = get_batch_ref_df(ref_df, batch_number, batchsize)
        batch_file = write_to_plink_batch(batch_ref_df, batchsize, batch_number, metadata)
        batch_files[batch_number] = batch_file

    outfile = @sprintf("%s.bed", metadata.outfile_prefix)
    merge_batch_files(batch_files, outfile, metadata.plink)
end


"""Returns the subset of ref_df containing only rows relevant to the current batch
"""
function get_batch_ref_df(ref_df, batch_number, batchsize)
    start_haplotype = ((batch_number-1)*batchsize)*2+1
    end_haplotype = start_haplotype+batchsize*2-1
    batch_ref_df = ref_df[(ref_df.H .<= end_haplotype).&(ref_df.H .>= start_haplotype), :]
    return batch_ref_df
end


"""Merge all the batch files together
"""
function merge_batch_files(batch_files, outfile, plink)

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
        run(`rm $file_prefix.bed`)
        run(`rm $file_prefix.bim`)
        run(`rm $file_prefix.fam`)
    end

    run(`rm $mergelist`)
end


"""Construct the data into the format for writing to plink output
"""
function get_genostr(batch_ref_df, batchsize, start_haplotype, metadata)

    I_hap = Dict(1=>Array{Int8}(undef, batchsize, metadata.nvariants), 2=>Array{Int8}(undef, batchsize, metadata.nvariants))
    
    @showprogress for genotype in 1:batchsize
        for hap in [1,2]
            # construct the synthetic haplotypes, using the coordinates from the reference dataframe
            I = Vector{Int8}(undef, metadata.nvariants)

            if hap == 1
                hap_df = sort!(batch_ref_df[batch_ref_df.H .== (start_haplotype + (genotype-1)*2), :])
            else
                hap_df = sort!(batch_ref_df[batch_ref_df.H .== (start_haplotype + (genotype-1)*2 + 1), :])
            end
            
            I_pos = 1
            for row in eachrow(hap_df)
                segment = metadata.H1[metadata.index_map[row.I], (row.S+1):row.E]
                length_of_segment = length(segment)
                I[I_pos:I_pos+length_of_segment-1] = segment
                I_pos += length_of_segment
            end
            
            I = vcat(I...)
            I_hap[hap][genotype,:] = I
        end
    end

    genostr = I1_batch + I2_batch
    return genostr
end


"""Writes the plink output for a single batch
"""
function write_to_plink_batch(batch_ref_df, batchsize, batch_number, metadata)

    properties = Dict("fid"=>[string("syn",x) for x in ((batch_number-1)*batchsize+1):((batch_number-1)*batchsize+batchsize)],
            "iid"=>[string("syn",x) for x in ((batch_number-1)*batchsize+1):((batch_number-1)*batchsize+batchsize)],
            "chromosome"=>[split(f, "\t")[1][4:end] for f in metadata.fixed_fields], 
            "sid"=>[split(f, "\t")[3] for f in metadata.fixed_fields],
            "bp_position"=>[split(f, "\t")[2] for f in metadata.fixed_fields], 
            "allele_1"=>[split(f, "\t")[5] for f in metadata.fixed_fields], 
            "allele_2"=>[split(f, "\t")[4] for f in metadata.fixed_fields])

    start_haplotype = ((batch_number-1)*batchsize)*2+1
    genostr = get_genostr(batch_ref_df, batchsize, start_haplotype, metadata)

    batch_file = @sprintf("%s_i.bed", metadata.outfile_prefix, (batch_number-1))
    bed_reader.to_bed(batch_file, genostr, properties=properties)

    return batch_file
end


"""Writes genotype data to VCF output

Warning: only use VCF format for smaller datasets, because it is much slower than PLINK
"""
function write_to_vcf(ref_df, metadata)

    outfile = @sprintf("%s.vcf", metadata.outfile_prefix)

    # the file format line is always required
    file_format = "##fileformat=VCFv4.1\n"

    # the first 9 columns are fixed
    fixed_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"

    # one column for each synthetic sample
    data_fields = string(join([string("syn",x) for x in 1:metadata.nsamples], "\t"), "\n") 

    # note that we skip metadata lines
    
    out_string = string(file_format, fixed_columns, data_fields)
    
    open(outfile, "w") do io

        # write the initial lines to the file
        write(io, out_string)
        
        # write the remaining lines, variant-by-variant
        @showprogress for snp_number in 0:metadata.nvariants-1
            
            batch_string = string(batch_string, metadata.fixed_fields[snp_number+1]) # append fixed fields at start of line
            snp_df = sort!(ref_df[(ref_df.S .<= snp_number).&(ref_df.E .> snp_number), :], :H)
            I1_map = [string(metadata.H1[metadata.index_map[i],snp_number+1]) for i in snp_df.I[1:2:end]]
            I2_map = [string(metadata.H2[metadata.index_map[i],snp_number+1]) for i in snp_df.I[1:2:end]]
            batch_string = string(batch_string, join(broadcast.(*, I1_map, "|", I2_map, "\t"))[1:end-1], "\n")

            write(io, out_string)
        
        end
    
    end
end

