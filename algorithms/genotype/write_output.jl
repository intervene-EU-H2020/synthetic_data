"""Writes genotype data to plink BED/BIM/FAM output
"""
function write_to_plink(ref_df, metadata)
end


"""Writes genotype data to VCF output
"""
function write_to_vcf(ref_df, metadata)
end

# TODO refactor functions below
function write_to_vcf(H1, H2, ind_map, ref_df, fixed_fields, n_snps, n_syn_samples, outfile, batchsize)
    file_format = "##fileformat=VCFv4.1\n" # the file format line is always required
    fixed_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" # the first 9 columns are fixed
    data_fields = string(join([string("syn",x) for x in 1:n_syn_samples], "\t"), "\n") # one column for each synthetic sample
    
    # TODO optional metadata lines
    
    batch_string = string(file_format, fixed_columns, data_fields)
    
    open(outfile, "w") do io
        
        for snp_number in 0:n_snps-1
            
            # write out the output after a specified number of batch iterations
            if snp_number%batchsize==0
                write(io, batch_string)
                batch_string = ""
            end
            
            batch_string = string(batch_string, fixed_fields[snp_number+1]) # append fixed fields at start of line
            snp_df = sort!(ref_df[(ref_df.S .<= snp_number).&(ref_df.E .> snp_number), :], :H)
            I1_map = [string(H1[ind_map[i],snp_number+1]) for i in snp_df.I[1:2:end]]
            I2_map = [string(H2[ind_map[i],snp_number+1]) for i in snp_df.I[1:2:end]]
            batch_string = string(batch_string, join(broadcast.(*, I1_map, "|", I2_map, "\t"))[1:end-1], "\n")
        
        end
        
        # write out the last output string
        if length(batch_string)>0
            write(io, batch_string)
        end
    
    end
    
end


function write_to_plink_batch(batch, batch_file, batchsize, H1, H2, ind_map, batch_ref_df, fixed_fields, n_snps)
    batch_number = batch - 1

    I1_batch = Array{Int8}(undef, batchsize, n_snps)
    I2_batch = Array{Int8}(undef, batchsize, n_snps)
    
    start_haplotype = (batch_number*batchsize)*2+1
    
    @showprogress for genotype in 1:batchsize
        # construct the synthetic haplotypes, using the coordinates from the reference dataframe
        I1 = Vector{Int8}(undef, n_snps)
        hap_df1 = sort!(batch_ref_df[batch_ref_df.H .== (start_haplotype + (genotype-1)*2), :])
        I1_pos = 1
        for row in eachrow(hap_df1)
            segment = H1[ind_map[row.I], (row.S+1):row.E]
            length_of_segment = length(segment)
            I1[I1_pos:I1_pos+length_of_segment-1] = segment
            I1_pos += length_of_segment
        end
        I1 = vcat(I1...)
        I1_batch[genotype,:] = I1
    
        I2 = Vector{Int8}(undef, n_snps)
        hap_df2 = sort!(batch_ref_df[batch_ref_df.H .== (start_haplotype + (genotype-1)*2 + 1), :])
        I2_pos = 1
        for row in eachrow(hap_df2)
            segment = H2[ind_map[row.I], (row.S+1):row.E]
            length_of_segment = length(segment)
            I2[I2_pos:I2_pos+length_of_segment-1] = segment
            I2_pos += length_of_segment
        end
        I2 = vcat(I2...)
        I2_batch[genotype,:] = I2
    end
    
    # write batches to output file
    genostr = I1_batch + I2_batch
    properties = Dict("fid"=>[string("syn",x) for x in (batch_number*batchsize+1):(batch_number*batchsize+batchsize)],
            "iid"=>[string("syn",x) for x in (batch_number*batchsize+1):(batch_number*batchsize+batchsize)],
            #"sex"=>[g+1 for g in rand(Bernoulli(0.5), batchsize)], # TODO this needs to be consistent across all chromosomes
            "chromosome"=>[split(f, "\t")[1][4:end] for f in fixed_fields], 
            "sid"=>[split(f, "\t")[3] for f in fixed_fields],
            "bp_position"=>[split(f, "\t")[2] for f in fixed_fields], 
            "allele_1"=>[split(f, "\t")[5] for f in fixed_fields], 
            "allele_2"=>[split(f, "\t")[4] for f in fixed_fields])
    bed_reader.to_bed(batch_file, genostr, properties=properties)
end

function write_to_plink(H1, H2, ind_map, ref_df, fixed_fields, n_snps, n_syn_samples, outfile, batchsize, plink)
    number_of_batches = Int(ceil(n_syn_samples/batchsize))

    batch_files = Vector{String}(undef, number_of_batches)
    
    # write data to PLINK for each batch
    for batch in 1:number_of_batches
        # correct batchsize for final batch
        if batch == number_of_batches
            if batchsize >= n_syn_samples
                batchsize = n_syn_samples
            else
                batchsize = n_syn_samples - batchsize*(number_of_batches-1)
            end
        end

        @info @sprintf("Writing PLINK for batch %i", batch)
        batch_file = @sprintf("%s_%i.bed", outfile, batch-1)
        batch_files[batch] = batch_file[1:end-4]
        start_haplotype = ((batch-1)*batchsize)*2+1
        end_haplotype = start_haplotype+batchsize*2-1
        batch_ref_df = ref_df[(ref_df.H .<= end_haplotype).&(ref_df.H .>= start_haplotype), :]
        write_to_plink_batch(batch, batch_file, batchsize, H1, H2, ind_map, batch_ref_df, fixed_fields, n_snps)
    end

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