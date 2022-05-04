using DelimitedFiles, LsqFit, DataFrames, CSV, Mmap, Impute


"""Using the vcftools software, create a VCF file retaining only the variants in a given list
"""
function extract_variants(vcftools, vcf_input, vcf_output_prefix, vcf_output, variant_list, remove_list)
    vcf_cmd = `$vcftools --gzvcf $vcf_input --snps $variant_list --remove $remove_list --out $vcf_output_prefix --recode`
    run(vcf_cmd)
    
    @assert isfile(vcf_output) "Error occurred with creation of VCF file"
end


"""Get a list of SNP identifiers in the datafile (basepair position format)
"""
function get_snpids(datafile)
    snpid = String[]
    for line in eachline(datafile)
        if !startswith(line, "#")
            push!(snpid, split(line)[3])
        end
    end
    return snpid
end


"""Create a dataframe with both the position Id and rs id
"""
function get_id_df(datafile, rsidfile)
    rs_df = CSV.read(rsidfile, DataFrame)
    rename!(rs_df,[:Id,:rsId])
    snpid = get_snpids(datafile)
    snp_df = DataFrame(Id=snpid)
    snp_df = innerjoin(rs_df, snp_df, on=:Id)
    @assert nrow(snp_df)==get_number_variants(datafile) "Not all variants have a matching RSid" 
    return snp_df
end


"""Create a map of variant positions to age of mutation
"""
function get_mutation_ages(datafile, mapfile, rsidfile, output)
    # open the mapfile and extract the relevant columns
    df = CSV.read(mapfile, DataFrame, header=4)
    select!(df, "VariantID", " AgeMean_Mut", " DataSource")
    rename!(df, [:rsId,:AgeMean_Mut,:DataSource])

    # merge on variants in preprocessed file
    snp_df = get_id_df(datafile, rsidfile)
    df = innerjoin(df, snp_df, on=:rsId)

    # there may be duplicate rows from multiple sources - in this case take the "Combined" estimate, followed by "TGP" or "SGDP" 
    df = combine(first, groupby(sort(df, :DataSource), :rsId))

    # impute the rest with the mean value
    final_df = leftjoin(snp_df, select!(df, Not(:Id)), on=:rsId)
    mean = sum(skipmissing(final_df.AgeMean_Mut))/count(!ismissing, final_df.AgeMean_Mut)
    replace!(final_df.AgeMean_Mut, missing => mean)
    replace!(final_df.DataSource, missing => "Imputed")

    # convert from years to generations 
    years_per_gen = 25
    final_df.AgeMean_Mut = round.(final_df.AgeMean_Mut/years_per_gen, digits=3)
    
    # save to output
    rename!(final_df, :Id => :Variant)
    final_df = sort_by_variant(final_df)
    df = DataFrame(Index = 1:nrow(final_df), Variant = final_df.Variant, Age = final_df.AgeMean_Mut)
    CSV.write(output, df)
end


"""Sort dataframe by position in "Variant" field
"""
function sort_by_variant(df)
    get_pos(variant) = [parse(Int64, x[2]) for x in split.(variant, ":")]
    transform!(df, :Variant => get_pos => :Pos)
    sort!(df, :Pos)
    select!(df, Not(:Pos))
    return df
end


"""Create a map of basepair distances to centimorgan distances
"""
function get_genetic_distances(datafile, mapfile, rsidfile, output)
    rsid_df = get_id_df(datafile, rsidfile)
    genetic_df = CSV.read(mapfile, DataFrame; header=["rsId", "pos", "map"])
    dist_df = leftjoin(rsid_df, genetic_df, on=:rsId)
    select!(dist_df, "Id", "map")
    rename!(dist_df, [:Variant,:Distance])
    dist_df = sort_by_variant(dist_df)
    # use linear interpolation on missing distances
    dist_df = Impute.interp(dist_df) |> Impute.locf() |> Impute.nocb()
    final_df = DataFrame(Index = 1:nrow(dist_df), Variant = dist_df.Variant, Distance = dist_df.Distance)
    CSV.write(output, final_df)
end


"""Create haplotype matrices from VCF files
"""
function convert_vcf_to_hap(datafile, hap1_output, hap2_output)
    num_variants = get_number_variants(datafile)
    num_samples = get_number_samples(datafile)
    
    s1 = open(hap1_output, "w+")
    s2 = open(hap2_output, "w+")
    
    # the first two entires written to the output are the matrix dimensions
    write(s1, num_variants)
    write(s2, num_variants)
    write(s1, num_samples)
    write(s2, num_samples)
    
    # write the snp array data in batches
    H1 = []
    H2 = []
    for line in eachline(datafile)
        if !startswith(line, "#")
            haplotype_data = split(line)[10:end]
            push!(H1, [parse(Int8, h[1]) for h in haplotype_data])
            push!(H2, [parse(Int8, h[3]) for h in haplotype_data])
        end
    end
    
    # write the snp array data for the final batch
    write(s1, vcat(H1...))
    write(s2, vcat(H2...))

    close(s1)
    close(s2)
end


"""Store metadata from VCF files
"""
function convert_vcf_to_meta(datafile, output)
    snp_info = []
    for line in eachline(datafile)
        if !startswith(line, "#")
            push!(snp_info, string(join(split(line)[1:9], "\t"), "\t"))
        end
    end
    writedlm(output, snp_info)
end


"""
Returns the number of samples in the VCF input
"""
function get_number_samples(datafile)
    num_samples = 0
    for line in eachline(datafile)
        if !startswith(line, "#")
            haplotype_data = split(line)[10:end]
            num_samples = size([parse(Int64, h[1]) for h in haplotype_data])[1]
            break
        end
    end
    return num_samples
end


"""
Returns the number of variants in the VCF input
"""
function get_number_variants(datafile)
    num_variants = 0
    for line in eachline(datafile)
        if !startswith(line, "#")
            num_variants += 1
        end
    end
    return num_variants
end
