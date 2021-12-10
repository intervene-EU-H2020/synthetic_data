using DelimitedFiles, LsqFit, DataFrames, CSV, Mmap


"""Using the vcftools software, create a VCF file retaining only the variants in a given list
"""
function extract_variants(vcftools, vcf_input, vcf_output_prefix, vcf_output, variant_list)
    vcf_cmd = `$vcftools --gzvcf $vcf_input --snps $variant_list --out $vcf_output_prefix --recode`
    run(vcf_cmd)
    
    @assert isfile(vcf_output) "Error occurred with creation of VCF file"
end


"""Create a map of basepair distances to centimorgan distances
"""
# TODO method for converting basepair to centimorgan distances is not accurate - needs to be replaced
function get_genetic_distances(datafile, mapfile, output)
    df = CSV.read(mapfile, DataFrame; header=["ID", "POS", "MAP"])

    xdata = [convert(Float64,x) for x in df.POS]
    ydata = [convert(Float64,y) for y in df.MAP]
    
    poly3(x, p) = p[1] .* xdata .^ 3 .+ p[2] .* xdata .^ 2 .+ p[3] .* xdata .+ p[4]
    p0 = [1.0, 0.0, 0.0, 0.0]
    fit = curve_fit(poly3, xdata, ydata, p0)
    
    # predict genetic distances in cM for basepair distances given in datefile
    bpdata = Float64[]
    snpid = String[]
    for line in eachline(datafile)
        if !startswith(line, "#")
            push!(bpdata, parse(Float64,split(line)[2]))
            push!(snpid, split(line)[3])
        end
    end
    
    cmpred = fit.param[1] .* bpdata .^ 3 .+ fit.param[2] .* bpdata .^ 2 .+ fit.param[3] .* bpdata .+ fit.param[4]
    if cmpred[1] < 0
        # shift negative values
        cmpred = cmpred .+ abs(cmpred[1])
    end
    
    df = DataFrame(Index = 1:length(snpid), Variant = snpid, Distance = cmpred)
    CSV.write(output, df)
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
