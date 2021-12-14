# TODO refactor and move these to evaluation directory


"""
Function that implements the LD decay metric for ABC

Arguments
- plink_file::String: string giving filepath to the plink file prefix containing genotype
"""
function LD_decay(plink_file::String)
    # thin snp list to speed up calculations
    # keep 20 snps every 10^6 bases and create a new bim
    bim_file = @sprintf("%s.bim", plink_file)
    snp_thin = @sprintf("%s_snp-thin_%s", parsed_args["outfile_prefix"], parsed_args["opt_superpop"])
    map_thin = parsed_args["mapthin"]
    plink = parsed_args["plink"]
    run(`$map_thin -n -b 20 $bim_file $snp_thin`)
    
    # compute r2 values between all pairs
    run(`$plink --bfile $plink_file --extract $snp_thin --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out $snp_thin`)
    
    # extract the distance between each pair and the r2 value
    ld_df = DataFrame(CSV.File(@sprintf("%s.ld", snp_thin), delim=' ', ignorerepeated=true))
    # convert bp to cm
    ld_df.CM_A = [bp_to_cm_map[x] for x in ld_df.BP_A]
    ld_df.CM_B = [bp_to_cm_map[x] for x in ld_df.BP_B]
    ld_df.dist = ld_df.CM_B - ld_df.CM_A
    
    # categorise distances into intervals of a fixed length and compute mean r2 within each interval
    interval = 0.1 
    ld_df.bin = cut(ld_df.dist, Vector(0:interval:maximum(ld_df.dist)), extend=true)
    gdf = groupby(ld_df, [:bin])
    ld_grp = combine(gdf, :R2 => mean)
    ld_grp.dist = Vector(0:interval:maximum(ld_df.dist))
    
    return [ld_grp.dist ld_grp.R2_mean]
end


function kinship_cross(plink_file, plink_file_ref)
    king = parsed_args["king"]

    output = @sprintf("%s_%s.cross.king", parsed_args["outfile_prefix"], parsed_args["opt_superpop"])
    input = @sprintf("%s.bed", plink_file)
    cross_input = @sprintf("%s.bed", plink_file_ref)
    
    run(`$king -b $input,$cross_input --related --prefix $output`)

    nsamples = 8000
    total_number_pairs = 3202 * nsamples # TODO don't hardcode these numbers

    cols = [:ID1, :ID2, :N_SNP, :HetHet, :IBS0, :HetConc, :HomIBS0, :Kinship, :IBD1Seg, :IBD2Seg, :PropIBD, :InfType]
    
    kin_cross = CSV.File(@sprintf("%s.kin0", output), delim='\t', select=cols) |> DataFrame 

    kin_cross_filtered = filter(
        [:ID1, :ID2] => (id1, id2) -> 
            (startswith(id2, "syn") && !startswith(id1, "syn")) ||
            (startswith(id1, "syn") && !startswith(id2, "syn")),
        kin_cross, view=true)
    
    duplicate = kin_cross_filtered[in(["Dup/MZ"]).(kin_cross_filtered.InfType), :]
    first_degree = kin_cross_filtered[in(["FS","PO"]).(kin_cross_filtered.InfType), :]
    second_degree = kin_cross_filtered[in(["2nd"]).(kin_cross_filtered.InfType), :]

    duplicate_pairs = size(duplicate)[1]/total_number_pairs
    first_degree_pairs = size(first_degree)[1]/total_number_pairs
    second_degree_pairs = size(second_degree)[1]/total_number_pairs
    total_related_pairs = duplicate_pairs + first_degree_pairs + second_degree_pairs
    
    duplicate_samples = length(Set(duplicate.ID1))/nsamples
    first_degree_samples = length(Set(first_degree.ID1))/nsamples
    second_degree_samples = length(Set(second_degree.ID1))/nsamples
    total_related_samples = duplicate_samples + first_degree_samples + second_degree_samples
    
    data = [duplicate_pairs, first_degree_pairs, second_degree_pairs, total_related_pairs, duplicate_samples, first_degree_samples, second_degree_samples, total_related_samples]
    print(data)
    
    return [1:length(data) data]
end