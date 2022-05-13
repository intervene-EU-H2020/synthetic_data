"""
Function that implements an efficient relatedness metric for ABC
"""
function kinship_cross(synfile, nsamples_synfile, reffile, king, out_prefix)
    output = @sprintf("%s.cross.king", out_prefix)
    input = @sprintf("%s.bed", synfile)
    cross_input = @sprintf("%s.bed", reffile)
    
    run(`$king -b $input,$cross_input --related --prefix $output`)

    kinship_file = @sprintf("%s.kin0", output)
    if isfile(kinship_file)
        kin_cross = CSV.File(kinship_file, delim='\t') |> DataFrame 
        
        kin_cross_filtered = filter(
            [:ID1, :ID2] => (id1, id2) -> 
                (startswith(id2, "syn") && !startswith(id1, "syn")) ||
                (startswith(id1, "syn") && !startswith(id2, "syn")),
            kin_cross, view=true)
        
        duplicate = kin_cross_filtered[kin_cross_filtered.Kinship .> 0.354, :]
        first_degree = kin_cross_filtered[0.177 .< kin_cross_filtered.Kinship .<= 0.354, :]
    end
        
    # compute proportion of samples in synthetic data that are close relatives of reference data
    duplicate_samples = isfile(kinship_file) ? length(Set(duplicate.ID1))/nsamples_synfile : 0
    first_degree_samples = isfile(kinship_file) ? length(Set(first_degree.ID1))/nsamples_synfile : 0
    total_related_samples = duplicate_samples + first_degree_samples
    
    data = [duplicate_samples, first_degree_samples, total_related_samples]
    
    return [1:length(data) data]
end


function kinship_syn(synfile, nsamples_synfile, king, out_prefix)
    output = @sprintf("%s.syn.king", out_prefix)
    input = @sprintf("%s.bed", synfile)

    run(`$king -b $input --related --prefix $output`)

    kinship_file = @sprintf("%s.kin0", output)
    if isfile(kinship_file)
        kin_data = CSV.File(kinship_file, delim='\t') |> DataFrame 
        duplicate = kin_data[kin_data.Kinship .> 0.354, :]
        first_degree = kin_data[0.177 .< kin_data.Kinship .<= 0.354, :]
    end
    
    # compute number of related pairs in synthetic samples
    duplicate_pairs = isfile(kinship_file) ? nrow(duplicate) : 0
    first_degree_pairs = isfile(kinship_file) ? nrow(first_degree) : 0
    total_related_pairs = duplicate_pairs + first_degree_pairs

    duplicate_samples = isfile(kinship_file) ? length(Set(hcat(duplicate.ID1, duplicate.ID2))) : 0
    first_degree_samples = isfile(kinship_file) ? length(Set(hcat(first_degree.ID1, first_degree.ID2))) : 0
    total_related_samples = duplicate_samples + first_degree_samples
    
    data = [duplicate_pairs, first_degree_pairs, total_related_pairs, duplicate_samples, first_degree_samples, total_related_samples]
    
    return [1:length(data) data]
end


function kinship_quick(synfile, nsamples_synfile, reffile, king, out_prefix)
    syn = kinship_syn(synfile, nsamples_synfile, king, out_prefix)
    df = DataFrame(metric=["syn_duplicate_pairs", "syn_1st_degree_pairs", "syn_combined_pairs", "syn_duplicate_samples", "syn_1st_degree_samples", "syn_combined_samples"], value=[syn[1,2], syn[2,2], syn[3,2], syn[4,2], syn[5,2], syn[6,2]])
    outfile = joinpath(dirname(out_prefix), "results-relatedness.csv")
    CSV.write(outfile, df)
    @info "Relatedness data saved at $outfile"
end