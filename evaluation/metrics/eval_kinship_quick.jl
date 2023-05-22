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
        second_degree = kin_cross_filtered[0.0884 .< kin_cross_filtered.Kinship .<= 0.177, :]
    end
        
    # compute proportion of samples in synthetic data that are close relatives of reference data
    duplicate_samples = isfile(kinship_file) ? length(Set(duplicate.ID1))/nsamples_synfile : 0
    first_degree_samples = isfile(kinship_file) ? length(Set(first_degree.ID1))/nsamples_synfile : 0
    second_degree_samples = isfile(kinship_file) ? length(Set(second_degree.ID1))/nsamples_synfile : 0
    total_related_samples = duplicate_samples + first_degree_samples + second_degree_samples
    
    data = [duplicate_samples, first_degree_samples, second_degree_samples, total_related_samples]
    
    return data
end