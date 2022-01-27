using CSV, DataFrames

"""Utility function for creating reference datasets
needed for evaluation and optimisation pipelines,
i.e. using a subset of ancestries from the original data
"""
function create_reference_dataset(vcfpath, poppanel, population_weights, plink, outdir)
    @info "Creating reference dataset"
    outfile = @sprintf("%s/reference", outdir)
    keeplist, nsamples = create_keepfile(poppanel, population_weights, outdir)
    # convert vcf to plink, keeping only the ancestries used in data generation
    run(`$plink --vcf $vcfpath --keep $keeplist --make-bed --out $outfile`)
    return outfile, nsamples
end


"""Helper function that returns list of populations used for data generation
"""
function get_population_list(population_weights)
    poplist = []
    for k in keys(population_weights)
        for pop in keys(population_weights[k])
            if pop ∉ poplist
                push!(poplist, pop)
            end
        end
    end
    return poplist
end


"""Create a file specifying which samples to keep in the reference data
"""
function create_keepfile(poppanel, population_weights, outdir)
    poplist = get_population_list(population_weights)
    poppanel_df = CSV.File(poppanel, normalizenames=true, select=["FamilyID", "SampleID", "Superpopulation"]) |> DataFrame
    keep_samples = filter(row -> row.Superpopulation ∈ poplist, poppanel_df)
    keep_samples = keep_samples[!, [:FamilyID, :SampleID]]
    keepfile = @sprintf("%s/keep.txt", outdir)
    CSV.write(keepfile, keep_samples, header=false, delim="\t")
    nsamples = nrow(keep_samples)
    return keepfile, nsamples
end
