using CSV, DataFrames, DelimitedFiles

"""Utility function for creating reference datasets
needed for evaluation and optimisation pipelines,
i.e. using a subset of ancestries from the original data
"""
function create_reference_dataset(vcfpath, poppanel, population_weights, plink, outdir, chromosome)
    @info "Creating reference dataset"
    outfile = @sprintf("%s/reference%s", outdir, chromosome)
    keeplist, nsamples = create_keepfile(poppanel, population_weights, outdir, outfile)
    # convert vcf to plink, keeping only the ancestries used in data generation
    run(`$plink --vcf $vcfpath --keep $keeplist --make-bed --out $outfile`)
    return outfile, nsamples
end


"""Splits the reference dataset into 2 datasets (e.g. for calculating cross relatedness)
"""
function create_cross_datasets(ref_prefix, plink)
    # make the cross dataset
    fam_df = DataFrame(CSV.File(@sprintf("%s.fam", ref_prefix), header=["FID","IID","F","M","S","P"]))
    fam_df.R = rand(size(fam_df)[1])
    sort!(fam_df, [:R])
    fam_df = fam_df[:, [:FID, :IID]]
    fam_file1 = @sprintf("%s_1.txt", ref_prefix)
    fam_file2 = @sprintf("%s_2.txt", ref_prefix)
    threshold = Int(ceil(size(fam_df)[1]/2))
    fam1_df = fam_df[1:threshold,:]
    fam2_df = fam_df[threshold+1:size(fam_df)[1],:]
    CSV.write(fam_file1, fam1_df, delim="\t")
    CSV.write(fam_file2, fam2_df, delim="\t")
    
    cross1_prefix = @sprintf("%s_1", ref_prefix)
    cross2_prefix = @sprintf("%s_2", ref_prefix)
    nsamples_cross1 = nrow(fam1_df)
    nsamples_cross2 = nrow(fam2_df)

    run(`$plink --bfile $ref_prefix --keep $fam_file1 --make-bed --out $cross1_prefix`)
    run(`$plink --bfile $ref_prefix --keep $fam_file2 -make-bed --out $cross2_prefix`)
    
    return cross1_prefix, cross2_prefix, nsamples_cross1, nsamples_cross2
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
function create_keepfile(poppanel, population_weights, outdir, outfile)
    # filter relevant populations
    poplist = get_population_list(population_weights)
    poppanel_df = CSV.File(poppanel, normalizenames=true, select=["FamilyID", "SampleID", "Superpopulation"]) |> DataFrame
    keep_samples_df = filter(row -> row.Superpopulation ∈ poplist, poppanel_df)
    # create keep file
    keep_samples = keep_samples_df[!, [:FamilyID, :SampleID]]
    keep_samples.FamilyID = keep_samples.SampleID # it seems the --keep command doesn't work if these columns aren't identical
    keepfile = @sprintf("%s/keep.txt", outdir)
    CSV.write(keepfile, keep_samples, header=false, delim="\t")
    nsamples = nrow(keep_samples)
    # create sample file
    samplefile = @sprintf("%s.sample", outfile)
    open(samplefile, "w") do io
        writedlm(io, keep_samples_df.Superpopulation)
    end
    return keepfile, nsamples
end


function create_bp_cm_ref(genetic_distfile)
    # create a mapping to convert bp to cm distances
    bp_to_cm_df = CSV.read(genetic_distfile, DataFrame)
    bp_to_cm_df.BP = [parse(Int64,split(x,":")[2]) for x in bp_to_cm_df.Variant]
    bp_to_cm_map = Dict(zip(bp_to_cm_df.BP,bp_to_cm_df.Distance))
    return bp_to_cm_map
end