using CSV
using DelimitedFiles
using Printf
using Plots
default(show=false) # to plot headless
ENV["GKSwstype"] = "100"

function main()
    ## Configure file paths
    plink_path = "software/plink" # relative to the root of the project

    real_data_vcffile = ARGS[1]
    real_data_prefix = ARGS[2]
    synt_data_prefix = ARGS[3]
    eval_dir = ARGS[4]

    ids_path = joinpath(eval_dir, "ld.extract")

    ## Extract 512 random variants
    ## Load IDs from real data
    @info "Loading IDs"
    ids = CSV.File(real_data_vcffile, delim='\t', normalizenames=true, 
        comment="##", select=[:ID]) # |> DataFrame
    nids = size(ids, 1) # Get number of IDs
    @printf("Loaded %d IDs", nids)

    nSelect = 512 # Number of IDs to select randomly

    # idx = rand(1:nids, nSelect) # Select random indexes

    startidx = rand(1:nids - nSelect)
    idx = startidx:startidx+nSelect # Select a consecutive segment

    selected_ids = ids[idx]
    @info @sprintf("Selected %d IDs", nSelect)
    CSV.write(ids_path, selected_ids, writeheader=false)

    # Run Plink with --exatract flag to compute LD of selected variants
    ld_real_file = "$eval_dir/ld_real.extracted"
    ld_command = `$plink_path --bfile $real_data_prefix --extract $ids_path --r2 square --out $ld_real_file`
    run(ld_command)

    ld_synt_file = "$eval_dir/ld_synt.extracted"
    ld_command = `$plink_path --bfile $synt_data_prefix --extract $ids_path --r2 square --out $ld_synt_file`
    run(ld_command)

    ld_real = readdlm(ld_real_file * ".ld", '\t', Float64, '\n')
    ld_synt = readdlm(ld_synt_file * ".ld", '\t', Float64, '\n')

    # Plot the results
    fig = plot(layout=(1, 2), size=(2*900, 2*450))
    heatmap!(fig, subplot=1, ld_real, 
        clims=(0.0, 1.0), border=:none, yflip=true,
        title="Real data", aspect_ratio=:equal)
    heatmap!(fig, subplot=2, ld_synt, 
        clims=(0.0, 1.0), border=:none, yflip=true,
        title="Synthetic data", aspect_ratio=:equal)

    outfile = joinpath(eval_dir, "ld_vis.png")
    savefig(fig, outfile)
    @info "LD Visualization saved to $outfile"
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
