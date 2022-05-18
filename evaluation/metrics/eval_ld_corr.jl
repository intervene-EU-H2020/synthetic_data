using CSV
using DelimitedFiles
using Printf
using Plots
default(show=false) # to plot headless
ENV["GKSwstype"] = "100"

function run_ld_corr(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
    ids_path = @sprintf("%s.ld.extract", eval_dir)

    # @info "Loading IDs"
    # ids = CSV.File(@sprintf("%s.bim",real_data_prefix), delim='\t', normalizenames=true, select=[:2]) # |> DataFrame
    # nids = size(ids, 1) # Get number of IDs
    # @printf("Loaded %d IDs", nids)

    # Select sample with only > 0.001 MAF
    @info "Computing MAF for reference set"
    maf_file_prefix = "$real_data_prefix.ld_maf"
    run(`$plink_path --bfile $real_data_prefix --freq --out $maf_file_prefix`)
    # Load computed MAF
    real_maf_file = "$maf_file_prefix.frq"
    maf_real = CSV.File(real_maf_file, delim=' ', ignorerepeated=true) |> DataFrame
    sel_snps = maf_real.MAF .> 0.001
    ids = maf_real.SNP[sel_snps]
    nids = size(ids, 1)

    nSelect = 512 # Number of IDs to select randomly

    # idx = rand(1:nids, nSelect) # Select random indexes

    startidx = rand(1:nids - nSelect)
    idx = startidx:startidx+nSelect # Select a consecutive segment

    selected_ids = ids[idx]
    @info @sprintf("Selected %d IDs from %d", nSelect, size(ids, 1))
    open(ids_path, "w") do io
        writedlm(io, selected_ids)
    end
    # CSV.write(ids_path, selected_ids, writeheader=false)

    # Run Plink with --exatract flag to compute LD of selected variants
    ld_real_file = "$eval_dir.ld_real.extracted"
    ld_command = `$plink_path --bfile $real_data_prefix --extract $ids_path --r2 square --out $ld_real_file`
    run(ld_command)
    
    ld_synt_file = "$eval_dir.ld_synt.extracted"
    ld_command = `$plink_path --bfile $synt_data_prefix --extract $ids_path --r2 square --out $ld_synt_file`
    run(ld_command)

    ld_real = readdlm(ld_real_file * ".ld", '\t', Float64, '\n')
    ld_synt = readdlm(ld_synt_file * ".ld", '\t', Float64, '\n')

    # Plot the results
    fig = Plots.plot(layout=(1, 2), size=(2*900, 2*450))
    heatmap!(fig, subplot=1, ld_real, 
        clims=(0.0, 1.0), border=:none, yflip=true,
        title="Real data", aspect_ratio=:equal)
    heatmap!(fig, subplot=2, ld_synt, 
        clims=(0.0, 1.0), border=:none, yflip=true,
        title="Synthetic data", aspect_ratio=:equal)

    outfile = @sprintf("%s.ld_vis.png", eval_dir)
    savefig(fig, outfile)
    @info "LD Visualization saved to $outfile"
end
