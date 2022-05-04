using CategoricalArrays, Distances


"""
Function that implements the LD decay metric for ABC
"""
function LD_decay(plink_file, plink, mapthin, out_prefix, bp_to_cm_map)
    # thin snp list to speed up calculations
    # keep 20 snps every 10^6 bases and create a new bim
    bim_file = @sprintf("%s.bim", plink_file)
    snp_thin = @sprintf("%s_snp_thin", out_prefix)
    run(`$mapthin -n -b 20 $bim_file $snp_thin`)
    
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
    ld_df.bin = cut(ld_df.dist, Vector(0:interval:maximum(ld_df.dist)), extend=missing)
    gdf = groupby(ld_df, [:bin])
    ld_grp = combine(gdf, :R2 => mean)
    ld_grp.dist = Vector(0:interval:maximum(ld_df.dist))
    
    return [ld_grp.dist[2:end] ld_grp.R2_mean[2:end]]
end


function run_ld_decay(synfile, realfile, plink, mapthin, out_prefix, bp_to_cm_map)
    ld_decay_real = LD_decay(realfile, plink, mapthin, out_prefix, bp_to_cm_map)
    ld_decay_syn = LD_decay(synfile, plink, mapthin, out_prefix, bp_to_cm_map)
    # store data points
    df = DataFrame(real_x=ld_decay_real[:,1], real_y=ld_decay_real[:,2], syn_x=ld_decay_syn[:,1], syn_y=ld_decay_syn[:,2])
    outfile = joinpath(dirname(out_prefix), "results-ld-decay.csv")
    CSV.write(outfile, df)
    @info "LD decay data saved at $outfile"
    # make plot
    outfile = joinpath(dirname(out_prefix), "results-ld-decay.png")
    fig = Plots.plot(size=(400, 400))
    plot!(fig, ld_decay_real[:,1], ld_decay_real[:,2], label="real", xaxis="Genetic distance (cm)", yaxis="LD estimate (r2)", title="LD decay")
    plot!(fig, ld_decay_syn[:,1], ld_decay_syn[:,2], label="synthetic")
    savefig(fig, outfile)
    @info "LD decay plot saved at $outfile"
    ld_distance = evaluate(Euclidean(), ld_decay_real[:,2], ld_decay_syn[:,2])
    @info @sprintf("LD decay distance between real and synthetic data is %f", ld_distance)
end
