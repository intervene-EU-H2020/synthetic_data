# TO plot kinship vs IBS Plots

using CSV
using DataFrames
using Printf
using Plots: default
using Plots.PlotMeasures
using StatsPlots
default(show=false) # to plot headless
ENV["GKSwstype"] = "100"


function plot_kinship_ibs(ibs_real, ibs_synt, results_dir)
    @info "Plotting Kinship vs IBS"
    
    # Randomly sample 15K pairs from each for plotting
    idx_real = rand(1:size(ibs_real, 1), 100000)
    idx_synt = rand(1:size(ibs_synt, 1), 100000)


    df_real = @view ibs_real[idx_real, :]
    df_synt = @view ibs_synt[idx_synt, :]

    plt = Plots.plot(layout=(1, 3), size=(3*400, 400), 
            left_margin = [5mm 0mm], bottom_margin = 10mm)

    @df df_real scatter!( plt, subplot=1,
        :F_IBS0, :Kinship, label="Real",
        xaxis=("Fraction of IBS=0"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:red, alpha=0.30
    )
    @df df_synt scatter!( plt, subplot=1,
        :F_IBS0, :Kinship, label="Synthetic",
        xaxis=("Fraction of IBS=0"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:blue, alpha=0.30
    )

    @df df_real scatter!( plt, subplot=2,
        :F_IBS1, :Kinship, label="Real",
        xaxis=("Fraction of IBS=1"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:red, alpha=0.30
    )
    @df df_synt scatter!( plt, subplot=2,
        :F_IBS1, :Kinship, label="Synthetic",
        xaxis=("Fraction of IBS=1"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:blue, alpha=0.30
    )

    @df df_real scatter!( plt, subplot=3,
        :F_IBS2, :Kinship, label="Real",
        xaxis=("Fraction of IBS=2"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:red, alpha=0.30
    )
    @df df_synt scatter!( plt, subplot=3,
        :F_IBS2, :Kinship, label="Synthetic",
        xaxis=("Fraction of IBS=2"), yaxis=("Kinship"),
        markersize=1, markerstrokewidth=1, color=:blue, alpha=0.30
    )

    outfile = joinpath(results_dir, "results-kinship-ibs.png")
    savefig(plt, outfile)
    @info "Output saved to $outfile"
end


function plot_kin_density(ibs_real, ibs_synt, ibs_cross_filtered, results_dir)
    @info "Plotting relatedness metrics density"
    plt = Plots.plot(layout=(3, 1), size=(800, 3*450))

    density!(plt, subplot=1, ibs_real[!, :Kinship], label="Real data")
    density!(plt, subplot=1, ibs_synt[!, :Kinship], label="Synthetic data")
    density!(plt, subplot=1, ibs_cross_filtered[!, :Kinship], label="Cross")
    plot!(plt, subplot=1, xlabel="Kinship coefficient", ylabel="Density", title="Kinship")

    density!(plt, subplot=2, ibs_real[!, :IBS], label="Real data")
    density!(plt, subplot=2, ibs_synt[!, :IBS], label="Synthetic data")
    density!(plt, subplot=2, ibs_cross_filtered[!, :IBS], label="Cross")
    plot!(plt, subplot=2, xlabel="IBS", ylabel="Density", title="IBS")

    density!(plt, subplot=3, ibs_real[!, :Dist], label="Real data")
    density!(plt, subplot=3, ibs_synt[!, :Dist], label="Synthetic data")
    density!(plt, subplot=3, ibs_cross_filtered[!, :Dist], label="Cross")
    plot!(plt, subplot=3, xlabel="Distance", ylabel="Density", title="Distance")

    outfile = joinpath(results_dir, "results-kinship-density.png")
    savefig(plt, outfile)
    @info "Output saved to $outfile"
end

function summarize_kinship(kinship_values)
    # Summarize kinsip distribution
    #   "an estimated kinship coefficient range >0.354, [0.177, 0.354], 
    #   [0.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ 
    #   twin, 1st-degree, 2nd-degree, and 3rd-degree relationships 
    #   respectively."
    # Ref: https://www.kingrelatedness.com/manual.shtml

    n = size(kinship_values, 1)

    n1 = sum(kinship_values .> 0.354)
    n2 = sum(0.177 .< kinship_values .<= 0.354)
    n3 = sum(0.0884 .< kinship_values .<= 0.177)
    n4 = sum(0.0442 .< kinship_values .<= 0.0884)
    n5 = sum(kinship_values .<= 0.0442)

    @assert n == n1 + n2 + n3 + n4 + n5

    @info "Kinship distribution"
    @info @sprintf("    MZ             : %9d    (%.8f)", n1, n1/n)
    @info @sprintf("    1-st Degree    : %9d    (%.8f)", n2, n2/n)
    @info @sprintf("    2-nd Degree    : %9d    (%.8f)", n3, n3/n)
    @info @sprintf("    3-rd Degree    : %9d    (%.8f)", n4, n4/n)
    @info @sprintf("    Unrelated      : %9d    (%.8f)", n5, n5/n)
    @info @sprintf("                     %9d", n)
end


function run_kinship(ibsfile_real, ibsfile_synt, ibsfile_cross)
    @info "Running Kinship evaluations"

    results_dir = dirname(ibsfile_synt)

    cols = [:ID1, :ID2, :N_SNP, :N_IBS0, :N_IBS1, :N_IBS2, :IBS, :Kinship, :Dist]

    ibs_synt = CSV.File(ibsfile_synt, delim='\t', select=cols) |> DataFrame
    ibs_real = CSV.File(ibsfile_real, delim='\t', select=cols) |> DataFrame
    ibs_cross = CSV.File(ibsfile_cross, delim='\t', select=cols) |> DataFrame 

    ibs_cross_filtered = filter(
        [:ID1, :ID2] => (id1, id2) -> 
            (startswith(id2, "syn") && !startswith(id1, "syn")) ||
            (startswith(id1, "syn") && !startswith(id2, "syn")),
        ibs_cross, view=true)

    # Compute fraction of IBS
    ibs_synt.F_IBS0 = ibs_synt.N_IBS0 ./ ibs_synt.N_SNP
    ibs_synt.F_IBS1 = ibs_synt.N_IBS1 ./ ibs_synt.N_SNP
    ibs_synt.F_IBS2 = ibs_synt.N_IBS2 ./ ibs_synt.N_SNP

    ibs_real.F_IBS0 = ibs_real.N_IBS0 ./ ibs_real.N_SNP
    ibs_real.F_IBS1 = ibs_real.N_IBS1 ./ ibs_real.N_SNP
    ibs_real.F_IBS2 = ibs_real.N_IBS2 ./ ibs_real.N_SNP

    plot_kinship_ibs(ibs_real, ibs_synt, results_dir)

    plot_kin_density(ibs_real, ibs_synt, ibs_cross_filtered, results_dir)

    @info "Kinship analysis for Real data"
    summarize_kinship(ibs_real.Kinship)
    @info "Kinship analysis for Synthetic data"
    summarize_kinship(ibs_synt.Kinship)
    @info "Kinship across both datasets"
    summarize_kinship(ibs_cross_filtered.Kinship)

end
