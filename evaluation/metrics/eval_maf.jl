"""
    Minimum Allele Frequency results and visualization

    **Divergence**: For each of the variant, a Bernoulli distribution is assumed.
    Divergence between individual variants are calculated and mean over all variants
    is reported as evaluation metric.


Requirements:
1. `*.frq` files using Plink
"""

using CSV
using DataFrames
using Statistics
using Plots: default
using StatsPlots
default(show=false) # to plot headless
ENV["GKSwstype"] = "100"


## KL divergence
compute_entropy(p::T, q::T) where T <: AbstractFloat = - p * log2(q) - (1 - p) * log2(1 - q)
compute_entropy(p::T) where T <: AbstractFloat = compute_entropy(p, p)
compute_kldiv(p::T, q::T) where T <: AbstractFloat = compute_entropy(p, q) - compute_entropy(p)


## Wasserstein distance
compute_wsd(p::T, q::T) where T <: AbstractFloat = 2 * abs(p - q)


## Computations
function compute_divergences(maf_data)

    kldiv_vec = compute_kldiv.(maf_data.MAF_REAL, maf_data.MAF_SYNT)
    mean_kl = mean(filter(!isnan, kldiv_vec))


    wsd_vec = compute_wsd.(maf_data.MAF_REAL, maf_data.MAF_SYNT)
    mean_wsd = mean(wsd_vec)

    return Dict(:kl => mean_kl, :wd => mean_wsd)
end


## main
function run_maf(real_maf_file, synt_maf_file)
    @info "Running MAF evaluations"

    ## Load MAF data
    @info "Loading MAF files"
    maf_real = CSV.File(real_maf_file, delim=' ', ignorerepeated=true) |> DataFrame
    maf_synt = CSV.File(synt_maf_file, delim=' ', ignorerepeated=true) |> DataFrame

    # Clean up MAFs with 0.0 (or very low)
    sel_snps = maf_real.MAF .> 0.001
    maf_real = maf_real[sel_snps, :]
    maf_synt = maf_synt[sel_snps, :]
    
    # Consolidate variants available in both the tables
    maf_data = innerjoin(maf_real[!, [:SNP, :MAF]], maf_synt[!, [:SNP, :MAF]], 
                    on = :SNP,
                    renamecols = "_REAL" => "_SYNT")
    n_variants = size(maf_data, 1)
    @info "MAF information loaded for $n_variants variants."
                    
    ## Divergence computations
    results = compute_divergences(maf_data)

    @info "-----   MAF Eval results   -----"
    @info "KL Divergence: ", results[:kl]
    @info "Wasserstein Divergence: ", results[:wd]
    @info "--------------------------------"

    ## Scatter plot
    outfile = joinpath(dirname(synt_maf_file), "results-maf-scatter.png")
    fig = Plots.plot(size=(400, 400))
    @df maf_data[(maf_data.MAF_REAL .> 0) .& (maf_data.MAF_SYNT .> 0), :]  scatter!(fig,
            :MAF_REAL, :MAF_SYNT,
            label=nothing, 
            xaxis=("MAF (Real)"), 
            yaxis=("MAF Synthetic"),
            markersize=6,
            markerstrokewidth=1,
            title="MAF Values"
    )
    savefig(fig, outfile)
    @info "MAF comparison plot save at $outfile"
end
