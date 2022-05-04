"""
    PCA based anlaysis of synthetic data

"""

using CSV
using DataFrames
using LinearAlgebra
using DelimitedFiles
using Plots
default(show=false) # to plot headless


function pca_alignment(real_data_prefix, synt_data_prefix)
    ## Load Allele counts
    real_acount = CSV.File(real_data_prefix * ".acount", normalizenames=true) |> DataFrame
    synt_acount = CSV.File(synt_data_prefix * ".acount", normalizenames=true) |> DataFrame

    ## Load projection vectors
    df_pc_real = CSV.File(real_data_prefix * ".eigenvec.allele", normalizenames=true) |> DataFrame
    df_pc_synt = CSV.File(synt_data_prefix * ".eigenvec.allele", normalizenames=true) |> DataFrame

    ## Load variance contribution of PCs
    pc_real_eigval = readdlm(real_data_prefix * ".eigenval"; header=false)
    pc_synt_eigval = readdlm(synt_data_prefix * ".eigenval"; header=false) 

    ##
    # Extract principal directions
    pc_real = Matrix(df_pc_real[:, 6:end])
    pc_synt = Matrix(df_pc_synt[:, 6:end])

    # Normalize the vectors to account for sample difference
    pc_real_norm = pc_real ./ .√(sum(pc_real .^ 2, dims=1))
    pc_synt_norm = pc_synt ./ .√(sum(pc_synt .^ 2, dims=1))

    # Normalize variance contributions
    var_real = pc_real_eigval ./ sum(pc_real_eigval)
    var_synt = pc_synt_eigval ./ sum(pc_synt_eigval)

    ## Iteratively find out which PC in synthetic data is closest match to real data
    n_pcs = size(pc_real_norm, 2) # number of PCs
    idx = collect(1:n_pcs) # Indices of synthetic PCs to check
    align_val = []
    for i ∈ 1:n_pcs
        # For i-th real PC, find the closest match
        _pcs = pc_synt_norm[:, idx]
        _align_vals = pc_real_norm[:, i]' * _pcs
        _align_vals = vec(abs.(_align_vals)) # Direction may be different
        v, j = findmax(_align_vals) # find the best pc align value and it's index
        deleteat!(idx, j) # remove this pc index from list
        append!(align_val, v) # add this pc alignment to the align value
    end

    # Compute weighted alignment according to variance captured -- ideally close to 1
    wgt_align = sum(var_real .* align_val)

    @info "Weighted alignment value = $wgt_align    (Ideal is 1.0)"
end


"""Function to plot PCA of real/synthetic datasets and compare superpopulations
"""
function pca_plot(real_data_prefix, synt_data_prefix, real_data_popfile, synt_data_popfile, eval_dir)
    # Load first 2 PCs of real and synthetic datasets
    real_eigvec = DataFrame(CSV.File(string(real_data_prefix, ".eigenvec")))
    synt_eigvec = DataFrame(CSV.File(string(synt_data_prefix, ".eigenvec")))

    # Load population files 
    real_pop = readdlm(real_data_popfile, header=false)
    synt_pop = readdlm(synt_data_popfile, header=false)
    num_real_pop = length(unique(real_pop))
    num_synt_pop = length(unique(synt_pop))

    # Assign colors for plotting
    # TODO just hardcoded this to quickly make an example, but should be generalised to arbitrary populations
    superpopulations = ["AFR", "AMR", "EAS", "EUR", "CSA", "MID"]
    col_map = Dict(zip(superpopulations, distinguishable_colors(length(superpopulations))))
    real_col_map = [haskey(col_map, x) ? col_map[x] : "blue" for x in real_pop]
    synt_col_map = [haskey(col_map, x) ? col_map[x] : "blue" for x in synt_pop]

    # Plot first 2 PCs of real and synthetic datasets, color by superpopulation
    fig = Plots.plot(layout=(1, 2), size=(2*900, 2*450))
    scatter!(fig, subplot=1, real_eigvec.PC1, real_eigvec.PC2, title="Real data", color=real_col_map, legend=nothing)
    scatter!(fig, subplot=2, synt_eigvec.PC1, synt_eigvec.PC2, title="Synthetic data", color=synt_col_map, legend=nothing)

    outfile = @sprintf("%s.pca_vis.png", eval_dir)
    savefig(fig, outfile)
    @info "PCA Visualization saved to $outfile"
end

"""Function to plot projection of synthetic samples on PCs of real data
"""
function plot_pca_proj(pcaproj_file, eval_dir)
    df = CSV.File(pcaproj_file) |> DataFrame
    # Split real and synthetics
    df_real = df[df.AFF .== 1, [:PC1, :PC2]]
    df_synt = df[df.AFF .== 2, [:PC1, :PC2]]

    # Plot
    fig = Plots.plot(size=(640, 640))
    scatter!(fig, df_real.PC1, df_real.PC2, color=:blue, label="Reference")
    scatter!(fig, df_synt.PC1, df_synt.PC2, color=:red, label="Synthetic")

    outfile = @sprintf("%s.pcaproj.png", eval_dir)
    savefig(fig, outfile)
    @info "PCA Projection Visualization saved to $outfile"
end


function run_pca(real_data_prefix_pca, synt_data_prefix_pca, pcaproj_file, real_data_popfile, synt_data_popfile, eval_dir)
    @info "Running PCA Evaluations"
    pca_alignment(real_data_prefix_pca, synt_data_prefix_pca)
    pca_plot(real_data_prefix_pca, synt_data_prefix_pca, real_data_popfile, synt_data_popfile, eval_dir)
    plot_pca_proj(pcaproj_file, eval_dir)
end
