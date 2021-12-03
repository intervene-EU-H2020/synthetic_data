"""
    PCA based anlaysis of synthetic data

"""

using CSV
using DataFrames
using LinearAlgebra
using DelimitedFiles




## File paths
## Configure file paths
data_dir = "./../../../data/1000G";



function main()
    @info "Running PCA Evaluations"

    real_data_prefix = ARGS[1]
    synt_data_prefix = ARGS[2]

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

    ## Compute how aligned each of the PCs are -- ideally this should be close to 1
    align_val = diag(pc_real_norm' * pc_synt_norm)

    # Compute weighted alignment according to variance captured -- ideally close to 1
    wgt_align = sum(var_real .* align_val)

    @info "Weighted alignment value = $wgt_align    (Ideal is 1.0)"
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
