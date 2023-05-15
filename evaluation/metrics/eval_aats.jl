"""
    Analyze quality of Synthetic data using Nearest Neighbor Adversarial Accuracy.

# Usage:
"""

using CSV
using DataFrames
using Printf
using ProgressMeter


"""
    Helper function to decide whether on cluster
"""
is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

function run_aats(ibsfile_cross)
    @info "Running AA_TS computations"
    
    results_dir = dirname(ibsfile_cross)

    cols = [:ID1, :ID2, :Dist]
    ibs_cross = CSV.File(ibsfile_cross, delim='\t', select=cols) |> DataFrame

    ## Split into real, synthetic and cross parts
    df_real = filter([:ID1, :ID2] => 
        (id1, id2) -> (!startswith(id1, "syn") && !startswith(id2, "syn")),
        ibs_cross, view=true
    )
    df_synt = filter([:ID1, :ID2] =>
        (id1, id2) -> (startswith(id1, "syn") && startswith(id2, "syn")),
        ibs_cross, view=true
    )
    df_cross = filter(
        [:ID1, :ID2] => (id1, id2) -> 
            (startswith(id2, "syn") && !startswith(id1, "syn")) ||
            (startswith(id1, "syn") && !startswith(id2, "syn")),
        ibs_cross
    )
    
    n1 = size(df_real, 1)
    n2 = size(df_synt, 1)
    n3 = size(df_cross, 1)
    n0 = size(ibs_cross, 1)
    @assert n1 + n2 + n3 == n0

    if startswith(df_cross[1, :ID1], "syn") # First column hold Synthetic samples
        rename!(df_cross, [:ID1 => :SYNT, :ID2 => :REAL])
    else
        rename!(df_cross, [:ID1 => :REAL, :ID2 => :SYNT])
    end

    ids_real = unique(df_cross.REAL)
    ids_synt = unique(df_cross.SYNT)

    # TODO: To speed up the search, we can also convert IDs from string to number
    #       speeding up the comparison
    # Find d_TT and d_TS
    d_TT = Array{Float64}(undef, size(ids_real, 1))
    d_TS = Array{Float64}(undef, size(ids_real, 1))
    # @showprogress "Computing on real dataset" for (idx, id) ∈ enumerate(ids_real)
    p = Progress(length(ids_real); desc="Computing on reference dataset", enabled=!is_logging(stderr))
    for (idx, id) ∈ enumerate(ids_real)
        d_TT[idx] = minimum(
            df_real[(df_real.ID1 .== id) .| (df_real.ID2 .== id), :Dist];
            init=10.0 # Some high value in case of no ID in file
        )
        d_TS[idx] = minimum(
            df_cross[df_cross.REAL .== id, :Dist];
            init=10.0
        )
        ProgressMeter.next!(p)
    end
    # Find d_ST and d_SS
    d_SS = Array{Float64}(undef, size(ids_synt, 1))
    d_ST = Array{Float64}(undef, size(ids_synt, 1))
    # @showprogress "Computing on synthetic dataset" for (idx, id) ∈ enumerate(ids_synt)
    p = Progress(length(ids_real); desc="Computing on synthetic dataset", enabled=!is_logging(stderr))
    for (idx, id) ∈ enumerate(ids_synt)
        d_SS[idx] = minimum(
            df_synt[(df_synt.ID1 .== id) .| (df_synt.ID2 .== id), :Dist];
            init=10.0
        )
        d_ST[idx] = minimum(
            df_cross[df_cross.SYNT .== id, :Dist];
            init=10.0
        )
        ProgressMeter.next!(p)
    end

    n_real = size(ids_real, 1)
    n_synt = size(ids_synt, 1)

    AA_T = 1 / n_real * sum(d_TS .> d_TT)
    AA_S = 1 / n_synt * sum(d_ST .> d_SS)
    AA_TS = 0.50 * (AA_T + AA_S)

    @info @sprintf("AA_TS results")
    @info @sprintf("    AA_T  = %.4f", AA_T)
    @info @sprintf("    AA_S  = %.4f", AA_S)
    @info @sprintf("    AA_TS = %.4f", AA_TS)
end
