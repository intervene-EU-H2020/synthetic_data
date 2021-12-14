"""Run likelihood-free inference for optimising model parameters
"""


"""
Simulator function for the optimisation pipeline 


Arguments
- var_params::Vector{Float64}: list of the form [Ne_s, rho_s]
"""
# TODO refactor code
function simulator_function(var_params::AbstractArray{Float64,1})
    Ne_s = var_params[1]
    rho_s = var_params[2]
    
    outfile_prefix = @sprintf("%s_%s", parsed_args["outfile_prefix"], parsed_args["opt_superpop"])
    
    # simulate a synthetic dataset with these parameter values
    generate_synthetic_data(parsed_args["hapfile1"], parsed_args["hapfile2"], ind_map, variant_dist, poplist, popdist, outfile_prefix, "plink", parsed_args["batchsize"], rho_s, Ne_s, parsed_args["sigma"], parsed_args["vcffile"], parsed_args["plink"])
    
    # calculate the summary statistic
    if parsed_args["opt_metric"]=="ld_decay"
        sumstat = LD_decay(outfile_prefix)
    elseif parsed_args["opt_metric"]=="kinship"
        sumstat = kinship_cross2(outfile_prefix, parsed_args["opt_reference"])
    elseif parsed_args["opt_metric"]=="combined"
        # concatenate multiple summary statistics
        sumstat_ld = LD_decay(outfile_prefix)
        sumstat_kin = kinship(outfile_prefix)
        sumstat = vcat(sumstat_ld, sumstat_kin)
    end
    
    return sumstat
end


function get_priors(options)
    # uses uniform priors
    Ne_lower = options["optimisation"]["priors"]["Ne"]["uniform_lower"]
    Ne_upper = options["optimisation"]["priors"]["Ne"]["uniform_upper"]
    rho_lower = options["optimisation"]["priors"]["rho"]["uniform_lower"]
    rho_upper = options["optimisation"]["priors"]["rho"]["uniform_upper"]
    return [Uniform(Ne_lower, Ne_upper), Uniform(rho_lower, rho_upper)] 
end


function get_calculation_data()
    # TODO 
    # setup any additional metadata required for the calculations
    if parsed_args["opt_metric"]=="ld_decay" || parsed_args["opt_metric"]=="combined"
        # create a mapping to convert bp to cm distances
        bp_to_cm_df = DataFrame(CSV.File(parsed_args["distfile"]))
        bp_to_cm_df.BP = [parse(Int64,split(x,":")[2]) for x in bp_to_cm_df.Variant]
        bp_to_cm_map = Dict(zip(bp_to_cm_df.BP,bp_to_cm_df.Distance))
    end
end


function get_reference_statistics()
    # TODO
end


function run_simulation_rejection_abc()
    # TODO
end


function run_emulation_rejection_abc()
    # TODO
end


function run_abc(options, fp)

    @info "Setting up optimisation pipeline"

    priors = get_priors(options)
    @info @sprintf("Using priors %s", priors)

    get_calculation_data()

    get_reference_statistics()

    if options["optimisation"]["simulation_rejection_ABC"]["run"]
        run_simulation_rejection_abc()
    elseif options["optimisation"]["emulation_rejection_ABC"]["run"]
        run_emulation_rejection_abc()
    end
end


"""Entry point to running the optimisation pipeline
"""
# TODO move the chromosome if/else block to the run_program.jl file
function run_optimisation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            run_abc(options, fp)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        run_abc(options, fp)
    end
end
