"""Run likelihood-free inference for optimising model parameters
"""

using GpABC

include("../evaluation/metrics/eval_ld_decay.jl")
include("../evaluation/metrics/eval_kinship_quick.jl")
include("../utils/reference_data.jl")
include("../algorithms/genotype/genotype_algorithm.jl")


"""Struct for storing metadata required by the ABC procedure
"""
struct ABCMetadata
    statistics::Vector # list of summary statistics included in the optimisation 
    ref_prefix::String # prefix for the reference dataset used to generate synthetic data
    eval_ref_prefix::String # prefix for other reference used to compute summary statistics
    out_prefix::String # prefix for the synthetic dataset
    superpopulation::String # abbreviation for superpopulation
    nsamples_ref::Integer # number of samples in reference dataset
    nsamples_eval_ref::Integer # number of samples for evaluation reference dataset
    nsamples_syn::Integer # number of samples in synthetic dataset
    bp_to_cm_map::Dict # a mapping to convert bp to cm distances
    plink::String # plink software path
    king::String # king software path
    mapthin::String # mapthin software path
end


"""
Simulator function for the optimisation pipeline 

Arguments
- var_params::Vector{Float64}: list of the form [Ne_s, rho_s, mu_s]
"""
function simulator_function(var_params::AbstractArray{Float64,1})
    Ne_s = var_params[1]
    rho_s = var_params[2]
    
    # simulate a synthetic dataset with these parameter values and calculate the summary statistic
    genomic_metadata.population_Nes[abc_metadata.superpopulation] = trunc(Int, Ne_s)
    genomic_metadata.population_rhos[abc_metadata.superpopulation] = rho_s
    create_synthetic_genotype_for_chromosome(genomic_metadata) 
    sumstat = calculate_summary_statistic(abc_metadata.out_prefix, abc_metadata.nsamples_syn, abc_metadata.ref_prefix)
    
    return sumstat
end


"""Calculate summary statistics for ABC

Set use_cross=true if calculating kinship for the reference dataset (because there is only one dataset)
"""
function calculate_summary_statistic(syn_prefix, nsamples_syn, ref_prefix, cross_prefix=NaN, nsamples_cross=NaN, use_cross=false)
    if "ld_decay" ∈ abc_metadata.statistics
        ld_decay = LD_decay(syn_prefix, abc_metadata.plink, abc_metadata.mapthin, syn_prefix, abc_metadata.bp_to_cm_map)
        sumstat = ld_decay[:,2]
    end
    if "kinship" ∈ abc_metadata.statistics
        if use_cross
            sumstat_kin = kinship_cross(cross_prefix, nsamples_cross, ref_prefix, abc_metadata.king, syn_prefix)
        else
            sumstat_kin = kinship_cross(syn_prefix, nsamples_syn, ref_prefix, abc_metadata.king, syn_prefix)
        end
        
        if "ld_decay" ∈ abc_metadata.statistics
            sumstat = vcat(sumstat, sumstat_kin)
        else
            sumstat = sumstat_kin
        end
    end
    sumstat = reshape(sumstat,length(sumstat),1)
    return sumstat
end


"""Create priors for ABC
"""
function get_priors(options)
    # uses uniform priors
    Ne_lower = options["optimisation"]["priors"]["Ne"]["uniform_lower"]
    Ne_upper = options["optimisation"]["priors"]["Ne"]["uniform_upper"]
    rho_lower = options["optimisation"]["priors"]["rho"]["uniform_lower"]
    rho_upper = options["optimisation"]["priors"]["rho"]["uniform_upper"]
    return [Uniform(Ne_lower, Ne_upper), Uniform(rho_lower, rho_upper)] 
end


"""Run the simulation rejection ABC
"""
function run_simulation_rejection_abc(sumstat_reference, simulator_function, priors, sim_options)
    abc_result = SimulatedABCRejection(sumstat_reference, simulator_function, priors, sim_options["threshold"], sim_options["n_particles"], max_iter=sim_options["max_iter"], write_progress=sim_options["write_progress"], progress_every=1)
    return abc_result
end


"""Run the emulation rejection ABC
"""
function run_emulation_rejection_abc(sumstat_reference, simulator_function, priors, emu_options)
    abc_result = EmulatedABCRejection(sumstat_reference, simulator_function, priors, emu_options["threshold"], emu_options["n_particles"], emu_options["n_design_points"], max_iter=emu_options["max_iter"], write_progress=emu_options["write_progress"], progress_every=1)
    return abc_result
end


"""Save results from ABC procedure
"""
function save_results(abc_result, output_prefix, simulation_type)
    plot_ref = Plots.plot(abc_result)
    savefig(plot_ref, @sprintf("%s_%s.png", output_prefix, simulation_type))
    
    Nes = abc_result.population[:,1]
    rhos = abc_result.population[:,2]
    ds = abc_result.distances
    df = DataFrame(Ne=Nes, ρ=rhos, d=ds)
    CSV.write(@sprintf("%s_%s.csv", output_prefix, simulation_type), df)
end


"""Calculates the ABC summary statistics for the reference dataset
"""
function get_reference_statistics(ref_prefix, nsamples_ref)
    if "kinship" ∈ abc_metadata.statistics
        cross1_prefix, cross2_prefix, nsamples_cross1, nsamples_cross2 = create_cross_datasets(ref_prefix, abc_metadata.plink)
        sumstat = calculate_summary_statistic(ref_prefix, nsamples_ref, cross2_prefix, cross1_prefix, nsamples_cross1, true)
    else
        cross_prefix=NaN
        sumstat = calculate_summary_statistic(ref_prefix, nsamples_ref, cross_prefix)
    end
    return sumstat
end


"""Create the struct containing metadata for the ABC procedure
"""
function get_abc_metadata(options, superpopulation, filepaths, chromosome)
    bp_to_cm_map = create_bp_cm_ref(filepaths.genetic_distfile)
    statistics = [k for k in ["ld_decay", "kinship"] if options["optimisation"]["summary_statistics"][k]]
    ref_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, chromosome)
    ref1_prefix, ref2_prefix, nsamples_ref1, nsamples_ref2 = create_cross_datasets(ref_prefix, filepaths.plink)
    out_prefix = filepaths.optimisation_output
    nsamples_syn = options["genotype_data"]["samples"]["default"]["nsamples"]
    return ABCMetadata(statistics, ref1_prefix, ref2_prefix, out_prefix, superpopulation, nsamples_ref1, nsamples_ref2, nsamples_syn, bp_to_cm_map, filepaths.plink, filepaths.king, filepaths.mapthin)
end


"""Runs the Approximate Bayesian Computation (ABC) procedure for the specified superpopulation
"""
function run_abc_for_superpopulation(options, filepaths, superpopulation, chromosome)
    @info "Setting up optimisation pipeline"

    priors  = get_priors(options)
    @info @sprintf("Using priors %s", priors)
    
    global genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)
    genomic_metadata.outfile_prefix = filepaths.optimisation_output
    global abc_metadata = get_abc_metadata(options, superpopulation, filepaths, chromosome)

    @info "Computing summary statistics for reference data"
    sumstat_reference = get_reference_statistics(abc_metadata.eval_ref_prefix, abc_metadata.nsamples_eval_ref)
    
    @info "Running optimisation pipeline"
    t = @elapsed begin
        if options["optimisation"]["simulation_rejection_ABC"]["run"]
            sim_options = options["optimisation"]["simulation_rejection_ABC"]
            abc_result = run_simulation_rejection_abc(sumstat_reference, simulator_function, priors, sim_options)
            save_results(abc_result, abc_metadata.out_prefix, "sim")
        elseif options["optimisation"]["emulation_rejection_ABC"]["run"]
            emu_options = options["optimisation"]["emulation_rejection_ABC"]
            abc_result = run_emulation_rejection_abc(sumstat_reference, simulator_function, priors, emu_options)
            save_results(abc_result, abc_metadata.out_prefix, "emu")
        end
    end

    @info @sprintf("Optimisation completed in %.2f minutes", t/60)
end


"""Entry point to running the optimisation pipeline
"""
function run_optimisation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            run_abc_for_superpopulation(options, fp, superpopulation, chromosome_i)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        run_abc_for_superpopulation(options, fp, superpopulation, chromosome)
    end
end
