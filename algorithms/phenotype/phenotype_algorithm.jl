"""Code for setting up execution of the phenotype 
algorithm from the Julia pipeline
"""

using DelimitedFiles, Printf


"""There are two ways of running the phenotype program, either:
1. Directly specify a list of causal SNPs; or
2. Specify the polygenicity and pleiotropy parameters 
"""
function create_parfile(phenotype_options, filepaths, genetics_prefix, out_prefix, chromosome)
    lines = []
    push!(lines, @sprintf("nPopulation %s", phenotype_options["nPopulation"]))
    push!(lines, @sprintf("nTrait %s", phenotype_options["nTrait"]))
    push!(lines, @sprintf("a %s", phenotype_options["a"]))
    push!(lines, @sprintf("b %s", phenotype_options["b"]))
    push!(lines, @sprintf("c %s", phenotype_options["c"]))
    push!(lines, @sprintf("nComponent %s", phenotype_options["nComponent"]))
    push!(lines, @sprintf("PropotionGeno %s", phenotype_options["PropotionGeno"]))
    push!(lines, @sprintf("PropotionCovar %s", phenotype_options["PropotionCovar"]))
    push!(lines, @sprintf("SampleList %s-%s.sample", genetics_prefix, chromosome))
    
    if phenotype_options["Causality"]["UseCausalList"]
        push!(lines, @sprintf("CausalList %s", filepaths.phenotype_causal_list))
    else
        push!(lines, @sprintf("Polygenicity %s", phenotype_options["Causality"]["Polygenicity"]))
        push!(lines, @sprintf("Pleiotropy %s", phenotype_options["Causality"]["Pleiotropy"]))
    end

    
    push!(lines, @sprintf("TraitCorr %s", phenotype_options["TraitCorr"]))
    push!(lines, @sprintf("PopulationCorr %s", phenotype_options["PopulationCorr"]))
    push!(lines, @sprintf("Prevalence %s", phenotype_options["Prevalence"]))
    push!(lines, @sprintf("Reference %s", filepaths.phenotype_reference))
    push!(lines, @sprintf("GenoFile %s", genetics_prefix))
    push!(lines, @sprintf("Output %s", out_prefix))
    push!(lines, @sprintf("CompWeight %s", phenotype_options["CompWeight"]))
    
    parfile = @sprintf("%s.parfile", out_prefix)

    open(parfile, "w") do io
        writedlm(io, lines)
    end
    
    return parfile
end


function synthetic_pheno(filepaths, options, genetics_prefix, out_prefix, seed, chromosome)
    # create the parfile and generation synthetic phenotypes
    parfile = create_parfile(options["phenotype_data"], filepaths, genetics_prefix, out_prefix, chromosome)
    phenoalg = filepaths.phenoalg
    run(`$phenoalg $parfile $seed`)
end


function create_synthetic_phenotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    seed = options["global_parameters"]["random_seed"]

    if chromosome == "all"
        # use filepaths for the 1st chromosome
        tmp_chromosome = 1
    else
        tmp_chromosome = chromosome
    end

    filepaths = parse_filepaths(options, tmp_chromosome, superpopulation)
    out_prefix = filepaths.synthetic_data_prefix[1:end-ndigits(tmp_chromosome)-1]

    if options["filepaths"]["phenotype"]["plink_override"] == "none"
        genetics_prefix = out_prefix
    else
        genetics_prefix = options["filepaths"]["phenotype"]["plink_override"]
    end

    synthetic_pheno(filepaths, options, genetics_prefix, out_prefix, seed, tmp_chromosome)
    
    @info @sprintf("Phenotype output is at %s.pheno{x}, where {x} is the phenotype number", out_prefix)

end
