"""Code for setting up execution of the phenotype 
algorithm from the Julia pipeline
"""

using DelimitedFiles, Printf


"""There are two ways of running the phenotype program, either:
1. Directly specify a list of causal SNPs; or
2. Specify the polygenicity and pleiotropy parameters 
"""
function create_parfile(phenotype_options, filepaths, traw_prefix, out_prefix)
    lines = []
    push!(lines, @sprintf("nPopulation %s", phenotype_options["nPopulation"]))
    push!(lines, @sprintf("nTrait %s", phenotype_options["nTrait"]))
    push!(lines, @sprintf("a %s", phenotype_options["a"]))
    push!(lines, @sprintf("b %s", phenotype_options["b"]))
    push!(lines, @sprintf("c %s", phenotype_options["c"]))
    push!(lines, @sprintf("nComponent %s", phenotype_options["nComponent"]))
    push!(lines, @sprintf("PropotionGeno %s", phenotype_options["PropotionGeno"]))
    push!(lines, @sprintf("PropotionCovar %s", phenotype_options["PropotionCovar"]))
    push!(lines, @sprintf("SampleList %s", filepaths.phenotype_sample_list))
    
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
    push!(lines, @sprintf("GenoFile %s", traw_prefix))
    push!(lines, @sprintf("Output %s", out_prefix))
    push!(lines, @sprintf("CompWeight %s", phenotype_options["CompWeight"]))
    
    parfile = @sprintf("%s.parfile", out_prefix)

    open(parfile, "w") do io
        writedlm(io, lines)
    end
    
    return parfile
end


function convert_genotype_data(input, output, plink, memory)
    # convert from .bed to .traw format
    run(`$plink --bfile $input --recode A-transpose --memory $memory --out $output`)
end


function synthetic_pheno(filepaths, options, traw_prefix, out_prefix, seed)
    # create the parfile and generation synthetic phenotypes
    parfile = create_parfile(options["phenotype_data"], filepaths, traw_prefix, out_prefix)
    phenoalg = filepaths.phenoalg
    run(`$phenoalg $parfile $seed`)
end


function create_synthetic_phenotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    seed = options["global_parameters"]["random_seed"]
    memory = options["global_parameters"]["memory"]
    
    if chromosome == "all"
        for chromosome_i in 1:22
            filepaths = parse_filepaths(options, chromosome_i, superpopulation)
            input = filepaths.synthetic_data_prefix
            output = string(filepaths.synthetic_data_traw_prefix,"-",chromosome_i)
            convert_genotype_data(input, output, filepaths.plink, memory)
        end
    else
        filepaths = parse_filepaths(options, chromosome_i, superpopulation)
        input = filepaths.synthetic_data_prefix
        output = string(filepaths.synthetic_data_traw_prefix,"-",chromosome_i)
        convert_genotype_data(input, output, filepaths.plink, memory)
    end

    traw_prefix = filepaths.synthetic_data_traw_prefix
    out_prefix = filepaths.synthetic_data_prefix 
    synthetic_pheno(filepaths, options, traw_prefix, out_prefix, seed)
    
    @info @sprintf("Phenotype output is at %s.pheno{x}, where {x} is the phenotype number", out_prefix)

end
