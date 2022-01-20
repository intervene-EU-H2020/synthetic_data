"""Code for setting up execution of the phenotype 
algorithm from the Julia pipeline
"""

using DelimitedFiles

function convert_genotype_data(syndata, plink)
    run(`$plink --bfile $syndata --recode A-transpose --out $syndata`)
end


function create_parfile(options, filepaths)
    lines = []
    push!(lines, @sprintf("nPopulation %s", options["phenotype_data"]["nPopulation"]))
    push!(lines, @sprintf("nTrait %s", options["phenotype_data"]["nTrait"]))
    push!(lines, @sprintf("a %s", options["phenotype_data"]["a"]))
    push!(lines, @sprintf("b %s", options["phenotype_data"]["b"]))
    push!(lines, @sprintf("c %s", options["phenotype_data"]["c"]))
    push!(lines, @sprintf("nComponent %s", options["phenotype_data"]["nComponent"]))
    push!(lines, @sprintf("PropotionGeno %s", options["phenotype_data"]["PropotionGeno"]))
    push!(lines, @sprintf("PropotionCovar %s", options["phenotype_data"]["PropotionCovar"]))
    push!(lines, @sprintf("Polygenicity %s", options["phenotype_data"]["Polygenicity"]))
    push!(lines, @sprintf("Pleiotropy %s", options["phenotype_data"]["Pleiotropy"]))
    push!(lines, @sprintf("TraitCorr %s", options["phenotype_data"]["TraitCorr"]))
    push!(lines, @sprintf("PopulationCorr %s", options["phenotype_data"]["PopulationCorr"]))
    push!(lines, @sprintf("CompWeight %s", options["phenotype_data"]["CompWeight"]))
    push!(lines, @sprintf("CausalList %s", filepaths.phenotype_causal_list))
    push!(lines, @sprintf("SampleList %s", filepaths.phenotype_sample_list))
    push!(lines, @sprintf("Reference %s", filepaths.phenotype_reference))
    push!(lines, @sprintf("GenoFile %s", filepaths.synthetic_data_prefix))
    push!(lines, @sprintf("Output %s", filepaths.synthetic_data_prefix))

    parfile = @sprintf("%s.parfile", filepaths.synthetic_data_prefix)

    open(parfile, "w") do io
        writedlm(io, lines)
    end

    return parfile
end


function create_synthetic_phenotype_for_chromosome(filepaths, options, seed)
    convert_genotype_data(filepaths.synthetic_data_prefix, filepaths.plink)
    parfile = create_parfile(options, filepaths)
    phenoalg = filepaths.phenoalg
    run(`$phenoalg $parfile $seed`)
end


function create_synthetic_phenotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    seed = options["global_parameters"]["random_seed"]

    # synthetic genotype files are generated chromosome-by-chromosome
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            create_synthetic_phenotype_for_chromosome(fp, options, seed)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        create_synthetic_phenotype_for_chromosome(fp, options, seed)
    end
end
