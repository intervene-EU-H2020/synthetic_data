"""Code for setting up execution of the phenotype 
algorithm from the Julia pipeline
"""

using DelimitedFiles

function convert_genotype_data(syndata, plink)
    run(`$plink --bfile $syndata --recode A-transpose --out $syndata`)
end


function create_parfile(options, filepaths)
    lines = []
    push!(lines, @sprintf("Heritability %s", options["phenotype_data"]["heritability"]))
    push!(lines, @sprintf("Prevalence %s", options["phenotype_data"]["prevalence"]))
    push!(lines, @sprintf("a %s", options["phenotype_data"]["a"]))
    push!(lines, @sprintf("b %s", options["phenotype_data"]["b"]))
    push!(lines, @sprintf("c %s", options["phenotype_data"]["c"]))
    push!(lines, @sprintf("nComponent %s", options["phenotype_data"]["n_component"]))
    push!(lines, @sprintf("CompWeight %s", options["phenotype_data"]["comp_weight"]))
    push!(lines, @sprintf("PopCovar %s", options["phenotype_data"]["pop_covar"]))
    push!(lines, @sprintf("Polygenicity %s", options["phenotype_data"]["polygenicity"]))
    push!(lines, @sprintf("CausalList %s", options["filepaths"]["causal_list"])) # TODO are filepaths chromosome or population specific?
    push!(lines, @sprintf("SampleList %s", options["filepaths"]["sample_list"]))
    push!(lines, @sprintf("Reference %s", options["filepaths"]["pheno_reference"]))
    push!(lines, @sprintf("GenoFile %s", @sprintf("%s.traw", filepaths.synthetic_data_prefix)))
    push!(lines, @sprintf("Output %s", filepaths.synthetic_data_prefix))

    parfile = @sprintf("%s.parfile", filepaths.synthetic_data_prefix)

    open(parfile, "w") do io
        writedlm(io, lines)
    end

    return parfile
end


function create_synthetic_phenotype_for_chromosome(filepaths, seed)
    convert_genotype_data(filepaths.synthetic_data_prefix, filepaths.plink)
    parfile = create_parfile(options, filepaths)
    phenoalg = filepaths.phenoalg # TODO add phenoalg to filepaths
    run(`$phenoalg $parfile $seed`)
end


# TODO is this chromosome specific?
function create_synthetic_phenotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    seed = options["global_parameters"]["random_seed"]

    # synthetic genotype files are generated chromosome-by-chromosome
    if chromosome == "all"
        for chromosome_i in 1:22
            fp = parse_filepaths(options, chromosome_i, superpopulation)
            create_synthetic_phenotype_for_chromosome(fp, seed)
        end
    else
        fp = parse_filepaths(options, chromosome, superpopulation)
        create_synthetic_phenotype_for_chromosome(fp, seed)
    end
end
