"""Code for setting up execution of the phenotype 
algorithm from the Julia pipeline
"""

using DelimitedFiles

function convert_genotype_data(syndata, plink, combine, synth_paths, traw_prefix, memory)
    if combine
        # combine all chromosomes into one .traw file
        mergefile = @sprintf("%s_merge.txt", syndata)
        open(mergefile, "w") do io
            writedlm(io, synth_paths)
        end
        run(`$plink --bfile $syndata --merge-list $mergefile --recode A-transpose --memory $memory --out $traw_prefix`)
    else
        run(`$plink --bfile $syndata --recode A-transpose --memory $memory --out $syndata`)
    end
end


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


function synthetic_pheno(filepaths, options, seed, combine, synth_paths)
    if options["filepaths"]["phenotype"]["traw_prefix"] != "none"
        # traw file already exists 
        traw_prefix = options["filepaths"]["phenotype"]["traw_prefix"]
        filepaths.phenotype_sample_list = @sprintf("%s.sample", traw_prefix)
        parfile = create_parfile(options["phenotype_data"], filepaths, traw_prefix, filepaths.synthetic_data_traw_prefix)
    else
        traw_prefix = filepaths.synthetic_data_traw_prefix
        memory = options["global_parameters"]["memory"]
        convert_genotype_data(filepaths.synthetic_data_prefix, filepaths.plink, combine, synth_paths, traw_prefix, memory)
        if combine
            # if merging all chromosomes, need to give the pheno file a different name
            parfile = create_parfile(options["phenotype_data"], filepaths, traw_prefix, traw_prefix)
        else
            parfile = create_parfile(options["phenotype_data"], filepaths, filepaths.synthetic_data_prefix, filepaths.synthetic_data_prefix)
        end
    end
    
    phenoalg = filepaths.phenoalg
    run(`$phenoalg $parfile $seed`)
end


function create_synthetic_phenotype(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    seed = options["global_parameters"]["random_seed"]

    # synthetic genotype files are generated chromosome-by-chromosome
    if chromosome == "all"
        # creating phenotype based on genotype for all chromosomes
        chromosome_i = 1
        fp = parse_filepaths(options, chromosome_i, superpopulation)
        all_synth_paths = []
        for chromosome_i in 1:22
            fp_tmp = parse_filepaths(options, chromosome_i, superpopulation)
            push!(all_synth_paths, fp_tmp.synthetic_data_prefix)
        end
        synthetic_pheno(fp, options, seed, true, all_synth_paths)
    else
        # creating phenotype based on genotype for a single chromosome
        fp = parse_filepaths(options, chromosome, superpopulation)
        synthetic_pheno(fp, options, seed, false, [])
    end
end
