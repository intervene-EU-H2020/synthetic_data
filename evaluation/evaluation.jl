"""Parses input options and runs evaluation pipeline
"""

function run_kinship_evaluation()
    # TODO
end


function run_aats_evaluation()
    # TODO
end


function run_ld_evaluation()
    # TODO
end


function run_maf_evaluation()
    # TODO
end


function run_pca_evaluation()
    # TODO 
end


function run_external_tools(options, reference, filepaths)
    plink = filepaths.plink
    plink2 = filepaths.plink2
    king = filepaths.king

    reffile_prefix = reference
    synfile_prefix = filepaths.synthetic_data_prefix
    reffile_bed = @sprintf("%s.bed", reference)
    synfile_bed = @sprintf("%s.bed", filepaths.synthetic_data_prefix)

    if options["evaluation"]["metrics"]["maf"]
        reffile_out = @sprintf("%s.ref.maf", filepaths.evaluation_output)
        synfile_out = @sprintf("%s.syn.maf", filepaths.evaluation_output)
        run(`$plink --bfile $reffile_prefix --freq --out $reffile_out`)
        run(`$plink --bfile $synfile_prefix --freq --out $synfile_out`)
    elseif options["evaluation"]["metrics"]["pca"]
        reffile_out = @sprintf("%s.ref.pca", filepaths.evaluation_output)
        synfile_out = @sprintf("%s.syn.pca", filepaths.evaluation_output)
        run(`$plink2 --bfile $reffile_prefix --freq counts --pca allele-wts --out $reffile_out`)
        run(`$plink2 --bfile $synfile_prefix --freq counts --pca allele-wts --out $synfile_out`)
    elseif options["evaluation"]["metrics"]["kinship"] || options["evaluation"]["metrics"]["aats"]
        reffile_out = @sprintf("%s.ref.king", filepaths.evaluation_output)
        synfile_out = @sprintf("%s.syn.king", filepaths.evaluation_output)
        crossfile_out = @sprintf("%s.cross.king", filepaths.evaluation_output)
        run(`$king -b $reffile_bed --ibs --prefix $reffile_out`)
        run(`$king -b $synfile_bed --ibs --prefix $synfile_out`)
        run(`$king -b $reffile_bed,$synfile_bed --ibs --prefix $crossfile_out`)
    end
end


function run_pipeline(options, chromosome, superpopulation)
    fp = parse_filepaths(options, chromosome, superpopulation)
    reference = get_reference_data(fp)

    run_external_tools(options, reference, filepaths)

    if options["evaluation"]["metrics"]["aats"]
        run_aats_evaluation()
    elseif options["evaluation"]["metrics"]["kinship"]
        run_kinship_evaluation()
    elseif options["evaluation"]["metrics"]["ld"]
        run_ld_evaluation()
    elseif options["evaluation"]["metrics"]["maf"]
        run_maf_evaluation()
    elseif options["evaluation"]["metrics"]["pca"]
        run_pca_evaluation()
    end
end


function get_reference_data(filepaths)
    # TODO this dataset needs to be created
    return filepaths.evaluation_reference
end


"""Entry point to running the evaluation pipeline
"""
function run_evaluation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        for chromosome_i in 1:22
            run_pipeline(options, chromosome_i, superpopulation)
        end
    else
        run_pipeline(options, chromosome, superpopulation)
    end
end
