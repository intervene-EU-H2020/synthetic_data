"""Executes the pipeline for evaluating synthetic data quality
"""

include("../utils/reference_data.jl")
include("metrics/eval_aats.jl")
include("metrics/eval_kinship_detail.jl")
include("metrics/eval_ld_corr.jl")
include("metrics/eval_ld_decay.jl")
include("metrics/eval_maf.jl")
include("metrics/eval_pca.jl")
include("metrics/eval_gwas.jl")

function run_kinship_evaluation(ibsfile_real, ibsfile_synt, ibsfile_cross)
    run_kinship(ibsfile_real, ibsfile_synt, ibsfile_cross)
end


function run_aats_evaluation(ibsfile_cross)
    run_aats(ibsfile_cross)
end


function run_ld_corr_evaluation(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
    run_ld_corr(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
end


function run_ld_decay_evaluation(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map)
    run_ld_decay(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map)
end


function run_maf_evaluation(real_maf_file, synt_maf_file)
    run_maf(real_maf_file, synt_maf_file)
end


function run_pca_evaluation(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_prefix, synt_data_prefix, eval_dir)
    real_data_pop = string(real_data_prefix, ".sample")
    synt_data_pop = string(synt_data_prefix, ".sample")
    run_pca(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_pop, synt_data_pop, eval_dir)
end


function run_gwas_evaluation(ntraits, plink2, covar, synt_data_prefix, eval_dir)
    run_gwas(ntraits, plink2, covar, synt_data_prefix, eval_dir)
end


"""Computations for MAF using PLINK
"""
function run_maf_tools(plink, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for MAF"
    reffile_out = @sprintf("%s.ref.maf", outdir)
    synfile_out = @sprintf("%s.syn.maf", outdir)
    run(`$plink --bfile $reffile_prefix --freq --out $reffile_out`)
    run(`$plink --bfile $synfile_prefix --freq --out $synfile_out`)
    real_maffile = @sprintf("%s.frq", reffile_out)
    syn_maffile = @sprintf("%s.frq", synfile_out)
    return real_maffile, syn_maffile
end


"""Computations for PCA using PLINK
"""
function run_pca_tools(plink2, king, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for PCA"
    reffile_out = @sprintf("%s.ref.pca", outdir)
    synfile_out = @sprintf("%s.syn.pca", outdir)
    run(`$plink2 --bfile $reffile_prefix --freq counts --pca allele-wts --out $reffile_out`)
    run(`$plink2 --bfile $synfile_prefix --freq counts --pca allele-wts --out $synfile_out`)

    # For PCA projection
    projfile_prefix = @sprintf("%s_proj", outdir)
    projfile_out = @sprintf("%spc.txt", projfile_prefix)
    run(`$king -b $reffile_prefix.bed,$synfile_prefix.bed --pca --projection --prefix $projfile_prefix`)

    return reffile_out, synfile_out, projfile_out
end


"""Computations for relatedness using KING
"""
function run_relatedness_tools(king, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for relatedness"
    reffile_bed = @sprintf("%s.bed", reffile_prefix)
    synfile_bed = @sprintf("%s.bed", synfile_prefix)
    reffile_out = @sprintf("%s.ref.king", outdir)
    synfile_out = @sprintf("%s.syn.king", outdir)
    crossfile_out = @sprintf("%s.cross.king", outdir)
    run(`$king -b $reffile_bed --ibs --prefix $reffile_out`)
    run(`$king -b $synfile_bed --ibs --prefix $synfile_out`)
    run(`$king -b $reffile_bed,$synfile_bed --ibs --prefix $crossfile_out`)
    real_ibsfile = @sprintf("%s.ibs0", reffile_out)
    syn_ibsfile = @sprintf("%s.ibs0", synfile_out)
    cross_ibsfile = @sprintf("%s.ibs0", crossfile_out)
    return real_ibsfile, syn_ibsfile, cross_ibsfile
end


"""Runs computations using external tools such as PLINK and KING, based on the selection of evaluation metrics
"""
function run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)
    external_files = Dict()
    if metrics["maf"]
        external_files["real_maffile"], external_files["syn_maffile"] = run_maf_tools(filepaths.plink, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["pca"] || metrics["gwas"]
        external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"] = run_pca_tools(
            filepaths.plink2, 
            filepaths.king,
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["kinship"] || metrics["aats"]
        external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"] = run_relatedness_tools(filepaths.king, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    return external_files
end


"""Executes evaluation for the metrics specified in the configuration file
"""
function run_pipeline(options, chromosome, superpopulation, metrics)
    compute_chromosome = chromosome=="all" ? 1 : chromosome
    filepaths = parse_filepaths(options, compute_chromosome, superpopulation)
    genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)

    reffile_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, compute_chromosome)
    synfile_prefix = chromosome=="all" ? filepaths.synthetic_data_traw_prefix : filepaths.synthetic_data_prefix

    external_files = run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)

    if metrics["aats"]
        run_aats_evaluation(external_files["cross_ibsfile"])
    end
    if metrics["kinship"]
        run_kinship_evaluation(external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"])
    end
    if metrics["ld_corr"]
        run_ld_corr_evaluation(reffile_prefix, synfile_prefix, filepaths.evaluation_output, filepaths.plink)
    end
    if metrics["ld_decay"]
        bp_to_cm_map = create_bp_cm_ref(filepaths.genetic_distfile)
        run_ld_decay_evaluation(synfile_prefix, reffile_prefix, filepaths.plink, filepaths.mapthin, filepaths.evaluation_output, bp_to_cm_map)
    end
    if metrics["maf"]
        run_maf_evaluation(external_files["real_maffile"],  external_files["syn_maffile"])
    end
    if metrics["pca"]
        run_pca_evaluation(
            external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"],
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["gwas"]
        run_gwas_evaluation(options["phenotype_data"]["nTrait"], filepaths.plink2, @sprintf("%s.eigenvec", external_files["syn_pcafile"]), synfile_prefix, filepaths.evaluation_output)
    end
end


"""Entry point to running the evaluation pipeline for genotype data

Note that the evaluation pipeline assumes that the synthetic data you want 
to evaluate has already been geneerated, using the setup specified in the
configuration file. It is therefore recommended to run the pipeline with the 
--genotype and --evaluation flags together, so that the program generates 
the data and then immediately evaluates it using the correct settings.
"""
function run_evaluation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)

    metrics = options["evaluation"]["metrics"]
    if chromosome == "all"
        if metrics["gwas"]
            # TODO gwas is the only metric implemented for multi-chromosome analysis
            new_metrics = copy(metrics)
            for key in keys(new_metrics)
                new_metrics[key] = false
            end
            new_metrics["gwas"] = true
            run_pipeline(options, chromosome, superpopulation, new_metrics)
        end
        
        for chromosome_i in 1:22
            metrics["gwas"] = false
            run_pipeline(options, chromosome_i, superpopulation, metrics)
        end
    else
        run_pipeline(options, chromosome, superpopulation, metrics)
    end
end
