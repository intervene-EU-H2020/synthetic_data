using CSV, DataFrames, StatsBase, DelimitedFiles

include("../evaluation/metrics/eval_gwas.jl")


"""
Split the geno/pheno datasets into train and test sets (50/50 ratio)
"""
function split_train_test(data_prefix, prspipe_outdir, plink, num_pheno)
    # select samples to retain in train/set sets
    samples_df = CSV.File(string(data_prefix, ".fam"), normalizenames=true, header = 0) |> DataFrame
    samples1 = sort(sample(1:nrow(samples_df), Int(ceil(nrow(samples_df)/2)); replace=false))
    samples2 = setdiff(1:nrow(samples_df), samples1)

    # create keepfiles for train/test sets
    traindir = string(prspipe_outdir, "/train")
    testdir = string(prspipe_outdir, "/test")
    for outdir in [traindir, testdir]
        if !isdir(outdir)
            mkpath(outdir)
        end
    end
    
    keepfile_train = string(traindir, "/keep.txt")
    keepfile_test = string(testdir, "/keep.txt")
    keeplist_train = samples_df[samples1, [:Column1, :Column2]]
    keeplist_test =  samples_df[samples2, [:Column1, :Column2]]
    CSV.write(keepfile_train, keeplist_train, header=false, delim="\t")
    CSV.write(keepfile_test, keeplist_test, header=false, delim="\t")
    
    # extract train/test genotypes 
    for chromosome in 1:22
        outfile_train = string(traindir, "/train_chr", chromosome)
        outfile_test = string(testdir, "/test_chr", chromosome)
        run(`$plink --bfile $data_prefix --keep $keepfile_train --make-bed --out $outfile_train`)
        run(`$plink --bfile $data_prefix --keep $keepfile_test --make-bed --out $outfile_test`)
    end
    
    # extract train/test phenotypes 
    for pheno in 1:num_pheno
        pheno_df = CSV.File(string(data_prefix, ".pheno", pheno), header = 1) |> DataFrame
        pheno_train = pheno_df[samples1, :]
        pheno_test = pheno_df[samples2, :]
        pheno_train_file = string(traindir, "/train.pheno", pheno)
        pheno_test_file = string(testdir, "/test.pheno", pheno)
        CSV.write(pheno_train_file, pheno_train, delim="\t")
        CSV.write(pheno_test_file, pheno_test, delim="\t")
    end
end


"""
Compute summary statistics on the train set
"""
function compute_summary_stats(prspipe_train_outdir, plink, plink2, num_pheno)
    # prepare pca input
    pcafile = string(prspipe_train_outdir, "/train.pca")
    mergelist = [string(prspipe_train_outdir, "/train_chr", x) for x in 1:22]
    mergefile = string(prspipe_train_outdir, "/merge.txt")
    open(mergefile, "w") do io
        writedlm(io, mergelist)
    end
    train_data_prefix = string(prspipe_train_outdir, "/train")
    first_train_data_prefix = mergelist[1]
    run(`$plink --bfile $first_train_data_prefix --merge-list $mergefile --make-bed --out $train_data_prefix`)
    run(`$plink2 --bfile $train_data_prefix --freq counts --pca allele-wts --out $pcafile`)
    run_gwas(num_pheno, plink2, string(pcafile, ".eigenvec"), train_data_prefix, string(prspipe_train_outdir, "/train"))
end


"""
Format summary statistics and test data (geno/pheno) for input to prspipe

The final files that should be transfered are 
test/test_chr{chr}.{bed/bim/fam} for the genotype,
test/test.pheno{num}.tsv for the phenotype, and 
train/train.trait{num}.PHENO1.glm.linear.final for the summary statistics
"""
function format_outputs(prspipe_train_outdir, prspipe_test_outdir, num_pheno)
    for pheno in 1:num_pheno
        # format the phenotype data (genodata data doesn't change)
        phenofile = string(prspipe_test_outdir, "/test.pheno", pheno)
        pheno_df = CSV.File(phenofile) |> DataFrame
        pheno_df.Sample1 = [split(x, "_")[1] for x in pheno_df.Sample]
        pheno_df.Sample2 = [split(x, "_")[2] for x in pheno_df.Sample]
        pheno_df = pheno_df[:, [:Sample1, :Sample2, :Phenotype]]
        CSV.write(string(phenofile, ".tsv"), pheno_df, writeheader=false, delim="\t")
        # format the gwas summary statistics
        gwasfile = string(prspipe_train_outdir, "/train.trait", pheno, ".PHENO1.glm.linear")
        gwas_df = CSV.File(gwasfile) |> DataFrame
        gwas_df = gwas_df[:, [1, 2, 3, 4, 5, 8, 9, 10, 14]]
        rename!(gwas_df, [:CHR, :POS, :SNP, :A2, :A1, :N, :BETA, :SE, :P])
        CSV.write(string(gwasfile, ".final"), gwas_df, writeheader=true, delim="\t")
    end
end


"""Prepares the data for input to INTERVENE's polygenic risk scoring pipeline
"""
function run_prspipe_prep(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)
    
    if chromosome == "all"
        chromosome_i = 1
        fp = parse_filepaths(options, chromosome_i, superpopulation)
        num_pheno = options["phenotype_data"]["nTrait"]
        split_train_test(fp.synthetic_data_traw_prefix, fp.prspipe_dir, fp.plink, num_pheno)
        compute_summary_stats(string(fp.prspipe_dir, "/train"), fp.plink, fp.plink2, num_pheno)
        format_outputs(string(fp.prspipe_dir, "/train"), string(fp.prspipe_dir, "/test"), num_pheno)
    else
        throw(error("Config error: the chromosome `all` option and data for all 22 chromosomes is required for prspipe"))
    end
end
