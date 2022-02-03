"""
    Run a fast GWAS and generate manhattan and qqplot
"""

using CSV
using DataFrames
using MendelPlots
using Printf


function run_gwas_tools(plink2, syngeno_prefix, synpheno_prefix, trait_idx, covar, outdir)
	# TODO still giving run errors, see e.g. https://docs.julialang.org/en/v1/manual/running-external-programs/
	syngeno_prefix_fam = @sprintf("%s.fam", syngeno_prefix)
	syngeno_prefix_bed = @sprintf("%s.bed", syngeno_prefix)
	syngeno_prefix_bim = @sprintf("%s.bim", syngeno_prefix)
	synpheno_prefix_pheno_trait_idx = @sprintf("%s.pheno%i", synpheno_prefix, trait_idx)
	syngeno_prefix_phe_trait_idx = @sprintf("%s.phe%i", syngeno_prefix, trait_idx)
	outdir_tmp = @sprintf("%s.tmp", outdir)
	outdur_gwas = @sprintf("%s.PHENO1.glm.linear", outdir)
	outdir_sumstat = @sprintf("%s.sumstat", outdir)

	@info  "Create plink phenotype file"
	run(`paste $syngeno_prefix_fam \< \(awk 'NR\>1\{print $5\}' $synpheno_prefix_pheno_trait_idx\) \| awk '\{print $1,$2,$3,$4,$5,$7\}' \> $syngeno_prefix_phe_trait_idx`)
	
	@info  "GWAS using plink 2"
	run(`$plink2 --bed $syngeno_prefix_bed --bim $syngeno_prefix_bim --fam $syngeno_prefix_phe_trait_idx --glm hide-covar --covar $covar --ci 0.95 --out $outdir`)
	
	@info  "Create summary statistics"
	run(`echo -e "CHR\tBP\tSNP\tA2\tA1\tA1\tTEST\tNMISS\tBETA\tSE\tL95\tU95\tSTAT\tP" \> $outdir_tmp`)
	run(pipeline(`awk 'NR\>1\{print $0\}' $outdur_gwas`, `grep -v NA`, `sed 's/ /\t/g' \>\> $outdir_tmp`))
	run(`cut -f 1-5,7- $outdir_tmp \> $outdir_sumstat \&\& rm $outdir_tmp`)
end


function plot_gwas(gwas_out_prefix, outdir)
	df = CSV.File(gwas_out_prefix * ".sumstat", normalizenames=true) |> DataFrame
	qq_out = @sprintf("%s.qq.png", outdir)
	man_out = @sprintf("%s.man.png", outdir)
	
	@info  "Plot qq plot"
	qq(df[:,"P"], dpi = 400, fontsize = 10pt, ystep = 5, outfile = qq_out)
	
	@info  "Plot manhattan plot"
	manhattan(df[:,"P"], df[:,"CHR"], df[:,"BP"], dpi = 400, fontsize = 10pt, ystep = 5, outfile = man_out)
end


"""Run GWAS for each phenotypic trait and plot results
"""
function run_gwas(ntraits, plink2, covar, synt_data_prefix, eval_dir)
	for trait_idx in 1:ntraits
		outdir = @sprintf("%s.trait%i", eval_dir, trait_idx)
		run_gwas_tools(plink2, synt_data_prefix, synt_data_prefix, trait_idx, covar, outdir)
		plot_gwas(outdir, outdir)
	end
end