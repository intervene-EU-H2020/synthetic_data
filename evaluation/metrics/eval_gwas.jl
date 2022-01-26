"""
    Run a fast GWAS and generate manhattan and qqplot
"""

using CSV
using DataFrames
using MendelPlots
using Printf

function run_gwas_tools(plink2, syngeno_prefix, synpheno_prefix, trait_idx, covar, outdir)
	@info  "Create plink phenotype file"
	run(`paste ${syngeno_prefix}.fam < (awk 'NR>1{print $5}' ${synpheno_prefix}.pheno${trait_idx}) | awk '{print $1,$2,$3,$4,$5,$7}' > $syngeno_prefix.phe${trait_idx}`)
	
	@info  "GWAS using plink 2"
	run(`$plink --bed ${syngeno_prefix}.bed --bim ${syngeno_prefix}.bim --fam ${syngeno_prefix}.phe${trait_idx} --glm hide-covar --covar ${covar} --ci 0.95 --out ${outdir}`)
	
	@info  "Create summary statistics"
	run(`echo -e "CHR\tBP\tSNP\tA2\tA1\tA1\tTEST\tNMISS\tBETA\tSE\tL95\tU95\tSTAT\tP" > ${outdir}.tmp`)
	run(`awk 'NR>1{print $0}' ${outdir}.PHENO1.glm.linear | grep -v NA | sed 's/ /\t/g' >> ${outdir}.tmp`)
	run(`cut -f 1-5,7- ${outdir}.tmp > ${outdir}.sumstat && rm ${outdir}.tmp`)
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
