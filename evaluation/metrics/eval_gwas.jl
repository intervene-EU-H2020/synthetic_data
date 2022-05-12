"""
    Run a fast GWAS and generate manhattan and qqplot
"""

using CSV
using DataFrames
using MendelPlots
using Printf


function run_gwas_tools(plink2, syngeno_prefix, synpheno_prefix, trait_idx, covar, outdir)
	@info  "Create plink phenotype file"
	fam = CSV.File(syngeno_prefix * ".fam", normalizenames=true, header = 0) |> DataFrame
	pheno = CSV.File(synpheno_prefix * ".pheno" * string(trait_idx), normalizenames=true) |> DataFrame
	fam[!,"Sample"] = string.(fam[:,1], "_",fam[:,2])
	
	PhenoFam = leftjoin(fam, pheno, on = :Sample)
	PhenoFam = PhenoFam[:, ["Column1","Column2","Column3","Column4","Column5","Phenotype_liability_"]]
	syngeno_prefix_phe_trait_idx = @sprintf("%s.phe%i", synpheno_prefix, trait_idx)
	CSV.write(syngeno_prefix_phe_trait_idx, DataFrame(PhenoFam), delim = "\t", header = false)
	
	@info  "GWAS using plink 2"
	syngeno_prefix_bed = @sprintf("%s.bed", syngeno_prefix)
	syngeno_prefix_bim = @sprintf("%s.bim", syngeno_prefix)
	run(`$plink2 --bed $syngeno_prefix_bed --bim $syngeno_prefix_bim --fam $syngeno_prefix_phe_trait_idx --glm hide-covar --covar $covar --ci 0.95 --out $outdir`)

	@info  "Create summary statistics"
	GWASout = CSV.File( outdir * ".PHENO1.glm.linear", normalizenames=true) |> DataFrame
	rename!(GWASout,[:CHR, :BP, :SNP, :A2, :A1, :A1_dup, :TEST, :NMISS, :BETA, :SE, :L95, :U95, :STAT, :P, :ERRCODE])
	select!(GWASout, Not(:A1_dup))

	CSV.write(outdir * ".sumstat", DataFrame(GWASout), delim = "\t")
end


function plot_gwas(gwas_out_prefix, outdir)
	# TODO appears to be a bug here
	df = CSV.File(gwas_out_prefix * ".sumstat", normalizenames=true) |> DataFrame
	df = df[df[:,"P"] .!= "NA", :] # drop nan
	df.P = parse.(Float64, df.P) # convert to float

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