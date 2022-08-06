include("../evaluation/metrics/eval_gwas.jl")

using ArgParse

"""Script for making gwas summary statistics
"""

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--plink2"
            help = "path to plink2 software tool"
            arg_type = String
            required = true
        "--geno_prefix"
            help = "path to synthetic dataset"
            arg_type = String
            required = true
        "--pheno_prefix"
            help = "path to phenotype files"
            arg_type = String
            required = true
        "--pca_prefix"
            help = "path to pca output"
            arg_type = String
            required = true
        "--outprefix"
            help = "location for output"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()

    trait_idx = 1
    pcafile = string(parsed_args["pca_prefix"], ".eigenvec")
    run_gwas_tools(parsed_args["plink2"], parsed_args["geno_prefix"], parsed_args["pheno_prefix"], trait_idx, pcafile, parsed_args["outprefix"])

end

main() 

