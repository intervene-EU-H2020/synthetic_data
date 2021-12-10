using ArgParse, YAML

include("preprocessing/preprocessing.jl")

"""Executes the program, according to which pipelines and configuration options are specified in the input

Note that any combination of pipelines can be run together, except the optimisation pipeline, which exits immediately after running

Also note that there is a specific ordering to the pipeline execution
"""
function run_program(pipelines, options)
    if pipelines["preprocessing"]
        @info "running the preprocessing pipeline"
        run_preprocessing(options)
    end

    if pipelines["optimisation"]
        @info "optimising model parameter values"
        # TODO
        exit(0)
    end

    if pipelines["genotype"]
        @info "generating synthetic genotype data"
        # TODO
    end

    if pipelines["phenotype"]
        @info "generating synthetic phenotype data"
        # TODO
    end

    if pipelines["evaluation"]
        @info "evaluating synthetic data quality"
        # TODO
    end
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "filepath to YAML configuration file"
            arg_type = String
            required = true
        "--preprocessing"
            help = "preprocesses raw data to use as input for synthetic data generation"
            action = :store_true
        "--genotype"
            help = "generate synthetic genotype data"
            action = :store_true
        "--phenotype"
            help = "generate synthetic phenotype data (requires synthetic genotype data input)"
            action = :store_true
        "--evaluation"
            help = "evaluate synthetic data quality (requires synthetic genotype data input)"
            action = :store_true
        "--optimisation"
            help = "run procedure for optimising model parameter values"
            action = :store_true
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()

    options = YAML.load_file(parsed_args["config"])

    pipelines = Dict("preprocessing" => parsed_args["preprocessing"], 
                     "genotype" => parsed_args["genotype"],
                     "phenotype" => parsed_args["phenotype"],
                     "evaluation" => parsed_args["evaluation"],
                     "optimisation" => parsed_args["optimisation"])

    println("Running pipelines:")
    for (arg,val) in pipelines
        println("  $arg  =>  $val")
    end

    run_program(pipelines, options)
    
end

main()
