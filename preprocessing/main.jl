using ArgParse, YAML

include("run_preprocessing.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            help = "path to the YAML configuration file"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()

    options = YAML.load_file(parsed_args["config"])
    run_preprocessing(options)
end


main()
