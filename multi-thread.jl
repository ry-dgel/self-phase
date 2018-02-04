@everywhere include("self-phase.jl")
using YAML
using Iterators

@everywhere function unpack(p, numSaves)
    fname = @sprintf("%.0fnm_%.0fμJ_%.0fbar_%.0ffs_%.1fm_%.0ffs^2_%.0fμm",
                     p["λ"]*1E9, p["Energy"]*1E6, p["Pin"], p["Tfwhm"]*1E15, p["zmax"], p["Chirp"], p["fiberD"]*1E6)
    derive_constants(p)
    E, zinit = initialize(fname, p, "resume" in ARGS, "keep" in ARGS)
    saveData(fname, E, 0, 0, zinit)
    simulate(E, p, zinit, fname, numSaves)
end

lists = YAML.load(open(ARGS[1]))
if haskey(lists, "lambda")
    merge!(lists, Dict("λ" => pop!(lists, "lambda")))
end
for key in keys(lists)
    if length(lists[key]) != 1
        lists[key] = collect(lists[key][1]:lists[key][3]:lists[key][2])
    else
        lists[key] = lists[key][1]
    end
end

param_tuples = collect(product(values(lists)...))

if "data" ∉ readdir()
    mkdir("data")
end
@everywhere cd("data")
numSaves = parse(Int, ARGS[2])
pmap(jawn -> unpack(Dict{Any,Any}(zip(keys(lists), jawn)), numSaves), param_tuples)
