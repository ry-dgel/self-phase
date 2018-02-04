include("self-phase.jl")
using YAML
using Iterators

function unpack(p)
    fname = @sprintf("%.0fnm_%.0fμJ_%.0fbar_%.0ffs_%.1fm_%.0ffs^2_%.0fμm",
                     p["λ"]*1E9, p["Energy"]*1E6, p["Pin"], p["Tfwhm"]*1E15, p["zmax"], p["Chirp"], p["fiberD"]*1E6)
    derive_constants(p)
    checkForData(fname, p)
    saveData(fname)
    simulate(E,p,zinit)
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

for jawn in param_tuples
    p = Dict{Any,Any}(zip(keys(lists), jawn))
    unpack(p)
end
