@everywhere include("self-phase.jl")
using YAML
using Iterators
#=
"""
To run, simply use:
```
julia multi-thread.jl <parameter file> <#saves>
```
e.g:
```
julia multi-thread.jl params.yaml 10
```
"""

"""
    unpack(p, numSaves)

given a dictionary of inital paremeters, and a number of times to save, performs
the entire simulation with those parameters.
"""
=#
@everywhere function unpack(p, numSaves)
    # Generate folder string
    fname = @sprintf("%.0fnm_%.0fμJ_%.2fbar_%.0ffs_%.1fm_%.0ffs^2_%.0fμm",
                     p["λ"]*1E9, p["Energy"]*1E6, p["Pin"], p["Tfwhm"]*1E15,
                     p["zmax"], p["Chirp"], p["fiberD"]*1E6)
    # Derive additional constants and spatial grids
    derive_constants(p)
    # Initialize electric field
    E, zinit = initialize(fname, p, "resume" in ARGS, "keep" in ARGS)
    # Save inital data
    saveData(fname, E, 0, zinit)
    # Simulate the whole thing
    simulate(E, p, zinit, fname, numSaves)
end

# Load paramters which are possible list
lists = YAML.load(open(ARGS[1]))
# Rename lambda to λ in dictionary.
if haskey(lists, "lambda")
    merge!(lists, Dict("λ" => pop!(lists, "lambda")))
end

# If paramter is a list of form [initial, final, stepsize] convert to list
# of all paramters in range. Otherwise keep it as a single value.
for key in keys(lists)
    if length(lists[key]) != 1
        lists[key] = collect(lists[key][1]:lists[key][3]:lists[key][2])
    else
        lists[key] = lists[key][1]
    end
end

# Generate all possible combinations of parameters.
param_tuples = collect(product(values(lists)...))

# Put all results into a single folder named data.
if "data" ∉ readdir()
    mkdir("data")
end
# Put all processes in same place.
@everywhere cd("data")

if length(ARGS) < 2
    println("No number of saves passed, assuming 2.")
    numSaves = 2
else
    numSaves = parse(Int, ARGS[2])
end

# Parallel map over all paramters to simulate!
pmap(jawn -> unpack(Dict{Any,Any}(zip(keys(lists), jawn)), numSaves), param_tuples)
