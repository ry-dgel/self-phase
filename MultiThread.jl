@everywhere push!(LOAD_PATH, pwd())
@everywhere using SelfPhase
@everywhere using Printf
using Distributed
using Base.Iterators: product
using ProgressMeter
using YAML
#=
To run, simply use:
```
julia multi-thread.jl <parameter file> <#saves>
```
e.g:
```
julia multi-thread.jl params.yaml 10
```
=#


@everywhere function write_stacktrace(fn, st)
    open(fn*"/error", "w") do f
        for e=st
          write(f,string(e)*"\n")
     end end
end
#=
    runSim(p, numSaves)

given a dictionary of inital paremeters, and a number of times to save, performs
the entire simulation with those parameters.
=#
@everywhere function runSim(p, numSaves)

    # Derive additional constants and spatial grids
    derive_constants(p)
    # Generate folder string
    fname = @sprintf("%.0fnm_%04.0fuJ_%.2fbar_%.2fbar_%.0ffs_%.1fm_%.0ffs^2_%.0fum",
                    p["λ"]*1E9, p["Energy"]*1E6, p["Pin"], p["Pout"], p["Tfwhm"]*1E15,
                    p["zmax"], p["Chirp"], p["fiberD"]*1E6)
    # Initialize electric field
    E, zinit = initialize(fname, p, "resume" in ARGS, "keep" in ARGS)
    # Save inital data
    #saveData(fname, E, 0, zinit)
    # Simulate the whole thing
    try
        simulate(E, p, zinit, fname, numSaves)
    catch ex
        if ex isa PropagationError
            warn("Error for $fname : $(ex.msg)")
        else
            warn("Unexpected error for $fname. Stacktrace written to 'error'" )
            write_stacktrace(fname, catch_stacktrace())
            #rethrow(ex)
        end
    end
end

# Load paramters which are possible list
lists = YAML.load(open(ARGS[1]))

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
param_tuples = vec(collect(product(values(lists)...)))

# Put all results into a single folder named data.
if "data" ∉ readdir()
    mkdir("data")
end
@everywhere cd("data") # Must be run on every thread

if length(ARGS) > 1
    numSaves = parse(Int, ARGS[2])
else
    println("No number of saves passed, assuming 2")
    numSaves = 2
end

print("Starting Simulations\n")
prog = Progress(length(param_tuples), 5, "Running Parallel Simulations...")
update!(prog, 0)
n = length(param_tuples)
channel = RemoteChannel(()->Channel{Bool}(n), 1)

@sync begin
    # this task prints the progress bar
    @async while take!(channel)
        next!(prog)
    end

    # this task does the computation
    @async begin
        @distributed for jawn in param_tuples
            put!(channel, true)
            runSim(Dict{Any, Any}(zip(keys(lists), jawn)), numSaves)
        end
        put!(channel, false) # this tells the printing task to finish
    end
end
