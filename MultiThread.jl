push!(LOAD_PATH, pwd())
@everywhere using SelfPhase
using Iterators
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

# Modifications required to pmap to allow for a global
# progress meter.
globalProgressMeters = Dict()
globalProgressValues = Dict()
globalPrintLock = Dict()

function Base.pmap(f::Function, p::Progress, values...; kwargs...)
    global globalProgressMeters
    global globalProgressValues
    global globalPrintLock

    id = randstring(50)
    globalProgressMeters[id] = p
    globalProgressValues[id] = 0
    globalPrintLock[id] = ReentrantLock()

    passcallback = false
    kwa = Dict(kwargs)
    if haskey(kwa,:passcallback)
      passcallback = true
      delete!(kwa,:passcallback)
    end

    out = pmap(values...; kwa...) do x...
        if passcallback
          v = f(n -> remotecall(updateProgressMeter, 1, id, n), x...)
        else
          v = f(x...)
          wait(remotecall(updateProgressMeter, 1, id, 1))
        end
        v
    end

    delete!(globalProgressMeters, id)
    out
end

#"This is remote-called by all the workers to update the progress."
@everywhere function updateProgressMeter(id,n)
    global globalProgressMeters
    global globalProgressValues
    global globalPrintLock

    lock(globalPrintLock[id])
    globalProgressValues[id] += n
    update!(globalProgressMeters[id] , globalProgressValues[id])
    unlock(globalPrintLock[id])
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
    fname = @sprintf("%.0fnm_%.0fuJ_%.2fbar_%.2fbar_%.0ffs_%.1fm_%.0ffs^2_%.0fum",
                     p["λ"]*1E9, p["Energy"]*1E6, p["Pin"], p["Pout"], p["Tfwhm"]*1E15,
                     p["zmax"], p["Chirp"], p["fiberD"]*1E6)
    # Initialize electric field
    E, zinit = initialize(fname, p, "resume" in ARGS, "keep" in ARGS)
    # Save inital data
    #saveData(fname, E, 0, zinit)
    # Simulate the whole thing
    simulate(E, p, zinit, fname, numSaves)
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
param_tuples = collect(product(values(lists)...))

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
# Parallel map over all paramters to simulate!
pmap(jawn -> runSim(Dict{Any,Any}(zip(keys(lists), jawn)), numSaves),
     prog,
     param_tuples)
ProgressMeter.finish!(prog)
