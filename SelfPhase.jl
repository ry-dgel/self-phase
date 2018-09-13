module SelfPhase
export derive_constants, simulate, initialize, saveData

using Dierckx
using SpecialFunctions.dawson
using Base.Filesystem
using ProgressMeter
using YAML

#=
 ██████  ██████  ███    ██ ███████ ████████  █████  ███    ██ ████████ ███████
██      ██    ██ ████   ██ ██         ██    ██   ██ ████   ██    ██    ██
██      ██    ██ ██ ██  ██ ███████    ██    ███████ ██ ██  ██    ██    ███████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██  ██ ██    ██         ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██ ██   ████    ██    ███████
=#
# Does this speed things up??
FFTW.set_num_threads(nprocs())
BLAS.set_num_threads(nprocs())

###################
# Fixed Constants #
###################
# Setup constants
const losses  = 0.13          # Absorption Coef.        1/m
const Ui_Ar   = 15.75*1.6E-19 # Ionization Energy of Ar J
const α       = 7E-13         # ???
const Zeff_Ar = 1.0           # Effective Charge of ionized Argon

# General Physics Constants
const c  = 299792458.0      # Speed of light     m/s
const ħ = 1.0545718E-34     # Reduced Planck     Js
const me = 9.10938356E-31   # Electron Mass      kg
const ee = 1.6021766208E-19 # Elementary Charge  C
const ϵ0 = 8.854187817E-12  # Vaccum Permitivity F/m
const Kb = 1.38064852E-23   # Boltzmann Constant J/K
const T  = 300.0            # ~Room Temperature  K

# Indices of refraction/dispersion for Argon 1/nm
const C1 = 0.012055
const C2 = 0.2075
const C3 = 91.012
const C4 = 0.0415
const C5 = 87.892
const C6 = 4.3330
const C7 = 214.02
C = [C1,C2,C3,C4,C5,C6,C7]

#Indices of Refraction for Fused Silica 1/nm
const Cfs1 = 0.6961663
const Cfs2 = 0.0684043
const Cfs3 = 0.4079426
const Cfs4 = 0.1162414
const Cfs5 = 0.8974794
const Cfs6 = 9.896161
Cfs = [Cfs1,Cfs2,Cfs3,Cfs4,Cfs5,Cfs6]

# Plasma Ionization Constants
const a0 = -185.8
const a1 =  11.16
const a2 =  4.763e-3
const a3 = -9.946e-4
const a4 = -9.722e-5
const a5 =  1.182e-5
const a6 =  4.63e-7
const a7 = -4.227e-8
as = [a0,a1,a2,a3,a4,a5,a6,a7]

#=
███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
█████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████
=#

#################
# File Handling #
#################
"""
    saveData(fname, E, ΔT_pulse, z)

Saves all the data to respective files within the folder given by fname
"""
function saveData(fname, E, ΔT_pulse, z)
    open(fname * "/E", "a") do f
        #write(f, @sprintf("\# z = %f\n", z))
        for pair in zip(real(E),imag(E))
            write(f, @sprintf("%.15e,%.15e\n",pair...))
        end
        write(f, "\n")
    end
    open(fname * "/z", "a") do f
        write(f, "$z\n")
    end
end

"""
    loadParams(fname)

Loads the parameters from fname, a yaml file, into a dictionary.
"""
function loadParams(fname)
    p = YAML.load(open(fname))

    # We use λ in the code, so replace the key by removing the value and reinserting it
    if haskey(p, "lambda")
        merge!(p, Dict("λ" => pop!(p, "lambda")))
    end

    # Derive additional constants
    derive_constants(p)
    return p
end

"""
    saveParams(fname, p)

Saves the important parameters from a dictionary into a file given by fname.
"""
function saveParams(fname, p)
    open("$fname/params", "w") do f
        write(f, @sprintf("Energy:    [%e]\n", p["Energy"]))
        write(f, @sprintf("Tfwhm:     [%e]\n", p["Tfwhm"]))
        write(f, @sprintf("lambda:    [%e]\n", p["λ"]))
        write(f, @sprintf("dz:        [%f]\n", p["dz"]))
        write(f, @sprintf("zmax:      [%f]\n", p["zmax"]))
        write(f, @sprintf("Nt:        [%d]\n", p["Nt"]))
        write(f, @sprintf("tmax:      [%e]\n", p["tmax"]))
        write(f, @sprintf("Pin:       [%f]\n", p["Pin"]))
        write(f, @sprintf("Pout:      [%f]\n", p["Pout"]))
        write(f, @sprintf("fiberD:    [%e]\n", p["fiberD"]))
        if haskey(p, "Chirp")
            write(f, @sprintf("Chirp:     [%d]\n",     p["Chirp"]))
        end
        if haskey(p, "TOD")
            write(f, @sprintf("TOD:       [%d]\n",     p["TOD"]))
        end
    end
    writedlm("$fname/t_vec", p["t_vec"])
end

"""
    initialize(fname, p, resume, keep)

Initializes a simulation given a folder name that may or may not already contain
data. Folder name should be relevant to the paramters in p. If resume is true,
will attempt to continue a simulaiton. If keep is true, will not overwrite old
data already present in folder given by fname.
"""
function initialize(fname, p, resume, keep)
    zinit = 0
    # Check for existing folder
    if fname ∈ readdir()
        # Resume simulations if desired
        if resume
            zinit = float(read("$fname/z"))
            E[:] = readcsv("$fname/E")[end][:]
        else
            # Don't overwrite data, append number to folder name.
            if keep
                i = 2
                while fname*"_($i)" ∈ readdir()
                    i += 1
                end
                fname = fname*"_($i)"
            else
                rm(fname, recursive=true)
            end
            # Otherwise make new folder and initialize data
            mkdir(fname)
            E = initField(p)
            saveParams(fname, p)
        end
    else
        # Otherwise make new folder and initialize data
        mkdir(fname)
        E = initField(p)
        saveParams(fname, p)
    end
    return E, zinit
end

#######################
# Handling Parameters #
#######################
"""
    derive_constants(p)

Adds several concrete values and vectors to the parameter dictionary.
"""
function derive_constants(p)
    if haskey(p, "lambda")
        merge!(p, Dict("λ" => pop!(p, "lambda")))
    end

    f      = c / p["λ"]                 # Pulse Frequency Hz
    ω      = 2*pi*f                       # Pulse Angular Frequency
    σ_t    = p["Tfwhm"]/sqrt(2*log(2))    # 1-sigma width of pulse
    Power  = sqrt(2/pi) * p["Energy"]/σ_t # Max power delivered by pulse

    dt     = p["tmax"]/(p["Nt"]-1)        # Time step
    points = (-p["Nt"]/2:1:p["Nt"]/2-1)
    t_vec  = (-p["Nt"]/2:1:p["Nt"]/2-1)*dt  # Time grid iterator

    ff     = points./p["tmax"]            # Frequency grid
    ωω     = (2*pi)*ff                    # Angular frequency grid

    # Peak centered Wavelength grid, in nm
    λ_tot  = 1E9 * c ./ (f + ff)
    λ_tot_micron = 1E6 * c ./ (f + ff)
    # Peak centered angular frequency grid
    ωω_tot = ω+ωω

    # Nonlinear index of refraction
    n_tot_0 = 1 + C1 * (C2 * (λ_tot_micron.^2) ./ (C3 * (λ_tot_micron.^2) -1) +
                        C4 * (λ_tot_micron.^2) ./ (C5 * (λ_tot_micron.^2) -1) +
                        C6 * (λ_tot_micron.^2) ./ (C7 * (λ_tot_micron.^2) -1))

    λ_min = 400
    λ_max = 1500
    λ_min = λ_tot[indmin(abs.(λ_tot-λ_min))]
    λ_max = λ_tot[indmin(abs.(λ_tot-λ_max))]

    ρ_crit = ω^2*me*ϵ0/ee^2
    k_Ar   = ceil(Ui_Ar / (ħ*ω))

    dp = Dict(
        "f"       => f,
        "ω"       => ω,
        "σ_t"     => σ_t,
        "Power"   => Power,

        "dt"      => dt,
        "points"  => points,
        "t_vec"   => t_vec,

        "ff"      => ff,
        "ωω"      => ωω,

        "λ_tot"   => λ_tot,
        "ωω_tot"  => ωω_tot,
        "λ_max"   => λ_max,
        "λ_min"   => λ_min,

        "n_tot_0" => n_tot_0,

        "ρ_crit"  => ρ_crit,
        "k_Ar"    => k_Ar
    )
    merge!(p,dp)
end

######################
# Simulation Helpers #
######################
"""
    initField(p)

Returns a vector representing an electric field defined in terms of the
peak power and spread of the intial gaussian pulse. Can also handle
Chirp and TOD. All paramters read from dictionary p.
"""
function initField(p)
    E             = exp.(-p["t_vec"].^2/p["σ_t"]^2)
    E0            = sqrt(2*p["Power"]/(pi*(p["fiberD"]/2)^2))
    E             = E0 .* E

    # Use get so as to have a default value of zero instead
    # of an exception in case of Chirp or TOD not being defined.
    Chirp_function = get(p,"Chirp", 0) * p["ωω"] .^ 2 +
                     get(p, "TOD", 0) .* p["ωω"] .^ 3
    if !iszero(Chirp_function)
	       E_TF           = fftshift(fft(fftshift(E))).*exp.(im * Chirp_function)
           E              = ifftshift(ifft(ifftshift(E_TF)))
    end

    return E
end

"""
    calc_pressure(p_in, p_out, z, zmax)

Computes the pressure at a point z along a linear pressure gradient.
"""
function calc_pressure(p_in, p_out, z, zmax) #tested
    #Calculate pressure along fiber
    return sqrt(p_in^2 + (z/zmax) * (p_out^2-p_in^2))
end

"""
    calc_duration(E, t1)

Computes the FHWM of a gaussian peak, doesn't appear to work properly.
"""
function calc_duration(E, t1) #tested
    center_pulse = sum(t1.*(abs2.(E)))/sum(abs2.(E))
    return 2 * sqrt(2*log(2)).*((sum((t1-center_pulse).^2.*abs2.(E))
                               / sum(abs.(E).^2)).^0.5)*1E15
end

"""
    prop_lin(p, E, deriv_t_2, losses, ft, ift)

Computes the linear evoltion of a field E. deriv_t_2 is computed in sim step.
ft and ift are pre-planned FFT matrices for forward and reverse FT. losses is
the absorption coeff of the fiber.
"""
function prop_lin(p, E, deriv_t_2, losses, ft, ift) #tested
    # Shift to frequency domain and compute linear propagation
    E_TF = fftshift(ft * (fftshift(E))) .* exp.(1im*(deriv_t_2)* p["dz"])
    # Shift back to time domain, compute losses and return
    return ifftshift(ift * (ifftshift(E_TF))) .* exp.(-losses/2 * p["dz"])
end

"""
    prop_non_lin(p, E, rrr, ρ, losses, kerr_response)

Computes the nonlinear evolution of a field E. rrr, ρ, losses and kerr_response
all computed in main.
"""
function prop_non_lin(p, E, rrr, ρ, losses, kerr_response) #tested
    return E.*exp.(rrr.*ρ*p["dz"] - losses + kerr_response)
end

"""
    smooth(values, radius)

Smooths a vector by averaging all values within radius of each element.
"""
function smooth(values, radius)
    smoothed_values = zeros(values)
    for i in eachindex(values)
        # Handles edges by smoothing over a smaller radius
        temp_radius = minimum([radius, i - 1, length(values) - i])
        smoothed_values[i] = mean(values[(i - temp_radius):(i + temp_radius)])
    end
    return smoothed_values
end

"""
    NL_response(p, E, γs)

Computes the electric field brought on by steepening.
"""
function NL_response(p, E, γs)
    #=
    Functionized version of the repeated code in steppening.m
    =#
    NL = sum([γs[i] * abs.(E).^(2*i) for i in eachindex(γs)]) * p["dz"]
    Temp = NL .* E
    return (im / p["ω"]) *
           ((circshift(Temp,1) - circshift(Temp,-1))/(2 * p["dt"]))
end

"""
    steepening(p, E, γs)

Computes the evoltion of field E due to steepening.
"""
function steepening(p, E, γs) #Tested
    Etemp = E + (0.5 * NL_response(p, E, γs))
    return E + NL_response(p, Etemp, γs)
end

"""
    calc_ks(p, n_tot)

Computes the kerr coefficients for a given frequency.
"""
function calc_ks(p, n_tot) #tested
    # Interpolation of n
    n_interp = Spline1D(p["ωω_tot"], n_tot,k=4)
    n_tot = evaluate(n_interp, p["ωω_tot"])
    dn    = derivative(n_interp, p["ωω_tot"], nu=1)
    d2n   = derivative(n_interp, p["ωω_tot"], nu=2)
    d3n   = derivative(n_interp, p["ωω_tot"], nu=3)
    d4n   = derivative(n_interp, p["ωω_tot"], nu=4)

    # Frequency Disperion
    k_tot = n_tot .* p["ωω_tot"]/c
    k_tot[p["λ_tot"] .< 245] = maximum(k_tot)

    k_first  = 1/c * (dn .* p["ωω_tot"] + n_tot)
    k_second = 1/c * (d2n .* p["ωω_tot"] + 2 * dn)
    k_third  = 1/c * (d3n .* p["ωω_tot"] + 3 * d2n)
    k_fourth = smooth(1/c * (d4n .* p["ωω_tot"] + 4 * d3n), 100)

    k  = k_tot[findfirst(x->x==p["ω"],p["ωω_tot"])]
    k1 = k_first[findfirst(x->x==p["ω"],p["ωω_tot"])]
    k2 = k_second[findfirst(x->x==p["ω"],p["ωω_tot"])]
    k3 = k_third[findfirst(x->x==p["ω"],p["ωω_tot"])]
    k4 = k_fourth[findfirst(x->x==p["ω"],p["ωω_tot"])]
    return k_tot, [k,k1,k2,k3,k4]
end

"""
    calc_ns(pressure, n, n_tot, λ_tot)

Computes the index of refraction for a given pressure along the entire frequency
range.
"""
function calc_ns(pressure, n, n_tot, λ_tot) #tested
    n2_800  = 1e-23   * pressure
    n4_800  =-3.7e-42 * pressure
    n6_800  = 4e-58   * pressure
    n8_800  =-1.7e-75 * pressure
    n10_800 = 8.8e-94 * pressure

    n_800 = n_tot[indmin(abs.(λ_tot - 800))]
    n2  = n2_800 * ((n^2 - 1)/(n_800^2-1))^4
    n4  = n4_800 * ((n^2 - 1)/(n_800^2-1))^6
    n6  = n6_800 * ((n^2 - 1)/(n_800^2-1))^8
    n8  = n8_800 * ((n^2 - 1)/(n_800^2-1))^10
    n10 = n10_800 * ((n^2 - 1)/(n_800^2-1))^12

    return [n,n2,n4,n6,n8,n10]
end

"""
    plasma_potential(E, p, Ui)

Computes the chance of ionization based off of PPT theory.
"""
function plasma_potential(E,p,Ui) #Tested
    """
    Derived From PPT Theory
    http://jetp.ac.ru/cgi-bin/dn/e_023_05_0924.pdf
    """
    Zeff = 1
    Uh = 13.5984*ee                     # Hydrogen Ionization Potential
    ω_au = 4.1E16                       # Ionization Potential (1/s Natural)
    E = abs.(E) * sqrt(2/(ϵ0*c))        # Renomrmalize E-field for proper units
    γ = p["ω"] .* sqrt(2*me*Ui)./(ee*E)
    Eh = ee^5*me^2/(ħ^4 * (4*π*ϵ0)^3)
    E0 = Eh * (Ui/Uh) ^ (3/2)
    A = zeros(γ)
    β = 2*γ ./ sqrt.(1+γ.^2)
    α = 2.*asinh.(γ) - β

    g=3./(2*γ).*((1+1./(2*γ.^2)).*asinh.(γ)-1./β)
    ν0 = Ui / (ħ*p["ω"])
    ν  = ν0 * (1 + 1./(2*γ.^2))
    kmin = minimum(floor.(ν) + 1)
    # These Values are gas dependent, see paper.
    l=0
    m=0
    n_star = Zeff * sqrt(Uh/Ui)

    C_nl2 = 2^(2*n_star) / (n_star*gamma(2*n_star))
    f = (2*l+1) * factorial(l+abs(m)) / (2^abs(m)) *
        factorial(abs(m)) * factorial(l-abs(m))

    A = sum([4/sqrt(3*π) * γ.^2 ./ (1+γ.^2) .* exp.(-α .* (z-ν)) .*
             dawson.(sqrt.(complex(abs.(β.*(z - ν))))) for z in kmin:kmin+2])

    potential = ω_au * C_nl2 * f * sqrt(6/π) *
                (Ui / (2 * Uh)) * A .*
                (2 * E0./(E .* sqrt.(1+γ.^2))).^(2 * n_star - abs(m) - 3/2) .*
                exp.(-2 * E0 * g ./ (3*E))

    # present in original code, seemingly not needed
    # potential[isnan.(potential)]=0

    return potential
end

"""
    plasma(p, α, ρ_at, Potentiel_Ar, E, coeff2)

Computes the plasma density along the pulse.
"""
function plasma(p, α, ρ_at, Potentiel_Ar, E, coeff2) #tested
    ρ_Ar = zeros(size(E))
    for i in 1:(p["Nt"]-1)
        ρ_Ar[i+1] = ρ_Ar[i] +
                    p["dt"] * (-α*ρ_Ar[i]^2 + Potentiel_Ar[i]*(ρ_at - ρ_Ar[i]) +
                               (coeff2 * abs(E[i])^2)*ρ_Ar[i])
    end
    return ρ_Ar
end

"""
    simulate(E, p, zinit, fname, num_saves)

Performs the full simuation of a pulse propagating. Saves results to fname. Will
save num_saves times with even propagation length between saves.
"""
function simulate(E, p, zinit, fname, num_saves)

    ##################
    # Initialisation #
    ##################
    steps = round(Int, p["zmax"]/p["dz"])-round(Int, zinit/p["dz"])

    # How often to save data.
    if num_saves > 2
        save_every = ceil((steps-1)/((num_saves)-2))
    else
        # Data is explicitely saved at initial and final points.
        save_every = Inf
    end

    #Plan FFT
    ft  = plan_fft(E)
    ift = plan_ifft(E)

    # Propagation variables
    z = zinit
    ρ = 0
    while z < p["zmax"]
        if (round(z/p["dz"])%save_every == 0)
            # Async data write
            saveData(fname, E,
                     calc_duration(E,p["t_vec"]*p["dt"]), z)
        end
        open(fname * "/PlasmaDensity", "a") do f
            write(f, @sprintf("%.5e\n", ρ))
        end

        E, ρ = simStep(E, p, z, ft, ift)

        # Update Distance
        z += p["dz"]
    end

    #Save final data
    saveData(fname, E, calc_duration(E,p["t_vec"]*p["dt"]), z)
end

"""
    simStep(E, p, z, ft, ift)

Performs a single split time step of the simulation.
"""
function simStep(E, p, z, ft, ift)
    # Calculate Pressure
    pressure_z = calc_pressure(p["Pin"], p["Pout"], z, p["zmax"])

    # Update of n, with cutoffs
    n_tot = sqrt.(complex(1+pressure_z * (p["n_tot_0"].^2 - 1)))
    n_max = n_tot[findfirst(λ->λ==p["λ_max"], p["λ_tot"])]
    n_tot[p["λ_tot"] .> p["λ_max"]] = n_max
    n_min = n_tot[findfirst(λ->λ==p["λ_min"], p["λ_tot"])]
    n_tot[p["λ_tot"] .< p["λ_min"]] = n_min

    k_tot, ks = calc_ks(p, n_tot)
    n  = n_tot[findfirst(x->x==p["ω"],p["ωω_tot"])]
    vg = 1/ks[2]

    # Argon Parameters
    ns = calc_ns(pressure_z, n, n_tot, p["λ_tot"])
    #β2 = pressure_z * ks[3]
    #β3 = pressure_z * ks[4]
    #β4 = pressure_z * ks[5]
    τ = 3.5E-13 / pressure_z
    ρ_at = pressure_z * 1E5 / (Kb * T)

    # Plasma Parameters
    σ_k = 2.81E-96 * p["Pout"]
    σ   = (ks[1]*ee^2) ./ (p["ω"] * me * ϵ0) .* τ./(1+(p["ω"] * τ).^2)
    #β_k = 10.^(-4 * p["k_Ar"]) .* p["k_Ar"] * ħ .* p["ω"] * ρ_at * 0.21 * σ_k
    rrr = -im * ks[1]./(2 * n[1]^2 * p["ρ_crit"]) - 0.5 * σ
    coeff2 = σ/Ui_Ar

    # Dispersion and Laplacian Operators
    dv_t_2_op = k_tot - ks[1] - ks[2] * p["ωω"]

    # Kerr factors
    γs = im * ns[2:end] * ks[1]/n[1]

    # Propagation
    E = prop_lin(p, E, -dv_t_2_op, losses, ft, ift)   #Linear
    E = steepening(p, E, γs)                          #Steepening

    # Plasma
    U_ion = plasma_potential(E, p, Ui_Ar)
    ρ = plasma(p, α, ρ_at, U_ion, E, coeff2)
    plasma_loss = U_ion ./ (2 * abs.(E).^2) * Ui_Ar .* (ρ_at - ρ) * p["dz"]
    plasma_loss[isnan.(plasma_loss)] = 0

    # Kerr and Plasma Propagation (NonLinear)
    kerr_response = -(γs[1] * abs.(E).^2 +
                      γs[2] * abs.(E).^4 +
                      γs[3] * abs.(E).^6 +
                      γs[4] * abs.(E).^8 +
                      γs[5] * abs.(E).^10
                     ) * p["dz"]

    E = prop_non_lin(p, E, rrr, ρ, plasma_loss, kerr_response)

    return E, maximum(ρ) * 1E-6
end

#module
end
