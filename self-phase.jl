using Revise
using Dierckx
using SpecialFunctions
using Base.Filesystem
using Plots

#=
 ██████  ██████  ███    ██ ███████ ████████  █████  ███    ██ ████████ ███████
██      ██    ██ ████   ██ ██         ██    ██   ██ ████   ██    ██    ██
██      ██    ██ ██ ██  ██ ███████    ██    ███████ ██ ██  ██    ██    ███████
██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██  ██ ██    ██         ██
 ██████  ██████  ██   ████ ███████    ██    ██   ██ ██   ████    ██    ███████
=#
#########################
# Simulation Parameters #
#########################

p = Dict(
    # Experimental Parameters
    "Energy" => 0.8E-3,   # Pulse Energy        J
    "Tfwhm"  => 130E-15,  # Pulse Width         s
    #"Chirp"  => 0,        # Pulse Chirp         s^2
    "λ"      => 800E-9,   # Pulse Wavelength    m
    # Numeric Parameters
    "dz"    => 2E-3,     # z-step
    "zmax"   => 0.10,     # Total length to sim m
    "Nt"     => 2*8192,   # Number of time steps
    "tmax"   => 1200E-15  # Maximum time
    )


###################
# Fixed Constants #
###################
const Fiber_D  = 270E-6     # Fiber Diameter      m
#const Fiber_L  = 100E-3     # Fiber Length        m
const Pressure = 1.4E3      # Capillary Pressure  Pa

const TOD      = 0          # 3rd Order Disp.     s^3
const losses   = 0.13       # Absorption Coef.    1/m

const Ui_Ar   = 15.75*1.6E-19 # Ionization Energy of Ar J
const α       = 7E-13
const Zeff_Ar = 1

λ_min = 400
λ_max = 1500

# General Physics Constants
const c  = 299792458        # Speed of light     m/s
const ħ = 1.0545718E-34     # Reduced Planck     Js
const me = 9.10938356E-31   # Electron Mass      kg
const ee = 1.6021766208E-19 # Elementary Charge  C
const ϵ0 = 8.854187817E-12  # Vaccum Permitivity F/m
const Kb = 1.38064852E-23   # Boltzmann Constant J/K
const T  = 300              # ~Room Temperature  K

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

#####################
# Derived Constants #
#####################
function derive_constants(p)
    f      = c / p["λ"]                 # Pulse Frequency Hz
    ω      = 2*pi*f                       # Pulse Angular Frequency
    σ_t    = p["Tfwhm"]/sqrt(2*log(2))    # 1-sigma width of pulse
    Power  = sqrt(2/pi) * p["Energy"]/σ_t # Max power delivered by pulse

    dt     = p["tmax"]/(p["Nt"]-1)        # Time step
    points = (-p["Nt"]/2:1:p["Nt"]/2-1)
    t_vec  = (-p["Nt"]/2:1:p["Nt"]/2)*dt  # Time grid iterator

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

        "n_tot_0" => n_tot_0
    )
    merge!(p,dp)
end
derive_constants(p)

λ_test = minimum(abs.(p["λ_tot"]-λ_min))
λ_min = filter(x -> abs(x - λ_min) == λ_test, p["λ_tot"])[1]
λ_test = minimum(abs.(p["λ_tot"]-λ_max))
λ_max = filter(x -> abs(x - λ_max) == λ_test, p["λ_tot"])[1]

const ρ_crit = p["ωω"].^2*me*ϵ0/ee^2
const k_Ar   = ceil.(Ui_Ar ./ (ħ*p["ωω"]))

#=
██    ██  █████  ██████  ██  █████  ██████  ██      ███████ ███████
██    ██ ██   ██ ██   ██ ██ ██   ██ ██   ██ ██      ██      ██
██    ██ ███████ ██████  ██ ███████ ██████  ██      █████   ███████
 ██  ██  ██   ██ ██   ██ ██ ██   ██ ██   ██ ██      ██           ██
  ████   ██   ██ ██   ██ ██ ██   ██ ██████  ███████ ███████ ███████
=#
##################
# Initialisation #
##################
# Electric Field
E             = exp.(-p["t_vec"].^2/p["σ_t"]^2)
E0            = sqrt(2*p["Power"]/(pi*Fiber_D^2))
E             = E0 .* E
spectrum_init = abs.(fftshift(fft(fftshift(E)))).^2
I_init        = sum(spectrum_init)

# Propagation variables
zpoints       = 1000
dist          = zeros(zpoints)
ρ_max         = zeros(zpoints)
I_max         = zeros(zpoints)
ΔT_pulse      = zeros(zpoints)
It_dist       = zeros(p["Nt"],zpoints)
spectrum_dist = zeros(p["Nt"],zpoints)

#Chirp_function = p["Chirp"] * p[entrance"ωω"] .^ 2 + TOD .* p["ωω"] .^ 3
Chirp_function = 0
E_TF           = fftshift(fft(fftshift(E))).*exp.(im * Chirp_function)
E              = ifftshift(ifft(ifftshift(E_TF)))

#Clear vars
Chirp_function = nothing; E_TF = nothing; gc()

z   = 0
JJJ = 0
ZZZ = 0

idxp = push!(collect(2:p["Nt"]), 1)
idxn = append!([p["Nt"]], collect(1:p["Nt"]-1))

pressure = []

#=
███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
█████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████
=#
function calc_pressure(p_in, p_out, z, zmax) #tested
    #=
    Calculate pressure along fiber
    param p_in:  Pressure at input of fiber
    param p_out: Pressure at exit of fiber
    param z:     Position to calculate pressure at
    param zmax:  Length of fiber

    return: Pressure at z
    =#
    return sqrt(p_in^2 + (z/zmax) * (p_out^2-p_in^2))
end

function calc_duration(E, t1) #tested
    center_pulse = sum(t1.*(abs.(E)).^2)/sum(abs.(E).^2)
    return 2 * sqrt(2*log(2)).*((sum((t1-center_pulse).^2.*abs.(E).^2)
                               / sum(abs.(E).^2)).^0.5)*1E15
end

function prop_lin(p, E, deriv_t_2, losses) #tested
    # Shift to frequency domain and compute linear propagation
    E_TF = fftshift(fft(fftshift(E))) .* exp.(1im*(deriv_t_2)* p["dz"])
    # Shift back to time domain, compute losses and return
    return ifftshift(ifft(ifftshift(E_TF))).*exp.(-losses/2 * p["dz"])
end

function prop_non_lin(p,E, rrr, ρ, losses, kerr_response) #tested
    return E.*exp(rrr*ρ*p["dz"] - losses + kerr_response)
end

function calc_compression(p, width, Cs) #tested
    λ_mu = p["λ_tot"]*1E-3
    λ_test = minimum(abs.(p["λ_tot"] - 600))
    λ_first = p["λ_tot"][findfirst(x -> abs(x - 600) == λ_test && x > 0,
                                   p["λ_tot"])]

    λ_test = minimum(abs.(p["λ_tot"] - 6000))
    λ_final = p["λ_tot"][findfirst(x -> abs(x - 6000) == λ_test && x > 0,
                                   p["λ_tot"])]
    n_fs = sqrt.(1 + 0im + Cs[1]*(λ_mu).^2./((λ_mu).^2-Cs[2]^2) +
                           Cs[3]*(λ_mu).^2./((λ_mu).^2-Cs[4]^2) +
                           Cs[5]*(λ_mu).^2./((λ_mu).^2-Cs[6]^2))

    n_begin = n_fs[findfirst(x->x==λ_first, p["λ_tot"])]
    n_fs[p["λ_tot"].<600] = n_begin

    n_final = n_fs[findfirst(x->x==λ_final, p["λ_tot"])]
    n_fs[p["λ_tot"].>6000] = n_final

    n_fs_spline = Spline1D(p["ωω_tot"], n_fs)
    n_fs_interpd = evaluate(n_fs_spline, p["ωω_tot"])
    dn_fs_interpd = derivative(n_fs_spline, p["ωω_tot"])

    k_FS=n_fs_interpd.*p["ωω_tot"]/c
    k_FS[find(x->x < 245, p["λ_tot"])] = maximum(k_FS)

    k_prime_FS=1 / c * (dn_fs_interpd.*p["ωω_tot"]+n_fs_interpd)

    k0_FS=k_FS[findfirst(x->x==p["ω"], p["ωω_tot"])]
    k1_FS=k_prime_FS[findfirst(x->x==p["ω"], p["ωω_tot"])]

    return [(k_FS-k0_FS-k1_FS.*p["ωω"]).*(w.*1e-3) for w in width]
end

function smooth(values, radius)
    smoothed_values = zeros(values)
    for i in eachindex(values)
        temp_radius = minimum([radius, i - 1, length(values) - i])
        smoothed_values[i] = mean(values[(i - temp_radius):(i + temp_radius)])
    end
    return smoothed_values
end

#Place this definition inside the next one?
function NL_response(p, E, idxp, idxn, γs) #Tested, fixed error in matlab
    #=
    Functionized version of the repeated code in steppening.m
    =#
    NL = sum([γs[i] * abs.(E).^(2*i) for i in eachindex(γs)]) * p["dz"]
    Temp = NL .* E
    return (im / p["ω"]) * ((Temp[idxp] - Temp[idxn])/(2 * p["dt"]))
end

#Make this operate in place?
function steepening(p, E, idxp, idxn, γs) #Tested
    Etemp = E + (0.5 * NL_response(p, E, idxp, idxn, γs))
    return E + NL_response(p, Etemp, idxp, idxn, γs)
end

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
    return [k,k1,k2,k3,k4]
end

function calc_ns(pressure, n, n_tot, λ_tot) #tested
    n2_800  = 1e-23   * pressure
    n4_800  =-3.7e-42 * pressure
    n6_800  = 4e-58   * pressure
    n8_800  =-1.7e-75 * pressure
    n10_800 = 8.8e-94 * pressure

    n_800 = n_tot[findfirst(x->x==minimum(abs.(λ_tot-800)), abs.(λ_tot - 800))]
    n2  = n2_800 * ((n^2 - 1)/(n_800^2-1))^4
    n4  = n4_800 * ((n^2 - 1)/(n_800^2-1))^6
    n8  = n8_800 * ((n^2 - 1)/(n_800^2-1))^10
    n10 = n10_800 * ((n^2 - 1)/(n_800^2-1))^12

    return [n,n2,n4,n8,n10]
end

function plasma_potential(E,ω,Zeff,Ui) #Tested
    """
    Derived From PPT Theory
    http://jetp.ac.ru/cgi-bin/dn/e_023_05_0924.pdf
    """
    Uh = 13.5984*ee                     # Hydrogen Ionization Potential
    ω_au = 4.1E16                       # Ionization Potential (1/s Natural)
    E = abs.(E) * sqrt(2/(ϵ0*c))              # Renomrmalize E-field for proper units
    γ = ω .* sqrt(2*me*Ui)./(ee*E)
    Eh = ee^5*me^2/(ħ^4 * (4*π*ϵ0)^3)
    E0 = Eh * (Ui/Uh) ^ (3/2)
    A = zeros(γ)
    β = 2*γ ./ sqrt.(1+γ.^2)
    α = 2.*asinh.(γ) - β

    g=3./(2*γ).*((1+1./(2*γ.^2)).*asinh.(γ)-1./β)
    ν0 = Ui / (ħ*ω)
    ν  = ν0 * (1 + 1./(2*γ.^2))
    kmin = minimum(floor.(ν) + 1)
    l=0
    m=0 #??Why
    n_star = Zeff * sqrt(Uh/Ui)

    C_nl2 = 2^(2*n_star) / (n_star*gamma(2*n_star))
    f = (2*l+1) * factorial(l+abs(m)) / (2^abs(m)) *
        factorial(abs(m)) * factorial(l-abs(m))

    A = sum([4/sqrt(3*π) * γ.^2 ./ (1+γ.^2) .* exp.(-α .* (z-ν)) .*
             dawson.(sqrt.(abs.(β.*(z-ν)))) for z in kmin:kmin+2])

    potential = ω_au * C_nl2 * f * sqrt(6/π) * (Ui / (2 * Uh)) *
                A .* (2 * E0./(E .* sqrt.(1+γ.^2))).^(2 * n_star - abs(m) - 3/2) .* exp.(-2 * E0 * g ./ (3*E))
    #potential[isnan.(potential)]=0
    return potential
end

function plasma(p, α, ρ_at, Potentiel_Ar, E, C2_Ar) #tested
    ρ_Ar = zeros(size(E))

    for i in 1:(p["Nt"]-1)
        ρ_Ar[i+1] = ρ_Ar[i] +
                    p["dt"] * (-α*ρ_Ar[i]^2+Potentiel_Ar[i]*(ρ_at - ρ_Ar[i]) +
                    (C2_Ar * abs(E[i])^2)*ρ_Ar[i])
    end
    return ρ_Ar
end

#=
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
=#
# Saving/Loading Files
#=
fname = "$(p["λ"])nm_$(p["Energy"])J_$(p["Tfwhm"])s"
if fname ∈ readdir()
    print("Found data for these parameters, load and continue? (y)/n: ")
    input = readline()
    if input != "n"
        #Load Data
    else
        print("Overwrite old data? y/(n): ")
        input = readline()
        if input != "y"
            i = 2
            while fname*"_($i)" ∈ readdir()
                i += 1
            end
            fname = fname*"_($i)"
        else
            rm(fname)
        end
        mkdir(fname)
    end
else
    mkdir(fname)
end
fparams = open("$fname/params", "w")
fE = open("$fname/E", "w")
=#

#Temporary Flag
run = false
if run
for iter in round(Int, zinit/p["dz"]):1:round(Int, p["zmax"]/p["dz"])
    z = iter * p["dz"]
    # Calculate Pressure
    pressure_z = calc_pressure(0.008, Pressure, z, p["zmax"])
    push!(pressure, pressure_z)

    # Update of n, with cutoffs
    n_tot = sqrt(1+pressure_z * (p["n_tot_0"].^2 - 1))
    n_max = n_tot[findfirst(λ->λ==λ_max, p["λ_tot"])]
    n_tot[p["λ_tot"] .> λ_max] = n_max
    n_min = n_tot[findfirst(λ->λ==λ_min, p["λ_tot"])]
    n_tot[p["λ_tot"] .< λ_min] = n_min

    ks = calc_ks()
    n  = n_tot[findfirst(x->x==p["ωω"],p["ωω_tot"])]
    vg = 1/ks[2]

    # Argon Parameters
    ns = calc_ns()
    β2 = pressure_z * ks[3]
    β3 = pressure_z * ks[4]
    β4 = pressure_z * ks[5]
    τ = 3.5E-13 / pressure_z
    ρ_at = pressure_z * 1E5 / (Kb * T)

    # Plasma Parameters
    σ_k = 2.81E-96 * pressure
    σ   = (ks[1]*ee^2) / (p["ωω"] * me * ϵ0) * τ/(1+(p["ωω"] * τ).^2)
    β_k = 10^(-4 * k_Ar) * k_Ar * ħ * p["ωω"] * ρ_at * 0.21 * σ_k
    rrr = -im * ks[1]/(2 * n[1]^2 * ρ_crit) - 0.5 * σ
    coeff2 = σ/Ui_Ar

    # Dispersion and Laplacian Operators
    dv_t_2_op = p["k_tot"] - ks[1] - ks[2] * p["ωω"]

    # Kerr factors
    γs = im * n[2:end] * ks[1]/n[1]

    # Propagation
    E = prop_lin(p, E, dv_t_2_op, losses)    #Linear
    E = steepening(p, E, idxp, idxn, γs)     #Steepening

    # Plasma
    U_ion = plasma_potential(abs(E).^2, c, p["ωω"])
    ρ = plasma(p, α, ρ_at, U_ion, E, coeff2, )
    plasma_loss = U_ion / (2 * abs.(E).^2) * Ui_Ar * (ρ_at - ρ) * p["dz"]
    plasma_loss[isnan.(plasma_loss)] = 0

    # Kerr and Plasma Propagation (NonLinear)
    kerr_response = -(γs[1]*(abs.(E)).^2 + γs[2]*(abs.(E)).^4 + γs[3] * abs.(E).^6 + γs[4]*abs.(E).^8 + γ[5].*abs.(E).^10) * p["dz"]
    E = prop_non_lin(E, rrr, ρ, p["dz"], plasma_loss, kerr_response)

    # Update Distances
    ZZZ += p["dz"]
    z   += p["dz"]
    distance = 100*ZZZ
    dist[1, JJJ+1] = z

    #Saving E-field
    writecsv(fE, E)
    I_max[iter+1]         = maximum(abs.(E).^2)
    It_dist[iter+1]       = abs.(E).^2
    spectrum_dist[iter+1] = abs.(fftshift(fft(fftshift(E)))).^2
    ΔT_pulse[iter+1]      = calc_duration(E,t)
    ρ_max[iter+1]         = maximum(ρ*1E-6)
    #TODO:Save intermediate data
    if (iter%save_every == 0)
    end
end
end
