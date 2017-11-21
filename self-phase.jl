using Revise
using Dierckx
using SpecialFunctions
#using Plots

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
# Experimental Parameters
const Energy   = 0.8E-3     # Pulse Energy       J
const Tfwhm    = 130E-15     # Pulse Width        s
const Chirp    = 0          # Pulse Chirp        s^2
const λ        = 800E-9     # Pulse Wavelength   m

const Fiber_D  = 270E-6     # Fiber Diameter     m
const Fiber_L  = 100E-3     # Fiber Length       m
const Pressure = 1.4E3      # Capillary Pressure Pa

const TOD      = 0          # 3rd Order Disp.    s^3
const losses   = 0.13       # Absorption Coef.   1/m

const Ui_Ar   = 15.75*1.6E-19 # Ionization Energy of Ar J
const α       = 7E-13
const Zeff_Ar = 1

# Numeric Parameters
const dz0=2E-3       # z-step 1
const dz1=2E-4       # z-step 2
const zmax = Fiber_L # Total z

const Nt=2*8192          # Number of time steps
const tmax=1200E-15      # Maximum time

const λ_min = 400
const λ_max = 1500
###################
# Fixed Constants #
###################
# General Physics Constants
const c  = 299792458        # Speed of light     m/s
const ħ = 1.0545718E-34     # Reduced Planck     Js
const me = 9.10938356E-31   # Electron Mass      kg
const ee = 1.6021766208E-19 # Elementary Charge  C
const ϵ0 = 8.854187817E-12  # Vaccum Permitivity F/m
const Kb = 1.38064852E-23   # Boltzmann Constant J/K
const T  = 300              # ~Room Temperature  K

# Indices of refraction/dispersion for Argon 1/nm
const C1 = 0.012055E3
const C2 = 0.2075E3
const C3 = 91.012E3
const C4 = 0.0415E3
const C5 = 87.892E3
const C6 = 4.3330E3
const C7 = 214.02E3
const C = [C1,C2,C3,C4,C5,C6,C7]

#Indices of Refraction for Fused Silica 1/nm
const Cfs1 = 0.6961663E3
const Cfs2 = 0.0684043E3
const Cfs3 = 0.4079426E3
const Cfs4 = 0.1162414E3
const Cfs5 = 0.8974794E3
const Cfs6 = 9.896161
const Cfs = [Cfs1,Cfs2,Cfs3,Cfs4,Cfs5,Cfs6]

# Plasma Ionization Constants
const a0 = -185.8
const a1 =  11.16
const a2 =  4.763e-3
const a3 = -9.946e-4
const a4 = -9.722e-5
const a5 =  1.182e-5
const a6 =  4.63e-7
const a7 = -4.227e-8

#####################
# Derived Constants #
#####################
const f = c/λ                           # Pulse Frequency Hz
const ω = 2*pi*f                        # Pulse Angular Frequency
const σ_t = Tfwhm/sqrt(2*log(2))        # 1-sigma width of pulse
const Power = sqrt(2/pi) * Energy/σ_t   # Max power delivered by pulse

const dt=tmax/(Nt-1)             # Time step
const points = (-Nt/2:1:Nt/2)
const t_vec = (-Nt/2:1:Nt/2)*dt  # Time grid iterator

const ff=points./tmax            # Frequency grid
const ωω=(2*pi)*ff               # Angular frequency grid

const λ_tot = 1E9 * c ./ (f + ff) # Peak centered Wavelength grid, in nm
const ωω_tot = ω+ωω               # Peak centered angular frequency grid

const n_tot_0 = 1+ C1 * (C2 * (λ_tot.^2) ./ (C3 * (λ_tot.^2) -1) +
                         C4 * (λ_tot.^2) ./ (C5 * (λ_tot.^2) -1) +
                         C6 * (λ_tot.^2) ./ (C7 * (λ_tot.^2) -1))

#= How is this to be translated?
l1min=λ_tot(abs(lambda_tot_rouge-lmin)==min(abs(lambda_tot_rouge-lmin)));
l1max=λ_tot(abs(lambda_tot_rouge-lmax)==min(abs(lambda_tot_rouge-lmax)));
This should do it I believe:=#
λ_test = minimum(abs.(λ_tot-λ_min))
const λ_min = filter(x -> abs(x - λ_min) == λ_test, λ_tot)[1]
λ_test = minimum(abs.(λ_tot-λ_max))
const λ_max = filter(x -> abs(x - λ_max) == λ_test, λ_tot)[1]

const ρ_crit = ω^2*me*ϵ0/ee^2
const k_Ar   = ceil.(Ui_Ar ./ (hb*ωω))

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
E                 = exp.(-t_vec.^2/σ_t^2)
E0                = sqrt(2*Power/(pi*Fiber_D^2))
E                 = E0 .* E
spectrum_entrance = abs.(fftshift(fft(fftshift(E)))).^2
I_entrance        = sum(spectrum_entrance)

# Propagation variables
zpoints       = 1000
dist          = zeros(zpoints)
ρ_max         = zeros(zpoints)
I_max         = zeros(zpoints)
ΔT_pulse      = zeros(zpoints)
It_dist       = zeros(Nt,zpoints)
spectrum_dist = zeros(Nt,zpoints)

Chirp_function = Chirp * ωω .^ 2 + TOD .* ωω .^ 3
E_TF           = fftshift(fft(fftshift(E))).*exp.(im * Chirp_function)
E              = ifftshift(ifft(ifftshift(E_TF)))

#Clear vars
Chirp_function = nothing; E_TF           = nothing; gc()

z   = 0
JJJ = 0
dz  = dz0
ZZZ = 0

idxp = 2:Nt
idxn = Nt:(Nt-1)

pressure = []

#=
███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
█████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████
=#
function calc_pressure(p_in, p_out, z, zmax)
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

function calc_duration(E, t1)
    center_pulse = sum(t1.*(abs.(E)).^2)/sum(abs.(E).^2)
    return 2 * sqrt(2*log(2)).*((sum((t1-center_pulse).^2.*abs.(E).^2)
                               / sum(abs.(E).^2)).^0.5)*1E15
end

function plasma(dt, α, ρ_at, Potentiel_Ar, E, C2_Ar, Nt)
    #Init
    ρ_Ar = zeros(size(E))

    for i in 1:(Nt-1)
        ρ_Ar[i+1] = ρ_Ar[i] +
                    dt * (-α*ρ_Ar[i]^2+Potentiel_Ar[i]*(ρ_at - ρ_Ar[i]) +
                    (C2_Ar * abs(E[i])^2)*ρ_Ar[i])
    end
    return ρ_Ar
end

function prop_lin(E, deriv_t_2, dz, losses)
    # Shift to frequency domain and compute linear propagation
    E_TF = fftshift(fft(fftshift(E))) .* exp.(1im*(deriv_t_2)*dz)
    # Shift back to time domain, compute losses and return
    return ifftshift(ifft(ifftshift(E_TF))).*exp.(-losses/2 * dz)
end

function prop_non_lin(E, rrr, ρ, dz, losses, kerr_response)
    return E.*exp(rrr*ρ*dz - losses + kerr_response)
end

function calc_compression(width, Cs)
    λ_test = minimum(abs.(λ_tot - 600))
    λ_begin = λ_tot[findfirst(x -> abs(x - 600) == λ_test & x, λ_tot)]

    λ_test = minimum(abs.(λ_tot - 6000))
    λ_final = λ_tot[findfirst(x -> abs(x - 6000) == λ_test && x > 0, λ_tot)]

    n_fs = sqrt.(1 + Cs[1]*(λ_tot).^2./((λ_tot).^2-Cs[2]^2) +
                     Cs[3]*(λ_tot).^2./((λ_tot).^2-Cs[4]^2) +
                     Cs[5]*(λ_tot).^2./((λ_tot).^2-Cs[6]^2))

    n_begin = n_fs[findfirst(x->x==λ_begin, λ_tot)]
    n_fs[λ_tot.<λ_begin] = n_begin

    n_final = n_fs[findfirst(x->x==λ_final, λ_tot)]
    n_fs[λ_tot.>λ_final] = n_final

    spline = Spline1D(ωω_tot, n_fs)
    n_fs_interpd = evaluate(n_fs_spline, ωω_tot)
    dn_fs_interpd = derivative(n_fs_spline, ωω_tot)
    k_FS=n_fs_interpd.*ωω_tot/c
    k_FS[find(x->x<245,λ_tot)] = maximum(k_FS)

    k_prime_FS=1 / c * (dn_fs_interpd.*ωω_tot+n_fs_interpd)

    k0_FS=k_FS[findfirst(x->x==ω,ωω_tot)]
    k1_FS=k_prime_FS[findfirst(x->x==ω,ωω_tot)]

    return [(k_FS-k0_FS-k1_FS.*ω).*(w.*1e-3) for w in width]
end

function smooth(values, radius)
    smoothed_values = zeros(values)
    for i in eachindex(values)
        temp_radius = minimum([radius, i - 1, length(values) - i])
        smoothed_values[i] = mean(values[(i - temp_radius):(i + temp_radius)])
    end
    return smoothed_values
end

function steepening(E, idxp, idxn, γs , ω_0, dt, dz)
    NL = (γ[1] * abs.(E)^2 + γ[2] * abs.(E)^4 + γ[3] * abs.(E)^6 +
          γ[4] * abs.(E)^8 + γ[5] * abs.(E)^10) * dz
    Temp = NL .* E
    k1 = (im / ω_0) * ((Temp[idxp] - Temp[idxn])/(2 * dt))
    E = E + (0.5 * k1)

    NL = (γ[1] * abs.(E)^2 + γ[2] * abs.(E)^4 + γ[3] * abs.(E)^6 +
          γ[4] * abs.(E)^8 + γ[5] * abs.(E)^10) * dz
    Temp = NL .* E
    k2 = (im / ω_0) * ((Temp[idxp] - Temp[idxn])/(2 * dt))

    return E + k2
end

function calc_ks(ωω, λ_tot, ωω_tot, n_tot)
    # Interpolation of n
    n_interp = Spline1D(ωω_tot, n_tot)
    n_tot = evaluate(n_interp, ωω_tot)
    dn    = derivative(n_interp, ωω_tot, nu=1)
    d2n   = derivative(n_interp, ωω_tot, nu=2)
    d3n   = derivative(n_interp, ωω_tot, nu=3)
    d4n   = derivative(n_interp, ωω_tot, nu=4)

    # Frequency Disperion
    k_tot = n_tot .* ωω_tot/c
    k_tot[λ_tot .< 245] = max(k_tot)

    k_first  = 1/c * (dn .* ωω_tot + n_tot)
    k_second = 1/c * (d2n .* ωω_tot + 2 * dn)
    k_third  = 1/c * (d3n .* ωω_tot + 3 * d2n)
    k_fourth = smooth(1/c * (d4n .* ωω_tot + 4 * d3n), 100)

    k  = k_tot[find_first(x->x==ωω,ωω_tot)]
    k1 = k_first[find_first(x->x==ωω,ωω_tot)]
    k2 = k_second[find_first(x->x==ωω,ωω_tot)]
    k3 = k_third[find_first(x->x==ωω,ωω_tot)]
    k4 = k_fourth[find_first(x->x==ωω,ωω_tot)]
    return [k,k1,k2,k3,k4]
end

function calc_ns(pressure, n, n_tot, λ_tot)
    n2_800  = 1e-23   * pression_z
    n4_800  =-3.7e-42 * pression_z
    n6_800  = 4e-58   * pression_z
    n8_800  =-1.7e-75 * pression_z
    n10_800 = 8.8e-94 * pression_z

    n_800 = n_tot[abs.(λ_tot - 800) == minimum(abs.(λ_tot-800))]
    n2  = n2_800 * ((n^2 - 1)/(n_800^2-1))^4
    n4  = n2_800 * ((n^2 - 1)/(n_800^2-1))^6
    n8  = n2_800 * ((n^2 - 1)/(n_800^2-1))^10
    n10 = n2_800 * ((n^2 - 1)/(n_800^2-1))^12

    return [n,n2,n4,n8,n10]
end

function plasma_potential(E,ω,Zeff,Ui)
    """
    Derived From PPT Theory
    http://jetp.ac.ru/cgi-bin/dn/e_023_05_0924.pdf
    """
    Uh = 13.5984*ee                     # Hydrogen Ionization Potential
    ω_au = 4.1E16                       # Ionization Potential (1/s Natural)
    γ = ω .* sqrt(2*me*Ui)./(ee*E)
    Eh = ee^5*me^2/(ħ^4 * (4*π*ϵ0)^3)
    E0 = Eh * (Ui/Uh) ^ (3/2)
    A = zeros(γ)
    β = 2*γ ./ sqrt.(1+gamma1.^2)
    α = 2.*asinh.(γ) - β

    g=3/(2*γ).*((1+1/(2*γ.^2)).*asinh.(γ)-1/beta)
    ν0 = Ui / (ħ*ω)
    ν  = ν0 * (1 + 1/(2*γ.^2))
    kmin = floor(nu) + 1

    l=0 #????????????????????????????????????
    m=0 #????????????????????????????????????
    n⋆ = Zeff * sqrt(Uh/Ui)

    C_nl2 = 2^(2*n⋆) / (n⋆*γ[2*n⋆]*γ[1])
    f = (2*l+1) * factorial(l+abs(m)) / (2^abs(m)) *
        factorial(abs(m)) * factorial(l-abs(m))

    for z in Kmin:(Kmin+2)
        A += 4/sqrt(3*π) * γ.^2 ./ (1+γ.^2) * exp(-α * (z-nu)) * dawson(sqrt(abs(beta*(z-nu))))
    end

    potential = ω_au * C_nl2 * f * sqrt(6/π) * (Ui / (2 * Uh)) * A .* (2 * E0/(E .* sqrt.(1+γ.^2))).^(2 * n⋆ - abs(m) - 3/2) .* exp.(-2 * E0 * g ./ (3*E))
    potential[isnan.(potential)]=0
end

#=
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
=#

#Temporary Flag
run = false
while z <= zmax && run
    # Calculate Pressure
    pressure_z = calc_pressure(0.008, Pressure, z, zmax)
    push!(pressure, pressure_z)

    # Update of n, with cutoffs
    n_tot = sqrt(1+pressure_z * (n_tot_0.^2 - 1))
    n_max = n_tot[findfirst(λ->λ==λ_max, λ_tot)]
    n_tot[λ_tot .> λ_max] = n_max
    n_min = n_tot[findfirst(λ->λ==λ_min, λ_tot)]
    n_tot[λ_tot .< λ_min] = n_min

    ks = calc_ks()
    n  = n_tot[find_first(x->x==ωω,ωω_tot)]
    vg = 1/ks[2]

    # Argon Parameters
    ns = calc_ns()
    β2 = pressure_z * k[3]
    β3 = pressure_z * k[4]
    β4 = pressure_z * k[5]
    τ = 3.5E-13 / pressure_z
    ρ_at = pressure_z * 1E5 / (kB * T)

    # Plasma Parameters
    σ_k = 2.81E-96 * pressure
    σ   = (k[1]*ee^2) / (ωω * m * ϵ0) * τ/(1+(ωω * τ).^2)
    β_k = 10^(-4 * k_Ar) * k_Ar * ħ * ωω * ρ * 0.21 * σ_k
    rrr = -im * k[1]/(2 * n[1]^2 * ρ_crit) - 0.5 * σ
    coeff2 = σ/Ui_Ar

    # Dispersion and Laplacian Operators
    dv_t_2_op = k_tot - k[1] - k[2] * ωω

    # Kerr factors
    γs = im * n[2:end] * k[0]/n[0]

    # Propagation
    E = prop_lin(E, dv_t_2_op, dz, losses)        #Linear
    E = steepening(E, idxp, idxn, γs, ωω, dt, dz) #Steepening

    # Plasma
    U_ion = PPT(abs(E).^2, c, ωω)
    ρ = plasma(dt, α, ρ_at, U_ion, E, coeff2, Nt)
    plasma_loss = U_ion / (2 * abs.(E).^2) * Ui_Ar * (ρ_at - ρ) * dz
    plasma_loss[isnan.(plasma_loss)] = 0

    # Kerr and Plasma Propagation (NonLinear)
    kerr_response = -(γs[1]*(abs.(E)).^2 + γs[2]*(abs.(E)).^4 + γs[3] * abs.(E).^6 + γs[4]*abs.(E).^8 + γ[5].*abs.(E).^10) * dz
    E = prop_non_lin(E, rrr, ρ, dz, plasma_loss, kerr_response)

    # Update Distances
    ZZZ += dz
    z   += dz
    distance = 100*ZZZ
    dist[1, JJJ+1] = z

    # Calculating Results
    I_max[JJJ+1]         = maximum(abs.(E).^2)
    It_dist[JJJ+1]       = abs.(E).^2
    spectrum_dist[JJJ+1] = abs.(fftshift(fft(fftshift(E)))).^2
    ΔT_pulse[JJJ+1]      = calc_duration(E,t)
    ρ_max[JJJ+1]         = maximum(ρ*1E-6)
    #TODO:Save intermediate data
    JJJ += 1
end
