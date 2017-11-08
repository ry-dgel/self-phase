using Revise
using Dierckx
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
const hb = 1.0545718E-34    # Reduced Planck     Js
const me = 9.10938356E-31   # Electron Mass      kg
const ee = 1.6021766208E-19 # Elementary Charge  C
const ϵ0 = 8.854187817E-12  # Vaccum Permitivity F/m
const Kb = 1.38064852E-23   # Boltzmann Constant J/K
const T  = 300              # ~Room Temperature  Ka generalized theory of flexa

# Indices of refraction/dispersion for Argon
const C1 = 0.012055
const C2 = 0.2075
const C3 = 91.012
const C4 = 0.0415
const C5 = 87.892
const C6 = 4.3330
const C7 = 214.02
const C = [C1,C2,C3,C4,C5,C6,C7]

#Indices of Refraction for Fused Silica
const Cfs1 = 0.6961663;
const Cfs2 = 0.0684043;
const Cfs3 = 0.4079426;
const Cfs4 = 0.1162414;
const Cfs5 = 0.8974794;
const Cfs6 = 9.896161;
const Cfs = [Cfs1,Cfs2,Cfs3,Cfs4,Cfs5, Cfs6]

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
const t_vec = (-Nt/2:1:Nt/2)*dt  # Iteration of timesteps

const ff=points./tmax
const ωω=(2*pi)*ff

const λ_tot = 1E9 * c ./ (f + ff)
const ωω_tot = ω+ωω

const n_tot_0 = 1+ C1 * (C2 * (λ_tot*1E3.^2) ./ (C3 * (λ_tot*1E3.^2) -1) +
                         C4 * (λ_tot*1E3.^2) ./ (C5 * (λ_tot*1E3.^2) -1) +
                         C6 * (λ_tot*1E3.^2) ./ (C7 * (λ_tot*1E3.^2) -1))

#= How is this to be translated?
l1min=λ_tot(abs(lambda_tot_rouge-lmin)==min(abs(lambda_tot_rouge-lmin)));
l1max=λ_tot(abs(lambda_tot_rouge-lmax)==min(abs(lambda_tot_rouge-lmax)));
This should do it I believe:=#
λ_test = minimum(abs.(λ_tot-λ_min))
const λ1_min = filter(x -> abs.(λ_tot - λ_min) == λ_test, λ_tot)
λ_test = minimum(abs.(λ_tot-λ_max))
const λ1_max = filter(x -> abs.(λ_tot - λ_max) == λ_test, λ_tot)

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

function prop_lin(deriv_t_2, dz1, E, losses)
    # Shift to frequency domain and compute linear propagation
    E_TF = fftshift(fft(fftshift(E))) .* exp.(1im*(deriv_t_2)*dz1)
    # Shift back to time domain, compute losses and return
    return ifftshift(ifft(ifftshift(E_TF))).*exp.(-losses/2 * dz1)
end

function prop_non_lin(E, rrr, ρ, dz, losses, kerr_response)
    return E.*exp(rrr*ρ*dz - losses + kerr_response)
end

function calc_compression(width, Cs)
    λ_test = minimum(abs.(λ_tot - 600))
    λ_begin = filter(x -> abs(x - 600) == λ_test && x > 0, λ_tot)[1]
    println(λ_begin)
    λ_test = minimum(abs.(λ_tot - 6000))
    λ_final = filter(x -> abs(x - 6000) == λ_test && x > 0, λ_tot)[1]
    println(λ_final)
    n_fs = sqrt.(1 + Cs[1].*(λ_tot).^2./((λ_tot).^2-Cs[2].^2) +
                     Cs[3].*(λ_tot).^2./((λ_tot).^2-Cs[4].^2) +
                     Cs[5].*(λ_tot).^2./((λ_tot).^2-Cs[6].^2))
    println(mean(n_fs))
    n_begin = n_fs[find(x->x==λ_begin, λ_tot)][1]
    n_fs[find(x->x<λ_begin, λ_tot)] = n_begin
    println(mean(n_fs))
    n_final = n_fs[find(x->x==λ_final, λ_tot)][1]
    n_fs[find(x->x>λ_final, λ_tot)] = n_final
    println(mean(n_fs))
    n_interp_fs = Spline1D(ωω_tot, n_fs)
    n_fs_tot = evaluate(n_interp_fs, ωω_tot) #Is this not the same as n_fs??
    dn_fs_tot = derivative(n_interp_fs, ωω_tot)
    k_FS=n_fs_tot.*ωω_tot/c
    k_FS[find(x->x<245,λ_tot)]=maximum(k_FS)

    k_prime_FS=1./c.*(dn_fs_tot.*ωω_tot+n_fs_tot)

    k0_FS=k_FS[findfirst(x->x==ω,ωω_tot)]
    k1_FS=k_prime_FS[findfirst(x->x==ω,ωω_tot)]

    return [(k_FS-k0_FS-k1_FS.*ω).*(w.*1e-3) for w in width]
end

function smooth(values, radius)
    smoothed_values = zeros(size(values))
    smoothed_values[1], smoothed_values[end] = values[1], values[end]
    for i in 2:(length(values) - 1)
        temp_radius = minimum([radius, i - 1, length(values) - i])
        smoothed_values[i] = mean(values[(i - temp_radius):(i + temp_radius)])
    end
    return smoothed_values
end

function steepening(E1, idxp, idxn, γ_2, γ_4, γ_6, γ_8, γ_10, ω_0, dt, dz)
    NL = (γ_2 * abs.(E1)^2 + γ_4 * abs.(E1)^4 + γ_6 * abs.(E1)^6 +
        γ_8 * abs.(E1)^8 + γ_10 * abs.(E1)^10) * dz
    Temp = NL .* E1
    k1 = (im / ω_0) * ((Temp[idxp] - Temp[idxn])/(2 * dt))
    E_temp = 0
end

#=
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
=#
#while z <= zmax
#    pressure = calc_pressure(0.008, Pressure, z, zmax)
#end
