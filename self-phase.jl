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
const Energy   = # Pulse Energy       J
const Tfwhm    = # Pulse Width        s
const Chirp    = # Pulse Chirp        s
const λ        = # Pulse Wavelength   m

const Fiber_D  = # Fiber Diameter     m
const Fiber_L  = # Fiber Length       m
const Pressure = # Capillary Pressure Pa

const TOD      = # ??
const losses   = # ??

const Ui_Ar   = 15.75*1.6E-19 # Ionization Energy of Ar J
const α       = 7E-13
const Zeff_Ar = 1

# Numeric Parameters
const dz0=2E-3       # z-step 1
const dz1=2E-4       # z-step 2
const zmax = Fiber_L # Total z

const Nt=2*8192          # Number of time steps
const tmax=1200E-15      # Maximum time

const lmin = 400
const lmax = 1500
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
const T  = 300              # ~Room Temperature  K

# Indices of refraction/dispersion for Argon
const C1 = 0.012055
const C2 = 0.2075
const C3 = 91.012
const C4 = 0.0415
const C5 = 87.892
const C6 = 4.3330
const C7 = 214.02
const C = [C1,C2,C3,C4,C5,C6,C7]

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
const t_vec = collect((1:Nt)*dt) # Vector of times

const ff=t_vec./tmax
const ωω=(2*pi)*FF

const λ_tot = 1E9 * c ./ (f + ff)
const ωω_tot = ω+ωω

const n_tot_0 = 1+ C1 * (C2 * (λ_tot*1E3.^2) ./ (C3 * (λ_tot*1E3.^2) -1) +
                         C4 * (λ_tot*1E3.^2) ./ (C5 * (λ_tot*1E3.^2) -1) +
                         C6 * (λ_tot*1E3.^2) ./ (C7 * (λ_tot*1E3.^2) -1))

#= How is this to be translated?
l1min=λ_tot(abs(lambda_tot_rouge-lmin)==min(abs(lambda_tot_rouge-lmin)));
l1max=λ_tot(abs(lambda_tot_rouge-lmax)==min(abs(lambda_tot_rouge-lmax)));
This should do it I believe:=#
const l1min = λ_tot[(abs(λ_tot - lmin) == min(abs(λ_tot-lmin))) + 1]
const l1max = λ_tot[(abs(λ_tot - lmax) == max(abs(λ_tot-lmax))) + 1]

const ρ_crit = ωω.^2*m*ϵ0/ee^2
const k_Ar   = ceil(Ui_Ar / (hbar*ωω))

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
E                 = exp(-t_vec.^2/σ_t^2)
E0                = sqrt(2*Power/(pi*Fiber_D^2))
E                 = E0 .* E
spectrum_entrance = abs(fftshift(fft(fftshift(E)))).^2
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
E_TF           = fftshift(fft(fftshift(E))).*exp(im * Chirp_function)
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
    return sqrt(p_in^2 + (z/zmax)*(p_out^2-p_in^2))
end

function calc_duration(E, t1)
    center_pulse = sum(t1.*(abs(E)).^2)/sum(abs(E).^2) #Remove . from /?
    return 2 * sqrt(2*log(2)).*((sum((t1-center_pulse).^2.*abs(E).^2)
                               / sum(abs(E).^2)).^0.5)*1E15
end

function plasma(dt, α, ρ_at, Potentiel_Ar, E, C2_Ar, Nt)
    #Init
    ρ_Ar = zeros(size(E))

    for i in 1:Nt-1
        ρ_Ar[i+1] = ρ_Ar[i]
                    + dt * (-α*ρ_Ar[i]^2+Potentiel_Ar[i]*(ρ_at - ρ_ar[i])
                            + (C2_Ar * abs(E[q])^2)*ρ_ar[i])
    end

    return ρ_Ar
end

function prop_lin(deriv_t_2, dz1, E, losses)
    # Shift to frequency domain and compute linear propagation
    E_TF = fftshift(fft(fftshift(E))) .* exp(1im*(deriv_t_2)*dz1)
    # Shift back to time domain, compute losses and return
    return ifftshift(ifft(ifftshift(E_TF))).*exp(-losses/2 * dz1)
end

function prop_non_lin(E, rrr, ρ, dz, losses, kerr_response)
    return E.*exp(rrr*ρ*dz - losses + kerr_response)

#=
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
=#
while z <= zmax
    pressure = calc_pressure(0.008, Pressure, z, zmax)
end
