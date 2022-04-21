#------------------------------------------------------------------------------------------------------------------------
# Implements solution to the diffusion equation in an infinite medium as given in [1].
# Solutions are given in the steady-state, frequency, and time domains for an isotroptic source.
#
# [1] Patterson et. al., "Time resolved reflectance and transmittance for the noninvasive measurement of tissue optical properties," 
#     Appl. Opt. 28, 2331-2336 (1989)
#-------------------------------------------------------------------------------------------------------------------------
"""
    fluence_DA_inf_CW(ρ, μa, μsp)

Compute the steady-state fluence in an infinite medium.

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Examples
```
julia> fluence_DA_inf_CW(1.0, 0.1, 10.0)
```
"""
function fluence_DA_inf_CW end

"""
    fluence_DA_inf_FD(ρ, μa, μsp; ω = 0.0, n_med = 1.0)

Compute the fluence for a frequency modulated source in an infinite medium. 

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword Arguments
- `ω`: modulation frequency (1/ns)
- `n_med`: medium's index of refraction

# Examples
```
julia> fluence_DA_inf_FD(1.0, 0.1, 10.0, ω = 1.0, n_med = 1.4)
```
"""
function fluence_DA_inf_FD end

"""
    fluence_DA_inf_TD(t, ρ, μa, μsp; n_med = 1.0)

Compute the time-domain fluence in an infinite medium.

# Arguments
- `t`: the time vector (ns). 
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword Arguments
- `n_med`= 1.0: medium's index of refraction

# Examples
```
julia> fluence_DA_inf_TD(0.1, 1.0, 0.1, 10.0, n_med = 1.4)
julia> fluence_DA_inf_TD(0.1:0.5:5.0, 1.0, 0.1, 10.0, n_med = 1.4)
```
"""
function fluence_DA_inf_TD end


#-------------------------------------
# Steady-State Fluence 
#-------------------------------------

function fluence_DA_inf_CW(ρ, μa, μsp)
    D = D_coeff(μsp)
    return exp(-sqrt(3 * μsp * μa) * ρ) / (4 * π * ρ * D)
end

#-------------------------------------
# Frequency-Domain Fluence 
#-------------------------------------

function fluence_DA_inf_FD(ρ, μa, μsp; ω = 0.0, n_med = 1.0)
    ν = ν_coeff(n_med)
    μa_complex = μa + ω * im / ν
    return fluence_DA_inf_CW(ρ, μa_complex, μsp)
end

#-------------------------------------
# Time-Domain Fluence 
#-------------------------------------

function fluence_DA_inf_TD(t, ρ, μa, μsp; n_med = 1.0)
    D, ν = D_coeff(μsp), ν_coeff(n_med)
    return map(t -> _kernel_fluence_DA_inf_TD(t, ρ, μa, ν, D), t)
end
@inline function _kernel_fluence_DA_inf_TD(t, ρ, μa, ν, D)
    tmp1 = 4 * D * ν * t
    ϕ = exp(-(ρ^2 / tmp1) - μa * ν * t)
    ϕ *= ν / ((tmp1 * π )^(3/2))
end
