#----------------------------------------------------------------------------------------------------------------------------------------
# Implements solutions to the diffusion equation for the fluence in the steady-state and time-domain in a semi-infinite medium using
# equations 1, 3, 5, and 6 from Kienle and Patterson [1].
#
# [1] Alwin Kienle and Michael S. Patterson, "Improved solutions of the steady-state and the time-resolved diffusion equations 
#     for reflectance from a semi-infinite turbid medium," J. Opt. Soc. Am. A 14, 246-254 (1997) 
#----------------------------------------------------------------------------------------------------------------------------------------
"""
    fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)

Compute the steady-state fluence in a semi-infinite geometry.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)

# Examples
```
julia> fluence_DA_semiinf_CW(1.0, 0.1, 10.0)
julia> fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0)
```
"""
function fluence_DA_semiinf_CW end

"""
    flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)

Compute the steady-state flux (D*∂ϕ(ρ)/∂z @ z = 0) from a semi-infinite medium.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction

# Examples
```
julia> flux_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)
```
"""
function flux_DA_semiinf_CW end

"""
fluence_DA_semiinf_FD(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)

Compute the frequency domain fluence in a semi-infinite geometry.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)
- `ω`: the modulation frequency (1/ns)

# Examples
```
julia> fluence_DA_semiinf_FD(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)
```
"""
function fluence_DA_semiinf_FD end

"""
    fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)

Compute the time-domain fluence in a semi-infinite medium.  

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)

# Examples
```
julia> fluence_DA_semiinf_TD(0.1:0.1:2.0, 1.0, 0.1, 10.0)
julia> fluence_DA_semiinf_TD(0.2:0.2:1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0)
```
"""
function fluence_DA_semiinf_TD end

"""
    flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)

Compute the time-domain flux (D*∂ϕ(t)/∂z @ z = 0) in a semi-infinite medium.

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction

# Examples
```
julia> flux_DA_semiinf_TD(0.1:0.1:2.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)
```
"""
function flux_DA_semiinf_TD end

#------------------------------
# Steady-State Fluence 
#------------------------------

function fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    μeff = sqrt(3 * μa * μsp)
    D, A = D_coeff(μsp), A_coeff(n_med / n_ext)
    z0, zb = z0_coeff(μsp), zb_coeff(A, D)
    return _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
end
function _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
    a = hypot(ρ, z - z0)
    b = hypot(ρ, z + 2 * zb + z0)

    # the following tries to compute ϕ = exp(-μeff * a) / a - exp(-μeff * b) / b more accurately
    # h also tries to compute h = b - a more accurately
    h = 4 * (z0 * z + z0 * zb + z * zb + zb^2) / (a + b)
    ϕ = exp(-b * μeff) * fma(b, expm1(h * μeff), h) / (b * a)
    return ϕ / (4 * π * D)
end
function _kernel_fluence_DA_semiinf_CW(μeff::Complex, ρ, z, z0, zb, D)
    a = hypot(ρ, z - z0)
    b = hypot(ρ, z + 2 * zb + z0)
    ϕ = exp(-μeff * a) / a - exp(-μeff * b) / b 
    return ϕ / (4 * π * D)
end

#------------------------------
# Steady-State Flux
#------------------------------

function flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    T = eltype(ρ)

    μeff = sqrt(3 * μa * μsp)
    D, A = D_coeff(μsp), A_coeff(n_med / n_ext)
    z0, zb = z0_coeff(μsp), zb_coeff(A, D)

    a = hypot(z0, ρ)
    b = hypot(z0 + 2*zb, ρ)
    ϕ = z0 * (μeff + inv(a)) * exp(-μeff * a) / a^2
    ϕ += (z0 + 2*zb) * (μeff + inv(b)) * exp(-μeff * b) / b^2

    return ϕ / (4 * T(π))
end

#------------------------------
# Frequency-Domain Fluence 
#------------------------------

function fluence_DA_semiinf_FD(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)
    D, A, ν = D_coeff(μsp), A_coeff(n_med / n_ext), ν_coeff(n_med)
    z0, zb = z0_coeff(μsp), zb_coeff(A, D)

    μa_complex = μa + ω * im / ν
    μeff = sqrt(3 * μa_complex * μsp)
    return _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
end

#------------------------------
# Time-Domain Fluence 
#------------------------------

function fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    T = eltype(ρ)
    D, A, ν = D_coeff(μsp), A_coeff(n_med / n_ext), ν_coeff(n_med)
    z0, zb = z0_coeff(μsp), zb_coeff(A, D)

    b = ρ^2 + (z + z0 + 2zb)^2
    h = 4 * (z * z0 + z * zb + zb^2 + z0 * zb)
    ϕ = ν * inv(T(pi)^(1.5))

    return map(t -> _kernel_fluence_DA_semiinf_TD(ϕ, b, h, D, ν, t, μa), t)
end
@inline function _kernel_fluence_DA_semiinf_TD(ϕ, b, h, D, ν, t, μa)
    tmp1 = inv(4 * D * ν * t)
    ϕ *= exp(-μa * ν * t) * tmp1^(1.5)
    ϕ *= exp(-b * tmp1) * expm1(tmp1 * h)
    return ϕ
end

#------------------------------
# Time-Domain Flux
#------------------------------

function flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    D, A, ν = D_coeff(μsp), A_coeff(n_med / n_ext), ν_coeff(n_med)
    z0, zb = z0_coeff(μsp), zb_coeff(A, D)

    a = hypot(z0, ρ)^2
    b = hypot(z0 + 2 * zb, ρ)^2

    return  map(t -> _kernel_flux_DA_semiinf_TD(D, ν, t, μa, a, b, z0, zb), t)
end

@inline function _kernel_flux_DA_semiinf_TD(D, ν, t, μa, a, b, z0, zb)
    tmp = inv(4 * D * ν * t)
    ϕ = z0 * exp(-a * tmp)
    ϕ += (z0 + 2 * zb) * exp(-b * tmp)
    ϕ *= exp(-μa * ν * t)
    ϕ /= 2 * (4 * D * ν * π)^(3/2) * t^(5/2)

    return ϕ 
end
