#----------------------------------------------------------------------------------------------------------------------------------------
# Implements solution to the diffusion equation for the fluence in the steady-state and time-domain in a semi-infinite medium [1][2][3].
#
# [1] Alwin Kienle and Michael S. Patterson, "Improved solutions of the steady-state and the time-resolved diffusion equations 
#     for reflectance from a semi-infinite turbid medium," J. Opt. Soc. Am. A 14, 246-254 (1997) 
# [2] Martelli, Fabrizio, et al. Light propagation through biological tissue and other diffusive media: theory, solutions and software. 
#     Vol. 10. No. 3.824746. Bellingham: SPIE press, 2010.
#----------------------------------------------------------------------------------------------------------------------------------------

""" Structure containing inputs for simulating the fluence under the diffusion approximation in the semi-infinite space."""
struct DAsemiinf_params{T <: AbstractFloat} <: DiffusionParameters
    ρ::T                                 # distance away from isotropic point source (cm)
    μa::T                                # absorption coefficient (cm⁻¹)
    μsp::T                               # reduced scattering coefficient (cm⁻¹)
    n_med::T                             # medium's index of refraction
    n_ext::T                             # external medium's index of refraction
    ω::T                                 # modulation frequency (1/ns)
    z::T                                 # z depth of detector within medium (cm)
    t::Union{T, AbstractVector{T}}       # time vector (ns)

    # do not need to provide these arguments (calculated from previous inputs)
    D::T                                 # Diffusion coefficient                        
    ν::T                                 # Speed of light in medium (cm/ns)
    z0::T                                # isotroptic source depth
    zb::T                                # Extrapolation length
    function DAsemiinf_params{T}(ρ::T, μa::T, μsp::T, n_med::T, n_ext::T, ω::T, z::T, t::Union{T, AbstractVector{T}}) where {T <: AbstractFloat}
        @assert ρ >= zero(T) "ρ must greater than or equal to 0"
        @assert μa >= zero(T) "μa must greater than or equal to 0"
        @assert μsp > zero(T) "μsp must greater than 0"
        @assert n_med > zero(T) "n_med must be positive"
        @assert n_ext > zero(T) "n_ext must be positive"
        @assert all(t .> zero(T)) "t must be positive"
        D = D_coeff(μsp)
        return new{T}(ρ, μa, μsp, n_med, n_ext, ω, z, t, D, ν_coeff(n_med), z0_coeff(μsp), zb_coeff(A_coeff(n_med / n_ext), D))
    end
end

"""
    Generator function for DAsemiinf_params structure.

Provides default parameters to use in the semi-infinite space fluence calculation in either the CW, FD, or TD.

# Keyword Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_med`: medium's index of refraction
- `n_ext`: external medium's index of refraction (air or detector)
- `ω`: modulation frequency (1/ns)
- `z`: the z-depth orthogonal from the boundary (cm)
- `t`: the time vector (ns). 

# Examples
```
julia> data = DAsemiinf_params() # return default parameters
julia> data = DAsemiinf_params(ρ = 1.5) # return ρ = 1.5 with the rest of the parameters given by defaults
julia> fluence_DA_semiinf_CW(data) # can then be used with corresponding functions
```
"""
function DAsemiinf_params(;
    ρ::T = 1.0,
    μa::T = 0.1,
    μsp::T = 10.0,
    n_med::T = 1.0,
    n_ext::T = 1.0,
    ω::T = 0.0,
    z::T = 0.0,
    t::Union{T, AbstractVector{T}} = 1.0
) where {T<:AbstractFloat}
    return DAsemiinf_params{T}(ρ, μa, μsp, n_med, n_ext, ω, z, t)
end

#------------------------------
# Steady-State Fluence 
#------------------------------

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
function fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    μeff = sqrt(3 * μa * μsp)
    D = D_coeff(μsp)
    (z0, zb) = (z0_coeff(μsp), zb_coeff(A_coeff(n_med / n_ext), D))
    return _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
end
function _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
    a = hypot(ρ, z - z0)
    b = hypot(ρ, z + 2 * zb + z0)

    # the following tries to compute ϕ = exp(-μeff * a) / a - exp(-μeff * b) / b more accurately
    # there is some loss of significance when a ≈ b and a > 5
    h = b - a
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
# Time-Domain Fluence 
#------------------------------

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
function fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t -> _kernel_fluence_DA_semiinf_TD(params.D, params.ν, t, μa, z, params.z0, ρ, params.zb), t)
end

"""
    fluence_DA_semiinf_TD(t, data::DiffusionParameters)

Wrapper to `fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)`

# Examples
```
julia> data = DAsemiinf_params() # use structure to generate inputs
julia> fluence_DA_semiinf_TD(data) # then call the function
```
"""
function fluence_DA_semiinf_TD(data::DiffusionParameters)
    return map(t -> _kernel_fluence_DA_semiinf_TD(data.D, data.ν, t, data.μa, data.z, data.z0, data.ρ, data.zb), data.t)
end

@inline function _kernel_fluence_DA_semiinf_TD(D, ν, t, μa, z, z0, ρ, zb)
    tmp1 = 4 * D * ν * t
    ϕ = ν * exp(-μa * ν * t) / (π * tmp1)^(3/2)
    ϕ *= (exp(-((z - z0)^2 + ρ^2) / tmp1) - exp(-((z + z0 + 2*zb)^2 + ρ^2) / tmp1))
    return ϕ
end

#------------------------------
# Frequency-Domain Fluence 
#------------------------------
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
function fluence_DA_semiinf_FD(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    μa_complex = μa + ω * im / params.ν
    μeff = sqrt(3 * μa_complex * μsp)
    return _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, params.z0, params.zb, params.D)
end

#------------------------------
# Time-Domain Flux
#------------------------------
"""
    flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)

Compute the time-domain flux (D*∂ϕ(t)/∂z @ z = 0) from a semi-infinite medium from Eqn. 36 Contini 97 (m = 0).

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
function flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    z3m = - params.z0
    z4m = 2 * params.zb + params.z0

    return map(t -> _kernel_flux_DA_semiinf_TD(params.D, params.ν, t, μa, z3m, z4m, ρ), t)
 end

 """
    flux_DA_semiinf_TD(t, data::DiffusionParameters)

Wrapper to `flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)`

# Examples
```
julia> data = DAsemiinf_params() # use structure to generate inputs
julia> flux_DA_semiinf_TD(data) # then call the function
```
"""
function flux_DA_semiinf_TD(data::DiffusionParameters)
    z3m = - data.z0
    z4m = 2 * data.zb + data.z0

    return map(t -> _kernel_flux_DA_semiinf_TD(data.D, data.ν, t, data.μa, z3m, z4m, data.ρ), data.t)
end

@inline function _kernel_flux_DA_semiinf_TD(D, ν, t, μa, z3m, z4m, ρ)
    tmp1 = 4 * D * ν * t
    Rt = -(ρ^2 / (tmp1))
    Rt -= μa * ν * t
    Rt = -exp(Rt) / (2 * (4 * π * D * ν)^(3/2) * sqrt(t * t * t * t * t))
    Rt *= (z3m * exp(-(z3m^2 / (tmp1))) - z4m * exp(-(z4m^2 / (tmp1))))

    return Rt
end

#------------------------------
# Steady-State Flux
#------------------------------
"""
    flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)

Compute the steady-state flux (D*∂ϕ(ρ)/∂z @ z = 0) from a semi-infinite medium.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction

# Examples
```
julia> flux_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)
```
"""
function flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    D = D_coeff(μsp)
    return D * ForwardDiff.derivative(z -> fluence_DA_semiinf_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, z = z), 0.0)
end

"""
    flux_DA_semiinf_CW(data::DiffusionParameters)

Wrapper to `flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)`

# Examples
```
julia> data = DAsemiinf_params(ρ = 1.0) # use structure to generate inputs
julia> flux_DA_semiinf_CW(data) # then call the function
```
"""
function flux_DA_semiinf_CW(data)
    return flux_DA_semiinf_CW(data.ρ, data.μa, data.μsp, n_ext = data.n_ext, n_med = data.n_med)
end
