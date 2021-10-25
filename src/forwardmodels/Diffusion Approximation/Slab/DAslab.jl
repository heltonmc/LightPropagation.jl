#----------------------------------------------------------------------------------------------------------------------------------------
 
# Implements solution to the diffusion equation for the fluence in the steady-state and time-domain in a slab [1][2].
#
# [1] Contini et. al. "Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory." 
#     Applied optics 36.19 (1997): 4587-4599.
# [2] Martelli et al. Light propagation through biological tissue and other diffusive media: theory, solutions and software. 
#     Vol. 10. No. 3.824746. Bellingham: SPIE press, 2010.
#----------------------------------------------------------------------------------------------------------------------------------------
""" Structure containing inputs for simulating the fluence under the diffusion approximation in the slab space."""
struct DAslab{T <: AbstractFloat} <: DiffusionParameters
    ρ::T                                 # distance away from isotropic point source (cm)
    μa::T                                # absorption coefficient (cm⁻¹)
    μsp::T                               # reduced scattering coefficient (cm⁻¹)
    n_med::T                             # medium's index of refraction
    n_ext::T                             # external medium's index of refraction
    ω::T                                 # modulation frequency (1/ns)
    z::T                                 # z depth of detector within medium (cm)
    t::Union{T, AbstractVector{T}}       # time vector (ns)
    s::T                                 # slab thickness (cm) in z direction
    xs::Int                              # the number of sources to compute in the series

    # do not need to provide these arguments (calculated from previous inputs)
    D::T                                 # Diffusion coefficient                        
    ν::T                                 # Speed of light in medium (cm/ns)
    z0::T                                # isotroptic source depth
    zb::T                                # Extrapolation length
    function DAslab{T}(ρ::T, μa::T, μsp::T, n_med::T, n_ext::T, ω::T, z::T, t::Union{T, AbstractVector{T}}, s::T, xs::Int) where {T <: AbstractFloat}
        @assert ρ >= zero(T) "ρ must greater than or equal to 0"
        @assert μa >= zero(T) "μa must greater than or equal to 0"
        @assert μsp > zero(T) "μsp must greater than 0"
        @assert s > zero(T) "s must greater than 0"
        @assert xs > 0 "xs must greater than 0"
        @assert all(t .> zero(T)) "t must be positive"
        D = D_coeff(μsp)
        return new{T}(ρ, μa, μsp, n_med, n_ext, ω, z, t, s, xs, D, ν_coeff(n_med), z0_coeff(μsp), zb_coeff(A_coeff(n_med / n_ext), D))
    end
end

"""
    Generator function for DAslab structure.

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
- `s`: the thickness (z) of the slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
```
julia> data = DAslab() # return default parameters
julia> data = DAslab(ρ = 1.5) # return ρ = 1.5 with the rest of the parameters given by defaults
julia> fluence_DA_slab_CW(data) # can then be used with corresponding functions
```
"""
function DAslab(;
    ρ::T = 1.0,
    μa::T = 0.1,
    μsp::T = 10.0,
    n_med::T = 1.0,
    n_ext::T = 1.0,
    ω::T = 0.0,
    z::T = 0.0,
    t::Union{T, AbstractVector{T}} = 1.0,
    s::T = 2.0,
    xs::Int = 10
) where {T<:AbstractFloat}
    return DAslab{T}(ρ, μa, μsp, n_med, n_ext, ω, z, t, s, xs)
end

#--------------------------------------
# Steady-State Fluence 
#--------------------------------------
"""
    fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)

Compute the steady-state fluence from a slab geometry (x, y -> inf, z -> finite).
The fluence is calculated from an infinite summation that converges rather rapidly.
Generally only 5 terms are needed (xs = 5), but more terms will be required if ρ >> s.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth within slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
```
julia> fluence_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0)
```
"""
function fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0, xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    μeff = sqrt(3 * μa * μsp)
    return _sum_fluence_DA_slab_CW(xs, s, params.zb, params.z0, ρ, z, μeff, params.D)
end

"""
    fluence_DA_slab_CW(data::DiffusionParameters)

Wrapper to `fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0, xs = 10)`

# Examples
```
julia> data = DAslab(ρ = 1.0, s = 3.0) # use structure to generate inputs
julia> fluence_DA_slab_CW(data) # then call the function
```
"""
function fluence_DA_slab_CW(data::DiffusionParameters)
    μeff = sqrt(3 * data.μa * data.μsp)
    return _sum_fluence_DA_slab_CW(data.xs, data.s, data.zb, data.z0, data.ρ, data.z, μeff, data.D)
end

function _sum_fluence_DA_slab_CW(xs, s, zb, z0, ρ, z, μeff, D)
    ϕ = zero(eltype(μeff))

    for m in xs:-1:1
        ϕ += _kernel_fluence_DA_slab_CW(-m, s, zb, z0, ρ, z, μeff)
        ϕ += _kernel_fluence_DA_slab_CW(m, s, zb, z0, ρ, z, μeff)
    end
    ϕ += _kernel_fluence_DA_slab_CW(0, s, zb, z0, ρ, z, μeff)

    return ϕ / (4 * π * D)
end

@inline function _kernel_fluence_DA_slab_CW(m, s, zb, z0, ρ, z, μeff)
    tmp1 = 2 * m * (s + 2 * zb)
    zmp = tmp1 + z0
    zmm = tmp1 - 2 * zb - z0

    tmp2 = sqrt(ρ^2 + (z - zmp)^2)
    tmp3 = sqrt(ρ^2 + (z - zmm)^2)

    ϕ = exp(-μeff * tmp2) / tmp2
    ϕ -= exp(-μeff * tmp3) / tmp3
    return ϕ
end

#--------------------------------------
# Time-Domain Fluence 
#--------------------------------------
"""
    fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)

Compute the time-domain fluence from a slab geometry (x, y -> inf, z -> finite).
The fluence is calculated from an infinite summation that converges rather rapidly.
Generally only 5 terms are needed (xs = 5), but more terms will be required if ρ >> s.
At long time values (low fluence values), the solution struggles with limited precision.
Using double floats or BigFloats is advised at these times.

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth coordinate (cm)
- `xs`: the number of sources to compute in the series

# Examples
julia> `fluence_DA_slab_TD(0.1, 1.0, 0.1, 10.0, s = 40.0)`
"""
function fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t-> _kernel_fluence_DA_slab_TD(params.D, params.ν, t, ρ, μa, s, params.zb, params.z0, z, xs), t)
end

"""
    fluence_DA_slab_TD(t, data::DiffusionParameters)

Wrapper to `fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)`

# Examples
```
julia> data = DAslab() # use structure to generate inputs
julia> fluence_DA_slab_TD(data) # then call the function
```
"""
function fluence_DA_slab_TD(data::DiffusionParameters)
    return map(t -> _kernel_fluence_DA_slab_TD(data.D, data.ν, t, data.ρ, data.μa, data.s, data.zb, data.z0, data.z, data.xs), data.t)
end

function _kernel_fluence_DA_slab_TD(D, ν, t, ρ, μa, s, zb, z0, z, xs)
    tmp1 = 4 * D * ν * t
    ϕ = ν * exp(-(ρ^2 / tmp1) - μa * ν * t)
    ϕ /= (π * tmp1)^(3/2)
    ϕ_tmp = zero(eltype(ϕ))

    for m in xs:-1:1
        ϕ_tmp += _images_fluence_DA_slab_TD(-m, s, zb, z0, z, tmp1)
        ϕ_tmp += _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    end
    ϕ_tmp += _images_fluence_DA_slab_TD(0, s, zb, z0, z, tmp1)
    ϕ *= ϕ_tmp
end

@inline function _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    tmp2 = 2 * m * (s + 2 * zb)
    zmp = tmp2 + z0
    zmm = tmp2 - 2 * zb - z0
    return exp(-((z - zmp)^2 / tmp1)) - exp(-((z - zmm)^2 / tmp1))
end

#--------------------------------------
# Time-Domain Flux
#--------------------------------------
"""
    flux_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 15)

    Compute the time-domain flux, D*∂ϕ(t)/∂z for z = 0 and -D*∂ϕ(t)/∂z for z = s from a slab geometry (x, y -> inf, z -> finite).
    If z != 0 or s will default to compute D*∂ϕ(t)/∂z for the z given. 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth coordinate (cm) must be equal
- `xs`: the number of sources to compute in the series

# Examples
julia> `flux_DA_slab_TD(1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0)`
"""
function flux_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0, xs = 15)
    if z == zero(eltype(z))
        return _refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s, xs)
    elseif z == s
        return _trans_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s, xs)
    else 
        return D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, z = dz, xs = xs), z)
    end
end

"""
    flux_DA_slab_TD(t, data::DiffusionParameters)

Wrapper to `flux_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0, xs = 15)`

# Examples
```
julia> data = DAslab(ρ = 1.0) # use structure to generate inputs
julia> flux_DA_slab_TD(data) # then call the function
```
"""
function flux_DA_slab_TD(data)
    return flux_DA_slab_TD(data.t, data.ρ, data.μa, data.μsp, n_ext = data.n_ext, n_med = data.n_med, z = data.z, s = data.s)
end

"""
    _refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs = -15:15)

Compute the time-domain flux or reflectance at the top surface from a slab geometry (x,y->inf, z-> finite). 

"""
function _refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s, xs)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t -> _kernel_refl_DA_slab_TD(ρ, params.D, params.ν, t, μa, xs, s, params.zb, params.z0), t)
end

@inline function _kernel_refl_DA_slab_TD(ρ, D, ν, t, μa, xs, s, zb, z0)
    Rt1 = -exp(-(ρ^2 / (4 * D * ν * t)) - μa * ν * t) / (2 * (4 * π * D * ν)^(3/2) * t^(5/2))
    Rt2 = zero(eltype(Rt1))
	for m in -xs:xs
	    z3m = -2 * m * s - 4 * m * zb - z0
	    z4m = -2 * m * s - (4m - 2) * zb + z0
	    Rt2 += z3m * exp(-(z3m^2 / (4 * D * ν * t))) - z4m * exp(-(z4m^2 / (4 * D * ν * t)))
    end	

	return Rt1 * Rt2  
end

"""
    trans_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs = -15:15)

Compute the time-domain transmittance (flux) from a slab geometry (x,y -> inf, z -> finite) with Eqn. 39 from Contini 97. 
"""
function _trans_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s, xs)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t -> _kernel_trans_DA_slab_TD(ρ, params.D, params.ν, t, μa, xs, s, params.zb, params.z0), t)
end

@inline function _kernel_trans_DA_slab_TD(ρ, D, ν, t, μa, xs, s, zb, z0)
    Rt1 = exp(-(ρ^2 / (4 * D * ν * t)) - μa * ν * t) / (2 * (4 * π * D * ν)^(3/2) * t^(5/2))
    Rt2 = zero(eltype(Rt1))
	for m in -xs:xs
	    z1m = (1 - 2 * m) * s - 4 * m * zb - z0
		z2m = (1 - 2 * m) * s - (4 * m - 2) * zb + z0
	    Rt2 += z1m * exp(-(z1m^2 / (4 * D * ν * t))) - z2m * exp(-(z2m^2 / (4 * D * ν * t)))
    end	

	return Rt1 * Rt2  
end

#--------------------------------------
# Steady-State Flux
#--------------------------------------
"""
    fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)

Compute the steady-state flux, D*∂ϕ(ρ)/∂z for z = 0 and -D*∂ϕ(ρ)/∂z for z = s from a slab geometry (x, y -> inf, z -> finite).
If z != 0 or s will default to compute D*∂ϕ(ρ)/∂z for the z given.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth within slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
```
julia> flux_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0)
```
"""
function flux_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    if z == zero(eltype(z))
        return params.D * ForwardDiff.derivative(dz -> fluence_DA_slab_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, s = s, z = dz, xs = xs), z)
    elseif z == s
        return -params.D * ForwardDiff.derivative(dz -> fluence_DA_slab_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, s = s, z = dz, xs = xs), z)
    else 
        return params.D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, z = dz, xs = xs), z)
    end
end

"""
    flux_DA_slab_CW(data::DiffusionParameters)

Wrapper to `flux_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)`

# Examples
```
julia> data = DAslab(ρ = 1.0) # use structure to generate inputs
julia> flux_DA_slab_CW(data) # then call the function
```
"""
function flux_DA_slab_CW(data)
    return flux_DA_slab_CW(data.ρ, data.μa, data.μsp, n_ext = data.n_ext, n_med = data.n_med, s = data.s, z = data.z, xs = data.xs)
end
