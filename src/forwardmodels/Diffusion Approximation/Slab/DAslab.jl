################################################################################################################################## 
# Implements solution to the diffusion equation for the fluence in the steady-state and time-domain in a slab [1][2].
#
# [1] Contini et. al. "Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory." 
#     Applied optics 36.19 (1997): 4587-4599.
# [2] Martelli et al. Light propagation through biological tissue and other diffusive media: theory, solutions and software. 
#     Vol. 10. No. 3.824746. Bellingham: SPIE press, 2010.
################################################################################################################################## 

#####################################
# Steady-State Fluence 
#####################################
"""
    fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)

Compute the steady-state fluence from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth within slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
julia> `fluence_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0)`
"""
function fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)
    D = D_coeff(μsp, μa)
    μeff = sqrt(3 * μa * μsp)
    A = A_coeff(n_med / n_ext)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    ϕ = zero(eltype(μeff))

    for m in -xs:xs
        ϕ += _kernel_fluence_DA_slab_CW(m, s, zb, z0, ρ, z, μeff)
    end

    return (ϕ) / (4 * π * D)
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

#####################################
# Time-Domain Fluence 
#####################################
"""
    fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)

Compute the time-domain fluence from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth coordinate (cm)
- `xs`: the number of sources to compute in the series

# Examples
julia> `fluence_DA_slab_TD(0.1, 1.0, 0.1, 10.0, s = 40.0)`
"""
function fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 10)
    D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

	if isa(t, AbstractFloat)
    	return _kernel_fluence_DA_slab_TD(D, ν, t, ρ, μa, s, zb, z0, z, xs)
    elseif isa(t, AbstractArray)
    	T = promote_type(eltype(ρ), eltype(D), eltype(s), eltype(z))
        ϕ = zeros(T, length(t))

        @inbounds Threads.@threads for ind in eachindex(t)
            ϕ[ind] = _kernel_fluence_DA_slab_TD(D, ν, t[ind], ρ, μa, s, zb, z0, z, xs)
        end
    	return ϕ
    end
end
@inline function _kernel_fluence_DA_slab_TD(D, ν, t, ρ, μa, s, zb, z0, z, xs)
	@assert t > zero(eltype(D))
    tmp1 = 4 * D * ν * t
    ϕ = ν * exp(-(ρ^2 / tmp1) - μa * ν * t)
    ϕ /= (π * tmp1)^(3/2)
    ϕ_tmp = zero(eltype(ϕ))

    for m in -xs:xs
    	ϕ_tmp += _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    end
    ϕ *= ϕ_tmp
end
@inline function _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    tmp2 = 2 * m * (s + 2 * zb)
    zmp = tmp2 + z0
    zmm = tmp2 - 2 * zb - z0
    return exp(-((z - zmp)^2 / tmp1)) - exp(-((z - zmm)^2 / tmp1))
end

#####################################
# Time-Domain Flux
#####################################
"""
    flux_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0, xs = 15)

    Compute the time-domain flux (D*∂ϕ(t)/∂z @ z = 0 or z = s) from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth coordinate (cm)
- `xs`: the number of sources to compute in the series

# Examples
julia> `flux_DA_slab_TD(1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0)`
"""
function flux_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, s = 1.0, xs = 15)
    if z == zero(eltype(z))
        return _refl_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, xs = xs)
    elseif z == s
        return _trans_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, xs = xs)
    else 
        return D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, z = dz, xs = xs), z)
    end
end

"""
    _refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs = -15:15)

Compute the time-domain flux or reflectance at the top surface from a slab geometry (x,y->inf, z-> finite). 

"""
function _refl_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, xs = 15)
	D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)
    
    if isa(t, AbstractFloat)
    	return _kernel_refl_DA_slab_TD(ρ, D, ν, t, μa, xs, s, zb, z0)
    elseif isa(t, AbstractVector)
    	T = promote_type(eltype(ρ), eltype(D), eltype(s))
        Rt = zeros(T, length(t))

        @inbounds Threads.@threads for ind in eachindex(t)
            Rt[ind] = _kernel_refl_DA_slab_TD(ρ, D, ν, t[ind], μa, xs, s, zb, z0)
        end

    	return Rt
    end
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
function _trans_DA_slab_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0, xs = 15)
	D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)
    
    if isa(t, AbstractFloat)
    	return _kernel_trans_DA_slab_TD(ρ, D, ν, t, μa, xs, s, zb, z0)
    elseif isa(t, AbstractVector)
    	T = promote_type(eltype(ρ), eltype(D), eltype(s))
        Rt = zeros(T, length(t))

        @inbounds Threads.@threads for ind in eachindex(t)
            Rt[ind] = _kernel_trans_DA_slab_TD(ρ, D, ν, t[ind], μa, xs, s, zb, z0)
        end

    	return Rt
    end
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

#####################################
# Steady-State Flux
#####################################
"""
fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)

    flux_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)

Compute the steady-state flux (D*∂ϕ(ρ)/∂z for z = 0 or s from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)
- `z`: the z-depth within slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
```jldoctest
julia> `flux_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0)`
```
"""
function flux_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0, xs = 10)
    D = D_coeff(μsp, μa)
    if z == zero(eltype(z))
        return D * ForwardDiff.derivative(dz -> fluence_DA_slab_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, s = s, z = dz, xs = xs), z)
    elseif z == s
        return -D * ForwardDiff.derivative(dz -> fluence_DA_slab_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, s = s, z = dz, xs = xs), z)
    else 
        return D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, s = s, z = dz, xs = xs), z)
    end
end