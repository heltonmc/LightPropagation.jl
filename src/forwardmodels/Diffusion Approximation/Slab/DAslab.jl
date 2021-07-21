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
    fluence_DA_slab_CW(ρ, μa, μsp; n_ext, n_med, s, z, xs)

Compute the steady-state fluence from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext::Float64`: the boundary's index of refraction (air or detector)
- `n_med::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)
- `z::Float64`: the z-depth within slab (cm)
- `xs`: the number of sources to compute in the series

# Examples
julia> `fluence_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0)`
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
    fluence_DA_slab_TD(t, ρ, μa, μsp; n_ext, n_med, s, z, xs)

Compute the time-domain fluence from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)
- `z::Float64`: the z-depth coordinate (cm)
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
    	ϕ = zeros(eltype(D), length(t))
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
# Time-Domain Reflectance (flux) 
#####################################
"""
    refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs)

Compute the time-domain reflectance (flux) from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ρ`: the source detector separation (cm⁻¹)
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)

# Examples
julia> `refl_DA_slab_TD(1.0, 1.0, 0.1, 10.0, 1.0, 1.0, 1.0)`
"""
function refl_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs = -15:15)
	D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)
    
    Rt1 = Array{eltype(ρ)}(undef, length(t))
    Rt2 = zeros(eltype(ρ), length(t))
	
    Threads.@threads for n in eachindex(t) 
        Rt1[n] = -exp(-(ρ^2 / (4 * D * ν * t[n])) - μa * ν * t[n]) / (2 * (4 * π * D * ν)^(3/2) * t[n]^(5/2))
	    for m in xs
	        z3m = -2 * m * s - 4 * m * zb - z0
	        z4m = -2 * m * s - (4m - 2) * zb + z0
	    	Rt2[n] += z3m * exp(-(z3m^2 / (4 * D * ν * t[n]))) - z4m * exp(-(z4m^2 / (4 * D * ν * t[n])))
    	end		
	    Rt1[n] *= Rt2[n]   
	end

	return Rt1
end

#####################################
# Time-Domain Transmittance (flux) 
#####################################
"""
    trans_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs)

Compute the time-domain transmittance (flux) from a slab geometry (x,y -> inf, z -> finite) with Eqn. 39 from Contini 97. 

# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ρ`: the source detector separation (cm⁻¹)
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction
- `s`: the thickness (z-depth) of the slab (cm)

# Examples
julia> `trans_DA_slab_TD(1.0, 1.0, 0.1, 10.0, 1.0, 1.0, 1.0)`
"""
function trans_DA_slab_TD(t, ρ, μa, μsp, n_ext, n_med, s; xs = -15:15)
	D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    Rt1 = Array{eltype(ρ)}(undef, length(t))
    Rt2 = zeros(eltype(ρ), length(t))

	Threads.@threads for n in eachindex(t) 
        Rt1[n] = exp(-(ρ^2 / (4 * D * ν * t[n])) - μa * ν * t[n]) / (2 * (4 * π * D * ν)^(3/2) * t[n]^(5/2))
        for m in xs
            z1m = (1 - 2 * m) * s - 4 * m * zb - z0
		    z2m = (1 - 2 * m) * s - (4 * m - 2) * zb + z0

		    Rt2[n] += z1m * exp(-(z1m^2 / (4 * D * ν * t[n]))) - z2m * exp(-(z2m^2 / (4 * D * ν * t[n])))
	    end
	    Rt1[n] *= Rt2[n]
    end
			
	return Rt1
end
