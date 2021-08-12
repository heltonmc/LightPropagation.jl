#########################################################################################################################################
# Implements solution to the diffusion equation for the fluence in the steady-state and time-domain in a semi-infinite medium [1][2][3].
#
# [1] Alwin Kienle and Michael S. Patterson, "Improved solutions of the steady-state and the time-resolved diffusion equations 
#     for reflectance from a semi-infinite turbid medium," J. Opt. Soc. Am. A 14, 246-254 (1997) 
# [2] Martelli, Fabrizio, et al. Light propagation through biological tissue and other diffusive media: theory, solutions and software. 
#     Vol. 10. No. 3.824746. Bellingham: SPIE press, 2010.
#########################################################################################################################################

#####################################
# Steady-State Fluence 
#####################################
"""
    fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)

Compute the steady-state fluence in a semi-infinite geometry according to Eqn. 3 of Kienle 1997. 

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)

# Examples
julia> `fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0)`
"""
function fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    D = D_coeff(μsp, μa)
    μeff = sqrt(3 * μa * μsp)
    A = A_coeff(n_med / n_ext)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    μeff = sqrt(3 * μa * μsp)

    return _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
end
@inline function _kernel_fluence_DA_semiinf_CW(μeff, ρ, z, z0, zb, D)
    ϕ = exp(-μeff * sqrt(ρ^2 + (z - z0)^2)) / (sqrt(ρ^2 + (z - z0)^2))
    ϕ -= exp(-μeff * sqrt(ρ^2 + (z + 2 * zb + z0)^2)) / (sqrt(ρ^2 + (z + 2 * zb + z0)^2))

    return ϕ / (4 * π * D)
end

#####################################
# Time-Domain Fluence 
#####################################
"""
    fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)

Compute the time-domain fluence in a semi-infinite medium (Eqn. 33 Contini). 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth in medium

# Examples
julia> `fluence_DA_semiinf_TD(0.1:0.1:1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0)`
"""
function fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)
    D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    if isa(t, AbstractFloat)
		return _kernel_fluence_DA_semiinf_TD(D, ν, t, μa, z, z0, ρ, zb)
	elseif isa(t, AbstractVector)
        T = promote_type(eltype(ρ), eltype(μa), eltype(μsp), eltype(z))
        ϕ = zeros(T, length(t))

        @inbounds Threads.@threads for ind in eachindex(t)
    		ϕ[ind] = _kernel_fluence_DA_semiinf_TD(D, ν, t[ind], μa, z, z0, ρ, zb)
    	end
    	return ϕ
	end
end
@inline function _kernel_fluence_DA_semiinf_TD(D, ν, t, μa, z, z0, ρ, zb)
    @assert t > zero(eltype(D))
    tmp1 = 4 * D * ν * t
    ϕ = ν * exp(-μa * ν * t) / (π * tmp1)^(3/2)
    ϕ *= (exp(-((z - z0)^2 + ρ^2) / tmp1) - exp(-((z + z0 + 2*zb)^2 + ρ^2) / tmp1))

    return ϕ
end
#####################################
# Frequency-Domain Fluence 
#####################################
"""
fluence_DA_semiinf_FD(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)

Compute the frequency domain fluence in a semi-infinite geometry.

# Arguments
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)
- `ω`: the modulation frequency (1/ns)

# Examples
julia> `fluence_DA_semiinf_FD(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)`
"""
function fluence_DA_semiinf_FD(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0, ω = 1.0)
    ν = ν_coeff(n_med)
    μa_complex = μa + ω * im / ν

    return fluence_DA_semiinf_CW(ρ, μa_complex, μsp, n_ext = n_ext, n_med = n_med, z = z)
end

#####################################
# Time-Domain Flux
#####################################
"""
    flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)

Compute the time-domain flux (D*∂ϕ(t)/∂z @ z = 0) from a semi-infinite medium from Eqn. 36 Contini 97 (m = 0).

# Arguments
- `t`: the time vector (ns). 
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_det`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction

# Examples
```jldoctest
julia> `flux_DA_semiinf_TD(1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)`
```
"""
function flux_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    z3m = - z0
    z4m = 2 * zb + z0

    if isa(t, AbstractFloat)
        return _kernel_flux_DA_semiinf_TD(D, ν, t, μa, z3m, z4m, ρ)
    elseif isa(t, AbstractArray)
        T = promote_type(eltype(ρ), eltype(μa), eltype(μsp))
        Rt = zeros(T, length(t))
        @inbounds Threads.@threads for ind in eachindex(t)
            Rt[ind] = _kernel_flux_DA_semiinf_TD(D, ν, t[ind], μa, z3m, z4m, ρ)
        end
        return Rt
    end
 end
@inline function _kernel_flux_DA_semiinf_TD(D, ν, t, μa, z3m, z4m, ρ)
    @assert t > zero(eltype(D))
    tmp1 = 4 * D * ν * t
    Rt = -(ρ^2 / (tmp1))
    Rt -= μa * ν * t
    Rt = -exp(Rt) / (2 * (4 * π * D * ν)^(3/2) * sqrt(t * t * t * t * t))
    Rt *= (z3m * exp(-(z3m^2 / (tmp1))) - z4m * exp(-(z4m^2 / (tmp1))))

    return Rt
end

#####################################
# Steady-State Flux
#####################################
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
```jldoctest
julia> `flux_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)`
```
"""
function flux_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
    D = D_coeff(μsp, μa)
    return D * ForwardDiff.derivative(z -> fluence_DA_semiinf_CW(ρ, μa, μsp, n_med = n_med, n_ext = n_ext, z = z), 0.0)
end
