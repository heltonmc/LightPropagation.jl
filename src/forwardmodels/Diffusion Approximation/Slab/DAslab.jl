################################################################################################################################## 
# Implements solution to the diffusion equation for the fluence in the steady-state and time-domain in a slab [1][2].
#
# [1] Contini, Daniele, Fabrizio Martelli, and Giovanni Zaccanti. "Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory." 
#     Applied optics 36.19 (1997): 4587-4599.
# [2] Martelli, Fabrizio, et al. Light propagation through biological tissue and other diffusive media: theory, solutions and software. Vol. 10. No. 3.824746. Bellingham: SPIE press, 2010.
################################################################################################################################## 

# Steady-State Fluence
@inline function _kernel_fluence_DA_slab_CW(m, s, zb, z0, ρ, z, μeff)
    tmp1 = 2 * m * (s + 2 * zb)
    zmp = tmp1 + z0
    zmm = tmp1 - 2 * zb - z0

    tmp2 = sqrt(ρ^2 + (z - zmp)^2)
    tmp3 = sqrt(ρ^2 + (z - zmm)^2)
    ϕ = exp(-μeff * tmp2) / tmp2
    ϕ -= exp(-μeff * tmp3) / tmp3
end
function fluence_DA_slab_CW(ρ, μa, μsp, n_det, n_med, s, z; rtol = 1e-20, maxiter = 1e5)
    D = D_coeff(μsp, μa)
    μeff = sqrt(3 * μa * μsp)
    A = get_afac(n_med / n_det)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    ϕ = _kernel_fluence_DA_slab_CW(0, s, zb, z0, ρ, z, μeff)
    ϕ_new = zero(eltype(ϕ))

    for m in 1:maxiter
	    tmp = _kernel_fluence_DA_slab_CW(m, s, zb, z0, ρ, z, μeff)
	    tmp += _kernel_fluence_DA_slab_CW(-m, s, zb, z0, ρ, z, μeff)
	    ϕ_new += tmp

	    if abs(tmp) / (abs(ϕ_new))  < rtol
            break
	    end
    end
             
    return (ϕ + ϕ_new) / (4 * π * D)
end
function fluence_DA_slab_CW(data)
    return fluence_DA_slab_CW(data.ρ, data.μa, data.μsp, data.n_det, data.n_med, data.s, data.z)
end
@doc """
    fluence_DA_slab_CW(ρ, μa, μsp, n_det, n_med, s, z; rtol, maxiter)

Compute the steady-state fluence from a slab geometry (x, y -> inf, z -> finite). 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_det::Float64`: the boundary's index of refraction (air or detector)
- `n_med::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)
- `z::Float64`: the z-depth within slab (cm)

# Examples
julia> fluence_DA_slab_CW(1.0, [0.1,10.], 1.0,1.0, 2.0, 0.)
""" fluence_DA_slab_CW

### Fluence TD ###
@inline function _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    tmp2 = 2 * m * (s + 2 * zb)
    zmp = tmp2 + z0
    zmm = tmp2 - 2 * zb - z0
    return exp(-((z - zmp)^2 / tmp1)) - exp(-((z - zmm)^2 / tmp1))
end
@inline function _kernel_fluence_DA_slab_TD(D, ν, t, ρ, μa, s, zb, z0, z, maxiter, rtol)
    tmp1 = 4 * D * ν * t
    ϕ = ν * exp(-(ρ^2 / tmp1) - μa * ν * t)
    ϕ /= (π * tmp1)^(3/2)
    ϕ_tmp = _images_fluence_DA_slab_TD(0, s, zb, z0, z, tmp1)

    for m in 1:maxiter
    	tmp = _images_fluence_DA_slab_TD(m, s, zb, z0, z, tmp1)
    	tmp += _images_fluence_DA_slab_TD(-m, s, zb, z0, z, tmp1)
    	ϕ_tmp += tmp
			
    	if abs(tmp) / (abs(ϕ_tmp))  < rtol
    		break
    	end
    end
    ϕ *= ϕ_tmp
end
function fluence_DA_slab_TD(t::AbstractFloat, ρ, μa, μsp, n_det, n_med, s, z; rtol = 1e-20, maxiter = 1e5)
    D = D_coeff(μsp, μa)
    A = get_afac(n_med / n_det)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    return _kernel_fluence_DA_slab_TD(D, ν, t, ρ, μa, s, zb, z0, z, maxiter, rtol)
end
function fluence_DA_slab_TD(t::AbstractArray, ρ, μa, μsp, n_det, n_med, s, z; rtol = 1e-20, maxiter = 1e5)
    D = D_coeff(μsp, μa)
    A = get_afac(n_med / n_det)
    ν = ν_coeff(n_med)

    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    ϕ = zeros(eltype(ρ), length(t))
    for ind in eachindex(t)
    	ϕ[ind] = _kernel_fluence_DA_slab_TD(D, ν, t[ind], ρ, μa, s, zb, z0, z, maxiter, rtol)
    end
    return ϕ
end
function fluence_DA_slab_TD(data)
    return fluence_DA_slab_TD(data.t, data.ρ, data.μa, data.μsp, data.n_det, data.n_med, data.s, data.z)
end
@doc """
    fluence_DA_slab_TD(t, ρ, μa, μsp, n_det, n_med, s, z; rtol, maxiter)

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

# Examples
julia> fluence_DA_slab_TD(0.1:0.1:1.0, 1.0, 0.1, 10.0, 1.0, 1.0, 40.0, 0.0)
""" fluence_DA_slab_TD


### Reflectance TD ###
"""
    refl_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain reflectance (flux) from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)

# Examples
julia> refl_DA_slab_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
"""
function refl_DA_slab_TD(t, β, ρ,ndet, nmed, s)
	n = nmed/ndet
	μa = β[1]
    μsp = β[2]
    D = 1/3μsp
    ν = 29.9792458/nmed
    xs = -10:10

    Rt1 = Array{eltype(ρ)}(undef, length(t))
	Rt2 = zeros(eltype(ρ), length(t))

	@inbounds begin

		if n == 1.0
			A= 1.0
		elseif n > 1.0
			A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
	 		5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
		else
			A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
		end


		zs = 1/μsp
		ze = 2A*D

			
        Threads.@threads for n in eachindex(t) 
            Rt1[n] = @. -exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

			for m in xs

				z3m = -2m*s - 4m*ze - zs
				z4m = -2m*s - (4m - 2)*ze + zs

	 			Rt2[n] += z3m*exp(-(z3m^2 / (4D*ν*t[n]))) - z4m*exp(-(z4m^2 / (4D*ν*t[n])))
			end
					
			Rt1[n] *= Rt2[n]

            if isnan(Rt1[n])
                Rt1[n] = 0
            end
		    end
		end

	return Rt1
end

### Transmittance TD ###
"""
    trans_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain transmittance from a slab geometry (x,y->inf, z-> finite) with Eqn. 39 from Contini 97. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)

# Examples
julia> trans_DA_slab_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
"""
function trans_DA_slab_TD(t, β, ρ, ndet, nmed, s)
	  n = nmed/ndet
      μa = β[1]
      μsp = β[2]
      D = 1/3μsp
      ν = 29.9792458/nmed
      xs = -10:10

      Rt1 = Array{eltype(ρ)}(undef, length(t))
	  Rt2 = zeros(eltype(ρ), length(t))

		@inbounds begin

			if n == 1.0
				A= 1.0
			elseif n > 1.0
				A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 		5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
			else
				A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
			end

			zs = 1/μsp
			ze = 2A*D
	
                Threads.@threads for n in eachindex(t) 
                    Rt1[n] = @. exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

					for m in xs

						z1m = (1 - 2*m)*s - 4m*ze - zs
						z2m = (1 - 2*m)*s - (4m - 2)*ze + zs

			 			Rt2[n] += z1m*exp(-(z1m^2 / (4D*ν*t[n]))) - z2m*exp(-(z2m^2 / (4D*ν*t[n])))
					end
					
					Rt1[n] *= Rt2[n]

                    if isnan(Rt1[n])
                        Rt1[n] = 0
                    end
				end
			end

	return Rt1
end


#= old version of fluence slab that is left for clarity
function fluence_DA_slab_CW(ρ, μa, μsp, n_det, n_med, s, z; xs = -10:10)
    D = D_coeff(μsp, μa)
	μeff = sqrt(3 * μa * μsp)
	A = get_afac(n_med / n_det)
	
	ϕ = zero(eltype(ρ))

	z0 = z0_coeff(μsp)
	zb = zb_coeff(A, D)

	for m in xs
		tmp1 = 2 * m * (s + 2 * zb)
        zmp = tmp1 + z0
		zmm = tmp1 - 2 * zb - z0

		tmp2 = sqrt(ρ^2 + (z - zmp)^2)
		tmp3 = sqrt(ρ^2 + (z - zmm)^2)
		ϕ += exp(-μeff * tmp2) / tmp2
		ϕ -= exp(-μeff * tmp3) / tmp3
	end
             
	return ϕ / (4 * π * D)
end

function fluence_DA_slab_TD(t, ρ, μa, μsp, n_det, n_med, s, z; xs = -10:10)
    D = D_coeff(μsp, μa)
	A = get_afac(n_med / n_det)
	ν = ν_coeff(n_med)

	z0 = z0_coeff(μsp)
	zb = zb_coeff(A, D)

	Rt1 = Vector{eltype(ρ)}(undef, length(t))
	Rt2 = zeros(eltype(Rt1), length(t))

	for ind in eachindex(t)
		tmp1 = 4 * D * ν * t[ind]
		Rt1[ind] = ν * exp(-(ρ^2 / tmp1) - μa * ν * t[ind])
		Rt1[ind] /= (π * tmp1)^(3/2)
		for m in xs
			tmp2 = 2 * m * (s + 2 * zb)
		    zmp = tmp2 + z0
		    zmm = tmp2 - 2*zb - z0

		    Rt2[ind] += exp(-((z - zmp)^2 / tmp1)) - exp(-((z - zmm)^2 / tmp1))
		end	
		Rt1[ind] *= Rt2[ind]
	end
	return Rt1
end
=#
