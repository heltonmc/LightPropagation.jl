##### Diffusion approximation for semi-infinite media #####

#=
@article{contini1997photon,
  title={Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory},
  author={Contini, Daniele and Martelli, Fabrizio and Zaccanti, Giovanni},
  journal={Applied optics},
  volume={36},
  number={19},
  pages={4587--4599},
  year={1997},
  publisher={Optical Society of America}
}

@article{kienle1997improved,
  title={Improved solutions of the steady-state and the time-resolved diffusion equations for reflectance from a semi-infinite turbid medium},
  author={Kienle, Alwin and Patterson, Michael S},
  journal={JOSA A},
  volume={14},
  number={1},
  pages={246--254},
  year={1997},
  publisher={Optical Society of America}
}
=#

### TD Reflectance ###
"""
    refl_DA_semiinf_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)

Compute the time-domain reflectance from a semi-infinite medium from Eqn. 36 Contini 97 (m = 0). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction

# Examples
```jldoctest
julia> refl_DA_semiinf_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0)
6-element Array{Float64,1}:
 0.0
 0.0001440103022493725
 1.446739954231315e-6
 2.7354735244571076e-8
 6.794070985483474e-10
 1.9657536202689858e-11
```
"""
function refl_DA_semiinf_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)
    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
  	ν::Float64 = 29.9792345/nmed

    Rt = Array{Float64}(undef, length(t))

	if n > 1.0
		  A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		  5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	elseif n < 1.0
		  A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	else 
		  A = 1.0
	end

	zs::Float64 = 1/μsp
	ze::Float64 = 2A*D

	z3m::Float64 = - zs
    z4m::Float64 = 2ze +zs

    Threads.@threads for n in eachindex(t)

        Rt[n] = -(ρ^2/(4D*ν*t[n]))
        Rt[n] = Rt[n] - μa*ν*t[n]
        Rt[n] = -exp(Rt[n])/(2*(4π*D*ν)^(3/2)*sqrt(t[n]*t[n]*t[n]*t[n]*t[n]))
        Rt[n] = Rt[n]*(z3m*exp(-(z3m^2/(4D*ν*t[n]))) - z4m*exp(-(z4m^2/(4D*ν*t[n]))))

        if isnan(Rt[n])
            Rt[n] = 0
        end
     
    end
 
     return Rt
end

### TD Fluence ###
"""
    fluence_DA_semiinf_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64, z::Float64)

Compute the time-domain fluence in a semi-infinite medium (Eqn. 33 Contini). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `z::Float64`: the z-depth in medium

# Examples
julia> fluence_DA_semiinf_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 1.0)
"""
function fluence_DA_semiinf_TD(t, β, ρ, ndet, nmed, z)
    n = nmed/ndet
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
  	ν = 29.9792345/nmed

    ϕ = similar(t)

	if n > 1.0
		  A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		  5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	elseif n < 1.0
		  A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	else 
		  A = 1.0
	end

	z0 = 1/μsp
	zb = 2A*D

    Threads.@threads for ind in eachindex(t)
        ϕ[ind] = ν*exp(-μa*ν*t[ind])/(4*π*D*ν*t[ind])^(3/2)
        ϕ[ind] *= (exp(-((z - z0)^2 + ρ^2)/(4*D*ν*t[ind])) - exp(-((z + z0 + 2*zb)^2 + ρ^2)/(4*D*ν*t[ind])))

        if isnan(ϕ[ind])
            ϕ[ind] = 0
        end
     
    end
 
     return ϕ
end

### CW Fluence ###
"""
    fluence_DA_semiinf_CW(ρ::Float64, β::Array{Float64,1}, ndet::Float64, nmed::Float64, z::Float64)

Compute the steady-state fluence in a semi-infinite geometry according to Eqn. 3 of Kienle 1997. 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `z::Float64`: the z-depth within slab (cm)

# Examples
julia> fluence_DA_semiinf_CW(1.0, [0.1,10.0], 1.0,1.0, 0.0)
"""
function fluence_DA_semiinf_CW(ρ::Float64, β::Array{Float64,1}, ndet::Float64, nmed::Float64, z::Float64)
      n = nmed/ndet
      μa = β[1]
      μsp = β[2]
      D = 1/3μsp
      μeff = sqrt(3*μa*μsp)

	  ϕ = 0.0

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

        ϕ += exp(-μeff*sqrt(ρ^2 + (z - zs)^2))/(sqrt(ρ^2 + (z - zs)^2))
        ϕ -= exp(-μeff*sqrt(ρ^2 + (z + 2*ze + zs)^2))/(sqrt(ρ^2 + (z + 2*ze + zs)^2))

        if isnan(ϕ)
            ϕ = 0
        end
             
	return ϕ/(4*π*D)
end


function fluence_DA_semiinf_FD(ρ::Float64, β::Array{Float64,1}, ndet::Float64, nmed::Float64, z::Float64, ω)
    n = nmed/ndet
    ν = 29.9792345/nmed
    μa = β[1] + ω*im/ν
    μsp = β[2]
    D = 1/3μsp
    μeff = sqrt(3*μa*μsp)

    ϕ = 0.0

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

      ϕ += exp(-μeff*sqrt(ρ^2 + (z - zs)^2))/(sqrt(ρ^2 + (z - zs)^2))
      ϕ -= exp(-μeff*sqrt(ρ^2 + (z + 2*ze + zs)^2))/(sqrt(ρ^2 + (z + 2*ze + zs)^2))

      if isnan(ϕ)
          ϕ = 0
      end
           
  return ϕ/(4*π*D)
end