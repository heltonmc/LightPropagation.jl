##### Diffusion approximation for Slab Geometry #####

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
@inproceedings{martelli2009light,
  title={Light propagation through biological tissue and other diffusive media: theory, solutions, and software},
  author={Martelli, Fabrizio and Del Bianco, Samuele and Ismaelli, Andrea},
  year={2009},
  organization={Society of Photo-Optical Instrumentation Engineers}
}
=#

### Fluence TD ###
"""
    fluence_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64, s::Float64, z::Float64)

Compute the time-domain fluence from a slab geometry (x, y->inf, z-> finite) with Eqn. 33 Contini 97. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)
- `z::Float64`: the z-depth coordinate (cm)

# Examples
julia> fluence_DA_slab_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0, 0.0)
"""
function fluence_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64, s::Float64, z::Float64)
	  n::Float64 = nmed/ndet
      μa::Float64 = β[1]
      μsp::Float64 = β[2]
      D::Float64 = 1/3μsp
      ν::Float64 = 29.9792345/nmed
      xs::UnitRange{Int64} = -10:10

      Rt1 = Array{Float64}(undef, length(t))
	  Rt2 = zeros(Float64, length(t))

		@inbounds begin

			if n == 1.0
				A= 1.0
			elseif n > 1.0
				A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 		5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
			else
				A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
			end

			zs::Float64 = 1/μsp
			ze::Float64 = 2A*D

                Threads.@threads for n in eachindex(t) 
                    Rt1[n] = @. ν*exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])/((4π*D*ν*t[n])^(3/2))

					for m in xs

						zmp = 2m.*(s + 2*ze) .+ zs
						zmm = 2m.*(s + 2*ze) .- 2*ze .- zs

			 			Rt2[n] += exp(-((z - zmp)^2 / (4D*ν*t[n]))) - exp(-((z - zmm)^2 / (4D*ν*t[n])))
                    end	

                    Rt1[n] *= Rt2[n]

                    if isnan(Rt1[n])
                        Rt1[n] = 0
                    end
             
				end
		end

	return Rt1
end

### Reflectance TD ###
"""
    refl_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain reflectance from a slab geometry (x,y->inf, z-> finite). 

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
function refl_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	  n::Float64 = nmed/ndet
      μa::Float64 = β[1]
      μsp::Float64 = β[2]
      D::Float64 = 1/3μsp
      ν::Float64 = 29.9792345/nmed
      xs::UnitRange{Int64} = -10:10

      Rt1 = Array{Float64}(undef, length(t))
	  Rt2 = zeros(Float64, length(t))

		@inbounds begin

			if n == 1.0
				A= 1.0
			elseif n > 1.0
				A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 		5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
			else
				A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
			end


			zs::Float64 = 1/μsp
			ze::Float64 = 2A*D

			
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
function trans_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	  n::Float64 = nmed/ndet
      μa::Float64 = β[1]
      μsp::Float64 = β[2]
      D::Float64 = 1/3μsp
      ν::Float64 = 29.9792345/nmed
      xs::UnitRange{Int64} = -10:10

      Rt1 = Array{Float64}(undef, length(t))
	  Rt2 = zeros(Float64, length(t))

		@inbounds begin

			if n == 1.0
				A= 1.0
			elseif n > 1.0
				A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 		5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
			else
				A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
			end

			zs::Float64 = 1/μsp
			ze::Float64 = 2A*D
	
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

### Fluence CW ###
"""
    fluence_DA_slab_CW(ρ::Float64, β::Array{Float64,1}, ndet::Float64, nmed::Float64, s::Float64, z::Float64)

Compute the steady-state fluence from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)
- `z::Float64`: the z-depth within slab (cm)

# Examples
julia> fluence_DA_slab_CW(1.0, [0.1,10.], 1.0,1.0, 2.0, 0.)
"""
function fluence_DA_slab_CW(ρ::Float64, β::Array{Float64,1}, ndet::Float64, nmed::Float64, s::Float64, z::Float64)
      
      n::Float64 = nmed/ndet
      μa::Float64 = β[1]
      μsp::Float64 = β[2]
      D::Float64 = 1/3μsp
      μeff::Float64 = sqrt(3*μa*μsp)
      xs::UnitRange{Int64} = -10:10

	  ϕ = 0.0

		if n == 1.0
			A= 1.0
		elseif n > 1.0
			A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
			5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
		else
			A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
		end

		zs::Float64 = 1/μsp
		ze::Float64 = 2A*D

        Threads.@threads for m in xs 

			zmp = 2m.*(s + 2*ze) .+ zs
			zmm = 2m.*(s + 2*ze) .- 2*ze .- zs

            ϕ += exp(-μeff*sqrt(ρ^2 + (z - zmp)^2))/(sqrt(ρ^2 + (z - zmp)^2))
            ϕ -= exp(-μeff*sqrt(ρ^2 + (z - zmm)^2))/(sqrt(ρ^2 + (z - zmm)^2))
        end	

        if isnan(ϕ)
            ϕ = 0
        end
             
	return ϕ/(4*π*D)
end