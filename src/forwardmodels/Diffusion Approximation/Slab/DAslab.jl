##### Diffusion approximation for Slab Geometry (TD-fluence) #####

#=
The reflectance for the slab is calculated from Equation 36 from Contini 1997 for the first 11 dipoles (-10:1:10).
The number of sources considered in the summation can be changed with xs. Larger values of m will be needed for large values of SDS and t.
=#

"""
    fluence_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain reflectance from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)



# Examples

julia> fluence_DA_slab_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)

"""
function fluence_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64, z::Float64)
	#s is slab thickness

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

##### Diffusion approximation for Slab Geometry (TD-reflectance) #####

#=
The reflectance for the slab is calculated from Equation 36 from Contini 1997 for the first 11 dipoles (-10:1:10).
The number of sources considered in the summation can be changed with xs. Larger values of m will be needed for large values of SDS and t.
=#

"""
    TPSF_DA_slab_refl(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain reflectance from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)



# Examples
```jldoctest
julia> TPSF_DA_slab_refl(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
6-element Array{Float64,1}:
 0.0
 0.00011885510267563147
 3.8254510325745536e-7
 1.5189689097582458e-9
 6.645341571314721e-12
 3.075337855389727e-14


```
"""
function TPSF_DA_slab_refl(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness

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


##### Diffusion approximation for Slab Geometry (TD-transmittance) #####

#=
The transmittance for the slab is calculated from Equation 39 from Contini 1997 for the first 11 dipoles (-10:1:10).
The number of sources considered in the summation can be changed with xs. Larger values of m will be needed for large values of SDS and t.
=#

"""
    trans_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

Compute the time-domain reflectance from a slab geometry (x,y->inf, z-> finite). 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `s::Float64`: the thickness (z-depth) of the slab (cm)



# Examples
```jldoctest
julia> trans_DA_slab_TD(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
6-element Array{Float64,1}:
 0.0
 0.00011885510267563147
 3.8254510325745536e-7
 1.5189689097582458e-9
 6.645341571314721e-12
 3.075337855389727e-14


```
"""
function trans_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness

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


##### Diffusion approximation for Slab Geometry (CW-fluence) #####

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