

##### Diffusion approximation for semi-infinite media  #####

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

The reflectance for the semi-infinite medium is calculated from Equation 36 from Contini 1997. 
The semi-infinite medium is described by retaining the first of the dipole sources (m=0) from Eqn. 36.
=#

#Reflectance semi-infinite

"""
    DA_semiinf(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)

Compute the time-domain reflectance from a semi-infinite medium. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction


# Examples
```jldoctest
julia> DA_semiinf(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0)
6-element Array{Float64,1}:
 0.0
 0.0001440103022493725
 1.446739954231315e-6
 2.7354735244571076e-8
 6.794070985483474e-10
 1.9657536202689858e-11
```
"""
function DA_semiinf(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)


    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed

    Rt1 = Array{Float64}(undef, length(t))
    Rt2 = Array{Float64}(undef, length(t))


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
		
   		Rt1[n] = -exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])
   		Rt1[n] = Rt1[n]/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

   		Rt2[n] = z3m*exp(-(z3m^2/(4D*ν*t[n]))) - z4m*exp(-(z4m^2/(4D*ν*t[n])))
		
	end
	
    

    replace!(Rt1.*Rt2, NaN => 0)

end



##### Diffusion approximation for Slab Geometry   #####
#The reflectance for the slab is calculated from Equation 36 from Contini 1997 for the first 11 dipoles (-10:1:10).
#The number of sources considered in the summation can be changed with xs. Larger values of m will be needed for large values of SDS and t.

"""
    DA_reflslab(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)

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
julia> DA_reflslab(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
6-element Array{Float64,1}:
 0.0
 0.00011885510267563147
 3.8254510325745536e-7
 1.5189689097582458e-9
 6.645341571314721e-12
 3.075337855389727e-14


```
"""
function DA_reflslab(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
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

		    Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))
			
				Threads.@threads for n in eachindex(t) 
					for m in xs

						z3m = -2m.*s .- 4m.*ze .- zs
						z4m = -2m.*s .- (4m .- 2).*ze .+ zs

			 			Rt2[n] += z3m*exp(-(z3m^2 / (4D*ν*t[n]))) - z4m*exp(-(z4m^2 / (4D*ν*t[n])))
					end
					
					
				end
			end

	return replace!(Rt1.*Rt2, NaN => 0)
end
##### Diffusion approximation for a turbid parallelepiped (reflectance)  #####



#=
@article{kienle2005light,
  title={Light diffusion through a turbid parallelepiped},
  author={Kienle, Alwin},
  journal={JOSA A},
  volume={22},
  number={9},
  pages={1883--1888},
  year={2005},
  publisher={Optical Society of America}
}

The fluence is given by equation 6 in Kienle 2005. The reflectance as calculated here uses Fick's law where Rt(rho, t) = D∂/∂zϕ(rho, z=0, t)
=#

"""
DA_reflparalpip(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})

Compute the time-domain reflectance from a parallelepiped [lx, ly, lz]. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `rd::Array{Float64,1}`: the location of the detector [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs::Array{Float64,1}`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L::Array{Float64,1}`: the dimenensions [lx, ly, lz] of the parallelepied



# Examples
```jldoctest
julia> DA_reflparalpip(0:1:5, [0.1,10.0], 1.0, 1.0, [2.5,5.0], [5.0,5.0], [10.0,10.0,10.0])
6-element Array{Float64,1}:
 0.0
 3.872453933271761e-5
 7.490933986486167e-7
 1.741880709155373e-8
 4.687789101530684e-10
 1.3825359487332909e-11


```
"""
function DA_reflparalpip(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})

	n::Float64 = nmed/ndet
	μa::Float64 = β[1]
	μsp::Float64 = β[2]
	D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x::Float64 = rd[1]
	y::Float64 = rd[2]
	xu::Float64 = rs[1]
	yu::Float64 = rs[2]
	lx::Float64 = L[1]
	ly::Float64 = L[2]
	lz::Float64 = L[3]


	Rt2 = zeros(Float64, length(t))
	Rt3 = zeros(Float64, length(t))
	Rt4 = zeros(Float64, length(t))

		if n == 1.0
			A= 1.0
		elseif n > 1.0
			A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
			 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
		else 
			A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
		end
	
		zo::Float64 = 1/μsp
		zb::Float64 = 2A*D

	  Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

  
	  Threads.@threads for n in eachindex(t) 
		  for m in xs
			    x1l = 2m.*lx .+ 4m.*zb .+ xu 
				x2l = 2m.*lx + (4m-2).*zb .- xu
				y1m = 2m.*ly .+ 4m.*zb .+ yu
				y2m = 2m.*ly .+ (4m.-2).*zb .- yu
				z1n = 2m.*lz .+ 4m.*zb .+ zo
				z2n = 2m.*lz .+ (4m.-2).*zb .- zo
		  
				Rt2[n] += exp(-(x-x1l)^2 / (4D*ν*t[n])) - exp(-(x-x2l)^2 / (4D*ν*t[n]))

				Rt3[n] += exp(-(y-y1m)^2 / (4D*ν*t[n])) - exp(-(y-y2m)^2 / (4D*ν*t[n]))

				Rt4[n] += z1n*exp(-z1n^2 /(4D*ν*t[n])) - z2n'*exp(-z2n^2 /(4D*ν*t[n]))

		  end
	  end
	return replace!(Rt1.*Rt2.*Rt3.*Rt4, NaN => 0)
end

