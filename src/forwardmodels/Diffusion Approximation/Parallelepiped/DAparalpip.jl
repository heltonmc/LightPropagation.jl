##### Diffusion approximation for a turbid parallelepiped #####

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
=#

### TD Fluence ###
"""
    fluence_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})

Compute the time-domain fluence in a parallelepiped [lx, ly, lz] with Eqn. 6 Kienle 05. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `ndet::Float64`: the boundary's index of refraction (air or detector)
- `nmed::Float64`: the sample medium's index of refraction
- `rd::Array{Float64,1}`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs::Array{Float64,1}`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L::Array{Float64,1}`: the dimensions [lx, ly, lz] of the parallelepied

# Examples
julia> fluence_DA_paralpip_TD(0:1:5, [0.1,10.0], 1.0, 1.0, [2.5,5.0,1.0], [5.0,5.0], [10.0,10.0,10.0])
"""
function fluence_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})
	n = nmed/ndet
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
	ν = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x = rd[1]
	y = rd[2]
	z = rd[3]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	Rt1 = Array{Float64}(undef, length(t))
	Rt2 = zeros(Float64, length(t))
	Rt3 = zeros(Float64, length(t))
	Rt4 = zeros(Float64, length(t))

	if n == 1.0
		A = 1.0
	elseif n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end

	zo = 1/μsp
    zb = 2A*D
    
    Threads.@threads for n in eachindex(t) 
        
        Rt1[n] = @. ν*exp(-μa*ν*t[n])/((4π*D*ν*t[n])^(3/2))

		  for m in xs
			    x1l = 2m.*lx .+ 4m.*zb .+ xu 
				x2l = 2m.*lx + (4m-2).*zb .- xu
				y1m = 2m.*ly .+ 4m.*zb .+ yu
				y2m = 2m.*ly .+ (4m.-2).*zb .- yu
				z1n = 2m.*lz .+ 4m.*zb .+ zo
				z2n = 2m.*lz .+ (4m.-2).*zb .- zo
		  
    	        Rt2[n] += exp(-(x-x1l)^2 / (4D*ν*t[n])) - exp(-(x-x2l)^2 / (4D*ν*t[n]))
	            Rt3[n] += exp(-(y-y1m)^2 / (4D*ν*t[n])) - exp(-(y-y2m)^2 / (4D*ν*t[n]))
	            Rt4[n] += exp(-(z-z1n)^2 /(4D*ν*t[n])) - exp(-(z-z2n)^2 /(4D*ν*t[n]))
          end

          Rt1[n] *= (Rt2[n]*Rt3[n]*Rt4[n])

        if isnan(Rt1[n])
            Rt1[n] = 0
        end

    end
      
	return Rt1
end

### TD Reflectance ###
"""
    refl_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})

Compute the time-domain reflectance from a parallelepiped [lx, ly, lz] applying Fick's law to Eqn. 6 Kienle 05. 

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
julia> refl_DA_paralpip_TD(0:1:5, [0.1,10.0], 1.0, 1.0, [2.5,5.0], [5.0,5.0], [10.0,10.0,10.0])
"""
function refl_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})
	n = nmed/ndet
	μa = β[1]
	μsp = β[2]
	D = 1/3μsp
	ν = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x = rd[1]
	y = rd[2]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
    lz = L[3]
    
    Rt1 = Array{Float64}(undef, length(t))
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
	
	zo = 1/μsp
	zb = 2A*D
  
    Threads.@threads for n in eachindex(t) 
        
        Rt1[n] = @. exp(-μa*ν*t[n])/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

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

        Rt1[n] *= (Rt2[n]*Rt3[n]*Rt4[n])

        if isnan(Rt1[n])
            Rt1[n] = 0
        end

    end
      
	return Rt1
end

### TD Transmittance ###
"""
    trans_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})

Compute the time-domain transmittance from a parallelepiped [lx, ly, lz] with Eqn. 8 Kienle 05. 

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
julia> trans_DA_paralpip_TD(0:1:5, [0.1,10.0], 1.0, 1.0, [2.5,5.0], [5.0,5.0], [10.0,10.0,10.0])
"""
function trans_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})
	n = nmed/ndet
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
	ν = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x = rd[1]
	y = rd[2]
	z = L[3]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	Rt1 = Array{Float64}(undef, length(t))
	Rt2 = zeros(Float64, length(t))
	Rt3 = zeros(Float64, length(t))
	Rt4 = zeros(Float64, length(t))

	if n == 1.0
		A = 1.0
	elseif n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end

	zo = 1/μsp
    zb = 2A*D
    
    Threads.@threads for n in eachindex(t) 
        
        Rt1[n] = @. exp(-μa*ν*t[n])/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

		  for m in xs
			    x1l = 2m.*lx .+ 4m.*zb .+ xu 
				x2l = 2m.*lx + (4m-2).*zb .- xu
				y1m = 2m.*ly .+ 4m.*zb .+ yu
				y2m = 2m.*ly .+ (4m.-2).*zb .- yu
				z1n = 2m.*lz .+ 4m.*zb .+ zo
				z2n = 2m.*lz .+ (4m.-2).*zb .- zo
		  
    	        Rt2[n] += exp(-(x-x1l)^2 / (4D*ν*t[n])) - exp(-(x-x2l)^2 / (4D*ν*t[n]))
	            Rt3[n] += exp(-(y-y1m)^2 / (4D*ν*t[n])) - exp(-(y-y2m)^2 / (4D*ν*t[n]))
	            Rt4[n] += (lz - z1n)*exp(-(lz-z1n)^2 /(4D*ν*t[n])) - (lz - z2n)*exp(-(lz-z2n)^2 /(4D*ν*t[n]))
          end

          Rt1[n] *= (Rt2[n]*Rt3[n]*Rt4[n])

        if isnan(Rt1[n])
            Rt1[n] = 0
        end

    end
      
	return Rt1
end
