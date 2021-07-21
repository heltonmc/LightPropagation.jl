#########################################################################################################################################
# Implements solution to the diffusion equation in the time-domain through a turbid parallelepiped [1]

# [1] Alwin Kienle, "Light diffusion through a turbid parallelepiped," J. Opt. Soc. Am. A 22, 1883-1888 (2005) 
#########################################################################################################################################

#####################################
# Time-Domain Fluence 
#####################################
"""
    fluence_DA_paralpip_TD(t, μa, μsp; n_ext, n_med, rd, rs, L, xs)

Compute the time-domain fluence in a parallelepiped [lx, ly, lz].

# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ρ`: the source detector separation (cm⁻¹)
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction
- `rd`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L`: the dimensions [lx, ly, lz] of the parallelepiped
- `xs`: the number of sources to compute in the series

# Examples
julia> fluence_DA_paralpip_TD(0.5, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 20)
"""
function fluence_DA_paralpip_TD(t, μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)
    D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
    ν = ν_coeff(n_med)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)

    x = rd[1]
	y = rd[2]
	z = rd[3]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	if isa(t, AbstractFloat)
		return _kernel_fluence_DA_paralpip_TD(D, ν, t, μa, zb, x, y, z, lx, ly, lz, xu, yu, z0, xs)
	elseif isa(t, AbstractArray)
		ϕ = zeros(eltype(D), length(t))
        @inbounds Threads.@threads for ind in eachindex(t)
    		ϕ[ind] = _kernel_fluence_DA_paralpip_TD(D, ν, t[ind], μa, zb, x, y, z, lx, ly, lz, xu, yu, z0, xs)
    	end
    	return ϕ
	end
end
@inline function _kernel_fluence_DA_paralpip_TD(D, ν, t, μa, zb, x, y, z, lx, ly, lz, xu, yu, z0, xs)
    tmp = 4 * D * ν * t
    ϕ1 = ν * exp(-μa * ν * t) / ((π * tmp)^(3/2))
	ϕ2 = zero(eltype(ϕ1)) 
	ϕ3 = zero(eltype(ϕ1)) 
	ϕ4 = zero(eltype(ϕ1)) 

    for m in -xs:xs
    	ϕ2, ϕ3, ϕ4 = _sources_DA_paralpip_TD!(ϕ2, ϕ3, ϕ4, m, zb, lx, ly, lz, xu, yu, z0, x, y, z, tmp)
    end
    
	return ϕ1 * ϕ2 * ϕ3 * ϕ4
end
@inline function _sources_DA_paralpip_TD!(ϕ2, ϕ3, ϕ4, m, zb, lx, ly, lz, xu, yu, z0, x, y, z, tmp)
	tmp1 = 4 * m * zb
	tmp2 = (4 * m - 2) * zb
	tmp3 = 2 * m

	x1l = tmp3 * lx + tmp1 + xu 
	x2l = tmp3 * lx + tmp2 - xu
	y1m = tmp3 * ly + tmp1 + yu
	y2m = tmp3 * ly + tmp2 - yu
	z1n = tmp3 * lz + tmp1 + z0
	z2n = tmp3 * lz + tmp2 - z0
		  
    ϕ2 += exp(-(x - x1l)^2 / tmp) - exp(-(x - x2l)^2 / tmp)
	ϕ3 += exp(-(y - y1m)^2 / tmp) - exp(-(y - y2m)^2 / tmp)
	ϕ4 += exp(-(z - z1n)^2 / tmp) - exp(-(z - z2n)^2 / tmp)
	
	return ϕ2, ϕ3, ϕ4
end


#####################################
# Steady-State Fluence 
#####################################
"""
    fluence_DA_paralpip_CW(μa, μsp; n_ext, n_med, rd, rs, L, xs)

Compute the steady-state fluence in a parallelepiped [lx, ly, lz].

# Arguments
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ρ`: the source detector separation (cm⁻¹)
- `ndet`: the boundary's index of refraction (air or detector)
- `nmed`: the sample medium's index of refraction
- `rd`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L`: the dimensions [lx, ly, lz] of the parallelepiped
- `xs`: the number of sources to compute in the series

# Examples
julia> fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 20)
"""
function fluence_DA_paralpip_CW(μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)
    D = D_coeff(μsp, μa)
    A = A_coeff(n_med / n_ext)
	
    z0 = z0_coeff(μsp)
    zb = zb_coeff(A, D)
	μeff = sqrt(3 * μa * μsp)


    x = rd[1]
	y = rd[2]
	z = rd[3]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	return _kernel_fluence_DA_paralpip_CW(μeff, D, zb, x, y, z, lx, ly, lz, xu, yu, z0, xs)
end
@inline function _kernel_fluence_DA_paralpip_CW(μeff, D, zb, x, y, z, lx, ly, lz, xu, yu, z0, xs)
    ϕ = zero(eltype(D))

    for l in -xs:xs, m in -xs:xs, n in -xs:xs
    	ϕ += _sources_DA_paralpip_CW(l, m, n, x, y, z, lx, ly, lz, xu, yu, z0, zb, μeff)
    end
    
	return ϕ / (4 * π * D)
end
@inline function _sources_DA_paralpip_CW(l, m, n, x, y, z, lx, ly, lz, xu, yu, z0, zb, μeff)
	
	x1l = 2 * l * lx + 4 * l * zb + xu 
	x2l = 2 * l * lx + (4 * l - 2) * zb - xu
	y1m = 2 * m * ly + 4 * m * zb + yu
	y2m = 2 * m * ly + (4 * m - 2) * zb - yu
	z1n = 2 * n * lz + 4 * n * zb + z0
	z2n = 2 * n * lz + (4 * n - 2) * zb - z0

	r1 = sqrt((x - x1l)^2 + (y - y1m)^2 + (z - z1n)^2)
	r2 = sqrt((x - x1l)^2 + (y - y1m)^2 + (z - z2n)^2)
	r3 = sqrt((x - x1l)^2 + (y - y2m)^2 + (z - z1n)^2)
	r4 = sqrt((x - x1l)^2 + (y - y2m)^2 + (z - z2n)^2)
	r5 = sqrt((x - x2l)^2 + (y - y1m)^2 + (z - z1n)^2)
	r6 = sqrt((x - x2l)^2 + (y - y1m)^2 + (z - z2n)^2)
	r7 = sqrt((x - x2l)^2 + (y - y2m)^2 + (z - z1n)^2)
	r8 = sqrt((x - x2l)^2 + (y - y2m)^2 + (z - z2n)^2)


	ϕ  = exp(-μeff * r1) / r1 - exp(-μeff * r2) / r2 - exp(-μeff * r3) / r3 + exp(-μeff * r4) / r4
	ϕ += -exp(-μeff * r5) / r5 + exp(-μeff * r6) / r6 + exp(-μeff * r7) / r7 - exp(-μeff * r8) / r8

	return ϕ
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
