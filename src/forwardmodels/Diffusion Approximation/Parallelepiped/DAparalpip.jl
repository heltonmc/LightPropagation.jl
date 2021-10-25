#----------------------------------------------------------------------------------------------------------------------------------------
# Implements solution to the diffusion equation in the steady-state and time-domain through a turbid parallelepiped [1]

# [1] Alwin Kienle, "Light diffusion through a turbid parallelepiped," J. Opt. Soc. Am. A 22, 1883-1888 (2005) 
#----------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------
# Steady-State Fluence 
#--------------------------------------
"""
    fluence_DA_paralpip_CW(μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)

Compute the steady-state fluence in a parallelepiped [lx, ly, lz].

# Arguments
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `rd`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L`: the dimensions [lx, ly, lz] of the parallelepiped
- `xs`: the number of sources to compute in the series

# Examples
julia> `fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 20)`
"""
function fluence_DA_paralpip_CW(μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    μeff = sqrt(3 * μa * μsp)

    return _kernel_fluence_DA_paralpip_CW(μeff, params.D, params.zb, rd, rs, L, params.z0, xs)
end

function _kernel_fluence_DA_paralpip_CW(μeff, D, zb, rd, rs, L, z0, xs)
    x = rd[1]
    y = rd[2]
    z = rd[3]
    xu = rs[1]
    yu = rs[2]
    lx = L[1]
    ly = L[2]
    lz = L[3]

    ϕ = zero(eltype(D))

    for l in -xs:xs, m in -xs:xs, n in -xs:xs
    	ϕ += _sources_DA_paralpip_CW(l, m, n, x, y, z, lx, ly, lz, xu, yu, z0, zb, μeff)
    end
    
    return ϕ / (4 * π * D)
end
function _sources_DA_paralpip_CW(l, m, n, x, y, z, lx, ly, lz, xu, yu, z0, zb, μeff)
	
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

#--------------------------------------
# Time-Domain Fluence 
#--------------------------------------
"""
    fluence_DA_paralpip_TD(t, μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)

Compute the time-domain fluence in a parallelepiped [lx, ly, lz].

# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `rd`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L`: the dimensions [lx, ly, lz] of the parallelepiped
- `xs`: the number of sources to compute in the series

# Examples
```
julia> fluence_DA_paralpip_TD(0.5, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 20)
```
"""
function fluence_DA_paralpip_TD(t, μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t -> _kernel_fluence_DA_paralpip_TD(params.D, params.ν, t, μa, params.zb, rd, rs, L, params.z0, xs), t)
end
function _kernel_fluence_DA_paralpip_TD(D, ν, t, μa, zb, rd, rs, L, z0, xs)
    x = rd[1]
    y = rd[2]
    z = rd[3]
    xu = rs[1]
    yu = rs[2]
    lx = L[1]
    ly = L[2]
    lz = L[3]

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

### TD Reflectance ###
"""
    flux_DA_paralpip_TD(t, μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = -10:10)

Compute the time-domain flux from a parallelepiped [lx, ly, lz] applying Fick's law to Eqn. 6 Kienle 05. 


# Arguments
- `t`: the time vector (ns). 
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `rd`: target location within medium [x,y,z] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `rs`: the location of the source [x,y] with origin x,y = 0 at corner of parallelepiped; z assumed to be 0
- `L`: the dimensions [lx, ly, lz] of the parallelepiped
- `xs`: the number of sources to compute in the series


# Examples
julia> `flux_DA_paralpip_TD(0.1:0.1:2.0, 0.1, 10.0)`
"""
function flux_DA_paralpip_TD(t, μa, μsp; n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10)
    params = DiffusionKernelParams(μsp, n_med, n_ext)
    return map(t -> _flux_DA_paralpip_TD(t, μa, params.D, params.ν, params.z0, params.zb, rd, rs, L, xs), t)
end

function _flux_DA_paralpip_TD(t, μa, D, ν, z0, zb, rd, rs, L, xs)
    x = rd[1]
    y = rd[2]
    xu = rs[1]
    yu = rs[2]
    lx = L[1]
    ly = L[2]
    lz = L[3]

    if rd[3] == zero(eltype(μa))
        f = _refl_DA_paralpip_TD
    elseif rd[3] == L[3]
        f =  _trans_DA_paralpip_TD
    else 
        @assert z == zero(eltype(z)) || z == s
    end

    Rt1 = exp(-μa * ν * t) / (2 * (4 * π * D * ν)^(3/2) * t^(5/2))
    Rt2 = zero(eltype(Rt1)) 
    Rt3 = zero(eltype(Rt1)) 
    Rt4 = zero(eltype(Rt1))
    tmp = 4 * D * ν * t

    for m in -xs:xs
        x1l = 2 * m * lx + 4 * m * zb + xu 
        x2l = 2 * m * lx + (4 * m - 2) * zb - xu
        y1m = 2 * m * ly + 4 * m * zb + yu
        y2m = 2 * m * ly + (4 * m - 2) * zb - yu
        z1n = 2 * m * lz + 4 * m * zb + z0
        z2n = 2 * m * lz + (4 * m - 2) * zb - z0
		  
        Rt2 += exp(-(x - x1l)^2 / tmp) - exp(-(x - x2l)^2 / tmp)
        Rt3 += exp(-(y - y1m)^2 / tmp) - exp(-(y - y2m)^2 / tmp)
        Rt4 += f(z1n, tmp, z2n, lz)
    end

    Rt1 *= Rt2 * Rt3 * Rt4
	return Rt1
end

_refl_DA_paralpip_TD(z1n, tmp, z2n, lz) = z1n * exp(-z1n^2 / tmp) - z2n * exp(-z2n^2 / tmp)
_trans_DA_paralpip_TD(z1n, tmp, z2n, lz) = (lz - z1n) * exp(-(lz - z1n)^2 / tmp) - (lz - z2n) * exp(-(lz - z2n)^2 / tmp)
