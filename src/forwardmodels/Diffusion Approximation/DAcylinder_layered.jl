#---------------------------------------------------------------------------------------------------------------------------------------- 
# Implements solution to the diffusion equation in a N-layered finite cylinder as given in [1].
# Solutions are given in the spatial, frequency, and time domains for a point source located in the middle of the cylinder top.
#
# [1] André Liemert and Alwin Kienle, "Light diffusion in a turbid cylinder. II. Layered case," Opt. Express 18, 9266-9279 (2010) 
#---------------------------------------------------------------------------------------------------------------------------------------- 

@with_kw struct Nlayer_cylinder{T <: Real} <: DiffusionParameters
    μsp::Vector{T} = [10.0, 10.0, 10.0, 10.0]               # reduced scattering coefficient (1/cm)
    μa::Vector{T} = [0.1, 0.1, 0.1, 0.1]                    # absorption coefficient (1/cm)
    n_ext::T = 1.0                                          # surrounding index of refraction
    n_med::Vector{T} = [1.0, 1.0, 1.0, 1.0]                 # layers index of refraction

    l::Vector{T} = [0.5, 0.8, 1.0, 5.0]                     # length of cylinder layers (cm)
    ρ::Union{T, AbstractVector{T}} = 1.0                    # source-detector separation (cm)
    a::T = 5.0                                              # radius of cylinder (cm)
    z::T = 0.0                                              # detector depth (cm)

    ω::T = 0.0                                              # modulation frequency
    bessels::AbstractVector{T} = besselroots[1:1000]        # roots of J0
end

#-------------------------------------------
# Steady-State Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots)

Compute the steady-state fluence in an N-layered cylinder. Source is assumed to be located on the top middle of the cylinder.
No checks on arguments are made except enforcing that the source is located in the top layer.
ρ, μa, and z should be >= 0.0.
μsp, n_ext, n_med, a should be > 0.0.
Each layer in l should be > 0.0.
In general, ρ should always be < a.
All arguments should have the same types and the length of μsp, μa, n_med, l should be equal.
It is advised to use the Nlayer_cylinder() structure format for inputs.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `z`: source depth within cylinder
- `besselroots`: roots of bessel function of first kind zero order J0(x) = 0

# Examples
julia> `fluence_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, besselroots)`
"""
function fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots)
    D = D_coeff.(μsp)
    N = length(D)
    A = A_coeff.(n_med / n_ext)
    z0 = z0_coeff(μsp[1])
    zb = zb_coeff.(A, D)
    @assert z0 < l[1]

    if z < l[1]
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, besselroots, _green_Nlaycylin_top, N) / (π * (a + zb[1])^2)
    elseif z > sum(l[1:end - 1])
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, besselroots, _green_Nlaycylin_bottom, N) / (π * (a + zb[1])^2)
    end
end
"""
    fluence_DA_Nlay_cylinder_CW(data)

Wrapper to fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `fluence_DA_Nlay_cylinder_CW(data)`
"""
function fluence_DA_Nlay_cylinder_CW(data)
    return fluence_DA_Nlay_cylinder_CW(data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.bessels)
end

#-------------------------------------------
# Steady-State Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots)

Compute the steady-state flux using Fick's law D[1]*∂ϕ(ρ)/∂z for z = 0 (reflectance) and -D[end]*∂ϕ(ρ)/∂z for z = sum(l) (transmittance) in an N-layered cylinder. 
Source is assumed to be located on the top middle of the cylinder. z must be equal to 0 or the total length sum(l) of cylinder.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `z`: source depth within cylinder: must be equal to 0 or sum(l)
- `besselroots`: roots of bessel function of first kind zero order J0(x) = 0

# Examples
julia> `flux_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, besselroots)`
"""
function flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, dz, besselroots), z)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, dz, besselroots), z)
    end
end
"""
    flux_DA_Nlay_cylinder_CW(data, besselroots)

Wrapper to flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `flux_DA_Nlay_cylinder_CW(data, besselroots)`
"""
function flux_DA_Nlay_cylinder_CW(data)
    return flux_DA_Nlay_cylinder_CW(data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.bessels)
end

#-------------------------------------------
# Time-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24, ILT = hyper_fixed)

Compute the time-domain fluence in an N-layered cylinder. Source is assumed to be located on the top middle of the cylinder.
See notes on fluence_DA_Nlay_cylinder_CW for input checks with the initial constraint that all values of t should be > 0.0 and ordered least to greatest.
Do not have any value of t be zero or negative.

It is best to use `hyper_fixed` as the inverse laplace transform if t consists of many time points. Utilize `hyperbola` for a single time point.
The value of N should be proportional to the dynamic range of the time-domain signal needed. For later times you will need a larger N.
The lowest fluence value you can compute will be no less than the machine precision you utilize. Float64 values are limited to fluence > ~2.2-16.

# Arguments
- `t`: time point or vector (ns)
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `z`: source depth within cylinder
- `besselroots`: roots of bessel function of first kind zero order J0(x) = 0
- `N`: number of Hankel-Laplace calculations
- `ILT`: inverse laplace transform function

# Examples
julia> `fluence_DA_Nlay_cylinder_TD(0.1:0.1:2.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, besselroots)`
"""
function fluence_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 18, ILT = hyperbola)
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, besselroots), t, N = N)
end
function fluence_DA_Nlay_cylinder_TD(t::AbstractArray, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24, ILT = hyper_fixed)
    T = promote_type(eltype(ρ), eltype(μa), eltype(μsp), eltype(z), eltype(l))
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, besselroots), t, N = N, T = T)
end
"""
    fluence_DA_Nlay_cylinder_TD(t, data; N = 24)

Wrapper to fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24, ILT = hyper_fixed) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `fluence_DA_Nlay_cylinder_TD(0.5:1.0:2.5, data)`
"""
function fluence_DA_Nlay_cylinder_TD(t, data; N = 24)
    return fluence_DA_Nlay_cylinder_TD(t, data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.bessels, N = N)
end

#-------------------------------------------
# Time-Domain Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24)

Compute the time-domain flux using Fick's law D[1]*∂ϕ(ρ, t)/∂z for z = 0 (reflectance) and -D[end]*∂ϕ(ρ, t)/∂z for z = sum(l) (transmittance) in an N-layered cylinder. 
Source is assumed to be located on the top middle of the cylinder. z must be equal to 0 or the total length sum(l) of cylinder.
Uses the hyperbola contour if t is an AbstractFloat and the fixed hyperbola contour if t is an AbstractVector.
∂ϕ(ρ, t)/∂z is calculated using forward mode auto-differentiation with ForwardDiff.jl

# Arguments
- `t`: time point or vector (ns)
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `z`: source depth within cylinder
- `besselroots`: roots of bessel function of first kind zero order J0(x) = 0
- `N`: number of Hankel-Laplace calculations

# Examples
julia> `fluence_DA_Nlay_cylinder_TD(0.1:0.1:2.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, besselroots)`
"""
function flux_DA_Nlay_cylinder_TD(t::AbstractVector, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)

    ## need to know the type of z to preallocate vectors coorectly in hyper_fixed for autodiff
    function _ILT(z::T) where {T} 
        return hyper_fixed(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, besselroots), t, N = N, T = T)
    end
    if z == zero(eltype(z))
        return D[1]*ForwardDiff.derivative(_ILT, z)
    elseif z == sum(l)
        return -D[end]*ForwardDiff.derivative(_ILT, z)
    end
end
function flux_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)  
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, besselroots), t, N = N), 0.0)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, besselroots), t, N = N), sum(l))
    end
end

"""
    flux_DA_Nlay_cylinder_TD(t, data; N = 24)

Wrapper to flux_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, besselroots; N = 24) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `flux_DA_Nlay_cylinder_TD(0.5:1.0:2.5, data)`
"""
function flux_DA_Nlay_cylinder_TD(t, data;  N = 24)
    return flux_DA_Nlay_cylinder_TD(t, data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.bessels, N = N)
end

#-------------------------------------------
# Frequency-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, besselroots)

Compute the frequency modulated fluence in an N-layered cylinder. Source is assumed to be located on the top middle of the cylinder.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `z`: source depth within cylinder
- `ω`: source modulation frequency (1/ns)
- `besselroots`: roots of bessel function of first kind zero order J0(x) = 0

# Examples
julia> `fluence_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, 1.0, besselroots)`
"""
function fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, besselroots)
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ ω * im ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex, μsp, n_ext, n_med, l, a, z, besselroots)
end
"""
    fluence_DA_Nlay_cylinder_FD(data)

Wrapper to fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, besselroots) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0, ω = 1.0)
julia> `fluence_DA_Nlay_cylinder_FD(data)`
"""
function fluence_DA_Nlay_cylinder_FD(data)
    return fluence_DA_Nlay_cylinder_FD(data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.ω, data.bessels)
end

#-------------------------------------------------------------------------------
# this function is the base for the laplace transform in the time-domain
#-------------------------------------------------------------------------------
function _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, besselroots)
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ s ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex, μsp, n_ext, n_med, l, a, z, besselroots)
end

#-------------------------------------------------------------------------------
# _kernel_fluence_DA_Nlay_cylinder provides the main kernel of the program.
# The inputs are similar to fluence_DA_Nlay_cylinder_CW.
# D is the diffusion coefficient and N is the number of layers.
# green is the Green's function for either the first or bottom layer (below).
#-------------------------------------------------------------------------------
function _kernel_fluence_DA_Nlay_cylinder(ρ::AbstractFloat, D, μa, a, zb, z, z0, l, n_med, besselroots, green, N)
    ϕ = zero(eltype(μa))
    ϕ_tmp = zero(eltype(μa))
    α = zeros(eltype(μa), N)

    for ind in eachindex(besselroots)
        ϕ_tmp = green(α, besselroots[ind] / (a + zb[1]), μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp *= besselj0(besselroots[ind] / (a + zb[1]) * ρ)
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ
end

function _kernel_fluence_DA_Nlay_cylinder(ρ::AbstractVector, D, μa, a, zb, z, z0, l, n_med, besselroots, green, N)
    ϕ = zeros(eltype(μa), length(ρ))
    ϕ_tmp = zero(eltype(μa))
    α = zeros(eltype(μa), N)

    for ind in eachindex(besselroots)
        tmp = besselroots[ind] / (a + zb[1])
        ϕ_tmp = green(α, tmp, μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        for ρ_ind in eachindex(ρ)
            ϕ[ρ_ind] += ϕ_tmp * besselj0(tmp * ρ[ρ_ind])
        end
    end

    return ϕ
end
#-------------------------------------------------------------------------------
# Calculates the Green's function in the first (top) and last (bottom) layer
# sinh and cosh have been expanded as exponentials.
# A common term has been factored out. This cancels for G1 but not GN (see _βγN_correction)
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
#-------------------------------------------------------------------------------
@inline function α_coeff!(α, μa, D, sn)
    @inbounds for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2)
    end
    return α
end
@inline function _green_Nlaycylin_top(α, sn, μa, D, z, z0, zb, l, n, N)
    α = α_coeff!(α, μa, D, sn)

    if N == 4
        β, γ = _get_βγ4(α, D, n, zb, l)
    elseif N == 3
        β, γ = _get_βγ3(α, D, n, zb, l)
    elseif N == 2
        β, γ = _get_βγ2(α, D, n, zb, l)
    elseif N > 4
        β, γ = _get_βγk(α, D, n, zb, l)
    end

    tmp1 = D[1] * α[1] * n[1]^2 * β
    tmp2 = D[2] * α[2] * n[2]^2 * γ
    tmp3 = exp(-2 * α[1] * (l[1] + zb[1]))

    g  = (exp(-α[1] * abs(z - z0)) - exp(-α[1] * (z + z0 + 2 * zb[1])))
    g1 = exp(α[1] * (z + z0 - 2 * l[1])) * (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * (z + zb[1])))
    g1 *= (tmp1 - tmp2)
    g1 /= (tmp1 * (1 + tmp3) + tmp2 * (1 - tmp3))

    return (g + g1) / (2 * D[1] * α[1])
end
@inline function _green_Nlaycylin_bottom(α, sn, μa, D, z, z0, zb, l, n, N)
    α = α_coeff!(α, μa, D, sn)

    if N == 4
        β, γ = _get_βγ4(α, D, n, zb, l)
        βγ_correction = _βγ4_correction(α, zb, l)  
    elseif N == 3
        β, γ = _get_βγ3(α, D, n, zb, l)
        βγ_correction = _βγ3_correction(α, zb, l)  
    elseif N == 2
        β, γ = _get_βγ2(α, D, n, zb, l)
        βγ_correction = _βγ2_correction(α, zb, l)  
    elseif N > 4
        β, γ = _get_βγk(α, D, n, zb, l)
        βγ_correction = _βγN_correction(α, zb, l)  
    end

    tmp = 1.0
    for ind in (N - 1):-1:2
        tmp *= D[ind] * α[ind] * n[ind]^2
    end

    tmp1 = exp(-2 * α[1] * (l[1] + zb[1]));
    gN = n[end]^2 * tmp * 2^(N - 1) / 2
    gN *= exp(α[1] * (z0 - l[1]) + α[end] * (sum(l) + zb[end] - z) - βγ_correction)
    gN /= D[1] * α[1] * n[1]^2 * β * (1 + tmp1) + D[2] * α[2] * n[2]^2 * γ * (1 - tmp1)
    gN *= (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[end] * (sum(l) + zb[end] - z)))

    return gN
 end

#-------------------------------------------------------------------------------
# Calculate β and γ coefficients with eqn. 17 in [1].
# sinh and cosh have been expanded as exponentials.
# A common term has been factored out. This cancels for G1 but not GN (see _βγN_correction)
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
#-------------------------------------------------------------------------------
@inline function _get_βγ2(α, D, n, zb, l)
    tmp1 = exp(-2 * α[2] * (l[2] + zb[2]))
    
    β = (1 - tmp1)
    γ = (1 + tmp1)
    
    return β, γ
end
@inline function _get_βγ3(α, D, n, zb, l)
    tmp1 = D[2] * α[2] * n[2]^2
    tmp2 = D[3] * α[3] * n[3]^2
    tmp3 = exp(-2 * α[2] * l[2])
    tmp4 = exp(-2 * α[3] * (l[3] + zb[2]))

    β  = tmp1 * (1 + tmp3) * (1 - tmp4)
    β += tmp2 * (1 - tmp3) * (1 + tmp4)

    γ  = tmp1 * (1 - tmp3) * (1 - tmp4)
    γ += tmp2 * (1 + tmp3) * (1 + tmp4)

    # used muladd but slower
    # β  = muladd((1 + tmp3) * (1 - tmp4), tmp1, tmp2 * (1 - tmp3) * (1 + tmp4))
    

    return β, γ
end
@inline function _get_βγ4(α, D, n, zb, l)
    tmp5 = D[2] * α[2] * n[2]^2
    tmp1 = D[3] * α[3] * n[3]^2
    tmp4 = D[4] * α[4] * n[4]^2

    tmp6 = exp(-2 * α[2] * l[2])
    tmp2 = exp(-2 * α[3] * l[3])
    tmp3 = exp(-2 * α[4] * (l[4] + zb[2]))

    β_4  = tmp1 * (1 + tmp2) * (1 - tmp3)
    β_4 += tmp4 * (1 - tmp2) * (1 + tmp3)

    γ_4  = tmp1 * (1 - tmp2) * (1 - tmp3)
    γ_4 += tmp4 * (1 + tmp2) * (1 + tmp3)

    β = (tmp5 * β_4 * (1 + tmp6) + tmp1 * γ_4 * (1 - tmp6))
    γ = (tmp5 * β_4 * (1 - tmp6) + tmp1 * γ_4 * (1 + tmp6))

    return β, γ
end
@inline function _get_βγk(α, D, n, zb, l)
    βN, γN = _get_βγN(α, D, n, zb, l)
  
    for k in length(α):-1:4
        tmp1 = D[k - 2] * α[k - 2] * n[k - 2]^2 * βN
        tmp2 = D[k - 1] * α[k - 1] * n[k - 1]^2 * γN
        tmp3 = exp(-2 * α[k - 2] * l[k - 2])

        # βN  =  tmp1 * (1 + tmp3)
        # βN +=  tmp2 * (1 - tmp3)
        # γN  =  tmp1 * (1 - tmp3)
        # γN +=  tmp2 * (1 + tmp3)


        # 2 ns optimization for median time of "fluence_DA_Nlay_cylinder_CW"
        a  =  muladd(tmp1, tmp3, tmp1)
        b  =  muladd(tmp2, -tmp3, tmp2)
        βN =  a + b
        c  =  muladd(tmp1, -tmp3, tmp1)
        d  =  muladd(tmp2, tmp3, tmp2)
        γN  =  c + d

    end

    return βN, γN
end
@inline function _get_βγN(α, D, n, zb, l)
    tmp1 = D[end - 1] * α[end - 1] * n[end - 1]^2
    tmp2 = D[end] * α[end] * n[end]^2
    tmp3 = exp(-2 * α[end - 1] * l[end - 1])
    tmp4 = exp(-2 * α[end] * (l[end] + zb[2]))
    
    βN  =   tmp1 * (1 + tmp3) *  (1 - tmp4)
    βN +=   tmp2 * (1 - tmp3) *  (1 + tmp4)

    γN =  tmp1 * (1 - tmp3) *  (1 - tmp4)
    γN +=  tmp2 * (1 + tmp3) *  (1 + tmp4)

    return βN, γN
end

#-------------------------------------------------------------------------------
# Calculate βγ correction for N layer green's function.
# This is done because a common term was factored out from general expressions
# They don't cancel for the N layer so have to divide them out.
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively
#-------------------------------------------------------------------------------
@inline function _βγ2_correction(α, zb, l)  
    return α[2] * (l[2] + zb[2])
end
function _βγ3_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * (l[3] + zb[2])
end
@inline function _βγ4_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * l[3] + α[4] * (l[4] + zb[2])
end
@inline function _βγ5_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * l[3] + α[4] * l[4] + α[5] * (l[5] + zb[2])
end
@inline function _βγN_correction(α, zb, l) 
    out = α[end] * (l[end] + zb[2])
    for ind in (length(α) - 1):-1:2
        out += α[ind] * l[ind]
    end
    return out
end
#-------------------------------------------------------------------------------
