#---------------------------------------------------------------------------------------------------------------------------------------- 
# Implements solution to the diffusion equation in a N-layered finite cylinder as given in [1].
# Solutions are given in the spatial, frequency, and time domains for a point source located in the middle of the cylinder top.
#
# [1] André Liemert and Alwin Kienle, "Light diffusion in a turbid cylinder. II. Layered case," Opt. Express 18, 9266-9279 (2010) 
#---------------------------------------------------------------------------------------------------------------------------------------- 

@with_kw struct Nlayer_cylinder{N, T <: Real} <: DiffusionParameters
    μsp::NTuple{N, T} = (10.0, 10.0, 10.0, 10.0)            # reduced scattering coefficient (1/cm)
    μa::NTuple{N, T} = (0.1, 0.1, 0.1, 0.1)                 # absorption coefficient (1/cm)
    n_ext::T = 1.0                                          # surrounding index of refraction
    n_med::NTuple{N, T} = (1.0, 1.0, 1.0, 1.0)              # layers index of refraction

    l::NTuple{N, T} = (0.5, 0.8, 1.0, 5.0)                  # length of cylinder layers (cm)
    ρ::T = 1.0                                              # source-detector separation (cm)
    a::T = 5.0                                              # radius of cylinder (cm)
    z::T = 0.0                                              # detector depth (cm)

    ω::T = 0.0                                              # modulation frequency
    N_J0Roots::Int = 1000                                   # Number of besselj0 roots in sum (N<=1e6)
end

#-------------------------------------------
# Steady-State Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)

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
- `N_J0Roots`: Number of roots of bessel function of first kind zero order J0(x) = 0. Must be less than 1e6

# Examples
julia> `fluence_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, 1000)`
"""
function fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)
    D = D_coeff.(μsp)
    N = length(D)
    A = A_coeff.(n_med ./ n_ext)
    z0 = z0_coeff(μsp[1])
    zb = zb_coeff.(A, D)
    n_med = @. D * n_med^2
    @assert z0 < l[1]
    
    roots = @view J0_ROOTS[1:N_J0Roots]
    if z < l[1]
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, roots, _green_Nlaycylin_top, N)
    elseif z > sum(l[1:end - 1])
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, roots, _green_Nlaycylin_bottom, N)
    end
end
"""
    fluence_DA_Nlay_cylinder_CW(data)

Wrapper to fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `fluence_DA_Nlay_cylinder_CW(data)`
"""
function fluence_DA_Nlay_cylinder_CW(data)
    return fluence_DA_Nlay_cylinder_CW(data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.N_J0Roots)
end

#-------------------------------------------
# Steady-State Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)

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
- `N_J0Roots`: number of roots of bessel function of first kind zero order J0(x) = 0. Must be less than 1e6.

# Examples
julia> `flux_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, 1000)`
"""
function flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, dz, N_J0Roots), z)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, dz, N_J0Roots), z)
    end
end
"""
    flux_DA_Nlay_cylinder_CW(data, N_J0Roots)

Wrapper to flux_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `flux_DA_Nlay_cylinder_CW(data, N_J0Roots)`
"""
function flux_DA_Nlay_cylinder_CW(data)
    return flux_DA_Nlay_cylinder_CW(data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.N_J0Roots)
end

#-------------------------------------------
# Time-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24, ILT = hyper_fixed)

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
- `N_J0Roots`: roots of bessel function of first kind zero order J0(x) = 0
- `N`: number of Hankel-Laplace calculations
- `ILT`: inverse laplace transform function

# Examples
julia> `fluence_DA_Nlay_cylinder_TD(0.1:0.1:2.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, 1000)`
"""
function fluence_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 18, ILT = hyperbola)
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, N_J0Roots), t, N = N)
end
function fluence_DA_Nlay_cylinder_TD(t::AbstractArray, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24, ILT = hyper_fixed)
    T = promote_type(eltype(ρ), eltype(μa), eltype(μsp), eltype(z), eltype(l))
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, N_J0Roots), t, N = N, T = T)
end
"""
    fluence_DA_Nlay_cylinder_TD(t, data; N = 24)

Wrapper to fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24, ILT = hyper_fixed) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `fluence_DA_Nlay_cylinder_TD(0.5:1.0:2.5, data)`
"""
function fluence_DA_Nlay_cylinder_TD(t, data; N = 24)
    return fluence_DA_Nlay_cylinder_TD(t, data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.N_J0Roots, N = N)
end

#-------------------------------------------
# Time-Domain Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24)

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
- `N_J0Roots`: roots of bessel function of first kind zero order J0(x) = 0
- `N`: number of Hankel-Laplace calculations

# Examples
julia> `fluence_DA_Nlay_cylinder_TD(0.1:0.1:2.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, N_J0Roots)`
"""
function flux_DA_Nlay_cylinder_TD(t::AbstractVector, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)

    ## need to know the type of z to preallocate vectors coorectly in hyper_fixed for autodiff
    function _ILT(z::T) where {T} 
        return hyper_fixed(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, N_J0Roots), t, N = N, T = T)
    end
    if z == zero(eltype(z))
        return D[1]*ForwardDiff.derivative(_ILT, z)
    elseif z == sum(l)
        return -D[end]*ForwardDiff.derivative(_ILT, z)
    end
end
function flux_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)  
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, N_J0Roots), t, N = N), 0.0)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, N_J0Roots), t, N = N), sum(l))
    end
end

"""
    flux_DA_Nlay_cylinder_TD(t, data; N = 24)

Wrapper to flux_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots; N = 24) with inputs given as a structure (data).

# Examples
julia> data = Nlayer_cylinder(a = 10.0, l = [1.0, 1.0, 1.0, 2.0], z = 5.0)
julia> `flux_DA_Nlay_cylinder_TD(0.5:1.0:2.5, data)`
"""
function flux_DA_Nlay_cylinder_TD(t, data;  N = 24)
    return flux_DA_Nlay_cylinder_TD(t, data.ρ, data.μa, data.μsp, data.n_ext, data.n_med, data.l, data.a, data.z, data.N_J0Roots, N = N)
end

#-------------------------------------------
# Frequency-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, N_J0Roots)

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
- `N_J0Roots`: roots of bessel function of first kind zero order J0(x) = 0

# Examples
julia> `fluence_DA_Nlay_cylinder_CW(1.0, [0.1, 0.1], [10.0, 10.0], 1.0, [1.0, 1.0], [4.5, 4.5], 10.0, 0.0, 1.0, 1000)`
"""
function fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, N_J0Roots)
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ ω * im ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex, μsp, n_ext, n_med, l, a, z, N_J0Roots)
end
"""
    fluence_DA_Nlay_cylinder_FD(data)

Wrapper to fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp, n_ext, n_med, l, a, z, ω, N_J0Roots) with inputs given as a structure (data).

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
function _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, N_J0Roots)
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ s ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex, μsp, n_ext, n_med, l, a, z, N_J0Roots)
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
    apzb = inv(a + zb[1])

    @inbounds for ind in eachindex(besselroots)
        ϕ_tmp = green(besselroots[ind] * apzb, μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp *= besselj0(besselroots[ind] * apzb * ρ)
        ϕ_tmp /= J1_J0ROOTS_2[ind] # replaces (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ / (π * (a + zb[1])^2)
end
function _kernel_fluence_DA_Nlay_cylinder(ρ::Tuple, D, μa, a, zb, z, z0, l, n_med, besselroots, green, N)
    ϕ = ρ .* zero(eltype(μa))
    ϕ_tmp = zero(eltype(μa))
    apzb = inv(a + zb[1])

    for ind in eachindex(besselroots)
        tmp = besselroots[ind] * apzb
        ϕ_tmp = green(tmp, μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp /= J1_J0ROOTS_2[ind] # replaces (besselj1(besselroots[ind]))^2
        ϕ = @. ϕ + ϕ_tmp * besselj0(tmp * ρ)
    end

    return ϕ ./ (π * (a + zb[1])^2)
end
#-------------------------------------------------------------------------------
# Calculates the Green's function in the first (top) and last (bottom) layer
# sinh and cosh have been expanded as exponentials.
# A common term has been factored out. This cancels for G1 but not GN (see _βγN_correction)
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
#-------------------------------------------------------------------------------
@inline function _green_Nlaycylin_top(sn, μa, D, z, z0, zb, l, n, N)
    α = @. sqrt(μa / D + sn^2)

    if N == 4
        β, γ = _get_βγ4(α, n, zb, l)
    elseif N == 3
        β, γ = _get_βγ3(α, n, zb, l)
    elseif N == 2
        β, γ = _get_βγ2(α, zb, l)
    elseif N > 4
        β, γ = _get_βγk(α, n, zb, l)
    end

    tmp1 = α[1] * n[1] * β
    tmp2 = α[2] * n[2] * γ
    tmp3 = exp(-2 * α[1] * (l[1] + zb[1]))

    g  = exp(-α[1] * abs(z - z0))
    g -= exp(-α[1] * (z + z0 + 2 * zb[1]))
    g1 = exp(α[1] * (z + z0 - 2 * l[1]))
    g1 *= (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * (z + zb[1])))
    g1 *= tmp1 - tmp2
    g1 /= muladd(tmp1, tmp3, tmp1) + muladd(tmp2, -tmp3, tmp2)

    return (g + g1) / (2 * D[1] * α[1])
end
@inline function _green_Nlaycylin_bottom(sn, μa, D, z, z0, zb, l, n, N)
    α = @. sqrt(μa / D + sn^2)

    if N == 4
        β, γ = _get_βγ4(α, n, zb, l)
        βγ_correction = _βγ4_correction(α, zb, l)  
    elseif N == 3
        β, γ = _get_βγ3(α, n, zb, l)
        βγ_correction = _βγ3_correction(α, zb, l)  
    elseif N == 2
        β, γ = _get_βγ2(α, zb, l)
        βγ_correction = _βγ2_correction(α, zb, l)  
    elseif N > 4
        β, γ = _get_βγk(α, n, zb, l)
        βγ_correction = _βγN_correction(α, zb, l)  
    end

    tmp = one(eltype(α))
    for ind in (N - 1):-1:2
        tmp *= α[ind] * n[ind]
    end

    tmp1 = exp(-2 * α[1] * (l[1] + zb[1]))
    gN = n[end] * tmp * 2^(N - 1) / 2 / D[end]
    gN *= exp(α[1] * (z0 - l[1]) + α[end] * (sum(l) + zb[end] - z) - βγ_correction)
    gN /= α[1] * n[1] * β * (1 + tmp1) + α[2] * n[2] * γ * (1 - tmp1)
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
@inline function _get_βγ2(α, zb, l)
    tmp1 = exp(-2 * α[2] * (l[2] + zb[2]))
    
    β = (1 - tmp1)
    γ = (1 + tmp1)
    
    return β, γ
end
@inline function _get_βγ3(α, n, zb, l)
    tmp1 = α[2] * n[2]
    tmp2 = α[3] * n[3]
    tmp3 = exp(-2 * α[2] * l[2])
    tmp4 = exp(-2 * α[3] * (l[3] + zb[2]))

    a = tmp1 * tmp3
    b = muladd(-tmp1, tmp4, tmp1)
    c = a * tmp4
    a = a - c
    β = b + a
    γ = b - a

    a = tmp2 * tmp3
    b = muladd(tmp2, tmp4, tmp2)
    c = a * tmp4
    a = a + c

    β += b - a
    γ += b + a
    
    return β, γ
end
@inline function _get_βγ4(α, n, zb, l)
    tmp5 = α[2] * n[2]
    tmp1 = α[3] * n[3]
    tmp4 = α[4] * n[4]

    tmp6 = exp(-2 * α[2] * l[2])
    tmp2 = exp(-2 * α[3] * l[3])
    tmp3 = exp(-2 * α[4] * (l[4] + zb[2]))

    a = tmp1 * tmp3
    b = tmp1 * tmp2
    c = a * tmp2
    a = tmp1 - a
    b = b - c
    β4 = a + b
    γ4 = a - b

    a = tmp4 * tmp2
    b = tmp4 * tmp3
    c = a * tmp3
    b = tmp4 + b
    c = a + c
    β4 += b - c
    γ4 += b + c
    
    a = tmp5 * β4
    b = a * tmp6
    c = tmp1 * γ4
    d = c * tmp6
    a = a + c
    b = b - d

    β = a + b
    γ = a - b

    return β, γ
end
@inline function _get_βγk(α, n, zb, l)
    βN, γN = _get_βγN(α, n, zb, l)
  
    for k in length(α):-1:4
        tmp1 = α[k - 2] * n[k - 2] * βN
        tmp2 = α[k - 1] * n[k - 1] * γN
        tmp3 = exp(-2 * α[k - 2] * l[k - 2])

        a = tmp1 * tmp3
        c = tmp1 + tmp2
        b = tmp2 * tmp3
        d = a - b

        βN = c + d
        γN = c - d
    end

    return βN, γN
end
@inline function _get_βγN(α, n, zb, l)
    tmp1 = α[end - 1] * n[end - 1]
    tmp2 = α[end] * n[end]
    tmp3 = exp(-2 * α[end - 1] * l[end - 1])
    tmp4 = exp(-2 * α[end] * (l[end] + zb[2]))
    
    a = tmp1 * tmp3
    b = tmp1 - tmp1 * tmp4
    c = a * tmp4
    a = a - c
    βN =  b + a
    γN =  b - a

    a = tmp2 * tmp3
    b = tmp2 + tmp2 * tmp4
    c = a * tmp4
    a = a + c
    βN += b - a
    γN += b + a

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
