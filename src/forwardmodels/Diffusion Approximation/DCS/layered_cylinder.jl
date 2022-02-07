"""
Nlayer_cylinder_DCS

Provides default parameters to use to compute the autocorrelation function g2 in a N-layered cylinder.
The arguments defined with a Vector should be the same length and same type.

# Keyword Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa::Vector{T}`: absorption coefficient (cm⁻¹)
- `μsp::Vector{T}`: reduced scattering coefficient (cm⁻¹)
- `n_med::Vector{T}`: medium's index of refraction
- `n_ext`: external medium's index of refraction (air or detector)
- `l`: the thicknesses of each layer (cm)
- `z`: the z-depth orthogonal from the boundary (cm)
- `a`: the radius of the cylinder (cm)
- `β`: constant in Siegert relation dependent on collection optics
- `BFi::Vector{T}`: Blood flow index ~αDb (cm²/s)
- `λ`: wavelength (nm)
- `N_J0Roots`: roots of J0

# Examples
```
julia> data = Nlayer_cylinder_DCS() # return default parameters
julia> data = Nlayer_cylinder_DCS(ρ = 1.5) # return ρ = 1.5 with the rest of the parameters given by defaults

# we can now define our correlation times τ
julia> τ = 10 .^(range(-10,stop=0,length=250))
julia> g2_DA_Nlay_cylinder_CW(τ, data)
```
"""
@with_kw struct Nlayer_cylinder_DCS{N, T <: Real}
    μsp::NTuple{N, T} = (10.0, 10.0)                    # reduced scattering coefficient (1/cm)
    μa::NTuple{N, T} = (0.1, 0.1)                       # absorption coefficient (1/cm)
    n_ext::T = 1.0                                      # surrounding index of refraction
    n_med::NTuple{N, T} = (1.0, 1.0)                    # layers index of refraction

    l::NTuple{N, T} = (0.5, 10.0)                       # length of cylinder layers (cm)
    ρ::T = 1.0                                          # source-detector separation (cm)
    a::T = 20.0                                         # radius of cylinder (cm)
    z::T = 0.0                                          # detector depth (cm)

    β::T = 1.0                                          # constant in Siegert relation dependent on collection optics
    BFi::NTuple{N, T} = (2.0e-8, 2.0e-8)                   # Blood flow index ~αDb (cm²/s)
    λ::T = 750.0                                        # wavelength (nm)

    N_J0Roots::Int = 600                                # Number of besselj0 roots in sum (N<=1e6)

    N_quad::Int = 50                                    # number of nodes in gauss-legendre integration
    N_laplace::Int = 12                                 # number of laplace evaluations in time-domain inversion
end

"""
    g1_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], n_ext = 1.0, n_med = [1.0, 1.0], l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 1:500)

Compute the electric field autocorrelation function function g1 in a layered cylinder (CW).

# Arguments
- `τ`: time lags (s)
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `BFi`: Blood flow index ~αDb (cm²/s)
- `β`: constant in Siegert relation dependent on collection optics
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: length of cylinder layers (cm)
- `a`: cylinder radius (cm)
- `z`: the z-depth orthogonal from the boundary (cm)
- `λ`: wavelength (nm)
- `N_J0Roots`: roots of J0 (bessel function of first kind order zero)
"""
function g1_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 500)
    g1 = similar(τ)

    k0 = @. 2 * π * n_med / (λ * 1e-7) # this should be in cm
    tmp = @. 2 * μsp * BFi * k0^2

    G0 = fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)

    Threads.@threads for ind in eachindex(τ)
        μa_dynamic = @. μa + tmp * τ[ind]
        G1 = fluence_DA_Nlay_cylinder_CW(ρ, μa_dynamic, μsp, n_ext, n_med, l, a, z, N_J0Roots)
        g1[ind] = G1 / G0
    end

    return g1
end

"""
    g1_DA_Nlay_cylinder_TD(τ::AbstractVector, t::AbstractFloat, ρ, μa, μsp; BFi = [2e-8, 2e-8], n_ext = 1.0, n_med = [1.0, 1.0], l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 1000, N_laplace = 8)
    g1_DA_Nlay_cylinder_TD(τ::AbstractVector, t::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], n_ext = 1.0, n_med = [1.0, 1.0], l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 1000, N_laplace = 8, N_quad = 100)

Compute the time-domain electric field autocorrelation function function g1 in a layered cylinder.

If t is a AbstractVector then g1(tau) will be calculated by ∫P(t)g1(tau, t)dt with integration bounds given by t[1] and t[end].
This is needed for time-gated applications.

# Arguments
- `τ`: time lags (s)
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `BFi`: Blood flow index ~αDb (cm²/s)
- `β`: constant in Siegert relation dependent on collection optics
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `l`: length of cylinder layers (cm)
- `a`: cylinder radius (cm)
- `z`: the z-depth orthogonal from the boundary (cm)
- `λ`: wavelength (nm)
- `N_J0Roots`: roots of J0 (bessel function of first kind order zero)
- `N_quad`: number of nodes in gauss-legendre integration
- `N_laplace`: number of laplace evaluations in inversion to time-domain
"""
function g1_DA_Nlay_cylinder_TD(τ::AbstractVector, t::AbstractFloat, ρ, μa, μsp; BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 1000, N_laplace = 8)
    g1 = similar(τ)

    k0 = @. 2 * π * n_med / (λ * 1e-7) # this should be in cm
    tmp = @. 2 * μsp * BFi * k0^2

    G0 = fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots, N = N_laplace)

    for ind in eachindex(τ)
        μa_dynamic = @. μa + tmp * τ[ind]
        G1 = fluence_DA_Nlay_cylinder_TD(t, ρ, μa_dynamic, μsp, n_ext, n_med, l, a, z, N_J0Roots, N = N_laplace)
        g1[ind] = G1 / G0
    end

    return g1
end

function g1_DA_Nlay_cylinder_TD(τ::AbstractVector, t::AbstractVector, ρ, μa, μsp; BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 1000, N_laplace = 8, N_quad = 100)
    g1 = similar(τ)

    k0 = @. 2 * π * n_med / (λ * 1e-7) # this should be in cm
    tmp = @. 2 * μsp * BFi * k0^2
    x, w = gausslegendre(N_quad)

    tpsf_norm = integrate_compute_y_first(t -> fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots, N = N_laplace), x, w, t[1], t[2])

    for ind in eachindex(τ)
        μa_dynamic = @. μa + tmp * τ[ind]
        g1[ind] = integrate_compute_y_first(t -> (fluence_DA_Nlay_cylinder_TD(t, ρ, μa_dynamic, μsp, n_ext, n_med, l, a, z, N_J0Roots, N = N_laplace) / tpsf_norm), x, w, t[1], t[2])
    end

    return g1
end

function g2_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 500, β = 1.0)
    g1 = g1_DA_Nlay_cylinder_CW(τ, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, N_J0Roots = N_J0Roots)
    for ind in eachindex(τ)
        g1[ind] = 1 + β * abs(g1[ind])^2
    end
    return g1
end
function g2_DA_Nlay_cylinder_CW(τ::AbstractVector, data::Nlayer_cylinder_DCS)
    return g2_DA_Nlay_cylinder_CW(τ, data.ρ, data.μa, data.μsp, BFi = data.BFi, n_ext = data.n_ext, n_med = data.n_med, l = data.l, a = data.a, λ = data.λ, N_J0Roots = data.N_J0Roots, β = data.β)
end

function g2_DA_Nlay_cylinder_TD(τ::AbstractVector, t::Union{AbstractVector, AbstractFloat}, ρ, μa, μsp; BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 500, N_laplace = 12, N_quad = 50, β = 1.0)
    if t isa AbstractVector
        g1 = g1_DA_Nlay_cylinder_TD(τ, t, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, N_J0Roots = N_J0Roots, N_laplace = N_laplace, N_quad = N_quad)
    elseif t isa AbstractFloat
        g1 = g1_DA_Nlay_cylinder_TD(τ, t, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, N_J0Roots = N_J0Roots, N_laplace = N_laplace)
    end

    for ind in eachindex(τ)
        g1[ind] = 1 + β * abs(g1[ind])^2
    end
    return g1
end
function g2_DA_Nlay_cylinder_TD(τ::AbstractVector, t, data::Nlayer_cylinder_DCS)
    return g2_DA_Nlay_cylinder_TD(τ, t, data.ρ, data.μa, data.μsp, BFi = data.BFi, n_ext = data.n_ext, n_med = data.n_med, l = data.l, a = data.a, λ = data.λ, N_J0Roots = data.N_J0Roots, N_laplace = data.N_laplace, N_quad = data.N_quad, β = data.β)
end
