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
- `bessels::AbstractVector{T}`: roots of J0

# Examples
```
julia> data = Nlayer_cylinder_DCS() # return default parameters
julia> data = Nlayer_cylinder_DCS(ρ = 1.5) # return ρ = 1.5 with the rest of the parameters given by defaults

# we can now define our correlation times τ
julia> τ = 10 .^(range(-10,stop=0,length=250))
julia> g2_DA_Nlay_cylinder_CW(τ, data)
```
"""
@with_kw struct Nlayer_cylinder_DCS{T <: Real}
    μsp::Vector{T} = [10.0, 10.0]                       # reduced scattering coefficient (1/cm)
    μa::Vector{T} = [0.1, 0.1]                          # absorption coefficient (1/cm)
    n_ext::T = 1.0                                      # surrounding index of refraction
    n_med::Vector{T} = [1.0, 1.0]                       # layers index of refraction

    l::Vector{T} = [0.5, 10.0]                          # length of cylinder layers (cm)
    ρ::Union{T, AbstractVector{T}} = 1.0                # source-detector separation (cm)
    a::T = 20.0                                         # radius of cylinder (cm)
    z::T = 0.0                                          # detector depth (cm)

    β::T = 1.0                                          # constant in Siegert relation dependent on collection optics
    BFi::Vector{T} = [2.0e-8, 2.0e-8]                   # Blood flow index ~αDb (cm²/s)
    λ::T = 750.0                                        # wavelength (nm)

    bessels::AbstractVector{T} = besselroots[1:1000]    # roots of J0
end

"""
g2_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], β = 1.0, n_ext = 1.0, n_med = 1.0, l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, bessels = besselroots[1:1000])

Compute the electric field autocorrelation function function g2 in a N-layered cylinder. μa, μsp, n_med and BFi must be vectors of the same type.

# Arguments
- `τ`: time lags (s)
- `ρ`: the source detector separation (cm⁻¹)
- `μa::Vector{T}`: absorption coefficient (cm⁻¹)
- `μsp::Vector{T}`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `n_med::Vector{T}`: medium's index of refraction
- `n_ext`: external medium's index of refraction (air or detector)
- `z`: the z-depth orthogonal from the boundary (cm)
- `l`: the thicknesses of each layer (cm)
- `a`: the radius of the cylinder (cm)
- `β`: constant in Siegert relation dependent on collection optics
- `BFi::Vector{T}`: Blood flow index ~αDb (cm²/s)
- `λ`: wavelength (nm)
- `bessels::AbstractVector{T}`: roots of J0

# Examples
```
julia> τ = 10 .^(range(-10,stop=0,length=250)) # need to define τ vector first
julia> g2_DA_Nlay_cylinder_CW(τ, 1.0, [0.1, 0.1], [10.0, 10.0], BFi = [2.1e-8, 3.0e-8]) # only define BFi
julia> g2_DA_Nlay_cylinder_CW(τ, 1.0, [0.1, 0.1], [10.0, 10.0], BFi = [2.1e-8, 3.0e-8], l = [0.8, 3.0]) # define BFi with layer thickness
```
"""
function g2_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], β = 1.0, n_ext = 1.0, n_med = [1.0, 1.0], l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, bessels = besselroots[1:1000])
    g2 = similar(τ)
    μa_dynamic = similar(μa)

    k0 = 2 * π * n_med / (λ * 1e-7) # this should be in cm
    tmp = @. 2 * μsp * BFi * k0^2

    G0 = fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, bessels)

    Threads.@threads for ind in eachindex(τ)
        μa_dynamic = μa + tmp * τ[ind]
        G1 = fluence_DA_Nlay_cylinder_CW(ρ, μa_dynamic, μsp, n_ext, n_med, l, a, z, bessels)
        g1 = G1 / G0
        g2[ind] = 1 + β * abs(g1)^2
    end

    return g2
end

"""
g2_DA_Nlay_cylinder_CW(τ::AbstractVector, data::Nlayer_cylinder_DCS)

Provides a wrapper to `g2_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], β = 1.0, n_ext = 1.0, n_med = 1.0, l = [1.0, 10.0], a = 20.0, z = 0.0, λ = 700.0, bessels = besselroots[1:1000])`.
Allows for arguments to be defined in a structure with `Nlayer_cylinder_DCS`. All vectors must be the same length (you will get a bounds error if not).

# Examples
```
julia> τ = 10 .^(range(-10,stop=0,length=250)) # need to define τ vector first
julia> data = Nlayer_cylinder_DCS(ρ = 1.5, BFi = [3.1e-8, 2.0e-7]) # return ρ = 1.5 and custom BFi
julia> g2_DA_Nlay_cylinder_CW(τ, data) # can then simulate g2 in semiinf
```
"""
function g2_DA_Nlay_cylinder_CW(τ::AbstractVector, data::Nlayer_cylinder_DCS)
    return g2_DA_Nlay_cylinder_CW(τ, data.ρ, data.μa, data.μsp, BFi = data.BFi, n_ext = data.n_ext, n_med = data.n_med, l = data.l, a = data.a, λ = data.λ, bessels = data.bessels)
end
