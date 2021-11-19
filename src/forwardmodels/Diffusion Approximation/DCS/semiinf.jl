"""
    g1_DA_semiinf_CW(τ::AbstractVector, ρ, μa, μsp; BFi = 2e-8, β = 1.0, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0)

Compute the electric field autocorrelation function function g1 in a semi-infinite geometry (CW).

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
- `z`: the z-depth orthogonal from the boundary (cm)
- `λ`: wavelength (nm)
"""
function g1_DA_semiinf_CW(τ::AbstractVector, ρ, μa, μsp; BFi = 2e-8, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0)
    g1 = similar(τ)
    k0 = 2 * n_med * π / (λ * 1e-7) # this should be in cm
    tmp = 2 * μsp * BFi * k0^2

    G0 = fluence_DA_semiinf_CW(ρ, μa, μsp, n_ext = n_ext, n_med = n_med, z = z)

    for ind in eachindex(τ)
        μa_dynamic = μa + tmp * τ[ind]
        G1 = fluence_DA_semiinf_CW(ρ, μa_dynamic, μsp, n_ext = n_ext, n_med = n_med, z = z)
        g1[ind] = G1 / G0
    end

    return g1
end

"""
    g1_DA_semiinf_TD(τ::AbstractVector, t::AbstractFloat, ρ, μa, μsp; BFi = 2e-8, β = 1.0, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0)
    g1_DA_semiinf_TD(τ::AbstractVector, t::Abstractector, ρ, μa, μsp; BFi = 2e-8, β = 1.0, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0)

Compute the time-domain electric field autocorrelation function function g1 in a semi-infinite geometry.

If t is a AbstractVector then g1(tau) will be calculated by ∫P(t)g1(tau, t)dt with integration bounds given by t[1] and t[end].
This is needed for time-gated applications.

# Arguments
- `τ`: time lags (s)
- `t`: photon time of flight (related to photon path length)
- `ρ`: the source detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Keyword arguments
- `BFi`: Blood flow index ~αDb (cm²/s)
- `β`: constant in Siegert relation dependent on collection optics
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction
- `z`: the z-depth orthogonal from the boundary (cm)
- `λ`: wavelength (nm)
- `N_quad`: number of nodes in gauss-legendre integration
"""
function g1_DA_semiinf_TD(τ::AbstractVector, t::AbstractFloat, ρ, μa, μsp; BFi = 2e-8, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0)
    g1 = similar(τ)
    k0 = 2 * n_med * π / (λ * 1e-7) # this should be in cm
    tmp = 2 * μsp * BFi * k0^2

    G0 = fluence_DA_semiinf_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, z = z)

    for ind in eachindex(τ)
        μa_dynamic = μa + tmp * τ[ind]
        G1 = fluence_DA_semiinf_TD(t, ρ, μa_dynamic, μsp, n_ext = n_ext, n_med = n_med, z = z)
        g1[ind] = G1 / G0
    end

    return g1
end
function g1_DA_semiinf_TD(τ::AbstractVector, t::AbstractVector, ρ, μa, μsp; BFi = 2e-8, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0, N_quad = 50)
    g1 = similar(τ)
    k0 = 2 * n_med * π / (λ * 1e-7) # this should be in cm
    tmp = 2 * μsp * BFi * k0^2
    x, w = gausslegendre(N_quad)

    tpsf_norm = integrate(t -> fluence_DA_semiinf_TD(t, ρ, μa, μsp, n_ext = n_ext, n_med = n_med, z = z), x, w, t[1], t[2])

    Threads.@threads for ind in eachindex(τ)
        μa_dynamic = μa + tmp * τ[ind]
        g1[ind] = integrate(t -> (fluence_DA_semiinf_TD(t, ρ, μa_dynamic, μsp, n_ext = n_ext, n_med = n_med, z = z) / tpsf_norm), x, w, t[1], t[2])
    end

    return g1
end

function g2_DA_semiinf_CW(τ::AbstractVector, ρ, μa, μsp; BFi = 2e-8, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0, β = 1.0)
    g1 = g1_DA_semiinf_CW(τ, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
    for ind in eachindex(τ)
        g1[ind] = 1 + β * abs(g1[ind])^2
    end
    return g1
end

function g2_DA_semiinf_TD(τ::AbstractVector, t::Union{AbstractVector, AbstractFloat}, ρ, μa, μsp; BFi = 2e-8, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700.0, β = 1.0, N_quad = 50)
    if t isa AbstractVector
        g1 = g1_DA_semiinf_TD(τ, t, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ, N_quad = N_quad)
    elseif t isa AbstractFloat
        g1 = g1_DA_semiinf_TD(τ, t, ρ, μa, μsp, BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
    end

    for ind in eachindex(τ)
        g1[ind] = 1 + β * abs(g1[ind])^2
    end
    return g1
end

function integrate(f::Function, x, w, lb, ub)
    out = zero(eltype(x))
    for m in eachindex(x)
        out += f(x[m] * (ub - lb) / 2 + (lb + ub) / 2) * w[m]
    end
    return out * (ub - lb) / 2
end

# use this function when it is more efficient to compute y points all together for layered media (Laplace contour is fixed for each time point)
function integrate_compute_y_first(f::Function, x, w, lb, ub)
    x = @. x * (ub - lb) / 2 + (lb + ub) / 2
    y = f(x) # the function should be able to take a vector and return a vector
    out = zero(eltype(x))
    for m in eachindex(x)
        out += y[m] * w[m]
    end
    return out * (ub - lb) / 2
end
