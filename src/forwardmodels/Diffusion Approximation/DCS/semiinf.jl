@with_kw struct DAsemiinf_DCS{T <: Real}
    μsp::T = 10.0                                       # reduced scattering coefficient (1/cm)
    μa::T = 0.1                                         # absorption coefficient (1/cm)
    n_ext::T = 1.0                                      # surrounding index of refraction
    n_med::T = 1.0                                      # layers index of refraction

    ρ::Union{T, AbstractVector{T}} = 1.0                # source-detector separation (cm)
    z::T = 0.0                                          # detector depth (cm)

    β::T = 1.0                                          # constant in Siegert relation dependent on collection optics
    BFi::T = 2.0e-8                                     # Blood flow index ~αDb (cm²/s)
    λ::T = 750.0                                        # wavelength (nm)
end

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


function integrate(f, x, w, lb, ub)
    out = zero(eltype(x))
    for m in eachindex(x)
        out += f(x[m] * (ub - lb) / 2 + (lb + ub) / 2) * w[m]
    end
    return out * (ub - lb) / 2
end
