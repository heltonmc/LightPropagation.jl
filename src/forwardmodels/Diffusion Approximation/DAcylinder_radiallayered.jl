

dIm(m, γ, ρ) = γ*besseli(m + 1, γ*ρ) + m*besseli(m, γ*ρ)/ρ
dKm(m, γ, ρ) = -γ*besselk(m + 1, γ*ρ) + m*besselk(m, γ*ρ)/ρ

function _get_αβ_start(m, γ, D, ρ, n, ind)

    α = D[ind - 1] * (n[ind - 1])^2 * dKm(m, γ[ind - 1], ρ[ind]) * besseli(m, γ[ind]*ρ[ind])
    α -= D[ind] * (n[ind])^2 * besselk(m, γ[ind - 1]*ρ[ind]) * dIm(m, γ[ind], ρ[ind]) 

    β = D[ind - 1] * (n[ind - 1])^2 * dIm(m, γ[ind - 1], ρ[ind]) * besseli(m, γ[ind]*ρ[ind])
    β -= D[ind] * (n[ind])^2 * besseli(m, γ[ind - 1]*ρ[ind]) * dIm(m, γ[ind], ρ[ind]) 

    return α, β
end

function _get_αn_βn(m, αnp1, βnp1, γ, D, ρ, n, ind)

    αn = D[ind] * (n[ind])^2 * dKm(m, γ[ind], ρ[ind + 1])
    αn *= (αnp1*besseli(m, γ[ind + 1]*ρ[ind + 1]) - βnp1*besselk(m, γ[ind + 1]*ρ[ind + 1]))
    αn_aux = D[ind + 1] * (n[ind + 1])^2 * besselk(m, γ[ind]*ρ[ind + 1])
    αn_aux *= (αnp1*dIm(m, γ[ind + 1], ρ[ind + 1]) - βnp1*dKm(m, γ[ind + 1], ρ[ind + 1]))
    αn -= αn_aux

    βn = D[ind] * (n[ind])^2 * dIm(m, γ[ind], ρ[ind + 1])
    βn *= (αnp1*besseli(m, γ[ind + 1]*ρ[ind + 1]) - βnp1*besselk(m, γ[ind + 1]*ρ[ind + 1]))
    βn_aux = D[ind + 1] * (n[ind + 1])^2 * besseli(m, γ[ind]*ρ[ind + 1])
    βn_aux *= (αnp1*dIm(m, γ[ind + 1], ρ[ind + 1]) - βnp1*dKm(m, γ[ind + 1], ρ[ind + 1]))
    βn -= βn_aux

    return αn, βn
end


function _get_α1_β1(m, γ, D, ρ, n)

    ind = length(γ)
    α, β = _get_αβ_start(m, γ, D, ρ, n, ind)
    ind -= 1
    for i in ind:-1:2
        println(α, β)
        α, β = _get_αn_βn(m, α, β, γ, D, ρ, n, i)
    end

    return α, β
end


function _get_A1_B1(m, γ, D, ρ, n, ρ0)

    α1, β1 = _get_α1_β1(m, γ, D, ρ, n)

    A1 = besselk(m, γ[1]*ρ[1])
    A1 *= (β1*besselk(m, γ[1]*ρ0) - α1*besseli(m, γ[1]*ρ0))
    A1 /= D[1]*(α1*besseli(m, γ[1]*ρ[1]) - β1*besselk(m, γ[1]*ρ[1]))

    B1 = β1/D[1]
    B1 /= (α1*besseli(m, γ[1]*ρ[1]) - β1*besselk(m, γ[1]*ρ[1]))
    B1 *= (besseli(m, γ[1]*ρ0)*besselk(m, γ[1]*ρ[1]) - besselk(m, γ[1]*ρ0)*besseli(m, γ[1]*ρ[1]))

    return A1, B1
end


function _integrand_modbessel(m, A1, B1, γ, ρ, ϕ)

    out = A1*besseli(m, γ[1]*ρ[1]) + B1*besselk(m, γ[1]*ρ[1])
    out *= cos(m*ϕ)

    return out
end

function sum_modbessel(k, μa, D, ρ, n, ρ0, ϕ)

    a = 0.0
    A1 = 0.0
    B1 = 0.0
    γ = @. sqrt(μa/D + k^2)

    bessel_index = -100:100
    for m in bessel_index
        A1, B1 = _get_A1_B1(m, γ, D, ρ, n, ρ0)
        a += _integrand_modbessel(m, A1, B1, γ, ρ, ϕ)
    end

    return a
end


function _greens_homogenous(μa, D, ρ, n, ρ0, ϕ, z)

    g = quadgk(k -> cos(k*z)*sum_modbessel(k, μa, D, ρ, n, ρ0, ϕ), 0, Inf)

    return g[1]/(2*π^2)
end
    

function _green_Nradialcyl(μa, D, ρ, n, ρ0, ϕ, z)

    r = sqrt(ρ[1]^2 + ρ0^2 - 2*ρ[1]*ρ0*cosh(ϕ) + z^2)

    g = exp(-r*sqrt(μa[1]/D[1]))/(4*π*D[1]*r)
    gh = _greens_homogenous(μa, D, ρ, n, ρ0, ϕ, z)

    return g + gh
end