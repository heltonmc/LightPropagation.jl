using SpecialFunctions: besselix, besselkx
using ArbNumerics: besseli, besselk, ArbFloat


dImx(m, γ, ρ) = γ*besselix(m + 1, γ*ρ) + m*besselix(m, γ*ρ)/ρ
dKmx(m, γ, ρ) = -γ*besselkx(m + 1, γ*ρ) + m*besselkx(m, γ*ρ)/ρ

dIm(m, γ, ρ) = γ*besseli(m + 1, γ*ρ) + m*besseli(m, γ*ρ)/ρ
dKm(m, γ, ρ) = -γ*besselk(m + 1, γ*ρ) + m*besselk(m, γ*ρ)/ρ

function _get_αβ_start(m, γ, D, ρ, n, ind)

    α = D[ind - 1] * (n[ind - 1])^2 * dKm(m, γ[ind - 1], ρ[ind]) * besseli(m, γ[ind]*ρ[ind])
    α -= D[ind] * (n[ind])^2 * besselk(m, γ[ind - 1]*ρ[ind]) * dIm(m, γ[ind], ρ[ind]) 

    β = D[ind - 1] * (n[ind - 1])^2 * dIm(m, γ[ind - 1], ρ[ind]) * besseli(m, γ[ind]*ρ[ind])
    β -= D[ind] * (n[ind])^2 * besseli(m, γ[ind - 1]*ρ[ind]) * dIm(m, γ[ind], ρ[ind]) 

    return α, β
end

##reduced form 
function _get_αβ_startx(m, γ, D, ρ, n, ind)

    α = D[ind - 1] * (n[ind - 1])^2 * dKmx(m, γ[ind - 1], ρ[ind]) * besselix(m, γ[ind]*ρ[ind])
    α -= D[ind] * (n[ind])^2 * besselkx(m, γ[ind - 1]*ρ[ind]) * dImx(m, γ[ind], ρ[ind]) 
    α *= exp(ρ[ind]*(γ[ind] - γ[ind - 1]))

    β = D[ind - 1] * (n[ind - 1])^2 * dImx(m, γ[ind - 1], ρ[ind]) * besselix(m, γ[ind]*ρ[ind])
    β -= D[ind] * (n[ind])^2 * besselix(m, γ[ind - 1]*ρ[ind]) * dImx(m, γ[ind], ρ[ind]) 
    β *= exp(ρ[ind]*(γ[ind] - γ[ind - 1]))

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
##reduced form
function _get_αn_βnx(m, αnp1, βnp1, γ, D, ρ, n, ind)

    αn = D[ind] * (n[ind])^2 * dKmx(m, γ[ind], ρ[ind + 1])
    αn *= (αnp1*besselix(m, γ[ind + 1]*ρ[ind + 1]) - βnp1*exp(-2*ρ[ind + 1]*γ[ind + 1])*besselkx(m, γ[ind + 1]*ρ[ind + 1]))
    αn_aux = D[ind + 1] * (n[ind + 1])^2 * besselkx(m, γ[ind]*ρ[ind + 1])
    αn_aux *= (αnp1*dImx(m, γ[ind + 1], ρ[ind + 1]) - βnp1*exp(-2*ρ[ind + 1]*γ[ind + 1])*dKmx(m, γ[ind + 1], ρ[ind + 1]))
    αn -= αn_aux
    αn *= exp(ρ[ind + 1]*(γ[ind + 1] - γ[ind]))

    βn = D[ind] * (n[ind])^2 * dImx(m, γ[ind], ρ[ind + 1])
    βn *= (αnp1*besselix(m, γ[ind + 1]*ρ[ind + 1]) - βnp1*exp(-2*ρ[ind + 1]*γ[ind + 1])*besselkx(m, γ[ind + 1]*ρ[ind + 1]))
    βn_aux = D[ind + 1] * (n[ind + 1])^2 * besselix(m, γ[ind]*ρ[ind + 1])
    βn_aux *= (αnp1*dImx(m, γ[ind + 1], ρ[ind + 1]) - βnp1*exp(-2*ρ[ind + 1]*γ[ind + 1])*dKmx(m, γ[ind + 1], ρ[ind + 1]))
    βn -= βn_aux
    βn *= exp(ρ[ind + 1]*(γ[ind + 1] + γ[ind]))


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
##reduced form
function _get_α1_β1x(m, γ, D, ρ, n)

    ind = length(γ)
    α, β = _get_αβ_startx(m, γ, D, ρ, n, ind)
    ind -= 1
    for i in ind:-1:2
        println(α, β)
        α, β = _get_αn_βnx(m, α, β, γ, D, ρ, n, i)
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

##reduced form (calculates A1 and B1 together)
function _integrand_modbesselx(m, γ, D, ρ, n, ρ0, ϕ)

    α1, β1 = _get_α1_β1x(m, γ, D, ρ, n)

    A1 = (β1*besselkx(m, γ[1]*ρ0)*exp(-2*γ[1]*ρ0) - α1*besselix(m, γ[1]*ρ0))
    A1 /= (α1*besselix(m, γ[1]*ρ[1]) - β1*besselkx(m, γ[1]*ρ[1])*exp(-2*γ[1]*ρ[1]))
    A1 *= besselkx(m, γ[1]*ρ[1])*besselix(m, γ[1]*ρ[1])/D[1]
    A1 *= exp(γ[1]*(ρ0 - ρ[1]))

    B1 = (besselix(m, γ[1]*ρ0)*besselkx(m, γ[1]*ρ[1])*exp(γ[1]*(ρ0 - ρ[1])) - besselkx(m, γ[1]*ρ0)*besselix(m, γ[1]*ρ[1])*exp(γ[1]*(-ρ0 + ρ[1])))
    B1 /= (α1*besselix(m, γ[1]*ρ[1]) - β1*besselkx(m, γ[1]*ρ[1])*exp(-2*γ[1]*ρ[1]))
    B1 *=  besselkx(m, γ[1]*ρ[1])*β1/D[1]
    B1 *= exp(-2*γ[1]*ρ[1])

    return A1*cos(m*ϕ)
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

##reduced form
function sum_modbesselx(k, μa, D, ρ, n, ρ0, ϕ)
    
    a = 0.0
    A1 = 0.0
    B1 = 0.0
    γ = @. sqrt(μa/D + k^2)

    bessel_index = 0:200
    for m in bessel_index
        a += (2 - (m == 0 ? 1 : 0)) * _integrand_modbesselx(m, γ, D, ρ, n, ρ0, ϕ)
    end

    return a
end
## switch to arbitrary precision
function sum_modbessel(k, μa, D, ρ, n, ρ0, ϕ)
    a = 0.0
    A1 = 0.0
    B1 = 0.0
    γ = @. sqrt(μa/D + k^2)
    γ1 = @. sqrt(ArbFloat.(μa)/D + k^2)
    if k < 10.0
        bessel_index = 0:500
    else
        bessel_index = 0:500
    end

    for m in bessel_index

        if besselix(m,minimum(γ)*minimum([ρ; ρ0])) < floatmin(eltype(μa))*50000 # besselkx(m,maximum(γ)*maximum([ρ; ρ0])) > typemax(eltype(μa))/2 #besselix(m, minimum(γ)*minimum([ρ; ρ0])) == 0
            #println(m); println("breaking")
            #break
            A1, B1 = _get_A1_B1(m, γ1, D, ρ, n, ρ0)
            #println(besseli(m, γ1[2]*ρ[2]))
            a += convert(eltype(μa), (2 - (m == 0 ? 1 : 0)) * _integrand_modbessel(m, A1, B1, γ1, ρ, ϕ))
            #println(a)
        else
            besselix(m,minimum(γ)*minimum([ρ; ρ0]))
            a += (2 - (m == 0 ? 1 : 0)) * _integrand_modbesselx(m, γ, D, ρ, n, ρ0, ϕ)
           # println(a)
        end

    end

    return a
end

function _greens_homogenous(μa, D, ρ, n, ρ0, ϕ, z)

    g = quadgk(k -> cos(k*z)*sum_modbessel(k, μa, D, ρ, n, ρ0, ϕ), 0, 400; maxevals=10^3)

    return g[1]/(2*π^2)
end
    

function _green_Nradialcyl(μa, D, ρ, n, ρ0, ϕ, z)

    r = sqrt(ρ[1]^2 + ρ0^2 - 2*ρ[1]*ρ0*cosh(ϕ) + z^2)

    g = exp(-r*sqrt(μa[1]/D[1]))/(4*π*D[1]*r)
    println(g)
    gh = _greens_homogenous(μa, D, ρ, n, ρ0, ϕ, z)

    return g + gh
end

function _greens_homogenous(μa, D, ρ, n, ρ0, ϕ, z)
    
    N = 190
    x, w = gausslegendre(N)
    ub = 150

    green = zeros(BigFloat, N)

    Threads.@threads for m in eachindex(x)
        green[m] = cos(x[m]*ub/2 +ub/2*z)*sum_modbessel(x[m]*ub/2 +ub/2, μa, D, ρ, n, ρ0, ϕ)*w[m]*ub
    end

    return sum(green)/(2*pi^2)
end


