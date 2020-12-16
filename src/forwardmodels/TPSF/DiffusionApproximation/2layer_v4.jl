using SpecialFunctions, QuadGK, ForwardDiff


function fluence_spatial_freq(z, s, μa, μsp, l)

    #get numerical instabillity at large s
    D = @. 1/(3(μsp + μa))
    α = @. sqrt((D*s^2 + μa)/D)
    A = 0.493
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])
    
    ## equation 11
    if z < z0

        # have to account for exp blowing up to inf
        if (α[1]*z0 - α[1]*z) > 15
            ϕ = exp(α[1]*(z0 - 2*l)) - exp(α[1]*(-z0 - 2*l - zb)) - exp(α[1]*(z0 - 2*l - 2*zb)) + exp(α[1]*(-z0 - 2*l - 4*zb))
            ϕ = ϕ*(D[1]*α[1] - D[2]*α[2])/(D[1]*α[1] + D[2]*α[2])
            ϕ = ϕ + exp(-α[1]*z0) - exp(-α[1]*(2*zb + z0))
            ϕ = ϕ/(2*D[1]*α[1])


        else
            ϕ = D[1]*α[1]*(1 + exp(-2*α[1]*(l - z))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l - z)))
            ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
            ϕ = ϕ*(1 + exp(-2*α[1]*(zb + z0)))/(2*D[1]*α[1])
            ϕ = ϕ*exp(α[1]*z0 - α[1]*z)
            ϕ = ϕ - ((exp(α[1]*(z0 - z))*(1 - exp(-2*α[1]*(z0 - z))))/(2*D[1]*α[1]))
        end

    ## equation 12
    elseif z > z0 && z < l
        ϕ = D[1]*α[1]*(1 + exp(-2*α[1]*(l - z))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l - z)))
        ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
        ϕ = ϕ*(1 + exp(-2*α[1]*(zb + z0)))/(2*D[1]*α[1])
        ϕ = ϕ*exp(α[1]*z0 - α[1]*z)


    #equation 13
    elseif z >= l

        ϕ = exp(α[1]*(z0 - l) + α[2]*(l - z))
        ϕ = ϕ*(1 - exp(-2*α[1]*(zb + z0)))
        ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) + D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
    end


    return ϕ
end

# compute equation 14 integrand
function compute_integral(ρ, z, s, μa, μsp, l)
    return fluence_spatial_freq(z, s, μa, μsp, l)*s*besselj0(ρ*s)
end


# numerically integrate equation 14
function fluence_steadystate(ρ, z, μa, μsp, l)
    return quadgk((s) -> compute_integral(ρ, z, s, μa, μsp, l), 0, Inf, rtol=1e-8)[1]/2π 
 end

# Eqn 16 to compute reflectance (use auto diff)
function refl_steadystate(ρ, μa, μsp, l)
    f(z) = fluence_steadystate(ρ, z, μa, μsp, l)
    g(z) = ForwardDiff.derivative(z -> f(z), z)
    return 0.118*fluence_steadystate(ρ, 0, μa, μsp, l) + 0.306/(3*(μa[1] + μsp[1]))*g(0)
end

#compare semi-infinite fluence
function ss_compare(ρ, z, μa, μsp)
    D = 1/(3(μsp + μa))
    A = 0.493
    z0 = 1/(μsp + μa)
    zb = 2*A*D
    μe = sqrt(3*μa*(μa + μsp))

    ϕ = exp(-μe*sqrt((z - z0)^2 + ρ^2))/sqrt((z - z0)^2 + ρ^2)
    ϕ = ϕ - exp(-μe*sqrt((z + z0 + 2*zb)^2 + ρ^2))/sqrt((z + z0 + 2*zb)^2 + ρ^2)
    ϕ = ϕ/(4*π*D)
end


function refl_steadystate1(ρ, μa, μsp, l)
    f(z) = ss_compare(ρ, z, μa, μsp)
    g(z) = ForwardDiff.derivative(z -> f(z), z)
    return 0.118*ss_compare(ρ, 0, μa, μsp) + 0.306/(3*(μa[1] + μsp[1]))*g(0)
end