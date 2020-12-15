
# equations 11-13 for spatial frequencies converted as cosh = e^x - e^-x /2
function fluence_spatial_freq(z, s, μa, μsp, l)



    D = @. 1/3μsp
    α = @. sqrt((D*s^2 + μa)/D)
    A = 1
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])



    if z < l

        ϕ = D[1]*α[1]*(1 + exp(-2*α[1]*(l - z))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l - z)))
        ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
        ϕ = ϕ*(1 + exp(-2*α[1]*(zb + z0)))/(2*D[1]*α[1])
        ϕ = ϕ*exp(α[1]*z0 - α[1]*z)

        if z < z0
            ϕ = ϕ - exp(α[1]*(z0 - z))*(1 - exp(-2*α[1]*(z0 - z)))/(2*D[1]*α[1])
        end

    else

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


    ϕ = quadgk((s) -> compute_integral(ρ, z, s, μa, μsp, l), 0, 10000)
    return ϕ[1]/2π
 
 end
 
 
function ss_compare(ρ, z, μa, μsp)

    D = 1/3μsp
    A = 1
    z0 = 1/(μsp + μa)
    zb = 2*A*D
    μe = sqrt(3*μa*(μa + μsp))

    ϕ = exp(-μe*sqrt((z - z0)^2 + ρ^2))/sqrt((z - z0)^2 + ρ^2)
    ϕ = ϕ - exp(-μe*sqrt((z + z0 + 2*zb)^2 + ρ^2))/sqrt((z + z0 + 2*zb)^2 + ρ^2)
    ϕ = ϕ/(4*π*D)
end





ϕ = sinh(α[1]*(zb + z0))*exp(α[2]*(l - z))
ϕ = ϕ/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))




ϕ = (D[1]*α[1])*cosh(α[1]*(l - z)) + D[2]*α[2]*sinh(α[1]*(l - z))/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))

ϕ = ϕ*sinh(α[1]*(zb + z0))/(D[1]*α[1])
