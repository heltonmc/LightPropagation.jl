
# steady-state solution Kienle 2 layer 1998

@. D = 1/3musp

function fluence_spatial_freq(z, s, μa, μsp, l)

    D = @. 1/3*μsp
    α = @. sqrt((D*s^2 + μa)/D)
    A = 1
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])

    ϕ = BigFloat(0.0)

    if z < l && z >= 0

        ϕ = 1/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))

        ϕ = ϕ*sinh(α[1]*(zb + z0))/(D[1]* α[1])
        ϕ = ϕ*(D[1]*α[1]*cosh(α[1]*(l - z)) + D[2]*α[2]*sinh(α[1]*(l - z)))

        if z < z0
            ϕ  = ϕ - sinh(α[1]*(z0 - z))/(D[1]*α[1])
        end

    else
        ϕ = sinh(α[1]*(zb + z0))*exp(α[2]*(l - z))
        ϕ = ϕ/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))
    end

    return ϕ
end


function fluence_steadystate(ρ, z, μa, μsp, l)


   ϕ = quadgk((s) -> compute_integral(z, s, μa, μsp, l, ρ), 0, 400)
   return ϕ[1]/2π

end



function compute_integral(z, s, μa, μsp, l, ρ)
    return fluence_spatial_freq(z, s, μa, μsp, l)*s*besselj0(ρ*s)
end

function fluence_spatial_freq(z, s, μa, μsp, l)
    # keep l under 10 because of bounds error
    D = @. 1/3μsp
    α = @. sqrt((D*s^2 + μa)/D)
    A = 1
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])

    if 2*cosh(α[1]*(l + zb)) == Inf # cosh(715) returns inf
        ϕ = fluence_spatial_freq_large(z, α, z0, l, D)
    else
        ϕ = fluence_spatial_freq_small(z, α, zb, z0, l, D)
    end
    
    return ϕ
end


function fluence_spatial_freq_small(z, α, zb, z0, l, D)

    if z < l

        ϕ = (cosh(α[1]*(l - z))/(D[2]*α[2]) + sinh(α[1]*(l - z))/(D[1]*α[1]))
        ϕ = ϕ/(cosh(α[1]*(l + zb))/(D[2]*α[2]) + sinh(α[1]*(l + zb))/(D[1]*α[1]))
        ϕ = ϕ*sinh(α[1]*(zb + z0))/(D[1]*α[1])

        if z < z0
            ϕ  = ϕ - sinh(α[1]*(z0 - z))/(D[1]*α[1])
        end

    else
        ϕ = sinh(α[1]*(zb + z0))/(D[1]*α[1]*D[2]*α[2])
        ϕ = ϕ/(cosh(α[1]*(l + zb))/(D[2]*α[2]) + sinh(α[1]*(l + zb))/(D[1]*α[1]))
        ϕ = ϕ*exp(α[2]*(l - z))
    end

    if isinf(ϕ)
        return NaN
    end
    
    return ϕ
end


function fluence_spatial_freq_large(z, α, z0, l, D)

    if z < z0
        ϕ = 0.0

    elseif z < l && z >= z0
        ϕ = exp(α[1]*z0 - α[1]*z)/(2*D[1]*α[1])

    else
        ϕ = exp(α[1]*z0 + α[2]*l - α[2]*z - α[1]*l)/(D[1]*α[1] + D[2]*α[2])
    end

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



function fluence_spatial_freq1(z, s, μa, μsp, l)



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






ϕ = ϕ - exp(α[1]*(z0 - z))*(1 - exp(-2*α[1]*(z0 - z)))/(2*D[1]*α[1])