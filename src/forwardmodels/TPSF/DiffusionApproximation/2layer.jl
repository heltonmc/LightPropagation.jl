
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
    D = @. 1/3*μsp
    α = @. sqrt((D*s^2 + μa)/D)
    A = 1
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])

    if isnan(fluence_spatial_freq_small(z, α, zb, z0, l, D))
        ϕ = fluence_spatial_freq_large(z, α, z0, l, D)
    
    else
        ϕ = fluence_spatial_freq_small(z, α, zb, z0, l, D)
    end

end


function fluence_spatial_freq_small(z, α, zb, z0, l, D)

    if z < l && z >= 0

        ϕ = sinh(α[1]*(zb + z0))/(D[1]*α[1])
        ϕ = ϕ*(D[1]*α[1]*cosh(α[1]*(l - z)) + D[2]*α[2]*sinh(α[1]*(l - z)))
        ϕ = ϕ/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))

        if z < z0
            ϕ  = ϕ - sinh(α[1]*(z0 - z))/(D[1]*α[1])
        end

    else
        ϕ = sinh(α[1]*(zb + z0))*exp(α[2]*(l - z))
        ϕ = ϕ/(D[1]*α[1]*cosh(α[1]*(l + zb)) + D[2]*α[2]*sinh(α[1]*(l + zb)))
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
