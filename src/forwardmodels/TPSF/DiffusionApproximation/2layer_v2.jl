
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
            ϕ = ϕ - ((exp(α[1]*(z0 - z))*(1 - exp(-2*α[1]*(z0 - z))))/(2*D[1]*α[1]))
        end

    else

        ϕ = exp(α[1]*(z0 - l) + α[2]*(l - z))
        ϕ = ϕ*(1 - exp(-2*α[1]*(zb + z0)))
        ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) + D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
    end

    if isnan(ϕ) || isinf(ϕ)
        ϕ = 0
    end

    return ϕ
end

# compute equation 14 integrand
function compute_integral(ρ, z, s, μa, μsp, l)
    return fluence_spatial_freq(z, s, μa, μsp, l)*s*besselj0(ρ*s)
end



# numerically integrate equation 14
function fluence_steadystate(ρ, z, μa, μsp, l)


    ϕ = quadgk((s) -> compute_integral(ρ, z, s, μa, μsp, l), 0, Inf, rtol=1e-8)
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




function fluence_spatial_freq(z, s, μa, μsp, l)



    D = @. 1/3μsp
    α = @. sqrt((D*s^2 + μa)/D)
    A = 1
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])
    tolr = 1e-50
    ϕ = []

    ## equation 11
    if z < z0

        # limit as neg exp gets large
        if exp(-2*α[1]*z0) < tolr
            ϕ = 0.0
        elseif exp(-2*l) < tolr
            ϕ = exp(α[1]*z0 - α[1]*z)/(2*D[1]*α[1])
            ϕ = ϕ*(exp(-2*α[1]*(z0 - z)) - exp(-2*α[1]*(zb + z0)))
        elseif exp(-2*l*α[1]) < tolr
            ϕ = 0.0
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
    elseif z > l

        ϕ = exp(α[1]*(z0 - l) + α[2]*(l - z))
        ϕ = ϕ*(1 - exp(-2*α[1]*(zb + z0)))
        ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) + D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
    end


    return ϕ
end





function fluence_spatial_freq(z, s, μa, μsp, l)
    #get numerical instabillity at large s 


    D = @. 1/(3(μsp + μa))
    α = @. sqrt((D*s^2 + μa)/D)
    A = 0.493
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])
    tolr = 1e-50
    ϕ = zeros(Float64, 1)

    ## equation 11
    if z < z0

        # have to account for exp blowing up to inf
        if exp(α[1]*z0 - α[1]*z) > 1e8
            ϕ = 0.0
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



function manual_integral(ρ, z, μa, μsp, l)

    xs = 0:0.001:300
    a = map((s) -> fluence_spatial_freq(z, s, μa, μsp, l), xs)
    b = a.*xs
    c  = b.*[besselj0(ρ*x) for x in xs]
    return c
end


# compute equation 14 integrand
function compute_integral(ρ, z, s, μa, μsp, l)
    return fluence_spatial_freq(z, s, μa, μsp, l)*s*besselj0(ρ*s)
end



# numerically integrate equation 14
function fluence_steadystate(ρ, z, μa, μsp, l)


    ϕ = quadgk((s) -> compute_integral(ρ, z, s, μa, μsp, l), 0, Inf, rtol=1e-8)
    return ϕ[1]/2π
 
 end