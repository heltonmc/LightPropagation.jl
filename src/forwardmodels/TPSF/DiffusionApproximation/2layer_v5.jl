using SpecialFunctions, QuadGK, ForwardDiff


##### Diffusion approximation for 2-layered medium with matched index of refraction #####


"""
    fluence_spatial_freq(z, s, μa, μsp, l)

Compute the spatial frequency components of the fluence for 2-layered model 

# Arguments
- `z`: depth (z=0 at the surface, z=l at interface of two layers). 
- `s`: spatial frequency component s^2 = sx^2 + sy^2
- `μa`: absorption coefficient of both layers [μa1, μa2]
- `μsp`: reduced scattering coefficient of both layers [μsp1, μsp2]
- `l`: thickness of first layer


# Examples
```jldoctest
julia> fluence_spatial_freq(0., 1,  [0.005,0.1], [10,100.], 1)
0.8647044203305079
```
"""
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
            ϕ = ϕ*(1 - exp(-2*α[1]*(zb + z0)))/(2*D[1]*α[1])
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

function fluence_spatial_complexfreq(z, s, μa, μsp, l, f, nmed)

    #get numerical instabillity at large s
    D = @. 1/(3(μsp + μa))
    ν = 29.9792345/nmed
    ω = 2*π*f
    α = @. sqrt((D*s^2 + μa + ω*im/ν)/D)
    A = 0.493
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])
    
    ## equation 11
    if z < z0

        # have to account for exp blowing up to inf
        if abs2(α[1]*z0 - α[1]*z) > 15
            ϕ = exp(α[1]*(z0 - 2*l)) - exp(α[1]*(-z0 - 2*l - zb)) - exp(α[1]*(z0 - 2*l - 2*zb)) + exp(α[1]*(-z0 - 2*l - 4*zb))
            ϕ = ϕ*(D[1]*α[1] - D[2]*α[2])/(D[1]*α[1] + D[2]*α[2])
            ϕ = ϕ + exp(-α[1]*z0) - exp(-α[1]*(2*zb + z0))
            ϕ = ϕ/(2*D[1]*α[1])


        else
            
            ϕ = D[1]*α[1]*(1 + exp(-2*α[1]*(l - z))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l - z)))
            ϕ = ϕ/(D[1]*α[1]*(1 + exp(-2*α[1]*(l + zb))) +  D[2]*α[2]*(1 - exp(-2*α[1]*(l + zb))))
            ϕ = ϕ*(1 - exp(-2*α[1]*(zb + z0)))/(2*D[1]*α[1])
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

function compute_complexintegral(ρ, z, s, μa, μsp, l, f, nmed)
    return fluence_spatial_complexfreq(z, s, μa, μsp, l, f, nmed)*s*besselj0(ρ*s)
end

function fluence_FD(ρ, z, μa, μsp, l, f, nmed)
    return (quadgk((s) -> compute_complexintegral(ρ, z, s, μa, μsp, l, f, nmed), 0, Inf, rtol=1e-6)[1])/2π 
end

function fluence_TD(ρ, z, μa, μsp, l, nmed)
    fr = 0.19531 .* Vector(1:512)
	Rt = Array{Complex{Float64},1}(undef, length(fr))
    map!((f) -> fluence_FD(ρ, z, μa, μsp, l, f, nmed), Rt, fr) 
end


# compute equation 14 integrand
function compute_integral(ρ, z, s, μa, μsp, l)
    return fluence_spatial_freq(z, s, μa, μsp, l)*s*besselj0(ρ*s)
end


# numerically integrate equation 14
function fluence_steadystate(ρ, z, μa, μsp, l)
    return (quadgk((s) -> compute_integral(ρ, z, s, μa, μsp, l), 0, Inf, rtol=1e-6)[1])/2π 
 end

# Eqn 16 to compute reflectance (use auto diff)
function refl_steadystate(ρ, μa, μsp, l)
    f(z) = fluence_steadystate(ρ, z, μa, μsp, l)
    g(z) = ForwardDiff.derivative(z -> f(z), z)
    return 0.118*fluence_steadystate(ρ, 0, μa, μsp, l) + 0.306/(3*(μa[1] + μsp[1]))*g(0)
end

# Eqn 3 Kienle 97 for spatial fluence semi-infinite (use as comparison for l -> ∞)
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


##reproduce Figure 2 Kienle 98 *seem to be slightly off by some factor
a = map((ρ) -> refl_steadystate(ρ, [0.02,0.01], [1.3,1.2], 2), 0.5:0.5:20)
b = map((ρ) -> refl_steadystate(ρ, [0.02,0.01], [1.3,0.7], 2), 0.5:0.5:20)

plot(0.5:0.5:20, a, yscale=:log10, lw = 2, label = "crosses", xlabel = "distance [mm]", ylabel = "reflectance")
plot!(0.5:0.5:20, b, yscale=:log10, lw = 2, label = "circles")

##reproduce Figure 3 Kienle 98
b = map((ρ) -> refl_steadystate(ρ, [0.005,0.01], [1.3,1.0], 6), 0.5:0.5:20)
a = map((ρ) -> refl_steadystate(ρ, [0.005,0.022], [1.3,1.0], 6), 0.5:0.5:20)

plot(0.5:0.5:20, a, yscale=:log10, lw = 2, label = "crosses", xlabel = "distance [mm]", ylabel = "reflectance")
plot!(0.5:0.5:20, b, yscale=:log10, lw = 2, label = "circles")


##compare to semi-infinite model from Kienle 97 as l gets large
a = map((ρ) -> fluence_steadystate(ρ, 0, [0.02,0.02], [1.3,1.3], 10), 0.5:0.5:20)
b = map((ρ) -> ss_compare(ρ, 0, 0.02, 1.3), 0.5:0.5:20)

plot(0.5:0.5:20, a, yscale=:log10, lw = 2, label = "2-layer", xlabel = "distance [mm]", ylabel = "reflectance")
plot!(0.5:0.5:20, b, yscale=:log10, lw = 2, label = "semi-inf", alpha = 0.7)

