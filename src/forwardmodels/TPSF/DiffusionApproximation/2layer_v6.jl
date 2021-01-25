function fluence_spatial_complexfreq(z, s, μa, μsp, l, f, nmed)

    #get numerical instabillity at large s
    D = @. 1/(3(μsp + μa))
    ν = 29.9792345/nmed
    ω = 2*π*f
    α = @. sqrt((D*s^2 + μa + ω*im/ν)/D)
    A = 1.0
    zb = 2*A*D[1]
    z0 = 1/(μsp[1] + μa[1])
    
    ## equation 11
    if z < z0

        # have to account for exp blowing up to inf
        if real(α[1]*z0 - α[1]*z) > 15
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


  
x, w = gausslegendre(N)
x = x.*ub/2 .+ ub/2
w = w.*ub/2

function FD_test(ρ, z, μa, μsp, l, f, nmed, x, w)
	
	FD(s) = compute_complexintegral(ρ, z, s, μa, μsp, l, f, nmed)

	return sum(FD.(x) .* w)/2pi
end

function compute_TDML(t, ρ, μa, μsp, l, ub, N)
	
	Rt = zeros(Float64, length(t))
	freqcomp = zeros(Float64, N)
	x, w = gausslegendre(N)
	x = x.*ub/2 .+ ub/2
	w = w.*ub/2
    
    f = 1
	
	FD(s, f) = compute_complexintegral(ρ, 0.0, s, μa, μsp, l, f, 1.0)
	Threads.@threads for m in eachindex(f)
		freqcomp[m] = FD_test(ρ, 0.0, μa, μsp, l, f[m], 1.0, x, w)
    end
    
	return freqcomp
end	