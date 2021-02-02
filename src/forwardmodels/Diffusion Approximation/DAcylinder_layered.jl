
function greenNcylin(μsp, μa, sn, lz)
    A = 1
    D = 1/3μsp

    α = sqrt(μa/D + sn^2)


    zb = 2*A*D
    z0 = 1/μsp    
    
    
    g = (exp(-α[1]*z0) - exp(-α[1]*(z0 + 2*zb[1]))) / (2*D*α)

    g1 = exp(α[1]*(z0 - 2*l[1]))*(1 - exp(-2*α[1]*(z0 + zb[1])))*(1 - exp(-2*α[1]*zb[1]))
    g1 *= (D[1]*α[1]*n[1]^2*β - D[2]*α[2]*n[2]^2*γ)
    g1 /= D[1]*α[1]
    g1 /= (D[1]*α[1]*n[1]^2*β*(1 + exp(-2*α[1]*(l[1] + zb[1]))) + D[2]*α[2]*n[2]^2*γ*(1 + exp(-2*α[1]*(l[1] + zb[1]))))

    g += g1

    return g
end


##### for N = 2, expand sinh and cosh in terms of exp and factor out common exp and combine
function green2cylin(μsp, μa, sn, l, n)

    D = @. 1/(3(μsp))
    α = @. sqrt(sn^2 + μa/D)

    A = 1.0
    zb = @. 2*A*D
    z0 = 1/(μsp[1])
    
    g = (exp(-α[1]*z0) - exp(-α[1]*(z0 + 2*zb[1]))) / (2*D[1]*α[1])

    g1 = exp(α[1]*(z0 - 2*l[1]))*(1 - exp(-2*α[1]*(z0 + zb[1])))*(1 - exp(-2*α[1]*zb[1]))
    g1 *= (D[1]*α[1]*n[1]^2*(1 - exp(-2*α[2]*(l[2] + zb[2]))) - D[2]*α[2]*n[2]^2*(1 + exp(-2*α[2]*(l[2] + zb[2]))))
    g1 /= D[1]*α[1]
    g1 /= (D[1]*α[1]*n[1]^2*(1 - exp(-2*α[2]*(l[2] + zb[2])))*(1 + exp(-2*α[1]*(l[1] + zb[1]))) + D[2]*α[2]*n[2]^2*(1 + exp(-2*α[2]*(l[2] + zb[2])))*(1 + exp(-2*α[1]*(l[1] + zb[1]))))

    return g + g1
end

## equation 22
function fluence_DA_2lay_cylinder_CW(μsp, μa, l, a, ρ, n)

    A = 1
    D = 1/3μsp[1]
    zb = 2*A*D[1]

    ϕ = 0.0

    for ind in eachindex(besselroots)
        ϕ += green2cylin(μsp, μa, besselroots[ind]/(a + zb), l, n)*besselj0(besselroots[ind]/(a + zb)*ρ)/(besselj1(besselroots[ind]))^2
    end

    return ϕ/(π*(a + zb)^2)
end


##### for N = 3, same concept


function green3cylin(μsp, μa, sn, l, n)
  
    D = @. 1/(3(μsp))
    α = @. sqrt(sn^2 + μa/D)

    A = 1.0
    zb = @. 2*A*D
    z0 = 1/(μsp[1])  

    β = D[2]*α[2]*n[2]^2*(1 + exp(-2*α[2]*l[2]))*(1 - exp(-2*α[3]*(l[3] + zb[2])))
    β += D[3]*α[3]*n[3]^2*(1 - exp(-2*α[2]*l[2]))*(1 + exp(-2*α[3]*(l[3] + zb[2])))

    γ = D[2]*α[2]*n[2]^2*(1 - exp(-2*α[2]*l[2]))*(1 - exp(-2*α[3]*(l[3] + zb[2])))
    γ += D[3]*α[3]*n[3]^2*(1 + exp(-2*α[2]*l[2]))*(1 + exp(-2*α[3]*(l[3] + zb[2])))

    
    
    g = (exp(-α[1]*z0) - exp(-α[1]*(z0 + 2*zb[1]))) / (2*D[1]*α[1])

    g1 = exp(α[1]*(z0 - 2*l[1]))*(1 - exp(-2*α[1]*(z0 + zb[1])))*(1 - exp(-2*α[1]*zb[1]))
    g1 *= (D[1]*α[1]*n[1]^2*β - D[2]*α[2]*n[2]^2*γ)
    g1 /= D[1]*α[1]
    g1 /= (D[1]*α[1]*n[1]^2*β*(1 + exp(-2*α[1]*(l[1] + zb[1]))) + D[2]*α[2]*n[2]^2*γ*(1 + exp(-2*α[1]*(l[1] + zb[1]))))


    return g + g1
end

## equation 22
function fluence_DA_3lay_cylinder_CW(μsp, μa, l, a, ρ, n)

    A = 1
    D = 1/3μsp[1]
    zb = 2*A*D[1]

    ϕ = 0.0

    for ind in eachindex(besselroots)
        ϕ += green3cylin(μsp, μa, besselroots[ind]/(a + zb), l, n)*besselj0(besselroots[ind]/(a + zb)*ρ)/(besselj1(besselroots[ind]))^2
    end

    return ϕ/(π*(a + zb)^2)
end



### for n = 4


function green4cylin(μsp, μa, sn, l, n)
  
    D = @. 1/(3(μsp))
    α = @. sqrt(sn^2 + μa/D)

    A = 1.0
    zb = @. 2*A*D
    z0 = 1/(μsp[1])  

    β_4 = D[3]*α[3]*n[3]^2*(1 + exp(-2*α[3]*l[3]))*(1 - exp(-2*α[4]*(l[4] + zb[2])))
    β_4 += D[4]*α[4]*n[4]^2*(1 - exp(-2*α[3]*l[3]))*(1 + exp(-2*α[4]*(l[4] + zb[2])))

    γ_4 = D[3]*α[3]*n[3]^2*(1 - exp(-2*α[3]*l[3]))*(1 - exp(-2*α[4]*(l[4] + zb[2])))
    γ_4 += D[4]*α[4]*n[4]^2*(1 + exp(-2*α[3]*l[3]))*(1 + exp(-2*α[4]*(l[4] + zb[2])))

    β = (D[2]*α[2]*n[2]^2*β_4*(1 + exp(-2*α[2]*l[2])) +  D[3]*α[3]*n[3]^2*γ_4*(1 - exp(-2*α[2]*l[2])))

    γ = (D[2]*α[2]*n[2]^2*β_4*(1 - exp(-2*α[2]*l[2])) +  D[3]*α[3]*n[3]^2*γ_4*(1 + exp(-2*α[2]*l[2])))
    
    
    g = (exp(-α[1]*z0) - exp(-α[1]*(z0 + 2*zb[1]))) / (2*D[1]*α[1])

    g1 = exp(α[1]*(z0 - 2*l[1]))*(1 - exp(-2*α[1]*(z0 + zb[1])))*(1 - exp(-2*α[1]*zb[1]))
    g1 *= (D[1]*α[1]*n[1]^2*β - D[2]*α[2]*n[2]^2*γ)
    g1 /= D[1]*α[1]
    g1 /= (D[1]*α[1]*n[1]^2*β*(1 + exp(-2*α[1]*(l[1] + zb[1]))) + D[2]*α[2]*n[2]^2*γ*(1 + exp(-2*α[1]*(l[1] + zb[1]))))


    return g + g1
end

function fluence_DA_4lay_cylinder_CW(μsp, μa, l, a, ρ, n)

    A = 1
    D = 1/3μsp[1]
    zb = 2*A*D[1]

    ϕ = 0.0

    for ind in eachindex(besselroots[1:1000])
        ϕ += green4cylin(μsp, μa, besselroots[ind]/(a + zb), l, n)*besselj0(besselroots[ind]/(a + zb)*ρ)/(besselj1(besselroots[ind]))^2
    end

    return ϕ/(π*(a + zb)^2)
end