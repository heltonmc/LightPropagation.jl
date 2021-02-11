# straightforward implementation of equation 37: this is not stable for large lz 
function green(μsp, μa, sn, lz)
    A = 1
    D = 1/3μsp
    ν = 30.0
    ω = 0.0

    α = sqrt(μa/D + sn^2 + im*ω/(D*ν))


    zb = 2*A*D
    z0 = 1/μsp    
    
    
    g = (exp(-α*z0) - exp(-α*(z0 + 2*zb))) / (2*D*α)
    g -= sinh(α*(z0 + zb))*sinh(α*(zb))*exp(-α*(lz + 2zb)) / (D*α*sinh(α*(lz + 2zb)))

    return g
end


function greencylin_FD(μsp, μa, sn, lz, ω)
    A = 1
    D = 1/3μsp
    ν = 29.9792345 #/nmed

    α = sqrt(μa/D + sn^2 + im*ω/(D*ν))


    zb = 2*A*D
    z0 = 1/μsp    
    
    
    g = (exp(-α*z0) - exp(-α*(z0 + 2*zb)))
    g -= exp(α*(z0 - 2*lz - 2*zb))*(1 - exp(-2*α*(z0 + zb)))*(1 - exp(-2*α*zb))/(1 - exp(-2*α*(lz + 2*zb)))

    return g / (2*D*α)
end

function greencylin_CW(μsp, μa, sn, lz)
    A = 1
    D = 1/3μsp

    α = sqrt(μa/D + sn^2)


    zb = 2*A*D
    z0 = 1/μsp    
    
    
    g = (exp(-α*z0) - exp(-α*(z0 + 2*zb)))
    g -= exp(α*(z0 - 2*lz - 2*zb))*(1 - exp(-2*α*(z0 + zb)))*(1 - exp(-2*α*zb))/(1 - exp(-2*α*(lz + 2*zb)))

    return g / (2*D*α)
end

using SpecialFunctions, Roots

besselj_zero_asymptotic(nu, n) = pi * (n - 1 + nu / 2 + 3//4)

#besselj_zero_asymptotic(n, A) = π*(4*n + 3) / (4*A) # with factor A and 0 order
besselj_zero(n, A; order=2) = Roots.fzero((x) -> SpecialFunctions.besselj0(A*x), besselj_zero_asymptotic(n, A); order=order)
sn = map(s -> besselj_zero(s, 1; order = 2), 1:100)

const besselroots = map(s -> besselj_zero(s, 1; order = 2), 0:10000)



# eqn 21 fluence of finite cylinder 
function fluence_DA_cylinder_CW(μsp, μa, lz, a, ρ)

    A = 1
    D = 1/3μsp
    zb = 2*A*D

    sn = map(s -> besselj_zero(s, (a + zb); order = 2), 0:10000)

    ϕ = map(s -> real.(greencylin_CW(μsp, μa, s, lz))*besselj0(s*ρ)/(besselj1((a + zb)*s))^2, sn)

    return sum(ϕ)/(π*(a + zb)^2)
end

function fluence_DA_cylinder_FD(μsp, μa, lz, a, ρ, ω)

    A = 1
    D = 1/3μsp
    zb = 2*A*D

    sn = map(s -> besselj_zero(s, (a + zb); order = 2), 0:10000)

    ϕ = map(s -> greencylin_FD(μsp, μa, s, lz, ω)*besselj0(s*ρ)/(besselj1((a + zb)*s))^2, sn)

    return sum(ϕ)/(π*(a + zb)^2)
end

#improved functions that use besselroots constant and loop over indices 
function fluence_DA_cylinder_FD1(μsp, μa, lz, a, ρ, ω)

    A = 1
    D = 1/3μsp
    zb = 2*A*D

    ϕ = 0.0 + 0.0im

    for n in eachindex(besselroots)
        ϕ += greencylin_FD(μsp, μa, besselroots[n]/(a + zb), lz, ω)*besselj0(besselroots[n]/(a + zb)*ρ)/(besselj1(besselroots[n]))^2
    end

    return ϕ/(π*(a + zb)^2)
end


function fluence_DA_cylinder_CW1(μsp, μa, lz, a, ρ)

    A = 1
    D = 1/3μsp
    zb = 2*A*D

    ϕ = 0.0

    for n in eachindex(besselroots)
        ϕ += greencylin_CW(μsp, μa, besselroots[n]/(a + zb), lz)*besselj0(besselroots[n]/(a + zb)*ρ)/(besselj1(besselroots[n]))^2
    end

    return ϕ/(π*(a + zb)^2)
end
#### TD

function fluence_DA_cylinder_TD(t, μa, μsp, lz, a, ρ)
    ub = 1200
	N = 1000

    Rt = zeros(Float64, length(t))
    freqcomp = zeros(Float64, N)
    x, w = gausslegendre(N)
	x = x.*ub/2 .+ ub/2
	w = w.*ub/2
    map!(f -> real.(fluence_DA_cylinder_FD(μsp, μa, lz, a, ρ, f)), freqcomp, x)
    
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m])))*freqcomp[m]*w[m]
        end
    end
    return Rt./94.182542865
end

function TD_test(t, μa, μsp, lz, a, ρ)
    ub = 1200
	N = 1000

    Rt = zeros(Float64, length(t))
    freqcomp = zeros(Float64, N)
    x, w = gausslegendre(N)
	x = x.*ub/2 .+ ub/2
	w = w.*ub/2
    map!(f -> real.(fluence_DA_cylinder_FD1(μsp, μa, lz, a, ρ, f)), freqcomp, x)
    
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m])))*freqcomp[m]*w[m]
        end
    end
    return Rt./94.182542865
end

t = 0.01:0.01:4.0
μsp = 10.0
μa = 0.1
ρ = 3.0
cylTD = TD_test(t, μa, μsp, 5.0, 10.0, ρ)
inds1 = cylTD .> 0
plot(t[inds1], cylTD[inds1], yscale =:log10)

siTD = fluence_DA_semiinf_TD(t, [μa, μsp], ρ, 1.0, 1.0, 0.0)
cylTD = TD_test(t, μa, μsp, 5.0, 10.0, ρ)

plot(t, siTD, yscale=:log10, label= "Semi-inf", lw=2)
plot!(t, cylTD, label= "Cylinder, r = 10cm", lw=2)


cylTD = TD_test(t, μa, μsp, 5.0, 5.0, 1.5)
plot!(t, cylTD, label= "Cylinder, r = 5cm", lw=2)


cylTD = TD_test(t, μa, μsp, 5.0, 3.0, 1.5)
plot!(t, cylTD, label= "Cylinder, r = 3cm", lw=2)


cylTD = TD_test(t, μa, μsp, 5.0, 2.0, 1.5)
plot!(t, cylTD, label= "Cylinder, r = 2cm", lw=2)



### paralepidp
fluence_DA_paralpip_TD(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})


t = 0.01:0.01:5.0
μsp = 10.0
μa = 0.2
ρ = 1.5

siTD = fluence_DA_semiinf_TD(t, [μa, μsp], ρ, 1.0, 1.0, 0.0)
parTD = fluence_DA_paralpip_TD(t, [μa, μsp], 1.0, 1.0, [8.5,10.0,0.0], [10.0,10.0], [20.0,20.0,20.0])


plot(t, siTD, yscale=:log10, label= "Semi-inf", lw=2)
plot!(t, parTD./59.96, label= "Cube, L = 20cm", lw=2)

parTD = fluence_DA_paralpip_TD(t, [μa, μsp], 1.0, 1.0, [3.5,5.0,0.0], [5.0,5.0], [10.0,10.0,10.0])
plot!(t, parTD./59.96, label= "Cube, L = 10cm", lw=2)


parTD = fluence_DA_paralpip_TD(t, [μa, μsp], 1.0, 1.0, [2.5,4.0,0.0], [4.0,4.0], [8.0,8.0,8.0])
plot!(t, parTD./59.96, label= "Cube, L = 8cm", lw=2)


parTD = fluence_DA_paralpip_TD(t, [μa, μsp], 1.0, 1.0, [1.0,2.5,0.0], [2.5,2.5], [5.0,5.0,5.0])
plot!(t, parTD./59.96, label= "Cube, L = 5cm", lw=2)


parTD = fluence_DA_paralpip_TD(t, [μa, μsp], 1.0, 1.0, [0.5,2.0,0.0], [2.0,2.0], [4.0,4.0,4.0])
plot!(t, parTD./59.96, label= "Cube, L = 4cm", lw=2)

TPSF_paralpip_8cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [3.0,4.0], [4.0,4.0], [8.0,8.0,8.0])
TPSF_paralpip_4cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [1.0,2.0], [2.0,2.0], [4.0,4.0,4.0])
