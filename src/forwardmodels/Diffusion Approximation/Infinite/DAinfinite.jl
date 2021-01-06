"""
    fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64)

Compute the time-domain fluence in an infinite medium with Eqn. 3 of Patterson. et al. 1989. 

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `nmed::Float64`: the sample medium's index of refraction

# Examples
julia> fluence_DA_inf_TD(0:1:5, [0.1,10.0], 1.0, 1.0)
"""
function fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64 = 1.0)
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
    ν::Float64 = 29.9792458/nmed

    ϕ = Array{Float64}(undef, length(t))

    Threads.@threads for n in eachindex(t)

        ϕ[n] = -(ρ^2/(4D*ν*t[n]))
        ϕ[n] = ϕ[n] - μa*ν*t[n]
        ϕ[n] = ν*exp(ϕ[n])/((4π*D*ν*t[n])^(3/2))

        if isnan(ϕ[n])
            ϕ[n] = 0
        end
 
    end

    return ϕ
end

"""
    fluence_DA_inf_CW(ρ::Float64, β::Array{Float64,1})

Compute the fluence for a steady-state source in an infinite medium. 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)

# Examples
julia> fluence_DA_inf_CW(1.0, [0.1,10.0])
"""
function fluence_DA_inf_CW(ρ::Float64, β::Array{Float64,1})
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp

    ϕ = exp(-sqrt(3*μsp*μa)*ρ)/(4*π*ρ*D)

    if isnan(ϕ)
        ϕ = 0
    end

    return ϕ
end

"""
    fluence_DA_inf_FD(ρ::Float64, β::Array{Float64,1})

Compute the fluence for a steady-state source in an infinite medium. 

# Arguments
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)

# Examples
julia> fluence_DA_inf_FD(1.0, [0.1,10.0])
"""
function fluence_DA_inf_FD(ρ, β, ω, nmed)
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
    ν = 29.9792458/nmed

    ϕ = exp(-ρ*sqrt(μa/D + im*ω/(ν*D)))/(4*π*ρ*D)

    if isnan(ϕ)
        ϕ = 0
    end

    return ϕ
end

N = 512
T = 10e-12
ω = 2*π*Vector(0:N-1)/(N*T)

ω = 0.0977*Vector(1:512)

map((ω) -> fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0), ω)
quadgk((ω) -> exp(im*ω*1)*fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0), -Inf, Inf)

function TD_test(t, N, ub)

    Rt = Array{Float64,1}(undef, length(t))
    
    x, w = gauss(N, 0, ub)

    map!(t -> sum(FD.(x, t) .* w), Rt, t)
   
   
return Rt./pi
end

FD(ω, t) = real.(exp(im*ω*t)*fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0))
x, w = gauss(1000, 0, 500)
sum(FD.(x, 1.0) .* w)

a = map(ω -> fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0)*exp(im*ω*t), [0,1,2, 3, 4, -3, -2, -1])


r = map(t -> (quadgk((ω) -> exp(im*ω*t)*fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0), -Inf, Inf))[1], 0.01:0.05:3)

map(t -> (quadgk((ω) -> 1/2pi*exp(im*ω*t)*fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0), 0, Inf))[1], 0.01:0.05:3)

map(t -> (quadgk((ω) -> exp(im*ω*t)*fluence_FD(1.0, 0., [0.1,0.1], [10.,10.], 1.0, ω, 1.0), 0, Inf))[1], 0.01:0.05:3)
r1 = map(t -> (quadgk((ω) -> exp(im*ω*t)*fluence_FD(1.0, 0., [0.1,0.1], [10.,10.], 1.0, ω, 1.0), -1000, 1000))[1], 1:1:3)

map(t -> (quadgk((ω) -> exp(im*ω*t)*fluence_DA_inf_FD(1.0, [0.1,10.0], ω, 1.0), 0, 1024))[1], 0.01:0.05:3)


function FD(ω)
    D = 1/30
    ν = 29.9792458 #speed of light in cm/ns

    ϕ = exp(-sqrt(0.1/D + im*ω/(ν*D)))
end

using QuadGK

tind = 0.1:0.1:3
r = map(t -> (quadgk((ω) -> exp(im*ω*t)*FD(ω), -Inf, Inf))[1], tind)
