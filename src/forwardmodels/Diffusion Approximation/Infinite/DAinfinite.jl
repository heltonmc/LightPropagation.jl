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
function fluence_DA_inf_FD(ρ, β, ω, nmed = 1.0)
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
    ν = 29.9792458/nmed

    ϕ = exp(-ρ*sqrt(μa/D + im*ω/(ν*D)))/(4*π*ρ*D)

    return ϕ
end

using FastGaussQuadrature

function computeTD_fromFD(t, β, ρ, ub, N)	

    Rt = zeros(Float64, length(t))
    
    x, w = gausslegendre(N)
	
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m]*ub/2 +ub/2)))*real(fluence_DA_inf_FD(ρ, [β[1],β[2]], x[m]*ub/2 +ub/2))*w[m]*ub/pi
        end
    end
   
   
return Rt
end


function computeTD_fromFD(t, β, ρ, ub, N)	

    Rt = zeros(Float64, length(t))
    
    x, w = gausslegendre(N)
	
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m]*ub/2 +ub/2)))*real(fluence_DA_inf_FD(ρ, [β[1],β[2]], x[m]*ub/2 +ub/2))*w[m]*ub/pi
        end
    end
   
   
return Rt
end


function computeTD_fromFD2(t, β, ρ, ub, N)	

    Rt = zeros(Float64, length(t))
    freqcomp = zeros(Float64, N)
    x, w = gausslegendre(N)


    Threads.@threads for m in eachindex(x)
        freqcomp[m] = real(fluence_DA_inf_FD(ρ, [β[1],β[2]], x[m]*ub/2 +ub/2))*w[m]*ub/pi
    end

	
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m]*ub/2 +ub/2)))*freqcomp[m]
        end
    end
   
   
return Rt
end