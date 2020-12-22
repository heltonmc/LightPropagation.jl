"""
fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64)

Compute the time-domain fluence in an infinite medium. 

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
    ν::Float64 = 29.9792345/nmed

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
