#####################################################################################################################################
# Implements solution to the diffusion equation in an infinite medium as given in [1].
# Solutions are given in the spatial, frequency, and time domains for a point source located in the middle of the cylinder top.
#
# [1] Patterson et. al., "Time resolved reflectance and transmittance for the noninvasive measurement of tissue optical properties," 
#     Appl. Opt. 28, 2331-2336 (1989)
#####################################################################################################################################

"""
    fluence_DA_inf_CW(ρ, μa, μsp)

Compute the steady-state fluence in an infinite medium. 

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Examples
julia> fluence_DA_inf_CW(1.0, 0.1, 10.0)
"""
function fluence_DA_inf_CW(ρ, μa, μsp)
    @assert ρ > zero(eltype(ρ)) "ρ must be greater than zero"
    D = D_coeff(μsp, μa)
    return exp(-sqrt(3 * μsp * μa) * ρ) / (4 * π * ρ * D)
end

"""
    fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64)

Compute the time-domain fluence in an infinite medium with Eqn. 3 of Patterson. et al. 1989. 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ω`: modulation frequency (1/ns)
- `n_med::Float64`: medium's index of refraction

# Examples
julia> fluence_DA_inf_TD(0.1:0.5:5.0, 1.0, 0.1, 10.0)
"""
function fluence_DA_inf_TD(t, ρ, μa, μsp, n_med = 1.0)
    @assert ρ > zero(eltype(ρ))
    D = D_coeff(μsp, μa)
    ν = ν_coeff(n_med)

    ϕ = Array{eltype(ρ)}(undef, length(t))

    for n in eachindex(t)
        @assert t[n] > zero(eltype(t)) "t must be greater than zero"
        tmp1 = 4 * D * ν * t[n]
        ϕ[n] = exp(-(ρ^2 / tmp1) - μa * ν * t[n])
        ϕ[n] *= ν / ((tmp1 * π )^(3 / 2))
    end

    return ϕ
end

"""
    fluence_DA_inf_FD(ρ, μa, μsp, ω, n_med)

Compute the fluence for a frequency modulated source in an infinite medium. 

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ω`: modulation frequency (1/ns)
- `n_med::Float64`: medium's index of refraction

# Examples
julia> fluence_DA_inf_FD(1.0, 0.1, 10.0, 1.0)
"""
function fluence_DA_inf_FD(ρ, μa, μsp, ω, n_med = 1.0)
    @assert ρ > zero(eltype(ρ)) "ρ must be greater than zero"
    D = D_coeff(μsp, μa)
    ν = ν_coeff(n_med)

    ϕ = exp(-ρ * sqrt(μa / D + im * ω / (ν * D))) / (4 * π * ρ * D)
    return ϕ
end