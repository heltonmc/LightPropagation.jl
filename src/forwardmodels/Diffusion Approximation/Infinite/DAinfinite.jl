#------------------------------------------------------------------------------------------------------------------------
# Implements solution to the diffusion equation in an infinite medium as given in [1].
# Solutions are given in the steady-state, frequency, and time domains for an isotroptic source.
#
# [1] Patterson et. al., "Time resolved reflectance and transmittance for the noninvasive measurement of tissue optical properties," 
#     Appl. Opt. 28, 2331-2336 (1989)
#-------------------------------------------------------------------------------------------------------------------------
""" Structure containing inputs for simulating the fluence under the diffusion approximation in the infinite space."""
struct DAinfinite{T <: AbstractFloat} <: DiffusionParameters
    ρ::T                                 # distance away from isotropic point source (cm)
    μa::T                                # absorption coefficient (cm⁻¹)
    μsp::T                               # reduced scattering coefficient (cm⁻¹)
    n_med::T                             # medium's index of refraction

    t::Union{T, AbstractVector{T}}       # time (ns)
    ω::T                                 # modulation frequency (1/ns)

    # do not need to provide these arguments (calculated from previous inputs)
    
    D::T                                 # Diffusion coefficient                        
    ν::T                                 # Speed of light in medium (cm/ns)
    function DAinfinite{T}(ρ::T, μa::T, μsp::T, n_med::T, t::Union{T, AbstractVector{T}}, ω::T) where {T <: AbstractFloat}
        @assert ρ > zero(T) "ρ must be positive"
        @assert μa >= zero(T) "μa must greater than or equal to 0"
        @assert μsp > zero(T) "μsp must greater than 0"
        @assert all(t .> zero(T)) "t must be greater than 0"
        return new{T}(ρ, μa, μsp, n_med, t, ω, D_coeff(μsp), ν_coeff(n_med))
    end
end

"""
    Generator function for DAinfinite structure.

Provides default parameters to use in the infinite space fluence calculation in either the CW, FD, or TD.

# Examples
```
julia> data = DAinfinite() # return default parameters
julia> data = DAinfinite(ρ = 1.5) # return ρ = 1.5 with the rest of the parameters given by defaults
julia> fluence_DA_inf_CW(data) # can then be used with corresponding functions
```
"""
function DAinfinite(;
    ρ::T = 1.0,
    μa::T = 0.1,
    μsp::T = 10.0,
    n_med::T = 1.0,
    ω::T = 0.0,
    t::Union{T, AbstractVector{T}} = 1.0
) where {T<:AbstractFloat}
    return DAinfinite{T}(ρ, μa, μsp, n_med, t, ω)
end

#-------------------------------------
# Steady-State Fluence 
#-------------------------------------

"""
    fluence_DA_inf_CW(ρ, μa, μsp)

Compute the steady-state fluence in an infinite medium. 

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Examples
```
julia> fluence_DA_inf_CW(1.0, 0.1, 10.0)
```
"""
function fluence_DA_inf_CW(ρ, μa, μsp; D = D_coeff(μsp))
    return exp(-sqrt(3 * μsp * μa) * ρ) / (4 * π * ρ * D)
end
"""
    fluence_DA_inf_CW(data::DiffusionParameters)

Wrapper to `fluence_DA_inf_CW(ρ, μa, μsp)`

# Examples
```
julia> data = DAinfinite(ρ = 1.0) # use structure to generate inputs
julia> fluence_DA_inf_CW(data) # then call the function
```
"""
function fluence_DA_inf_CW(data::DiffusionParameters)
    return fluence_DA_inf_CW(data.ρ, data.μa, data.μsp, D = data.D)
end

#-------------------------------------
# Frequency-Domain Fluence 
#-------------------------------------
"""
    fluence_DA_inf_FD(ρ, μa, μsp, ω; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))

Compute the fluence for a frequency modulated source in an infinite medium. 

# Arguments
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)
- `ω`: modulation frequency (1/ns)

# Optional Arguments
- `n_med`=1.0: medium's index of refraction
- `ν`: speed of light in the medium cm/ns
- `D`: Diffusion coefficient

# Examples
```
julia> fluence_DA_inf_FD(1.0, 0.1, 10.0, 1.0, n_med = 1.4)
```
"""
function fluence_DA_inf_FD(ρ, μa, μsp, ω; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))
    μa_complex = μa + ω * im / ν
    return fluence_DA_inf_CW(ρ, μa_complex, μsp, D = D)
end
"""
    fluence_DA_inf_FD(data::DiffusionParameters)

Wrapper to `fluence_DA_inf_FD(ρ, μa, μsp, ω; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))`

# Examples
```
julia> data = DAinfinite(ω = 1.0) # use structure to generate inputs
julia> fluence_DA_inf_FD(data) # then call the function
```
"""
function fluence_DA_inf_FD(data::DiffusionParameters)
    return fluence_DA_inf_FD(data.ρ, data.μa, data.μsp, data.ω, ν = data.ν, D = data.D)
end

#-------------------------------------
# Time-Domain Fluence 
#-------------------------------------
"""
    fluence_DA_inf_TD(t, ρ, μa, μsp; nmed = 1.0)

Compute the time-domain fluence in an infinite medium with Eqn. 3 of Patterson. et al. 1989. 

# Arguments
- `t`: the time vector (ns). 
- `ρ`: source-detector separation (cm⁻¹)
- `μa`: absorption coefficient (cm⁻¹)
- `μsp`: reduced scattering coefficient (cm⁻¹)

# Optional Arguments
- `n_med`=1.0: medium's index of refraction
- `ν`: speed of light in the medium cm/ns
- `D`: Diffusion coefficient

# Examples
```
julia> fluence_DA_inf_TD(0.1:0.5:5.0, 1.0, 0.1, 10.0, n_med = 1.4)
```
"""
function fluence_DA_inf_TD(t::AbstractVector, ρ, μa, μsp; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))
    T = promote_type(eltype(ρ), eltype(μa), eltype(μsp))
    ϕ = zeros(T, length(t))

    @inbounds Threads.@threads for ind in eachindex(t)
        ϕ[ind] = _kernel_fluence_DA_inf_TD(t[ind], ρ, μa, ν, D)
    end

    return ϕ
end
function fluence_DA_inf_TD(t::AbstractFloat, ρ, μa, μsp; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))
    return _kernel_fluence_DA_inf_TD(t, ρ, μa, ν, D)
end
"""
    fluence_DA_inf_TD(data::DiffusionParameters)

Wrapper to `fluence_DA_inf_TD(ρ, μa, μsp, ω; n_med = 1.0, ν = ν_coeff(n_med), D = D_coeff(μsp))`

# Examples
```
julia> data = DAinfinite(t = 0.01:0.1:2.0) # use structure to generate inputs
julia> fluence_DA_inf_TD(data) # then call the function
```
"""
function fluence_DA_inf_TD(data::DiffusionParameters)
    return fluence_DA_inf_TD(data.t, data.ρ, data.μa, data.μsp, n_med = data.n_med, ν = data.ν, D = data.D)
end

@inline function _kernel_fluence_DA_inf_TD(t, ρ, μa, ν, D)
    tmp1 = 4 * D * ν * t
    ϕ = exp(-(ρ^2 / tmp1) - μa * ν * t)
    ϕ *= ν / ((tmp1 * π )^(3/2))
end
