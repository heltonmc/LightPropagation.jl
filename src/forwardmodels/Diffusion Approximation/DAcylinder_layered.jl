#---------------------------------------------------------------------------------------------------------------------------------------- 
# Implements solution to the diffusion equation in a N-layered finite cylinder as given in [1].
# Solutions are given in the spatial, frequency, and time domains for a point source located in the middle of the cylinder top.
#
# [1] André Liemert and Alwin Kienle, "Light diffusion in a turbid cylinder. II. Layered case," Opt. Express 18, 9266-9279 (2010) 
#---------------------------------------------------------------------------------------------------------------------------------------- 

#-------------------------------------------
# Steady-State Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))

Compute the steady-state fluence in an N-layered cylinder.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficients in each layer (cm⁻¹)
- `μsp`: reduced scattering coefficient in each layer (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction in each layer
- `z`: the z-depth orthogonal from the boundary (cm) within the cylinder
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `MaxIter`: the maximum number of terms to consider in the infinite sum
- `atol`: the infinite sum will break after this absolute tolerance is met

The source must be located in the first layer (l[1] > 1/μsp[1]). Other arguments are not checked but should be restricted to:
- ρ, μa, and z >= 0.0
- μsp, n_ext, n_med, a > 0.0
- l .> 0.0
- ρ < a
- length(μa) == length(μsp) == length(n_med) == length(l)

μa, μsp, n_med, l should be tuples of the same length (e.g., (0.1, 0.1)) The input parameters should be of the same type, but will work with mixed types. 
The routine can be accurate until the machine precision usedin the calculation. 
Therefore, atol should be >= eps(eltype(μsp)). Larger values of μsp[1] will require a larger number of terms in the summation.
It is recommended to increase MaxIter if simulating higher scattering coefficients. It is also recommended to keep the cylinder radius as small as 
possible to increase convergence rate.  


# Examples
```
julia> fluence_DA_Nlay_cylinder_CW(1.0, (0.1, 0.1), (10.0, 10.0))
julia> fluence_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1), (12.0, 10.0), l = (10.0, 10.0), MaxIter=1000, atol=1.0e-8)
# to simulate a 3 layer media
julia> fluence_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=1000, atol=1.0e-8)
```
"""
function fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))
    D = D_coeff.(μsp)
    N = length(D)
    A = A_coeff.(n_med ./ n_ext)
    z0 = z0_coeff(μsp[1])
    zb = zb_coeff.(A, D)
    n_med = @. D * n_med^2
    @assert z0 < l[1] "The source must be located in the first layer (l[1] > 1/μsp[1])"
    @assert length(μa) == length(μsp) == length(n_med) == length(l) "μa, μsp, n_med, l should be tuples of the same length"
    
    if z < l[1]
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, MaxIter, _green_Nlaycylin_top, N, atol)
    elseif z > sum(l[1:end - 1])
        return _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, MaxIter, _green_Nlaycylin_bottom, N, atol)
    end
end
function fluence_DA_Nlay_cylinder_CW_approx(ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))
    @assert length(μa) == length(μsp) == length(n_med) == length(l) "μa, μsp, n_med, l should be tuples of the same length"
    D = D_coeff.(μsp)
    N = length(D)
    A = A_coeff.(n_med ./ n_ext)
    z0 = z0_coeff(μsp[1])
    zb = zb_coeff.(A, D)
    D_nmed2 = @. D * n_med^2
    @assert z0 < l[1] "The source must be located in the first layer (l[1] > 1/μsp[1])"
    abs(z - z0) > 0.5 && @warn "This approximation may yield inaccurate results. Consider using the exact version"
    
    if z <= l[1]
        ϕ = _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, D_nmed2, MaxIter, _green_approx_z_z0, N, atol)
    else
        throw(DomainError(z, "This approximation should only be used to calculate the fluence in the first layer when z ≈ z0"))
    end
    ϕ = @. ϕ + fluence_DA_semiinf_CW(ρ, μa[1], μsp[1], n_ext = n_ext, n_med = n_med[1], z = z)
    return ϕ
end

#-------------------------------------------
# Steady-State Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))

Compute the steady-state flux using Fick's law D[1]*∂ϕ(ρ)/∂z for z = 0 (reflectance) and -D[end]*∂ϕ(ρ)/∂z for z = sum(l) (transmittance) in an N-layered cylinder.
z must be equal to 0 or the total length sum(l) of cylinder.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficients in each layer (cm⁻¹)
- `μsp`: reduced scattering coefficient in each layer (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction in each layer
- `z`: the z-depth orthogonal from the boundary (cm) within the cylinder
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `MaxIter`: the maximum number of terms to consider in the infinite sum
- `atol`: the infinite sum will break after this absolute tolerance is met

The source must be located in the first layer (l[1] > 1/μsp[1]). Other arguments are not checked but should be restricted to:
- ρ, μa, and z >= 0.0
- μsp, n_ext, n_med, a > 0.0
- l .> 0.0
- ρ < a
- length(μa) == length(μsp) == length(n_med) == length(l)

μa, μsp, n_med, l should be tuples of the same length (e.g., (0.1, 0.1)) The input parameters should be of the same type, but will work with mixed types. 
The routine can be accurate until the machine precision usedin the calculation. 
Therefore, atol should be >= eps(eltype(μsp)). Larger values of μsp[1] will require a larger number of terms in the summation.
It is recommended to increase MaxIter if simulating higher scattering coefficients. It is also recommended to keep the cylinder radius as small as 
possible to increase convergence rate.  


# Examples
```
julia> flux_DA_Nlay_cylinder_CW(1.0, (0.1, 0.1), (10.0, 10.0))
julia> flux_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1), (12.0, 10.0), l = (10.0, 10.0), MaxIter=1000, atol=1.0e-8)
# to simulate a 3 layer media
julia> flux_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=1000, atol=1.0e-8)
```
"""
function flux_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=n_ext, n_med=n_med, l=l, a=a, z=dz, MaxIter=MaxIter, atol=atol), z)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp; n_ext=n_ext, n_med=n_med, l=l, a=a, z=dz, MaxIter=MaxIter, atol=atol), z)
    end
end

#-------------------------------------------
# Time-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N, ILT)

Compute the time-domain fluence in an N-layered cylinder.

# Arguments
- `t`: time point or vector of time values (ns)
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficients in each layer (cm⁻¹)
- `μsp`: reduced scattering coefficient in each layer (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction in each layer
- `z`: the z-depth orthogonal from the boundary (cm) within the cylinder
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `MaxIter`: the maximum number of terms to consider in the infinite sum
- `atol`: the infinite sum will break after this absolute tolerance is met

The source must be located in the first layer (l[1] > 1/μsp[1]). Other arguments are not checked but should be restricted to:
- ρ, μa, and z >= 0.0
- μsp, n_ext, n_med, a > 0.0
- l .> 0.0
- ρ < a
- length(μa) == length(μsp) == length(n_med) == length(l)
- t .> 0.0

t can be a scaler or a vector but all values of t should be greater than zero and if a vector should be ordered. If t is a scaler `hyperbola` is used for 
the inverse Laplace transform where `hyper_fixed` is used if it is a vector. The value of N should be proportional to the dynamic range of the time-domain signal needed.
If t[end]/t[1] is large then a higher value of N will be required. In general, a larger N will be needed to reconstruct very late times. 
μa, μsp, n_med, l should be tuples of the same length (e.g., (0.1, 0.1)) The input parameters should be of the same type, but will work with mixed types. 
The routine can be accurate until the machine precision used in the calculation. 
Therefore, atol should be >= eps(eltype(μsp)). Larger values of μsp[1] will require a larger number of terms in the summation.
It is recommended to increase MaxIter if simulating higher scattering coefficients. It is also recommended to keep the cylinder radius as small as 
possible to increase convergence rate.

The Laplace transform implements threaded parallelism. It is recommended to start Julia with multiple threads (check with Threads.nthreads() in the REPL)

# Examples
```
# at a single time point
julia> fluence_DA_Nlay_cylinder_TD(1.0, 1.0, (0.1, 0.1), (10.0, 10.0))
# for several time points
julia> fluence_DA_Nlay_cylinder_TD(0.2:0.2:2.0, 1.0, (0.1, 0.1), (10.0, 10.0))
julia> fluence_DA_Nlay_cylinder_TD([0.5, 0.8, 1.2], 1.0, (0.1, 0.1), (10.0, 10.0))
julia> fluence_DA_Nlay_cylinder_TD(0.1:0.3:5.0, 1.0, (0.1, 0.2), (12.0, 10.0), N = 48)
# to simulate a 3 layer media
julia> fluence_DA_Nlay_cylinder_TD(0.1:0.2:1.2, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=1000, atol=1.0e-8)
```
"""
function fluence_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N = 18, ILT = hyperbola)
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, MaxIter, atol), t, N = N)
end
function fluence_DA_Nlay_cylinder_TD(t::AbstractArray, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N = 24, ILT = hyper_fixed)
    T = promote_type(eltype(ρ), eltype(μa), eltype(μsp), eltype(z), eltype(l))
    return ILT(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, MaxIter, atol), t, N = N, T = T)
end

#-------------------------------------------
# Time-Domain Flux 
#-------------------------------------------
"""
    flux_DA_Nlay_cylinder_TD(t, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N, ILT)

Compute the time-domain flux using Fick's law D[1]*∂ϕ(ρ, t)/∂z for z = 0 (reflectance) and -D[end]*∂ϕ(ρ, t)/∂z for z = sum(l) (transmittance) in an N-layered cylinder.
∂ϕ(ρ, t)/∂z is calculated using forward mode auto-differentiation with ForwardDiff.jl

# Arguments
- `t`: time point or vector of time values (ns)
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficients in each layer (cm⁻¹)
- `μsp`: reduced scattering coefficient in each layer (cm⁻¹)

# Keyword arguments
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction in each layer
- `z`: the z-depth orthogonal from the boundary (cm) within the cylinder
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `MaxIter`: the maximum number of terms to consider in the infinite sum
- `atol`: the infinite sum will break after this absolute tolerance is met

The source must be located in the first layer (l[1] > 1/μsp[1]). Other arguments are not checked but should be restricted to:
- ρ, μa, and z >= 0.0
- μsp, n_ext, n_med, a > 0.0
- l .> 0.0
- ρ < a
- length(μa) == length(μsp) == length(n_med) == length(l)
- t .> 0.0

t can be a scaler or a vector but all values of t should be greater than zero and if a vector should be ordered. If t is a scaler `hyperbola` is used for 
the inverse Laplace transform where `hyper_fixed` is used if it is a vector. The value of N should be proportional to the dynamic range of the time-domain signal needed.
If t[end]/t[1] is large then a higher value of N will be required. In general, a larger N will be needed to reconstruct very late times. 
μa, μsp, n_med, l should be tuples of the same length (e.g., (0.1, 0.1)) The input parameters should be of the same type, but will work with mixed types. 
The routine can be accurate until the machine precision used in the calculation. 
Therefore, atol should be >= eps(eltype(μsp)). Larger values of μsp[1] will require a larger number of terms in the summation.
It is recommended to increase MaxIter if simulating higher scattering coefficients. It is also recommended to keep the cylinder radius as small as 
possible to increase convergence rate.

The Laplace transform implements threaded parallelism. It is recommended to start Julia with multiple threads (check with Threads.nthreads() in the REPL)

# Examples
```
# at a single time point
julia> flux_DA_Nlay_cylinder_TD(1.0, 1.0, (0.1, 0.1), (10.0, 10.0))
# for several time points
julia> flux_DA_Nlay_cylinder_TD(0.2:0.2:2.0, 1.0, (0.1, 0.1), (10.0, 10.0))
julia> flux_DA_Nlay_cylinder_TD([0.5, 0.8, 1.2], 1.0, (0.1, 0.1), (10.0, 10.0))
julia> flux_DA_Nlay_cylinder_TD(0.1:0.3:5.0, 1.0, (0.1, 0.2), (12.0, 10.0), N = 48)
# to simulate a 3 layer media
julia> flux_DA_Nlay_cylinder_TD(0.1:0.2:1.2, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=1000, atol=1.0e-8)
```
"""
function flux_DA_Nlay_cylinder_TD(t::AbstractVector, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)

    ## need to know the type of z to preallocate vectors coorectly in hyper_fixed for autodiff
    function _ILT(z::T) where {T} 
        return hyper_fixed(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, MaxIter, atol), t, N = N, T = T)
    end
    if z == zero(eltype(z))
        return D[1]*ForwardDiff.derivative(_ILT, z)
    elseif z == sum(l)
        return -D[end]*ForwardDiff.derivative(_ILT, z)
    end
end
function flux_DA_Nlay_cylinder_TD(t::AbstractFloat, ρ, μa, μsp; n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64), N = 24)
    @assert z == zero(eltype(z)) || z == sum(l)
    D = D_coeff.(μsp)  
    if z == zero(eltype(z))
        return D[1] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, MaxIter, atol), t, N = N), 0.0)
    elseif z == sum(l)
        return -D[end] * ForwardDiff.derivative(dz -> hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, dz, s, MaxIter, atol), t, N = N), sum(l))
    end
end

#-------------------------------------------
# Frequency-Domain Fluence 
#-------------------------------------------
"""
    fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp; ω=1.0, n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))

Compute the frequency modulated fluence in an N-layered cylinder.

# Arguments
- `ρ`: source-detector separation in cylindrical coordinates (distance from middle z-axis of cylinder) (cm⁻¹)
- `μa`: absorption coefficients in each layer (cm⁻¹)
- `μsp`: reduced scattering coefficient in each layer (cm⁻¹)

# Keyword arguments
- `ω`: source modulation frequency (1/ns)
- `n_ext`: the boundary's index of refraction (air or detector)
- `n_med`: the sample medium's index of refraction in each layer
- `z`: the z-depth orthogonal from the boundary (cm) within the cylinder
- `l`: layer thicknesses (cm)
- `a`: cylinder radius (cm)
- `MaxIter`: the maximum number of terms to consider in the infinite sum
- `atol`: the infinite sum will break after this absolute tolerance is met

The source must be located in the first layer (l[1] > 1/μsp[1]). Other arguments are not checked but should be restricted to:
- ρ, μa, and z >= 0.0
- μsp, n_ext, n_med, a > 0.0
- l .> 0.0
- ρ < a
- length(μa) == length(μsp) == length(n_med) == length(l)

μa, μsp, n_med, l should be tuples of the same length (e.g., (0.1, 0.1)) The input parameters should be of the same type, but will work with mixed types. 
The routine can be accurate until the machine precision usedin the calculation. 
Therefore, atol should be >= eps(eltype(μsp)). Larger values of μsp[1] will require a larger number of terms in the summation.
It is recommended to increase MaxIter if simulating higher scattering coefficients. It is also recommended to keep the cylinder radius as small as 
possible to increase convergence rate.  

# Examples
```
julia> fluence_DA_Nlay_cylinder_FD(1.0, (0.1, 0.1), (10.0, 10.0))
julia> fluence_DA_Nlay_cylinder_FD(1.0, (0.2, 0.1), (12.0, 10.0), l = (10.0, 10.0), MaxIter=1000, atol=1.0e-8)
# to simulate a 3 layer media
julia> fluence_DA_Nlay_cylinder_FD(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=1000, atol=1.0e-8)
```
"""
function fluence_DA_Nlay_cylinder_FD(ρ, μa, μsp; ω=1.0, n_ext=1.0, n_med=(1.0, 1.0), l=(1.0, 5.0), a=10.0, z=0.0, MaxIter=10000, atol=eps(Float64))
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ ω * im ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex, μsp; n_ext=n_ext, n_med=n_med, l=l, a=a, z=z, MaxIter=MaxIter, atol=atol)
end

#-------------------------------------------------------------------------------
# this function is the base for the laplace transform in the time-domain
#-------------------------------------------------------------------------------
function _fluence_DA_Nlay_cylinder_Laplace(ρ, μa, μsp, n_ext, n_med, l, a, z, s, MaxIter, atol)
    ν = ν_coeff.(n_med)
    μa_complex = μa .+ s ./ ν

    return fluence_DA_Nlay_cylinder_CW(ρ, μa_complex,  μsp; n_ext=n_ext, n_med=n_med, l=l, a=a, z=z, MaxIter=MaxIter, atol=atol)
end

#-------------------------------------------------------------------------------
# _kernel_fluence_DA_Nlay_cylinder provides the main kernel of the program.
# The inputs are similar to fluence_DA_Nlay_cylinder_CW.
# D is the diffusion coefficient and N is the number of layers.
# green is the Green's function for either the first or bottom layer (below).
#-------------------------------------------------------------------------------
function _kernel_fluence_DA_Nlay_cylinder(ρ, D, μa, a, zb, z, z0, l, n_med, MaxIter, green, N, atol)
    ϕ = ρ .* zero(eltype(μa))
    ϕ_tmp = zero(eltype(μa))
    apzb = inv(a + zb[1])
    #T = eltype(ρ)
    
    local ind
    @inbounds for outer ind in 1:MaxIter#eachindex(besselroots)
        tmp = J0_ROOTS[ind] * apzb
        ϕ_tmp = green(tmp, μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp /= J1_J0ROOTS_2[ind] # replaces (besselj1(besselroots[ind]))^2
        ϕ = @. ϕ + ϕ_tmp * besselj0(tmp * ρ)
        abs(ϕ_tmp) < atol && break
    end
    #ind == length(1000) && @warn "Failed to converge to desired tolerance: Increase MaxIter"
    return ϕ ./ (π * (a + zb[1])^2)
end
function _kernel_fluence_DA_Nlay_cylinder(ρ::BigFloat, D, μa, a, zb, z, z0, l, n_med, MaxIter, green, N, atol)
    ϕ = ρ .* zero(eltype(μa))
    ϕ_tmp = zero(eltype(μa))
    apzb = inv(a + zb[1])
    #T = eltype(ρ)
    
    local ind
    @inbounds for outer ind in 1:MaxIter#eachindex(besselroots)
        tmp = J0_ROOTSbig[ind] * apzb
        ϕ_tmp = green(tmp, μa, D, z, z0, zb, l, n_med, N)
        ϕ_tmp /= J1_J0ROOTS_2big[ind] # replaces (besselj1(besselroots[ind]))^2
        ϕ = @. ϕ + ϕ_tmp * besselj0(tmp * ρ)
        abs(ϕ_tmp) < atol && break
    end
    #ind == length(1000) && @warn "Failed to converge to desired tolerance: Increase MaxIter"
    return ϕ ./ (π * (a + zb[1])^2)
end

#-------------------------------------------------------------------------------
# Calculates the Green's function in the first (top) and last (bottom) layer
# sinh and cosh have been expanded as exponentials.
# A common term has been factored out. This cancels for G1 but not GN (see _βγN_correction)
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
#-------------------------------------------------------------------------------
@inline function _green_Nlaycylin_top(sn, μa, D, z, z0, zb, l, n, N)
    α = @. sqrt(μa / D + sn^2)

    if N == 4
        β, γ = _get_βγ4(α, n, zb, l)
    elseif N == 3
        β, γ = _get_βγ3(α, n, zb, l)
    elseif N == 2
        β, γ = _get_βγ2(α, zb, l)
    elseif N > 4
        β, γ = _get_βγk(α, n, zb, l)
    end

    tmp1 = α[1] * n[1] * β
    tmp2 = α[2] * n[2] * γ
    tmp3 = exp(-2 * α[1] * (l[1] + zb[1]))

    g  = exp(-α[1] * abs(z - z0))
    g -= exp(-α[1] * (z + z0 + 2 * zb[1]))
    g1 = exp(α[1] * (z + z0 - 2 * l[1]))
    g1 *= (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * (z + zb[1])))
    g1 *= tmp1 - tmp2
    g1 /= muladd(tmp1, tmp3, tmp1) + muladd(tmp2, -tmp3, tmp2)

    return (g + g1) / (2 * D[1] * α[1])
end
@inline function _green_Nlaycylin_bottom(sn, μa, D, z, z0, zb, l, n, N)
    α = @. sqrt(μa / D + sn^2)

    if N == 4
        β, γ = _get_βγ4(α, n, zb, l)
        βγ_correction = _βγ4_correction(α, zb, l)  
    elseif N == 3
        β, γ = _get_βγ3(α, n, zb, l)
        βγ_correction = _βγ3_correction(α, zb, l)  
    elseif N == 2
        β, γ = _get_βγ2(α, zb, l)
        βγ_correction = _βγ2_correction(α, zb, l)  
    elseif N > 4
        β, γ = _get_βγk(α, n, zb, l)
        βγ_correction = _βγN_correction(α, zb, l)  
    end

    tmp = one(eltype(α))
    for ind in (N - 1):-1:2
        tmp *= α[ind] * n[ind]
    end

    tmp1 = exp(-2 * α[1] * (l[1] + zb[1]))
    gN = n[end] * tmp * 2^(N - 1) / 2 / D[end]
    gN *= exp(α[1] * (z0 - l[1]) + α[end] * (sum(l) + zb[end] - z) - βγ_correction)
    gN /= α[1] * n[1] * β * (1 + tmp1) + α[2] * n[2] * γ * (1 - tmp1)
    gN *= (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[end] * (sum(l) + zb[end] - z)))

    return gN
 end
 @inline function _green_approx_z_z0(sn, μa, D, z, z0, zb, l, n, N)
    α = @. sqrt(μa / D + sn^2)

    if N == 4
        β, γ = _get_βγ4(α, n, zb, l)
    elseif N == 3
        β, γ = _get_βγ3(α, n, zb, l)
    elseif N == 2
        β, γ = _get_βγ2(α, zb, l)
    elseif N > 4
        β, γ = _get_βγk(α, n, zb, l)
    end

    tmp1 = α[1] * n[1] * β
    tmp2 = α[2] * n[2] * γ
    tmp3 = exp(-2 * α[1] * (l[1] + zb[1]))

    g1 = exp(α[1] * (z + z0 - 2 * l[1]))
    g1 *= (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * (z + zb[1])))
    g1 *= tmp1 - tmp2
    g1 /= muladd(tmp1, tmp3, tmp1) + muladd(tmp2, -tmp3, tmp2)

    return g1 / (2 * D[1] * α[1])
end

#-------------------------------------------------------------------------------
# Calculate β and γ coefficients with eqn. 17 in [1].
# sinh and cosh have been expanded as exponentials.
# A common term has been factored out. This cancels for G1 but not GN (see _βγN_correction)
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
#-------------------------------------------------------------------------------
@inline function _get_βγ2(α, zb, l)
    tmp1 = exp(-2 * α[2] * (l[2] + zb[2]))
    
    β = (1 - tmp1)
    γ = (1 + tmp1)
    
    return β, γ
end
@inline function _get_βγ3(α, n, zb, l)
    tmp1 = α[2] * n[2]
    tmp2 = α[3] * n[3]
    tmp3 = exp(-2 * α[2] * l[2])
    tmp4 = exp(-2 * α[3] * (l[3] + zb[2]))

    a = tmp1 * tmp3
    b = muladd(-tmp1, tmp4, tmp1)
    c = a * tmp4
    a = a - c
    β = b + a
    γ = b - a

    a = tmp2 * tmp3
    b = muladd(tmp2, tmp4, tmp2)
    c = a * tmp4
    a = a + c

    β += b - a
    γ += b + a
    
    return β, γ
end
@inline function _get_βγ4(α, n, zb, l)
    tmp5 = α[2] * n[2]
    tmp1 = α[3] * n[3]
    tmp4 = α[4] * n[4]

    tmp6 = exp(-2 * α[2] * l[2])
    tmp2 = exp(-2 * α[3] * l[3])
    tmp3 = exp(-2 * α[4] * (l[4] + zb[2]))

    a = tmp1 * tmp3
    b = tmp1 * tmp2
    c = a * tmp2
    a = tmp1 - a
    b = b - c
    β4 = a + b
    γ4 = a - b

    a = tmp4 * tmp2
    b = tmp4 * tmp3
    c = a * tmp3
    b = tmp4 + b
    c = a + c
    β4 += b - c
    γ4 += b + c
    
    a = tmp5 * β4
    b = a * tmp6
    c = tmp1 * γ4
    d = c * tmp6
    a = a + c
    b = b - d

    β = a + b
    γ = a - b

    return β, γ
end
@inline function _get_βγk(α, n, zb, l)
    βN, γN = _get_βγN(α, n, zb, l)
  
    for k in length(α):-1:4
        tmp1 = α[k - 2] * n[k - 2] * βN
        tmp2 = α[k - 1] * n[k - 1] * γN
        tmp3 = exp(-2 * α[k - 2] * l[k - 2])

        a = tmp1 * tmp3
        c = tmp1 + tmp2
        b = tmp2 * tmp3
        d = a - b

        βN = c + d
        γN = c - d
    end

    return βN, γN
end
@inline function _get_βγN(α, n, zb, l)
    tmp1 = α[end - 1] * n[end - 1]
    tmp2 = α[end] * n[end]
    tmp3 = exp(-2 * α[end - 1] * l[end - 1])
    tmp4 = exp(-2 * α[end] * (l[end] + zb[2]))
    
    a = tmp1 * tmp3
    b = tmp1 - tmp1 * tmp4
    c = a * tmp4
    a = a - c
    βN =  b + a
    γN =  b - a

    a = tmp2 * tmp3
    b = tmp2 + tmp2 * tmp4
    c = a * tmp4
    a = a + c
    βN += b - a
    γN += b + a

    return βN, γN
end

#-------------------------------------------------------------------------------
# Calculate βγ correction for N layer green's function.
# This is done because a common term was factored out from general expressions
# They don't cancel for the N layer so have to divide them out.
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively
#-------------------------------------------------------------------------------
@inline function _βγ2_correction(α, zb, l)  
    return α[2] * (l[2] + zb[2])
end
function _βγ3_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * (l[3] + zb[2])
end
@inline function _βγ4_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * l[3] + α[4] * (l[4] + zb[2])
end
@inline function _βγ5_correction(α, zb, l)  
    return α[2] * l[2] + α[3] * l[3] + α[4] * l[4] + α[5] * (l[5] + zb[2])
end
@inline function _βγN_correction(α, zb, l) 
    out = α[end] * (l[end] + zb[2])
    for ind in (length(α) - 1):-1:2
        out += α[ind] * l[ind]
    end
    return out
end
#-------------------------------------------------------------------------------
