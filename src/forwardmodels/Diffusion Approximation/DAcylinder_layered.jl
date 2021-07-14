################################################################################################################################## 
# Implements solution to the diffusion equation in a N-layered finite cylinder as given in [1].
# Solutions are given in the spatial, frequency, and time domains for a point source located in the middle of the cylinder top.
#
# [1] André Liemert and Alwin Kienle, "Light diffusion in a turbid cylinder. II. Layered case," Opt. Express 18, 9266-9279 (2010) 
################################################################################################################################## 

@with_kw struct Nlayer_cylinder{T <: AbstractFloat}
    μsp::Array{T,1} = [10.0, 10.0, 10.0, 10.0]  # reduced scattering coefficient (1/cm)
    μa::Array{T,1} = [0.1, 0.1, 0.1, 0.1]       # absorption coefficient (1/cm)
    n_ext::T = 1.0                              # surrounding index of refraction
    n_med::Array{T,1} = [1.0, 1.0, 1.0, 1.0]    # layers index of refraction

    l::Array{T,1} = [0.5, 0.8, 1.0, 5.0]        # length of cylinder layers (cm)
    ρ::T = 1.0                                  # source-detector separation (cm)
    a::T = 5.0                                  # radius of cylinder (cm)
    z::T = 0.0                                  # detector depth (cm)

    ω::T = 0.0                                  # modulation frequency
end

# Calculate β and γ coefficients with eqn. 17 in [1].
# sinh and cosh have been expanded as exponentials and factored.
# For N = 2, 3, 4 coefficients are explicitly calculated.
# For N > 4, β and γ are calculated recursively using eqn. 17 & 18.
function _get_βγ2(α, D, n, zb, l)
    tmp1 = exp(-2 * α[2] * (l[2] + zb[2]))
    
    β = (1 - tmp1)
    γ = (1 + tmp1)
    
    return β, γ
end
function _get_βγ3(α, D, n, zb, l)
    tmp1 = D[2] * α[2] * n[2]^2
    tmp2 = D[3] * α[3] * n[3]^2
    tmp3 = exp(-2 * α[2] * l[2])
    tmp4 = exp(-2 * α[3] * (l[3] + zb[2]))

    β  = tmp1 * (1 + tmp3) * (1 - tmp4)
    β += tmp2 * (1 - tmp3) * (1 + tmp4)

    γ  = tmp1 * (1 - tmp3) * (1 - tmp4)
    γ += tmp2 * (1 + tmp3) * (1 + tmp4)

    return β, γ
end
function _get_βγ4(α, D, n, zb, l)
    tmp5 = D[2] * α[2] * n[2]^2
    tmp1 = D[3] * α[3] * n[3]^2
    tmp4 = D[4] * α[4] * n[4]^2

    tmp6 = exp(-2 * α[2] * l[2])
    tmp2 = exp(-2 * α[3] * l[3])
    tmp3 = exp(-2 * α[4] * (l[4] + zb[2]))

    β_4  = tmp1 * (1 + tmp2) * (1 - tmp3)
    β_4 += tmp4 * (1 - tmp2) * (1 + tmp3)

    γ_4  = tmp1 * (1 - tmp2) * (1 - tmp3)
    γ_4 += tmp4 * (1 + tmp2) * (1 + tmp3)

    β = (tmp5 * β_4 * (1 + tmp6) + tmp1 * γ_4 * (1 - tmp6))
    γ = (tmp5 * β_4 * (1 - tmp6) + tmp1 * γ_4 * (1 + tmp6))

    return β, γ
end

# Calculates the Green's function G1 in the first layer using eqn. 15 from [1].
function _green_Nlaycylin(α, D, z, z0, zb, l, n)
    if isequal(length(α), 4)
        β, γ = _get_βγ4(α, D, n, zb, l)
    elseif isequal(length(α), 3)
        β, γ = _get_βγ3(α, D, n, zb, l)
    elseif isequal(length(α), 2)
        β, γ = _get_βγ2(α, D, n, zb, l)
    end

    tmp1 = D[1] * α[1] * n[1]^2 * β
    tmp2 = D[2] * α[2] * n[2]^2 * γ
    tmp3 = exp(-2 * α[1] * (l[1] + zb[1]))

    g  = (exp(-α[1] * abs(z - z0)) - exp(-α[1] * (z + z0 + 2 * zb[1])))
    g1 = exp(α[1] * (z + z0 - 2 * l[1])) * (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * (z + zb[1])))
    g1 *= (tmp1 - tmp2)
    g1 /= (tmp1 * (1 + tmp3) + tmp2 * (1 - tmp3))

    return (g + g1) / (2 * D[1] * α[1])
end

# Calculates α as given in eqn. 27 for each domain.
# In spatial (CW) domain, ω = 0. In Hankel-Laplace space, iω -> s.
function _green_Nlaycylin_CW(α, sn, μa, D, z, z0, zb, l, n)
    @inbounds for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2)
    end
    return _green_Nlaycylin(α, D, z, z0, zb, l, n)
end
function _green_Nlaycylin_FD(sn, μa, D, z, z0, zb, l, n, ν, ω)
    α = @. sqrt(μa / D + sn^2 + im * ω / (D * ν))
    return _green_Nlaycylin(α, D, z, z0, zb, l, n)
end
function _green_Nlaycylin_Laplace(α, sn, μa, D, z, z0, zb, l, n, ν, s)
    for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2 + s / (D[ind] * ν[ind]))
    end
    return _green_Nlaycylin(α, D, z, z0, zb, l, n)
end

# Calculates the fluence in the first layer (G1) due to a point source that is incident onto the center top of the cylinder with eqn. [22].
function fluence_DA_Nlay_cylinder_CW(data, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(eltype(data.ρ))
    ϕ_tmp = zero(eltype(ϕ))
    α = zeros(eltype(data.ρ), length(data.μa))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_CW(α, besselroots[ind] / (data.a + zb[1]), data.μa, D, data.z, z0, zb, data.l, data.n_med)
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ)
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ / (π * (data.a + zb[1])^2)
end
function fluence_DA_Nlay_cylinder_FD(data, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(Complex{eltype(data.ρ)})
    ϕ_tmp = zero(eltype(ϕ))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_FD(besselroots[ind] / (data.a + zb[1]), data.μa, D, data.z, z0, zb, data.l, data.n_med, ν, ω) 
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ) 
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end
  
    return ϕ / (π * (data.a + zb[1])^2)
end
function _fluence_DA_Nlay_cylinder_Laplace(s, data, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(eltype(s))
    ϕ_tmp = zero(eltype(ϕ))

    α = zeros(eltype(s), length(data.μa))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_Laplace(α, besselroots[ind] / (data.a + zb[1]), data.μa, D, data.z, z0, zb, data.l, data.n_med, ν, s) 
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ) 
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ / (π * (data.a + zb[1])^2)
end

# Calculates the fluence in the time-domain with the Laplace transform.
# This solution is not explicitly shown in [1].
function fluence_DA_Nlay_cylinder_TD(t::AbstractFloat, data; bessels = besselroots, N = 16)
    Rt = zeros(eltype(t), length(t))
    Rt = hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(s, data, bessels), t, N = N)

    return Rt
end
function fluence_DA_Nlay_cylinder_TD(t::AbstractArray, data; bessels = besselroots, N = 28)
    Rt = zeros(eltype(t), length(t))
    Rt = hyper_fixed(s -> _fluence_DA_Nlay_cylinder_Laplace(s, data, bessels), t, N = N)

    return Rt
end
