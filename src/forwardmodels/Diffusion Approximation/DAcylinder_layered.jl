@with_kw struct Nlayer_cylinder{T <: AbstractFloat}
    μsp::Array{T,1} = [10.0, 10.0, 10.0, 10.0]
    μa::Array{T,1} = [0.1, 0.1, 0.1, 0.1] 
    n_ext::T = 1.0                              # surrounding index of refraction
    n_med::Array{T,1} = [1.0, 1.0, 1.0, 1.0]    # layers index of refraction

    #source, detector
    l::Array{T,1} = [0.5, 0.8, 1.0, 5.0]        # length of cylinder layers
    ρ::T = 1.0                                  # source-detector separation
    a::T = 5.0                                  # radius of cylinder

    ω::T = 0.0                                  # modulation frequency
end

# Calculate β and γ coefficients 
function _get_βγ2(α, D, n, zb, l)
    β = (1 - exp(-2 * α[2] * (l[2] + zb[2])))

    γ = (1 + exp(-2 * α[2] * (l[2] + zb[2])))
    
    return β, γ
end
function _get_βγ3(α, D, n, zb, l)
    β  = D[2] * α[2] * n[2]^2 * (1 + exp(-2 * α[2] * l[2])) * (1 - exp(-2 * α[3] * (l[3] + zb[2])))
    β += D[3] * α[3] * n[3]^2 * (1 - exp(-2 * α[2] * l[2])) * (1 + exp(-2 * α[3] * (l[3] + zb[2])))

    γ  = D[2] * α[2] * n[2]^2 * (1 - exp(-2 * α[2] * l[2])) * (1 - exp(-2 * α[3] * (l[3] + zb[2])))
    γ += D[3] * α[3] * n[3]^2 * (1 + exp(-2 * α[2] * l[2])) * (1 + exp(-2 * α[3] * (l[3] + zb[2])))

    return β, γ
end
function _get_βγ4(α, D, n, zb, l)
    β_4  = D[3] * α[3] * n[3]^2 * (1 + exp(-2 * α[3] * l[3])) * (1 - exp(-2 * α[4] * (l[4] + zb[2])))
    β_4 += D[4] * α[4] * n[4]^2 * (1 - exp(-2 * α[3] * l[3])) * (1 + exp(-2 * α[4] * (l[4] + zb[2])))

    γ_4  = D[3] * α[3] * n[3]^2 * (1 - exp(-2 * α[3] * l[3])) * (1 - exp(-2 * α[4] * (l[4] + zb[2])))
    γ_4 += D[4] * α[4] * n[4]^2 * (1 + exp(-2 * α[3] * l[3])) * (1 + exp(-2 * α[4] * (l[4] + zb[2])))

    β = (D[2] * α[2] * n[2]^2 * β_4 * (1 + exp(-2 * α[2] * l[2])) + D[3] * α[3] * n[3]^2 * γ_4 * (1 - exp(-2 * α[2] * l[2])))
    γ = (D[2] * α[2] * n[2]^2 * β_4 * (1 - exp(-2 * α[2] * l[2])) + D[3] * α[3] * n[3]^2 * γ_4 * (1 + exp(-2 * α[2] * l[2])))

    return β, γ
end

function _green_Nlaycylin(α, D, z0, zb, l, n)
    if isequal(length(α), 4)
        β, γ = _get_βγ4(α, D, n, zb, l)
    elseif isequal(length(α), 3)
        β, γ = _get_βγ3(α, D, n, zb, l)
    elseif isequal(length(α), 2)
        β, γ = _get_βγ2(α, D, n, zb, l)
    end

    g   = (exp(-α[1] * z0) - exp(-α[1] * (z0 + 2 * zb[1]))) / (2 * D[1] * α[1])
    g1  = exp(α[1] * (z0 - 2 * l[1])) * (1 - exp(-2 * α[1] * (z0 + zb[1]))) * (1 - exp(-2 * α[1] * zb[1]))
    g1 *= (D[1] * α[1] * n[1]^2 * β - D[2] * α[2] * n[2]^2 * γ)
    g1 /= 2 * D[1] * α[1]
    g1 /= (D[1] * α[1] * n[1]^2 * β * (1 + exp(-2 * α[1] * (l[1] + zb[1]))) + D[2] * α[2] * n[2]^2 * γ * (1 - exp(-2 * α[1] * (l[1] + zb[1]))))

    return g + g1
end

function _green_Nlaycylin_CW(α, sn, μa, D, z0, zb, l, n)
    @inbounds for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2)
    end
    return _green_Nlaycylin(α, D, z0, zb, l, n)
end
function _green_Nlaycylin_FD(sn, μa, D, z0, zb, l, n, ν, ω)
    α = @. sqrt(μa / D + sn^2 + im * ω / (D * ν))
    return _green_Nlaycylin(α, D, z0, zb, l, n)
end
function _green_Nlaycylin_Laplace(α, sn, μa, D, z0, zb, l, n, ν, s)
    for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2 + s / (D[ind] * ν[ind]))
    end
    return _green_Nlaycylin(α, D, z0, zb, l, n)
end

function fluence_DA_Nlay_cylinder_CW(data::Nlayer_cylinder, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(eltype(data.ρ))
    ϕ_tmp = zero(eltype(ϕ))
    α = zeros(eltype(data.ρ), length(data.μa))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_CW(α, besselroots[ind]/(data.a + zb[1]), data.μa, D, z0, zb, data.l, data.n_med)
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ)
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ / (π * (data.a + zb[1])^2)
end
function fluence_DA_Nlay_cylinder_FD(data::Nlayer_cylinder, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(Complex{eltype(data.ρ)})
    ϕ_tmp = zero(eltype(ϕ))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_FD(besselroots[ind]/(data.a + zb[1]), data.μa, D, z0, zb, data.l, data.n_med, ν, ω) 
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ) 
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end
  
    return ϕ / (π * (data.a + zb[1])^2)
end
function _fluence_DA_Nlay_cylinder_Laplace(s, data::Nlayer_cylinder, besselroots)
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(eltype(s))
    ϕ_tmp = zero(eltype(ϕ))

    α = zeros(eltype(s), length(data.μa))

    for ind in eachindex(besselroots)
        ϕ_tmp = _green_Nlaycylin_Laplace(α, besselroots[ind]/(data.a + zb[1]), data.μa, D, z0, zb, data.l, data.n_med, ν, s) 
        ϕ_tmp *= besselj0(besselroots[ind] / (data.a + zb[1]) * data.ρ) 
        ϕ_tmp /= (besselj1(besselroots[ind]))^2
        ϕ += ϕ_tmp
    end

    return ϕ / (π * (data.a + zb[1])^2)
end

function fluence_DA_Nlay_cylinder_TD(t, N, data::Nlayer_cylinder, besselroots)
    Rt = zeros(eltype(t), length(t))
    Rt = LT_hyperbola(s -> _fluence_DA_Nlay_cylinder_Laplace(s, data, besselroots), N, t)

    return Rt
end

function fluence_DA_Nlay_cylinder_TD(t, data::Nlayer_cylinder, besselroots; N = 28)
    Rt = zeros(eltype(t), length(t))
    Rt = LT_hyper_fixed(s -> _fluence_DA_Nlay_cylinder_Laplace(s, data, besselroots), N, t)

    return Rt
end
