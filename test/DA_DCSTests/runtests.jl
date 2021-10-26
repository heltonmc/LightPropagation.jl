module DA_DCSTests

using Test
using LightPropagation: g2_DA_semiinf_CW, g2_DA_Nlay_cylinder_CW
using LightPropagation: DAsemiinf_DCS, Nlayer_cylinder_DCS, besselroots

τ = 10 .^(range(-10,stop=0,length=250))
ρ = 1.0; μa = 0.1; μsp = 10.0; n_med = 1.0; n_ext = 1.0
β = 1.0; λ = 700.0; z = 0.0
BFi = 2.0e-8

data = DAsemiinf_DCS(ρ = ρ, μa = μa, μsp = μsp, n_med = n_med, n_ext = n_ext, β = β, λ = λ, BFi = BFi, z = z)

  
@test g2_DA_semiinf_CW(τ, data) ≈ g2_DA_semiinf_CW(τ, ρ, μa, μsp; BFi = BFi, β = β, n_ext = n_ext, n_med = n_med, z = z, λ = λ)

si = g2_DA_semiinf_CW(τ, data)

### test layered against semi-inf

μa = [0.1, 0.1]; μsp = [10.0, 10.0]; n_med = [1.0, 1.0]; n_ext = 1.0
BFi = [2.0e-8, 2.0e-8]; l = [1.0, 10.0]; a = 25.0

data = Nlayer_cylinder_DCS(ρ = ρ, μa = μa, μsp = μsp, n_med = n_med, n_ext = n_ext, β = β, λ = λ, BFi = BFi, z = z, a = a, l = l, bessels = besselroots[1:2000])


@test g2_DA_Nlay_cylinder_CW(τ, data) ≈ g2_DA_Nlay_cylinder_CW(τ, ρ, μa, μsp; BFi = BFi, β = β, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, bessels = besselroots[1:2000])

layered = g2_DA_Nlay_cylinder_CW(τ, data)

@test si ≈ layered

end
