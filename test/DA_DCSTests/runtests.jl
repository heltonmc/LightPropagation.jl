module DA_DCSTests

using Test
using LightPropagation

τ = 10 .^(range(-10,stop=0,length=250))


# test that the sum of all pathlengths (solution given in time-domain) matches CW solution
@test isapprox(g1_DA_semiinf_TD(τ,[0.00001, 10.0], 1.0, 0.1, 10.0, N_quad = 200), g1_DA_semiinf_CW(τ, 1.0, 0.1, 10.0))

# check other optical properties and index of refraction
@test g1_DA_semiinf_TD(τ,[0.00001, 10.0], 0.5, 0.04, 14.2, n_med = 1.2, n_ext = 1.1, N_quad = 300) ≈ g1_DA_semiinf_CW(τ, 0.5, 0.04, 14.2, n_med = 1.2, n_ext = 1.1)

# test if a very narrow gate integrates close to the single point time domain solution 
@test isapprox(g1_DA_semiinf_TD(τ,[1.101, 1.102], 1.0, 0.1, 10.0), g1_DA_semiinf_TD(τ,1.1015, 1.0, 0.1, 10.0), rtol = 1e-7)

#check corresponding g2 functions (these call g1 functions)
@test g2_DA_semiinf_TD(τ,[0.00001, 10.0], 1.8, 0.32, 4.1, BFi = 1.2e-7, n_ext = 1.8, n_med = 1.4, z = 0.0, λ = 820.1, β = 0.88, N_quad = 200) ≈ g2_DA_semiinf_CW(τ, 1.8, 0.32, 4.1, BFi = 1.2e-7, n_ext = 1.8, n_med = 1.4, z = 0.0, λ = 820.1, β = 0.88)
@test isapprox(g2_DA_semiinf_TD(τ,[2.101, 2.102], 0.4, 0.01, 12.0), g2_DA_semiinf_TD(τ,2.1015, 0.4, 0.01, 12.0), rtol = 1e-7)


ρ = 1.0; μa = 0.1; μsp = 10.0; n_med = 1.0; n_ext = 1.0
β = 1.0; λ = 700.0; z = 0.0
BFi = 2.0e-8

si = g2_DA_semiinf_CW(τ, ρ, μa, μsp, BFi = BFi)

### test layered against semi-inf

μa = [0.1, 0.1]; μsp = [10.0, 10.0]; n_med = [1.0, 1.0]; n_ext = 1.0
BFi = [2.0e-8, 2.0e-8]; l = [1.0, 10.0]; a = 25.0

data = Nlayer_cylinder_DCS(ρ = ρ, μa = μa, μsp = μsp, n_med = n_med, n_ext = n_ext, β = β, λ = λ, BFi = BFi, z = z, a = a, l = l, bessels = besselroots[1:2000])


@test g2_DA_Nlay_cylinder_CW(τ, data) ≈ g2_DA_Nlay_cylinder_CW(τ, ρ, μa, μsp; BFi = BFi, β = β, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, bessels = besselroots[1:2000])

layered = g2_DA_Nlay_cylinder_CW(τ, data)

@test si ≈ layered

end
