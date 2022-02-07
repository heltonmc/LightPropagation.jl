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

### test layered against semi-inf

ρ = 1.0; μa = 0.1; μsp = 10.0; n_med = 1.0; n_ext = 1.0
β = 1.0; λ = 700.0; z = 0.0
BFi = 2.0e-8

si = g2_DA_semiinf_CW(τ, ρ, μa, μsp, BFi = BFi)
μa = (0.1, 0.1); μsp = (10.0, 10.0); n_med = (1.0, 1.0); n_ext = 1.0
BFi = (2.0e-8, 2.0e-8); l = (1.0, 10.0); a = 25.0

@test si ≈ g2_DA_Nlay_cylinder_CW(τ, ρ, μa, μsp; BFi = BFi, β = β, n_ext = n_ext, n_med = n_med, l = l, a = a, z = z, λ = λ, N_J0Roots = 2000)



# test layered TD solution for single time 

τ = 10 .^(range(-10,stop=0,length=250))
ρ = 1.0; mua = 0.1; musp = 10.0; BFi = 2.0e-8; n_ext = 1.0; n_med = 1.0; z = 0.0; λ = 700.0
t = 1.0
si = g1_DA_semiinf_TD(τ, t, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
layered = g1_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 10, N_laplace = 12)

@test si ≈ layered

# test for different optical properties and SDS (limited by degeneration to semi-infinite solution have adjusted relative error)
τ = 10 .^(range(-10,stop=0,length=250))
ρ = 2.4; mua = 0.18; musp = 20.1; BFi = 1.2e-6; n_ext = 1.2; n_med = 1.1; z = 0.0; λ = 740.0
t = 2.1
si = g1_DA_semiinf_TD(τ, t, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
layered = g1_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 20, N_laplace = 12)

@test isapprox(si, layered, rtol = 1e-5)

# test layered TD solution for 


# test range of time vector matches 
τ = 10 .^(range(-10,stop=0,length=250))
ρ = 1.0; mua = 0.1; musp = 10.0; BFi = 2.0e-8; n_ext = 1.0; n_med = 1.0; z = 0.0; λ = 700.0
t = [1.0, 1.5]
si = g1_DA_semiinf_TD(τ, t, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
layered = g1_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 200, N_laplace = 8, N_quad = 100)
@test isapprox(si, layered, rtol = 1e-5)

# test range of time vector matches different optical properties
τ = 10 .^(range(-10,stop=0,length=250))
ρ = 2.2; mua = 0.42; musp = 18.2; BFi = 8.2e-8; n_ext = 1.2; n_med = 1.4; z = 0.0; λ = 700.0
t = [1.0, 2.6]
si = g1_DA_semiinf_TD(τ, t, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)
layered = g1_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 200, N_laplace = 16, N_quad = 100)
@test isapprox(si, layered, rtol = 1e-5)

# test that the integral over whole time window range matches cw solution
τ = 10 .^(range(-10,stop=0,length=250))
ρ = 1.0; mua = 0.2; musp = 10.0; BFi = 3.0e-8; n_ext = 1.2; n_med = 1.4; z = 0.0; λ = 700.0
t = [1e-5, 10.0]
layered = g1_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 700, N_laplace = 54, N_quad = 100)
si = g1_DA_semiinf_CW(τ, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)

@test isapprox(si, layered, rtol = 1e-5)

# test that g2 is working

τ = 10 .^(range(-10,stop=0,length=250))
ρ = 1.0; mua = 0.2; musp = 10.0; BFi = 3.0e-8; n_ext = 1.2; n_med = 1.4; z = 0.0; λ = 700.0
t = [1e-5, 10.0]
layered = g2_DA_Nlay_cylinder_TD(τ, t, ρ, (mua, mua), (musp, musp); BFi = (BFi, BFi), n_ext = n_ext, n_med = (n_med, n_med), z = z, λ = λ, N_J0Roots = 700, N_laplace = 54, N_quad = 100)
si = g2_DA_semiinf_CW(τ, ρ, mua, musp; BFi = BFi, n_ext = n_ext, n_med = n_med, z = z, λ = λ)

@test isapprox(si, layered, rtol = 1e-5)

# test Structures
data = Nlayer_cylinder_DCS(BFi = (2e-8, 2e-8), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 10.0), a = 20.0, z = 0.0, λ = 700.0, N_J0Roots = 500, β = 1.0)
@test g2_DA_Nlay_cylinder_CW(τ, data) ≈ g2_DA_Nlay_cylinder_CW(τ, data.ρ, data.μa, data.μsp)
@test g2_DA_Nlay_cylinder_TD(τ,1.0, data) ≈ g2_DA_Nlay_cylinder_TD(τ,1.0, data.ρ, data.μa, data.μsp)

end
