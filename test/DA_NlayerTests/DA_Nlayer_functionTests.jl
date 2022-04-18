module DA_Nlayer_functionTests

using Test
using LightPropagation

# CW comparison of N-layer cylinder to semi-infinite
# note that for large layer thicknesses and radii many more besselroots are needed

### test at multiple SDS

# ρ = 0.0 cm
ρ, μa, μsp = (0.0, 0.1, 10.0)
@test fluence_DA_semiinf_CW(ρ, μa, μsp) ≈ fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0))
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-6), atol=1.0e-6)
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-3), atol=1.0e-3)

# ρ = 1.0 cm
ρ, μa, μsp = (1.0, 0.1, 10.0)
@test fluence_DA_semiinf_CW(ρ, μa, μsp) ≈ fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0))
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-6), atol=1.0e-6)
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-3), atol=1.0e-3)

# ρ = 4.0 cm
ρ, μa, μsp = (4.0, 0.1, 10.0)
@test fluence_DA_semiinf_CW(ρ, μa, μsp) ≈ fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0))
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-6), atol=1.0e-6)
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), atol=1.0e-3), atol=1.0e-3)

# ρ = 8.0 cm; for long SDS you need more besselroots or decrease the dimensions of cylinder
ρ, μa, μsp = (8.0, 0.1, 10.0)
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), a=20.0), rtol=1.0e-6)
@test isapprox(fluence_DA_semiinf_CW(ρ, μa, μsp), fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), a=20.0, atol=1.0e-8), atol=1.0e-9)


### test that 2, 3, 4, 8 layers equal with same optical properties
ρ, μa, μsp = (1.0, 0.1, 10.0)
cyl_2 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (4.0, 8.0), n_med=(1.0, 1.0), a=10.0)
cyl_3 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa), (μsp, μsp, μsp), l = (4.0, 4.0, 4.0), n_med=(1.0, 1.0, 1.0), a=10.0)
cyl_4 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa, μa), (μsp, μsp, μsp, μsp), l = (4.0, 4.0, 4.0, 4.0), n_med=(1.0, 1.0, 1.0, 1.0), a=10.0)
cyl_8 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa, μa, μa, μa, μa, μa), (μsp, μsp, μsp, μsp, μsp, μsp, μsp, μsp), l = (4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0), n_med=(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), a=10.0)

@test cyl_2 ≈ cyl_3 ≈ cyl_4 ≈ cyl_8

### test high scatterring
# SDS = 1.0, μsp = 60.0, once again need more besselroots if you want higher scattering coefficients
ρ, μa, μsp = (1.0, 0.1, 60.0)
@test fluence_DA_semiinf_CW(ρ, μa, μsp) ≈ fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa, μa), (μsp, μsp, μsp, μsp), l = (2.0, 2.0, 2.0, 2.0), n_med=(1.0, 1.0, 1.0, 1.0), a=10.0)



### test for varying index of refractions
ρ, μa, μsp = (1.0, 0.1, 10.0)
cyl_2 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (4.0, 8.0), n_med=(1.5, 1.5), n_ext=1.4, a=10.0)
cyl_3 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa), (μsp, μsp, μsp), l = (4.0, 4.0, 4.0), n_med=(1.5, 1.5, 1.5), n_ext=1.4, a=10.0)
cyl_4 = fluence_DA_Nlay_cylinder_CW(ρ, (μa, μa, μa, μa), (μsp, μsp, μsp, μsp), l = (4.0, 4.0, 4.0, 4.0), n_med=(1.5, 1.5, 1.5, 1.5), n_ext=1.4, a=10.0)
@test cyl_2 ≈ cyl_3 ≈ cyl_4 ≈ fluence_DA_semiinf_CW(ρ, μa, μsp, n_med = 1.5, n_ext = 1.4)

#### tests for CW fluence in the bottom layer GN
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1), (10.0, 10.0), l = (1.0, 1.0), z = 2.0) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(1.0, (0.1, 0.1), (10.0, 10.0), l = (1.0, 1.0), z = 2.0) ≈ fluence_DA_slab_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(1.6, (0.21, 0.21), (12.1, 12.1), l = (1.0, 1.0), z = 2.0) ≈ fluence_DA_slab_CW(1.6, 0.21, 12.1; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(2.2, (0.31, 0.31, 0.31), (12.1, 12.1, 12.1), l = (1.0, 1.0, 1.0), n_med = (1.0, 1.0, 1.0), z = 3.0) ≈ fluence_DA_slab_CW(2.2, 0.31, 12.1; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(2.8, (0.41, 0.41, 0.41, 0.41), (13.1, 13.1, 13.1, 13.1), l = (1.0, 1.0, 1.0, 1.0), n_med = (1.0, 1.0, 1.0, 1.0), z = 4.0) ≈ fluence_DA_slab_CW(2.8, 0.41, 13.1; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1), (10.0, 10.0), l = (1.5, 1.5), z = 3.0) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1, 0.1, 0.1, 0.1, 0.1), (10.0, 10.0, 10.0, 10.0, 10.0, 10.0), n_med = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), l = (0.5, 0.5, 0.5, 0.5, 0.5, 0.5), z = 3.0) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1, 0.1, 0.1, 0.1, 0.1), (10.0, 10.0, 10.0, 10.0, 10.0, 10.0), n_med = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), l = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), z = 6.0) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 6.0, z = 6.0, xs = 50)

### test that scalar and vector ρ matches
μsp = (10.0, 10.0, 10.0, 10.0); μa = (0.1, 0.1, 0.1, 0.1); n_ext = 1.0;
n_med = (1.0, 1.0, 1.0, 1.0); l = (0.5, 0.8, 1.0, 5.0); a = 5.0; z = 0.0
                                        
ρ1 = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
@test [fluence_DA_Nlay_cylinder_CW(ρ1, μa, μsp, n_med=n_med, l=l, a=a)...] ≈ map(ρ -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_med=n_med, l=l, a=a), 1.0:0.1:2.0)

### test approximations against higher precision

# test for ρ = 1.0 cm
#fluence_DA_Nlay_cylinder_CW(big"1.0", BigFloat.((0.2, 0.1, 0.2)), BigFloat.((12.0, 10.0, 11.0)), l = BigFloat.((1.0, 1.2, 4.0)), n_med = BigFloat.((1.0, 1.0, 1.0)), z=big"0.0", MaxIter=100000, atol=1.0e-60, a = big"30.0") ...
# = 0.01059246710599644995407690834032245361228582286358680220137244538434879662408047
exact = fluence_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = fluence_DA_Nlay_cylinder_CW_approx(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ 0.01059246710599644995407690834032245361228582286358680220137244538434879662408047 
@test approx ≈ 0.01059246710599644995407690834032245361228582286358680220137244538434879662408047

# test for ρ = 3.5 cm
#fluence_DA_Nlay_cylinder_CW(big"3.5", BigFloat.((0.01, 0.2, 0.35)), BigFloat.((6.0, 15.0, 11.0)), l = BigFloat.((0.5, 1.2, 4.0)), n_med = BigFloat.((1.0, 1.0, 1.0)), z=big"0.0", MaxIter=100000, atol=1.0e-60, a = big"30.0") ...
# = 2.37808800608361746480719102898448922930259853687958126602277737309411249063373e-05
exact = fluence_DA_Nlay_cylinder_CW(3.5, (0.01, 0.2, 0.35), (6.0, 15.0, 11.0), l = (0.5, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = fluence_DA_Nlay_cylinder_CW_approx(3.5, (0.01, 0.2, 0.35), (6.0, 15.0, 11.0), l = (0.5, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ 2.37808800608361746480719102898448922930259853687958126602277737309411249063373e-05
@test approx ≈ 2.37808800608361746480719102898448922930259853687958126602277737309411249063373e-05

# test for ρ = 10.0 cm
# fluence_DA_Nlay_cylinder_CW(big"10.0", BigFloat.((0.2, 0.1, 0.2)), BigFloat.((12.0, 10.0, 11.0)), l = BigFloat.((1.0, 1.2, 4.0)), n_med = BigFloat.((1.0, 1.0, 1.0)), z=big"0.0", MaxIter=100000, atol=1.0e-60, a = big"30.0") ...
# = 5.130032007256222345294315223040138886305142618609630358861997082580778793718697e-13
#exact = fluence_DA_Nlay_cylinder_CW(10.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = fluence_DA_Nlay_cylinder_CW_approx(10.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
#@test exact ≈ 5.130032007256222345294315223040138886305142618609630358861997082580778793718697e-13 this test fails because of numerical errors
@test approx ≈ 5.130032007256222345294315223040138886305142618609630358861997082580778793718697e-13 # the approximation is actually more accurate

##################
# Time-Domain
##################

## test reflectance in first layer
# test long times at SDS = 0.5 cm
t = range(0.03, 8.0, length = 60)

ρ, μa, μsp = (0.5, 0.1, 10.0)
si = fluence_DA_semiinf_TD(t, ρ, μa, μsp)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 72, MaxIter=500, atol=1.0e-16)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 96, MaxIter=500, atol=1.0e-16)
@test a72 ≈ si
@test a96 ≈ si

# test long times at SDS = 3.0
ρ, μa, μsp = (3.0, 0.1, 10.0)
si = fluence_DA_semiinf_TD(t, ρ, μa, μsp)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 96)
@test a72 ≈ si
@test a96 ≈ si


# test high scattering at SDS = 1.5 cm
ρ, μa, μsp = (1.5, 0.1, 60.0)
si = fluence_DA_semiinf_TD(t, ρ, μa, μsp)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 72, MaxIter=500, atol=1.0e-12)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(4.0, 5.0), N = 96, MaxIter=5000)
@test a72 ≈ si
@test a96 ≈ si

## in general the error in the TD is like 1e-14..

## test transmittance in bottom layer
ρ, μa, μsp = (0.5, 0.1, 10.0)
slab = fluence_DA_slab_TD(t, ρ, μa, μsp; s = 2.0, z = 2.0, xs = 50)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.0, 1.0), z = 2.0, N = 72, MaxIter=100)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.0, 1.0), z = 2.0, N = 96, MaxIter=100)
@test a72 ≈ slab
@test a96 ≈ slab

ρ, μa, μsp = (1.5, 0.1, 10.0)
slab = fluence_DA_slab_TD(t, ρ, μa, μsp; s = 3.0, z = 3.0, xs = 50)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.5, 1.5), z = 3.0, N = 72, MaxIter=100)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.5, 1.5), z = 3.0, N = 96, MaxIter=100)
@test a72 ≈ slab
@test a96 ≈ slab


ρ, μa, μsp = (0.5, 0.1, 10.0)
slab = fluence_DA_slab_TD(t, ρ, μa, μsp; s = 3.0, z = 3.0, xs = 50, n_ext = 1.3, n_med = 1.2)
a72 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.5, 1.5), z = 3.0, n_ext = 1.3, n_med = (1.2, 1.2), N = 72, MaxIter=100)
a96 = fluence_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l=(1.5, 1.5), z = 3.0, n_ext = 1.3, n_med = (1.2, 1.2), N = 96, MaxIter=100)
@test a72 ≈ slab
@test a96 ≈ slab

# approx
# TD
t = 1.0
exact = fluence_DA_Nlay_cylinder_TD(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = fluence_DA_Nlay_cylinder_TD_approx(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ approx

t = 0.5:0.1:1.5
exact = fluence_DA_Nlay_cylinder_TD(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = fluence_DA_Nlay_cylinder_TD_approx(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ approx


## frequency domain
ρ, μa, μsp = (0.0, 0.1, 10.0)
@test LightPropagation.fluence_DA_semiinf_FD(ρ, μa, μsp, ω = 1.4) ≈ LightPropagation.fluence_DA_Nlay_cylinder_FD(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0), ω = 1.4)


### Flux tests

# CW 
ρ, μa, μsp = (1.0, 0.1, 10.0)
@test flux_DA_semiinf_CW(ρ, μa, μsp) ≈ flux_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0))
@test flux_DA_slab_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 15) ≈ flux_DA_Nlay_cylinder_CW(ρ, (μa, μa), (μsp, μsp), l = (2.0, 2.0), z = 4.0)

exact = flux_DA_Nlay_cylinder_CW(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = flux_DA_Nlay_cylinder_CW_approx(1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ approx

# TD
t = 1.0
ρ, μa, μsp = (1.0, 0.1, 10.0)
@test flux_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0)) ≈ flux_DA_semiinf_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)
@test flux_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l = (2.0, 2.0), z = 4.0) ≈ flux_DA_slab_TD(t, 1.0, 0.1, 10.0; s = 4.0, z = 4.0, xs = 15)
exact = flux_DA_Nlay_cylinder_TD(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = flux_DA_Nlay_cylinder_TD_approx(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ approx

t = 0.5:0.1:1.5
@test flux_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l = (10.0, 10.0)) ≈ flux_DA_semiinf_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)
@test flux_DA_Nlay_cylinder_TD(t, ρ, (μa, μa), (μsp, μsp), l = (2.0, 2.0), z = 4.0) ≈ flux_DA_slab_TD(t, 1.0, 0.1, 10.0; s = 4.0, z = 4.0, xs = 15)
exact = flux_DA_Nlay_cylinder_TD(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
approx = flux_DA_Nlay_cylinder_TD_approx(t, 1.0, (0.2, 0.1, 0.2), (12.0, 10.0, 11.0), l = (1.0, 1.2, 4.0), n_med = (1.0, 1.0, 1.0), MaxIter=100000, atol=1.0e-60, a = 30.0, z = 0.0)
@test exact ≈ approx

end # module
