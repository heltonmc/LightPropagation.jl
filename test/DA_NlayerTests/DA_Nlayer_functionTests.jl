module DA_Nlayer_functionTests

using Test
using LightPropagation

# CW comparison of N-layer cylinder to semi-infinite
# note that for large layer thicknesses and radii many more besselroots are needed

### test at multiple SDS

# ρ = 0.0 cm
cylinder_data = Nlayer_cylinder(ρ = 0.0, a = 100.0, l = (10.0, 10.0, 10.0, 10.0), N_J0Roots = 10000)
@test fluence_DA_semiinf_CW(0.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data)

# ρ = 1.0 cm
cylinder_data = Nlayer_cylinder(ρ = 1.0, a = 100.0, l = (10.0, 10.0, 10.0, 10.0), N_J0Roots = 10000)
@test fluence_DA_semiinf_CW(1.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data)

# ρ = 4.0 cm
cylinder_data = Nlayer_cylinder(ρ = 4.0, a = 100.0, l = (10.0, 10.0, 10.0, 10.0), N_J0Roots = 10000)
@test fluence_DA_semiinf_CW(4.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data)

# ρ = 8.0 cm; for long SDS you need more besselroots or decrease the dimensions of cylinder
cylinder_data = Nlayer_cylinder(ρ = 8.0, a = 100.0, l = (10.0, 10.0, 10.0, 10.0),N_J0Roots = 10000)
@test isapprox(fluence_DA_semiinf_CW(8.0, 0.1, 10.0), fluence_DA_Nlay_cylinder_CW(cylinder_data), rtol=1e-6)


### test that 2, 3, 4 layers equal with same optical properties
cylinder_4 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1, 0.1), l = (10.0, 10.0, 10.0, 10.0), n_med = (1.0, 1.0, 1.0, 1.0), N_J0Roots = 10000)
cylinder_3 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1), l = (10.0, 10.0, 10.0), n_med = (1.0, 1.0, 1.0), N_J0Roots = 10000)
cylinder_2 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (10.0, 10.0), n_med = (1.0, 1.0), N_J0Roots = 10000)
@test fluence_DA_Nlay_cylinder_CW(cylinder_4) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_3)  ≈ fluence_DA_Nlay_cylinder_CW(cylinder_2)

### test that arbitrary layers are equal
cylinder_5 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1, 0.1, 0.1), l = (2.0, 2.0, 2.0, 5.0, 5.0), n_med = (1.0, 1.0, 1.0, 1.0, 1.0), N_J0Roots = 10000)
cylinder_8 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), l = (2.0, 2.0, 2.0, 5.0, 5.0, 3.0, 3.0, 3.0), n_med = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), N_J0Roots = 10000)

@test fluence_DA_Nlay_cylinder_CW(cylinder_4) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_5) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_8)
### test high scatterring
# SDS = 1.0, μsp = 60.0, once again need more besselroots if you want higher scattering coefficients
cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (60.0, 60.0, 60.0, 60.0), μa = (0.1, 0.1, 0.1, 0.1), l = (2.0, 2.0, 2.0, 2.0), n_med = (1.0, 1.0, 1.0, 1.0), a = 20.0, N_J0Roots = 10000)
@test fluence_DA_semiinf_CW(1.0, 0.1, 60.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data)

### test for varying index of refractions
cyl_2 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (10.0, 10.0), n_med = (1.5, 1.5), n_ext = 1.4, a = 20.0, N_J0Roots = 10000)
cyl_3 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1), l = (10.0, 10.0, 10.0), n_med = (1.5, 1.5, 1.5), n_ext = 1.4, a = 20.0, N_J0Roots = 10000)
cyl_4 = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0, 10.0, 10.0), μa = (0.1, 0.1, 0.1, 0.1), l = (10.0, 10.0, 10.0, 10.0), n_med = (1.5, 1.5, 1.5, 1.5), n_ext = 1.4, a = 20.0, N_J0Roots = 10000)
@test fluence_DA_Nlay_cylinder_CW(cyl_4) ≈ fluence_DA_Nlay_cylinder_CW(cyl_3) ≈ fluence_DA_Nlay_cylinder_CW(cyl_2) ≈ fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_med = 1.5, n_ext = 1.4)

#### tests for CW fluence in the bottom layer GN
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1), (10.0, 10.0), 1.0, (1.0, 1.0), (1.0, 1.0), 10.0, 2.0, 10000) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(1.0, (0.1, 0.1), (10.0, 10.0), 1.0, (1.0, 1.0), (1.0, 1.0), 10.0, 2.0, 10000) ≈ fluence_DA_slab_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(1.6, (0.21, 0.21), (12.1, 12.1), 1.0, (1.0, 1.0), (1.0, 1.0), 10.0, 2.0, 10000) ≈ fluence_DA_slab_CW(1.6, 0.21, 12.1; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(2.2, (0.31, 0.31, 0.31), (12.1, 12.1, 12.1), 1.0, (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), 10.0, 3.0, 10000) ≈ fluence_DA_slab_CW(2.2, 0.31, 12.1; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(2.8, (0.41, 0.41, 0.41, 0.41), (13.1, 13.1, 13.1, 13.1), 1.0, (1.0, 1.0, 1.0, 1.0), (1.0, 1.0, 1.0, 1.0), 10.0, 4.0, 10000) ≈ fluence_DA_slab_CW(2.8, 0.41, 13.1; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 50)

@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1), (10.0, 10.0), 1.0, (1.0, 1.0), (1.5, 1.5), 10.0, 3.0, 10000) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)

@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1, 0.1, 0.1, 0.1, 0.1), (10.0, 10.0, 10.0, 10.0, 10.0, 10.0), 1.0, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), (0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 10.0, 3.0, 10000) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50)
@test fluence_DA_Nlay_cylinder_CW(0.0, (0.1, 0.1, 0.1, 0.1, 0.1, 0.1), (10.0, 10.0, 10.0, 10.0, 10.0, 10.0), 1.0, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 10.0, 6.0, 10000) ≈ fluence_DA_slab_CW(0.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 6.0, z = 6.0, xs = 50)

### test that scalar and vector ρ matches (they call different routines)
μsp = (10.0, 10.0, 10.0, 10.0); μa = (0.1, 0.1, 0.1, 0.1); n_ext = 1.0;
n_med = (1.0, 1.0, 1.0, 1.0); l = (0.5, 0.8, 1.0, 5.0); a = 5.0; z = 0.0                                         
N_J0Roots = 1000                       
ρ1 = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
@test [fluence_DA_Nlay_cylinder_CW(ρ1, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots)...] ≈ map(ρ -> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots), 1.0:0.1:2.0)

#########
# Time-Domain
##########

## test reflectance in first layer
# test long times at SDS = 0.5 cm
t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 0.5, N_J0Roots = 600)

si = fluence_DA_semiinf_TD(t, 0.5, 0.1, 10.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ si
@test a96 ≈ si

# test long times at SDS = 3.0
t = range(0.04, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 3.0, l = (4.0, 4.0, 4.0, 5.0), N_J0Roots = 600)

si = fluence_DA_semiinf_TD(t, 3.0, 0.1, 10.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ si
@test a96 ≈ si

# test high scattering at SDS = 1.5 cm
t = range(0.04, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 1.5, l = (4.0, 4.0, 4.0, 5.0), μsp = (60.0, 60.0, 60.0, 60.0), N_J0Roots = 600)

si = fluence_DA_semiinf_TD(t, 1.5, 0.1, 60.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ si
@test a96 ≈ si

## in general the error in the TD is like 1e-14..

## test transmittance in bottom layer
t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(μa = (0.1, 0.1), μsp = (10.0, 10.0), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.0, 1.0), a = 8.0, ρ = 0.5, z = 2.0, N_J0Roots = 600)

slab = fluence_DA_slab_TD(t, 0.5, 0.1, 10.0; s = 2.0, z = 2.0, xs = 50)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ slab
@test a96 ≈ slab

t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(μa = (0.1, 0.1), μsp = (10.0, 10.0), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.5, 1.5), a = 8.0, ρ = 0.5, z = 3.0, N_J0Roots = 600)

slab = fluence_DA_slab_TD(t, 0.5, 0.1, 10.0; s = 3.0, z = 3.0, xs = 50)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ slab
@test a96 ≈ slab

t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(μa = (0.1, 0.1), μsp = (10.0, 10.0), n_ext = 1.0, n_med = (1.0, 1.0), l = (1.5, 1.5), a = 8.0, ρ = 1.5, z = 3.0, N_J0Roots = 600)

slab = fluence_DA_slab_TD(t, 1.5, 0.1, 10.0; s = 3.0, z = 3.0, xs = 50)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ slab
@test a96 ≈ slab

t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(μa = (0.1, 0.1), μsp = (10.0, 10.0), n_ext = 1.3, n_med = (1.2, 1.2), l = (1.0, 1.0), a = 8.0, ρ = 0.5, z = 2.0, N_J0Roots = 600)

slab = fluence_DA_slab_TD(t, 0.5, 0.1, 10.0; s = 2.0, z = 2.0, xs = 50, n_ext = 1.3, n_med = 1.2)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, N = 96)

@test a72 ≈ slab
@test a96 ≈ slab

### Flux tests

# CW 
cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (10.0, 10.0), n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_CW(cylinder_data) ≈ flux_DA_semiinf_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)

cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (2.0, 2.0), z = 4.0, n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_CW(cylinder_data) ≈ flux_DA_slab_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 15)

# TD
t = 1.0
cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (10.0, 10.0), n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_TD(t, cylinder_data) ≈ flux_DA_semiinf_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)

cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (2.0, 2.0), z = 4.0, n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_TD(t, cylinder_data) ≈ flux_DA_slab_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 15)

t = 0.5:0.1:1.5
cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (10.0, 10.0), n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_TD(t, cylinder_data) ≈ flux_DA_semiinf_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)

cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = (10.0, 10.0), μa = (0.1, 0.1), l = (2.0, 2.0), z = 4.0, n_med = (1.0, 1.0), n_ext = 1.0, a = 40.0, N_J0Roots = 10000)
@test flux_DA_Nlay_cylinder_TD(t, cylinder_data) ≈ flux_DA_slab_TD(t, 1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 15)

end # module
