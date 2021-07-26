module DA_Nlayer_functionTests

using Test
using LightPropagation

# CW comparison of N-layer cylinder to semi-infinite
# note that for large layer thicknesses and radii many more besselroots are needed

### test at multiple SDS

# ρ = 0.0 cm
cylinder_data = Nlayer_cylinder(ρ = 0.0, a = 100.0, l = [10.0, 10.0, 10.0, 10.0])
@test fluence_DA_semiinf_CW(0.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data, besselroots)

# ρ = 1.0 cm
cylinder_data = Nlayer_cylinder(ρ = 1.0, a = 100.0, l = [10.0, 10.0, 10.0, 10.0])
@test fluence_DA_semiinf_CW(1.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data, besselroots)

# ρ = 4.0 cm
cylinder_data = Nlayer_cylinder(ρ = 4.0, a = 100.0, l = [10.0, 10.0, 10.0, 10.0])
@test fluence_DA_semiinf_CW(4.0, 0.1, 10.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data, besselroots)

# ρ = 8.0 cm; for long SDS you need more besselroots or decrease the dimensions of cylinder
cylinder_data = Nlayer_cylinder(ρ = 8.0, a = 100.0, l = [10.0, 10.0, 10.0, 10.0])
@test isapprox(fluence_DA_semiinf_CW(8.0, 0.1, 10.0), fluence_DA_Nlay_cylinder_CW(cylinder_data, besselroots), rtol=1e-6)


### test that 2, 3, 4 layers equal with same optical properties
cylinder_4 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1, 0.1], l = [10.0, 10.0, 10.0, 10.0], n_med = [1.0, 1.0, 1.0, 1.0])
cylinder_3 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1], l = [10.0, 10.0, 10.0], n_med = [1.0, 1.0, 1.0])
cylinder_2 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0], μa = [0.1, 0.1], l = [10.0, 10.0], n_med = [1.0, 1.0])
@test fluence_DA_Nlay_cylinder_CW(cylinder_4, besselroots) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_3, besselroots)  ≈ fluence_DA_Nlay_cylinder_CW(cylinder_2, besselroots)

### test that arbitrary layers are equal
cylinder_5 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1, 0.1, 0.1], l = [2.0, 2.0, 2.0, 5.0, 5.0], n_med = [1.0, 1.0, 1.0, 1.0, 1.0])
cylinder_8 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], l = [2.0, 2.0, 2.0, 5.0, 5.0, 3.0, 3.0, 3.0], n_med = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

@test fluence_DA_Nlay_cylinder_CW(cylinder_4, besselroots) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_5, besselroots) ≈  fluence_DA_Nlay_cylinder_CW(cylinder_8, besselroots)
### test high scatterring
# SDS = 1.0, μsp = 60.0, once again need more besselroots if you want higher scattering coefficients
cylinder_data = Nlayer_cylinder(ρ = 1.0, μsp = [60.0, 60.0, 60.0, 60.0], μa = [0.1, 0.1, 0.1, 0.1], l = [2.0, 2.0, 2.0, 2.0], n_med = [1.0, 1.0, 1.0, 1.0], a = 20.0)
@test fluence_DA_semiinf_CW(1.0, 0.1, 60.0) ≈ fluence_DA_Nlay_cylinder_CW(cylinder_data, besselroots)

### test for varying index of refractions
cyl_2 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0], μa = [0.1, 0.1], l = [10.0, 10.0], n_med = [1.5, 1.5], n_ext = 1.4, a = 20.0)
cyl_3 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1], l = [10.0, 10.0, 10.0], n_med = [1.5, 1.5, 1.5], n_ext = 1.4, a = 20.0)
cyl_4 = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 10.0, 10.0, 10.0], μa = [0.1, 0.1, 0.1, 0.1], l = [10.0, 10.0, 10.0, 10.0], n_med = [1.5, 1.5, 1.5, 1.5], n_ext = 1.4, a = 20.0)
@test fluence_DA_Nlay_cylinder_CW(cyl_4, besselroots) ≈ fluence_DA_Nlay_cylinder_CW(cyl_3, besselroots) ≈ fluence_DA_Nlay_cylinder_CW(cyl_2, besselroots) ≈ fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_med = 1.5, n_ext = 1.4)


# TD

# test long times at SDS = 0.5 cm
t = range(0.03, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 0.5)

si = fluence_DA_semiinf_TD(t, [0.1, 10.0], 0.5, 1.0, 1.0, 0.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 96)

@test a72 ≈ si
@test a96 ≈ si

# test long times at SDS = 3.0
t = range(0.04, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 3.0, l = [4.0, 4.0, 4.0, 5.0])

si = fluence_DA_semiinf_TD(t, [0.1, 10.0], 3.0, 1.0, 1.0, 0.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 96)

@test a72 ≈ si
@test a96 ≈ si

# test high scattering at SDS = 1.5 cm
t = range(0.04, 8.0, length = 60)
cylinder_data = Nlayer_cylinder(a = 8.0, ρ = 1.5, l = [4.0, 4.0, 4.0, 5.0], μsp = [60.0, 60.0, 60.0, 60.0])

si = fluence_DA_semiinf_TD(t, [0.1, 60.0], 1.5, 1.0, 1.0, 0.0)

a72 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 72)
a96 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 96)

@test a72 ≈ si
@test a96 ≈ si

## in general the error in the TD is like 1e-14..

end # module



fluence([10.0, 10.0, 10.0, 10.0], [0.1, 0.1, 0.1, 0.1], 1.0, [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], 0.0, 10.0, 4.0, besselroots)
8.077937297438846e-5  == fluence_DA_slab_CW(0.0, big(0.1), 10.0; n_ext = 1.0, n_med = 1.0, s = 4.0, z = 4.0, xs = 50)

fluence_DA_slab_CW(0.0, big(0.1), 10.0; n_ext = 1.0, n_med = 1.0, s = 3.0, z = 3.0, xs = 50) == fluence([10.0, 10.0], [0.1, 0.1], 1.0, [1.0, 1.0], [1.5, 1.5], 0.0, 10.0, 3.0, besselroots)
fluence_DA_slab_CW(0.0, big(0.1), 10.0; n_ext = 1.0, n_med = 1.0, s = 2.0, z = 2.0, xs = 50) == fluence([10.0, 10.0], [0.1, 0.1], 1.0, [1.0, 1.0], [1.0, 1.0], 0.0, 10.0, 2.0, besselroots)
fluence_DA_slab_CW(0.0, big(0.1), 10.0; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 1.0, xs = 50) == fluence([10.0, 10.0], [0.1, 0.1], 1.0, [1.0, 1.0], [0.5, 0.5], 0.0, 10.0, 1.0, besselroots)
fluence_DA_slab_CW(1.0, big(0.1), 10.0; n_ext = 1.0, n_med = 1.0, s = 1.0, z = 1.0, xs = 50) == fluence([10.0, 10.0], [0.1, 0.1], 1.0, [1.0, 1.0], [0.5, 0.5], 1.0, 10.0, 1.0, besselroots)
### add some more tests for N layers in transmittance and reflectance in the first layer ( like N = 12)