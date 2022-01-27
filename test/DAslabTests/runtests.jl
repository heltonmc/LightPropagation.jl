module DAslabTests

using Test
using ForwardDiff
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD
using LightPropagation: fluence_DA_slab_CW, fluence_DA_slab_TD
using LightPropagation: flux_DA_slab_CW, flux_DA_slab_TD, flux_DA_semiinf_TD, flux_DA_semiinf_CW


# Test slabs against semi-inf for sufficiently thick slabs

# CW
@test fluence_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 20.0, z = 0.0) == 0.024035998288332954
@test fluence_DA_slab_CW(2.4, 0.6, 18.0, n_ext = 1.3, n_med = 1.2, s = 15.0, z = 1.0) == 1.77784095787581e-7
@test fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 20.0)
@test fluence_DA_semiinf_CW(4.0, 0.01, 50.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_CW(4.0, 0.01, 50.0,n_ext = 1.0, n_med = 1.0, z = 0.0, s = 20.0)

# TD
t = 0.1:0.5:8
@test fluence_DA_semiinf_TD(t, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 40.0)
@test fluence_DA_semiinf_TD(t, 1.0, 0.01, 70.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_TD(t, 1.0, 0.01, 70.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 40.0)

@test fluence_DA_semiinf_TD(t, 3.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_TD(t, 3.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 40.0)
@test fluence_DA_semiinf_TD(t, 3.0, 0.01, 70.0, n_ext = 1.0, n_med = 1.0, z = 0.0) ≈ fluence_DA_slab_TD(t, 3.0, 0.01, 70.0, n_ext = 1.0, n_med = 1.0, z = 0.0, s = 40.0)
@test fluence_DA_semiinf_TD(t, 3.0, 0.01, 70.0, n_ext = 1.4, n_med = 1.2, z = 0.1) ≈ fluence_DA_slab_TD(t, 3.0, 0.01, 70.0, n_ext = 1.4, n_med = 1.2, z = 0.1, s = 40.0)

## check flux

# reflectance

D = 1 / (3 * 10.0)
@test (D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 1.0, z = dz), 0.0)) ≈ flux_DA_slab_TD(1.0, 1.0, 0.1, 10.0,  n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0)

D = 1 / (3 * 20.0)
@test (D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(1.0, 1.3, 0.3, 20.0, n_ext = 1.0, n_med = 1.0, s = 1.0, z = dz), 0.0)) ≈ flux_DA_slab_TD(1.0, 1.3, 0.3, 20.0,  n_ext = 1.0, n_med = 1.0, s = 1.0, z = 0.0)

# transmittance
D = 1 / (3 * 20.0)
@test (-D * ForwardDiff.derivative(dz -> fluence_DA_slab_TD(1.0, 1.3, 0.3, 20.0, n_ext = 1.0, n_med = 1.0, s = 1.0, z = dz), 1.0)) ≈ flux_DA_slab_TD(1.0, 1.3, 0.3, 20.0,  n_ext = 1.0, n_med = 1.0, s = 1.0, z = 1.0)

## CW
@test flux_DA_slab_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0, s = 40.0, z = 0.0, xs = 15) ≈ flux_DA_semiinf_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)
@test flux_DA_slab_CW(1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, s = 2.0, z = 0.0,xs = 10) ≈ flux_DA_semiinf_CW(1.0, 0.1, 10.0; n_ext = 1.0, n_med = 1.0)

end 