module DAslabTests

using Test
using ForwardDiff
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD
using LightPropagation: fluence_DA_slab_CW, fluence_DA_slab_TD
using LightPropagation: refl_DA_slab_TD, trans_DA_slab_TD


# Test slabs against semi-inf for sufficiently thick slabs

# CW
@test fluence_DA_slab_CW(1.0, 0.1, 10.0, 1.0, 1.0, 20.0, 0.0) == 0.024035998288332954
@test fluence_DA_slab_CW(2.4, 0.6, 18.0, 1.3, 1.2, 15.0, 1.0) == 1.77784095787581e-7
@test fluence_DA_semiinf_CW(1.0, 0.1, 10.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_CW(1.0, 0.1, 10.0, 1.0, 1.0, 20.0, 0.0)
@test fluence_DA_semiinf_CW(4.0, 0.01, 50.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_CW(4.0, 0.01, 50.0, 1.0, 1.0, 20.0, 0.0)

# TD
t = 0.1:0.5:8
@test fluence_DA_semiinf_TD(t, 1.0, 0.1, 10.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, 1.0, 1.0, 40.0, 0.0)
@test fluence_DA_semiinf_TD(t, 1.0, 0.01, 70.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, 1.0, 0.01, 70.0, 1.0, 1.0, 40.0, 0.0)

@test fluence_DA_semiinf_TD(t, 3.0, 0.1, 10.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, 3.0, 0.1, 10.0, 1.0, 1.0, 40.0, 0.0)
@test fluence_DA_semiinf_TD(t, 3.0, 0.01, 70.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, 3.0, 0.01, 70.0, 1.0, 1.0, 40.0, 0.0)

## check flux

# reflectance
D = 1 / (3 * 10.0)
@test (D * ForwardDiff.derivative(z -> fluence_DA_slab_TD(1.0, 1.0, 0.1, 10.0, 1.0, 1.0, 1.0, z), 0.0)) ≈ refl_DA_slab_TD(1.0, [0.1, 10.0], 1.0, 1.0, 1.0, 1.0)[1]

D = 1 / (3 * 20.0)
@test (D * ForwardDiff.derivative(z -> fluence_DA_slab_TD(1.5, 1.0, 0.4, 20.0, 1.0, 1.0, 1.0, z), 0.0)) ≈ refl_DA_slab_TD(1.5, [0.4, 20.0], 1.0, 1.0, 1.0, 1.0)[1]

# transmittance

D = 1 / (3 * 20.0)
@test (-D * ForwardDiff.derivative(z -> fluence_DA_slab_TD(1.5, 1.0, 0.4, 20.0, 1.0, 1.0, 1.0, z), 1.0)) ≈ trans_DA_slab_TD(1.5, [0.4, 20.0], 1.0, 1.0, 1.0, 1.0)[1]

end 