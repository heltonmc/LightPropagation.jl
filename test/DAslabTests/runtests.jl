module DAslabTests

using Test
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD
using LightPropagation: fluence_DA_slab_CW, fluence_DA_slab_TD

# Test slabs against semi-inf for sufficiently thick slabs

# CW
@test fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0) ≈ fluence_DA_slab_CW(1.0, [0.1, 10.0], 1.0, 1.0, 20.0, 0.0)
@test fluence_DA_semiinf_CW(4.0, [0.01, 50.0], 1.0, 1.0, 0.0) ≈ fluence_DA_slab_CW(4.0, [0.01, 50.0], 1.0, 1.0, 20.0, 0.0)

# TD
t = 0.1:0.5:8
@test fluence_DA_semiinf_TD(t, [0.1, 10.0], 1.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, [0.1, 10.0], 1.0, 1.0, 1.0, 40.0, 0.0)
@test fluence_DA_semiinf_TD(t, [0.01, 70.0], 1.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, [0.01, 70.0], 1.0, 1.0, 1.0, 40.0, 0.0)

@test fluence_DA_semiinf_TD(t, [0.1, 10.0], 3.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, [0.1, 10.0], 3.0, 1.0, 1.0, 40.0, 0.0)
@test fluence_DA_semiinf_TD(t, [0.01, 70.0], 3.0, 1.0, 1.0, 0.0) ≈ fluence_DA_slab_TD(t, [0.01, 70.0], 3.0, 1.0, 1.0, 40.0, 0.0)


end # module