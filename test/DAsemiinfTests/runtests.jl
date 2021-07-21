module DAsemiinfTests

using Test
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD


# Hard code DA semiinf solution to compare other geometries

# CW
@test fluence_DA_semiinf_CW(1.0, 0.1, 10.0, n_med = 1.0, n_ext = 1.0, z = 0.0) ≈ 0.024035998288332954
@test fluence_DA_semiinf_CW(4.0, 0.01, 50.0, n_med = 1.0, n_ext = 1.0, z = 0.0) ≈ 7.2878997230756655e-6
@test fluence_DA_semiinf_CW(3.0, 0.1, 10.0, n_med = 1.0, n_ext = 1.0, z = 100.0) ≈ 6.938564742712919e-78


# TD
d = [0.000288659079311426, 2.896685181887841e-6, 5.474977667638328e-8, 1.3595622802163355e-9]
@test fluence_DA_semiinf_TD(1.0:4.0, 1.0, 0.1, 10.0, n_med = 1.0, n_ext = 1.0, z = 0.0) ≈ d

end # module