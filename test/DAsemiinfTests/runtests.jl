module DAsemiinfTests

using Test
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD


# Hard code DA semiinf solution to compare other geometries

# CW
@test fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0) ≈ 0.02403599828833302
@test fluence_DA_semiinf_CW(4.0, [0.01, 50.0], 1.0, 1.0, 0.0) ≈ 7.2878997230756655e-6
@test fluence_DA_semiinf_CW(3.0, [0.1, 10.0], 1.0, 1.0, 100.0) ≈ 6.938564742712919e-78


# TD
d = [ 0.00028865954060618826, 2.896693225195505e-6, 5.4749991457247e-8, 1.3595691610567748e-9]
@test fluence_DA_semiinf_TD(1.0:4.0, [0.1, 10.0], 1.0, 1.0, 1.0, 0.0) ≈ d

end # module