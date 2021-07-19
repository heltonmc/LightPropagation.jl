module DAinfTests

using Test
using LightPropagation: fluence_DA_inf_CW, fluence_DA_inf_TD

# Hard code DA inf solution to compare other geometries

# CW
@test fluence_DA_inf_CW(0.1, 0.1, 10.0) ≈ 20.076563644369305
@test fluence_DA_inf_CW(0.1, 0.2, 15.0) ≈ 26.52859839465854
@test fluence_DA_inf_CW(60.0, 0.2, 15.0) ≈ 4.007233568620553e-80

# TD
d = [5.87809115624927e-14, 0.0019712791801801076, 2.10256922566133e-6, 3.038559267148949e-9]
@test  fluence_DA_inf_TD(0.01:1.0:4.0, 1.0, 0.2, 15.0) ≈ d

end # module