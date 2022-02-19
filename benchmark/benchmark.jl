using BenchmarkTools
using LightPropagation

suite = BenchmarkGroup()

suite["Layered_fluence_CW"] = BenchmarkGroup()
suite["Layered_fluence_CW"]["2-layer"] = @benchmarkable fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) setup=(ρ = 1.0; μa = (0.1, 0.2); μsp = (10.0, 11.0); n_ext = 1.0; n_med = (1.4, 1.4); l = (1.0, 10.0); a = 10.0; z = 0.0; N_J0Roots = 1000)
suite["Layered_fluence_CW"]["3-layer"] = @benchmarkable fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) setup=(ρ = 1.0; μa = (0.1, 0.2, 0.15); μsp = (10.0, 11.0, 12.5); n_ext = 1.0; n_med = (1.4, 1.4, 1.4); l = (1.0, 2.0, 10.0); a = 10.0; z = 0.0; N_J0Roots = 1000)
suite["Layered_fluence_CW"]["4-layer"] = @benchmarkable fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) setup=(ρ = 1.0; μa = (0.1, 0.2, 0.15, 0.2); μsp = (10.0, 11.0, 12.5, 12.0); n_ext = 1.0; n_med = (1.4, 1.4, 1.4, 1.4); l = (1.0, 2.0, 2.0, 10.0); a = 10.0; z = 0.0; N_J0Roots = 1000)
suite["Layered_fluence_CW"]["10-layer"] = @benchmarkable fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, N_J0Roots) setup=(ρ = 1.0; μa = (0.1, 0.2, 0.15, 0.2, 0.1, 0.2, 0.15, 0.2, 0.1, 0.2); μsp = (10.0, 11.0, 12.5, 12.0, 10.0, 11.0, 12.5, 12.0, 10.0, 11.0); n_ext = 1.0; n_med = (1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4); l = (1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0, 2.0, 2.0); a = 10.0; z = 0.0; N_J0Roots = 1000)

suite["semiinf_fluence_CW"] = BenchmarkGroup()
suite["semiinf_fluence_CW"] = @benchmarkable fluence_DA_semiinf_CW(ρ, μa, μsp) setup=(ρ = 1.0*rand(); μa = 0.1*rand(); μsp = 10.0*rand())

tune!(suite)