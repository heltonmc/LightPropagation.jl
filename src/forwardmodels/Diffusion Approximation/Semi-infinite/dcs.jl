function g2_DA_semiinf_CW(τ::AbstractVector, ρ, μa, μsp; BFi = 2e-8, β = 1.0, n_ext = 1.0, n_med = 1.0, z = 0.0, λ = 700)
    g2 = similar(τ)

    k0 = 2 * π / (λ * 1e-7) # this should be in cm
    tmp = 2 * μsp * BFi * k0^2

    G0 = fluence_DA_semiinf_CW(ρ, μa, μsp, n_ext = n_ext, n_med = n_med, z = z)

    for ind in eachindex(τ)
        μa_dynamic = μa + tmp * τ[ind]
        G1 = fluence_DA_semiinf_CW(ρ, μa_dynamic, μsp, n_ext = n_ext, n_med = n_med, z = z)
        g1 = G1 / G0
        g2[ind] = 1 + β * abs(g1)^2
    end

    return g2
end

function g2_DA_Nlay_cylinder_CW(τ::AbstractVector, ρ, μa, μsp; BFi = [2e-8, 2e-8], β = 1.0, n_ext = 1.0, n_med = 1.0, l = [1.0, 10.0], a = 20.0, z = 0.0, λ =700.0, bessels = besselroots[1:1000])
    g2 = similar(τ)
    μa_dynamic = similar(μa)

    k0 = 2 * π / (λ * 1e-7) # this should be in cm
    tmp = @. 2 * μsp * BFi * k0^2

    G0 = fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_ext, n_med, l, a, z, bessels)

    Threads.@threads for ind in eachindex(τ)

        μa_dynamic = μa + tmp * τ[ind]

        G1 = fluence_DA_Nlay_cylinder_CW(ρ, μa_dynamic, μsp, n_ext, n_med, l, a, z, bessels)
    
        g1 = G1 / G0
        g2[ind] = 1 + β * abs(g1)^2
    end

    return g2
end

τ = 10 .^(range(-10,stop=0,length=250))
BFi = 2e-8

si = g2_DA_semiinf_CW(τ,1.0,0.1, 10.0, BFi = BFi)

layered  = g2_DA_Nlay_cylinder_CW(τ, 1.0, [0.1, 0.1], [10.0, 10.0]; n_ext = 1.0, n_med = [1.0, 1.0], l = [0.5, 10.0], a = 20.0, z = 0.0, bessels = besselroots[1:1000])
layered2  = g2_DA_Nlay_cylinder_CW(τ, 1.0, [0.1, 0.1], [10.0, 10.0]; BFi = [2e-8, 5e-8], n_ext = 1.0, n_med = [1.0, 1.0], l = [0.5, 10.0], a = 20.0, z = 0.0, bessels = besselroots[1:1000])
layered3  = g2_DA_Nlay_cylinder_CW(τ, 1.0, [0.1, 0.1], [10.0, 10.0]; BFi = [5e-8, 2e-8], n_ext = 1.0, n_med = [1.0, 1.0], l = [0.5, 10.0], a = 20.0, z = 0.0, bessels = besselroots[1:1000])

scatter(log10.(τ), si, label = "Semi-inf - BFi = 2e-8 cm²/s", ylabel = "g2(τ)", xlabel = "log(τ (s))")
plot!(log10.(τ), layered, label = "Layered - BFi = [2e-8, 2e-8] cm²/s", lw = 1.5, linecolor =:black)
plot!(log10.(τ), layered2, label = "Layered - BFi = [2e-8, 5e-8] cm²/s", linestyle=:dash, lw = 1.5, linecolor =:red)
plot!(log10.(τ), layered3, label = "Layered - BFi = [5e-8, 2e-8] cm²/s", linestyle=:dot, lw = 1.5, linecolor =:blue4)

