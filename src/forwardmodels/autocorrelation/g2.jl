function g2_semiinf(τ, B, ρ; nmed = 1.33, ndet = 1.5, μsp = 10.0, μa = 0.1, λ = 700)



    αDB::Float64 = B[1]
    β::Float64 = B[2]
    n::Float64 = nmed/ndet

    D::Float64 = 1/3μsp



    if n > 1.0
        A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
        5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
    elseif n < 1.0
        A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
    else 
        A = 1.0
    end

    z0::Float64 = 1/(μsp + μa)
    zb::Float64 = 2A*D
    k0::Float64 = (n*2*pi)/(λ*10^-7)

    r1 = sqrt(ρ^2 +z0^2)
    r2 = sqrt(ρ^2 + (z0 +2zb)^2)

    C1 = 3*μsp*μa
    C2 = 6*μsp^2*k0^2*αDB

    A = zeros(length(τ))

    Threads.@threads for n in eachindex(τ)
        Kτ = sqrt(C1 + C2*τ[n])

        A[n] = exp(-r1*Kτ)/r1 - exp(-r2.*Kτ)/r2
    end

    g1 = (A)./(A[1])
    g2 = 1 .+ β*g1.^2
end


function test(τ, B::Array{Float64,1}, ρ::Float64; nmed::Float64 = 1.33, ndet::Float64 = 1.5, μsp::Float64 = 10.0, μa::Float64 = 0.1, λ::Float64 = 700.0)
 

    αDB::Float64 = B[1]
    β::Float64 = B[2]
    n::Float64 = nmed/ndet

    D::Float64 = 1/3μsp



    if n > 1.0
        A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
        5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
    elseif n < 1.0
        A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
    else 
        A = 1.0
    end

    z0::Float64 = 1/(μsp + μa)
    zb::Float64 = 2A*D
    k0::Float64 = (n*2*pi)/(λ*10^-7)

    r1 = sqrt(ρ^2 +z0^2)
    r2 = sqrt(ρ^2 + (z0 +2zb)^2)

    C1 = 3*μsp*μa
    C2 = 6*μsp^2*k0^2*αDB

    A = zeros(length(τ))

    Threads.@threads for n in eachindex(τ)
        Kτ = sqrt(C1 + C2*τ[n])

        A[n] = exp(-r1*Kτ)/r1 - exp(-r2.*Kτ)/r2
    end

    g1 = (A)./(A[1])
    g2 = 1 .+ β*g1.^2
end

beta = 0.5
BFi = 2e-8
τ = 1e-8:1e-6:1e-2
test(τ, [BFi, beta], 1.)



t = 0:0.01:10

ydata =  test(τ, [BFi, beta], 1.) + 0.05*randn(length(τ))


function curvetestff(ydata, τ)

    fit = curve_fit((τ, β) -> test(τ, β, 1.0), τ, ydata, [5e-9, 0.3]; autodiff=:finiteforward)
end


#=
@benchmark test(τ, [BFi, beta], 1.)
BenchmarkTools.Trial: 
  memory estimate:  704.81 KiB
  allocs estimate:  19576
  --------------
  minimum time:     96.266 μs (0.00% GC)
  median time:      101.922 μs (0.00% GC)
  mean time:        146.992 μs (23.95% GC)
  maximum time:     8.491 ms (97.80% GC)
  --------------
  samples:          10000
  evals/sample:     1
