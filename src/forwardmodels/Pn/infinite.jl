function compute_roots(μa, μa1, μa2, μa3)
    b = 35 * μa2 * μa3 + 28 * μa * μa3 + 27 * μa * μa1
    c = 105 * μa * μa1 * μa2 * μa3

    γ1 = (-b + sqrt(b^2 - 36 * c)) / 18
    γ2 = (-b - sqrt(b^2 - 36 * c)) / 18
    return γ1, γ2
end

char_poly(γ, μa1, μa2, μa3) = (105 * μa1 * μa2 * μa3) + (28 * μa3 * γ) + (27 * μa1 * γ)

q_poly(γ, μa, μa1, μa2, μa3) = (105 * μa * μa1 * μa2 * μa3) + (35 * μa2 * μa3 * γ) + (28 * μa * μa3 * γ) + (27 * μa * μa1 * γ) + (9 * γ^2)

function green(ρ, μa, μs, g)

    μa1 = μa + (1 - g) * μs #/ (1 - g)
    μa2 = μa + (1 - g^2) * μs #/ (1 - g)
    μa3 = μa + (1 - g^3) * μs #/ (1 - g)


    γ1, γ2 = compute_roots(μa, μa1, μa2, μa3)
    #@show q_poly(γ1, μa, μa1, μa2, μa3)
    #@show q_poly(γ2, μa, μa1, μa2, μa3)

    v1 = sqrt(-γ1)
    v2 = sqrt(-γ2)
    @show γ1, γ2
    #@show v1, v2

    P1 = char_poly(γ1, μa1, μa2, μa3)
    P2 = char_poly(γ2, μa1, μa2, μa3)

    A1 = P1 / 9 # 3550.392247461887
    #A1 = 3 * μa1

    A2 = P2 / (9 * (γ2 - γ1))
   
    @show v1
    ϕ = A1 * exp(-v1 * ρ) + A2 * exp(-v2 * ρ)

    return ϕ / (4 * π * ρ)
end

function green1(ρ, μa, μs, g)
    μa1 = μa + (1 - g) * μs

    @show sqrt(3*μa1*μa)
    ϕ = exp(-sqrt(3*μa1*μa) * ρ) * 3 *μa1
    return ϕ / (4 * π * ρ)
end



exp(-sqrt(3 * μsp * μa) * ρ) / (4 * π * ρ * D)
phi_0 = phi_0 + 1./(2*pi*mua(1)*r).*exp(-r./D(k,k))./D(k,k).^2.*V(1,k).^2;