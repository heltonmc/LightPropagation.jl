module DA_NlayerkernelTests

using Test
using Parameters, SpecialFunctions

abstract type DiffusionParameters end

include(joinpath(dirname(@__FILE__), "..", "..","src/forwardmodels/Diffusion Approximation/DAcylinder_layered.jl"))

## Test the derivations with original

@inline function α_coeff!(α, μa, D, sn)
    @inbounds for ind in 1:length(μa)
        α[ind] = sqrt(μa[ind] / D[ind] + sn^2)
    end
    return α
end

# greens function for 2 layer cylinder as presented by Liemert
function green_2layer(α, sn, μa, D, z, z0, zb, l, n)
    α = α_coeff!(α, μa, D, sn)


    β3 = sinh(α[2]*(l[2] + zb[2]))
    γ3 = cosh(α[2]*(l[2] + zb[2]))

    g1 = sinh(α[1]*(z0 + zb[1]))*sinh(α[1]*(z + zb[1])) / (D[1]*α[1]*exp(α[1]*(l[1] + zb[1])))
    g1 *= (D[1]*α[1]*n[1]^2*β3 - D[2]*α[2]*n[2]^2*γ3) / (D[1]*α[1]*n[1]^2*β3*cosh(α[1]*(l[1] + zb[1])) + (D[2]*α[2]*n[2]^2*γ3*sinh(α[1]*(l[1] + zb[1]))))

    g = (exp(-α[1]*abs(z - z0)) - exp(-α[1]*(z + z0 + 2zb[1]))) / (2*D[1]*α[1])

    return g + g1
end
# greens function for 3 layer cylinder as presented by Liemert
function green_3layer(α, sn, μa, D, z, z0, zb, l, n)
    α = α_coeff!(α, μa, D, sn)

    β3 = D[2]*α[2]*n[2]^2*cosh(α[2]*l[2])*sinh(α[3]*(l[3] + zb[2])) + D[3]*α[3]*n[3]^2*sinh(α[2]*l[2])*cosh(α[3]*(l[3] + zb[2]))
    γ3 = D[2]*α[2]*n[2]^2*sinh(α[2]*l[2])*sinh(α[3]*(l[3] + zb[2])) + D[3]*α[3]*n[3]^2*cosh(α[2]*l[2])*cosh(α[3]*(l[3] + zb[2]))

    g1 = sinh(α[1]*(z0 + zb[1]))*sinh(α[1]*(z + zb[1])) / (D[1]*α[1]*exp(α[1]*(l[1] + zb[1])))
    g1 *= (D[1]*α[1]*n[1]^2*β3 - D[2]*α[2]*n[2]^2*γ3) / (D[1]*α[1]*n[1]^2*β3*cosh(α[1]*(l[1] + zb[1])) + (D[2]*α[2]*n[2]^2*γ3*sinh(α[1]*(l[1] + zb[1]))))

    g = (exp(-α[1]*abs(z - z0)) - exp(-α[1]*(z + z0 + 2zb[1]))) / (2*D[1]*α[1])

    return g + g1
end

### test the 2 layer case
α = [1.2, 1.5]
sn = 2.2
D = [0.2, 0.5]
μa = [0.1, 0.14]
z0 = 0.2
zb = [0.2, 0.3]
l = [1.0, 2.0]
n = [1.0, 1.0]
n2 = @. D * n^2
z = 0.0

@test _green_Nlaycylin_top(sn, μa, D, z, z0, zb, l, n2, 2) ≈ green_2layer(α, sn, μa, D, z, z0, zb, l, n)

### test the 3 layer case
α = [1.2, 1.5, 1.3]
sn = 5.8
D = [0.2, 0.5, 0.7]
μa = [0.1, 0.14, 0.2]
z0 = 0.2
zb = [0.2, 0.3, 0.3]
l = [1.0, 2.0, 1.2]
n = [1.0, 1.0, 1.1]
n2 = @. D * n^2
z = 0.0

@test _green_Nlaycylin_top(sn, μa, D, z, z0, zb, l, n2, 3) ≈ green_3layer(α, sn, μa, D, z, z0, zb, l, n)

end # module
