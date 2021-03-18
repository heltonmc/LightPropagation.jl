module DA_NlayerkernelTests

using Test
using Parameters, SpecialFunctions

include(joinpath(dirname(@__FILE__), "..", "..","src/forwardmodels/Diffusion Approximation/DAcylinder_layered.jl"))

## Test the derivations with original

# greens function for 2 layer cylinder as presented by Liemert
function green_2layer(α, D, z0, zb, l, n, z)

    β3 = sinh(α[2]*(l[2] + zb[2]))
    γ3 = cosh(α[2]*(l[2] + zb[2]))

    g1 = sinh(α[1]*(z0 + zb[1]))*sinh(α[1]*(z + zb[1])) / (D[1]*α[1]*exp(α[1]*(l[1] + zb[1])))
    g1 *= (D[1]*α[1]*n[1]^2*β3 - D[2]*α[2]*n[2]^2*γ3) / (D[1]*α[1]*n[1]^2*β3*cosh(α[1]*(l[1] + zb[1])) + (D[2]*α[2]*n[2]^2*γ3*sinh(α[1]*(l[1] + zb[1]))))

    g = (exp(-α[1]*abs(z - z0)) - exp(-α[1]*(z + z0 + 2zb[1]))) / (2*D[1]*α[1])

    return g + g1
end
# greens function for 3 layer cylinder as presented by Liemert
function green_3layer(α, D, z0, zb, l, n, z)

    β3 = D[2]*α[2]*n[2]^2*cosh(α[2]*l[2])*sinh(α[3]*(l[3] + zb[2])) + D[3]*α[3]*n[3]^2*sinh(α[2]*l[2])*cosh(α[3]*(l[3] + zb[2]))
    γ3 = D[2]*α[2]*n[2]^2*sinh(α[2]*l[2])*sinh(α[3]*(l[3] + zb[2])) + D[3]*α[3]*n[3]^2*cosh(α[2]*l[2])*cosh(α[3]*(l[3] + zb[2]))

    g1 = sinh(α[1]*(z0 + zb[1]))*sinh(α[1]*(z + zb[1])) / (D[1]*α[1]*exp(α[1]*(l[1] + zb[1])))
    g1 *= (D[1]*α[1]*n[1]^2*β3 - D[2]*α[2]*n[2]^2*γ3) / (D[1]*α[1]*n[1]^2*β3*cosh(α[1]*(l[1] + zb[1])) + (D[2]*α[2]*n[2]^2*γ3*sinh(α[1]*(l[1] + zb[1]))))

    g = (exp(-α[1]*abs(z - z0)) - exp(-α[1]*(z + z0 + 2zb[1]))) / (2*D[1]*α[1])

    return g + g1
end

### test the 2 layer case
α = [1.2, 1.5]
D = [0.2, 0.5]
z0 = 0.2
zb = [0.2, 0.3]
l = [1.0, 2.0]
n = [1.0, 1.0]
z = 0.0

@test _green_Nlaycylin(α, D, z0, zb, l, n) ≈ green_2layer(α, D, z0, zb, l, n, z)

### test the 3 layer case
α = [1.2, 1.5, 1.3]
D = [0.2, 0.5, 0.7]
z0 = 0.2
zb = [0.2, 0.3, 0.3]
l = [1.0, 2.0, 1.2]
n = [1.0, 1.0, 1.1]
z = 0.0

@test _green_Nlaycylin(α, D, z0, zb, l, n) ≈ green_3layer(α, D, z0, zb, l, n, z)


end # module