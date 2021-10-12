"""
    abstract type DiffusionParameters

Abstract type for all input structures to diffusion functions.
"""
abstract type DiffusionParameters end

"""
    A_coeff(n)

Function to approximate the A coefficient (mismatch between external and internal index of refractions) given in Eq. 25 from Ref. [1]. 
n is the ratio between the refractive index of the diffusing medium and that of the external medium (n_med / n_ext).
Uses the polynomial fit given in Appendix A [1]. Used to calculate the extrapolation length.

References:
[1] Daniele Contini, Fabrizio Martelli, and Giovanni Zaccanti, 
    "Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory," 
    Appl. Opt. 36, 4587-4599 (1997) 
"""
function A_coeff(n)
    if n > 1.0
        A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
        5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
    elseif n < 1.0
        A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
    else 
        A = 1.0
    end
    return A
end

"""
    D_coeff(μsp)

Returns the diffusion coefficient as a function of the reduced scattering coefficient.
Note that this is sometimes defined with the absorption coefficient D = 1 / (3 * (μsp + μa)).
"""
D_coeff(μsp) = 1 / (3 * μsp)

"""
    ν_coeff(n_med)

Computes the speed of light in the medium (cm/ns). 
This is a good place to change the units of the returning function along with optical properties if interested in mm or ps.
"""
ν_coeff(n_med) = 29.9792458 / n_med

"""
    z0_coeff(μsp)

The location of the isotropic point source within the medium assuming some perpindicular incidence of light beam.
Note that this is sometimes defined with the absorption coefficient z0 = 1 / (μsp + μa)
"""
z0_coeff(μsp) = 1 / μsp

"""
    zb_coeff(A, D)

The extrapolation length calculated with the A coefficient (accounts for mismatch index of refraction at boundary)
and D (diffusion coefficient).
"""
zb_coeff(A, D) = 2 * A * D
