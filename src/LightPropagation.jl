
module LightPropagation

using Parameters
import SpecialFunctions
import Bessels
using JLD
using ForwardDiff
using FastGaussQuadrature: gausslegendre

besselj0(x::T) where T <: Union{Float32, Float64} = Bessels.besselj0(x)
besselj1(x::T) where T <: Union{Float32, Float64} = Bessels.besselj1(x)

besselj0(x) = SpecialFunctions.besselj0(x)
besselj1(x) = SpecialFunctions.besselj1(x)

### Diffusion Theory

# Infinite
export fluence_DA_inf_CW
export fluence_DA_inf_TD

# Semi-infinite
export fluence_DA_semiinf_TD
export fluence_DA_semiinf_CW

export flux_DA_semiinf_CW
export flux_DA_semiinf_TD

# Slab
export fluence_DA_slab_CW
export fluence_DA_slab_TD

export flux_DA_slab_CW
export flux_DA_slab_TD

# Parallelepiped
export fluence_DA_paralpip_TD
export fluence_DA_paralpip_CW

# Layered cylinder
export fluence_DA_Nlay_cylinder_CW, fluence_DA_Nlay_cylinder_CW_approx
export fluence_DA_Nlay_cylinder_TD, fluence_DA_Nlay_cylinder_TD_approx

export flux_DA_Nlay_cylinder_CW, flux_DA_Nlay_cylinder_CW_approx
export flux_DA_Nlay_cylinder_TD, flux_DA_Nlay_cylinder_TD_approx

# g2 for DCS
export g1_DA_semiinf_CW, g1_DA_semiinf_TD
export g2_DA_semiinf_CW, g2_DA_semiinf_TD

export g1_DA_Nlay_cylinder_CW, g1_DA_Nlay_cylinder_TD
export g2_DA_Nlay_cylinder_CW, g2_DA_Nlay_cylinder_TD

# Structures
export Nlayer_cylinder
export DAsemiinf_DCS, Nlayer_cylinder_DCS

# Abstract types
export DiffusionParameters

# constants
export J0_ROOTS
export J1_J0ROOTS_2

#export getfit
#export load_asc_data
#export fitDTOF, DTOF_fitparams

include("forwardmodels/Diffusion Approximation/diffusionparameters.jl")

include("forwardmodels/Diffusion Approximation/Infinite/DAinfinite.jl")
include("forwardmodels/Diffusion Approximation/Semi-infinite/DAsemiinf.jl")
include("forwardmodels/Diffusion Approximation/Slab/DAslab.jl")
include("forwardmodels/Diffusion Approximation/Parallelepiped/DAparalpip.jl")

include("forwardmodels/Diffusion Approximation/transforms.jl")
include("forwardmodels/Diffusion Approximation/DAcylinder_layered.jl")

include("forwardmodels/Diffusion Approximation/DCS/semiinf.jl")
include("forwardmodels/Diffusion Approximation/DCS/layered_cylinder.jl")

const J0_ROOTS = load(joinpath(@__DIR__,"..", "utils/besselroots/J0_ROOTS.jld"))["besselroots"]

# computes (besselj1(besselroots[ind]))^2
const J1_J0ROOTS_2 = load(joinpath(@__DIR__,"..", "utils/besselroots/J1_J0ROOTS_2.jld"))["J1_J0ROOTS_2"]

file_big = joinpath(@__DIR__,"..", "utils/besselroots/besselzeroroots_big.jld")

if isfile(file_big)
    const J0_ROOTSbig = load(file_big)["big_besselroots"]
    const J1_J0ROOTS_2big = (besselj1.(J0_ROOTSbig)).^2
end

end
