
module LightPropagation

using Parameters
using SpecialFunctions
using JLD
using ForwardDiff
using FastGaussQuadrature: gausslegendre

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
export fluence_DA_Nlay_cylinder_CW
export fluence_DA_Nlay_cylinder_TD

export flux_DA_Nlay_cylinder_CW
export flux_DA_Nlay_cylinder_TD

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
    
end
