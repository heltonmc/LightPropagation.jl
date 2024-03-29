
module LightPropagation

using DelimitedFiles
using LsqFit
using DSP
using FFTW
using Parameters
using SpecialFunctions
using JLD
using ForwardDiff

export TPSF_DA_semiinf_refl
export TPSF_DA_slab_refl
export TPSF_DA_paralpip_refl

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

# Structures
export Nlayer_cylinder

# constants
export besselroots

export getfit
export load_asc_data
export fitDTOF, DTOF_fitparams

include("forwardmodels/Diffusion Approximation/diffusionparameters.jl")

include("forwardmodels/Diffusion Approximation/Infinite/DAinfinite.jl")
include("forwardmodels/Diffusion Approximation/Semi-infinite/DAsemiinf.jl")
include("forwardmodels/Diffusion Approximation/Slab/DAslab.jl")
include("forwardmodels/Diffusion Approximation/Parallelepiped/DAparalpip.jl")

#include("inversefitting/FitStructures/fitstructures.jl")
#include("inversefitting/FitFunctions/fitfunction.jl")
#include("inversefitting/FitFunctions/fitmodels.jl")

include("forwardmodels/Diffusion Approximation/transforms.jl")
include("forwardmodels/Diffusion Approximation/DAcylinder_layered.jl")

const besselroots = load(joinpath(@__DIR__,"..", "src/forwardmodels/Diffusion Approximation/besselzeroroots.jld"))["besselroots"]
       
end
