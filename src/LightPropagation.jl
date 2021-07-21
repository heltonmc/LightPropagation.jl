
module LightPropagation

using DelimitedFiles
using LsqFit
using DSP
using FFTW
using Parameters
using SpecialFunctions
using JLD

export TPSF_DA_semiinf_refl
export TPSF_DA_slab_refl
export TPSF_DA_paralpip_refl

export fluence_DA_inf_CW
export fluence_DA_inf_TD
export fluence_DA_semiinf_TD
export fluence_DA_semiinf_CW
export fluence_DA_slab_CW
export fluence_DA_slab_TD
export fluence_DA_paralpip_TD
export fluence_DA_paralpip_CW

export Nlayer_cylinder
export fluence_DA_Nlay_cylinder_CW
export fluence_DA_Nlay_cylinder_TD
export besselroots

export getfit

export load_asc_data

export fitDTOF, DTOF_fitparams

include("data/loadData/loaddata.jl")

include("forwardmodels/Diffusion Approximation/diffusionparameters.jl")

include("forwardmodels/Diffusion Approximation/Infinite/DAinfinite.jl")
include("forwardmodels/Diffusion Approximation/Semi-infinite/DAsemiinf.jl")
include("forwardmodels/Diffusion Approximation/Slab/DAslab.jl")
include("forwardmodels/Diffusion Approximation/Parallelepiped/DAparalpip.jl")

include("inversefitting/FitStructures/fitstructures.jl")
include("inversefitting/FitFunctions/fitfunction.jl")
include("inversefitting/FitFunctions/fitmodels.jl")

include("forwardmodels/Diffusion Approximation/transforms.jl")
include("forwardmodels/Diffusion Approximation/DAcylinder_layered.jl")

const besselroots = load(joinpath(@__DIR__,"..", "src/forwardmodels/Diffusion Approximation/besselzeroroots.jld"))["besselroots"]
       
end
