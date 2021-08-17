# LightPropagation.jl

![Repo status](https://www.repostatus.org/badges/latest/active.svg?style=flat-square)  [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://heltonmc.github.io/LightPropagation.jl/stable/) [![Join the chat at https://gitter.im/LightPropagation/community?source=orgpage](https://badges.gitter.im/LightPropagation/community.svg)](https://gitter.im/LightPropagation/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) ![GitHub tag (latest SemVer pre-release)](https://img.shields.io/github/v/tag/heltonmc/LightPropagation.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square) ![MIT license](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)

## Overview

LightPropagation.jl provides a set of tools to model and analyze the propagation of light in turbid media written in the [Julia programming language](https://julialang.org/). The primary goals of this project are to create robust analytical solutions of light transport that are easy to use, fast, accurate, and well tested. Current development is focused on robust implementations of solutions to the Diffusion Equation. 

The library currently supports simulating the fluence and flux (using Fick's law) in the following geometries under the diffusion approximation: **(a) Infinite, (b) Semi-Infinite, (c) Slab, (d) Parallelepiped, (e) Cylinder, (f) N-Layered Cylinder** for the Steady-State (CW), Frequency-Domain (FD), and Time-Domains (TD). We also provide a fitting routine for time-domain measurements using the standard Levenberg-Marquardt algorithm with consideration of the Instrument Response Function (IRF). 

## Installation

Install Julia by [downloading](https://julialang.org/downloads/) the latest version from the offical site and following the [platform specific installations](https://julialang.org/downloads/platform/). 

Then, use Julia's  built-in package manager (accessed by pressing `]` in the Julia REPL command prompt) to add the package and build all the required dependencies.

```julia
julia>]
(@v1.6) pkg> add "https://github.com/heltonmc/LightPropagation.jl.git"

julia> using LightPropagation
```

LightPropagation.jl can be updated to the latest tagged release from the package manager by typing

```julia
(@v1.6) pkg> update LightPropagation
```
LightPropagation.jl is only tested on Julia v1.5 and later. The package is under rapid development and breaking changes to the user API do occur, so update with care. Please open an issue if something breaks or doesn't seem right!

## Quick Start

To see a list of all exported models check out the [API in the documentation](https://heltonmc.github.io/LightPropagation.jl/stable/API/). A list of all exported names can also be found using the `names` function in the Julia REPL:
```julia
julia> names(LightPropagation)
26-element Vector{Symbol}:
 ....
 ....
 :fluence_DA_Nlay_cylinder_CW
 :fluence_DA_Nlay_cylinder_TD
 :fluence_DA_inf_CW
 :fluence_DA_inf_TD
 :fluence_DA_paralpip_CW
 :fluence_DA_paralpip_TD
 :fluence_DA_semiinf_CW
 :fluence_DA_semiinf_TD
 :fluence_DA_slab_CW
 :fluence_DA_slab_TD
 :flux_DA_Nlay_cylinder_CW
 :flux_DA_Nlay_cylinder_TD
 :flux_DA_semiinf_CW
 :flux_DA_semiinf_TD
 :flux_DA_slab_CW
 :flux_DA_slab_TD
 ....
 ....
```

LightPropagation.jl follows a naming pattern like PhysicalQuantity_Approximation_Geometry_MeasurementDomain. For example, `fluence_DA_semiinf_CW` models the fluence under the diffusion approximation for the semi-infinite geometry in the steady-state (continuous wave) domain.

To see how to use a specific model we can type `?` in the Julia REPL followed by the name of the model:
```julia
julia>?
help?> fluence_DA_semiinf_TD
search: fluence_DA_semiinf_TD fluence_DA_semiinf_CW

  fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, z = 0.0)

  Compute the time-domain fluence in a semi-infinite medium (Eqn. 33 Contini).

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  t: the time vector (ns).

    •  ρ: the source detector separation (cm⁻¹)

    •  μa: absorption coefficient (cm⁻¹)

    •  μsp: reduced scattering coefficient (cm⁻¹)

    •  n_ext: the boundary's index of refraction (air or detector)

    •  n_med: the sample medium's index of refraction

    •  z: the z-depth in medium

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> fluence_DA_semiinf_TD(0.1:0.1:1.0, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, z = 0.0)
```
This will give the function header, a description of the model, description of arguments, and an example on how to execute it. More detailed information on usage and examples can be found in the [documentation](https://heltonmc.github.io/LightPropagation.jl/stable/). If any problem is observed in the documentation or running a model please file an issue!

### Contributing to LightPropagation.jl

One of the primary goals of LightPropagation.jl was to create robust implementations of analytical models of light transport in turbid media. Ideally, such models would be applicable to a wide variety of research problems and would be highly tested, accurate, and performant. 

If you are interested in contributing in any way we would love to have your help no matter the size of contribution. If a model is not giving expected results, not performing the way you would expect, not matching your implementation please let us know! Some good first issues would be increasing performance of existing functions, adding more tests, or creating better documentation.

Let us know by [opening an issue](https://github.com/heltonmc/LightPropagation.jl/issues/new) if you would like to work on a new feature or if you are new to open-source and want to find a cool little project or issue to work on that fits your interests! We are more than happy to help you along the way.

### Citing

We hope you fork and modify the existing code base in anyway for your application. If you make an enhancement to the code please either try to become a contributor or let us know where we could be better! If you use one of the models in your research please cite the original authors who derived the equations.

### Contact

The best place to discuss usage, errors, bugs, features/requests, and potential contributions are in the issues and discussion forums here so everyone can benefit and participate. For other questions, please contact [Michael Helton](mailto:heltonmc@umich.edu).
