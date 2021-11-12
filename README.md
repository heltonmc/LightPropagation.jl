# LightPropagation.jl

![Repo status](https://www.repostatus.org/badges/latest/active.svg?style=flat-square)  [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://heltonmc.github.io/LightPropagation.jl/stable/)![GitHub tag (latest SemVer pre-release)](https://img.shields.io/github/v/tag/heltonmc/LightPropagation.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square) ![MIT license](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)

## Overview

Numerical implementations of analytical solutions to the radiative transport equation focused on modeling light propagation through turbid media like biological tissue.

Primary goals: provide fast, well-tested numerical recipes for use in inverse problems in biomedical optics research (though solutions can be used for a variety of applications in turbid media)

Current development is focused on solutions to the diffusion equation (and diffuse correlation equation). The library supports simulating the fluence and flux (using Fick's law) in the following geometries under the diffusion approximation: **(a) Infinite, (b) Semi-Infinite, (c) Slab, (d) Parallelepiped, (e) Cylinder, (f) N-Layered Cylinder** for the Steady-State (CW), Frequency-Domain (FD), and Time-Domains (TD). 

We also provide fast interfaces for inverse fitting of time-domain measurements with consideration of the Instrument Response Function.

## Installation

Install Julia by [downloading](https://julialang.org/downloads/) the latest version from the offical site. This package requires Julia versions >=1.5.

**LightPropagation.jl** is a [registered package](https://juliahub.com/ui/Packages/LightPropagation/Wheva/) in the Julia package manager.
Just add the package to your environment using Julia's built-in package manager (accessed by pressing `]` in the Julia REPL command prompt) to add the package and build all the required dependencies.

```julia
julia> ] add LightPropagation

julia> using LightPropagation
```

This will add the latest tagged release to your environment. LightPropagation.jl can be updated to the latest tagged release from the package manager with

```julia
julia> ] update LightPropagation
```

You can also add the package directly from GitHub to get the latest changes between releases:
```
julia> ] add "https://github.com/heltonmc/LightPropagation.jl.git"
```

LightPropagation.jl is only tested on Julia v1.5 and later. The package is under rapid development and breaking changes to the user API do occur, so update with care. Please open an issue if something breaks or doesn't seem right!

## Quick Start

To see a list of all exported models check out the [API in the documentation](https://heltonmc.github.io/LightPropagation.jl/stable/API/). 
LightPropagation.jl follows a naming pattern like PhysicalQuantity_Approximation_Geometry_Domain. For example, `fluence_DA_semiinf_CW` models the fluence under the diffusion approximation for the semi-infinite geometry in the steady-state (continuous wave) domain.

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

This will give the function header, a description of the model, description of arguments, and an example on how to execute it. 
More detailed information on usage and examples can be found in the [documentation](https://heltonmc.github.io/LightPropagation.jl/stable/). If any problem is observed in the documentation or running a model please file an issue!

### How to support and contribute

You can support the project by actively using it and [raising issues](https://github.com/heltonmc/LightPropagation.jl/issues/new). User feedback on API is especially appreciated as we hope to get to a more stable (v1.0) release. If you like the project you can also leave a star. 

Feature requests are especially important as we want this to be in anway helpful in your research or personal enjoyment. For example, if you need more example scripts for loading that weird TCSPC data... or you want us to support a new model... let us know! Opening an issue or a [discussion](https://github.com/heltonmc/LightPropagation.jl/discussions/new) will be the best place to do this.

Contributions via pull request are also greatly appreciated and we would love to have your help no matter the size of contribution. If new to Julia or Git, we are happy to help you through the process. Let us know by [opening an issue](https://github.com/heltonmc/LightPropagation.jl/issues/new) if you would like to work on a new feature or if you are new to open-source and want to find a cool little project or issue to work on that fits your interests! We are more than happy to help you along the way. If a model is not giving expected results, not performing the way you would expect, not matching your implementation please let us know! Some good first issues would be increasing performance of existing functions, adding more tests, or creating better documentation.

### Citing

We hope you fork and modify the existing code base in anyway for your application. If you make an enhancement to the code please either try to become a contributor or let us know where we could be better! If you use one of the models in your research please cite the original authors who derived the equations.

### Contact

The best place to discuss usage, errors, bugs, features/requests, and potential contributions are in the issues and discussion forums here so everyone can benefit and participate. For other questions, please contact [Michael Helton](mailto:heltonmc@umich.edu).
