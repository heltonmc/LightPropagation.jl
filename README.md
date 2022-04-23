<h1 align="center">
LightPropagation.jl
</h1>
<p align="center">Fast and accurate numerical routines for simulating light transport in scattering media.
  
 <div align="center">

![Repo status](https://www.repostatus.org/badges/latest/active.svg?style=flat-square)  [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://heltonmc.github.io/LightPropagation.jl/stable/)  ![GitHub tag (latest SemVer pre-release)](https://img.shields.io/github/v/tag/heltonmc/LightPropagation.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square)  ![MIT license](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)  [![Slack](https://img.shields.io/badge/chat-slack-e01e5a)](https://join.slack.com/t/lightpropagationjl/shared_invite/zt-z4r47k5w-cco_rRK5EmKQqq3pdqykdw)
 
</div>  
  
## Overview

Numerical implementations of analytical solutions to the radiative transport equation focused on modeling light propagation through turbid media like biological tissue.

Primary goals: provide fast, well-tested numerical recipes for use in inverse problems in biomedical optics research (though solutions can be used for a variety of applications in turbid media)

Current development is focused on solutions to the diffusion equation (and correlation equation). The library supports simulating the fluence and flux (using Fick's law) in the following geometries under the diffusion approximation: **(a) Infinite, (b) Semi-Infinite, (c) Slab, (d) Parallelepiped, (e) N-Layered Cylinder** for the Steady-State (CW), Frequency-Domain (FD), and Time-Domains (TD). 

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

**Operating Systems:** LightPropagation.jl is tested on the [latest versions](https://github.com/actions/virtual-environments) of MacOS, Windows, and Ubuntu using v1.6 and the latest stable release of Julia. Please see Julia's [supported platforms](https://julialang.org/downloads/#supported_platforms) form more information. On Julia versions >1.8.0, LightPropagation.jl will pass all tests on macOS ARM (M-series Processor).

## Quick Start

To see a list of all exported models check out the [API in the documentation](https://heltonmc.github.io/LightPropagation.jl/stable/API/). 
LightPropagation.jl follows a naming pattern like PhysicalQuantity_Approximation_Geometry_Domain. For example, `fluence_DA_semiinf_CW` models the fluence under the diffusion approximation for the semi-infinite geometry in the steady-state (continuous wave) domain.

```julia
julia> ρ = 1.0; μa = 0.1; μsp = 10.0
julia> fluence_DA_semiinf_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
0.024035998288332864

# s is the slab thickness
julia> fluence_DA_slab_CW(ρ, μa, μsp; n_ext = 1.0, n_med = 1.0, s = 1.0)
0.021728999570708202
```
To see how to use a specific model we can type `?` in the Julia REPL followed by the name of the model (`?fluence_DA_slab_CW`) and press enter.

For simulations in the time-domain we can initially define a vector (or scalar) of time values...
```julia
julia> t = 0.1:0.4:2.1

julia> fluence_DA_semiinf_TD(t, ρ, μa, μsp; n_ext = 1.0, n_med = 1.0)
6-element Vector{Float64}:
 0.13275269870908793
 0.005646876427088415
 0.0004926251803161931
 6.468625297894211e-5
 1.0448483663828698e-5
 1.9116131901703835e-6
 ```
 
 ### Layered media
 In layered media we can simply give the optical coefficients for each layer...
 ```julia
 julia> μa = (0.1, 0.2); μsp = (10.0, 12.0)
(10.0, 12.0)

julia> fluence_DA_Nlay_cylinder_CW(ρ, μa, μsp, n_med=(1.0, 1.0))
0.023755152650609134

julia> fluence_DA_Nlay_cylinder_TD(t, ρ, μa, μsp, n_med=(1.0, 1.0))
6-element Vector{Float64}:
 0.13275618056673205
 0.005134284906594011
 0.00028378694581441963
 1.9158665807291887e-5
 1.4243249932846215e-6
 1.1205939970902316e-7
```


More detailed information on usage and examples can be found in the [documentation](https://heltonmc.github.io/LightPropagation.jl/stable/). If any problem is observed in the documentation or running a model please file an issue!

### How to support and contribute

You can support the project by actively using it and [raising issues](https://github.com/heltonmc/LightPropagation.jl/issues/new). User feedback on API is especially appreciated as we hope to get to a more stable (v1.0) release. Leaving a star if you like the project is also appreciated and motivating!

Feature requests are especially important as we want this to be in anway helpful in your research or personal enjoyment. For example, if you need more example scripts for loading that weird TCSPC data... or you want us to support a new model... let us know! Opening an issue or a [discussion](https://github.com/heltonmc/LightPropagation.jl/discussions/new) will be the best place to do this.

Contributions via pull request are also greatly appreciated and we would love to have your help no matter the size of contribution. If new to Julia or Git, we are happy to help you through the process. Let us know by [opening an issue](https://github.com/heltonmc/LightPropagation.jl/issues/new) if you would like to work on a new feature or if you are new to open-source and want to find a cool little project or issue to work on that fits your interests! We are more than happy to help you along the way. If a model is not giving expected results, not performing the way you would expect, not matching your implementation please let us know! Some good first issues would be increasing performance of existing functions, adding more tests, or creating better documentation.

### Citing

We hope you fork and modify the existing code base in anyway for your application. If you make an enhancement to the code please either try to become a contributor or let us know where we could be better! If you use one of the models in your research please cite the original authors who derived the equations. They are listed at the top of the source file for each function.

### Contact

The best place to discuss usage, errors, bugs, features/requests, and potential contributions are in the issues and discussion forums here so everyone can benefit and participate. For more casual discussion, [join us on Slack](https://join.slack.com/t/lightpropagationjl/shared_invite/zt-z4r47k5w-cco_rRK5EmKQqq3pdqykdw). For other questions, please contact [Michael Helton](mailto:heltonmc@umich.edu).
