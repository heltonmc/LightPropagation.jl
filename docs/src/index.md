# LightPropagation.jl

## Introduction

LighPropagation provides a set of tools for modeling light transport in turbid media written in the
[Julia programming language](https://julialang.org/).
The main motivation behind the development of this library is to provide an easy-to-use, open-source framework that highlights and standardizes prominent analytical techniques to model light transport in complex heterogenous structures. The radiative transport equation (RTE) and its approximations to model propagation of particles in random media are frequently applied in astrophysics, nuclear physics, biophotonics, heat transfer, computer graphics, and climate research. This package focuses on their use in biomedical optics with a focus on performance and standardization.

The library has two main components split into (a) the forward modelling of light transport and (b) the associated inverse problem. Currently, photon migration in the continuous-wave (CW), frequency-domain (FD)and time-domain (TD) under the Diffusion Approximation (DA) are described in homogenous media for semi-infinite, slab, and parallelepiped geometries. Laterally infinite multi-layered solutions are also described. Least squares fitting based on the Levenberg-Marquardt algorithm are also available to fit time-resolved measurements after convolution with the Instrument Response Function (IRF). 

## Julia educational resources

A basic knowledge of the Julia programming language is needed to use the LightPropagation package.
Here, one can find a list of resources to get started with this programming language.

* Official webpage [docs.julialang.org](https://docs.julialang.org/)
* Official list of learning resources [julialang.org/learning](https://julialang.org/learning/)
* Official [YouTubechannel] (https://www.youtube.com/c/TheJuliaLanguage)

## Manual

```@contents
Pages = [
  "DA_slab_semiinfgeom.md",
  ]
```

