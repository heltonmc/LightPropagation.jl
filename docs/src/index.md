# LightPropagation.jl

## Introduction

LighPropagation provides a set of tools for modeling light transport in turbid media written in the
[Julia programming language](https://julialang.org/).
The main motivation behind the development of this library is to provide an easy-to-use, open-source framework that highlights and standardizes prominent analytical techniques to model light transport in different geometries. The radiative transport equation (RTE) and its approximations to model propagation of particles in random media are frequently applied in astrophysics, nuclear physics, biophotonics, heat transfer, computer graphics, and climate research. This package focuses on their use in biomedical applications with a focus on performance and accuracy.

The library has two main components split into (a) the forward modelling of light transport and (b) the associated inverse problem. Currently, photon migration in the continuous-wave (CW), frequency-domain (FD) and time-domain (TD) under the Diffusion Approximation (DA) are described in homogenous media for infinite, semi-infinite, slab, and parallelepiped geometries and for heterogenous media in the case of layered cylinders. Least squares fitting based on the Levenberg-Marquardt algorithm are also available to fit time-resolved measurements after convolution with the Instrument Response Function (IRF). 

## Julia educational resources

If you are interested in learning Julia (definitely encouraged!) here are some resources. In general, only a basic knowledge of the language is needed to use LightPropagation.jl.

* Official [documentation](https://docs.julialang.org/)
* [List of learning resources](https://julialang.org/learning/)
* [YouTube] (https://www.youtube.com/c/TheJuliaLanguage)

## How to get help?

If you are interested in using LightPropagation.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions. Do [open an issue](https://github.com/heltonmc/LightPropagation.jl/issues/new) or [pull request](https://github.com/heltonmc/LightPropagation.jl/pulls/new) if you have questions, suggestions or solutions. [Start a new discussion](https://github.com/heltonmc/LightPropagation.jl/discussions/new) for a more informal discussion!




