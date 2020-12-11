# LightPropagation

## Purpose

LightPropagation provides a set of tools to model and analyze the propagation of light in turbid media written in the [Julia programming language](https://julialang.org/).
The current library supports simulating the time-resolved reflectance in the forward direction under the diffusion approximation for 3 different geometries:
- __Semi-infinite__ 
- __Slabs__
- __Parallelepipeds__

## Usage

The easiest way to install julia is by [downloading julia](https://julialang.org/downloads/) from the offical site and following the [platform specific installations](https://julialang.org/downloads/platform/). 

Launch Julia and open Julia's package manager by typing `]` in the REPL. You should see the command line change from `julia>` to `(@v1.5) pkg>`. (The @v1.x will display your current version) Once the package has cloned, backspace to bring back the `julia>` in your REPL and type `using LightPropagation` as shown below.

```julia
(@v1.5) pkg> add "https://github.com/heltonmc/LightPropagation.git"

julia> using LightPropagation

```
## Tutorials

You should now have access to all of the exported functions of the package.

#### Forward Simulation

To view a functions inputs and methods type `?` in the REPL and then the name of the function. Let's first look at simulating a temporal point spread function (TPSF) for a semi-infinite medium.

```julia
help?> TPSF_DA_semiinf_refl
search: TPSF_DA_semiinf_refl

  TPSF_DA_semiinf_refl(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)

  Compute the time-domain reflectance from a semi-infinite medium. 

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •    t: the time vector (ns). 

    •    β::Array{Float64,1}: the optical properties μa, μs' (cm⁻¹)

    •    ρ::Float64: the source detector separation (cm⁻¹)

    •    ndet::Float64: the boundary's index of refraction (air or
        detector)

    •    nmed::Float64: the sample medium's index of refraction

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> TPSF_DA_semiinf_refl(0:1:5, [0.1,10.0], 1.0, 1.0, 1.0)
  6-element Array{Float64,1}:
   0.0
   0.0001440103022493725
   1.446739954231315e-6
   2.7354735244571076e-8
   6.794070985483474e-10
   1.9657536202689858e-11
```
This is a convenient way to see a brief description of the function and its arguments. Let's call this function and set it equal to a variable called TPSF. We can also plot this with the below example. 

```julia

julia> using Plots

julia> t = 0:0.01:10
julia> TPSF = TPSF_DA_semiinf_refl(t, [0.1,10.0], 1.0, 1.0, 1.0)
julia> ind = TPSF .> 0

julia> plot(t[ind], TPSF[ind], yscale=:log10, lw = 2, ylabel="Counts", xlabel="time (ns)", label="TPSF")

```
There are several other forward models in different geometries. The naming scheme of the forward models follow a physicalquantity_approximation_geometry_measurementtype scheme. So something that simulates the reflected temporal point spread function using the diffusion approximation for a semi-infinite geometry would be `TPSF_DA_semiinf_refl`. Try simulating the TPSF for a slab and parralelepiped geometry by using `TPSF_DA_slab_refl` and `TPSF_DA_paralpip_refl` functions.
## Documentation

## LightPropagation community

Currently none.

## Contributing to LightPropagation

LightPropgation is a collaborative project open to contributions and discussion.

## How to cite LightPropagation

In order to give credit to the `LightPropagation` contributors, we ask you to cite the reference below in any publication in which you have made use of `LightPropagation` packages:

## Contact

Please contact the project adminstrators, [Michael Helton](mailto:heltonmc@umich.edu), for further questions on usage and contributing. 

