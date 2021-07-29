# LightPropagation.jl


| **Documentation** |
|:------------ |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://heltonmc.github.io/LightPropagation.jl/dev/) |
| [![Join the chat at https://gitter.im/LightPropagation/community?source=orgpage](https://badges.gitter.im/LightPropagation/community.svg)](https://gitter.im/LightPropagation/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)|


## Purpose

LightPropagation provides a set of tools to model and analyze the propagation of light in turbid media written in the [Julia programming language](https://julialang.org/).
The current library supports simulating the time-resolved reflectance in the forward direction under the diffusion approximation in the following geometries:
- __Semi-infinite__ 
- __Slabs__
- __Parallelepipeds__
- __Cylinders__
- __N-layer Cylinder__

## Usage

The easiest way to install julia is by [downloading julia](https://julialang.org/downloads/) from the offical site and following the [platform specific installations](https://julialang.org/downloads/platform/). 

Launch Julia and open Julia's package manager by typing `]` in the REPL. You should see the command line change from `julia>` to `(@v1.5) pkg>`. (The @v1.x will display your current version) Once the package has cloned, backspace to bring back the `julia>` in your REPL and type `using LightPropagation` as shown below.

```julia
(@v1.5) pkg> add "https://github.com/heltonmc/LightPropagation.git"

julia> using LightPropagation

```
## Tutorials

You should now have access to all of the exported functions of the package.

### N-layered media

The n-layered solutions take the structure `Nlayer_cylinder` as an input argument. We can define the input structure as follows with default arguments:

```julia
julia> data = Nlayer_cylinder()
Nlayer_cylinder{Float64}
  μsp: Array{Float64}((4,)) [10.0, 10.0, 10.0, 10.0]
  μa: Array{Float64}((4,)) [0.1, 0.1, 0.1, 0.1]
  n_ext: Float64 1.0
  n_med: Array{Float64}((4,)) [1.0, 1.0, 1.0, 1.0]
  l: Array{Float64}((4,)) [0.5, 0.8, 1.0, 5.0]
  ρ: Float64 1.0
  a: Float64 5.0
  ω: Float64 0.0  
```
`ρ` is the source-detector separation in cylindrical coordinates,  `a` is the cylinder radius, `l` is the thickness of each layer.
If we wanted to explicitly define our optical properties for a 2-layered media we can define our input structure using keyword arguments:
```julia
julia> data = Nlayer_cylinder(ρ = 1.0, μsp = [10.0, 11.1], μa = [0.1, 0.18], l = [1.2, 10.0], n_med = [1.0, 1.0])
Nlayer_cylinder{Float64}
  μsp: Array{Float64}((2,)) [10.0, 11.1]
  μa: Array{Float64}((2,)) [0.1, 0.18]
  n_ext: Float64 1.0
  n_med: Array{Float64}((2,)) [1.0, 1.0]
  l: Array{Float64}((2,)) [1.2, 10.0]
  ρ: Float64 1.0
  a: Float64 5.0
  ω: Float64 0.0
```

#### Steady-State
To simulate the fluence in the steady-state we can simply run:
```julia
julia> data = Nlayer_cylinder()
julia> fluence_DA_Nlay_cylinder_CW(data, besselroots[1:600])
0.02403599777515162
```
The values default to the same optical properties in both layers to allow for comparison to homogoenous media. For example, to compare to the semi-infinite approximation you can run:
```julia
julia> fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0) # inputs are (ρ, [μa, μsp], n_ext, n_med, z)
0.024035998288332954
```
A few notes: `besselroots` contains the first 1000000 roots of `J0`. Generally `<1000` roots are needed to reach adequate convergence. In the above example you can see that the semi-inf and layered solution have `rel_error ~ 1e-7`. You can increase the size of the cylinder radius and layer thickness to better approximate semi-infinite. Generally, the number of roots needed to reach convergence is directly related to the inputs `a`, `l`, and `μsp`. 
It is recommended to keep `a` and `l` reasonable for your applications. `a = 10` and `l = [1.0, 10.0]` will approximate an infinite bottom layer and laterally infinite sides well. If you have large scattering coefficients you will need to increase number of roots. All units are in `cm`.
Due to floating point arithmetic and the limited precision of the besselroots, you can expect absolute errors on the order of `~1e-14` once converged. This means that to simulate a fluence value less than `1e-15`, you will need to use higher precision calculations. Though, for practical measurements such low values aren't usually detectable.
#### Time-domain
To simulate the fluence in the time-domain we can simply run:
```julia
julia> t = range(0.03, 5.0, length = 120) # nanoseconds
julia> phi_lay = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600])
## to compare to semi-infinite 
julia> SI = fluence_DA_semiinf_TD(t, [0.1, 10.0], 1.0, 1.0, 1.0, 0.0)

## you can then plot like
using Plots
julia> plot(t, log10.(abs.(phi_lay)), label="layered solution")
julia> plot!(t, log10.(SI), label="semi-inf")
```

So what is happening after `4ns`? The layered solution can take two keyword arguments: `besselroots` which is the number of roots to consider in the sum and `N` which is the number of Hankel-Laplace calculations in the inversion from the Hankel-Laplace space to the time-domain. The number of calculations needed depends on the dynamic range of the time-domain signal you want. Try different values of N like below and see how the signal is reconstructed:
```julia
phi_lay12 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 12)
phi_lay52 = fluence_DA_Nlay_cylinder_TD(t, cylinder_data, bessels = besselroots[1:600], N = 52)
```
Again, we are limited by the precision of the besselroots in our calculation. Utilizing double precision we are not able to generate fluence values below the machine precision `~2.2e-16`. If you want values lower than this you can use quad precision calculation to generate lower values. This occurs in the time-domain much more frequently than the spatial domain at long times and high absorption values.

#Performance notes#: For the best performance you will need to start julia in your terminal with multiple threads. See https://docs.julialang.org/en/v1/manual/multi-threading/ for start up help. 

### Calculation in higher precision
As mentioned previously, we are limited to absolute errors on the order of machine precision. So we can't expect to simulate any fluence values below eps. To simulate in arbitrary precision (down to at least absolute errors ~`1e-77` we need to load bessel roots in arbitrary precision. Unzip the file `besselzeroroots_big.jld.zip` located in the `src/forwardmodels/Diffusion Approximation/` folder. You will then need to load the file like:
```julia
using JLD
bessel_arb = load("besselzeroroots_big.jld")["big_besselroots"]
cylinder_data = Nlayer_cylinder{BigFloat}(ρ = 15.0, μsp = [20.0,20.0,20.0,20.0],μa = [0.3,0.3,0.3,0.3], a = 25.0, ω = 0.0, l = [1.0, 2.0, 3.0, 10.0]) # convert inputs to arbitrary floats
julia> cyl = fluence_DA_Nlay_cylinder_CW(cylinder_data, bessel_arb[1:12000])
1.166981416744208493461258014108689524878588065968633120335354050360868699344118e-31
julia> fluence_DA_semiinf_CW(15.0, [0.3, 20.0], 1.0, 1.0, 0.0) ### check against semi-infinite solution
1.1669489849622033e-31
```
Using arbitrary floats requires signficinatly more time to compute. Using Float128 values gives absolute errors at eps `~1.92592994438723e-34` while being much faster. Example below:
```julia
using Quadmath # gives us Float128 values
cylinder_data = Nlayer_cylinder{Float128}(ρ = 15.0, μsp = [20.0,20.0,20.0,20.0],μa = [0.3,0.3,0.3,0.3], a = 25.0, ω = 0.0, l = [1.0, 2.0, 3.0, 10.0]) # convert inputs to Float128 
julia> cyl = fluence_DA_Nlay_cylinder_CW(cylinder_data, Float128.(bessel_arb[1:12000]))
```

(OLD VERSION/DEPRECATED)
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

## LightPropagation community

Join the [gitter](https://gitter.im/LightPropagation/community) chat to ask questions and interact with the LightPropagation.jl community.

## Contributing to LightPropagation

LightPropgation is a collaborative project open to contributions and discussion.

## Contact

Please contact the project adminstrators, [Michael Helton](mailto:heltonmc@umich.edu), for further questions on usage and contributing. 

