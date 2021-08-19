# _When is semi-infinite semi-infinite?_

### Learning outcomes:
- Simulate the time-domain fluence under the diffusion approximation for **semi-finite, slab, and parralelepiped** geometries
- Utilize the non-linear curve fitting functions to extract optical properties from simulated experimental data
- Gain physical intuition of approximations and how to apply them

!!! info "Disclaimer"
    The following is not a recommendation on when to apply certain approximations as that would depend on experimental setup and application. Please use this as an introduction to the software.

##### Loading LightPropagation.jl
Let's begin by loading the package into the current environment. This can be done as below or as shown on the main page.

```julia
using LightPropagation
```

## Forward Models

The functions we are using to simulate the TPSF for the semi-infinite, slab, and parallelpiped geometries are:
- **fluence\_DA\_semiinf\_TD**
- **fluence\_DA\_slab_TD**
- **fluence\_DA\_paralpip_TD**

### Semi-infinite

Let's simulate the fluence for a semi-infinite medium under the following conditions in the time-domain:
- t = 0.01:0.01:5 (ns)
- μa = 0.1 -- absorption coefficient (1/cm)
- μsp = 10.0 -- reduced scattering coefficient (1/cm)
- ρ = 1.0 -- source-detector separation (cm)
- n_ext = 1.0 -- index of refraction of external medium
- n_med = 1.0 -- index of refraction of medium

```julia
t = 0.01:0.01:5
semiinf = fluence_DA_semiinf_TD(t, 1.0, 0.1, 10.0, n_ext = 1.0, n_med = 1.0)
```	
### Slab

Let's simulate the time-domain fluence for a slab (bounded in the z-direction) for two thicknesses (2cm and 8cm) with the same optical properties and SDS as simulated for the semi-infinite medium. The inputs are the same as the semi-infinite geometry with an additional input for the slab thickness.

```julia
slab_2cm = fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, s = 2.0) # n_ext and n_med default to 1.0, z defaults to 0.0
slab_8cm = fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, s = 8.0)
```

### Parallelepiped

The parralelepiped geometry is bounded in the x, y, and z direction. Check the difference in input arguments by typing `? fluence_DA_paralpip_TD` into the REPL. You can see that the inputs are slightly different than the slab and semi-infinite geometry. Instead of listing the SDS we need to specify the location of the source and detector in x, y, z coordinates. Here, rd is the x, y, z location of the source on the top surface, and rs is the location of the detector on the top source. These are measured from the corner. Let's have two cubes one with L = 8cm and L = 4cm with the source located in the center so that rs = [4.0, 4.0] and [2.0, 2.0] respectively. We want the detector to be 1cm from the source to match SDS in previous examples. Below they are defined as rd = [3.0, 4.0, 0.0] and [1.0, 2.0, 0.0].        

```julia
paralpip_8cm = fluence_DA_paralpip_TD(t, 0.1, 10.0, rd = [3.0, 4.0, 0.0], rs = [4.0,4.0], L = [8.0,8.0,8.0])
paralpip_4cm = fluence_DA_paralpip_TD(t, 0.1, 10.0, rd = [1.0, 2.0, 0.0], rs = [2.0, 2.0], L = [4.0,4.0,4.0])
```

#### Compare all three geometries

Now let's compare and visualize all 5 of the curves we have simulated as shown below:

```julia
using Plots
plot(t, semiinf, yscale=:log10, lw = 10, ylabel = "fluence (1/cm²)", xlabel = "time (ns)", label = "semiinf")
plot!(t, slab_2cm, lw = 2, label = "lab_2cm")
plot!(t, slab_8cm, lw = 4, label = "slab_8cm", c="black")
plot!(t, paralpip_8cm, lw = 2, label = "paralpip_8cm", alpha = 0.9)
plot!(t, paralpip_4cm, lw = 2, label = "paralpip_4cm", alpha = 0.9)
```

![modelcomparison](./assets/modelcomparisonTD.png)

Great. Now change the input parameters for each curve (optical properties and SDS) to investigate different domains. Here, for this set of optical properties and SDS the  semi-infinite approximation matches exactly the slab solution when the slab thickness is 8cm (for t < 5ns), though it would be a poor approximation for thin samples (2cm). The parallelepiped geometry also matches the semi-infinite geometry when its sides are greater than 8cm, though there are some differences at longer time scales. In practice, experimental systems would not have enough dynamic range to fit over such large changes in intensity.

#### Key Takeaways
- Good agreement at the earliest times (<1ns) where the boundary effects are minimized. 
- For confined geometries, we see there are significant differences compared to the semi-infinite approximations.
- If the boundaries are significantly away from the source the semi-infinite approximation appears to be a good approximation.
- If the underlying media is a parralelepiped, would a semi-infinite geometry under or overestimate the optical properties? Hint: Changes in absorption mainly affect the later times, while scattering dominates early arrival times.

## Inverse Problem

Now that we have visually observed the differences between the three geometries, let's investigate and try to quantify the error when using the models in the inverse problem to extract optical properties.

The first thing we must do is simulate some experimental data to fit. In practice, samples are not infinite so let's use the parralelepiped geometry as the ground 'truth' and let's fit all three models to that curve. You can add noise (follows a poisson distribution) if you would like, but let's assume no noise for now.

```julia
ydata = fluence_DA_paralpip_TD(t, 0.1, 10.0, rd = [0.5, 1.5, 0.0], rs = [1.5, 1.5], L = [3.0, 3.0, 3.0])
```

##### Setup Inverse problem

The first thing we must do is set up the input data and models. See the Inverse Fit tutorial for more information!

```julia
julia> t_arr = 0.02:0.01:4.0
julia> ydata = log10.(fluence_DA_paralpip_TD(t_arr, 0.1, 10.0, rd = [0.5, 1.5, 0.0], rs = [1.5, 1.5], L = [3.0, 3.0, 3.0]))
julia> ub = [1.0, 50.0] # upper bound for [μa, μsp]
julia> lb = [0.01, 3.0] # lower bound for [μa, μsp]
julia> p0 = [0.12, 18.2] # some initial guess for our algorithm

julia> model_semiinf(t, β) = log10.(fluence_DA_semiinf_TD(t, ρ, β[1], β[2]))
julia> model_slab(t, β) = log10.(fluence_DA_slab_TD(t, ρ, β[1], β[2], s = 3.0))
julia> model_paral(t, β) = log10.(fluence_DA_paralpip_TD(t, β[1], β[2], rd = [0.5, 1.5, 0.0], rs = [1.5, 1.5], L = [3.0, 3.0, 3.0]))

julia> fit_semiinf = curve_fit(model_semiinf, t_arr, ydata, p0, lower=lb, upper=ub)
julia> fit_slab = curve_fit(model_slab, t_arr, ydata, p0, lower=lb, upper=ub)
julia> fit_paral = curve_fit(model_paral, t_arr, ydata, p0, lower=lb, upper=ub)
```

```julia
scatter(t_arr, ydata, ms = 8, label = "ydata", c="black", ylabel = "fluence (1/cm²)", xlabel = "time (ns)")
plot!(t_arr, model_paral(t_arr, fit_paral.param), lw = 3, alpha = 1, label = "paralpip model", c="blue")
plot!(t_arr, model_semiinf(t_arr, fit_semiinf.param), lw = 3, alpha = 1,label = "semiinf model", c="red")
plot!(t_arr, model_slab(t_arr, fit_slab.param),lw = 3, alpha = 1, label = "slab model", c="orange")
```
![fitcomparison](./assets/fitcomparisonTD.png)

To check the extracted optical properties we can call the param field of our fit result like:
```julia
julia> fit_paral.param
2-element Vector{Float64}:
  0.1000000000000203
 10.000000000002833

julia> fit_slab.param
2-element Vector{Float64}:
 0.13805160496964075
 8.432187697059568

julia> fit_semiinf.param
2-element Vector{Float64}:
  0.14792469689680052
 10.222734757205624
```
