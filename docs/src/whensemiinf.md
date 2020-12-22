# _When is semi-infinite semi-infinite?_

### Learning outcomes:
- Simulate the time-domain reflectance under the diffusion approximation for **semi-finite, slab, and parralelepiped** geometries
- Utilize the non-linear curve fitting functions to extract optical properties from simulated experimental data
- Gain physical intuition of approximations and how to apply them

!!! info "Disclaimer"
    The following is not a recommendation on when to apply certain approximations as that would depend on     experimental setup and application. Please use this as an introduction to the software.

##### Loading LightPropagation.jl
Let's begin by loading the package into the current environment. This can be done as below or as shown on the main page.

```julia
using Pkg
Pkg.add(url = "https://github.com/heltonmc/LightPropagation.git")
using LightPropagation
```

## Forward Models

The functions we are using to simulate the TPSF for the semi-infinite, slab, and parallelpiped geometries are:
- **TPSF\_DA\_semiinf\_refl**
- **TPSF\_DA\_slab_refl**
- **TPSF\_DA\_paralpip_refl**

!!! tip "Available Functions"
    To view all the available functions is is often useful to check the exported functions in the main file. For example if you go to ~/LightPropagation/src/LightPropagation.jl in the repository you will find `export TPSF_DA_semiinf_refl` lines that indicate available functions.

### Semi-infinite

Let's simulate the TPSF for a semi-infinite medium under the following conditions:
- t = 0:0.01:5 (ns)
- β = [0.1, 10.0] -- β[1] and β[2] are the absorption and reduced scattering coefficient (1/cm)
- ρ = 1.0 -- source-detector separation (cm)
- ndet = 1.0 -- index of refraction of detector
- nmed = 1.0 -- index of refraction of medium

```julia
t = 0.01:0.01:5
TPSF_semiinf = TPSF_DA_semiinf_refl(t, [0.1,10.0], 1.0, 1.0, 1.0)
```	
### Slab

Let's simulate the TPSF for a slab (bounded in the z-direction) for two thicknesses (2cm and 8cm) with the same optical properties and SDS as simulated for the semi-infinite medium. The inputs are the same as the semi-infinite geometry with an additional input for the slab thickness.

```julia
TPSF_slab_2cm = TPSF_DA_slab_refl(t, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
TPSF_slab_8cm = TPSF_DA_slab_refl(t, [0.1,10.0], 1.0, 1.0, 1.0, 8.0)
```

### Parallelepiped

The parralelepiped geometry is bounded in the x, y, and z direction. Check the difference in input arguments by typing `? TPSF_DA_paralpip_refl` into the REPL. You can see that the inputs are slightly different than the slab and semi-infinite geometry. Instead of listing the SDS we need to specify the location of the source and detector in x, y coordinates. Here, rd is the x, y location of the source on the top surface, and rs is the location of the detector on the top source. These are measured from the corner. Let's have two cubes one with L = 8cm and L = 4cm with the source located in the center so that rs = [4.0, 4.0] and [2.0, 2.0] respectively. We want the detector to be 1cm from the source to match SDS in previous examples. Below they are defined as rd = [3.0, 4.0] and [1.0, 2.0].        

```julia
TPSF_paralpip_8cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [3.0,4.0], [4.0,4.0], [8.0,8.0,8.0])
TPSF_paralpip_4cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [1.0,2.0], [2.0,2.0], [4.0,4.0,4.0])
```

#### Compare all three geometries

Now let's compare and visualize all 5 of the curves we have simulated as shown below:

```julia
using Plots
plot(t, TPSF_semiinf, yscale=:log10, lw = 10, ylabel="Counts", xlabel="time (ns)", label="TPSF_semiinf")
plot!(t, TPSF_slab_2cm, lw = 2, label = "TPSF_slab_2cm")
plot!(t, TPSF_slab_8cm, lw = 4, label = "TPSF_slab_8cm", c="black")
plot!(t, TPSF_paralpip_8cm, lw = 2, label = "TPSF_paralpip_8cm", alpha = 0.9)
plot!(t, TPSF_paralpip_4cm, lw = 2, label = "TPSF_paralpip_4cm", alpha = 0.9)
```

![modelcomparison](./assets/modelcomparison.png)

Great. Now change the input parameters for each curve (optical properties and SDS) to investigate different domains. Here, for this set of optical properties and SDS the  semi-infinite approximation matches exactly the slab solution when the slab thickness is 8cm (for t < 5ns), though it would be a poor approximation for thin samples (2cm). The parallelepiped geometry also matches the semi-infinite geometry when its sides are greater than 8cm, though there are some differences at longer time scales. In practice, experimental systems would not have enough dynamic range to fit over such large changes in intensity.

#### Key Takeaways
- Good agreement at the earliest times (<1ns) where the boundary effects are minimized. 
- For confined geometries, (e.g. TPSF\_paralpip\_4cm and TPSF\_slab\_2cm) we see there are significant differences compared to the semi-infinite approximations.
- If the boundaries are significantly away from the source the semi-infinite approximation appears to be a good approximation.
- If the underlying media is a parralelepiped, would a semi-infinite geometry under or overestimate the optical properties? Hint: Changes in absorption mainly affect the later times, while scattering dominates early arrival times.

## Inverse Problem

Now that we have visually observed the differences between the three geometries, let's investigate and try to quantify the error when using the models in the inverse problem to extract optical properties.

The first thing we must do is simulate some experimental data to fit. In practice, samples are not infinite so let's use the parralelepiped geometry as the ground 'truth' and let's fit all three models to that curve. You can add noise (follows a poisson distribution) if you would like, but let's assume no noise for now.

```julia
ydata = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [0.5,1.5], [1.5,1.5], [3.0,3.0,3.0])
```

##### Setup Inverse problem

The first thing we must do is set up the input data and model parameters. The input data is your experimental known parameters. Here, these will be two separate data structures with different fields. There will be required fields that you must enter and fields that will default if you don't define them. 

You must specify the time scale, the experimental DTOF, and the IRF. Depending on your geometry, you would then specify the SDS, index of refraction of medium and detector, the slab thickness, or the parallelepiped dimensions. Since we know the dimensions of the cube (assuming these are known), and because we are testing the model with all three approximations we need to specify each field. If you were just using a semi-infinite approximation you would only have to specify (t, DTOF, IRF, and ρ). 

Before we can define our data structure we must define the IRF. We must also define our model parameters which take the analytical model you are fitting with, the fit ranges, the initial guess for the inverse fitting procedures as well as upper and lower bounds. If you don't define these they will use the default values. (The model normalizes the photon counts, so we have normalized the ydata)

```julia

#Define the IRF as a delta function length of ydata
IRF = [1; zeros(length(ydata) - 1)]

# Define a structure called data with the specified fields 
data = fitDTOF(t = t, DTOF = ydata./maximum(ydata), IRF = IRF, ρ = 1.0, nmed = 1.0, ndet = 1.0, s = 3.0, rd = [0.5,1.5], rs = [1.5,1.5], L = [3.0,3.0,3.0])

# Define the fitting parameters called m for parallepiped geometry
m = DTOF_fitparams(model = TPSF_DA_paralpip_refl, risefactor = 0.01, tailfactor = 1e-10,initparams = [0.05, 20.0])

# Define model with semi-inf  approximation
m1 = DTOF_fitparams(model = TPSF_DA_semiinf_refl, risefactor = 0.01, tailfactor = 1e-10, initparams = [0.05, 20.0])

#Define model with slab approximation
m2 = DTOF_fitparams(model = TPSF_DA_slab_refl, risefactor = 0.01, tailfactor = 1e-10, initparams = [0.05, 20.0])

#Fit data with the three models
f, tfit, yraw, yfit  = getfit(data, m)
f1, tfit1, yraw1, yfit1  = getfit(data, m1)
f2, tfit2, yraw2, yfit2  = getfit(data, m2)
```
As shown the `getfit` function returns four separate variables. f will be the fit structure where you can call the fitting parameters and residuals, t will be the time scale the fit happened on, and the yraw data that was fit to, and the resulting yfit.

Let's plot these below to see how they look.

```julia
scatter(tfit, yraw, m="*", ms = 8, label = "ydata", c="black")
plot!(tfit, yfit, lw = 3, alpha = 1, label = "paralpip model", c="blue")
plot!(tfit1, yfit1, lw = 3, alpha = 1,label = "semiinf model", c="red")
plot!(tfit2, yfit2,lw = 3, alpha = 1, label = "slab model", c="orange")
plot!(ylabel="Counts", xlabel="time (ns)")
```

![fitcomparison](./assets/fitcomparison.png)

To check the extracted optical properties we can call the param field of our fit result like:

```julia
f.param # [0.0999999 9.999999]
f1.param # [0.158778 13.51499]
f2.param # [0.153302 12.61908]
```
We can see that when using the semi-infinite approximation (f1) results in 59% error in the extracted absoprtion coefficient and 35% error in the extracted reduced scattering coefficient. The error when using the slab approximation (f2) is 53% and 26% respectively. Utilizing the parallelepiped geometry returned the exact properties. Remember that we are fitting over a large range where the difference at long time scales can be significant. Most of the time fitting only occurs at 1 or 0.1% on the falling tail. Change the fitting parameters in your model to see how it can improve or decrease in accuracy depending on fit range.

_Note_: If you are just wanting to fit using a semi-infinite or slab geometry you do not need to specify the data structure as defined earlier. Whatever model you define in the parameter structure will call the appropriate inputs. For example if you just wanted to fit a semi-infinite model you can define data and m as below:

```julia
data = fitDTOF(t = t, DTOF = ydata./maximum(ydata), IRF = IRF, ρ = 1.0, nmed = 1.0, ndet = 1.0)
m = DTOF_fitparams(model = TPSF_DA_semiinf_refl, risefactor = 0.01, tailfactor = 1e-10,initparams = [0.05, 20.0])
```