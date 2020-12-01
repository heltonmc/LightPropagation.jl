# Solutions of the Diffusion Equation for the Slab Geometry

The following describes the light propagation through turbid media bounded by parallel planes. 

The described derivations follow the methods from [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587).[^1]

## Diffusion Equation

```math
(\frac{1}{\nu}\frac{\partial }{\partial t} - D\nabla^2 + \mu_a)\Phi(\vec{r}, t) = Q(\vec{r}, t)
```

Where (``Q(\vec{r}, t)``) is the isotropic source term and D is the diffusion coefficient, (''D = \frac{1}{3\mu_s'}``

### Solution of the Diffusion Equation for homogeneous media in a slab geometry 

The diffusion equation is a partial-differential equation that will require boundary conditions to solve for a specific geometry. Here we will utilize the
extrapolated boundary conditions as described in [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587).[^1] that assumes that the flux is equal to 0
on an extrapolated surface at a distance of (``2AD``). After the boundary conditions have been determined, utilizng the method of images will allow us to reconstruct the fluence
inside the medium. 

Utilizing Equation 33 from [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587),[^1] the time-dependent Green's function for the fluence rate at (``\vec{r}``) can be described by:


```math
\Phi(\vec{r}, t) = \frac{\nu}{(4\pi D \nu t)^\frac{3}{2}}
exp(- \frac{\rho ^2}{4 D \nu t} - \mu_a \nu t) \times
\sum_{m=-\infty}^{m=+\infty} \{exp(-\frac{(z-z_m^+)^2}{4 D \nu t}) - 
exp(-\frac{(z-z_m^-)^2}{4 D \nu t})\}
```

Where:
- ``z_m^+ = 2m(s+ 2z_e)``
- ``z_m^- = 2m(s +2z_e) - 2z_e - z_s``
- ``m = 0, \pm 1, \pm 2, ...... \pm \infty``

To obtain solutions for the time-dependent transmittance and reflectance we can utilize Fick's law where:


```math
R (\rho, t) = D\frac{\partial}{\partial z} \Phi(\rho, z = 0, t)
```
and 

```math
T (\rho, t) = - D\frac{\partial}{\partial z} \Phi(\rho, z = s, t)
```

Which yields Equation 36 from [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587)[^1] for the time-dependent reflectance on the surface:
```math
R(\rho, t) = - \frac{exp(-\mu_a \nu t - \frac{\rho^2}{4 D \nu t})}{2(4\pi D \nu)^\frac{3}{2}t^\frac{5}{2}}
\times
\sum_{m=-\infty}^{m=+\infty} [z_3_,_mexp(-\frac{z_3_,_m^2}{4 D \nu t}) - 
z_4_,_mexp(-\frac{z_4_,_m^2}{4 D \nu t})]
```
and Equation 39 from [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587).[^1] for the time-dependent transmittance on at the distance z=s where s is the thickness of the slab:
```math
T(\rho, t) = \frac{exp(-\mu_a \nu t - \frac{\rho^2}{4 D \nu t})}{2(4\pi D \nu)^\frac{3}{2}t^\frac{5}{2}}
\times
\sum_{m=-\infty}^{m=+\infty} [z_1_,_mexp(-\frac{z_1_,_m^2}{4 D \nu t}) - 
z_2_,_mexp(-\frac{z_2_,_m^2}{4 D \nu t})]
```

Where:
- ``z_1_,_m = s(1-2m) - 4mz_e - z_o``
- ``z_2_,_m = s(1-2m) - (4m-2)z_e - z_o``
- ``z_3_,_m = -2ms - 4mz_e - z_o``
- ``z_4_,_m = -2ms - (4m-2)z_e - z_o``


To obtain solution
```julia
λ0 = 1 #doesn't matter since everything is normalized to λ0
k0 = 2π/λ0
kin = 3k0
θ_i = 0.0 #incident wave is left->right
pw = PlaneWave(θ_i)
N = 260
P = 10
shapes = [rounded_star(0.1λ0, 0.05λ0, 5, N)]
ids = [1] # the particle at centers[1,:] has the parametrization shapes[ids[1]]
centers = [0.0 0.0] # our particle is centered at the origin
φs = [0.0] #zero rotation angle
sp = ScatteringProblem(shapes, ids, centers, φs)
```

[^1]: Daniele Contini, Fabrizio Martelli, and Giovanni Zaccanti, "Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory," Appl. Opt. 36, 4587-4599 (1997) 
