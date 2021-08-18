# Getting Started

Please use the installation section to download and install Julia and LightPropagation.jl. Once downloaded we can open the Julia REPL and type...
```julia
using LightPropagation
```
The initial download and compilation could take same time but should only take a couple seconds afterwards. 

## Forward Simulation

Let's simulate the steady-state (continuous-wave) fluence from an isotropic point source in an infinite medium using `fluence_DA_inf_CW`.
```julia
julia> fluence_DA_inf_CW(1.0, 0.1, 10.0)
0.4223682678487986
```

This is a bit confusing because we don't know what `1.0, 0.1, 10.0` represent. To get more information about a function's inputs and methods type `?` in the REPL and then the name of the function...

```julia
julia>?
help?> fluence_DA_inf_CW
search: fluence_DA_inf_CW fluence_DA_semiinf_CW fluence_DA_inf_TD

  fluence_DA_inf_CW(ρ, μa, μsp)

  Compute the steady-state fluence in an infinite medium.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  ρ: source-detector separation (cm)

    •  μa: absorption coefficient (cm⁻¹)

    •  μsp: reduced scattering coefficient (cm⁻¹)

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> fluence_DA_inf_CW(1.0, 0.1, 10.0)
```

This is a convenient way to see a brief description of the function and its arguments. So now in our previous example `1.0` corresponds to `ρ` which is the source-detector separation in centimeters. `0.1` corresponds to `μa` which is the absorption coefficient and `10.0` corresponds to `μsp` which is the reduced scattering coefficient.

If we wanted to simulate this for many `ρ` values we can do...

```julia
julia> ρ = 0.1:0.1:1.0
0.1:0.1:1.0

julia> fluence_DA_inf_CW.(ρ, 0.1, 10.0)
10-element Vector{Float64}:
 20.076563644369305
  8.441844991553953
  4.732864855014603
  2.985130735997168
  2.0083126892986254
  1.4074341205286183
  1.0145168743692066
  0.7465266519850676
  0.5580470079778079
  0.4223682678487986
```

Notice the `.` after the function call to `fluence_DA_inf_CW` which implies broadcasting the function over each `ρ` value. We could also use `map`...
```julia
julia> map(ρ -> fluence_DA_inf_CW(ρ, 0.1, 10.0), 0.1:0.1:1.0)
10-element Vector{Float64}:
 20.076563644369305
  8.441844991553953
  4.732864855014603
  2.985130735997168
  2.0083126892986254
  1.4074341205286183
  1.0145168743692066
  0.7465266519850676
  0.5580470079778079
  0.4223682678487986
```

We can also do this for a single `ρ` value for several different permutations of `μsp` and `μa` with both approaches...

```julia
julia> fluence_DA_inf_CW.(1.0, [0.1, 0.2], [10.0, 12.1])
2-element Vector{Float64}:
 0.4223682678487986
 0.1952166673254022

 julia> map((μa, μsp) -> fluence_DA_inf_CW(1.0, μa, μsp), [0.1, 0.2], [10.0, 12.1])
2-element Vector{Float64}:
 0.4223682678487986
 0.1952166673254022
```

The `x -> f(x...)` creates an anonymous function that can be advantageous when passing functions.

Sometimes we want a computation to be done with higher precision. Take a slightly unrealistic example where we want the fluence way away from the source for large scattering coefficients......
```julia
julia> fluence_DA_inf_CW(100.0, 1.0, 10000.0)
0.0
```
Well this is slightly unexpected because we wouldn't expect this to be quite `0.0`. This is due to using `Float64` numbers which have a finite minimum `floatmin(Float64) = 2.2250738585072014e-308`. We can do this calculation in higher precision by converting one of the inputs to higher precision...
```julia
julia> fluence_DA_inf_CW(100.0, 1.0, big(10000.0))
1.502554931959911174198461266529831255882290472488206080785039554131741898622508e-7521

julia> floatmin(BigFloat)
8.50969131174083613912978790962048280567755996982969624908264897850135431080301e-1388255822130839284
```

The previous example works by promoting the types of the other inputs and doing calculations in the higher precision. It is recommended that all inputs be of the same type to make sure some of the calculation isn't done in a lower precision...

Perhaps we were interested in something slightly more interesting such as the time-domain fluence in a semi-infinite medium...
```julia
julia> fluence_DA_semiinf_TD(1.0, 1.0, 0.1, 10.0)
0.000288659079311426
```

More parameters to decipher... let's check the docs!
```julia
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

So now in addition to `ρ`, `μa`, `μsp` we have the `t` time input. Some of the parameters follow `;` symbol meaning they are keyword arguments. We can either explicitly define them with their name...
```julia
julia> fluence_DA_semiinf_TD(0.1:0.3:1.0, 1.0, 0.1, 10.0, n_ext = 1.0)
4-element Vector{Float64}:
 0.13275269870908793
 0.011700696523308796
 0.0015494334532460989
 0.000288659079311426
```
Or if they aren't defined they utilize their default value. Notice that we don't need the `.` here because methods in the time-domain are constructed to more efficiently compute values at many time points. Therefore, there are methods for both `t::AbstractFloat` and `t::AbstractVector`.

## Simple Plotting

We can compare semi-infinite to slab models in the time domain like the below example...
```julia
julia> t = 0.01:0.01:4.0; ρ = 1.0; μa = 0.1; μsp = 10.0
julia> semiinf = fluence_DA_semiinf_TD(t, ρ, μa, μsp)
julia> slab = fluence_DA_slab_TD(t, ρ, μa, μsp, s = 1.0, xs = 20) # s is the slab thickness

julia> using Plots # make sure to add the Plots package
julia> plot(t, log10.(semiinf), lw = 2, ylabel = "fluence (1/cm^2)", xlabel = "time (ns)", label = "Semi-infinite")
julia> plot!(t, log10.(slab), lw = 2, label = "Slab")
```

There are several other forward models in different geometries. The naming scheme of the forward models follow a PhysicalQuantity\_Approximation\_Geometry\_MeasurementType scheme. Check out the index below:

## Index

```@index
Pages = ["API.md"]
```
