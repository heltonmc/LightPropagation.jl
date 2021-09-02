### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ c9aa684a-43a8-11eb-3400-b383901a7f05
begin
	using Pkg
	Pkg.add(url = "https://github.com/heltonmc/LightPropagation.git")
	using LightPropagation
end


# ╔═╡ 4399d234-43a5-11eb-2389-21beb59e41fd
md"# _When is semi-infinite semi-infinite?_

### Learning outcomes:
- Simulate the time-domain reflectance under the diffusion approximation for **semi-finite, slab, and parralelepiped** geometries
- Utilize the non-linear curve fitting functions to extract optical properties from simulated experimental data
- Gain physical intuition of approximations and how to apply them
"

# ╔═╡ 31c8aab0-43a7-11eb-3fea-417ceb3ba0b5
md"
!!! info \"Disclaimer\"
	The following is not a recommendation on when to apply certain approximations as       these will depend on experimental setup and application. Please use this as an         introduction to the software. 
"

# ╔═╡ 117d72ca-43a9-11eb-1394-e584e55e76e9
md"
##### Loading LightPropagation.jl
Let's begin by loading the package into the current environment. This can be done as below or as shown on the main page.
"

# ╔═╡ f7a02712-43a8-11eb-1cff-eb83d99a0989
md"
## Forward Models

The functions we are using to simulate the TPSF for the semi-infinite, slab, and parallelpiped geometries are:
- **TPSF\_DA\_semiinf\_refl**
- **TPSF\_DA\_slab_refl**
- **TPSF\_DA\_paralpip_refl**

"

# ╔═╡ f70fa3ee-43a9-11eb-1599-c9d55ca04787
md"
!!! tip \"Available Functions\"
	To view all the available functions it is often useful to check the exported 		functions in the main file. For example if you go to                                   ~/LightPropagation/src/LightPropagation.jl in the repository you will find             `export TPSF_DA_semiinf_refl` lines that indicate available functions.
"

# ╔═╡ acd4eee0-43ab-11eb-14f7-f95c2291d580
md"
Let's simulate the TPSF for a semi-infinite medium under the following conditions:
- t = 0:0.01:5 (ns)
- β = [0.1, 10.0] -- β[1] and β[2] are the absorption and reduced scattering coefficient (1/cm)
- ρ = 1.0 -- source-detector separation (cm)
- ndet = 1.0 -- index of refraction of detector
- nmed = 1.0 -- index of refraction of medium
"

# ╔═╡ 27dee8a8-43ab-11eb-33a2-81af1f08c66d
begin
	t = 0.01:0.01:5
	TPSF_semiinf = TPSF_DA_semiinf_refl(t, [0.1,10.0], 1.0, 1.0, 1.0)
end

# ╔═╡ a32eb024-43ab-11eb-09fa-9dc34018f56f
begin
	using Plots
	plot(t, TPSF_semiinf, yscale=:log10, lw = 2, ylabel="Counts", 							xlabel="time (ns)", label="TPSF_semiinf")
end

# ╔═╡ 5864c6da-43ad-11eb-117a-b7270071791d
md"
##### Comparison of semi-infinite and slab
Now that we know how to run a function and plot it, let's compare the semi-infinite approximation to a slab (bounded in the z-direction) for two thickness (2cm and 8cm)
"

# ╔═╡ 2a29c586-43ad-11eb-007f-45911e0d0ca5
begin
	TPSF_slab_2cm = TPSF_DA_slab_refl(t, [0.1,10.0], 1.0, 1.0, 1.0, 2.0)
	TPSF_slab_8cm = TPSF_DA_slab_refl(t, [0.1,10.0], 1.0, 1.0, 1.0, 8.0)
end

# ╔═╡ d5282152-43ae-11eb-3381-5bb1db31fece
begin
	plot(t, TPSF_semiinf, yscale=:log10, lw = 6, ylabel="Counts", 			 			xlabel="time (ns)", label="TPSF_semiinf") #make lw thicker for clarity
	plot!(t, TPSF_slab_2cm, lw = 2, label = "TPSF_slab_2cm")
	plot!(t, TPSF_slab_8cm, lw = 2, label = "TPSF_slab_8cm", c="black")
end

# ╔═╡ 5bd4961c-43b0-11eb-0d78-c3850f2f643d
md"
Great. Now change the input parameters for each curve (optical properties and SDS) to investigate different domains. Here, for this set of optical properties and SDS the  semi-infinite approximation matches exactly the slab solution when the slab thickness is 8cm (for t<5ns), though it would be a poor approximation for thin samples (2cm).
"

# ╔═╡ e98e6dae-43b3-11eb-2e8a-01012c7befe2
md"
#### Compare all three geometries
Let's introduce the parralelepiped geometry where it is bounded in the x, y, and z direction. Remember the inputs look like this:                           

TPSF\_DA\_paralpip\_refl(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})  

You can also check arguments by typing ? in the REPL and typing in the function name as shown in main page. Here, rd is the x, y location of the source on the top surface, and rs is the location of the detector on the top source. These are measured from the corner. Let's have two cubes one with L = 8cm and L = 4cm with the source located in the center so that rs = [4.0, 4.0] and [2.0, 2.0] respectively. We want the detector to be 1cm from the source to match SDS in previous examples. Below they are defined as rd = [3.0, 4.0] and [1.0, 2.0].
"

# ╔═╡ d40bc0b0-43b0-11eb-1f81-597d8d2a017b
begin
	TPSF_paralpip_8cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [3.0,4.0], 			[4.0,4.0], [8.0,8.0,8.0])
	TPSF_paralpip_4cm = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [1.0,2.0], 			[2.0,2.0], [4.0,4.0,4.0])
end

# ╔═╡ 0cb1c518-43b6-11eb-3bac-052d89f819ae
begin
	plot(t, TPSF_semiinf, yscale=:log10, lw = 10, ylabel="Counts", 			 			xlabel="time (ns)", label="TPSF_semiinf") #make lw thicker for clarity
	plot!(t, TPSF_slab_2cm, lw = 2, label = "TPSF_slab_2cm")
	plot!(t, TPSF_slab_8cm, lw = 4, label = "TPSF_slab_8cm", c="black")
	plot!(t, TPSF_paralpip_8cm, lw = 2, label = "TPSF_paralpip_8cm", alpha = 0.9)
	plot!(t, TPSF_paralpip_4cm, lw = 2, label = "TPSF_paralpip_4cm", alpha = 0.9)

end

# ╔═╡ 9a7d26ee-4479-11eb-1eb0-f79e3c07fb28
pwd()

# ╔═╡ 46710c9e-43b8-11eb-0172-65d7f2171f56
md"
#### Key Takeaways
- Good agreement at the earliest times (<1ns) where the boundary effects are minimized. 
- For confined geometries, (e.g. TPSF\_paralpip\_4cm and TPSF\_slab\_2cm) we see there are significant differences compared to the semi-infinite approximations.
- If the boundaries are significantly away from the source the semi-infinite approximation appears to be a good approximation.
"

# ╔═╡ 3cb5c3e2-43b9-11eb-2299-7503e415f13c
md"
## Inverse Problem

Now that we have visually observed the differences between the three geometries, let's investigate and try to quantify the error in the extracted optical properties when using different approximations.

The first thing we must do is simulate some experimental data. In practice, samples are not infinite so let's use the parralelepiped geometry as the ground 'truth' and let's fit all three models to that data. 
"

# ╔═╡ d2ee9b40-43eb-11eb-1f68-7f94589ba3f4
ydata = TPSF_DA_paralpip_refl(t, [0.1,10.0], 1.0, 1.0, [0.5,1.5], [1.5,1.5], [3.0,3.0,3.0])

# ╔═╡ 2d46a2d4-43b7-11eb-10fa-03d55dba74e8
md"
##### Setup Inverse problem

The first thing we must do is set up the input data and model parameters. These use two separate structures. The input data has the following structure: 
```julia
@with_kw struct fitDTOF
   
    #required inputs
    t::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    DTOF::Array{Float64,1}
    IRF::Array{Float64,1}

    nmed::Float64 = 1.5 
	ndet::Float64 = 1.51
	lambda::Float64 = 750.0
     
    #semi-infinite and slab
    ρ::Union{Float64, Missing} = missing

    #slab
    s::Union{Float64, Missing} = missing

    #parallelepiped
    rd::Union{Array{Float64,1}, Missing} = missing
    rs::Union{Array{Float64,1}, Missing} = missing
    L::Union{Array{Float64,1}, Missing} = missing
end
```
"

# ╔═╡ 27b1a654-43ec-11eb-1b4f-9f165debb017
md"
We can see that we must also describe and IRF. For now, let's assume the IRF is as defined below.
"

# ╔═╡ 6f0f914c-43c3-11eb-2561-f7c8cae30e99
IRF = [1; zeros(length(ydata) - 1)]

# ╔═╡ 369062ea-43c3-11eb-3143-a562cfe3afe5
data = fitDTOF(t = t, DTOF = ydata./maximum(ydata), IRF = IRF, ρ = 1.0, nmed = 1.0, ndet = 1.0, s = 3.0, rd = [0.5,1.5], rs = [1.5,1.5], L = [3.0,3.0,3.0])

# ╔═╡ 495b4f02-43c5-11eb-3ea6-6d35ebc81514
m = DTOF_fitparams(model = TPSF_DA_paralpip_refl, risefactor = 0.01, tailfactor = 1e-10,initparams = [0.05, 20.0])

# ╔═╡ 6508edfe-43c5-11eb-2d67-4fcc34883150
f, tfit, yraw, yfit  = getfit(data, m)

# ╔═╡ 50348fb8-43c8-11eb-2bfc-6102f19bcf0e
begin
	m1 = DTOF_fitparams(model = TPSF_DA_semiinf_refl, risefactor = 0.01, tailfactor = 1e-10, initparams = [0.05, 20.0])
	m2 = DTOF_fitparams(model = TPSF_DA_slab_refl, risefactor = 0.01, tailfactor = 		1e-10, initparams = [0.05, 20.0])
end


# ╔═╡ 6fa2a252-43c8-11eb-2ffd-73dfebb5454e
begin
	f1, tfit1, yraw1, yfit1  = getfit(data, m1)
	f2, tfit2, yraw2, yfit2  = getfit(data, m2)
end

# ╔═╡ 8790e19e-43c8-11eb-11d4-e5b36da13513
begin
	scatter(tfit, yraw, m="*", ms = 8, label = "ydata", c="black")
	plot!(tfit, yfit, lw = 3, alpha = 1, label = "paralpip model", c="blue")
	plot!(tfit1, yfit1, lw = 3, alpha = 1,label = "semiinf model", c="red")
	plot!(tfit2, yfit2,lw = 3, alpha = 1, label = "slab model", c="orange")
	plot!(ylabel="Counts", xlabel="time (ns)")
	savefig("fitcomparison.png")
end

# ╔═╡ a8e788d4-43c8-11eb-3c20-b94cd8884937
f.param

# ╔═╡ af161984-43c8-11eb-0715-079cc75433ec
f1.param

# ╔═╡ b231ab40-43c8-11eb-0819-1d16af7aafb4
f2.param

# ╔═╡ c298b8e6-4488-11eb-2b48-b9a7fcc20740
abs.(f.param .- f2.param)./f.param

# ╔═╡ Cell order:
# ╟─4399d234-43a5-11eb-2389-21beb59e41fd
# ╠═31c8aab0-43a7-11eb-3fea-417ceb3ba0b5
# ╟─117d72ca-43a9-11eb-1394-e584e55e76e9
# ╠═c9aa684a-43a8-11eb-3400-b383901a7f05
# ╟─f7a02712-43a8-11eb-1cff-eb83d99a0989
# ╟─f70fa3ee-43a9-11eb-1599-c9d55ca04787
# ╟─acd4eee0-43ab-11eb-14f7-f95c2291d580
# ╠═27dee8a8-43ab-11eb-33a2-81af1f08c66d
# ╠═a32eb024-43ab-11eb-09fa-9dc34018f56f
# ╟─5864c6da-43ad-11eb-117a-b7270071791d
# ╠═2a29c586-43ad-11eb-007f-45911e0d0ca5
# ╠═d5282152-43ae-11eb-3381-5bb1db31fece
# ╟─5bd4961c-43b0-11eb-0d78-c3850f2f643d
# ╟─e98e6dae-43b3-11eb-2e8a-01012c7befe2
# ╠═d40bc0b0-43b0-11eb-1f81-597d8d2a017b
# ╠═0cb1c518-43b6-11eb-3bac-052d89f819ae
# ╠═9a7d26ee-4479-11eb-1eb0-f79e3c07fb28
# ╟─46710c9e-43b8-11eb-0172-65d7f2171f56
# ╟─3cb5c3e2-43b9-11eb-2299-7503e415f13c
# ╠═d2ee9b40-43eb-11eb-1f68-7f94589ba3f4
# ╟─2d46a2d4-43b7-11eb-10fa-03d55dba74e8
# ╟─27b1a654-43ec-11eb-1b4f-9f165debb017
# ╠═6f0f914c-43c3-11eb-2561-f7c8cae30e99
# ╠═369062ea-43c3-11eb-3143-a562cfe3afe5
# ╠═495b4f02-43c5-11eb-3ea6-6d35ebc81514
# ╠═6508edfe-43c5-11eb-2d67-4fcc34883150
# ╠═50348fb8-43c8-11eb-2bfc-6102f19bcf0e
# ╠═6fa2a252-43c8-11eb-2ffd-73dfebb5454e
# ╠═8790e19e-43c8-11eb-11d4-e5b36da13513
# ╠═a8e788d4-43c8-11eb-3c20-b94cd8884937
# ╠═af161984-43c8-11eb-0715-079cc75433ec
# ╠═b231ab40-43c8-11eb-0819-1d16af7aafb4
# ╠═c298b8e6-4488-11eb-2b48-b9a7fcc20740
