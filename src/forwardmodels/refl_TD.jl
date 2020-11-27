

##### Diffusion approximation for semi-infinite media  #####

#=

@article{contini1997photon,
  title={Photon migration through a turbid slab described by a model based on diffusion approximation. I. Theory},
  author={Contini, Daniele and Martelli, Fabrizio and Zaccanti, Giovanni},
  journal={Applied optics},
  volume={36},
  number={19},
  pages={4587--4599},
  year={1997},
  publisher={Optical Society of America}
}

The reflectance for the semi-infinite medium is calculated from Equation 36 from Contini 1997. 
The semi-infinite medium is described by retaining the first of the dipole sources (m=0) from Eqn. 36.
=#

#Reflectance semi-infinite
function DA_semiinf(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)


    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed

    Rt1 = Array{Float64}(undef, length(t))
    Rt2 = Array{Float64}(undef, length(t))


	if n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	elseif n < 1.0
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	else 
		A = 1.0
	end

	zs::Float64 = 1/μsp
	ze::Float64 = 2A*D

	z3m::Float64 = - zs
    z4m::Float64 = 2ze +zs

   	Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
   	Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

   	Rt2 = @. z3m*exp(-(z3m^2/(4D*ν*t))) - z4m*exp(-(z4m^2/(4D*ν*t)))

   	Rt = @. Rt1*Rt2
    

   	replace!(Rt, NaN => 0)
    return Rt

end

##### Diffusion approximation for Slab Geometry   #####
#The reflectance for the slab is calculated from Equation 36 from Contini 1997 for the first 11 dipoles (-10:1:10).
#The number of sources considered in the summation can be changed with xs. Larger values of m will be needed for large values of SDS and t.
   
function DA_reflslab(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness

	n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	z3m = Array{Float64}(undef, length(xs))
	z4m = Array{Float64}(undef, length(xs))


    Rt1 = Array{Float64}(undef, length(t))
    Rt2 = Array{Float64}(undef, length(t), length(z3m))
	
	if n == 1.0
		A= 1.0
	elseif n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end

	zs::Float64 = 1/μsp
	ze::Float64 = 2A*D
	
    z3m = Float64[-2m.*s .- 4m.*ze .- zs for m in xs]
	z4m = Float64[-2m.*s .- (4m .- 2).*ze .+ zs for m in xs]
	
    
    Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
   	Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))


	Rt2 = @. z3m'*exp(-(z3m'^2 / (4D*ν*t))) - z4m'*exp(-(z4m'^2 / (4D*ν*t)))

    Rt = Rt1.*sum(Rt2, dims=2)
	replace!(Rt, NaN => 0)

	return Rt
end


##### Diffusion approximation for a turbid parallelepiped (reflectance)  #####



#=
@article{kienle2005light,
  title={Light diffusion through a turbid parallelepiped},
  author={Kienle, Alwin},
  journal={JOSA A},
  volume={22},
  number={9},
  pages={1883--1888},
  year={2005},
  publisher={Optical Society of America}
}

The fluence is given by equation 6 in Kienle 2005. The reflectance as calculated here uses Fick's law where Rt(rho, t) = D∂/∂zϕ(rho, z=0, t)
=#

function DT_refl_paralpip(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})


	# rd is location of detector [x,y]
	# rs is location of source [x,y]
	# L is dimensions of parallelpiped [Lx, Ly, Lz]
	n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x::Float64 = rd[1]
	y::Float64 = rd[2]
	xu::Float64 = rs[1]
	yu::Float64 = rs[2]
	lx::Float64 = L[1]
	ly::Float64 = L[2]
	lz::Float64 = L[3]

	x1l = Array{Float64}(undef, length(xs))
	x2l = Array{Float64}(undef, length(xs))
	y1m = Array{Float64}(undef, length(xs))
	y2m = Array{Float64}(undef, length(xs))
	z1n = Array{Float64}(undef, length(xs))
	z2n = Array{Float64}(undef, length(xs))

	Rt1 = Array{Float64}(undef, length(t))
	Rt2 = Array{Float64}(undef, length(t), length(xs))
	Rt3 = Array{Float64}(undef, length(t), length(xs))
	Rt4 = Array{Float64}(undef, length(t), length(xs))

	if n == 1.0
		A= 1.0
	elseif n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end
	
	zo::Float64 = 1/μsp
	zb::Float64 = 2A*D

	x1l = Float64[2l*lx + 4l*zb + xu for l in xs]
	x2l = Float64[2l*lx + (4l-2)*zb - xu for l in xs]
	y1m = Float64[2m*ly + 4m*zb + yu for m in xs]
	y2m = Float64[2m*ly + (4m-2)*zb - yu for m in xs]
	z1n = Float64[2n*lz + 4n*zb + zo for n in xs]
	z2n = Float64[2n*lz + (4n-2)*zb - zo for n in xs]


	Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

	Rt2 = @. exp(-(x-x1l')^2 / (4D*ν*t)) - exp(-(x-x2l')^2 / (4D*ν*t))

	Rt3 = @. exp(-(y-y1m')^2 / (4D*ν*t)) - exp(-(y-y2m')^2 / (4D*ν*t))

	Rt4 = @. z1n'*exp(-z1n'^2 /(4D*ν*t)) - z2n'*exp(-z2n'^2 /(4D*ν*t))

	Rt = Rt1.*sum(Rt2, dims=2).*sum(Rt3, dims=2).*sum(Rt4, dims=2)
	replace!(Rt, NaN =>0)
	return Rt
end


