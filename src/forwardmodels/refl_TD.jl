

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
    	zₛ::Float64 = 1/μsp
    	zₑ::Float64 = 2A*D
    	ν::Float64 = 29.9792345/nmed
	A::Float64

    	z₃ₘ::Float64 = - zₛ
	z₄ₘ::Float64 = 2zₑ +zₛ

	Rt1 = Array{Float64}(length(t))
	Rt2 = Array{Float64}(length(t))
	Rt = Array{Float64}(length(t))


	if n == 1.0
		A = 1.0
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end


   	 Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
   	 Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

   	 Rt2 = @. z₃ₘ*exp(-(z₃ₘ^2/(4D*ν*t))) - z₄ₘ*exp(-(z₄ₘ^2/(4D*ν*t)))

   	 Rt = @. Rt1*Rt2
    

   	 replace!(Rt, NaN => 0)
    return Rt

end

##### Diffusion approximation for Slab Geometry   #####
   
function DA_reflslab(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness
	n::Float64 = nmed/ndet
	μa::Float64 = β[1]
	μsp::Float64 = β[2]
   	D::Float64 = 1/3μsp
    	zₛ::Float64 = 1/μsp
	zₑ::Float64 = 2A*D
	A::Float64
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	z₃ₘ::Array{Float64}(length(xs))
	z₄ₘ::Array{Float64}(length(xs))

	Rt1::Array{Float64}(length(t))
	Rt2::Array{Float64}(length(t), length(z₃ₘ))
	Rt::Array{Float64}(length(t))
	
	if n == 1
		A= 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end
	
	x₁ₗ = Float64[2l*lx + 4l*zb + xu for l in xs]

    	z₃ₘ = Float64[-2m.*s .- 4m.*zₑ .- zₛ for m in xs]
	z₄ₘ = Float64[-2m.*s .- (4m .- 2).*zₑ .+ zₛ for m in xs]
	
    
    	Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
   	Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))


	Rt2 = @. z₃ₘ'*exp(-(z₃ₘ'^2 / (4D*ν*t))) - z₄ₘ'*exp(-(z₄ₘ'^2 / (4D*ν*t)))

    	Rt = @. Rt1*sum(Rt2, dims=2)
	replace!(Rt, NaN => 0)

	return Rt
end



function DT_transslab(t, β, ρ, s)

    nmed = 1.4
	ndet = 1.4
	n = nmed/ndet

    μa = β[1]
    μsp = β[2]
    # s = slab thickness
	if n == 1
		A = 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end


    D = 1/3μsp
    zₛ = 1/μsp
    zₑ = 2A*D

	ν = 29.9792345/nmed

    m = Vector(-10:1:10)
    z₁ₘ = (1 .- 2m).*s .- 4m.*zₑ .- zₛ
    z₂ₘ = (1 .- 2m).*s .- (4m .- 2).*zₑ .+ zₛ
    


    Rt1 = @. exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
    Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

    
    Rt2 = @. z₁ₘ'*exp(-(z₁ₘ'^2 / (4D*ν*t))) - z₂ₘ'*exp(-(z₂ₘ'^2 / (4D*ν*t)))
    Rt21 = sum(Rt2, dims=2)
    
    

    Rt = @. Rt1*Rt21
    

    replace!(Rt, NaN => 0)
    return Rt
end

A1(n) = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3 + 1.171382n^4 # for n<1
A2(n) = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7 #for n>1

# seems to work for range 0.5<n<1.8 




function DT_refl_paralpip(t, β, rd, rs, L)
	# rd is location of detector [x,y]
	# rs is location of source [x,y]
	# L is dimensions of parallelpiped [Lx, Ly, Lz]



	nmed = 1.4
	ndet = 1.4
	μa = β[1]
	μsp = β[2]
	x = rd[1]
	y = rd[2]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	if (nmed == ndet)
		A = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = 2/(1-R0) - 1
	end

	D = 1/3μsp
	zo = 1/μsp
	zb = 2A*D

	ν = 29.9792345/nmed


	xs = -10:10

	x₁ₗ = Float64[2l*lx + 4l*zb + xu for l in xs]
	x₂ₗ = Float64[2l*lx + (4l-2)*zb - xu for l in xs]
	y₁ₘ = Float64[2m*ly + 4m*zb + yu for m in xs]
	y₂ₘ = Float64[2m*ly + (4m-2)*zb - yu for m in xs]
	z₁ₙ = Float64[2n*lz + 4n*zb + zo for n in xs]
	z₂ₙ = Float64[2n*lz + (4n-2)*zb - zo for n in xs]


	Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

	Rt2 = @. exp(-(x-x₁ₗ')^2 / (4D*ν*t)) - exp(-(x-x₂ₗ')^2 / (4D*ν*t))

	Rt3 = @. exp(-(y-y₁ₘ')^2 / (4D*ν*t)) - exp(-(y-y₂ₘ')^2 / (4D*ν*t))

	Rt4 = @. z₁ₙ'*exp(-z₁ₙ'^2 /(4D*ν*t)) - z₂ₙ'*exp(-z₂ₙ'^2 /(4D*ν*t))

	Rt = Rt1.*sum(Rt2, dims=2).*sum(Rt3, dims=2).*sum(Rt4, dims=2)
	replace!(Rt, NaN =>0)
	return Rt
end




function DT_trans_paralpip(t, β, rd, rs, L)
	# rd is location of detector [x,y]
	# rs is location of source [x,y]
	# L is dimensions of parallelpiped [Lx, Ly, Lz]



	nmed = 1.4
	ndet = 1.4
	μa = β[1]
	μsp = β[2]
	x = rd[1]
	y = rd[2]
	z = rd[3]
	xu = rs[1]
	yu = rs[2]
	lx = L[1]
	ly = L[2]
	lz = L[3]

	if (nmed == ndet)
		A = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = 2/(1-R0) - 1
	end

	D = 1/3μsp
	zo = 1/μsp
	zb = 2A*D

	ν = 29.9792345/nmed


	xs = -10:10

	x₁ₗ = Float64[2l*lx + 4l*zb + xu for l in xs]
	x₂ₗ = Float64[2l*lx + (4l-2)*zb - xu for l in xs]
	y₁ₘ = Float64[2m*ly + 4m*zb + yu for m in xs]
	y₂ₘ = Float64[2m*ly + (4m-2)*zb - yu for m in xs]
	z₁ₙ = Float64[2n*lz + 4n*zb + zo for n in xs]
	z₂ₙ = Float64[2n*lz + (4n-2)*zb - zo for n in xs]


	Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

	Rt2 = @. exp(-(x-x₁ₗ')^2 / (4D*ν*t)) - exp(-(x-x₂ₗ')^2 / (4D*ν*t))

	Rt3 = @. exp(-(y-y₁ₘ')^2 / (4D*ν*t)) - exp(-(y-y₂ₘ')^2 / (4D*ν*t))

	Rt4 = @. (lz - z₁ₙ')*exp(-(lz-z₁ₙ')^2 /(4D*ν*t)) - (lz - z₂ₙ')*exp(-(lz-z₂ₙ')^2 /(4D*ν*t))

	Rt = Rt1.*sum(Rt2, dims=2).*sum(Rt3, dims=2).*sum(Rt4, dims=2)
	replace!(Rt, NaN =>0)
	return Rt
end

