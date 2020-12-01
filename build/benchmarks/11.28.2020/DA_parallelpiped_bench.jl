#Similar changes to the parallelpiped as were made to the slab on 11.26.2020. This results in new version being >2.9x faster

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

#=
@benchmark DT_refl_paralpip(0:0.01:10, [0.1,10.0], 1.0, 1.0, [2.5,5.0], [5.0,5.0], [10.0,10.0,10.0])
BenchmarkTools.Trial: 
  memory estimate:  1.01 MiB
  allocs estimate:  46
  --------------
  minimum time:     1.903 ms (0.00% GC)
  median time:      2.141 ms (0.00% GC)
  mean time:        2.227 ms (2.61% GC)
  maximum time:     5.770 ms (0.00% GC)
  --------------
  samples:          2245
  evals/sample:     1
  =#

  function DT_refl_paralpip1(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})


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


		Rt2 = zeros(Float64, length(t))
		Rt3 = zeros(Float64, length(t))
		Rt4 = zeros(Float64, length(t))

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

		Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

	
		Threads.@threads for n in eachindex(t) 
			for m in xs
				x1l = 2m.*lx .+ 4m.*zb .+ xu 
		  		x2l = 2m.*lx + (4m-2).*zb .- xu
		  		y1m = 2m.*ly .+ 4m.*zb .+ yu
		  		y2m = 2m.*ly .+ (4m.-2).*zb .- yu
		  		z1n = 2m.*lz .+ 4m.*zb .+ zo
		  		z2n = 2m.*lz .+ (4m.-2).*zb .- zo
			
  				Rt2[n] += exp(-(x-x1l)^2 / (4D*ν*t[n])) - exp(-(x-x2l)^2 / (4D*ν*t[n]))

  				Rt3[n] += exp(-(y-y1m)^2 / (4D*ν*t[n])) - exp(-(y-y2m)^2 / (4D*ν*t[n]))

  				Rt4[n] += z1n*exp(-z1n^2 /(4D*ν*t[n])) - z2n'*exp(-z2n^2 /(4D*ν*t[n]))

			end
		end
  	return replace!(Rt1.*Rt2.*Rt3.*Rt4, NaN => 0)
  end
#=
  @benchmark DT_refl_paralpip1(0:0.01:10, [0.1,10.0], 1.0, 1.0, [2.5,5.0], [5.0,5.0], [10.0,10.0,10.0])
  BenchmarkTools.Trial: 
  memory estimate:  43.34 KiB
  allocs estimate:  30
  --------------
  minimum time:     681.021 μs (0.00% GC)
  median time:      733.641 μs (0.00% GC)
  mean time:        990.380 μs (0.41% GC)
  maximum time:     9.611 ms (49.77% GC)
  --------------
  samples:          5046
  evals/sample:     1
=#