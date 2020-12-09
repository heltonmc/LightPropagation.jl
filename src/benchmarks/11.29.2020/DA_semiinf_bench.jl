# Benchmark tests comparing the semiinf solution on 1 and four threads.
# DA_semiinf has almost no time dependence on the # of threads used 
# DA_semiinf1 has been edited for custom loops that can utilize multi-threading
# This results in DA_semiinf1 being 1.6x faster 
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
#=
1 THREAD
BenchmarkTools.Trial: 
  memory estimate:  48.16 KiB
  allocs estimate:  9
  --------------
  minimum time:     78.548 μs (0.00% GC)
  median time:      81.226 μs (0.00% GC)
  mean time:        89.842 μs (5.96% GC)
  maximum time:     6.547 ms (98.55% GC)
  --------------
  samples:          10000
  evals/sample:     1

4 THREADS
BenchmarkTools.Trial: 
  memory estimate:  48.16 KiB
  allocs estimate:  9
  --------------
  minimum time:     78.612 μs (0.00% GC)
  median time:      81.431 μs (0.00% GC)
  mean time:        87.015 μs (5.90% GC)
  maximum time:     6.618 ms (98.17% GC)
  --------------
  samples:          10000
  evals/sample:     1

12 THREADS (Michaels Desktop CPU: Intel 8700k)
BenchmarkTools.Trial: 
  memory estimate:  48.16 KiB
  allocs estimate:  9
  --------------
  minimum time:     46.689 μs (0.00% GC)
  median time:      48.945 μs (0.00% GC)
  mean time:        50.164 μs (0.79% GC)
  maximum time:     613.391 μs (84.43% GC)
  --------------
  samples:          10000
  evals/sample:     1

=#

function DA_semiinf1(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)


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

	Threads.@threads for n in eachindex(t)
		
   		Rt1[n] = -exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])
   		Rt1[n] = Rt1[n]/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

   		Rt2[n] = z3m*exp(-(z3m^2/(4D*ν*t[n]))) - z4m*exp(-(z4m^2/(4D*ν*t[n])))
		
	end
	
    

    replace!(Rt1.*Rt2, NaN => 0)

end



#=

1 THREAD
BenchmarkTools.Trial: 
  memory estimate:  24.16 KiB
  allocs estimate:  6
  --------------
  minimum time:     62.370 μs (0.00% GC)
  median time:      64.190 μs (0.00% GC)
  mean time:        68.621 μs (4.76% GC)
  maximum time:     6.965 ms (98.77% GC)
  --------------
  samples:          10000
  evals/sample:     1


4 THREADS
BenchmarkTools.Trial: 
  memory estimate:  27.05 KiB
  allocs estimate:  27
  --------------
  minimum time:     43.496 μs (0.00% GC)
  median time:      51.549 μs (0.00% GC)
  mean time:        98.719 μs (2.47% GC)
  maximum time:     8.069 ms (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1

  12 THREADS : (Michael's Desktop)


=#





function DT_model3(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)

    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
  	ν::Float64 = 29.9792345/nmed

    Rt = Array{Float64}(undef, length(t))


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

    Threads.@threads for n in eachindex(t)
		
        Rt[n] = -exp(-(ρ^2/(4D*ν*t[n])) - μa*ν*t[n])
        Rt[n] = Rt[n]/(2*(4π*D*ν)^(3/2)*t[n]^(5/2))

        Rt[n] = Rt[n]*(z3m*exp(-(z3m^2/(4D*ν*t[n]))) - z4m*exp(-(z4m^2/(4D*ν*t[n]))))

        if Rt[n] == NaN
            Rt[n] = 0
        end
     
    end
 
     return Rt

end

#= 12 threads
julia> @benchmark DT_model3(0:0.01:10, [0.1,10.], 1.0, 1.0, 1.0)
BenchmarkTools.Trial: 
  memory estimate:  16.95 KiB
  allocs estimate:  62
  --------------
  minimum time:     9.657 μs (0.00% GC)
  median time:      11.596 μs (0.00% GC)
  mean time:        13.072 μs (3.51% GC)
  maximum time:     1.661 ms (96.79% GC)
  --------------
  samples:          10000
  evals/sample:     1
