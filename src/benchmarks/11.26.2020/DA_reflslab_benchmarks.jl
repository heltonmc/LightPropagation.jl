#=
DA_reflslab is the first iteration to compute the TD reflectance from a slab. 
This code is very 'MATLABY'... as in you generally want to avoid loops in MATLAB and vectorize everything.
Though it is slow, I find this very composable and readable (I feel there should always be a trade off)
With that said the more optimized version is ~2.5x faster with the ability to scale much better.
The main performance increase came from completely devectorizing the code and iteratively looping through 
the summation and time vector.
Some improvement also came by imposing the @inbounds macro which skips the checking bounds. 
The primary advantage of creating loops and devectorizing is it allows multi-threading. 
The Threads.@threads runs the for loops on multiple threads on your cpu. This ran on 4 threads but as mentioned
	this operation is easily scalable with more threads in CPU. 
*** Julia's default is 1 thread where MATLAB I believe defaults to largest available. You have to initiate Julia
running multiple threads....
=#

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


function DA_reflslab1(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness

	  n::Float64 = nmed/ndet
      μa::Float64 = β[1]
      μsp::Float64 = β[2]
      D::Float64 = 1/3μsp
      ν::Float64 = 29.9792345/nmed
      xs::UnitRange{Int64} = -10:10



      Rt1 = Array{Float64}(undef, length(t))
	  Rt2 = zeros(Float64, length(t))

		@inbounds begin

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

		    Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))
			
				Threads.@threads for n in eachindex(t) 
					for m in xs

						z3m = -2m.*s .- 4m.*ze .- zs
						z4m = -2m.*s .- (4m .- 2).*ze .+ zs

			 			Rt2[n] += z3m*exp(-(z3m^2 / (4D*ν*t[n]))) - z4m*exp(-(z4m^2 / (4D*ν*t[n])))
					end
					
					
				end
			end



	return replace!(Rt1.*Rt2, NaN => 0)

end
