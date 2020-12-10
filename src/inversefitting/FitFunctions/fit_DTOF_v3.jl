using DSP, DelimitedFiles, LsqFit

function DT_model(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)
	#optimzed semi-infinite model

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
		
	end
	
    

   replace!(Rt, NaN => 0)
   return Rt./maximum(Rt)

end


function get_fit_window(counts, risefactor, tailfactor)

    maxval = maximum(counts)
    ind1::Int = findfirst(x -> x > maxval*risefactor, counts)
    ind2::Int =  findlast(x -> x >maxval*tailfactor, counts)
    
    return ind1, ind2
end


function pads(u::Array{Float64,1})
    outsize = 2 .* length(u) .- 1
    nffts = nextfastfft(outsize)
    out = zeros(eltype(u), nffts)
end

function get_fftplan(u::Array{Float64,1}, out)
    upad = zeros(eltype(u), length(out))
    copyto!(upad,u)
    pl = plan_rfft(upad)
    upadRfft = pl * upad
    return pl, upadRfft
end


function conv_DT(t, β::Array{Float64,1}, data::input_data, out, pl,upadRfft;
    vpad = zeros(eltype(data.DTOF), length(out)), 
    RtDT = Array{Float64}(undef, length(data.IRF)),
    vpadRfft = zeros(eltype(upadRfft), length(upadRfft))
    )


        RtDT = DT_model(data.t, β, data.ρ, data.nmed, data.ndet) 
        copyto!(vpad, RtDT)

        vpadRfft = pl * vpad
        vpadRfft .*= upadRfft
        out = irfft(vpadRfft, length(out))

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end



function getfit(input_data, model_params)
    ## get FFT padding and plans
    out::Array{eltype(input_data.IRF)} = pads(input_data.IRF)
    pl, upadRfft = get_fftplan(input_data.IRF, out)

    #get fit windows
    ind1, ind2 = get_fit_window(input_data.DTOF, model_params.risefactor, model_params.tailfactor)


    curve_fit((t, β) -> conv_DT(t, β, input_data, out, pl, upadRfft), input_data.t[ind1:ind2], log.(input_data.DTOF[ind1:ind2]),
    model_params.initparams, lower=model_params.lb, upper=model_params.ub)

end

#= 12/3/2020
@benchmark getfit(data,modelp)
BenchmarkTools.Trial: 
  memory estimate:  30.68 MiB
  allocs estimate:  9561
  --------------
  minimum time:     15.015 ms (0.00% GC)
  median time:      16.626 ms (5.95% GC)
  mean time:        18.319 ms (3.91% GC)
  maximum time:     64.462 ms (1.59% GC)
  --------------
  samples:          273
  evals/sample:     1
=#