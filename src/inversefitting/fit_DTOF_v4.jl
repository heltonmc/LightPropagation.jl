sing DSP, DelimitedFiles, LsqFit

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
    out = zeros(eltype(u), nextfastfft(outsize))
end

function get_fftplan(u::Array{Float64,1}, out)
    upad = zeros(eltype(u), length(out))
    copyto!(upad,u)
    pl = plan_rfft(upad)
    upadRfft = pl * upad
    return pl, upadRfft
end


function conv_DT(t, β::Array{Float64,1}, data::input_data, out, fftIRF, pl, pli;
    RtDT = Array{Float64}(undef, length(data.IRF))
    )


        RtDT = DT_model(data.t, β, data.ρ, data.nmed, data.ndet) 
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end



function getfit(input_data, model_params)
    ## get FFT padding and plans
    out::Array{eltype(input_data.IRF)} = pads(input_data.IRF)
    copyto!(out, input_data.IRF)

    pl = plan_rfft(out; flags = FFTW.MEASURE)
    fftIRF = pl * out
    pli = plan_irfft(fftIRF, length(out))



    #get fit windows
    ind1, ind2 = get_fit_window(input_data.DTOF, model_params.risefactor, model_params.tailfactor)


    curve_fit((t, β) -> conv_DT(t, β, input_data, out, fftIRF, pl, pli), input_data.t[ind1:ind2], log.(input_data.DTOF[ind1:ind2]),
    model_params.initparams, lower=model_params.lb, upper=model_params.ub)

end

#=
@benchmark getfit(data,modelp)
BenchmarkTools.Trial: 
  memory estimate:  21.43 MiB
  allocs estimate:  5717
  --------------
  minimum time:     8.924 ms (0.00% GC)
  median time:      9.845 ms (0.00% GC)
  mean time:        12.204 ms (4.42% GC)
  maximum time:     54.299 ms (3.45% GC)
  --------------
  samples:          410
  evals/sample:     1

this is if you use only 1024 points
    data1 = input_data(tDTOF[1:1024], RtDTOF[500:1524], RtIRF[500:1524], rand(1024), 2.56, 1.5, 1.51, 700, false)

julia> @benchmark getfit(data1,modelp)
BenchmarkTools.Trial: 
  memory estimate:  6.13 MiB
  allocs estimate:  5501
  --------------
  minimum time:     3.154 ms (0.00% GC)
  median time:      3.261 ms (0.00% GC)
  mean time:        5.035 ms (3.86% GC)
  maximum time:     88.410 ms (0.00% GC)
  --------------
  samples:          993
  evals/sample:     1
