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
    out = zeros(eltype(u), nextfastfft(outsize))
    copyto!(out, u)

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
    out = pads(input_data.IRF)

    pl = plan_rfft(out)
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


using DelimitedFiles, LsqFit, DSP, FFTW
filename1 = "/home/heltonmc/Desktop/DTOF.asc"
filename = "/home/heltonmc/Desktop/IRF.asc"

  data_array = readdlm(filename,skipstart=10)
data_array1 = readdlm(filename1,skipstart=10)

RtIRF = data_array[1:end-1,2]
RtIRF = RtIRF./maximum(RtIRF)
tIRF = data_array[1:end-1,1]

RtDTOF = data_array1[1:end-1,2]
RtDTOF = RtDTOF./maximum(RtDTOF)
tDTOF = data_array1[1:end-1,1]


struct input_data
    t::Array{Float64,1}
	DTOF::Array{Float64,1}
	IRF::Array{Float64,1}
	DTOFerrors::Array{Float64,1}
	ρ::Float64
	nmed::Float64
	ndet::Float64
	lambda::Float64
	normpeaks::Bool
end

struct model_params
	#yerrors::Array{Float64,1}          #uncertainties in dependent variable
    model::Function                    #model to be fitted against
	initparams::Array{Float64,1}
	lb::Array{Float64,1}
	ub::Array{Float64,1}
	risefactor::Float64  			#percent of peak to fit on rising edge
	tailfactor::Float64				#perfecnt of peak to fit on failing tail
	dof::Int64
end


struct fitresult
    #Fit inputs
    xdata::Array{Float64,1}     #independent variable
    ydata::Array{Float64,1}     #dependent variable
    yerrors::Array{Float64,1}   #uncertainties in dependent variable
	model::Function             #model to be fitted against

    #Fit results
    residuals::Array{Float64,1}  	 #weighted residuals
    params::Array{Float64,1}    		#best fit parameters for model
	marginerror::Array{Float64,1}        #product of standard error and the critical value of each parameter at a certain significance level (default is 5%) from t-distribution.
    perrors::Array{Float64,1}            #uncertainties in best fit parameters
    dof::Int64                  		#number of degrees of freedom in fit
    chisq::Float64             		 #chi^2 value of fit
    probability::Float64      		  #probability of chi^2 value
end




data = input_data(tDTOF, RtDTOF, RtIRF, rand(length(RtDTOF)), 2.56, 1.5, 1.51, 700, false)

modelp = model_params(DT_model, [0.2,10], [0.001,1], [1,20], 0.9, 0.1, 2)
