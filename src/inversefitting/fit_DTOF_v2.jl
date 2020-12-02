using LsqFit, Distributions, DSP, DelimitedFiles

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



function DT_model(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)
	#optimzed semi-infinite model

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
	
    

   Rt = replace!(Rt1.*Rt2, NaN => 0)
   return Rt./maximum(Rt)

end


function conv_DT(
    t, 
    β::Array{Float64,1},
    input_data, 
    vpad, RtDT,
    pl,
    vpadRfft,
    upadRfft,
    u_vpadRfft,
    convpad
)

    RtDT = DT_model(input_data.t, β, input_data.ρ, input_data.nmed, input_data.ndet) 
    copyto!(vpad, RtDT)

    vpadRfft = pl * vpad
    u_vpadRfft = upadRfft .* vpadRfft
    convpad = irfft(u_vpadRfft, length(upad))
    convpad = convpad./maximum(convpad)

    tidx = findfirst(x -> x == t[1], data.t)

    convpad = convpad[tidx:tidx+length(t)-1]
    return log.(convpad)
end




function getfit3(input_data, model_params; alpha=0.1)



## get FFT padding and plans
    outsize = 2 * length(input_data.IRF) .- 1
    nffts = nextfastfft(outsize)

    upad = zeros(Float64, nffts)
    vpad = zeros(Float64, nffts)
    convpad = zeros(Float64, nffts)

    copyto!(upad, input_data.IRF)
    pl = plan_rfft(upad)

    upadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
    upadRfft = pl * upad

    RtDT = Array{Float64}(undef, length(input_data.IRF))


    vpadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
    u_vpadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))

#get fitting window
    counts = input_data.DTOF
    maxvalue, maxindex = findmax(counts)

    ind1 = findfirst(x -> x > maxvalue*model_params.risefactor, counts)
    ind2 =  findlast(x -> x >maxvalue*model_params.tailfactor, counts)

    weights = 1 ./ input_data.DTOFerrors

    fit = curve_fit((t, β) -> conv_DT(t, β, input_data, vpad, RtDT, pl, vpadRfft, upadRfft, u_vpadRfft, convpad), input_data.t[ind1:ind2], log.(input_data.DTOF[ind1:ind2]),
    model_params.initparams, lower=model_params.lb, upper=model_params.ub)

    perrors = estimate_errors(fit,alpha)
    marginerror = margin_error(fit,alpha)

    chisq = sum(fit.resid.^2)
    chisqdist = Distributions.Chisq(2)
    probability = Distributions.ccdf(chisqdist,chisq)
    return fitresult(input_data.t[ind1:ind2],log.(input_data.DTOF[ind1:ind2]),input_data.DTOFerrors, DT_model, fit.resid,fit.param,marginerror,perrors,2,chisq,probability)
end