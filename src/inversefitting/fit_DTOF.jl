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



function DT_model(t,β,ρ,nmed,ndet)
	#diffusion theory model for homogenous semi-infinite layer
	μa = β[1]
	μsp = β[2]


	if (nmed == ndet)
		Afac = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = 2/(1-R0) - 1
	end

	z0 = 1/(μa+μsp)
	D = 1/(3*(μa + μsp))
	zb = 2*Afac*D

	v = 29.9792345/nmed # speed of light in medium

	Rt1 = @. v*exp(-μa*v*t)
	Rt1 = @. Rt1/((4*π*D*v*t)^1.5)
	Rt1 = @. Rt1*(exp(-(z0^2 + ρ^2)/(4*D*v*t)) - exp(-((2*zb + z0)^2 + ρ^2)/(4*D*v*t)))

	Rt2 = @. 3*exp(-μa*v*t)/(2*((4*π*D*v)^1.5)*(t^2.5))
	Rt2 = @. Rt2*(z0*exp(-((z0^2 + ρ^2)/(4*D*v*t))) + (2*zb + z0)*exp(-((2*zb+z0)^2 + ρ^2)/(4*D*v*t)))

	Rt = @. abs(Rt1)+abs(Rt2)
	replace!(Rt, NaN => 0)
	m = findmax(Rt)
	Rt = Rt./m[1]
end


function conv_DT(t, β, data)

	RtDT = DT_model(data.t, β, data.ρ, data.nmed, data.ndet)

	convDT = conv(data.IRF, RtDT)
	convDT = convDT./maximum(convDT)

	tidx = findfirst(x -> x == t[1], data.t)

	if data.normpeaks == true
		peakindDTOF, peakindconvDT = findmax(data.DTOF)[2], findmax(convDT)[2]
		if peakindDTOF > peakindconvDT
			convDT = [zeros(peakindDTOF-peakindconvDT); convDT]
		elseif peakindDTOF < peakindconvDT
			convDT = convDT[peakindconvDT-peakindDTOF:end]
		end
	end


	convDT = convDT[tidx:tidx+length(t)-1]
	return log.(convDT)
end




function getfit(input_data, model_params; alpha=0.1)

	#get fitting window
	counts = input_data.DTOF
	maxvalue, maxindex = findmax(counts)

	ind1 = findfirst(x -> x > maxvalue*model_params.risefactor, counts)
	ind2 =  findlast(x -> x >maxvalue*model_params.tailfactor, counts)

	weights = 1 ./ input_data.DTOFerrors

	fit = curve_fit((t, β) -> conv_DT(t, β, input_data), input_data.t[ind1:ind2], log.(input_data.DTOF[ind1:ind2]), model_params.initparams, lower=model_params.lb, upper=model_params.ub)

	perrors = estimate_errors(fit,alpha)
	marginerror = margin_error(fit,alpha)

	chisq = sum(fit.resid.^2)
	chisqdist = Distributions.Chisq(2)
	probability = Distributions.ccdf(chisqdist,chisq)
	return fitresult(input_data.t[ind1:ind2],log.(input_data.DTOF[ind1:ind2]),input_data.DTOFerrors, DT_model, fit.resid,fit.param,marginerror,perrors,2,chisq,probability)
end


function plotfit(fit::fitresult1, data::input_data)
	scatter(fit.xdata, fit.ydata, color="black", label="DTOF", markersize=3, alpha=0.8)

	#scatter(data.t[600:900], log.(data.DTOF[600:900]), color="black", label="DTOF", markersize=3, alpha=0.8)
	#plot!(data.t[600:900], conv_DT(data.t[600:900], fit.params, data), color="red", lw=3, label="DT-fit", alpha=0.8)

	plot!(fit.xdata, conv_DT(fit.xdata, fit.params, data), color="red", lw=3, label="DT-fit", alpha=0.8)
	xlabel!("time (ns)")
	ylabel!("Counts")

end

function plotfit_n(fit::fitresult1, data::input_data)

	xind1 = findmax(fit.ydata)[2]
	xind2 = findmax(conv_DT(fit.xdata, fit.params, data))[2]
	scatter(data.t[600:900].- fit.xdata[xind1], log.(data.DTOF[600:900]), color="black", label="DTOF", markersize=3, alpha=0.8)
	plot!(data.t[600:900].- fit.xdata[xind2], conv_DT(data.t[600:900], fit.params, data), color="red", lw=3, label="DT-fit", alpha=0.8)

	#scatter(fit.xdata .- fit.xdata[xind1], fit.ydata, color="black", label="DTOF", markersize=3, alpha=0.8)
	#plot!(fit.xdata .- fit.xdata[xind2], conv_DT(fit.xdata, fit.params, data), color="red", lw=3, label="DT-fit", alpha=0.8)
	xlabel!("time (ns)")
	ylabel!("Counts")
end



function loadData(filename)
    #loads data from asc
    data_array = readdlm(filename,skipstart=10)
    t = data_array[1:end-1,1]
    Rt = data_array[1:end-1,2]
    (time = convert(Array{Float64},t), counts = convert(Array{Float64},Rt))
end

