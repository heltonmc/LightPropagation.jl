

outsize = 2 * length(data.IRF) .- 1
nffts = nextfastfft(outsize)

upad = zeros(Float64, nffts)
vpad = zeros(Float64, nffts)
convpad = zeros(Float64, nffts)

copyto!(upad, data.IRF)
pl = plan_rfft(upad)

upadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
upadRfft = pl * upad

RtDT = Array{Float64}(undef, length(data.IRF))


vpadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
u_vpadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))

#copyto!(vpad, RtDT)


function testing1(t, pl, upadRfft, vpad, RtDT, u_vpadRfft, vpadRfft)

    RtDT = DT_model(t, [0.1,10.], 1., 1., 1.)

    copyto!(vpad, RtDT)


    vpadRfft = pl * vpad

    u_vpadRfft = upadRfft .* vpadRfft

    irfft(u_vpadRfft, length(upad))
end





function conv_DT(t, β::Array{Float64,1}, data::input_data, 
                 vpad, RtDT, pl, vpadRfft, upadRfft, u_vpadRfft,
                 convpad
                 )

    RtDT = DT_model(data.t, β, data.ρ, data.nmed, data.ndet) 
    copyto!(vpad, RtDT)

    vpadRfft = pl * vpad
    u_vpadRfft = upadRfft .* vpadRfft
    convpad = irfft(u_vpadRfft, length(upad))
	convpad = convpad./maximum(convpad)

	tidx = findfirst(x -> x == t[1], data.t)

	convpad = convpad[tidx:tidx+length(t)-1]
	return log.(convpad)
end




function getfit(input_data, model_params; alpha=0.1)



    ## get FFT padding and plans
    outsize = 2 * length(data.IRF) .- 1
    nffts = nextfastfft(outsize)

    upad = zeros(Float64, nffts)
    vpad = zeros(Float64, nffts)
    convpad = zeros(Float64, nffts)

    copyto!(upad, data.IRF)
    pl = plan_rfft(upad)

    upadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
    upadRfft = pl * upad

    RtDT = Array{Float64}(undef, length(data.IRF))


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




function test3(input_data, model_params)


## get FFT padding and plans
outsize = 2 * length(input_data.IRF) .- 1
nffts = nextfastfft(outsize)

upad = zeros(Float64, nffts)
vpad = zeros(Float64, nffts)
convpad = zeros(Float64, nffts)

copyto!(upad, input_data.IRF)


end