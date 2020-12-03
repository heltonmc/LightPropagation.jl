

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

    ind1, ind2 = get_fit_window(input_data.DTOF, model_params.risefactor, model_params.tailfactor)

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


end[[[[[[[[[]]]]]]]]]


function get_fit_window(counts, risefactor, tailfactor)

    maxval = maximum(counts)
    ind1 = findfirst(x -> x > maxval*risefactor, counts)
    ind2 =  findlast(x -> x >maxval*tailfactor, counts)
    
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



out = pads(u)
pl, upadRfft = get_fftplan(u, out)




function conv_DT(t, β::Array{Float64,1},
    data::input_data, 
    out, 
    pl,
    upadRfft;
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





function getfit(input_data, model_params; alpha=0.1)
    ## get FFT padding and plans
    out = pads(input_data.IRF)
    pl, upadRfft = get_fftplan(input_data.IRF, out)

    #get fit windows
    ind1, ind2 = get_fit_window(input_data.DTOF, model_params.risefactor, model_params.tailfactor)


    fit = curve_fit((t, β) -> conv_DT(t, β, input_data, out, pl, upadRfft), input_data.t[ind1:ind2], log.(input_data.DTOF[ind1:ind2]),
     model_params.initparams, lower=model_params.lb, upper=model_params.ub)

end




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
