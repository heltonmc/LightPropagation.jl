@with_kw struct fitDTOF
   
    #required inputs
    t::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    DTOF::Array{Float64,1}
    IRF::Array{Float64,1}

    nmed::Float64 = 1.5 
	ndet::Float64 = 1.51
	lambda::Float64 = 750.0
     
    #semi-infinite and slab
    ρ::Union{Float64, Missing} = missing

    #slab
    s::Union{Float64, Missing} = missing

    #parallelepiped
    rd::Union{Array{Float64,1}, Missing} = missing
    rs::Union{Array{Float64,1}, Missing} = missing
    L::Union{Array{Float64,1}, Missing} = missing
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

function compute_plans(IRFpad)
    pl = plan_rfft(IRFpad)
    fftIRF = pl * IRFpad
    pli = plan_irfft(fftIRF, length(IRFpad))

    return pl, fftIRF, pli
end



function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{TPSF_DA_semiinf_refl};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )


        RtDT = TPSF_DA_semiinf_refl(data.t, β, data.ρ, data.ndet, data.nmed) 
        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end


function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{TPSF_DA_slab_refl};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )

        RtDT = TPSF_DA_slab_refl(data.t, β, data.ρ, data.ndet, data.nmed, data.s) 
        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end


function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{TPSF_DA_paralpip_refl};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )

        RtDT = TPSF_DA_paralpip_refl(data.t, β, data.ndet, data.nmed, data.rd, data.rs, data.L) 
        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end


function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{fluence_DA_semiinf_TD};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )

        RtDT = fluence_DA_semiinf_TD(data.t, β, data.ρ, data.ndet, data.nmed, 0.0)
        RtDT =RtDT./(2*get_afac(data.nmed/data.ndet))

        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end

function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{fluence_DA_semiinf_TD};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )

        RtDT = fluence_DA_semiinf_TD(data.t, β, data.ρ, data.ndet, data.nmed, 0.0)
        RtDT =RtDT./(2*get_afac(data.nmed/data.ndet))

        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end

function _conv_DT(t, β::Array{Float64,1}, data::fitDTOF, out, fftIRF, pl, pli, mod::Val{fluence_DA_semiinf_TD};
    RtDT = Array{Float64}(undef, length(data.IRF))
    )

        RtDT = fluence_DA_semiinf_TD(data.t, β, data.ρ, data.ndet, data.nmed, 0.0)
        RtDT =RtDT./(2*get_afac(data.nmed/data.ndet))

        RtDT = RtDT./maximum(RtDT)
        copyto!(out, RtDT)

        fftRtDT = pl * out
        fftRtDT .*= fftIRF
        out = pli * fftRtDT

        out = out./maximum(out)

	    tidx = findfirst(x -> x == t[1], data.t)

	    out = out[tidx:tidx+length(t)-1]
	return log.(out)
end



function getfit(data::fitDTOF, model_params::DTOF_fitparams)
     ## get FFT padding and plans
     out = pads(data.IRF)

     pl, fftIRF, pli = compute_plans(out)



    #get fit windows
    ind1, ind2 = get_fit_window(data.DTOF, model_params.risefactor, model_params.tailfactor)
    

    fit = curve_fit((t, β) -> _conv_DT(t, β, data, out, fftIRF, pl, pli, Val(model_params.model)), data.t[ind1:ind2], log.(data.DTOF[ind1:ind2]),
        model_params.initparams, lower=model_params.lb, upper=model_params.ub)#; autodiff=:finiteforward)


    return fit, data.t[ind1:ind2], log.(data.DTOF[ind1:ind2]), _conv_DT(data.t[ind1:ind2], fit.param, data, out, fftIRF, pl, pli, Val(model_params.model))

end




using LsqFit, DSP, DelimitedFiles, Plots

struct DTOF_data
    ### Experimental measurements
    # length(t) == length(DTOF) == length(IRF)
    # t should be the same for DTOF and IRF measurements
    t
    DTOF
    IRF
end
Base.@kwdef struct DTOF_fit_params
    model::Function = flux_DA_semiinf_TD       # model to be fitted against
    initparams = [0.2, 10.0]                   # initial starting point of LsqFit
    lb = [0.001, 1.0]                          # lowerbound of fit for [mua, musp]
    ub = [1.0, 20.0]                           # upperbound of fit for [mua, musp]
    risefactor  = 0.9                          # percent of peak to fit on rising edge
    tailfactor  = 0.1                          # percent of peak to fit on falling tail
end

function get_fit_window(counts, risefactor, tailfactor)
    maxval = maximum(counts)
    ind1 = findfirst(x -> x > maxval*risefactor, counts)
    ind2 =  findlast(x -> x >maxval*tailfactor, counts)
    return ind1, ind2
end

function convolve_model(t, β, data, ind1, ind2)

        TPSF = data.model(data.t, β)


        convRt = conv(RtDT, data.IRF)
        convRt = abs.(convRt./maximum(convRt))

    return log10.(convRt[ind1:ind2])
end

function getfit(data, model_params)
    ind1, ind2 = get_fit_window(data.DTOF, model_params.risefactor, model_params.tailfactor)

    fit = curve_fit(
                    (t, β) -> conv_DT(t, β, data, ind1, ind2), 
                    data.t[ind1:ind2],
                    log10.(data.DTOF[ind1:ind2]),
                    model_params.initparams,
                    lower=model_params.lb, 
                    upper=model_params.ub
                    )

    return fit, data.t[ind1:ind2], log10.(data.DTOF[ind1:ind2]), conv_DT(data.t[ind1:ind2], fit.param, data, ind1, ind2)
end

function load_data(IRFfilename, DTOFfilename; ρ = 1.5, n_med = 1.35, n_det = 1.5)
    IRF = readdlm(IRFfilename, skipstart = 10)
    DTOF = readdlm(DTOFfilename, skipstart = 10)

    t = IRF[1:end-1, 1]
    @assert t == DTOF[1:end-1, 1] # ensure time scale of DTOF and IRF are same

    IRF = IRF[1:end-1, 2]
    IRF = IRF ./ maximum(IRF)
    IRF = IRF[100:end]
    IRF = [IRF; fill(eps(Float64),99)] # need to shift time-scale closer to zero because we don't want to use negative t values

    DTOF = DTOF[1:end-1, 2]
    DTOF = DTOF ./ maximum(DTOF)

    return DTOF_fit_inputs(t = t, IRF = IRF, DTOF = DTOF, ρ = ρ, n_med = n_med, n_det = n_det)
end
### code to run 

IRFfilename = "IRF_mirror/TRcurve_10mm_SDS_700nm_10secInt_1.asc"
DTOFfilename = "40mL_IL/TRcurve_15mm_SDS_700nm_30secInt_concen9_1.asc"

data = load_data(IRFfilename, DTOFfilename, ρ = 1.5, n_med = 1.35, n_det = 1.5)

params = DTOF_fit_params(model = flux_DA_semiinf_TD, initparams = [0.06, 12.1, 0.3],
                        risefactor = 0.4,tailfactor = 0.005, lb = [0.0001, 3.0, 0.0], 
                        ub = [0.8, 18.1, 0.5])

result, times, ydata, yfit = getfit(data, params)

plot(times, ydata)
plot!(times, yfit)
result.param






using LsqFit, DSP, DelimitedFiles, Plots

Base.@kwdef struct DTOF_fit_inputs{T}
    ### Experimental measurements from TCSPC
    t::Vector{T}
    DTOF::Vector{T}
    IRF::Vector{T}
end

Base.@kwdef struct DTOF_fitparams
    initparams = [0.2, 10.0, -0.1]             # initial starting point of LsqFit
    lb = [0.001, 1.0, -0.5]                    # lowerbound of fit for [mua, musp, tshift]
    ub = [1.0, 20.0, 0.0]                      # upperbound of fit for [mua, musp, tshift]
    risefactor  = 0.9                          # percent of peak to fit on rising edge
    tailfactor  = 0.1                          # percent of peak to fit on falling tail
end

function load_data(IRFfilename, DTOFfilename)
    IRF = readdlm(IRFfilename, skipstart = 10)
    DTOF = readdlm(DTOFfilename, skipstart = 10)

    t = Float64.(IRF[2:end-1, 1])
    @assert t == DTOF[2:end-1, 1] # ensure time scale of DTOF and IRF are same

    IRF = IRF[2:end-1, 2]
    IRF = IRF ./ maximum(IRF)
    
    DTOF = DTOF[2:end-1, 2]
    DTOF = DTOF ./ maximum(DTOF)

    return DTOF_fit_inputs(t = t, IRF = IRF, DTOF = DTOF)
end

function get_fit_window(counts, risefactor, tailfactor)
    maxval = maximum(counts)
    ind1 = findfirst(x -> x > maxval*risefactor, counts)
    ind2 =  findlast(x -> x >maxval*tailfactor, counts)
    return ind1, ind2
end
function conv_DT(t, β, model, ind1, ind2)
    RtDT = model(t, β)
    convRt = conv(RtDT, data.IRF)
    convRt = abs.(convRt./maximum(convRt))
    return log10.(convRt[ind1:ind2])
end

function getfit(model, data, model_params)
    ind1, ind2 = get_fit_window(data.DTOF, model_params.risefactor, model_params.tailfactor)

    fit = curve_fit((t, β) -> conv_DT(t, β, model, ind1, ind2), data.t[ind1:ind2], log10.(data.DTOF[ind1:ind2]),
                    model_params.initparams, lower=model_params.lb, upper=model_params.ub)

    return fit, data.t[ind1:ind2], log10.(data.DTOF[ind1:ind2]), conv_DT(data.t[ind1:ind2], model, fit.param, ind1, ind2)
end

### code to run 

IRFfilename = "IRF_mirror/TRcurve_10mm_SDS_700nm_10secInt_1.asc"
DTOFfilename = "40mL_IL/TRcurve_15mm_SDS_700nm_30secInt_concen9_1.asc"

data = load_data(IRFfilename, DTOFfilename)

params = DTOF_fitparams(
    initparams = [0.06, 12.1, 0.3],
    risefactor = 0.4, tailfactor = 0.005, 
    lb = [0.0001, 3.0, 0.0], 
    ub = [0.8, 18.1, 0.5]
)


model(t, β) = flux_DA_semiinf_TD(t, 1.5, β[1], β[2], n_ext = 1.4, n_med = 1.3)



result, times, ydata, yfit = getfit(model, data, params)

plot(times, ydata)
plot!(times, yfit)
result.param