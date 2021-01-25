

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