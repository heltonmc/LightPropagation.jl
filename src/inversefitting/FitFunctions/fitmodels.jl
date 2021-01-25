function get_afac(n::Float64)
    if n > 1.0
        A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
        5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
    elseif n < 1.0
        A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
    else 
        A = 1.0
    end
    return A
end

### Semi-infinite models ###

function _get_fit_model!(Rt::Array{Float64}, β::Array{Float64,1}, data::fitDTOF, mod::Val{refl_DA_semiinf_TD})
    Rt = refl_DA_semiinf_TD(data.t, β, data.ndet, data.nmed)
end

# fit DTOF to fluence with EBPC conditions
function _get_fit_model!(Rt::Array{Float64}, β::Array{Float64,1}, data::fitDTOF, mod::Val{fluence_DA_semiinf_TD})
    Rt = fluence_DA_semiinf_TD(data.t, β, data.ndet, data.nmed, 0.0)
    return Rt./(2*get_afac(data.nmed/data.ndet))
end

### Slab models ###

#=
function _get_fit_model!(Rt::Array{Float64}, β::Array{Float64,1}, data::fitDTOF, mod::Val{refl_DA_slab_TD})
    Rt = refl_DA_slab_TD(data.t, β, data.ndet, data.nmed, data.s)
end


### Parallelepiped models ###

function _get_fit_model!(Rt::Array{Float64}, β::Array{Float64,1}, data::fitDTOF, mod::Val{refl_DA_paralpip_TD})
    Rt = refl_DA_paralpip_TD(data.t, β, data.ndet, data.nmed, data.rd, data.rs, data.L)
end

