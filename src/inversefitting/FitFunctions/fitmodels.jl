

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
=#
