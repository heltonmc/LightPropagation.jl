@with_kw struct fitDTOF
   
    #required inputs
    t::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    DTOF::Array{Float64,1}
    IRF::Array{Float64,1}
     
    #semi-infinite and slab
    ρ::Union{Float64, Missing} = missing

    #slab
    s::Union{Float64, Missing} = missing

    #parallelepiped
    rd::Union{Array{Float64,1}, Missing} = missing
    rs::Union{Array{Float64,1}, Missing} = missing
end



@with_kw struct fitg2
   
    #required inputs
    τ::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    g2::Array{Float64,1}
     
    #semi-infinite
    ρ::Union{Float64, Missing}
end
