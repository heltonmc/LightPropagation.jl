
using DSP, BenchmarkTools

function model(t, β::Array{Float64,1}, ρ::Float64)

    Rt1 = Array{Float64}(undef, length(t))
    Rt2 = Array{Float64}(undef, length(t))
	A::Float64 = - 1/β[2]
    B::Float64 = 2/β[2] + 1/β[2]

	Threads.@threads for n in eachindex(t)
		
   		Rt1[n] = -exp(-(ρ^2/(t[n])) - β[1]*t[n])
   		Rt1[n] = Rt1[n]/(t[n]^(5/2))

   		Rt2[n] = A*exp(-(A^2/(t[n]))) - B*exp(-(B^2/(t[n])))
		
	end
	
    Rt = replace!(Rt1.*Rt2, NaN => 0)
    return Rt./maximum(Rt)
 
 end
 

 function conv_Rt(t, β::Array{Float64,1},  ρ::Float64, It)

	Rt = Array{Float64}(undef, length(t))
	convRt = Array{Float64}(undef, 2*length(t) - 1)

	Rt = model(t, β, ρ)

	convRT = conv(It, Rt)
    convRt = convRt./maximum(convRt)
    
end






function _conv_kern_fft!(out, u, v, su, sv, outsize, nffts)
    upad = _zeropad(u, nffts)
    vpad = _zeropad(v, nffts)
    p! = plan_fft!(upad)
    p! * upad # Operates in place on upad
    p! * vpad
    upad .*= vpad
    ifft!(upad)
    copyto!(out,
            CartesianIndices(out),
            upad,
            CartesianIndices(UnitRange.(1, outsize)))
end


@inline function _zeropad!(
    padded::AbstractVector,
    u::AbstractVector,
    padded_axes = axes(padded),
    data_dest::Tuple = (first(padded_axes[1]),),
    data_region = CartesianIndices(u),
)
    datasize = length(data_region)
    # Use axes to accommodate arrays that do not start at index 1
    data_first_i = first(data_region)[1]
    dest_first_i = data_dest[1]
    copyto!(padded, dest_first_i, u, data_first_i, datasize)
    padded[first(padded_axes[1]):dest_first_i - 1] .= 0
    padded[dest_first_i + datasize : end] .= 0

    padded
end

@inline function _zeropad!(
    padded::AbstractArray,
    u::AbstractArray,
    padded_axes = axes(padded),
    data_dest::Tuple = first.(padded_axes),
    data_region = CartesianIndices(u),
)
    fill!(padded, zero(eltype(padded)))
    dest_axes = UnitRange.(data_dest, data_dest .+ size(data_region) .- 1)
    dest_region = CartesianIndices(dest_axes)
    copyto!(padded, dest_region, u, data_region)

    padded
end




IRFpad = zeros(2*length(IRF) -1)
copyto!(IRFpad, IRF)

p! = plan_fft!(IRFpad)
p! * IRFpad




su = size(u)
sv = size(v)
outsize = su .+ sv .-1
axesu = axes(u)
axesv = axes(v)
out_offsets = first.(axesu) .+ first.(axesv)
out_axes = UnitRange.(out_offsets, out_offsets .+outsize .- 1)
out = similar(u, out_axes)

nffts = nextfastfft(outsize)

padded = similar(u, nffts)
padded_axes = axes(padded)
data_region = CartesianIndices(u)

fill!(padded, zero(eltype(padded)))
data_dest = first.(padded_axes)
dest_axes = UnitRange.(data_dest, data_dest .+ size(data_region) .- 1)
dest_region = CartesianIndices(dest_axes)
copyto!(padded, dest_region, u, data_region)

p! = plan_fft!(convert(Array{Complex{Float64}}, upad))
p! * upad
p! * vpad
upad .*= vpad
ifft!(upad)



function pad_fft(u, v)


    su = size(u)
    sv = size(v)
    outsize = su .+ sv .-1
    nffts = nextfastfft(outsize)


    upad = similar(u, nffts)
    vpad = similar(v, nffts)
    fill!(upad, zero(eltype(upad)))
    fill!(vpad, zero(eltype(vpad)))

    copyto!(upad, u)
    copyto!(vpad, v)

    p! = plan_fft!(convert(Array{Complex{Float64}}, upad))
    upad = p! * upad



function testing(upad, vpad, p!)

    
    vpad = p! * vpad
    vpad .*= upad
    real(ifft!(vpad))
end




outsize = 2 * length(IRF) .- 1
nffts = nextfastfft(outsize)

upad = zeros(Float64, nffts)
vpad = zeros(Float64, nffts)
convpad = zeros(Float64, nffts)

copyto!(upad, IRF)
copyto!(vpad, RtDT)

pl = plan_rfft(upad)

function testing1(pl, nffts, upad, vpad)
upadRfft = zeros(Complex{Float64}, 1 .+ div(nffts, 2, RoundDown))
vpadRfft = zeros(Complex{Float64}, 1 .+ div(nffts, 2, RoundDown))
u_vpadRfft = zeros(Complex{Float64}, 1 .+ div(nffts, 2, RoundDown))

upadRfft = pl * upad
vpadRfft = pl * vpad

u_vpadRfft = upadRfft .* vpadRfft

irfft(u_vpadRfft, length(upad))
end








outsize = 2 * length(IRF) .- 1
nffts = nextfastfft(outsize)

upad = zeros(Float64, nffts)
vpad = zeros(Float64, nffts)
convpad = zeros(Float64, nffts)

copyto!(upad, IRF)
pl = plan_rfft(upad)

upadRfft = zeros(Complex{Float64}, 1 .+ div(length(upad), 2, RoundDown))
upadRfft = pl * upad

RtDT = Array{Float64}(undef, length(IRF))


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



function conv_DT(t, IRF)

	RtDT = Array{Float64}(undef, length(t))
	convDT = Array{Float64}(undef, 2*length(t) - 1)

    RtDT = DT_model(t, [0.1,10.], 1., 1., 1.)

	convDT = conv(IRF, RtDT)
	convDT = convDT./maximum(convDT)
end
