#=
Here are two approaches to deal with convolutions in the inverse fitting approach. 
Taking the convolution of the IRF and the model appears to be computationally expensive.
Each iteration in the least squares fit has to compare the new fit with the data after convolving the model w/IRF.

Here I am trying to take advantage of the fact that the fft(IRF) is fixed for each step in the iteration and there is no reason to repeatedly calculate it.
The challenge is that the conv function from the FFTW library is very optimzed as it is a wrapper of a C/C++ library that is highly optimized.
Even still, there is an opportunity to plan out the FFT beforehand using the fact that the u and v in conv(u, v) are of the same length for our problem. 
This allows to calculate the linear operator by the FFT needed to compute fft for a given set of dimensions.
In other words we can calculate a plan_fft (call this p) operator and then calculate the fft by u*p and v*p. If u is the IRF this value is constant for each iteration in the least squares.

=#


model(t, p) = p[1] * exp.(-p[2] *t)

tdata = 1:0.01:10

ydata = model(tdata, [1.0, 2.0]) + randn(length(tdata))

IRF = [Vector(1:500); Vector(500:-1:100)]

ydata = conv(ydata, IRF)


function convFIT(tdata, p, IRF)

    Rt = model(tdata, p)

    Rt = conv(Rt, IRF)
end



function fit(tdata, ydata, IRF)
    curve_fit((tdata, p) -> convFIT(tdata, p, IRF), tdata, ydata, [0.5, 0.5])
end


function convFIT1(tdata, p, pl, IRF_fft)
    Rt = model(tdata, p)
    dims = 2*length(Rt) -1

    Rt_fft = convert(Array{Complex{Float64},1}, [Rt ; zeros(Float64, dims-length(Rt))])
    r = similar(Rt_fft)
    i = similar(Rt_fft)
    c = similar(Rt_fft)


    mul!(r, pl, Rt_fft)
    mul!(i, pl, IRF_fft)

   ldiv!(c, pl, (r.*i))
   return real(c)
end





function fit2(tdata, ydata, IRF)

    dims = 2*length(IRF) -1

    IRF_fft = convert(Array{Complex{Float64},1}, [IRF ; zeros(Float64, dims-length(IRF))])

    pl = plan_fft(IRF_fft)


    curve_fit((tdata, p) -> convFIT1(tdata, p, pl, IRF_fft), tdata, ydata, [0.5, 0.5])
end

#=

############# Run these commands after running previous lines ################
julia> @benchmark f = fit(tdata,ydata,IRF)
BenchmarkTools.Trial: 
  memory estimate:  132.05 MiB
  allocs estimate:  837915
  --------------
  minimum time:     792.865 ms (0.64% GC)
  median time:      801.574 ms (0.64% GC)
  mean time:        801.359 ms (0.64% GC)
  maximum time:     807.842 ms (0.63% GC)
  --------------
  samples:          7
  evals/sample:     1


 julia> @benchmark f2 = fit2(tdata,ydata,IRF)
  BenchmarkTools.Trial: 
    memory estimate:  255.31 MiB
    allocs estimate:  1578872
    --------------
    minimum time:     212.654 ms (4.39% GC)
    median time:      214.032 ms (4.38% GC)
    mean time:        218.092 ms (4.41% GC)
    maximum time:     242.318 ms (3.91% GC)
    --------------
    samples:          23
    evals/sample:     1
  =#



  su = size(IRF)
  outsize = 2 .* su .-1
 # nffts = nextfastfft(outsize)


  upad = similar(IRF, outsize)
  vpad = similar(IRF, outsize)
  fill!(upad, zero(eltype(upad)))
  fill!(vpad, zero(eltype(vpad)))

  copyto!(upad, IRF)

  p! = plan_fft!(convert(Array{Complex{Float64}}, upad))
  upad = p! * upad




function conv_DT(t, data::input_data, upad, p!, vpad)

	RtDT = Array{Float64}(undef, length(0:0.01:10))
    RtDT = DT_model(0:0.01:10, [0.01,10.0], 1., 1., 1.)

    copyto!(vpad, RtDT)


    vpad = p! * vpad
    vpad .*= upad
    real(ifft!(vpad))

	Rt = Rt./maximum(Rt)

	tidx = findfirst(x -> x == t[1], data.t)

	convDT = convDT[tidx:tidx+length(t)-1]
	return log.(convDT)
end


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


#=@benchmark conv_DT(data.t[ind1:ind2], [0.1,10.], data, vpad, RtDT, pl, vpadRfft, upadRfft, u_vpadRfft, convpad)
BenchmarkTools.Trial: 
memory estimate:  401.22 KiB
allocs estimate:  134
--------------
minimum time:     194.223 μs (0.00% GC)
median time:      199.783 μs (0.00% GC)
mean time:        241.692 μs (2.87% GC)
maximum time:     5.035 ms (0.00% GC)
--------------
samples:          10000
evals/sample:     1


old version
@benchmark conv_DT(data.t[ind1:ind2], [0.1,10.], data)
BenchmarkTools.Trial: 
  memory estimate:  628.81 KiB
  allocs estimate:  184
  --------------
  minimum time:     296.018 μs (0.00% GC)
  median time:      308.628 μs (0.00% GC)
  mean time:        423.330 μs (2.68% GC)
  maximum time:     7.641 ms (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
=#



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