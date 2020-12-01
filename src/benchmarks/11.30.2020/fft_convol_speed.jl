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