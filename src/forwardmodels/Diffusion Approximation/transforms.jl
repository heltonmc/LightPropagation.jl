# Compute Post-Widder coefficients Vk
function _PWcoeffs(N::Integer)
    if isodd(N)
        N += 1
        @warn "N must be even... increment N += 1"
    end
    v = zeros(BigFloat, N)
    aux = zero(eltype(v))
    for k in 1:N
        for j in fld(k + 1, 2):minimum(Int, [k, div(N, 2)])
            aux = big(j)^(div(N, 2)) * factorial(big(2 * j))
            aux /= factorial(big(div(N, 2) - j)) * factorial(big(j)) * factorial(big(j - 1))
            aux /= factorial(big(k - j)) * factorial(big(2 * j - k))
            v[k] += aux
        end
        v[k] *= (-1)^(k + div(N, 2))
    end
    return v
end

# compute the laplace transform using Post-Widder for a single time pt
function postwid(func::Function, t::AbstractFloat; v = _PWcoeffs(18))
    N = length(v)
    a = zero(eltype(v))
    for k in 1:N
        bk = convert(eltype(v), k)
        a += v[k] * func(bk * log(convert(eltype(a), 2)) / t)
    end
    return a * log(convert(eltype(a), 2)) / t
end

# for an array of t using multithreading
function postwid(f::Function, t::AbstractArray; v = _PWcoeffs(18))
    N = length(v)
    a = zeros(eltype(v), length(t))
    Threads.@threads for ind in eachindex(t)
        a[ind] = postwid(f, t[ind], v = v)
    end
    return a
end


#=
####### How to use
function fluence_DA_inf_FD(ρ, β, s, nmed = 1.0)
    μa = β[1]
    μsp = β[2]
    D = 1/3μsp
    ν = 29.9792458/nmed

    ϕ = exp(-ρ*sqrt(μa/D + s/(ν*D)))/(4*π*ρ*D)

    return ϕ
end

function fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64 = 1.0)
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
    ν::Float64 = 29.9792458/nmed

    ϕ = Array{Float64}(undef, length(t))

    Threads.@threads for n in eachindex(t)

        ϕ[n] = -(ρ^2/(4D*ν*t[n]))
        ϕ[n] = ϕ[n] - μa*ν*t[n]
        ϕ[n] = ν*exp(ϕ[n])/((4π*D*ν*t[n])^(3/2))

        if isnan(ϕ[n])
            ϕ[n] = 0
        end
 
    end

    return ϕ
end


function computeTD_fromFD2(t, β, ρ, ub, N)	

    Rt = zeros(Float64, length(t))
    freqcomp = zeros(Float64, N)
    x, w = gausslegendre(N)


    Threads.@threads for m in eachindex(x)
        freqcomp[m] = real(fluence_DA_inf_FD(ρ, [β[1],β[2]], x[m]*ub/2 +ub/2))*w[m]*ub/pi
    end

	
    Threads.@threads for n in eachindex(t)
        for m in eachindex(x)
            Rt[n] += real(exp(im*t[n]*(x[m]*ub/2 +ub/2)))*freqcomp[m]
        end
    end
   
   
return Rt
end

#=
N = 18
v = _PWcoeffs(N)
t = 0.01:0.01:5
lt = LT_postwid(s -> fluence_DA_inf_FD(1.0, [0.1, 10.0], s),v, t)
Rt = fluence_DA_inf_TD(t, [0.1, 10.0], 1.0)

plot(t, Rt, yscale=:log10, label = "TD analytical 'true'")
plot!(t, abs.(lt), label = "LT Post-Widder")
=#

=#
# Implement Laplace transfrom along a hyperbola contour

## adaptive contour for each time point (can't precompute f(ω) or sk)

function s(θ, N, t)
    μ = 4.492075287 * N / t
    ϕ = 1.172104229 
    return μ + im * μ * sinh(θ + im * ϕ)
end
# derivitive of hyperbola contour
function ds(θ, N, t)
    μ = 4.492075287 * N / t
    ϕ = 1.172104229
    return im * μ * cosh(θ + im * ϕ)
end

# compute the laplace transform along a hyperbola contour fixed time point
function hyperbola(f::Function, t::AbstractFloat; N = 16)
    a =  zero(Complex{eltype(t)})
    N = convert(eltype(t), N)
    h = 1.081792140 / N
    for k in 0:N-1
        sk = s((k + 1/2) * h, N, t)
        dsk = ds((k + 1/2) * h, N, t)
        a += f(sk) * exp(sk * t) * dsk
    end
    return imag(a) * h / pi
end

# loop through an entire time array with multithreads
function hyperbola(f::Function, t::AbstractArray; N = 16)
    out = similar(t)
    Threads.@threads for ind in eachindex(t)
        out[ind] = hyperbola(f, t[ind], N = N)
    end
    return out
end

#### fixed integration path (can precompute)
# get the fixed integration components

# get the fixed integration components
function hyper_coef(N, t::AbstractArray; ϕ = 1.09)
    A = acosh(((π - 2 * ϕ) * t[end] / t[1] + 4 * ϕ - π) / ((4 * ϕ - π) * sin(ϕ)))
    μ = (4 * π * ϕ - π^2) * N / t[end] / A
    h = A / N
    return μ, h
end

# compute the fixed hyperbola contour
s_fixed(θ, μ; ϕ = 1.09) = μ + im * μ * sinh(θ + im * ϕ)
ds_fixed(θ, μ; ϕ = 1.09) = im * μ * cosh(θ + im * ϕ)

# compute the function values over the fixed contour nodes
function fixed_sk(f::Function, N, t::AbstractArray)
    μ, h = hyper_coef(N, t)
    a = zeros(Complex{eltype(h)}, Int(N))
    sk = similar(a)
    Threads.@threads for k in 0:Int(N)-1
        sk[k+1] = s_fixed((k + 1/2) * h, μ)
        dsk = ds_fixed((k + 1/2) * h, μ)
        a[k+1] = f(sk[k + 1]) * dsk
    end
    return a, sk, h
end

function hyper_fixed_points(a, sk, h, t::AbstractFloat)
    b = zero(eltype(a))
    for ind in eachindex(sk)
        b += a[ind] * exp(sk[ind] * t)
    end
    return imag(b) * h / π
end



function hyper_fixed(f::Function, t::AbstractArray; N = 24)
    N = convert(eltype(t), N)
    a, sk, h = fixed_sk(f, N, t)
    out = zeros(eltype(h), length(t))
    Threads.@threads for ind in eachindex(t)
        out[ind] = hyper_fixed_points(a, sk, h, t[ind])
    end
    return out
end

# run as
# LT_hyper_fixed(s -> fluence_DA_inf_FD(1.0, [0.1, 10.0], s),28, t)

