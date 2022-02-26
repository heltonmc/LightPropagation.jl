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

function hyperbola(f::Function, t::AbstractArray; N = 16)
    out = fill(f(t[1]), length(t))
    Threads.@threads for ind in 2:length(t)
        out[ind] = hyperbola(f, t[ind], N = N)
    end
    return out
end

#### fixed integration path (can precompute)
# get the fixed integration components

# get the fixed integration components
function hyper_coef(N, t::AbstractVector; ϕ = 1.09)
    A = acosh(((π - 2 * ϕ) * t[end] / t[1] + 4 * ϕ - π) / ((4 * ϕ - π) * sin(ϕ)))
    μ = (4 * π * ϕ - π^2) * N / t[end] / A
    h = A / N
    return μ, h
end

# compute the fixed hyperbola contour
s_fixed(θ, μ; ϕ = 1.09) = μ + im * μ * sinh(θ + im * ϕ)
ds_fixed(θ, μ; ϕ = 1.09) = im * μ * cosh(θ + im * ϕ)

# compute the function values over the fixed contour nodes
function fixed_sk(f::Function, N, t::AbstractVector, T)
    μ, h = hyper_coef(N, t)
    a = zeros(Complex{eltype(T)}, N)
    sk = zeros(Complex{eltype(t)}, N)
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

function hyper_fixed(f::Function, t::AbstractVector; N = 24, T = eltype(t))
    a, sk, h = fixed_sk(f, N, t, T)
    out = zeros(T, length(t))
    Threads.@threads for ind in eachindex(t)
        out[ind] = hyper_fixed_points(a, sk, h, t[ind])
    end
    return out
end


#=

# compute the fixed hyperbola contour
s_fixed(θ, μ; ϕ = 1.09) = μ + im * μ * sinh(θ + im * ϕ)
ds_fixed(θ, μ; ϕ = 1.09) = im * μ * cosh(θ + im * ϕ)

# compute the function values over the fixed contour nodes
function fixed_sk_array(f::Function, N, M, t::AbstractVector, T)
    μ, h = hyper_coef(N, t)
    a = zeros(Complex{eltype(T)}, N, M)
    sk = zeros(Complex{eltype(t)}, N)
    for k in 0:Int(N)-1
        sk[k+1] = s_fixed((k + 1/2) * h, μ)
        dsk = ds_fixed((k + 1/2) * h, μ)
        a[k+1, :] = f(sk[k + 1]) * dsk
    end
    return a, sk, h
end

function hyper_fixed_points_array(a, sk, h, t::AbstractFloat, M)
    b = zeros(eltype(a), M)
    for ind in eachindex(sk)
        b += a[ind, :] * exp(sk[ind] * t)
    end
    return imag(b) * h / π
end

function hyper_fixed_array(f::Function, t::AbstractVector; N = 24, M = 2, T = eltype(t))
    a, sk, h = fixed_sk_array(f, N, M, t, T)
    out = zeros(T, length(t), M)
    for ind in eachindex(t)
        out[ind, :] = hyper_fixed_points_array(a, sk, h, t[ind], M)
    end
    return out
end
=#