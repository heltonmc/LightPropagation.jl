"""
    j_fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64)

Compute the jacobian of fluence_DA_inf_TD to be used in Levenberg-Marquardt.

# Arguments
- `t`: the time vector (ns). 
- `β::Array{Float64,1}`: the optical properties [μa, μs'] (cm⁻¹)
- `ρ::Float64`: the source detector separation (cm⁻¹)
- `nmed::Float64`: the sample medium's index of refraction

# Examples
julia> j_fluence_DA_inf_TD(1:1:5, [0.1,10.0], 1.0, 1.0)
"""
function j_fluence_DA_inf_TD(t, β::Array{Float64,1}, ρ::Float64, nmed::Float64 = 1.0)
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
    ν::Float64 = 29.9792345/nmed

    jϕ = Array{Float64}(undef, length(t), 2)

    Threads.@threads for n in eachindex(t)

        aux = exp(-(ρ^2/(4*D*ν*t[n])) - μa*ν*t[n])

        jϕ[n, 1] = -3^(3/2)*t[n]*ν^2*μsp^(3/2)*aux/(8*(π*t[n]*ν)^(3/2))
        jϕ[n, 2] = 3^(5/2)*sqrt(μsp)*(2*t[n]*ν - ρ^2*μsp)/(32*(t[n]*π*ν)^(3/2)*t[n])*aux

    end

    return jϕ
end

t = 0.01:0.01:3

ydata =  fluence_DA_inf_TD(t, [0.01,1.], 1.0) # + 0.00001*randn(length(t))

function curvetest(ydata, t)

    fit = curve_fit((t, β) -> fluence_DA_inf_TD(t, β, 1.0), t, ydata, [0.1, 60.])
end

function curvetestj(ydata, t)
    fit = curve_fit((t, β) -> fluence_DA_inf_TD(t, β, 1.0), (t, β) -> j_fluence_DA_inf_TD(t, β, 1.0), t, ydata, [0.1, 60.])
end

function curvetestff(ydata, t)
    fit = curve_fit((t, β) -> fluence_DA_inf_TD(t, β, 1.0), t, ydata, [0.1, 60.]; autodiff=:finiteforward)
end


function j_refl_DA_slab_TD(t, β::Array{Float64,1}, ρ::Float64, ndet::Float64, nmed::Float64)
    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
  	ν::Float64 = 29.9792345/nmed
    Rt = Array{Float64}(undef, length(t))

	if n > 1.0
		  A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		  5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	elseif n < 1.0
		  A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	else 
		  A = 1.0
	end

    jϕ = Array{Float64}(undef, length(t), 2)

    Threads.@threads for n in eachindex(t)

        jϕ[n, 1] = (3*μsp)^(3/2)*exp(-t[n]*ν*μa - 3*ρ*μsp/(4*t[n]*ν))
        jϕ[n, 1] *= (-(4*A + 3)*exp(-(4*A + 3)^2/(12*t[n]*ν*μsp))/(3*μsp) - exp(-3/(4*t[n]*ν*μsp))/μsp)
        jϕ[n, 1] /= 16*π^(3/2)*t[n]^(3/2)*sqrt(ν)

        jϕ[n, 2] = ((27*ρ*μsp^2 - 18*t[n]*ν*μsp - 27)*exp((4*A + 3)^2/(12*t[n]*ν*μsp)))
        jϕ[n, 2] += ((36*A + 27)*ρ*μsp^2 + (-24*A - 18)*t[n]*ν*μsp - 64*A^3 - 144*A^2 - 108*A - 27)*exp(3/(4*t[n]*ν*μsp))
        jϕ[n, 2] *= -exp(-3*ρ*μsp/(4*t[n]*ν) - (4*A + 3)^2/(12*ν*t[n]*μsp) - 3/(4*t[n]*ν*μsp) - t[n]*ν*μa)
        jϕ[n, 2] /= 64*sqrt(3)*(π*μsp)^(3/2)*t[n]^(7/2)*ν^(5/2)
    end
    return jϕ

end




t = 0.01:0.01:3

ydata =  refl_DA_slab_TD(t, [0.01,1.], 1.0, 1.0, 1.0) # + 0.00001*randn(length(t))

function curvetest(ydata, t)

    fit = curve_fit((t, β) -> refl_DA_slab_TD(t, β, 1.0, 1.0, 1.0), t, ydata, [0.1, 20.])
end

#this is much slower.... jacobian takes >3x longer to calculate 
function curvetestj(ydata, t)
    fit = curve_fit((t, β) -> refl_DA_slab_TD(t, β, 1.0, 1.0, 1.0), (t, β) -> j_refl_DA_slab_TD(t, β, 1.0, 1.0, 1.0), t, ydata, [0.1, 20.])
end

function curvetestff(ydata, t)
    fit = curve_fit((t, β) -> refl_DA_slab_TD(t, β, 1.0, 1.0, 1.0), t, ydata, [0.1, 20.]; autodiff=:finiteforward)
end