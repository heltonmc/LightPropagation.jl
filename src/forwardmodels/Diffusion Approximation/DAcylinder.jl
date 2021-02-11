
besselroots = load("besselzeroroots.jld")["besselroots"]

@with_kw struct cylinder_inputs
    μsp::Float64 = 10.0
    μa::Float64 = 0.1 
    n_ext::Float64 = 1.0  #surrounding index of refraction
    n_med::Float64 = 1.4 # layers index of refraction

    #source, detector
    lz::Float64 = 5.0 # length of cylinder
    ρ::Float64 = 1.0 # 
    a::Float64 = 5.0 # radius of cylinder
    z::Float64 = 0.0

    ω::Float64 = 0.0 #
end

# look at get_afac for some reason gives type instabillity
function diffusionparams(μsp, n_med, n_ext)
    ## Diffusion parameters
    D = 1/3μsp
    ν = 29.9792345/n_med
    A = 1.0 #get_afac(n_ext/n_med) # need to calculate reflection between layers and surrounding medium
    zb = 2*A*D
    z0 = 1/(μsp)

    return D, ν, A, zb, z0
end

function _greencylin(α, D, z, z0, zb, lz)
    g = (exp(-α*abs(z - z0)) - exp(-α*(z + z0 + 2*zb)))
    g -= exp(α*(z0 - 2*lz - 2*zb))*(1 - exp(-2*α*(z0 + zb)))*(1 - exp(-2*α*zb))/(1 - exp(-2*α*(lz + 2*zb)))
    return g / (2*D*α)
end

function _greencylin_CW(sn, μa, D, z, z0, zb, lz)
    α = sqrt(μa/D + sn^2)
    return _greencylin(α, D, z, z0, zb, lz)
end

function _greencylin_FD(sn, μa, D, z, z0, zb, lz, ν, ω)
    α = sqrt(μa/D + sn^2 + im*ω/(D*ν))
    return _greencylin(α, D, z, z0, zb, lz)
end

function _greencylin_Laplace(sn, μa, D, z, z0, zb, lz, ν, s)
    α = sqrt(μa/D + sn^2 + s/(D*ν))
    return _greencylin(α, D, z, z0, zb, lz)
end

function fluence_DA_cylinder_CW(data::cylinder_inputs, besselroots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0
    for bessel in besselroots
        ϕ += _greencylin_CW(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz)*besselj0(bessel/(data.a + zb)*data.ρ)/(besselj1(bessel))^2
    end

    return ϕ / (π*(data.a + zb)^2)
end


function fluence_DA_cylinder_FD(data::cylinder_inputs, besselroots)

    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0 + 0.0im

    for n in eachindex(besselroots)
        ϕ += _greencylin_FD(besselroots[n]/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz, ν, data.ω)*besselj0(besselroots[n]/(data.a + zb)*data.ρ)/(besselj1(besselroots[n]))^2
    end

    return ϕ / (π*(data.a + zb)^2)
end

### Compute Time Domain Signal

# calulcate fluence in laplace space: s -> im*ω
function _fluence_DA_cylinder_Laplace(s, data::cylinder_inputs, besselroots)

    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0

    for bessel in besselroots
        ϕ += _greencylin_Laplace(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz, ν, s)*besselj0(bessel/(data.a + zb)*data.ρ)/(besselj1(bessel))^2
    end

    return ϕ/(π*(data.a + zb)^2)
end

# use adaptive hyperbola contour for each time point Laplace Transform
# much more accurate but slower. N = 8 is usually sufficient
# need larger N for very short times (< 0.1 ns) and very long times (> 3ns)
function fluence_DA_cylinder_TD(t, N, data::cylinder_inputs, besselroots)
    Rt = zeros(Float64, length(t))
    Rt = LT_hyperbola(s -> _fluence_DA_cylinder_Laplace(s, data, besselroots), N, t)

    return Rt./59.958469
end

# use fixed contour for every time point
# signifanctly faster (10-100x) than adaptive contour but less accurate
# utilizing this is usually ok but further testing is needed 
# N = 28 is usuall pretty accurate
# Need about 600 Bessel Zeros
# the 59.9 is dividing it by some constant to match the fluence in semi-inf. Unsure why this is needed (not needed in the steady-state: match exactly)
function fluence_DA_cylinder_TD(t, N, data::cylinder_inputs, besselroots)
    Rt = zeros(Float64, length(t))
    Rt = LT_hyper_fixed(s -> _fluence_DA_cylinder_Laplace(s, data, besselroots), N, t)

    return Rt./59.958469
end


### for a finite flat beam incident onto center of cylinder top
# have data structure contain pw (radius of beam)

function fluence_DA_cylinder_CW_beam(data::cylinder_inputs, besselroots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0
    sn = 0.0
    for n in eachindex(besselroots)
        sn = besselroots[n]/(data.a + zb)
        ϕ += _greencylin_CW(sn, data.μa, D, data.z, z0, zb, data.lz)*besselj0(sn*data.ρ)*besselj1(sn*data.pw)/(sn*(besselj1(sn*(data.a + zb))^2))
    end

    return 2*ϕ/(data.pw*π*(data.a + zb)^2)
end
