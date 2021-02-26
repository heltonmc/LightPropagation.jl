
besselroots = load("besselzeroroots.jld")["besselroots"]

@with_kw struct cylinder_inputs1
    μsp::Float64 = 10.0
    μa::Float64 = 0.1 
    n_ext::Float64 = 1.0  #surrounding index of refraction
    n_med::Float64 = 1.4 # layers index of refraction

    #source, detector
    lz::Float64 = 5.0 # length of cylinder
    ρ::Float64 = 1.0 # 
    a::Float64 = 5.0 # radius of cylinder
    z::Float64 = 0.0
    pw::Float64 = 0.1

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
    g -= exp(α*(z0 + z - 2*lz - 2*zb))*(1 - exp(-2*α*(z0 + zb)))*(1 - exp(-2*α*(zb + z)))/(1 - exp(-2*α*(lz + 2*zb)))
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

function fluence_DA_cylinder_CW_beam(data::cylinder_inputs1, besselroots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0
    sn = 0.0
    for n in eachindex(besselroots)
        sn = besselroots[n]/(data.a + zb)
        ϕ += _greencylin_CW(sn, data.μa, D, data.z, z0, zb, data.lz)*besselj0(sn*data.ρ)*besselj1(sn*data.pw)/(sn*(besselj1(sn*(data.a + zb))^2))
    end

    return 2*ϕ/(data.pw*π*(data.a + zb)^2)
end

function _fluence_DA_cylinder_Laplace_beam(s, data::cylinder_inputs1, besselroots)

    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)
    ϕ = 0.0

    sn = 0.0
    for n in eachindex(besselroots)
        sn = besselroots[n]/(data.a + zb)
        ϕ += _greencylin_Laplace(sn, data.μa, D, data.z, z0, zb, data.lz, ν, s)*besselj0(sn*data.ρ)*besselj1(sn*data.pw)/(sn*(besselj1(sn*(data.a + zb))^2))
    end

    return 2*ϕ/(data.pw*π*(data.a + zb)^2)
end

function fluence_DA_cylinder_TD_beam(t, N, data::cylinder_inputs1, besselroots)
    Rt = zeros(Float64, length(t))
    Rt = LT_hyper_fixed(s -> _fluence_DA_cylinder_Laplace_beam(s, data, besselroots), N, t)

    return Rt./59.958469
end

### arbitrary position


@with_kw struct cylinder_inputs_arb3
    μsp::Float64 = 10.0
    μa::Float64 = 0.1 
    n_ext::Float64 = 1.0  #surrounding index of refraction
    n_med::Float64 = 1.4 # layers index of refraction

    #source, detector
    lz::Float64 = 5.0 # length of cylinder
    ρ::Float64 = 5.0 # 
    ρ0::Float64 = 5.0
    a::Float64 = 5.0 # radius of cylinder
    z::Float64 = 2.0
    z0::Float64 = 3.0
    
    θ = 0.0

    ω::Float64 = 0.0 #
end


function fluence_DA_cylinder_CW(data::cylinder_inputs_arb2, besselj_roots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = 0.0
    ϕ1 = 0.0

    for m = 0:100
        ϕ1 = 0.0
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            ϕ1 += _greencylin_CW(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz)*besselj(m, bessel/(data.a + zb)*data.ρ)*besselj(m, bessel/(data.a + zb)*data.ρ0)/(besselj(m + 1, bessel))^2

        end
        ϕ1 *= ((2 - (m == 0 ? 1 : 0))*cos(m*data.θ))

        ϕ += ϕ1
    end

    return ϕ / (π*(data.a + zb)^2)
end



function fluence_DA_cylinder_CW(data::cylinder_inputs_arb3, besselj_roots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ρ0 = data.ρ0 -  1/data.μsp[1]
    z0 = data.z0

    ϕ = 0.0
    ϕ1 = 0.0

    for m = 0:1200
        ϕ1 = 0.0
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            ϕ1 += _greencylin_CW(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz)*besselj(m, bessel/(data.a + zb)*data.ρ)*besselj(m, bessel/(data.a + zb)*ρ0)/(besselj(m + 1, bessel))^2

        end
        #println(ϕ1)
        ϕ1 *= ((2 - (m == 0 ? 1 : 0))*cos(m*data.θ))

        ϕ += ϕ1
        #println(ϕ)
    end

    return ϕ / (π*(data.a + zb)^2)
end


## work in process: need to test best convergence
## can't break for order m when theta is non-zero. You get oscillating term
##appears to be a bug when z = z0. Slow convergence or something off? Seems when theta = 0 everything is ok...
function fluence_DA_cylinder_CW(data::cylinder_inputs_arb3, besselj_roots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ρ0 = data.ρ0 -  1/data.μsp[1]
    z0 = data.z0

    ϕ = 0.0
    ϕ1 = 0.0

    root_old = 0.0
    order_old = 0.0
    count = 0

    TOL = 1e-12

    for m = -1600:1600#0:8000
        ϕ1 = 0.0
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            root_old = ϕ1
            ϕ1 += _greencylin_CW(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz)*besselj(m, bessel/(data.a + zb)*data.ρ)*besselj(m, bessel/(data.a + zb)*ρ0)/(besselj(m + 1, bessel))^2

            if (abs(ϕ1 - root_old))/ϕ1 < TOL
                #println(bessel)
                break
            end
        end
        
     

        ϕ1 *= #=((2 - (m == 0 ? 1 : 0))*=#cos(m*data.θ)
        order_old = ϕ1

        ϕ += ϕ1
    end

    return ϕ / (π*(data.a + zb)^2)
end




function _fluence_DA_cylinder_Laplace(s, data::cylinder_inputs_arb3, besselj_roots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ρ0 = data.ρ0 -  1/data.μsp[1]
    z0 = data.z0

    ϕ = 0.0
    ϕ1 = 0.0

    for m = 0:600
        ϕ1 = 0.0
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            ϕ1 += _greencylin_Laplace(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz, ν, s)*besselj(m, bessel/(data.a + zb)*data.ρ)*besselj(m, bessel/(data.a + zb)*ρ0)/(besselj(m + 1, bessel))^2

        end
        ϕ1 *= ((2 - (m == 0 ? 1 : 0))*cos(m*data.θ))

        ϕ += ϕ1
    end

    return ϕ / (π*(data.a + zb)^2)
end
function fluence_DA_cylinder_TD(t, N, data::cylinder_inputs_arb3, besselroots)
    Rt = zeros(Float64, length(t))
    Rt = LT_hyper_fixed(s -> _fluence_DA_cylinder_Laplace(s, data, besselroots), N, t)

    return Rt./59.958469
end
#=
Calculate Bessel Roots

using GSL

s = 1:1200 # number of roots to calculate
m = 0:1300 # order of besselj
besselj_roots = [sf_bessel_zero_Jnu_e(m, s).val for s in s, m in m]

=#