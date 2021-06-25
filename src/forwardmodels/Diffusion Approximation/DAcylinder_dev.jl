
@with_kw struct cylinder_inputs{T <: AbstractFloat}
    # medium optical properties
    μsp::T = 10.0      # Reduced scattering coefficient in cm-1
    μa::T = 0.1        # Absorption coefficient in cm-1
    n_ext::T = 1.0     # Index of refraction exterior
    n_med::T = 1.4     # Index of refraction inside cylinder

    #source, detector configuration
    lz::T = 5.0        # length of cylinder (cm)
    ρ::T = 1.0         # Detector ρ (cm) coordinate in cylindrical coordinates r = (ρ, ϕ, z)
    ρ0::T = 0.0        # Source ρ (cm) coordinate in cylindrical coordinates r = (ρ, ϕ, z)
    a::T = 5.0         # Radius of cylinder (cm)
    z::T = 0.0         # Detector z (cm) coordinate in cylindrical coordinates r = (ρ, ϕ, z)
    z0::T = 0.0        # Source z (cm) coordinate in cylindrical coordinates r = (ρ, ϕ, z)
    pw::T = 0.0        # Radius (cm) of circular flat beam

    θ::T = 0.0         # Angle difference between source and detector in cylin coordinates r = (ρ, ϕ, z)

    ω::T = 0.0         # Modulation frequency for frequency domain
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
    if z == z0
        g = - exp(-α*(z + z0 + 2*zb))
    else
        g = (exp(-α*abs(z - z0)) - exp(-α*(z + z0 + 2*zb)))
    end
    g -= exp(α*(z0 + z - 2*lz - 2*zb))*(1 - exp(-2*α*(z0 + zb)))*(1 - exp(-2*α*(zb + z)))/(1 - exp(-2*α*(lz + 2*zb)))

    return g / (2*D*α)
end
#=
function _greencylin(α, D, z, z0, zb, lz)
    
    g = -exp(-α*(2*z + 2*zb)) / (2*D*α)
    g -= exp(2*α*(z - lz - zb)) / (2*D*α) * 2*(1 - exp(-2*α*(z + zb))) / (1 - exp(-2*α*(lz + 2*zb)))

    return g
end
=#
function _greencylin_CW(sn, μa, D, z, z0, zb, lz)
    α = sqrt(μa/D + sn^2)
    return _greencylin(α, D, z, z0, zb, lz)
end

### for a finite flat beam incident onto center of cylinder top
# have data structure contain pw (radius of beam)

function fluence_DA_cylinder_CW_beam(data, besselroots)
   
    D, ν, A, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)

    ϕ = zero(eltype(data.ρ))
    sn = zero(eltype(data.ρ))

    for n in eachindex(besselroots)
        sn = besselroots[n]/(data.a + zb)
        ϕ += _greencylin_CW(sn, data.μa, D, data.z, z0, zb, data.lz)*besselj0(sn*data.ρ)*besselj1(sn*data.pw)/(sn*(besselj1(sn*(data.a + zb))^2))
    end

    return 2*ϕ/(data.pw*π*(data.a + zb)^2)
end

function fluence_DA_cylinder_CW_arb(data, D, zb, z0, besselj_roots)
    ϕ = zero(eltype(data.ρ))
    ϕ1 = zero(eltype(data.ρ))
    for m = 0:size(besselj_roots)[2] - 1 # order number
        ϕ1 = 0.0
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            ϕ1 += _greencylin_CW(bessel/(data.a + zb), data.μa, D, data.z, z0, zb, data.lz)*besselj(m, bessel/(data.a + zb)*data.ρ)*besselj(m, bessel/(data.a + zb)*data.ρ0)/(besselj(m + 1, bessel))^2
        end
        ϕ1 *= ((2 - (m == 0 ? 1 : 0))*cos(m*data.θ))
        ϕ += ϕ1
    end

    return ϕ / (π*(data.a + zb)^2)
end


function source_location(data, z0)
    # source located parallel to barrel of cylinder
    if data.ρ0 == data.a
        if data.z0 == zero(eltype(data.z0)) || data.z0 == data.lz # check if source is on corners of cylinder
            error("unclear source location on edge of cylinder")
        end
        ρ0 = data.ρ0 -  z0
        z0 = data.z0
    # source located on top of barrel
    elseif data.z0 == zero(eltype(data.z0)) && data.ρ0 < data.a # can't be on edge either
        ρ0 = data.ρ0
        z0 = z0 # this z0 is calculated with diffusionparams that defaults to 1/μsp
    elseif data.z0 > zero(eltype(data.z0)) && data.z0 < data.lz && data.ρ0 < data.a
        @warn "... assuming you are manually placing source inside cylinder and not modifying location"
        ρ0 = data.ρ0
        z0 = data.z0
    else # throw error if doesn't fit any conditions above. It could be that they have placed source term on bottom surface instead of top. Make z0 = 0 
        error("Unsure where source term is located please check ... could be that source term is located on bottom surface (z0 = lz) instead of top")
    end

    return ρ0, z0
end

function _fluencecylin_center!(ϕ, green::Function, a, D, μa, z, z0, zb, lz, ρ, besselroots)
    for bessel in besselroots # only need zero order bessel roots J0
        ϕ += green(bessel/(a + zb), μa, D, z, z0, zb, lz)*besselj0(bessel/(a + zb)*ρ)/(besselj1(bessel))^2
    end
    ϕ /= (π*(a + zb)^2)

    return ϕ
end

function _fluencecylin_arb!(ϕ, green::Function, a, D, μa, z, z0, zb, lz, ρ, ρ0, θ, besselj_roots)
    ϕ1 = zero(eltype(ρ))
    TOL = 1e-12

    for m = 0:size(besselj_roots)[2] - 1 # order number
        ϕ1 = zero(eltype(ρ))
        for bessel in besselj_roots[:, abs(m) + 1] # 1 column is 0 order, bessel roots symmetric J(-nu, x) = (-1)^2 J(nu, x)
            root_old = ϕ1
            ϕ1 += green(bessel/(a + zb), μa, D, z, z0, zb, lz)*besselj(m, bessel/(a + zb)*ρ)*besselj(m, bessel/(a + zb)*ρ0)/(besselj(m + 1, bessel))^2

            if (abs(ϕ1 - root_old))/ϕ1 < TOL
               # println(bessel)
               # break
            end
        end
        ϕ1 *= ((2 - (m == 0 ? 1 : 0))*cos(m*θ))
        ϕ += ϕ1
    end

    return ϕ / (π*(data.a + zb)^2)
end


function fluence_DA_cylinder_CW(data, besselj_roots)
    D, _, _, zb, z0 = diffusionparams(data.μsp, data.n_med, data.n_ext)
    ρ0, z0 = source_location(data, z0)
   
    ϕ = zero(eltype(data.ρ))
    ## pencil beam incident onto center of cylinder top
    if data.ρ0 == zero(eltype(data.ρ0)) && data.z0 == zero(eltype(data.z0))
        ϕ = _fluencecylin_center!(ϕ, _greencylin_CW, data.a, D, data.μa, data.z, z0, zb, data.lz, data.ρ, besselj_roots[:, 1]) # only need zero order bessel roots J0
    elseif data.pw > zero(eltype(data.pw))
        if data.ρ0 > zero(elytpe(data.pw)) || data.z0 > zero(eltype(data.z0))
            error("flat beam has to be located on cylinder top")
        end
        ϕ = fluence_DA_cylinder_CW_beam(data, besselj_roots[:, 1])
    else
        ϕ = _fluencecylin_arb!(ϕ, _greencylin_CW, data.a, D, data.μa, data.z, z0, zb, data.lz, data.ρ, ρ0, data.θ, besselj_roots)
        if z0 == data.z
            r = sqrt(data.ρ^2 + ρ0^2 - 2*data.ρ*ρ0*cos(data.θ))
            Gp = exp(-r*sqrt(data.μa/D))/(4*π*D*r)
            ϕ += Gp
        end
    end
    return ϕ
end


#=
Calculate Bessel Roots

using GSL

s = 1:1200 # number of roots to calculate
m = 0:1300 # order of besselj
besselj_roots = [sf_bessel_zero_Jnu_e(m, s).val for s in s, m in m]

s = 1:10000
b = [sf_bessel_zero_Jnu_e(m, s).val for s in s]

=#
