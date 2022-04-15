module DAinfTests

using Test
using LightPropagation

struct DAsemiinf_test{T <: AbstractFloat} <: DiffusionParameters
    ρ::T                                 # distance away from isotropic point source (cm)
    μa::T                                # absorption coefficient (cm⁻¹)
    μsp::T                               # reduced scattering coefficient (cm⁻¹)
    n_med::T                             # medium's index of refraction
    n_ext::T                             # external medium's index of refraction
    ω::T                                 # modulation frequency (1/ns)
    z::T                                 # z depth of detector within medium (cm)
    t::Union{T, AbstractVector{T}}       # time vector (ns)

    # do not need to provide these arguments (calculated from previous inputs)
    D::T                                 # Diffusion coefficient                        
    ν::T                                 # Speed of light in medium (cm/ns)
    z0::T                                # isotroptic source depth
    zb::T                                # Extrapolation length
end

# this manually changes the definition of the isotropic source depth to bury deep within a semi-infinite medium
# this approximates an infinite medium when we are close to source to compare the infinite and semi-infinite solutions

#SI_test = DAsemiinf_test(0.0, 0.1, 10.0, 1.0, 1.0, 0.0, 20.0, 0.1:0.1:2.5, 1.0 / (3.0 * 10.0), 29.9792458, 19.0, 2.0 / (3.0 * 10.0))
#@test fluence_DA_semiinf_CW(SI_test) ≈ fluence_DA_inf_CW(1.0, 0.1, 10.0)

#t = 0.1:0.1:2.5
#@test fluence_DA_semiinf_TD(SI_test) ≈ fluence_DA_inf_TD(t, 1.0, 0.1, 10.0)

# Hard code DA inf solution to compare other geometries

# CW
@test fluence_DA_inf_CW(0.1, 0.1, 10.0) ≈ 20.076563644369305
@test fluence_DA_inf_CW(0.1, 0.2, 15.0) ≈ 26.52859839465854
@test fluence_DA_inf_CW(60.0, 0.2, 15.0) ≈ 4.007233568620553e-80

# TD
d = [5.87809115624927e-14, 0.0019712791801801076, 2.10256922566133e-6, 3.038559267148949e-9]
@test fluence_DA_inf_TD(0.01:1.0:4.0, 1.0, 0.2, 15.0) ≈ d

end # module
