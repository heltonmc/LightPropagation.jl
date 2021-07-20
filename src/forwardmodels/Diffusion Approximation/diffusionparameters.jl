function get_afac(n)
    if n > 1.0
        A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
        5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
    elseif n < 1.0
        A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
    else 
        A = 1.0
    end
    return A
end

@inline D_coeff(μsp, μa) = 1 / (3 * μsp)
@inline ν_coeff(n_med) = 29.9792458 / n_med
@inline z0_coeff(μsp) = 1 / μsp
@inline zb_coeff(A, D) = 2 * A * D

function diffusionparams(μsp, n_med, n_ext)
    ## Diffusion parameters
    D = @. 1/3μsp
    ν = @. 29.9792345/n_med
    A = @. get_afac(n_ext/n_med) # need to calculate reflection between layers and surrounding medium
    zb = @. 2*A*D
    z0 = 1/(μsp[1])

    return D, ν, A, zb, z0
end