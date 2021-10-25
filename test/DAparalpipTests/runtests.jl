module DAparalpipTests

using Test
using LightPropagation: fluence_DA_semiinf_CW, fluence_DA_semiinf_TD
using LightPropagation: fluence_DA_slab_CW, fluence_DA_slab_TD
using LightPropagation: fluence_DA_paralpip_CW, fluence_DA_paralpip_TD
using LightPropagation: flux_DA_semiinf_TD, flux_DA_paralpip_TD, flux_DA_slab_TD

# Test parallelepied against semi-inf and slabs

# CW

# test against semi-inf
@test fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 10) ≈ 0.024035998288332954
@test fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_CW(1.0, 0.1, 10.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_CW(0.23, 25.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_CW(1.0, 0.23, 25.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_CW(0.4, 42.0, n_ext = 1.0, n_med = 1.0, rd = [3.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_CW(2.0, 0.4, 42.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_CW(0.1, 30.0, n_ext = 1.2, n_med = 1.4, rd = [3.5, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0], xs = 20) ≈  fluence_DA_semiinf_CW(1.5, 0.1, 30.0,  n_med = 1.4, n_ext = 1.2, z = 0.0)  

# test against slab
@test fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_slab_CW(1.0, 0.1, 10.0, s = 10.0, n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_CW(0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,0.5]) ≈  fluence_DA_slab_CW(1.0, 0.1, 10.0, s = 0.5, n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_CW(0.2, 21.0, n_ext = 1.1, n_med = 1.3, rd = [3.8, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,0.2]) ≈  fluence_DA_slab_CW(1.2, 0.2, 21.0, s = 0.2, n_med = 1.3, n_ext = 1.1, z = 0.0)  
@test fluence_DA_paralpip_CW(0.2, 21.0, n_ext = 1.1, n_med = 1.3, rd = [3.8, 5.0, 0.5], rs = [5.0,5.0], L = [10.0,10.0,1.0]) ≈  fluence_DA_slab_CW(1.2, 0.2, 21.0, s = 1.0, n_med = 1.3, n_ext = 1.1, z = 0.5)  

# TD 
t = 0.1:0.5:8

# compare against semi-inf
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 10) ≈  fluence_DA_semiinf_TD(t, 1.0, 0.1, 10.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_TD(t, 0.1, 70.0, n_ext = 1.0, n_med = 1.0, rd = [2.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_TD(t, 3.0, 0.1, 70.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  

# compare against slabs 
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 10.0], xs = 10) ≈  fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, s = 10.0, n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 0.5], xs = 10) ≈  fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, s = 0.5, n_med = 1.0, n_ext = 1.0, z = 0.0)  

@test fluence_DA_paralpip_TD(t, 0.2, 26.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0, 5.0], L = [10.0, 10.0, 0.2], xs = 10) ≈  fluence_DA_slab_TD(t, 1.0, 0.2, 26.0, s = 0.2, n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_TD(t, 0.4, 32.0, n_ext = 1.0, n_med = 1.0, rd = [3.5, 5.0, 0.1], rs = [5.0, 5.0], L = [10.0, 10.0, 0.2], xs = 10) ≈  fluence_DA_slab_TD(t, 1.5, 0.4, 32.0, s = 0.2, n_med = 1.0, n_ext = 1.0, z = 0.1)

# flux against semi-inf (reflectance)
@test flux_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [24.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 10) ≈  flux_DA_semiinf_TD(t, 1.0, 0.1, 10.0,  n_med = 1.0, n_ext = 1.0)  
@test flux_DA_paralpip_TD(t, 0.2, 12.0, n_ext = 1.2, n_med = 1.0, rd = [23.0, 25.0, 0.0], rs = [25.0,25.0], L = [50.0,50.0,50.0], xs = 10) ≈  flux_DA_semiinf_TD(t, 2.0, 0.2, 12.0,  n_med = 1.0, n_ext = 1.2)  

# flux against slabs (transmittance)
@test flux_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [5.0, 5.0, 1.0], rs = [5.0, 5.0], L = [10.0, 10.0, 1.0], xs = 10) ≈  flux_DA_slab_TD(t, 0.0, 0.1, 10.0, s = 1.0, n_med = 1.0, n_ext = 1.0, z = 1.0)  

end # module
