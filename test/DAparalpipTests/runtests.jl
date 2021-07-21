module DAparalpipTests

using Test
using LightPropagation: fluence_DA_semiinf_TD
using LightPropagation: fluence_DA_slab_TD
using LightPropagation: fluence_DA_paralpip_TD

# Test parallelepied against semi-inf and slabs

# TD 
t = 0.1:0.5:8

# compare against semi-inf
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, n_ext = 1.0, n_med = 1.0, rd = [4.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_TD(t, 1.0, 0.1, 10.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  
@test fluence_DA_paralpip_TD(t, 0.1, 70.0, n_ext = 1.0, n_med = 1.0, rd = [2.0, 5.0, 0.0], rs = [5.0,5.0], L = [10.0,10.0,10.0]) ≈  fluence_DA_semiinf_TD(t, 3.0, 0.1, 70.0,  n_med = 1.0, n_ext = 1.0, z = 0.0)  

# compare against slabs 
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, 1.0, 1.0, [4.0,5.0,0.0], [5.0,5.0], [10.0,10.0,10.0])  ≈ fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, 1.0, 1.0, 10.0, 0.0)
@test fluence_DA_paralpip_TD(t, 0.1, 10.0, 1.0, 1.0, [4.0,5.0,0.0], [5.0,5.0], [10.0,10.0,0.5])  ≈ fluence_DA_slab_TD(t, 1.0, 0.1, 10.0, 1.0, 1.0, 0.5, 0.0)

end # module