module DAcylinderTests

using Test
using LightPropagation: fluence_DA_semiinf_CW
using LightPropagation: fluence_DA_cylinder_CW


## How to calculate Bessel Roots
using GSL
s = 1:10000 # number of roots to calculate
m = 0:10000 # order of besselj
besselj_roots = [sf_bessel_zero_Jnu_e(m, s).val for s in s, m in m]

##############################################################################################################################
# Test cylinder against semi-infinite approximation when beam is located on top middle
##############################################################################################################################
## for these tests only need to consider first 600 bessel roots (only calculating 0 order)

### test for μsp = 10.0, μa = 0.1, ρ = 1.0
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.0, ρ0 = 0.0, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, :]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-7)

### test for μsp = 10.0, μa = 0.1, ρ = 2.0
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 2.0, ρ0 = 0.0, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, :]), fluence_DA_semiinf_CW(2.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-5)

### test for μsp = 10.0, μa = 0.1, ρ = 0.3
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 0.3, ρ0 = 0.0, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, :]), fluence_DA_semiinf_CW(0.3, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-5)

## test for μsp = 10.0, μa = 0.1, ρ = 1.0 and big lz = 20.0 (may be a convergence issue when lz is large)
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 20.0, ρ = 1.0, ρ0 = 0.0, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, :]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-7)

## test for μsp = 40.0, μa = 0.3, ρ = 1.0. Warning - you need to consider larger amount of roots here for large scattering. Using 1200
data = cylinder_inputs(μsp = 40.0, μa = 0.3, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.0, ρ0 = 0.0, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:1200, :]), fluence_DA_semiinf_CW(1.0, [0.3, 40.0], 1.0, 1.0, 0.0), rtol = 1e-5)

##############################################################################################################################
# Test cylinder against semi-infinite approximation when beam is located on top but source and detector not in middle (lack symmetry)
##############################################################################################################################

### We are going to test for source and detector located off the middle but still test for when |r - r0| = 1.
# So the semi-inf SDS still is 1.0 but the values in cylindrical coordinates will be different.
# test for μsp = 10.0, μa = 0.1, ρ = 0.4, ρ0 = 0.7, θ = 2.24593 which should give us a side length of triangle of 1.0
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 0.4, ρ0 = 0.7, 
                       a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 2.24593, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:60]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-5)

# move slightly off axis so ρ = 1.1, and ρ = 0.1
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.1, ρ0 = 0.1, 
                              a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:60]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-7)

# move slightly off axis so ρ = 1.1, and ρ = 0.1 and increase μsp = 40.0, μa = 0.3, need to increase number of roots
data = cylinder_inputs(μsp = 40.0, μa = 0.3, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.1, ρ0 = 0.1, 
                              a = 5.0, z= 0.0, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:1200, 1:100]), fluence_DA_semiinf_CW(1.0, [0.3, 40.0], 1.0, 1.0, 0.0), rtol = 1e-5)

## General comments: If the source/detector are located on the top the number of roots seem to have the highest effect on accuracy.
## As we increase the scattering we need more roots. If the source is not in middle we need to consider summation over orders as well.
## As you increase the scattering you will need to consider a larger amount of roots and a larger number of orders.
## Though when they are on the top,  ~100 seems to be sufficient.

##############################################################################################################################
# Test cylinder against semi-infinite approximation when beam is located on barrel and detector is on axis (ϕ = 0, changing z)
##############################################################################################################################
## these tests should represent the reflectance geometry. Measuring reflectance down the axis of the finger

### We are going to test for source and detector located off the middle but still test for when |r - r0| = 1.
# So the semi-inf SDS still is 1.0 but the values in cylindrical coordinates will be different.
# test for μsp = 10.0, μa = 0.1, ρ = 30.0, ρ0 = 30.0, θ = 0.0, lz = 5.0, z = 5.0, z0 = 4.0 
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:400]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-2)
# as you can see a curved geometry is much more noticeable here we have to drop the relative error to two significant digits. So within 99% of each other.

## let's move a little bit away from the source (z = 3.5, z = 5.5)
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 3.5, z0 = 5.5, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:400]), fluence_DA_semiinf_CW(2.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-1)
## had to drop it down but within 98% accuracy here.... this scenario looks to converge rather quickly so I think the differences are differences in geometries

# test for μsp = 40.0, μa = 0.3, ρ = 30.0, ρ0 = 30.0, θ = 0.0, lz = 5.0, z = 5.0, z0 = 4.0 
data = cylinder_inputs(μsp = 40.0, μa = 0.3, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:400]), fluence_DA_semiinf_CW(1.0, [0.3, 40.0], 1.0, 1.0, 0.0), rtol = 1e-2)
## one thing to note is the convergence primarily depends on |z - z0|. If the beam is located on the top this difference depends on the scattering coefficient.
## Which means that the number of roots to consider when the beam is on top will depend on how high the scattering coefficient is.
## When the beam is located on the barrel |z - z0| no longer depends on the scattering coefficient so convergence will depend less on scattering coefficient.
## The issue though is that |z - z0| is generally small on the barrel or worse equal to 0 which will not converge. If |z - zo| > 0.1 we are generally ok.
## Should probably think of erroring if |z - zo| < 0.01 and considering the special case when z = z0....

##############################################################################################################################
# Test cylinder against semi-infinite approximation when beam is located on barrel and detector has non-zero θ and z != z0
##############################################################################################################################
# this probably does not represent an experimental measurement nor semi-infinite well but it is to separate special case when z = z0

data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.5, pw = 0.0, θ = 0.028916732287767113, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:800]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-1)
# this is within 98% of each other which is pretty expected. 

data = cylinder_inputs(μsp = 40.0, μa = 0.3, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.5, pw = 0.0, θ = 0.028916732287767113, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:800]), fluence_DA_semiinf_CW(1.0, [0.3, 40.0], 1.0, 1.0, 0.0), rtol = 1e-1)

## let's test when z is close to z0; z-z0 = 0.1.
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.9, pw = 0.0, θ = 0.033223172543273956, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:1600]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-1)
# the test still passes but we need to consider signficantly higher orders to reach convergence
## let's test when z is close to z0; z-z0 = 0.01.
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 4.99, pw = 0.0, θ = 0.033388890726077015, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:1800, 1:1800]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-1)
# This can still pass but it is very slow to converge. 

##############################################################################################################################
# Test cylinder against semi-infinite approximation when beam is located on barrel and z = z0
##############################################################################################################################
# this is the hardest case for convergence. As we saw in the previous tests moving closer to the same z axis increases order number we need
# now analytically this also is a special case where one of the terms becomes one but as a suggestion from 
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 10.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 5.0, z0 = 5.0, pw = 0.0, θ = 0.03339056045285304, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:800]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.0), rtol = 1e-1)
# can not get good convergence for this case. Have to consider a ridiculous number of roots and orders 10,000 which takes hours to simulate

# place source manually inside cylinder in the middle top but source and detector located at 0.1.
# This should simulate well the semi-infinite case but detected at 0.1. This converges rapidly and very accurate.
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.0, ρ0 = 0.0, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:4]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-7)

data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 1.1, ρ0 = 0.1, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:12]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-7)
# this also passes but need to consider enough orders

data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 3.0, ρ0 = 2.0, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:50]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-7)
# this also passes but need to consider higher orders

# let's put the source on top z0 = 0.0 and it will calculate the source depth and be the same as previous manually putting in source
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 3.0, ρ0 = 2.0, 
                              a = 10.0, z= 0.1, z0 = 0.0, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:50]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-7)
# this also passes but need to consider higher orders

## test at same depth at 1/μsp but let's move off axis with source and detector

data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 0.4, ρ0 = 0.7, 
                       a = 5.0, z= 0.1, z0 = 0.0, pw = 0.0, θ = 2.24593, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:60]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-5)

## please not for these tests are carefully considered you can't arbitrarily change scattering and absorption properties here.
## The source depth depends on the optical properties as well as zb. If you change optical properties be aware of how that affects source location.


data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 30.0, ρ0 = 30.0, 
                       a = 30.0, z= 1.0, z0 = 1.0, pw = 0.0, θ = 0.03339056045285304, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:1000, 1:1000]), fluence_DA_semiinf_CW(1.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-5)


data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 6.1, ρ0 = 0.1, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:120]), fluence_DA_semiinf_CW(6.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-3)
# this also passes but need to consider enough orders

# this passes
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 6.1, ρ0 = 4.1, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:120]), fluence_DA_semiinf_CW(2.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-3)

# but this doesn't...
data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 5.0, ρ = 9.8, ρ0 = 7.8, 
                              a = 10.0, z= 0.1, z0 = 0.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:120]), fluence_DA_semiinf_CW(2.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-3)

### test 

end # module


data = cylinder_inputs(μsp = 10.0, μa = 0.1, n_ext= 1.0, n_med = 1.0, lz = 2.0, ρ = 1.0, ρ0 = 1.0, 
                              a = 1.0, z= 1.1, z0 = 1.1, pw = 0.0, θ = 0.0, ω = 0.0);
@test isapprox(fluence_DA_cylinder_CW(data, besselj_roots[1:600, 1:120]), fluence_DA_semiinf_CW(2.0, [0.1, 10.0], 1.0, 1.0, 0.1), rtol = 1e-3)
