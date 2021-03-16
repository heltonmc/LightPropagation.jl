# Light Diffusion in a N-layered turbid cylinder

The salient features and their implementation of [Liemert and Kienle's](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-18-9-9266&id=198273)[^1] solution to the diffusion equation for the N-layered finite cylinder are discussed. 


#### Helpful Background Reading:
- Light diffusion in a turbed cylinder. I. Homogoeneous case[^2]
- Light diffusion in a turbed cylinder. II. Layered case[^1]
- Noninvasive determination of the optical properties of two-layered turbid media[^3]
- Light diffusion in N-layered turbid media: steady-state domain[^4]
- Light diffusion in N-layered turbid media: frequency and time domains[^5]


#### Nomenclature:
We follow the solution for a point beam incident onto the center top of the first layer of N-layered cylinders with radii ``a``. In cylindrical cooridnates ``\vec{r}=(\rho, \phi, z) `` the source position can be defined as ``\vec{r_0}=(0,0,z_0)`` where ``z_0`` is the assumed position of an isotropic source at a distance ``z_0 = \frac{1}{\mu_{s1}'}``. We ignore the absorption term but if desired can be redefined in the `diffusionparams` function. The thickness, refractive index, reduced scattering and absorption coefficients of layer k are denoted by ``l_k, n_k, \mu_{s_k}', \mu_{a_k}``, respectively. The extrapolated boundary condition is used with the extrapolation length ``z_{b_k} = 2 A_k D_k`` where ``A_k`` reflection factor calculated from the angle-averaged probability for reflection at a boundary between layer ``k`` and the surrounding medium. We calculate ``A_k`` from the polonomial approximation listed in [Contini 1997](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-36-19-4587) but can be calculated in anyway by altering the `get_afac` function. Finally, the diffusion coefficient in layer k is ``D_k = 1/(3\mu_s')``.

#### Solution
The solution for the fluence in the first layer ``\Phi_1(\vec{r}, \omega)`` for a point source incident onto the center of the cylinder top can then be described by Equation 22 in Ref 1. [^1]

```math
\Phi_1(\vec{r}) = \frac{1}{\pi a'^2} \sum_{n=1}^{\infty} G_1(s_n, z) J_0(s_n \rho)J_1^{-2}(a' s_n)
```

Where ``a' = a + z_{b_k}``, ``\rho`` is the source-detector separation, and ``J_m`` is the Bessel function of the first kind and order ``m``. More interestingly, ``s_n`` are the roots of ``J_m`` satisfying ``J_m(a' s_n) = 0,\qquad n = 1, 2, ...,``. For an incident source on the center top, the summation over all orders of ``J_m`` is no longer needed so we only need to consider ``m=0``. ``s_n`` then becomes the roots to Bessel function of the first kind and zero order. It should be noted that both ``J_0`` and ``J_1`` appear in the solution but this is because the general solution is in the form of ``J_m`` and ``J_{m+1}``. For those familiar with the trickier case presented in [^3] and then extended in [^4] and [^5] this is main contribution to the accuracy, stability, and speed of the solutions. For those unfamilar, the decisive point is the numerically correct and efficient 2-D inverse Fourier transform of ``\Phi(z, s)`` must be calculated to obtain the fluence in real space ``\Phi(\vec{r})`` (directly taking language from Liemert here). I would say this is most easily implemented by taking advantage of rotational symmetry and using the 1-D inverse Hankel Transform. In contrast, the nodes of the integral are not fixed like in the cylindrical case presented here, but is solved by adaptive Gauss integration. If the nodes were consistent, they could be precomputed and used for each case but the amount of nodes and their location depend on the forward parameters in the model. All this to say, precalculating the roots and storing them is highly efficient and significantly more accurate and stable. The roots are divided by ``a'`` for each new condition.

#### Green's function
Finally, the Green's function in the first layer ``G_1(s_n, z)`` is:
```math
G_1(s_n, z) = \frac{e^{-\alpha_1 |z - z_0|} - e^{-\alpha_1 (z + z_0 + 2z_{b_1})}}{2D_1 \alpha_1} + \frac{\sinh{[\alpha_1(z_0 + z_{b_1})]} \sinh{[\alpha_1(z + z_{b_1})]}}{D_1\alpha_1 e^{\alpha_1(l_1 + z_{b_1})}} \\
\times \frac{D_1\alpha_1n_1^2\beta_3 - D_2\alpha_2n_2^2\gamma_3}{D_1\alpha_1n_1^2\beta_3\cosh{\alpha_1(l_1 + z_{b_1})} + D_2\alpha_2n_2^2\gamma_3\sinh{\alpha_1(l_1 + z_{b_1})}}
```

Where ``\alpha_k`` has the form in the steady-state:
```math
\alpha_k = \sqrt{s_n^2 + \frac{\mu_{a_k}}{D}}
```
Whereas a complex term is added in the Frequency Domain:
```math
\alpha_k = \sqrt{s_n^2 + \frac{\mu_{a_k}}{D} + \frac{i\omega}{Dc}}
```
As you can notice the first term of ``G_1(s_n, z)`` contains exponentially decaying functions which is great as you can see that ``\alpha`` has to be summed for very large numbers. Generally ``z, z_b, z_0`` are around 0.1 so this is easy to handle computationally. On the other hand, terms like ``\sinh{\alpha l}`` and ``\cosh{\alpha l}`` go to infinity very quickly. It is much better to define (though messier) the Green's function in terms of exponetials using ``\sinh{x} = (e^x - e^{-x})/2`` and ``\cosh{x} = (e^x + e^{-x})/2``. After some algebraic manipulation we can then define ``G_1(s_n, z)`` as:
```math
G_1(s_n, z) = \frac{e^{-\alpha_1 |z - z_0|} - e^{-\alpha_1 (z + z_0 + 2z_{b_1})}}{2D_1 \alpha_1} + \frac{e^{\alpha_1(z_0 + z - 2l_1)}[1 - e^{-2\alpha_1(z_0 + z_{b_1})}]}{2D_1\alpha_1} \\
\times \frac{D_1\alpha_1n_1^2\beta_3 - D_2\alpha_2n_2^2\gamma_3}{D_1\alpha_1n_1^2\beta_3[1 + e^{-2\alpha_1(l_1 + z_{b_1})}] + D_2\alpha_2n_2^2\gamma_3[1 - e^{-2\alpha_1(l_1 + z_{b_1})}]}
```
Great. So now we have expressed the Green's function in terms of exponentially decaying functions instead of dividing by exponentially growing functions so we are confident that we won't get overflow errors as ``\alpha`` goes to infinity. We are now ready to discuss the last pesky thing which are the ``\gamma_3`` and ``\beta_3`` factors. These are generally obtained from recurrence relations in [^1] by equation 17 and 18. As presented they also contain exponentially growing functions. Where before we only had ``\sinh{\alpha_1 l_1}`` now we have exponential functions of all the layers. This is especially challenging if we consider the bottom layer to be large. In general the same procedure we used prior should be used but we can't explicitly write them in terms of decaying functions because these factors actually do exponentially increase. What we can do is take advantage of the way they present in ``G_1(s_n, z)`` which is ``(\beta_3 - \gamma_3)/(\beta_3 + \gamma_3)``. We need to find a common exponentially growing factor that presents in both coefficients that we can then cancel.

After we factor a common term out we can then define the coefficients for an N-layered cylinder as :

``N = 2``:
```math
\beta_3 = 1 - e^{-2 \alpha_2 (l_2 + z_{b_2})} \\
\gamma_3 = 1 + e^{-2 \alpha_2 (l_2 + z_{b_2})}
```

``N = 3``:
```math
\beta_3 = D_2 \alpha_2 n_2^2 (1 + e^{-2\alpha_2 l_2})(1 - e^{-2 \alpha_3 (l_3 + z_{b_2})}) + D_3 \alpha_3 n_3^2 (1 - e^{-2\alpha_2 l_2})(1 + e^{-2 \alpha_3 (l_3 + z_{b_2})}) \\
\gamma_3 = D_2 \alpha_2 n_2^2 (1 - e^{-2\alpha_2 l_2})(1 - e^{-2 \alpha_3 (l_3 + z_{b_2})}) + D_3 \alpha_3 n_3^2 (1 + e^{-2\alpha_2 l_2})(1 + e^{-2 \alpha_3 (l_3 + z_{b_2})})
```

``N = 4``:
```math
\beta_4 = D_3 \alpha_3 n_3^2 (1 + e^{-2\alpha_3 l_3}) (1 - e^{-2\alpha_4 (l_4 + z_{b_2}}) + D_4 \alpha_4 n_4^2 (1 - e^{-2\alpha_3 l_3}) (1 + e^{-2\alpha_4 (l_4 + z_{b_2}}) \\
\gamma_4 = D_3 \alpha_3 n_3^2 (1 - e^{-2\alpha_3 l_3}) (1 - e^{-2\alpha_4 (l_4 + z_{b_2}}) + D_4 \alpha_4 n_4^2 (1 + e^{-2\alpha_3 l_3}) (1 + e^{-2\alpha_4 (l_4 + z_{b_2}}) \\
```

```math
\beta_3 = D_2 \alpha_2 n_2^2 \beta_4 (1 + e^{-2\alpha_2 l_2}) + D_3 \alpha_3 n_3^2 \gamma_4 (1 - e^{-2\alpha_2 l_2}) \\
\gamma_3 = D_2 \alpha_2 n_2^2 \beta_4 (1 - e^{-2\alpha_2 l_2}) + D_3 \alpha_3 n_3^2 \gamma_4 (1 + e^{-2\alpha_2 l_2})
```

``N > 4``: Use start values
```math
\beta_N = D_{N-1} \alpha_{N-1} n_{n-1}^2 (1 + e^{-2\alpha_{N-1}l_{N-1}})(1 - e^{-2\alpha_{N-1}(l_{N-1} + z_{b_2})}) + D_{N} \alpha_{N} n_{n}^2 (1 - e^{-2\alpha_{N-1}l_{N-1}})(1 + e^{-2\alpha_{N}(l_{N} + z_{b_2})}) \\
\gamma_N = D_{N-1} \alpha_{N-1} n_{n-1}^2 (1 - e^{-2\alpha_{N-1}l_{N-1}})(1 - e^{-2\alpha_{N-1}(l_{N-1} + z_{b_2})}) + D_{N} \alpha_{N} n_{n}^2 (1 + e^{-2\alpha_{N-1}l_{N-1}})(1 + e^{-2\alpha_{N}(l_{N} + z_{b_2})})
```
With downward recurrence relations: 
```math
\beta_{k - 1} = D_{k - 2} \alpha_{k - 2} n_{k - 2}^2 (1 + e^{-2\alpha_{k - 2}l_{k - 2}}) \beta_k + D_{k - 1} \alpha_{k - 1} n_{k - 1}^2 (1 - e^{-2\alpha_{k - 2}l_{k - 2}}) \gamma_k \\
\gamma_{k - 1} = D_{k - 2} \alpha_{k - 2} n_{k - 2}^2 (1 - e^{-2\alpha_{k - 2}l_{k - 2}}) \beta_k + D_{k - 1} \alpha_{k - 1} n_{k - 1}^2 (1 + e^{-2\alpha_{k - 2}l_{k - 2}}) \gamma_k
```

[^1]: Liemert, André, and Alwin Kienle. "Light diffusion in a turbid cylinder. II. Layered case." Optics Express 18.9 (2010): 9266-9279.
[^2]: Liemert, André, and Alwin Kienle. "Light diffusion in a turbid cylinder. I. Homogeneous case." Optics Express 18.9 (2010): 9456-9473.
[^3]: Kienle, Alwin, et al. "Noninvasive determination of the optical properties of two-layered turbid media." Applied optics 37.4 (1998): 779-791.
[^4]: Liemert, André, and Alwin Kienle. "Light diffusion in N-layered turbid media: steady-state domain." Journal of biomedical optics 15.2 (2010): 025003.
[^5]: Liemert, André, and Alwin Kienle. "Light diffusion in N-layered turbid media: frequency and time domains." Journal of biomedical optics 15.2 (2010): 025002.