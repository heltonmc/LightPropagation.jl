module LightPropagation

greet() = print("Hello World!")

end # module



function ffplan(x)
       pl = plan_rfft(x; flags=FFTW.ESTIMATE)
end

function ffplan1(x)
       pl = plan_rfft(x; flags=FFTW.MEASURE) #seems to be good choice
end

function ffplan2(x)
       pl = plan_rfft(x; flags=FFTW.PATIENT)
end

function ffplan3(x)
       pl = plan_rfft(x; flags=FFTW.EXHAUSTIVE)
end

@benchmark ffplan(rand(1024))

  median time:      77.641 μs (0.00% GC) (1024)
  median time:      120.184 μs (0.00% GC) (4095)





@benchmark ffplan1(rand(1024))

 median time:      39.838 μs (0.00% GC) (1024)
 median time:      128.125 μs (0.00% GC) (4095)



@benchmark ffplan2(rand(1024))

 median time:      30.988 μs (0.00% GC)
 median time:      125.710 μs (0.00% GC)



@benchmark ffplan3(rand(1024))

  median time:      31.056 μs (0.00% GC)
  median time:      125.198 μs (0.00% GC)


pl = ffplan(rand(4095))
pl1 = ffplan1(rand(4095))
pl2 = ffplan2(rand(4095))
pl3 = ffplan3(rand(4095))




function multpl(pl, u)
       pl * u
end

#identical for 1024

function I_0=ILT(mua,mus,g,r,N)
%Hello!...

%N:   order of approximation
%mua: absorption coefficient
%mus: scattering coefficient
%g :  anisotropic factor
%r:   distance to the isotropic source
%**************************************************************************

c = 2.99792458/1.4*1e11;                      %speed of light
D_k = 1/(3*(mus+mua*1));                      %diffusion coefficient
R_s = 50;                                     %radius of the sphere
l=0:N;
sigma(l+1) = (1-g.^l).*mus./(1-g);            %Henyey-Greenstein function
sigma(1) = 0;
k_end = 400;                                  %number of dircrete wave numbers
t_end = 120*1e-11;                            %time in seconds
dt = 0.5*1e-12;
n_end = round(t_end/dt);                      %number time values

eig_k=zeros(N+1,k_end);
Res=zeros(N+1,k_end);

for k=1:k_end
    ek = k*pi/R_s;
    A = diag(1i*ek*(1:N)./sqrt((1:2:2*N-1).*(3:2:2*N+1)),+1) + diag(sigma(1:N+1)) + diag(1i*ek*(1:N)./sqrt((1:2:2*N-1).*(3:2:2*N+1)),-1);
    [V,D] = eig(A);                             %eigenvalue decomposition
    b = V\eye(N+1,1);                           %first column vector of the inverse matrix to U
    eig_k(1:N+1,k) = diag(D);
    Res(1:N+1,k) = k.*sin(r*ek)*V(1,1:N+1).*transpose(b(1:N+1));        %Eq. 25 of the manusript

end

t =zeros(1,n_end);
I_0 =zeros(1,n_end);
G_DE=zeros(1,n_end);
for u = 1:n_end;
    t(u) = -0.25e-12 + u*dt;
    e_Matrix = exp(-eig_k*c*t(u));
    I_help = Res.*e_Matrix;
    I_0(u) = 1e-9*c/(2*r*R_s^2)*real(sum(I_help(:))).*exp(-mua*c*t(u)).*Heaviside(t(u)-r/c+1*dt);
    G_DE(u) = 1e-9*c*(1-(r./c./t(u)).^2).^0.125./(4*pi*c.*t(u)/3./mus).^1.5.*exp(-(mua+mus).*c.*t(u)).*exp(mus.*c.*t(u).*(1-(r./c./t(u)).^2).^0.75).*sqrt(1+2.026./(mus.*c.*t(u).*(1-(r./c./t(u)).^2).^0.75)).*Heaviside(t(u)-r/c-1*dt);
end
semilogy(1e9*t,I_0*1e-3,'r');









function testcase(x, a, b, c, d, e)

    ϕ = sinh(x*(a + b))/(d*x)
    ϕ = ϕ*BigFloat((d*x*cosh(x*(c - e)) + d*x*sinh(x*(c - e))))
    ϕ = ϕ/BigFloat((d*x*cosh(x*(c + a)) + d*x*sinh(x*(c + a))) - sinh(x*(b - e)))
end



musp = 10
mua = 0.1
z0 = 1/(musp + mua)
D = 1/3musp
zb = 2*D
l = 4.0

    a = zb
    b = z0
    c = l
    d = D
    e = z
