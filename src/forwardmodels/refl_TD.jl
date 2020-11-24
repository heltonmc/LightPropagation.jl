function refl_DT(t,β)
    μa = β[1]
    μsp = β[2]
    n = 1.35
    ρ = 1
    c = 29.9792458 # Speed of light cm/ns
    v = c/n # Speed of light in medium
    D = 1/(3μsp) # Diffusion Coefficient
    Rt = @. (v/(4π*D*v*t)^1.5)*exp(-(ρ^2)/(4D*v*t)-μa*v*t)
    replace!(Rt, NaN => 0)
    m = findmax(Rt)
    Rt = Rt./m[1]

end

function refl_DT1(t,β,ρ)

	nmed = 1.35
	ndet = 1.45
	μa = β[1]
	μsp = β[2]
	#ρ = β[3]


	if (nmed == ndet)
		Afac = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else #nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		Afac = 2/(1-R0) - 1
	end

	z0 = 1/(μa+μsp)
	D = 1/(3*(μa + μsp))
	zb = 2*Afac*D

	v = 29.9792345/nmed # speed of light in medium

	Rt1 = @. v*exp(-μa*v*t)
	Rt1 = @. Rt1/((4*π*D*v*t)^1.5)
	Rt1 = @. Rt1*(exp(-(z0^2 + ρ^2)/(4*D*v*t)) - exp(-((2*zb + z0)^2 + ρ^2)/(4*D*v*t)))

	Rt2 = @. 3*exp(-μa*v*t)/(2*((4*π*D*v)^1.5)*(t^2.5))
	Rt2 = @. Rt2*(z0*exp(-((z0^2 + ρ^2)/(4*D*v*t))) + (2*zb + z0)*exp(-((2*zb+z0)^2 + ρ^2)/(4*D*v*t)))

	Rt = @. abs(Rt1)+abs(Rt2)
	replace!(Rt, NaN => 0)
	m = findmax(Rt)
	#Rt = Rt./m[1]
	return Rt
end



function DT_semiinfinite(t, β, ρ)

	nmed = 1.4
	ndet = 1.4
	n = nmed/ndet

    μₐ = β[1]
    μₛₚ = β[2]
    # s = slab thickness
	if n == 1
		A = 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end
	

    D = 1/3μₛₚ
    zₛ = 1/μₛₚ
    zₑ = 2A*D

    ν = 29.9792345/nmed


    z₃ₘ = - zₛ
    z₄ₘ = 2zₑ +zₛ



    Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μₐ*ν*t)
    Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

    Rt2 = @. z₃ₘ*exp(-(z₃ₘ^2/(4D*ν*t))) - z₄ₘ*exp(-(z₄ₘ^2/(4D*ν*t)))

    Rt = @. Rt1*Rt2
    

    replace!(Rt, NaN => 0)
    return Rt

end



   
function DT_reflslab(t, β, ρ, s)

	nmed = 1.4
	ndet = 1.4
	n = nmed/ndet

    μₐ = β[1]
    μₛₚ = β[2]
    # s = slab thickness
	if n == 1
		A = 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end
	

    D = 1/3μₛₚ
    zₛ = 1/μₛₚ
    zₑ = 2A*D

    ν = 29.9792345/nmed

    m = Vector(-10:1:10)
    z₃ₘ = -2m.*s .- 4m.*zₑ .- zₛ
    z₄ₘ = -2m.*s .- (4m .- 2).*zₑ .+ zₛ
    
    Rt1 = @. -exp(-(ρ^2/(4D*ν*t)) - μₐ*ν*t)
    Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))


	#Rt2 = Array{Float64, 2}(undef, length(t), length(z₃ₘ))
	Rt2 = @. z₃ₘ'*exp(-(z₃ₘ'^2 / (4D*ν*t))) - z₄ₘ'*exp(-(z₄ₘ'^2 / (4D*ν*t)))


    Rt21 = sum(Rt2, dims=2)

    

    Rt = @. Rt1*Rt21
    

	replace!(Rt, NaN => 0)

	
	return Rt

end


function DT_transslab(t, β, ρ, s)

    nmed = 1.4
	ndet = 1.4
	n = nmed/ndet

    μₐ = β[1]
    μₛₚ = β[2]
    # s = slab thickness
	if n == 1
		A = 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end


    D = 1/3μₛₚ
    zₛ = 1/μₛₚ
    zₑ = 2A*D

	ν = 29.9792345/nmed

    m = Vector(-10:1:10)
    z₁ₘ = (1 .- 2m).*s .- 4m.*zₑ .- zₛ
    z₂ₘ = (1 .- 2m).*s .- (4m .- 2).*zₑ .+ zₛ
    


    Rt1 = @. exp(-(ρ^2/(4D*ν*t)) - μₐ*ν*t)
    Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

    
    Rt2 = @. z₁ₘ'*exp(-(z₁ₘ'^2 / (4D*ν*t))) - z₂ₘ'*exp(-(z₂ₘ'^2 / (4D*ν*t)))
    Rt21 = sum(Rt2, dims=2)
    
    

    Rt = @. Rt1*Rt21
    

    replace!(Rt, NaN => 0)
    return Rt
end

A1(n) = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3 + 1.171382n^4 # for n<1
A2(n) = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7 #for n>1

# seems to work for range 0.5<n<1.8 