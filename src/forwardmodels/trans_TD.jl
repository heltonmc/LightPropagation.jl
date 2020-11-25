
function DT_transslab(t, β::Array{Float64,1}, ρ::Float64,ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness
    n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10


	z1m = Array{Float64}(undef, length(xs))
	z2m = Array{Float64}(undef, length(xs))


    Rt1 = Array{Float64}(undef, length(t))
    Rt2 = Array{Float64}(undef, length(t), length(z1m))
	Rt = Array{Float64}(undef, length(t))
    
	if n == 1
		A = 1
	elseif n > 1
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 + 
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end


    zs::Float64 = 1/μsp
    ze::Float64 = 2A*D

    z1m = Float64[(1 .- 2m).*s .- 4m.*ze .- zs for m in xs]
    z2m = Float64[(1 .- 2m).*s .- (4m .- 2).*ze .+ zs for m in xs]
    
    Rt1 = @. exp(-(ρ^2/(4D*ν*t)) - μa*ν*t)
    Rt1 = @. Rt1/(2*(4π*D*ν)^(3/2)*t^(5/2))

    Rt2 = @. z1m'*exp(-(z1m'^2 / (4D*ν*t))) - z2m'*exp(-(z2m'^2 / (4D*ν*t)))
    
    Rt = Rt1.*sum(Rt2, dims=2)
	replace!(Rt, NaN => 0)
	
    return Rt
end





function DT_trans_paralpip(t, β::Array{Float64,1}, ndet::Float64, nmed::Float64, rd::Array{Float64,1}, rs::Array{Float64,1}, L::Array{Float64,1})
	# rd is location of detector [x,y]
	# rs is location of source [x,y]
	# L is dimensions of parallelpiped [Lx, Ly, Lz]

	n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	x::Float64 = rd[1]
	y::Float64 = rd[2]
	z::Float64 = L[3]
	xu::Float64 = rs[1]
	yu::Float64 = rs[2]
	lx::Float64 = L[1]
	ly::Float64 = L[2]
	lz::Float64 = L[3]

	x1l = Array{Float64}(undef, length(xs))
	x2l = Array{Float64}(undef, length(xs))
	y1m = Array{Float64}(undef, length(xs))
	y2m = Array{Float64}(undef, length(xs))
	z1n = Array{Float64}(undef, length(xs))
	z2n = Array{Float64}(undef, length(xs))

	Rt1 = Array{Float64}(undef, length(t))
	Rt2 = Array{Float64}(undef, length(t), length(xs))
	Rt3 = Array{Float64}(undef, length(t), length(xs))
	Rt4 = Array{Float64}(undef, length(t), length(xs))

	Rt = Array{Float64}(undef, length(t))

	if (nmed == ndet)
		A = 1
	elseif nmed > ndet
		costhetac = sqrt(1 - (ndet/nmed).^2)
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = (2/(1-R0) -1 + abs(costhetac.^3))./(1-abs(costhetac.^2))
	else nmed < ndet
		R0 = ((nmed/ndet-1)/(nmed/ndet +1)).^2
		A = 2/(1-R0) - 1
	end

	zo::Float64 = 1/μsp
	zb::Float64 = 2A*D


	x1l = Float64[2l*lx + 4l*zb + xu for l in xs]
	x2l = Float64[2l*lx + (4l-2)*zb - xu for l in xs]
	y1m = Float64[2m*ly + 4m*zb + yu for m in xs]
	y2m = Float64[2m*ly + (4m-2)*zb - yu for m in xs]
	z1n = Float64[2n*lz + 4n*zb + zo for n in xs]
	z2n = Float64[2n*lz + (4n-2)*zb - zo for n in xs]


	Rt1 = @. exp(-μa*ν*t)/(2*(4π*D*ν)^(3/2)*t^(5/2))

	Rt2 = @. exp(-(x-x1l')^2 / (4D*ν*t)) - exp(-(x-x2l')^2 / (4D*ν*t))

	Rt3 = @. exp(-(y-y1m')^2 / (4D*ν*t)) - exp(-(y-y2m')^2 / (4D*ν*t))

	Rt4 = @. (lz - z1n')*exp(-(lz-z1n')^2 /(4D*ν*t)) - (lz - z2n')*exp(-(lz-z2n')^2 /(4D*ν*t))

	Rt = Rt1.*sum(Rt2, dims=2).*sum(Rt3, dims=2).*sum(Rt4, dims=2)
	replace!(Rt, NaN =>0)
	return Rt
end

