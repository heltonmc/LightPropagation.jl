##### Steady State Transmittance from a Slab #####
# Taken from Eqn. 46 from Contini 1997

function DA_transslabCW(ρ::Float64, β::Array{Float64,1},ndet::Float64, nmed::Float64, s::Float64)
	#s is slab thickness

	n::Float64 = nmed/ndet
    μa::Float64 = β[1]
    μsp::Float64 = β[2]
    D::Float64 = 1/3μsp
	ν::Float64 = 29.9792345/nmed
	xs::UnitRange{Int64} = -10:10

	z1m = Array{Float64}(undef, length(xs))
	z2m = Array{Float64}(undef, length(xs))


    Rt1 = Array{Float64}(undef, length(xs))
    Rt2 = Array{Float64}(undef, length(xs))

	if n == 1.0
		A= 1.0
	elseif n > 1.0
		A = 504.332889 - 2641.00214n + 5923.699064n^2 - 7376.355814n^3 +
		 5507.53041n^4 - 2463.357945n^5 + 610.956547n^6 - 64.8047n^7
	else 
		A = 3.084635 - 6.531194n + 8.357854n^2 - 5.082751n^3
	end

	zs::Float64 = 1/μsp
    ze::Float64 = 2A*D

    z1m = Float64[(1 .- 2m).*s .- 4m.*ze .- zs for m in xs]
    z2m = Float64[(1 .- 2m).*s .- (4m .- 2).*ze .+ zs for m in xs]
    
    
    R1 = @. (z1m*(ρ^2 + z1m^2)^(-3/2))*(1+(μa*(ρ^2 +z1m^2))^(1/2))*exp(-(μa*(ρ^2 +z1m^2)/D)^(1/2))
    R2 = @. (z2m*(ρ^2 + z2m^2)^(-3/2))*(1+(μa*(ρ^2 +z2m^2))^(1/2))*exp(-(μa*(ρ^2 +z2m^2)/D)^(1/2))

    Rt = 1/(4π)*(sum(R1) - sum(R2))
    
	return Rt
end
