using LightPropagation, MAT

#=
Available mua values from MC sims:
μa = 0.001, 0.025500, 0.05, 0.075, 0.18125, 0.2875, 0.39375, 0.5

μsp = 3.0, 6.0, 9.0, 12.0, 15.0, 18.0

sds_cm = 0.025  0.15  0.35  0.55  0.75  0.95  1.15  1.35  1.55  1.75  1.95
=#

#MC sims have nmed = 1.33, ndet = 1.45

function getdata(μa, μsp, ρ)

    #open Monte-Carlo sims from .MAT file
    file = matopen("MCdata.mat")
    dset = read(file)
    dset = dset["MCdata"]

    #find indices that much mua, musp given
    muaind = findall(x -> x == μa, round.(dset["mua"], digits=5))
    muspind = findall(x -> x == μsp, dset["musp"])

    #find unique index for mua, musp pair and sds
    ind = intersect(muaind, muspind)[1][2]
    sdsind = findall(x -> x == ρ, dset["sds_cm"][1])[1][2]

    #extract appropriate counts and time
    Rt = dset["Rt"][ind][:,sdsind]
    t  = 0:0.01:3.0

    return t, Rt
end

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

# model inputs
μa = 0.18125
μsp = 9.0
ρ = 1.15
z = 0.0
ndet = 1.45
nmed = 1.33

A = get_afac(nmed/ndet)

t, RtMC = getdata(μa, μsp, ρ)

RtDA = fluence_DA_semiinf_TD(t, [μa,μsp], ρ, ndet, nmed, z)./(2A)

ind = RtDA .> 0

scatter(t[ind], RtMC[ind], ms = 3, label="MC")
plot!(t[ind], RtDA[ind], lw=3, yaxis=:log10, label="DA")


# Let's double check we can fit our MC data with DT 
## setup inverse structure

ydata = RtMC # data to be fit to is Monte-Carlo sim

#setup IRF as delta function length of ydata

IRF = [1; zeros(length(ydata) - 1)]


# Define a structure called data with the specified fields 
data = fitDTOF(t = t, DTOF = ydata./maximum(ydata), IRF = IRF, ρ = ρ, nmed = nmed, ndet = ndet)

# Define the fitting parameters called m for parallepiped geometry
m = DTOF_fitparams(model = fluence_DA_semiinf_TD, risefactor = 0.9, tailfactor = 1e-5,initparams = [0.1, 10.0])


f, tfit, yraw, yfit  = getfit(data, m)

#scatter(tfit, yraw, m="o", ms = 8, label = "ydata", c="black")
#plot!(tfit, yfit, lw = 3, alpha = 1, label = "DT model", c="red")

f.param

###make plot of figure and residuals
using CairoMakie, ColorSchemes, Colors

fig = Figure(resolution = (600, 400))
ax1 = Axis(fig, xlabel = "time (ns)", ylabel = "log(R(ρ, t))",
    xgridvisible = false, 
    ygridvisible = false
)

ax2 = Axis(fig, bbox = BBox(400, 540, 260, 350),
    xticklabelsize = 12, yticklabelsize = 12,showgrid=false,

)


pts1 = CairoMakie.scatter!(ax1, tfit, yraw, color =:black, markersize = 8)
line1 = lines!(ax1, tfit, yfit, color = :red, linewidth = 2)
limits!(ax1, -0.05, 2.1, -13, 1)

# inset
lines!(ax2, tfit, exp.(yraw) .- exp.(yfit), color = :red)
limits!(ax2, -0.1,2.1,-0.05,0.06)
ax2.yticks = [-0.05,0,0.05]
ax2.xticks = [0,1,2]

fig[1, 1] = ax1
#save("test.png", fig, px_per_unit = 1)

#### calculate the mua, musp, and chi2 for different fit ranges governed by rise and tail factor

function fitrange_error(data)

    risefact = [0.01, 0.1, 0.5, 0.9]
    tailfact = [1e-7, 1e-5, 1e-3, 1e-1]
    mua_mat = zeros(length(risefact), length(tailfact))
    musp_mat = similar(mua_mat)
    χ2 = similar(mua_mat)


    for l in eachindex(risefact)
        for n in eachindex(tailfact)

            m = DTOF_fitparams(model = fluence_DA_semiinf_TD, risefactor = risefact[l], tailfactor = tailfact[n], initparams = [0.16, 10.0])
            f, tfit, yraw, yfit  = getfit(data, m)

            mua_mat[l, n] = f.param[1]
            musp_mat[l, n] = f.param[2] 

            χ2[l, n] = sum((exp.(yraw) .- exp.(yfit)).^2 ./ exp.(yraw)) # outputted fit values are log
        end
    end

    return mua_mat, musp_mat, χ2
end

#plot the heat maps

using Plots

mua_mat, musp_mat, χ2 = fitrange_error(data)

muaerr = (mua_mat .- μa)./μa.*100
musperr = (musp_mat .- μsp)./μsp.*100

Plots.heatmap(log.(abs.(muaerr)), c=:thermal, cbar=:false)
Plots.heatmap(log.(abs.(musperr)), c=:thermal, cbar=:false)
Plots.heatmap(χ2, c=:thermal)




##### effect of error in the IRF


ydata = RtMC # data to be fit to is Monte-Carlo sim

using DSP

function uncertainty_IRF(ydata)

    #setup IRF that peaks at index 20
    IRFind = 20
    IRF = zeros(length(ydata))
    IRF[20] = 1

    #preallocate IRF shift to use in loop and optical tuple results
    IRFshift = similar(IRF)
    mua = []
    musp = []

    # convolve ydata with IRF
    convR = conv(IRF, ydata)[1:length(ydata)]

    m = DTOF_fitparams(model = fluence_DA_semiinf_TD, risefactor = 0.9, tailfactor = 1e-5,initparams = [0.1, 10.0])


    for indshift in IRFind-15:IRFind+9
        IRFshift = zeros(length(ydata))
        IRFshift[indshift] = 1

        data = fitDTOF(t = t, DTOF = convR./maximum(convR), IRF = IRFshift, ρ = ρ, nmed = nmed, ndet = ndet)
        f, tfit, yraw, yfit  = getfit(data, m)
        push!(mua, f.param[1])
        push!(musp, f.param[2])

    end

    return mua, musp
end


mua1, musp1= uncertainty_IRF(ydata)
