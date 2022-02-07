module LightPropagationRunTests

using Test

@time @testset "DAsemiinf" begin include("DAsemiinfTests/runtests.jl") end

@time @testset "DAslab" begin include("DAslabTests/runtests.jl") end

@time @testset "DAparalpip" begin include("DAparalpipTests/runtests.jl") end

@time @testset "DA_Nlayer" begin include("DA_NlayerTests/runtests.jl") end

@time @testset "DAinf" begin include("DAinfTests/runtests.jl") end

@time @testset "DAdcs" begin include("DA_DCSTests/runtests.jl") end

end
