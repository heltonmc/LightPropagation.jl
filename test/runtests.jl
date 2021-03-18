module LightPropagationRunTests

using Test

@time @testset "DAsemiinf" begin include("DAsemiinfTests/runtests.jl") end

@time @testset "DA_Nlayer" begin include("DA_NlayerTests/runtests.jl") end


end