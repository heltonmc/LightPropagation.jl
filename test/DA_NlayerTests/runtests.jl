module DA_NlayerTests

using Test

@testset "DA_Nlayer_functionTests" begin
    include("DA_Nlayer_functionTests.jl")
end

@testset "DA_Nlayer_kernelTests" begin
    include("DA_Nlayer_kernelTests.jl")
end

end # module