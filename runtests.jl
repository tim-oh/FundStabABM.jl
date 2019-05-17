using Test

@testset "All tests" begin
    include("test/types_test.jl")
    include("test/params_test.jl")
    include("test/functions_test.jl")
    include("tst.jl")
end
