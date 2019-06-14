using Test, FundStabABM, Random, StatsBase, StatsPlots, Distributions
using ProgressMeter

@testset "All tests" begin
    include("./types_test.jl")
    include("./params_test.jl")
    include("./functions_test.jl")
#    include("tst.jl")
end
