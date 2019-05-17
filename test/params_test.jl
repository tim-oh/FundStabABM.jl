include("../src/params.jl")
using Test, .Params

@testset "Parameters" begin
    @test Params.big_n >= Params.big_k >= 1 # No fewer investors than funds
    @test Params.big_t >= 1
    @test -1< Params.drift <1
    @test Params.marketvol >= 0
    @test Params.perfwindow[1] <= Params.perfwindow[2]
    @test Params.marketstartval > 0 # Positive starting value required given the multiplicative formula
    @test Params.marketvol >= 0
    @test Params.perfwindow[1] <= Params.perfwindow[2]
    @test Params.betamean > 0 # The average stock should correlate positively with the index
    @test Params.betastd >= 0
    @test Params.stockvolrange[2] >= Params.stockvolrange[1] >= 0
    @test Params.capitalrange[2] >= Params.capitalrange[1] >= 0

    ## TODO: @test Params.stockimpactrange

end
