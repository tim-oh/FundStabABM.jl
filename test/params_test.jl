include("../src/params.jl")
using Test, .Params

@testset "Parameters" begin
    @test Params.bign >= Params.bigk >= 1 # No fewer investors than funds
    @test Params.bigt >= 1
    @test -1< Params.drift <1
    @test Params.marketvol >= 0
    @test Params.perfwindow[1] <= Params.perfwindow[2]
    @test Params.marketstartval > 0
    # Positive starting value required given the multiplicative formula
    @test Params.marketvol >= 0
    @test Params.perfwindow[1] <= Params.perfwindow[2]
    @test Params.betamean > 0 # The average stock should correlate positively with the index
    @test Params.betastd >= 0
    @test Params.stockvolrange[2] >= Params.stockvolrange[1] >= 0
    @test Params.invcaprange[2] >= Params.invcaprange[1] >= 0
    @test Params.stockstartval >= 0
    @test Params.horizonrange[1] <= Params.horizonrange[2] < Params.bigt
    @test -1 < Params.thresholdmean < 1
    @test Params.thresholdstd >= 0
    @test 1 <= Params.portfsizerange[1] <=Params.portfsizerange[end] <=Params.bigm

    ## TODO: @test Params.stockimpactrange

end
