using FundStabABM.Params

@testset "Parameters" begin
    @test Params.bign >= Params.bigk >= 1 # No fewer investors than funds
    @test Params.bigt >= 1
    @test -1< Params.drift <1
    @test Params.marketvol >= 0
    @test 0 < Params.perfwindow[1] <= Params.perfwindow[end] < Params.bigt
    @test Params.marketstartval > 0
    # Positive starting value required given the multiplicative formula
    @test Params.marketvol >= 0
    @test Params.betamean > 0 # Average stock positively correlated with index
    @test Params.betastd >= 0
    @test 0 <= Params.stockvolrange[1] <= Params.stockvolrange[end]
    @test 0 <=Params.invcaprange[1] <= Params.invcaprange[2]
    @test Params.stockstartval >= 0
    @test -1 < Params.thresholdmean < 1
    @test Params.thresholdstd >= 0
    @test 1 <=Params.portfsizerange[1]<=Params.portfsizerange[end] <=Params.bigm
    @test 0 <= Params.impactrange[1] <= Params.impactrange[2]

end
