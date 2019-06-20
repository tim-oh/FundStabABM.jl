using FundStabABM.Params

@testset "Parameters" begin
    @test Params.default().bign >= Params.default().bigk >= 1
    @test Params.default().bigt >= 1
    @test -1< Params.default().drift <1
    @test Params.default().marketvol >= 0
    @test 0 < Params.default().perfwindow[1] <=
        Params.default().perfwindow[end] < Params.default().bigt
    @test Params.default().marketstartval > 0
    @test Params.default().marketvol >= 0
    @test Params.default().betamean > 0
    @test Params.default().betastd >= 0
    @test 0 <= Params.default().stockvolrange[1] <=
        Params.default().stockvolrange[end]
    @test 0 <=Params.default().invcaprange[1] <= Params.default().invcaprange[2]
    @test Params.default().stockstartval >= 0
    @test -1 < Params.default().thresholdmean < 1
    @test Params.default().thresholdstd >= 0
    @test 1 <=Params.default().portfsizerange[1] <=
        Params.default().portfsizerange[end] <=Params.default().bigm
    @test 0 <= Params.default().impactrange[1] <=Params.default().impactrange[2]

end
