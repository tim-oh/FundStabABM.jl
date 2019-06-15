using Random

using FundStabABM.Types
using FundStabABM.Params
using FundStabABM.Func


@testset "Market Type" begin

    market = Types.MarketIndex(vcat(
        Params.marketstartval,
        zeros(Params.bigt-1)))

    @test isa(market, Types.MarketIndex)

    @test typeof(market.value) == Array{Float64, 1}

    @test size(market.value) == (Params.bigt, )

    @test market.value[1] == Params.marketstartval
    @test all(market.value .>= 0)
end


@testset "Stocks Type" begin

    stocks = Types.Equity(
        zeros(Params.bigm, Params.bigt),
        zeros(Params.bigm),
        zeros(Params.bigm),
        zeros(Params.bigm),
        zeros(Params.bigm, Params.bigt))

    @test isa(stocks, Types.Equity)

    @test typeof(stocks.value) == Array{Float64, 2}
    @test typeof(stocks.vol) == Array{Float64, 1}
    @test typeof(stocks.beta) == Array{Float64, 1}
    @test typeof(stocks.impact) == Array{Float64, 1}
    @test typeof(stocks.volume) == Array{Float64, 2}

    @test size(stocks.value) == (Params.bigm, Params.bigt)
    @test size(stocks.beta) == (Params.bigm,)
    @test size(stocks.vol) == (Params.bigm,)
    @test size(stocks.impact) == (Params.bigm,)
    @test size(stocks.volume) == (Params.bigm, Params.bigt)

    # NOTE: No test of stock.beta now
    @test all(stocks.value .>= 0)
    @test all(stocks.vol .>= 0)
    @test all(stocks.impact .>= 0)
    @test all(stocks.volume .>= 0)
end

@testset "Fund Type" begin

    funds = Types.EquityFund(
    zeros(Params.bigk, Params.bigm),
    zeros(Params.bigk, Params.bign),
    zeros(Params.bigk, Params.bigt))

    @test isa(funds, Types.EquityFund)

    @test typeof(funds.holdings) == Array{Float64, 2}
    @test typeof(funds.stakes) == Array{Float64, 2}
    @test typeof(funds.value) == Array{Float64, 2}

    @test all(funds.holdings .>= 0)
    @test all(funds.stakes .>= 0)
    @test all(funds.value .>= 0)

    @test size(funds.holdings) == (Params.bigk, Params.bigm)
    @test size(funds.stakes) == (Params.bigk, Params.bign)
    @test size(funds.value) == (Params.bigk, Params.bigt)
end


@testset "Investor Type" begin

    investors = Types.RetailInvestor(
    zeros(Params.bign, Params.bigk + 1),
    zeros(Params.bign),
    zeros(Params.bign))

    @test isa(investors, Types.RetailInvestor)

    @test typeof(investors.assets) == Array{Float64, 2}
    @test typeof(investors.horizon) == Array{Int64, 1}
    @test typeof(investors.threshold) == Array{Float64, 1}

    @test size(investors.assets) == (Params.bign, Params.bigk + 1,)
    @test size(investors.horizon) == (Params.bign,)
    @test size(investors.threshold) == (Params.bign,)

    @test all(investors.assets .>= 0)
    @test all(investors.horizon .>= 0)
    @test all(1 .>= investors.threshold .>= -1)
end

@testset "Order Types" begin

    buyorder = Types.BuyMarketOrder(
    Array{Float64}(undef, 0, Params.bigm),
    Array{Int64}(undef, 0))

    sellorder = Types.SellMarketOrder(
    Array{Float64}(undef, 0, Params.bigm),
    Array{Int64}(undef, 0))

    @test isa(buyorder, Types.BuyMarketOrder)
    @test isa(sellorder, Types.SellMarketOrder)

    @test typeof(buyorder.values) == Array{Float64, 2}
    @test typeof(buyorder.funds) == Array{Int64, 1}
    @test typeof(sellorder.values) == Array{Float64, 2}
    @test typeof(sellorder.investors) == Array{Int64, 1}

    @test size(buyorder.values) == (0, Params.bigm)
    @test size(buyorder.funds) == (0, )
    @test size(sellorder.values) == (0, Params.bigm)
    @test size(sellorder.investors) == (0, )

    @test all(buyorder.values .>= 0)
    @test all(buyorder.funds .>= 0)
    @test all(sellorder.values .<= 0)
    @test all(sellorder.investors .>= 0)
end

#  Market makes perhaps not necessary
# @testset "Market Maker Type" begin
#
#     mktmaker = Types.MarketMaker(zeros(Params.bign, Params.bigm))
#
#     @test isa(mktmaker, Types.MarketMaker)
#
#     @test typeof(mktmaker.orderbook) == Array{Float64, 2}
#
#     @test size(mktmaker.orderbook) == (Params.bign, Params.bigm)
#
#end
