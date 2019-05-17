using Test, Random
include("../src/types.jl")
include("../src/params.jl")
include("../src/functions.jl")

using .Types
using .Params
using .Func


@testset "Market Type" begin
    market = Types.MarketIndex(vcat(
        Params.marketstartval,
        zeros(Params.big_t-1)))

    @test isa(market, Types.MarketIndex)

    @test typeof(market.value) == Array{Float64, 1}

    @test size(market.value) == (Params.big_t, )

    @test market.value[1] == Params.marketstartval
    @test all(market.value .>= 0)

    @test all(Func.marketinit!(market.value, Params.perfwindow[2]) .>= 0)
end


@testset "Stocks Type" begin
    stocks = Types.Equity(
        zeros(Params.big_m, Params.big_t),
        zeros(Params.big_m),
        zeros(Params.big_m))

    @test isa(stocks, Types.Equity)

    @test typeof(stocks.value) == Array{Float64, 2}
    @test typeof(stocks.vol) == Array{Float64, 1}
    @test typeof(stocks.beta) == Array{Float64, 1}

    @test size(stocks.value) == (Params.big_m, Params.big_t)
    @test size(stocks.beta) == (Params.big_m,)
    @test size(stocks.vol) == (Params.big_m,)

    @test all(stocks.value .>= 0)
    @test all(stocks.vol .>= 0)
    @test all(stocks.beta .>= 0)
    @test_broken all(stocks.impact .>= 0)
end

@testset "Fund Type" begin
    funds = Types.EquityFund(
    zeros(Params.big_k, Params.big_m),
    zeros(Params.big_k, Params.big_n),
    zeros(Params.big_k, Params.big_t))

    @test isa(funds, Types.EquityFund)

    @test typeof(funds.holdings) == Array{Float64, 2}
    @test typeof(funds.stakes) == Array{Float64, 2}
    @test typeof(funds.value) == Array{Float64, 2}

    @test all(funds.holdings .>= 0)
    @test all(funds.stakes .>= 0)
    @test all(funds.value .>= 0)

    @test size(funds.holdings) == (Params.big_k, Params.big_m)
    @test size(funds.stakes) == (Params.big_k, Params.big_n)
    @test size(funds.value) == (Params.big_k, Params.big_t)
end


@testset "Investor Type" begin
    investors = Types.RetailInvestor(
    zeros(Params.big_n, Params.big_k + 1),
    zeros(Params.big_n),
    zeros(Params.big_n))

    @test isa(investors, Types.RetailInvestor)

    @test typeof(investors.assets) == Array{Float64, 2}
    @test typeof(investors.horizon) == Array{Float64, 1}
    @test typeof(investors.threshold) == Array{Float64, 1}

    @test size(investors.assets) == (Params.big_n, Params.big_k + 1,)
    @test size(investors.horizon) == (Params.big_n,)
    @test size(investors.threshold) == (Params.big_n,)

    @test all(investors.assets .>= 0)
    @test all(investors.horizon .>= 0)
    @test all(investors.threshold .>= 0)
end

@testset "Market Maker Type" begin
    @test_broken isa(mktmaker, MarketMaker)
end
