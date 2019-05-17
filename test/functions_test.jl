using Test, Random, LinearAlgebra

include("../src/types.jl")
include("../src/params.jl")

using .Types

#
const big_m = 10 # Number of stocks
const big_n = 10 # Number of investors
const big_t = 20 # Number of time periods
const big_k = 5 # Number of funds
const marketstartval = 100 # Market index starting value
const drift = 0.05 # Market index drift
const marketvol = 0.1 # Market volatility (std?)
const perfwindow = (1, 3) # Performance window for investor
const betamean = 1 # Average stock beta
const betastd = 0.3 # Stock beta dispersion
const stockvolrange = (0, 0.2)
const capitalrange = (10, 10) # Investors' range of initial capital



@testset "Initialisation Functions" begin

    # TODO: Set seed. Don't the stocks have to be created first?
    @test Func.betainit()[1] == randn(MersenneTwister(1)) * Params.betastd + Params.betamean
    # Should I write the function so that it CAN take more arguments?

    # TODO: Set seed, remove [1]
    @test Func.stockvolinit()[1] == Params.stockvolrange[1] + rand(MersenneTwister(1)) * (Params.stockvolrange[2] - Params.stockvolrange[1])

    market = Types.MarketIndex(vcat(
        Params.marketstartval,
        zeros(Params.big_t-1)))
    market.value[2] = 102
    market.value[3] = 104
    stocks = Types.Equity(
        zeros(Params.big_m, Params.big_t),
        zeros(Params.big_m),
        zeros(Params.big_m))

    stocks.vol .= Func.stockvolinit()
    # TODO: Set seed

    @test Func.stockvalueinit!(stocks.value, Params.perfwindow[2], market.value, stocks.beta, stocks.vol)[1,2] == ((1 .+ (market.value[2] - market.value[1]) / market.value[1] .* stocks.beta) .* stocks.value[2] .+ stocks.vol .* randn(MersenneTwister(1)) .* stocks.value[2])[1]
    # ENTER VALUES INSTEAD


    @test_broken sum(funds.stakes == 1)

    ## TODO: @test stockimpactinit

    investors = Types.RetailInvestor(
    zeros(Params.big_n, Params.big_k + 1),
    zeros(Params.big_n),
    zeros(Params.big_n))

    # TODO: @test Func.invfundinit!(investors.assets) so that assets are non-zero on the diagonal for first k entries
    # TODO: test that last column is empty
    # TODO: test that each row has excatly one non-zero item
    # TODO: test that each column has at least one non-zero item
    # TODO: Change function so that the capital is generated as a vector. Then check that the row sum of the assets array is equal to this vector.
    # QUESTION: how do you have separate tests for A) directly after initialisation, and B) following periods. I really want three kinds (pre-initialisation, post-initialisation, during the runs)

    ## TODO: @test invhorizoninit
    ## TODO: @test invthreshinit

    ## TODO: @test fundholdinit
    ## TODO: @test fundstakeinit
    ## TODO: @test fundvalinit
end

@testset "Price Functions" begin

    @test Func.marketmove(100, MersenneTwister(1999), 0.05, 0.1) == 100 * (1 + 0.05) + randn(MersenneTwister(1999)) * 0.1

    market = Types.MarketIndex(vcat(
        Params.marketstartval,
        zeros(Params.big_t-1)))
    market.value[2] = 102

    stocks = Types.Equity(
        zeros(Params.big_m, Params.big_t),
        zeros(Params.big_m),
        zeros(Params.big_m))


    @test Func.stockmove(2, market.value, stocks.value[2], stocks.beta, stocks.vol, MersenneTwister(2000))[1] == (stocks.value[1] * (1 + ((market.value[2] - market.value[1]) / market.value[1]) * stocks.beta[1])) + stocks.vol[1] * randn(MersenneTwister(2000)) * stocks.value[1]

end
