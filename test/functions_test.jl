using Test, Random, LinearAlgebra

include("../src/types.jl")
using .Types

# Small parameter values for testing, like those in /src/.jl
const bigk = 3 # Number of funds, 3
const bign = 4 # Number of investors, 4
const bigm = 5 # Number of stocks, 5
const bigt = 6 # Number of time periods, 6
const marketstartval = 100 # Market index starting value, 100
const drift = 0.05 # Market index drift, 0.05
const marketvol = 0.1 # Market volatility (std?), 0.1
const perfwindow = (1, 3) # Performance window for investor (1,3)
const betamean = 1 # Average stock beta, 1
const betastd = 0.3 # Stock beta dispersion, 0.3
const stockstartval = 100
const stockvolrange = (0.1, 0.5) # Range of stock volatility (0.1, 0.5)
const invcaprange = (10, 1000) # Investors' range of initial capital, (10, 1000)
const horizonrange = (1,5)
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds
const portfsizerange = 1:5 # Range of number of stocks in fund portfolio (1,5)
rng = MersenneTwister(1)



@testset "Initialisation Functions" begin

    # Generate market value history
    market = Types.MarketIndex(
        zeros(bigt))
    @test all(Func.marketinit!(market.value, marketstartval, perfwindow[2])
    .>= 0)

    # Generate stock betas
    Random.seed!(1)
    @test Func.betainit(bigm, betastd, betamean) ==
    1 .+ 0.3 .* randn(MersenneTwister(1), 5)

    # Generate stock volatilities
    Random.seed!(2)
    @test Func.stockvolinit(bigm, stockvolrange) ==
    0.1 .+ (0.5 - 0.1) .* rand(MersenneTwister(2), 5)

    # Generate stocks' value history
    market.value[1] = 100
    market.value[2] = 101
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm))
    stocks.beta .= [0, 1, 2, 3, 4]
    stocks.vol .= [0.01, 0.01, 0.01, 0.01, 0.01]
    Random.seed!(3)
    @test Func.stockvalueinit!(stocks, stockstartval, perfwindow[2],
    market.value)[:,1:2] ==
    hcat(ones(5)*100,
    [100, 101, 102, 103, 104] +
    [1.1915557734285787, -2.5197330871745263, 2.0748097755419757,
    -0.9732501200602984, -0.1016067940589428])

    @test_broken sum(funds.stakes == 1)

    #TODO: stocks integration tests

    # TODO: @test stockimpactinit
    # TODO: funds integration tests

    investors = Types.RetailInvestor(
    zeros(bign, bigk + 1),
    zeros(bign),
    zeros(bign))

    Random.seed!(5)
    @test Func.invhorizoninit!(investors.horizon, horizonrange) ==
    [5, 2, 4, 4]

    Random.seed!(6)
    @test Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd) ==
    [-0.06570018342662688, 0.00665192247499981,
    0.06774218278688489, 0.002180862789142308]


    # The first bigk investors put their money in that of the first bigk funds that matches their own index
    Random.seed!(7)
    @test LinearAlgebra.diag(Func.invassetinit!(investors.assets,
    invcaprange, bigk)[1:bigk, 1:bigk]) == [318, 111, 692]
    @test all(investors.assets[:, end] .== 0)
    # Sum of investment matches equals implicit draw of capital vector
    @test sum(investors.assets) == 2013
    # One investment per investor
    @test sum(investors.assets .> 0, dims=2) == ones(Int64, (bign,1))
    # At least one investment per fund
    @test all(sum(investors.assets[:,1:end-1] .> 0, dims=1) .>= 1)

    #TODO: investors integration tests

    # QUESTION: how do you have separate tests for A) directly after initialisation, and B) following periods. I really want three kinds (pre-initialisation, post-initialisation, during the runs)

    funds = Types.EquityFund(
    zeros(bigk, bigm),
    zeros(bigk, bign),
    zeros(bigk, bigt))

    @test Func.fundvalinit!(funds.value, investors.assets) ==
    hcat([318.0, 1003.0, 692.0], zeros(3,5))

    @test Func.fundstakeinit!(funds.stakes, investors.assets) ==
    [1.0 0.0 0.0 0.0;
    0.0 0.1106679960119641 0.0 0.8893320039880359;
    0.0 0.0 1.0 0.0]
    @test sum(funds.stakes, dims=2) == ones(bigk,1)

    Random.seed!(8)
    @test Func.fundholdinit!(
    funds.holdings, portfsizerange, funds.value[:, 1], stocks.value[:, 1]) ≈
    [0 0 0 318/100 0; 1003/500 0 3009/500 1003/500 0;
    692/200 0 0 692/200 0]

    # TODO: @test funds.value[:,1] == funds.holdings .* stocks.prices[1]

end

@testset "Price Functions" begin

    Random.seed!(4)
    @test Func.marketmove(100, 0.05, 0.1) ==
    100 * (1 + 0.05) + 0.1 * randn(MersenneTwister(4))

    market = Types.MarketIndex(vcat(
        marketstartval,
        zeros(bigt-1)))
    market.value[2] = 102

    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm))


    @test Func.stockmove(2, market.value, stocks.value[2], stocks.beta, stocks.vol)[1] == (stocks.value[1] * (1 + ((market.value[2] - market.value[1]) / market.value[1]) * stocks.beta[1])) + stocks.vol[1] * randn(MersenneTwister(2000)) * stocks.value[1]

end
