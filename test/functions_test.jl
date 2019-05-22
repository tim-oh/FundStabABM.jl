using Test, Random, LinearAlgebra

# TODO: Use Reexport to hopefully get rid of WARNIN: replacing module ...
include("../src/types.jl")
using .Types

# Small parameter values for testing, like those in /src/.jl
const bigk = 3 # Number of funds, 3
const bign = 4 # Number of investors, 4
const bigm = 5 # Number of stocks, 5
const bigt = 6 # Number of time periods, 6
const mktstartval = 100 # Market index starting value, 100
const drift = 0.05 # Market index drift, 0.05
const marketvol = 0.1 # Market volatility (std?), 0.1
const perfwindow = 1:3 # Performance window for investor (1,3)
# QUESTION: should this be a continuous range or is discrete okay?
const betamean = 1 # Average stock beta, 1
const betastd = 0.3 # Stock beta dispersion, 0.3
const stockstartval = 100 # Starting value of stocks
# TODO: turn into range
const stockvolrange = range(0.01, stop=0.1, step=0.01)
# Range of stock volatility (0.1, 0.5). Consider if it should be continuous.
const invcaprange = (10,1000) # Investors' range of initial capital, (10, 1000)
# TODO: turn into proper range, though requires major recalculation OR
# different specs of the tests so that they're independent
const horizonrange = 1:5
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds
const portfsizerange = 1:5 # Range of number of stocks in fund portfolio (1,5)
const impactrange = 0.001:0.001:0.01 # Stock price impact per currency unit
rng = MersenneTwister(1)


@testset "Initialisation Functions" begin

    # Generate market value history
    market = Types.MarketIndex(
        zeros(bigt))

    Random.seed!(0)
    @test all(Func.marketinit!(market.value, mktstartval, perfwindow[end],
    drift, marketvol)
    .>= 0)
    # TODO: Replace inequality with specific values

    # TODO: stocks integration tests

    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm))

    # Generate stock betas
    Random.seed!(1)
    @test Func.betainit!(stocks.beta, bigm, betastd, betamean) ==
    1 .+ 0.3 .* randn(MersenneTwister(1), 5)

    # Generate stock volatilities
    Random.seed!(2)
    @test Func.stockvolinit!(stocks.vol, stockvolrange, bigm) ==
    vec([0.02 0.1 0.07 0.02 0.05])

    # QUESTION: should these be ! functions at all?

#   Test for a continuous stock vol range
#    @test Func.stockvolinit(bigm, stockvolrange) ==
#    0.1 .+ (0.5 - 0.1) .* rand(MersenneTwister(2), 5)

    Random.seed!(2)
    stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)
    # QUESTION: why do I have to assign it explicitly here and not elsewhere?

    # Generate stocks' value history
    Random.seed!(3)
    @test Func.stockvalueinit!(stocks, stockstartval, perfwindow[end],
    market.value)[:,1:2] ≈
    hcat(ones(5,1) .* 100,
    [105.51989943374585, 105.64929533253186, 104.15928328717519,
    105.05203007235151, 103.79227677999667] +
    [2.3831115468571573, -25.197330871745265, 14.523668428793831,
    -1.946500240120597, -0.508033970294714])

    Random.seed!(4)
    @test Func.stockimpactinit!(stocks.impact, impactrange, perfwindow[end]) ==
    [0.008, 0.004, 0.007, 0.001, 0.004]

    Random.seed!(4)
    Func.stockimpactinit!(stocks.impact, impactrange, perfwindow[end])

    # TODO: funds integration tests

    investors = Types.RetailInvestor(
    zeros(bign, bigk + 1),
    zeros(bign),
    zeros(bign))

    Random.seed!(29)
    @test Func.invhorizoninit!(investors.horizon, horizonrange) ==
    [3, 1, 3, 5]

    Random.seed!(6)
    @test Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd) ==
    [-0.06570018342662688, 0.00665192247499981,
    0.06774218278688489, 0.002180862789142308]

    # The first bigk investors put their money in that of the first bigk funds
    # that matches their own index
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

    # QUESTION: how do you have separate tests for A) directly after
    #initialisation, and B) following periods. I really want three kinds
    #(pre-initialisation, post-initialisation, during the runs)

    funds = Types.EquityFund(
    zeros(bigk, bigm),
    zeros(bigk, bign),
    zeros(bigk, bigt))

    @test Func.fundcapitalinit!(funds.value, investors.assets) ==
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

    @test funds.value[:,1] ≈ sum(funds.holdings .* stocks.value[1], dims=2)

    @test Func.fundvalinit!(
    funds.value, funds.holdings, stocks.value, perfwindow[end]) ≈
    hcat([318, 1003, 692],
    [327.8755848664943 ,1137.5171362972462, 730.0895512124052],
    [351.09084866012824, 1187.1302928457408, 764.4838542466905],
    zeros(3, 3))

    Random.seed!(333333333333333333)
    @test Func.drawreviewers(bign) == [3]

    Random.seed!(333333333333333333)
    reviewers = Func.drawreviewers(bign)
    funds.value[3, 4] = 730

    Random.seed!(29)
    investors.horizon .= Func.invhorizoninit!(investors.horizon, horizonrange)
    Random.seed!(6)
    investors.threshold .= Func.invthreshinit!(
    investors.threshold, thresholdmean, thresholdstd)
    @test Func.perfreview(4, reviewers, investors, funds.value) ==
    vcat(Array{Int64}(undef, 0, 2), [3 3])

    divestments = vcat([3 3], [4 2])
    @test Func.liquidate!(funds.holdings, funds.stakes, divestments) ==
    (vcat(
    [0.0 0.0 0.0 3.18 0.0],
    [0.22199999999999995 0.0 0.6659999999999997 0.22199999999999995 0.0],
    [0.0 0.0 0.0 0.0 0.0]),
    vcat(
    [-3.46 -0.0 -0.0 -3.46 -0.0 3],
    [-1.7839999999999998 -0.0 -5.351999999999999 -1.7839999999999998 -0.0 4]),
    vcat(
    [1.0 0.0 0.0 0.0],
    [0.0 1.0 0.0 0.0],
    [0.0 0.0 0.0 0.0]))

    # TODO: write test that checks that the sum of the sales orders and the
    # adjusted holdings matches the old holdings
    # Failed attempt at test that sellorder plus remaining holdings equals
    # initial holdings
    #@test Func.liquidate!(
    #funds.holdings, funds.stakes, divestments)[1][2, :] +
    #Func.liquidate!(
    #funds.holdings, funds.stakes, divestments)[2][2, 1:end-1] ==
    #funds.holdings[2, :]

    @test Func.marketmake!(stocks.value, sellorders) ==
    ([106.68540785714376 111.63061385709165 104.43603238737349
    109.79195727243169 107.73685990276314],
    vcat([3 749.011683348331],
    [4 945.1372647283855]))
#    NEW_STOCK_VALS, (INVESTOR=sellorder[end], CASH)
    # (1) NEW_STOCK_VALS: 1 + (SUM(SELLORDERS) * IMPACT_FACTOR)
    # (2) CASH: -1 .* QUANTITY OF SALES ORDER * NEW_STOCK_VALS

end

@testset "Price Functions" begin

    Random.seed!(4)
    @test Func.marketmove(100, 0.05, 0.1) ==
    100 * (1 + 0.05) + 0.1 * randn(MersenneTwister(4))

    market = Types.MarketIndex(vcat(
        mktstartval,
        zeros(bigt-1)))
    market.value[2] = 102

    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm))

    @test Func.stockmove(
    2, market.value, stocks.value[2], stocks.beta, stocks.vol)[1] ==
    ((1 + ((market.value[2]-market.value[1])/market.value[1]) * stocks.beta[1])
     * stocks.value[1]) + stocks.vol[1] * randn(MersenneTwister(2000)) *
      stocks.value[1]

end

@testset "Investor Behaviour" begin



#return history:
#horizon: 3
#threshold: 0.06774218278688489
#fund: 3
end

market = Types.MarketIndex(
    zeros(bigt))

Random.seed!(0)
Func.marketinit!(market.value, mktstartval, perfwindow[end],
drift, marketvol)

stocks = Types.Equity(
    zeros(bigm, bigt),
    zeros(bigm),
    zeros(bigm),
    zeros(bigm))

Random.seed!(1)
Func.betainit!(stocks.beta, bigm, betastd, betamean)

Random.seed!(2)
Func.stockvolinit!(stocks.vol, stockvolrange, bigm)

Random.seed!(3)
Func.stockvalueinit!(stocks, stockstartval, perfwindow[end],
market.value)

investors = Types.RetailInvestor(
zeros(bign, bigk + 1),
zeros(bign),
zeros(bign))

Random.seed!(29)
investors.horizon .= Func.invhorizoninit!(investors.horizon, horizonrange)

Random.seed!(6)
investors.threshold .= Func.invthreshinit!(investors.threshold, thresholdmean,
thresholdstd)

# The first bigk investors put their money in that of the first bigk funds
# that matches their own index
Random.seed!(7)
Func.invassetinit!(investors.assets, invcaprange, bigk)

funds = Types.EquityFund(
    zeros(bigk, bigm),
    zeros(bigk, bign),
    zeros(bigk, bigt))

Func.fundcapitalinit!(funds.value, investors.assets)

Func.fundstakeinit!(funds.stakes, investors.assets)

Random.seed!(8)
Func.fundholdinit!(
funds.holdings, portfsizerange, funds.value[:, 1], stocks.value[:, 1])

Func.fundvalinit!(
funds.value, funds.holdings, stocks.value, perfwindow[end])


Random.seed!(333333333333333333)
reviewers = Func.drawreviewers(bign)
funds.value[3, 4] = 730
Func.perfreview(4, reviewers, investors, funds.value)

# NOTE: Something is going wrong, as the code is not generating a long enough
# history of stock prices and fund values - it should be up to the top end of
# horizonrange
