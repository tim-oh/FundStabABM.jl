using Test, Random, LinearAlgebra

include("../src/types.jl")
using .Types

# Small parameter values for testing, same list as in /src/params.jl
const bigk = 3 # Number of funds, 3
const bign = 4 # Number of investors, 4
const bigm = 5 # Number of stocks, 5
const bigt = 6 # Number of time periods, 6
const mktstartval = 100 # Market index starting value, 100
const drift = 0.05 # Market index drift, 0.05
const marketvol = 0.1 # Market volatility (std?), 0.1
const perfwindow = 1:3 # Performance window for investor (1,3)
const betamean = 1 # Average stock beta, 1
const betastd = 0.3 # Stock beta dispersion, 0.3
const stockstartval = 100 # Starting value of stocks
const stockvolrange = range(0.01, stop=0.1, step=0.01)
# Range of stock volatility (0.1, 0.5). Consider if it should be continuous.
const invcaprange = (10,1000) # Investors' range of initial capital, (10, 1000)
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

    # Test: generated market index history is non-negative
    Random.seed!(0)
    @test all(Func.marketinit!(market.value, mktstartval, perfwindow[end],
    drift, marketvol)
    .>= 0)
    # TODO: Replace inequality with specific values, this one is in params_test

    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm))

    # Test: generation of stock betas
    Random.seed!(1)
    @test Func.betainit!(stocks.beta, bigm, betastd, betamean) ==
    1 .+ 0.3 .* randn(MersenneTwister(1), 5)

    # Test: generation of stock volatilities
    Random.seed!(2)
    @test Func.stockvolinit!(stocks.vol, stockvolrange, bigm) ==
    vec([0.02 0.1 0.07 0.02 0.05])

    Random.seed!(2)
    stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)
    # QUESTION: why do I have to assign it explicitly here and not elsewhere?

    # Test: generation of stock values histories
    Random.seed!(3)
    @test Func.stockvalueinit!(stocks, stockstartval, perfwindow[end],
    market.value)[:,1:2] ≈
    hcat(ones(5,1) .* 100,
    [105.51989943374585, 105.64929533253186, 104.15928328717519,
    105.05203007235151, 103.79227677999667] +
    [2.3831115468571573, -25.197330871745265, 14.523668428793831,
    -1.946500240120597, -0.508033970294714])

    # Test: generation of stock impact parameters
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

    # Test: generation of investor performance evaluation horizons
    Random.seed!(29)
    @test Func.invhorizoninit!(investors.horizon, horizonrange) ==
    [3, 1, 3, 5]

    # Test: generation of investor performance thresholds
    Random.seed!(6)
    @test Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd) ==
    [-0.06570018342662688, 0.00665192247499981,
    0.06774218278688489, 0.002180862789142308]

    # Test: first bigk investors put their money in funds with matching indices
    Random.seed!(7)
    @test LinearAlgebra.diag(Func.invassetinit!(investors.assets,
    invcaprange, bigk)[1:bigk, 1:bigk]) == [318, 111, 692]

    # Test: no investor holds cash at the beginning
    @test all(investors.assets[:, end] .== 0)

    # Test: Sum of investments equals funds' combined capital
    @test sum(investors.assets) == 2013

    # Test: One investment per investor
    @test sum(investors.assets .> 0, dims=2) == ones(Int64, (bign,1))

    # Test: At least one investor per fund
    @test all(sum(investors.assets[:,1:end-1] .> 0, dims=1) .>= 1)

    #TODO: investors integration tests

    # QUESTION: how do you have separate tests for A) directly after
    #initialisation, and B) following periods. I really want three kinds
    #(pre-initialisation, post-initialisation, during the runs)

    funds = Types.EquityFund(
    zeros(bigk, bigm),
    zeros(bigk, bign),
    zeros(bigk, bigt))

    # Test: Generation of funds' initial capital (== value of their assets)
    @test Func.fundcapitalinit!(funds.value, investors.assets) ==
    hcat([318.0, 1003.0, 692.0], zeros(3,5))

    # Test: Generation of investors' stakes in funds
    @test Func.fundstakeinit!(funds.stakes, investors.assets) ==
    [1.0 0.0 0.0 0.0;
    0.0 0.1106679960119641 0.0 0.8893320039880359;
    0.0 0.0 1.0 0.0]

    # Test: Investors' stakes sum to one for each fund
    @test sum(funds.stakes, dims=2) == ones(bigk,1)

    # Test: Generation of fund holdings
    Random.seed!(8)
    @test Func.fundholdinit!(
    funds.holdings, portfsizerange, funds.value[:, 1], stocks.value[:, 1]) ≈
    [0 0 0 318/100 0; 1003/500 0 3009/500 1003/500 0;
    692/200 0 0 692/200 0]

    # Test: Initial fund values equal the value of their holdings
    @test funds.value[:,1] ≈ sum(funds.holdings .* stocks.value[1], dims=2)

    # Test: Generation of fund value history
    @test Func.fundvalinit!(
    funds.value, funds.holdings, stocks.value, perfwindow[end]) ≈
    hcat([318, 1003, 692],
    [327.8755848664943 ,1137.5171362972462, 730.0895512124052],
    [351.09084866012824, 1187.1302928457408, 764.4838542466905],
    zeros(3, 3))

end # testset "Initialisation Functions"


@testset "Agent Behaviours" begin

    market = Types.MarketIndex(
        zeros(bigt))
    Random.seed!(0)
    Func.marketinit!(market.value, mktstartval, perfwindow[end], drift, marketvol)
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm))
    Random.seed!(1)
    Func.betainit!(stocks.beta, bigm, betastd, betamean)
    Random.seed!(2)
    stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)
    Random.seed!(3)
    Func.stockvalueinit!(stocks, stockstartval, perfwindow[end], market.value)
    Random.seed!(4)
    Func.stockimpactinit!(stocks.impact, impactrange, perfwindow[end])
    investors = Types.RetailInvestor(
    zeros(bign, bigk + 1),
    zeros(bign),
    zeros(bign))
    Random.seed!(29)
    investors.horizon .= Func.invhorizoninit!(investors.horizon, horizonrange)
    Random.seed!(6)
    investors.threshold .= Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd)
    Random.seed!(7)
    Func.invassetinit!(investors.assets, invcaprange, bigk)
    funds = Types.EquityFund(
        zeros(bigk, bigm),
        zeros(bigk, bign),
        zeros(bigk, bigt))
    Func.fundcapitalinit!(funds.value, investors.assets)
    Func.fundstakeinit!(funds.stakes, investors.assets)
    Random.seed!(8)
    Func.fundholdinit!(funds.holdings, portfsizerange, funds.value[:, 1],
    stocks.value[:, 1])
    Func.fundvalinit!(
    funds.value, funds.holdings, stocks.value, perfwindow[end])

    # Test: Random selection of investors that conduct a performance review
    Random.seed!(333333333333333333)
    @test Func.drawreviewers(bign) == [3]

    # Test: Fund-investor pairs for divestment following performance review
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

    # Test: Sell order resulting from divestment
    divestments = vcat([3 3], [4 2])
    @test Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)[1].values ≈
    vcat(
    [-382.47934595588924  -0.0    -0.0   -382.0045082904202  -0.0],
    [-197.20899618910363  -0.0  -661.5797843033668 -196.96417421712303  -0.0])

    # Test: Investor that will receive cash following divestment
    @test Func.liquidate!(
    funds.holdings, funds.stakes,stocks.value[:,3], divestments)[1].investors ==
    [3,4]

    # Test: Fund ownership of shares following divestment
    @test Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)[2] ==
    vcat(
    [0.0 0.0 0.0 3.18 0.0],
    [0.22199999999999995 0.0 0.6659999999999997 0.22199999999999995 0.0],
    [0.0 0.0 0.0 0.0 0.0])

    # Test: Investor stakes in funds following divestment
    @test Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)[3] ==
    vcat(
    [1.0 0.0 0.0 0.0],
    [0.0 1.0 0.0 0.0],
    [0.0 0.0 0.0 0.0])

    # TODO: write test that checks that the sum of the sales orders and the
    # adjusted holdings matches the old holdings
    # Failed attempt at test that sellorder plus remaining holdings equals
    # initial holdings
    #@test Func.liquidate!(
    #funds.holdings, funds.stakes, divestments)[1][2, :] +
    #Func.liquidate!(
    #funds.holdings, funds.stakes, divestments)[2][2, 1:end-1] ==
    #funds.holdings[2, :]

    # Test: Market impact resulting from divestment-driven sales
    divestments = vcat([3 3], [4 2])
    stocksvaltst = hcat(ones(5,1) .* 100,
    [107.903010980603, 80.4519644607866, 118.68295171596901,
    103.10552983223091, 103.28424280970195],
    [110.54316357113026,85.81763860422608,123.6135620896176,
    110.40592725161265, 111.56926365234682], zeros(5,3))
    stocks.value .= stocksvaltst
    sellorders = Func.Types.SellMarketOrder(vcat(
    [-3.46 -0.0 -0.0 -3.46 -0.0], [-1.784  -0.0  -5.352  -1.784  -0.0]), [3, 4])
    @test Func.marketmake!(stocks.value, stocks.impact, 3, sellorders)[1] ≈
    hcat(ones(5,1).*100,
    [ 107.903010980603, 80.4519644607866, 118.68295171596901,
    103.10552983223091, 103.28424280970195],
    [ 105.9056567729942, 85.81763860422608, 118.98250359949218,
    109.82695856910519, 111.56926365234682], zeros(5,3))

    # Test: Investor-cash pair resulting from divestment
    stocks.value .= stocksvaltst
    @test Func.marketmake!(stocks.value, stocks.impact, 3, sellorders)[2] ==
    vcat([3.0 746.434849083664], [4.0 1021.6613450347874])

    # Test: Updating of investor assets after divested funds are disbursed
    cashout = vcat([3.0 746.434849083664], [4.0 1021.6613450347874])
    divestments = vcat([3 3], [4 2])
    @test Func.disburse!(investors.assets, divestments, cashout) ==
    vcat([318.0 0.0 0.0 0.0],
        [0.0 111.0 0.0 0.0],
        [0.0 0.0 0.0 746.434849083664],
        [0.0 0.0 0.0 1021.6613450347874])

    # Test: Fund holdings after divestment
    resultstuple = Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)
    @test funds.holdings == resultstuple[2]

    # Test: No more dead funds than investors holding cash
    Func.liquidate!(funds.holdings, funds.stakes,stocks.value[:,3], divestments)
    cashout = vcat([3.0 746.434849083664], [4.0 1021.6613450347874])
    divestments = vcat([3 3], [4 2])
    Func.disburse!(investors.assets, divestments, cashout)
    @test length(findall(vec(sum(funds.holdings, dims=2) .== 0))) <=
    length(findall(vec(investors.assets[:, end] .> 0)))

    # Test: buyorders
    @test_broken Func.respawn!(funds, investors, 3, stocks.value)[1] == 0

    # Test: investors.assets
    @test_broken Func.respawn!(funds, investors, 3, stocks.value)[2] == 0

    # Test: funds.stakes
    @test_broken Func.respawn!(funds, investors, 3, stocks.value)[3] == 0

end # testset "Agent Behaviours"


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

end # testset "Price Functions"

@testset "Integration Tests" begin

end # testset "Integration Tests"





market = Types.MarketIndex(
    zeros(bigt))

Random.seed!(0)
Func.marketinit!(market.value, mktstartval, perfwindow[end], drift, marketvol)

stocks = Types.Equity(
    zeros(bigm, bigt),
    zeros(bigm),
    zeros(bigm),
    zeros(bigm))

Random.seed!(1)
Func.betainit!(stocks.beta, bigm, betastd, betamean)

Random.seed!(2)
stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)

Random.seed!(3)
Func.stockvalueinit!(stocks, stockstartval, perfwindow[end], market.value)

Random.seed!(4)
Func.stockimpactinit!(stocks.impact, impactrange, perfwindow[end])

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
Func.fundholdinit!(funds.holdings, portfsizerange, funds.value[:, 1],
stocks.value[:, 1])

Func.fundvalinit!(
funds.value, funds.holdings, stocks.value, perfwindow[end])


Random.seed!(333333333333333333)
reviewers = Func.drawreviewers(bign)
funds.value[3, 4] = 730
Func.perfreview(4, reviewers, investors, funds.value)
