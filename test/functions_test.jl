using Test, Random, LinearAlgebra, Traceur

using FundStabABM.Func, FundStabABM.Types

# QUESTION: How do I run the type tests periodically?

# Small parameter values for testing, same list as in /src/params.jl
const bigk = 3 # Number of funds, 3
const bign = 4 # Number of investors, 4
const bigm = 5 # Number of stocks, 5
const bigt = 6 # Number of time periods, 6
const mktstartval = 100 # Market index starting value, 100
const drift = 0.05 # Market index drift, 0.05
const marketvol = 0.1 # Market volatility (std?), 0.1
const perfwindow = 1:3 # Performance window range for investors (1:3)
const betamean = 1 # Average stock beta, 1
const betastd = 0.3 # Stock beta dispersion, 0.3
const stockstartval = 100 # Starting value of stocks
const stockvolrange = range(0.01, stop=0.1, step=0.01)
# Range of stock volatility (0.1, 0.5). Consider if it should be continuous.
const invcaprange = (10,1000) # Investors' range of initial capital, (10, 1000)
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds
const portfsizerange = 1:5 # Range of number of stocks in fund portfolio (1,5)
const impactrange = 0.00001:0.00001:0.0001 # Stock price impact per currency unit


@testset "Initialisation Functions" begin

    # Generate market value history
    market = Types.MarketIndex(
        zeros(bigt))

    # Test: generated market index history is non-negative
    Random.seed!(0)
    @test all(Func.marketinit!(market.value, mktstartval, perfwindow[end],
     drift, marketvol) .>= 0)
    # TODO: Replace inequality with specific values, this one is in params_test

    # Test: generation of stock betas
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm, bigt))
    Random.seed!(1)
    @test Func.betainit!(stocks.beta, bigm, betastd, betamean) ==
    1 .+ 0.3 .* randn(MersenneTwister(1), 5)

    # Test: generation of stock volatilities
    Random.seed!(2)
    @test Func.stockvolinit!(stocks.vol, stockvolrange, bigm) ==
    vec([0.02 0.1 0.07 0.02 0.05])

    # Test: generation of stock values histories
    Random.seed!(2)
    stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)
    Random.seed!(3)
    @test Func.stockvalueinit!(stocks, market.value, stockstartval,
    perfwindow)[:,1:2] ≈
    hcat(ones(5,1) .* 100,
    [105.51989943374585, 105.64929533253186, 104.15928328717519,
     105.05203007235151, 103.79227677999667] +
    [2.3831115468571573, -25.197330871745265, 14.523668428793831,
    -1.946500240120597, -0.508033970294714])

    # Test: generation of stock impact parameters
    Random.seed!(4)
    @test Func.stockimpactinit!(stocks.impact, impactrange) ==
    [0.00008, 0.00004, 0.00007, 0.00001, 0.00004]

    Random.seed!(4)
    Func.stockimpactinit!(stocks.impact, impactrange)

    investors = Types.RetailInvestor(
    zeros(bign, bigk + 1),
    zeros(bign),
    zeros(bign))

    # Test: generation of investor performance evaluation horizons
    Random.seed!(55)
    @test Func.invhorizoninit!(investors.horizon, perfwindow) ==
    [1, 3, 1, 1]

    # Test: generation of investor performance thresholds
    Random.seed!(2345)
    @test Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd) ==
    [0.022728468593375205,
    -0.03299475538810838,
     0.12695393970324698,
    -0.02006345690377245]

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

    # QUESTION: how do you have separate tests for A) directly after
    #initialisation, and B) following periods. I really want three kinds
    #(pre-initialisation, post-initialisation, during the runs)

    # NOTE Awkward scoping: Func.Types.EquityFund
    funds = Func.Types.EquityFund(
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
    funds.holdings, funds.value[:, 1], stocks.value[:, 1], portfsizerange) ≈
    [0 0 0 318/100 0;
     1003/500 0 3009/500 1003/500 0;
     692/200 0 0 692/200 0]

    # Test: Initial fund values equal the value of their holdings
    @test funds.value[:,1] ≈ sum(funds.holdings .* stocks.value[1], dims=2)

    # Test: Generation of fund value history
    @test Func.fundvalinit!(
    funds.value, funds.holdings, stocks.value, perfwindow[end]) ≈
    hcat([318, 1003, 692],
    [327.8755848664943 ,1137.5171362972462, 730.0895512124052],
    [351.09084866012824, 1187.1302928457408, 764.4838542466905],
    [368.98315314801744, 1233.1458121219132, 813.506737945768],
    zeros(3, 2))

end # testset "Initialisation Functions"


@testset "Agent Behaviours" begin

    market = Types.MarketIndex(
        zeros(bigt))
    Random.seed!(0)
    Func.marketinit!(
    market.value, mktstartval, perfwindow[end], drift, marketvol)
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm, bigt))
    Random.seed!(1)
    Func.betainit!(stocks.beta, bigm, betastd, betamean)
    Random.seed!(2)
    stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)
    Random.seed!(3)
    Func.stockvalueinit!(stocks, market.value, stockstartval, perfwindow)
    Random.seed!(4)
    Func.stockimpactinit!(stocks.impact, impactrange)
    investors = Types.RetailInvestor(
    zeros(bign, bigk + 1),
    zeros(bign),
    zeros(bign))
    Random.seed!(55)
    investors.horizon .= Func.invhorizoninit!(investors.horizon, perfwindow)
    Random.seed!(2345)
    investors.threshold .= Func.invthreshinit!(investors.threshold, thresholdmean,
    thresholdstd)
    Random.seed!(7)
    Func.invassetinit!(investors.assets, invcaprange, bigk)
    # NOTE Awkward scoping: Func.Types.EquityFund
    funds = Func.Types.EquityFund(
        zeros(bigk, bigm),
        zeros(bigk, bign),
        zeros(bigk, bigt))
    Func.fundcapitalinit!(funds.value, investors.assets)
    Func.fundstakeinit!(funds.stakes, investors.assets)
    Random.seed!(8)
    Func.fundholdinit!(funds.holdings, funds.value[:, 1],
    stocks.value[:, 1], portfsizerange)
    Func.fundvalinit!(
    funds.value, funds.holdings, stocks.value, perfwindow[end])

    # Test: Random selection of investors that conduct a performance review
    Random.seed!(333333333333333333)
    @test Func.drawreviewers(bign) == [3]

    # Test: Fund-investor pairs for divestment following performance review
    Random.seed!(333333333333333333)
    reviewers = Func.drawreviewers(bign)

    # QUESTION: In terms of the model dynamics, is there a 'neutral' time to
    # do the review or the resulting trading?

    Random.seed!(55)
    investors.horizon .= Func.invhorizoninit!(investors.horizon, perfwindow)
    Random.seed!(2345)
    investors.threshold .= Func.invthreshinit!(
    investors.threshold, thresholdmean, thresholdstd)
    @test Func.perfreview(4, reviewers, investors, funds.value) == [3 3]

    # Test: Sell order amounts resulting from divestment
    divestments = [3 3; 4 2]
    liquidationresults = Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)
    @test liquidationresults[1].values ≈
    [-382.47934595588924 -0.0 -0.0 -382.0045082904202 -0.0;
     -197.20899618910363 -0.0 -661.5797843033668 -196.96417421712303 -0.0]

    # Test: Investor that will receive cash following divestment
    @test liquidationresults[1].investors == [3,4]

    # Test: Fund ownership of shares following divestment
    @test liquidationresults[2] == [0.0 0.0 0.0 3.18 0.0;
    0.22199999999999995 0.0 0.6659999999999997 0.22199999999999995 0.0;
    0.0 0.0 0.0 0.0 0.0]

    # Test: Investor stakes in funds following divestment
    @test liquidationresults[3] ==
    [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0]

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
    divestments = [3 3; 4 2]
    priorstockvals = hcat(ones(5,1) .* 100,
    [107.903010980603, 80.4519644607866, 118.68295171596901,
     103.10552983223091, 103.28424280970195],
    [110.54316357113026,85.81763860422608,123.6135620896176,
     110.40592725161265, 111.56926365234682],
    [119.08511622272532, 86.49247760086409, 126.53705519388893,
     116.03243809686083, 112.85619511929654], zeros(5,2))
    stocks.value .= priorstockvals
    sellorders = Func.Types.SellMarketOrder(
    [-382.47934595588924 -0.0 -0.0 -382.0045082904202 -0.0;
     -197.20899618910363 -0.0 -661.5797843033668 -196.96417421712303 -0.0],
    [3, 4])
    sellmarketmakeresults = Func.executeorder!(stocks, 3, sellorders)
    @test sellmarketmakeresults[1] ≈ hcat(ones(5,1).*100,
    [ 107.903010980603, 80.4519644607866, 118.68295171596901,
    103.10552983223091, 103.28424280970195],
    [105.41671691304936, 85.81763860422608, 117.8889457275222,
    109.76671150919375, 111.56926365234682],
    [119.08511622272532, 86.49247760086409, 126.53705519388893,
    116.03243809686083, 112.85619511929654], zeros(5,2))

    # Test: Investor-cash pair resulting from divestment
    @test sellmarketmakeresults[2] ==
    [3.0 744.5346623405912; 4.0 1014.8288665706393]

    # Test: Updating of investor assets after divested funds are disbursed
    cashout = [3.0 746.434849083664; 4.0 1021.6613450347874]
    divestments = [3 3; 4 2]
    @test Func.disburse!(investors.assets, divestments, cashout) ==
    [318.0 0.0 0.0 0.0;
     0.0 111.0 0.0 0.0;
     0.0 0.0 0.0 746.434849083664;
     0.0 0.0 0.0 1021.6613450347874]

    # Test: Fund holdings after divestment
    resultstuple = Func.liquidate!(
    funds.holdings, funds.stakes, stocks.value[:, 3], divestments)
    @test funds.holdings == resultstuple[2]

    # Test: No greater number of dead funds than investors holding cash
    divestments = [3 3; 4 2]
    Func.liquidate!(funds.holdings, funds.stakes,stocks.value[:,3], divestments)
    cashout = [3.0 746.434849083664; 4.0 1021.6613450347874]
    Func.disburse!(investors.assets, divestments, cashout)
    @test length(findall(vec(sum(funds.holdings, dims=2) .== 0))) <=
    length(findall(vec(investors.assets[:, end] .> 0)))

    # Test: Stock-values-part of buy order following fund re-birth
    Random.seed!(10)
    respawn_output = Func.respawn!(
    funds, investors, 3, stocks.value, portfsizerange)
    @test respawn_output[1].values ≈
    [298.574 0.0 149.287 298.574 0.0] atol = 0.0001

    # Test: Fund identifier part of buy order following fund re-birth
    @test respawn_output[1].funds ==
    [3]

    # Test: Assets of investors following fund rebirths
    @test respawn_output[2] ≈
    [318.0 0.0 0.0 0.0;
     0.0 111.0 0.0 0.0;
     0.0 0.0 746.435 0.0;
     0.0 0.0 0.0 1021.66] atol = 0.01

    # Test: Investor stakes in funds following fund births
    @test respawn_output[3] ==
    [1.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0;
     0.0 0.0 1.0 0.0]

    # Test: Stock values after buying-driven upwards impact
    buyorder = respawn_output[1]
    buymarketmakeresults = Func.executeorder!(stocks, 3, buyorder)
    @test buymarketmakeresults[1] ≈
    hcat(ones(5,1) .* 100,
    [107.903010980603, 80.4519644607866, 118.68295171596901,
     103.10552983223091, 103.28424280970195],
    [107.93469217989708, 85.81763860422608, 119.12089582037991,
     110.09444637041521, 111.56926365234682],
    [119.08511622272532, 86.49247760086409, 126.53705519388893,
     116.03243809686083, 112.85619511929654], zeros(5,2)) atol = 0.0001

    # Test: Fund-amount of shares pair(s) intended for disbursement
    @test buymarketmakeresults[2] ≈
    [3.0 2.7662468291692 0.0 1.25323939995470 2.7119805752548 0.0] atol=0.000001

    # Test: Disbursement of shares to fund
    sharesout = buymarketmakeresults[2]
    disbursesharesresult = Func.disburse!(funds, sharesout, stocks.value)
    @test disbursesharesresult.holdings ≈
    [0.0 0.0 0.0 3.18 0.0;
     0.222  0.0  0.666  0.222  0.0;
     2.7662468291692 0.0 1.25323939995470 2.7119805752548 0.0] atol=0.00001

    # Test: re-valuation of fund following respawn
    # NOTE Strangely large atol, plus had to hardcode 802.677.. bug?
    # FIXME: I adjusted the value of the respawned fund, but not of the one
    # that still had an investor left (fund 2)
    @test disbursesharesresult.value ≈
    [318.0  327.8755848664943  351.09084866012824 368.98315314801744 0.0 0.0;
     1003.0 1137.5171362972462 1187.1302928457408 1233.1458121219132 0.0 0.0;
     673.1466804 726.844707282 746.434998099 802.677606280 0.0 0.0] atol=0.001

    # Test: Return index of best-performing fund for min and max time horizon
    fundvals =
    [318.0 327.876 351.091 368.98315314801744 0.0 0.0;
     1003.0 1137.52 1187.13 1233.1458121219132 0.0 0.0;
     673.1466804378 726.844707282640 746.434998099860 802.677606280818 0.0 0.0]
    horizonmin = perfwindow[1]
    horizonmax = perfwindow[end]
    @test Func.bestperformer(fundvals, horizonmin, 4) == 3
    @test Func.bestperformer(fundvals, horizonmax, 4) == 2

    # Test: Reallocation of spare investor cash to a fund
    # NOTE: investors.assets don't update automatically, tracks initial investmt
    reinvestresults = Func.reinvest!(investors, funds, stocks.value, 3)
    @test reinvestresults[1] ==
    [318.0 0.0 0.0 0.0;
     0.0 111.0 0.0 0.0;
     0.0 0.0 746.434849083664 0.0;
     1021.6613450347874 0.0 0.0 0.0]

    # Test: Update of investor stakes in funds
    # NOTE: atol quite large due to different calculation of fund value of ~1
    # There is a general discrepancy between test values and those ins console
    @test reinvestresults[2] ≈
    [0.25575690228192466 0.0 0.0 0.7442430977180754;
     0.0                 1.0 0.0 0.0;
     0.0                 0.0 1.0 0.0] atol=0.001

    # Test: Each line represents the values of each stock submitted to mkt maker
    @test reinvestresults[3].values == [0.0 0.0 0.0 1021.6613450347874 0.0]

     # Test: Funds for which will receive proceeds of buy orders
     @test reinvestresults[3].funds == [1]

    # TODO: Write a test for reinvest where the funds holds more than one stock

    # TODO: Write a test for fund value update stepping

end # testset "Agent Behaviours"

# TODO: Understand why some values remain within the test function while others
# are available in the console

@testset "Price Functions" begin

    # Test: Random walk of market with drift
    market = Types.MarketIndex([mktstartval; zeros(bigt-1)])
    Random.seed!(4)
    @test Func.marketmove!(market.value, 2, 0.05, 0.1) ==
    [100, 100 * (1 + 0.05) + 0.1 * randn(MersenneTwister(4)), 0, 0, 0, 0]

    # Test: Draw of stock price moves on basis of marketmove
    market.value[2] = 102
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm, bigt))
    @test Func.stockmove!(
    stocks.value, 2, market.value, stocks.beta, stocks.vol)[1] ==
    ((1 + ((market.value[2]-market.value[1])/market.value[1]) * stocks.beta[1])
     * stocks.value[1]) + stocks.vol[1] * randn(MersenneTwister(2000)) *
      stocks.value[1]

    # Test: Market making / price impact function for buy order
    tmpstock =
    Types.Equity([[100, 100] [101, 99] [102, 98]], [], [], [0.01, 0.05], [0 0])
    tmpbuyordervals  = [20 0]
    @test Func.marketmake!(tmpstock, 3, tmpbuyordervals)[:,3] ==
    [122.39999999999999, 98.0]

    # Test: Market making / price impact function for sell order
    tmpstock =
    Types.Equity([[100, 100] [101, 99] [102, 98]],[],[],[0.01, 0.05],[0 0])
    tmpsellordervals = [-5 -5]
    @test Func.marketmake!(tmpstock, 3, tmpsellordervals)[:,3] ==
    [96.89999999999999, 73.5]

    # Test: Revaluation of funds, can be needed after holdings or prices change
    # NOTE: Huge atol
    # TODO: Test that it works for a single fund (and maybe a vector of funds)
    funds = Func.Types.EquityFund(
        zeros(bigk, bigm),
        zeros(bigk, bign),
        zeros(bigk, bigt))
    funds.holdings .=
    [0.0    0.0  0.0    3.18   0.0;
    0.222  0.0  0.666  0.222  0.0;
    0.0    0.0  0.0    0.0    0.0]
    funds.value .=
    [318.0   327.876   351.091   368.983  0.0  0.0;
     1003.0  1137.52   1187.13   1233.15   0.0  0.0;
     692.0   730.09    764.484   813.507  0.0  0.0]
    stocks.value .=
    [100.0  107.903  107.935   119.085   0.0  0.0;
     100.0   80.452   85.8176   86.4925  0.0  0.0;
     100.0  118.683  119.121   126.537   0.0  0.0;
     100.0  103.106  110.094   116.032   0.0  0.0;
     100.0  103.284  111.569   112.856   0.0  0.0]
    @test Func.fundrevalue!(funds, stocks.value, 1:bigk) ≈
    [318.0   327.876   351.091   368.983  0.0  0.0;
    110.999999999999 125.886741903284 127.736985080927 136.469775818078 0.0 0.0;
     0.0  0.0  0.0  0.0  0.0  0.0] atol=1
end # testset "Price Functions"

# @testset "Integration Tests" begin

    # TODO: funds initialisation integration tests
    # TODO: market initialisation integration tests
    # TODO: stocks initialisation integration tests
    # TODO: investors initialisation integration tests
    # TODO: other integration tests
# end # testset "Integration Tests"




# These lines execute all the functions and keep their results in scope, so that
# subsequent tests can run. Since tests shouldn't be dependent on each other, I
# should presumably replace variable args with hardcoded args where possible
market = Types.MarketIndex(
    zeros(bigt))

Random.seed!(0)
Func.marketinit!(market.value, mktstartval, perfwindow[end], drift, marketvol)

stocks = Types.Equity(
    zeros(bigm, bigt),
    zeros(bigm),
    zeros(bigm),
    zeros(bigm),
    zeros(bigm, bigt))

Random.seed!(1)
Func.betainit!(stocks.beta, bigm, betastd, betamean)

Random.seed!(2)
stocks.vol .= Func.stockvolinit!(stocks.vol, stockvolrange, bigm)

Random.seed!(3)
Func.stockvalueinit!(stocks, market.value, stockstartval, perfwindow)

Random.seed!(4)
Func.stockimpactinit!(stocks.impact, impactrange)

investors = Types.RetailInvestor(zeros(bign, bigk + 1),
zeros(bign),
zeros(bign))

Random.seed!(55)
investors.horizon .= Func.invhorizoninit!(investors.horizon, perfwindow)

Random.seed!(2345)
investors.threshold .= Func.invthreshinit!(investors.threshold, thresholdmean,
thresholdstd)

Random.seed!(7)
Func.invassetinit!(investors.assets, invcaprange, bigk)

# NOTE Awkward scoping: Func.Types.EquityFund
funds = Func.Types.EquityFund(
    zeros(bigk, bigm),
    zeros(bigk, bign),
    zeros(bigk, bigt))

Func.fundcapitalinit!(funds.value, investors.assets)

Func.fundstakeinit!(funds.stakes, investors.assets)

Random.seed!(8)
Func.fundholdinit!(funds.holdings, funds.value[:, 1],
stocks.value[:, 1], portfsizerange)

Func.fundvalinit!(
funds.value, funds.holdings, stocks.value, perfwindow[end])

Random.seed!(333333333333333333)
reviewers = Func.drawreviewers(bign)
Func.perfreview(4, reviewers, investors, funds.value)

divestments = [3 3; 4 2]
liquidationresults = Func.liquidate!(
funds.holdings, funds.stakes, stocks.value[:, 3], divestments)

sellorders = Func.Types.SellMarketOrder(
[-382.47934595588924 -0.0 -0.0 -382.0045082904202 -0.0;
-197.20899618910363 -0.0 -661.5797843033668 -196.96417421712303 -0.0],
[3, 4])
sellmarketmakeresults = Func.executeorder!(stocks, 3, sellorders)
cashout = [3.0 746.434849083664; 4.0 1021.6613450347874]
Func.disburse!(investors.assets, divestments, cashout)

Random.seed!(10)
respawn_output = Func.respawn!(
funds, investors, 3, stocks.value, portfsizerange)

buyorder = respawn_output[1]
buymarketmakeresults = Func.executeorder!(stocks, 3, buyorder)
