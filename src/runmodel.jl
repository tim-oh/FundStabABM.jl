using Test, Random, LinearAlgebra, StatsBase, StatsPlots, Distributions

include("../src/types.jl")
include("../src/functions.jl")
include("../src/params.jl")
using .Types, .Func, .Params

# Set up types
market = Types.MarketIndex(
    zeros(Params.bigt))
stocks = Types.Equity(
    zeros(Params.bigm, Params.bigt),
    zeros(Params.bigm),
    zeros(Params.bigm),
    zeros(Params.bigm),
    zeros(Params.bigm, Params.bigt))
investors = Types.RetailInvestor(
    zeros(Params.bign, Params.bigk + 1),
    zeros(Params.bign),
    zeros(Params.bign))
funds = Func.Types.EquityFund(
    zeros(Params.bigk, Params.bigm),
    zeros(Params.bigk, Params.bign),
    zeros(Params.bigk, Params.bigt))

Func.initialise(market, stocks, investors, funds)
#Func.boundstest(market, stocks, investors, funds)
# TODO: Think about how I should have tested the initialisation
# TODO: Set test=true flag for testing bounds and other such things
#  or is that bad practice?

#
# @testset "Variable Bounds" begin
#
#     @test all(stocks.value .>= 0)
#     @test all(stocks.vol .>= 0)
#     @test all(stocks.beta .>= 0)
#     @test all(stocks.impact .>= 0)
#
#     @test all(funds.holdings .>= 0)
#     @test all(funds.stakes .>= 0)
#     @test all(funds.value .>= 0)
#
#     @test all(investors.assets .>= 0)
#     @test all(investors.horizon .>= 0)
#     @test all(1 .>= investors.threshold .>= -1)
#
# end

# TODO: Create logging set-up that records parameters, plots and outputs
# TODO: Think about how I should have tested the running cycle
# TODO: Write functions that collect activity of a) selling out and b) buying in
# TODO: Tests for cases with empty reviewers, divestments
# QUESTION: Need 'if reviewers not empty'/'divestments not empty' branches?
# Looping over time after the end of the  perfwindow[end]+1 initialisation phase

Func.modelrun(market, stocks, investors, funds)

@testset "Variable Bounds" begin

    @test all(stocks.value .>= 0)
    @test all(stocks.vol .>= 0)
    @test all(stocks.impact .>= 0)

    @test all(funds.holdings .>= 0)
    @test all(funds.stakes .>= 0)
    @test all(funds.value .>= 0)

    @test all(investors.assets .>= 0)
    @test all(investors.horizon .>= 0)
    @test all(1 .>= investors.threshold .>= -1)

end # testset

pyplot()
StatsPlots.plot(1:1000, transpose(stocks.value)[:, 1:100], legend=:none,
 title="History of stock returns")
StatsPlots.plot!(1:1000, market.value, lw=3, lc=:black, label="Market index")
png("returnshistory")

# Calculate returns
marketreturns = Func.assetreturns(market.value)
stockreturns = Func.assetreturns(stocks.value)

# Subtract mean returns
demeanedmarketreturns = Func.demean(marketreturns)
demeanedstockreturns = Func.demean(stockreturns)

# Plotting a histogram overlaid with a normal distribution
asset = demeanedmarketreturns
StatsPlots.histogram(asset, bins=100, title="Market returns vs normal distribution", label="Market returns", legend=:topright)
x = findmin(asset)[1]:0.001:findmax(asset)[1]
d = fit(Normal, asset)
StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
png("marketreturnsvsnormal")

# Plotting a histogram overlaid with a normal distribution
choice = 1
asset = demeanedstockreturns[choice, :]
StatsPlots.histogram(asset, bins=100, title="Example of stock returns", label="Random stock's returns", legend=:topright)
x = findmin(asset)[1]:0.001:findmax(asset)[1]
d = fit(Normal, asset)
StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
png("stockreturnsvsnormal")

pacfcoeffs = pacf(demeanedmarketreturns, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
StatsPlots.bar(pacfcoeffs, title="Market autocorrelation", label="Partial correlation coefficients")
StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)),], lc=:red, label="critical values", legend=:bottomright)
png("marketautocorrelation")

examplestock = demeanedstockreturns[1,: ]
pacfcoeffs = pacf(demeanedstockreturns[1,:], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
StatsPlots.bar(pacfcoeffs, title="Stock autocorrelation (random stock example)", label="Partial autocorrelation coefficients, demeaned returns", legend=:bottomright)
StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(examplestock, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(examplestock, 1)),], lc=:red, label="Critical values")
png("stockautocorrelation")

demeanedmarketreturnssquared = demeanedmarketreturns.^2
pacfcoeffs = pacf(demeanedmarketreturnssquared, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
StatsPlots.bar(pacfcoeffs, title="Volatility clustering, market", label="Partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)),], lc=:red, label="Critical values")
png("marketvolaclustering")

demeanedstockreturnssquared = demeanedstockreturns.^2
pacfcoeffs = pacf(demeanedstockreturnssquared[1,:], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
StatsPlots.bar(pacfcoeffs, title="Volatility clustering, random stock", label="partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(marketreturns, 1)),], lc=:red, label="Critical values")
png("stockvolaclustering")

# Calculate excess kurtosis for all stocks and the market
stockkurtoses = Func.calckurtoses(demeanedstockreturns)
StatsPlots.bar(sort!(stockkurtoses, rev=true), title="Kurtoses of stock returns", label="Kurtoses")
png("kurtoses")

cutoff = 5 # Percentile of absolute returns that determines 'large' return
ratios = Func.lossgainratio(demeanedstockreturns, cutoff)
plot(ratios, title="Gain-loss asymmetry in stock returns", label="# large losses / # large gains")
plot!(0:0.01:200, ones(length(0:0.01:200)) * mean(ratios), label="Mean ratio")
png("lossgainratio")

corrs = Func.volavolumecorr(stocks.volume, demeanedstockreturns)
bar(sort(corrs, rev=true), title="Volume-volatility correlation of stocks", label="Volume - absolute return correlation coefficients")
png("volavolumecorr")
