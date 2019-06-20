using Test, Random, LinearAlgebra, StatsBase, StatsPlots, Distributions
using Traceur, Parameters

using .Types, .Func, .Params

function runmodel(;params=Params.default(), fundselector="probabilistic")
    @unpack bigk, bigm, bign, bigt, impactrange, perfwindow, plotpath = params
    println("Run for ", bigt, " periods")
    println("Stock impact range:", impactrange)
    # Set up types
    market = Types.MarketIndex(
        zeros(bigt))
    stocks = Types.Equity(
        zeros(bigm, bigt),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm),
        zeros(bigm, bigt))
    investors = Types.RetailInvestor(
        zeros(bign, bigk + 1),
        zeros(bign),
        zeros(bign))
    funds = Func.Types.EquityFund(
        zeros(bigk, bigm),
        zeros(bigk, bign),
        zeros(bigk, bigt))

    Func.initialise(market, stocks, investors, funds, params)
    #Func.boundstest(market, stocks, investors, funds)
    Func.modelrun(market, stocks, investors, funds, params, fundselector)
    Func.boundstest(market, stocks, investors, funds)
    Func.plot_stylisedfacts(market.value, stocks.value, stocks.volume, params)
    # save("examplereturns.jld", "stockvals", stocks.value,
    #     "marketvals", market.value)
end

# TODO: Think about how I should have tested the initialisation
# TODO: Set test=true flag for testing bounds and other such things
#  or is that bad practice?
# TODO: Create logging set-up that records parameters, plots and outputs
# TODO: Think about how I should have tested the running cycle
# TODO: Write functions that collect activity of a) selling out and b) buying in
# TODO: Tests for cases with empty reviewers, divestments
