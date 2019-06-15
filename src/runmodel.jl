using Test, Random, LinearAlgebra, StatsBase, StatsPlots, Distributions, Traceur

using .Types, .Func, .Params

function runmodel()

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

    # TODO: Create logging set-up that records parameters, plots and outputs
    # TODO: Think about how I should have tested the running cycle
    # TODO: Write functions that collect activity of a) selling out and b) buying in
    # TODO: Tests for cases with empty reviewers, divestments
    # QUESTION: Need 'if reviewers not empty'/'divestments not empty' branches?
    # Looping over time after the end of the  perfwindow[end]+1 initialisation phase

    Func.modelrun(market, stocks, investors, funds)
    Func.boundstest(market, stocks, investors, funds)
    #Func.plot_stylisedfacts(market.value, stocks.value, stocks.volume)
end
