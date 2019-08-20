using Test, Random, LinearAlgebra, StatsBase, StatsPlots, Distributions
using Traceur, Parameters, Logging, Dates

using .Types, .Func, .Params

function runmodel(params=Params.default(); fundselector="probabilistic",
        doplot=false, boundstest=false)
    @unpack bigk, bigm, bign, bigt, impactrange, perfwindow, plotpath = params
    println("Run for ", bigt, " periods.")
    println("Sizes: ", bign, " investors, ", bigm, " funds, ", bigk, " stocks.")
    println("Stock impact range:", impactrange, ".")
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

    if boundstest
        Func.boundstest(market, stocks, investors, funds)
    end

    Func.modelrun(market, stocks, investors, funds, params, fundselector)

    if boundstest
        Func.boundstest(market, stocks, investors, funds)
    end

    stylefacts = Func.calc_stylisedfacts(market.value, stocks.value,
        stocks.volume, params)
    println(typeof(stylefacts))

    if doplot
        Func.plot_stylisedfacts(market.value, stocks.value, stylefacts, params)
        io = open(joinpath(plotpath, "parameters.txt"), "w+")
        logger = SimpleLogger(io)
        with_logger(logger) do
            @info "Parameters" Dates.today() Dates.Time(Dates.now()) fundselector params "zScoreKurtoses" stylefacts[10] "zScoreLossGain" stylefacts[11] "zScoreVolumeVola" stylefacts[12] "resultKStest" stylefacts[13]
        end
        close(io)
    end
end # runmodel()

# TODO: Think about how I should have tested the initialisation
# TODO: Think about how I should have tested the running cycle
# TODO: Tests for cases with empty reviewers, divestments
