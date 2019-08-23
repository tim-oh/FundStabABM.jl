using Test, Random, LinearAlgebra, StatsBase, StatsPlots, Distributions
using Traceur, Parameters, Logging, Dates, GLM

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

    liquidationlog = Func.modelrun(market, stocks, investors, funds, params,
        fundselector)
    println(liquidationlog)

    if boundstest
        Func.boundstest(market, stocks, investors, funds)
    end

    logit = Func.logisticregression(liquidationlog)
    println("\nLogistic regression results with whitened indep vars\n", logit)

    stylefacts = Func.calc_stylisedfacts(market.value, stocks.value,
        stocks.volume, params)

    if doplot
        Func.plot_stylisedfacts(market.value, stocks.value, stylefacts, params,
            liquidationlog)
        io = open(joinpath(plotpath, "parameters.txt"), "w+")
        logger = SimpleLogger(io)
        with_logger(logger) do
            @info "Parameters" Dates.today() Dates.Time(Dates.now()) fundselector params
            @info "zScoreKurtoses" stylefacts[10]
            @info "zScoreLossGain" stylefacts[11]
            @info "zScoreVolumeVola" stylefacts[12]
            @info "resultKStest" stylefacts[13]
            @info "logit" logit
            @info "liquidationslog" liquidationlog
        end
        close(io)
    end
end # runmodel()

# TODO: Think about how I should have tested the initialisation
# TODO: Think about how I should have tested the running cycle
# TODO: Tests for cases with empty reviewers, divestments
