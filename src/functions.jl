module Func
using Random, Test, StatsBase, ProgressMeter, StatsPlots, Distributions

include("types.jl")
include("params.jl")
import .Types, .Params

# NOTE Used to have default parameters: drift=Params.drift, marketvol=Params.marketvol. Consider these for the functions in general, so that function arguments are minimal outside of testing.

# TODO: Add description to each function, according to some standard format
# that should probably inclucde input, output and purpose, plus comments on
# assumptions or non-obvious points.

function marketmove!(
    marketvals,
    t::Int,
    drift=Params.drift,
    marketvol=Params.marketvol)

    marketvals[t] = marketvals[t-1]*(1 + drift) + marketvol * randn()

    return marketvals
end

function marketinit!(
    marketvals,
    marketstartval=Params.marketstartval,
    perfwindow=Params.perfwindow,
    drift=Params.drift,
    marketvol=Params.marketvol)

    horizon = perfwindow[end]
    marketvals[1] = marketstartval
    for t in 2:horizon+1
        Func.marketmove!(marketvals, t, drift, marketvol)
    end

    return marketvals
end

function stockmove!(
    stockvals,
    t::Int,
    mktvals,
    betas,
    stockvolas)

    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    stockvals[:, t] = ((1 .+ mktreturn .* betas) .* stockvals[:, t-1]) +
        stockvolas .* randn(length(stockvolas)) .* stockvals[:, t-1]

    return stockvals
end

function betainit!(
    betas,
    bigm=Params.bigm,
    betastd=Params.betastd,
    betamean=Params.betamean)

    betas .= randn(bigm) .* betastd .+ betamean

    return betas
end

# Function for a discrete volatility range
function stockvolinit!(
    stockvol,
    volrange=Params.stockvolrange,
    bigm=Params.bigm)

    stockvol .= rand(volrange, bigm)

    return stockvol
end

function stockvalueinit!(
    stocks,
    marketvals,
    stockstartval=Params.stockstartval,
    perfwindow=Params.perfwindow)

    horizon = perfwindow[end]
    stocks.value[:,1] .= stockstartval
    for t in 2:(horizon + 1)
        stocks.value .= Func.stockmove!(
        stocks.value, t, marketvals, stocks.beta, stocks.vol)
    end

    return stocks.value
end

function stockimpactinit!(
    impact,
    impactrange=Params.impactrange)

    impact .= rand(impactrange, length(impact))

    return impact
end

function invhorizoninit!(
    horizons,
    perfwindow=Params.perfwindow)
    horizons .= rand(perfwindow, length(horizons))

    return horizons
end

function invthreshinit!(
    thresholds,
    mean=Params.thresholdmean,
    std=Params.thresholdstd)

    thresholds .= randn(length(thresholds)) * std .+ mean

    return thresholds
end

# TODO: Sort out the capital range to match other items, as in rand(range, n)
function invassetinit!(
    invassets,
    caprange=Params.invcaprange,
    bigk=Params.bigk)

    ninvs = size(invassets, 1)
    capital = rand(caprange[1]:caprange[2], ninvs)
    for inv in 1:bigk
        invassets[inv, inv] = capital[inv]
    end
    for inv in (bigk+1):ninvs
        invassets[inv, rand(1:bigk)] = capital[inv]
    end

    return invassets
end

function fundcapitalinit!(
    fundvals,
    investorassets)

    fundvals[:,1] .= sum(investorassets, dims=1)[1:end-1]

    return fundvals
end

# TODO: Test of this function didn't throw error when ncols was misspecified
# as ncols=size(investorassets,1)- 1, i.e. wrong axis. Check why, write another
function fundstakeinit!(
    fundstakes,
    investorassets)

    ncols = size(investorassets, 2) - 1
    for col in 1:ncols
        fundstakes[col,:] = investorassets[:,col] ./ sum(investorassets[:,col])
    end

    return fundstakes
end

function fundholdinit!(
    holdings,
    initialcapital,
    stockvals,
    portfsizerange=Params.portfsizerange)

    kfunds = size(holdings, 1)
    mstocks = size(holdings, 2)
    nofm = rand(portfsizerange, kfunds) # Number of stocks in each portfolio
    for k in 1:kfunds # Loop over funds
        selection = rand(1:mstocks, nofm[k]) # Randomly select w/ replacement
        for stock in selection # Loop over stocks
            holdings[k, stock] +=
            (initialcapital[k] / length(selection)) / stockvals[stock]
        end
    end

    return holdings
end

function fundvalinit!(
    fundvals,
    holdings,
    stockvals,
    perfwindow=Params.perfwindow)

    horizon = perfwindow[end]
    kfunds = size(fundvals, 1)
    for t in 2:(horizon + 1) # t=1 is initialised by fundcapitalinit!
        for k in 1:kfunds
            fundvals[k, t] = sum(holdings[k, :] .* stockvals[:, t])
        end
    end

    return fundvals
end

function drawreviewers(
    bign=Params.bign,
    reviewprob=Params.reviewprobability)

    reviewers = rand(bign) .< reviewprob

    return findall(reviewers)
end

# QUESTION: Is it sufficient to compare fund values or do we need to compare the
# value of investors' stakes over time?

# TODO: Test that only a single investors.asset is > 0
function perfreview(
    t,
    reviewers,
    investors,
    fundvals)

    divestments = zeros(Int64, 1,2)
    for rev in reviewers
        horizon = investors.horizon[rev]
        fnds = findall(investors.assets[rev, 1:(end-1)] .> 0)
        for fnd in fnds
            valchange =
            (fundvals[fnd,t] - fundvals[fnd, t-horizon]) / fundvals[fnd,horizon]
            if valchange < investors.threshold[rev]
                reviewerfundpair = [convert(Int64, rev) convert(Int64, fnd)]
                if sum(divestments[1,:]) == 0
                    divestments = reviewerfundpair
                else
                    divestments = [divestments; reviewerfundpair]
                end
            end
        end
    end

    return divestments
end

# TODO: Write test that checks that sellorders and buyorders have the same
# number of rows in .values and .investors/.funds
# QUESTION: Should I implement size and value checks for all intermediate outpts
# NOTE: Passes one column of stock values only

function liquidate!(
    holdings,
    stakes,
    stockvals,
    divestments)

    sellorders = Types.SellMarketOrder(
    Array{Float64}(undef, 0, size(holdings, 2)), Array{Float64}(undef, 0))
    for row in 1:size(divestments, 1)
        investor = divestments[row, 1]
        fund = divestments[row, 2]

        # Note the minus
        salevals = -stakes[fund, investor] .* (holdings[fund, :] .* stockvals)'
        sellorders.values = vcat(sellorders.values, salevals)
        sellorders.investors = vcat(sellorders.investors, investor)

        remainder = (1 - stakes[fund, investor])
        holdings[fund, :] = holdings[fund, :] .* remainder

        stakes[fund, investor] = 0
        if sum(stakes[fund,: ]) > 0
            stakes[fund, :] = stakes[fund, :] ./ sum(stakes[fund,: ])
        end
    end

    return sellorders, holdings, stakes
end

# NOTE: Truncate minimum stock value at 0.000000001
# NOTE: Could move volume calculation into marketmake
function marketmake!(
    stocks,
    t,
    ordervals)

    netimpact = sum(ordervals, dims=1)' .* stocks.impact
    stocks.value[:, t] .= vec((1 .+ netimpact) .* stocks.value[:, t])

    nearzeroidx = findall(stocks.value[:, t] .< 0.000000001)
    stocks.value[nearzeroidx, t] .= 0.000000001

    return stocks.value
end

function executeorder!(
    stocks,
    t,
    orders::Types.SellMarketOrder)

    cashout = Array{Float64}(undef, 0, 2)
    oldstockvals = copy(stocks.value)
    stocks.value .= marketmake!(stocks, t, orders.values)
    priceratio = stocks.value[:, t] ./ oldstockvals[:, t]
    for order in 1:size(orders.values, 1)
        investor = orders.investors[order]
        amount = -sum(orders.values[order, :] .* priceratio)
        cashout = vcat(cashout, [investor amount])
    end

    return stocks.value, cashout
end

# NOTE: Small issue here that sell orders happen first, so sellers get a bad
# price, followed by buyers who get a good price, ~midmarket as their order
# neutralises the prior decline. However, the amount they buy is also diminished
function executeorder!(
    stocks,
    t,
    orders::Types.BuyMarketOrder)

    mstocks = size(stocks.value, 1)
    sharesout = Array{Float64}(undef, 0, mstocks + 1)
    stocks.value .= marketmake!(stocks, t, orders.values)
    for order in 1:size(orders.values, 1)
        fund = orders.funds[order]
        shares = (orders.values[order, :] ./ stocks.value[:, t])'
        sharesout = vcat(sharesout, [fund shares])
    end

    return stocks.value, sharesout
end

# QUESTION: Structure of disburse! methods is not symmetric in args... fix?
# Disburse cash to investors following  sell order / divestment
function disburse!(
    invassets,
    divestments,
    cashout)

    for row in 1:size(divestments,1)
        inv = divestments[row, 1]
        fund = divestments[row, 2]
        #@assert(inv == cashout[row, 1]) # Ensure divestment matches cashout
        invassets[inv, fund] = 0
        invassets[inv, end] = cashout[row, 2]
    end

    return invassets
end

# Disburse shares to funds following buy order/investment, update value history
# QUESTION: Why doesn't the value of funds drop as investors divest?
function disburse!(
    funds::Types.EquityFund,
    sharesout,
    stockvals)

    for row in 1:size(sharesout, 1)
        fund = convert(Int64, sharesout[row, 1])
        funds.holdings[fund, :] += sharesout[row, 2:end]
        fundrevalue!(funds, stockvals, [fund])
    end

    return funds
end

function fundrevalue!(
    funds,
    stockvals,
    targets=1:Params.bigk)

    for k in targets
        funds.value[k, :] .= vec(funds.holdings[k, :]' * stockvals)
    end

    return funds.value
end

function respawn!(
    funds,
    investors,
    t,
    stockvals,
    portfsizerange=Params.portfsizerange)

    buyorders = Types.BuyMarketOrder(
    Array{Float64}(undef, 0, size(funds.holdings, 2)), Array{Float64}(undef, 0))
    spawn = findall(vec(sum(funds.holdings, dims=2) .== 0))
    neggs = length(spawn)
    potentialinvs = findall(vec(investors.assets[:, end] .> 0))
    anchorinvs = Random.shuffle(potentialinvs)[1:neggs]
    for idx in 1:neggs
        inv = anchorinvs[idx]
        egg = spawn[idx]
        capital = investors.assets[inv, end]
        # Anchor investor notes investment in the fund
        investors.assets[inv, egg] = capital
        investors.assets[inv, end] = 0
        # Re-use fund holdings init fctn to draw new stocks and generate order
        buyorder = (fundholdinit!(funds.holdings[egg, :]',
        capital, stockvals[:, t], portfsizerange)' .* stockvals[:, t])'
        buyorders.values = vcat(buyorders.values, buyorder)
        buyorders.funds = vcat(buyorders.funds, egg)
        # Set anchor investor's stake
        funds.stakes[egg, inv] = 1
    end

    return buyorders, investors.assets, funds.stakes
end

# NOTE: Bit of an issue: as many reinvestments as dead funds are randomly drawn,
# which reduces the flow-performance relationship, especially if many funds die

function reinvest!(
    investors,
    funds,
    stockvals,
    t,
    selection="best")

    kfunds = size(investors.assets, 2) - 1

    buyorders = Types.BuyMarketOrder(
    Array{Float64}(undef, 0, size(funds.holdings, 2)), Array{Float64}(undef, 0))

    reinvestors = findall(vec(investors.assets[:, end] .> 0))
    for reinv in reinvestors
        injection = investors.assets[reinv, end]
        if selection == "best"
            # Choose best-performing fund
            choice = bestperformer(funds.value, investors.horizon[reinv], t)
        elseif selection == "goodenough"
            # Choose any fund above return threshold
            choice = goodenoughfund(
            funds.value, investors.horizon[reinv], investors.threshold[reinv],t)
        else
            println("Please specify a valid fund selection method in reinvest!")
        end
        # Move cash into the fund
        investors.assets[reinv, choice] = injection
        investors.assets[reinv, end] = 0

        fundval = funds.holdings[choice, :]' * stockvals[:,t]
        # QUESTION: Should assert hold true, i.e. funds.value be unchanged?
        #@assert fundval == funds.value[choice, t]
        stakevals = funds.stakes[choice, :] .* fundval
        stakevals[reinv] = injection
        newcapital = fundval + injection
        funds.stakes[choice, :] .= stakevals ./ newcapital

        portfoliovalues = funds.holdings[choice, :] .* stockvals[:, t]
        buyorder = (portfoliovalues ./ fundval)' .* injection
        buyorders.values = vcat(buyorders.values, buyorder)
        buyorders.funds = vcat(buyorders.funds, choice)
    end

    return investors.assets, funds.stakes, buyorders
end

function bestperformer(
    fundvals,
    horizon,
    t)

    # Return between investor's horizon and now/t
    horizonreturns =
    (fundvals[:, t] - fundvals[:, t-horizon]) ./ fundvals[:, t-horizon]
    # Index of best-performing fund
    _, bestfund = findmax(horizonreturns)

    return bestfund
end

# TODO: Test for goodenoughfund function
# NOTE: Choose best performer if no funds meet performance threshold
function goodenoughfund(
    fundvals,
    horizon,
    threshold,
    t)

    # Return between investor's horizon and now/t
    horizonreturns =
    (fundvals[:, t] - fundvals[:, t-horizon]) ./ fundvals[:, t-horizon]
    # Index of best-performing fund
    candidates = findall(horizonreturns .> threshold)
    if size(candidates,1) == 0
        _, goodenough = findmax(horizonreturns)
    else
        goodenough = rand(1:size(fundvals,1))
    end

    return goodenough
end

function initialise(market, stocks, investors, funds)
    Func.marketinit!(market.value)
    Func.betainit!(stocks.beta)
    Func.stockvolinit!(stocks.vol)
    Func.stockvalueinit!(stocks, market.value)
    Func.stockimpactinit!(stocks.impact)
    Func.invhorizoninit!(investors.horizon)
    Func.invthreshinit!(investors.threshold)
    Func.invassetinit!(investors.assets)
    Func.fundcapitalinit!(funds.value, investors.assets)
    Func.fundstakeinit!(funds.stakes, investors.assets)
    Func.fundholdinit!(funds.holdings, funds.value[:, 1], stocks.value)
    Func.fundvalinit!(funds.value, funds.holdings, stocks.value)
end

function boundstest(market, stocks, investors, funds)
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
    end
end

function modelrun(market, stocks, investors, funds)
    @showprogress for t in (Params.perfwindow[end]+2):Params.bigt

        # Market index and assets move
        Func.marketmove!(market.value, t)
        # Func.boundstest(market, stocks, investors, funds)
        Func.stockmove!(stocks.value, t, market.value, stocks.beta, stocks.vol)
        # Func.boundstest(market, stocks, investors, funds)
        Func.fundrevalue!(funds, stocks.value)
        # Func.boundstest(market, stocks, investors, funds)
        # print("Time: ", t)
        # QUESTION: Consider whether I should pass t instead of stocks.value[:, t]
        # Investors (may_ sell out of underperforming funds
        # TODO: Test sellorders.investors == cashout[1,:], similar for buyorders
        # TODO: No cashout of value 0
        # QUESTION: Enforce integers for fund/investor IDs in sellorders, cashout..?
        # QUESTION: Is a separate trackvolume! function appropriate?

        # Investors' performance review
        reviewers = Func.drawreviewers()
        divestments = Func.perfreview(t, reviewers, investors, funds.value)

        # Run liquidation/reinvestment cycle if there are divestments
        if sum(divestments[1,:]) > 0

            sellorders, _, _ = Func.liquidate!(funds.holdings, funds.stakes,
            stocks.value[:, t], divestments)

            Func.trackvolume!(stocks.volume, sellorders, t)

            _, cashout = Func.executeorder!(stocks, t, sellorders)
            Func.disburse!(investors.assets, divestments, cashout)

            # TODO: Write test that check correct outocme of the whole process,
            # step by step

            # Generate a new fund for each dead one, if there are dead ones
            if any(sum(funds.stakes, dims=2) .== 0)

                buyorders, _, _ = Func.respawn!(funds, investors, t, stocks.value)

                Func.trackvolume!(stocks.volume, sellorders, t)

                _, sharesout = Func.executeorder!(stocks, t, buyorders)

                Func.disburse!(funds, sharesout, stocks.value)

            end # respawn loop

            # Re-value funds prior to reinvestment decisions
            Func.fundrevalue!(funds, stocks.value)

            # If there are still investors with spare cash, they reinvest it
            if any(investors.assets[:, end] .> 0)

                _, _, buyorders = Func.reinvest!(
                investors, funds, stocks.value, t, "goodenough")

                Func.trackvolume!(stocks.volume, sellorders, t)

                _, sharesout = Func.executeorder!(stocks, t, buyorders)

                Func.disburse!(funds, sharesout, stocks.value)
            end # reinvestment loop

        end # Divestment loop
    end # Model run loop
end # runmodel function


# TODO: Write tests for both assetreturns functions and vecreturns
# NOTE: Duplication in assetreturns and demean, i.e. targets for refactoring
function vecreturns(vector::Array{Float64,1})
    returnsvec = (vector[2:end] - vector[1:(end-1)]) ./ vector[1:(end-1)]
    return returnsvec
end

function assetreturns(vector::Array{Float64,1})
    returnsvec = vecreturns(vector)
    return returnsvec
end

function assetreturns(array::Array{Float64,2})
    nassets = size(array, 1)
    nperiods = size(array, 2)
    returnsarray = zeros(nassets, nperiods-1)
    for asset in 1:nassets
        returnsarray[asset,:] .= vecreturns(array[asset,:])
    end
    return returnsarray
end

function subtractmean(returns::Array{Float64,1})
    demeaned = returns .- mean(returns)
    return demeaned
end

# TODO: Test demean functions
function demean(returns::Array{Float64,1})
    demeanedreturns = subtractmean(returns)
    return demeanedreturns
end

function demean(returns::Array{Float64,2})
    nassets = size(returns, 1)
    nperiods = size(returns, 2)
    demeanedreturns = zeros(nassets, nperiods)
    for asset in 1:nassets
        demeanedreturns[asset, :] .= subtractmean(returns[asset, :])
    end
    return demeanedreturns
end
function calckurtoses(returns, startidx=Params.perfwindow[end]+1)
    nassets = size(returns, 1)
    kurtoses = zeros(nassets)
    for asset in 1:nassets
        kurtoses[asset] = kurtosis(returns[asset, startidx:end])
    end
    return kurtoses
end

# TODO: Test for trackvolume!
# NOTE: Stylised facts should only be calculated for period after initialisation
function trackvolume!(volume, order::Types.BuyMarketOrder, t)
    volume[:, t] += sum(order.values, dims=1)'
    return volume
end

function trackvolume!(volume, order::Types.SellMarketOrder, t)
    volume[:, t] -= sum(order.values, dims=1)'
    return volume
end

# TODO: Test for volavolumecorr
function volavolumecorr(volumes, returns, startidx=Params.perfwindow[end]+1)
    nassets = size(volumes, 1)
    correlation = zeros(nassets)
    for asset in 1:nassets
        correlation[asset] =
        StatsBase.cor(volumes[asset, (startidx+1):end],
        broadcast(abs, returns[asset, startidx:end]))
    end
    return correlation
end

# TODO: Test for lossgainratio
function lossgainratio(returns, cutoff, startidx=Params.perfwindow[end]+1)
    nassets = size(returns, 1)
    ratios = zeros(nassets)
    for asset in 1:nassets
        stock = returns[asset, startidx:end]
        absolutevals = broadcast(abs, stock)
        largeamount = percentile(absolutevals, cutoff)
        upmoves = sum(stock .> largeamount)
        downmoves = sum(stock .< -largeamount)
        ratios[asset] = downmoves/upmoves
    end
    return ratios
end

# TODO: move plots to separate file
# TODO: save plots in plots file

function plot_stylisedfacts(
    marketval, stocksval, stocksvolume, returnspctile=5)
    marketreturns = Func.assetreturns(marketval)
    stockreturns = Func.assetreturns(stocksval)
    demeanedmarketreturns = Func.demean(marketreturns)
    demeanedstockreturns = Func.demean(stockreturns)
    demeanedmarketreturnssquared = demeanedmarketreturns.^2
    demeanedstockreturnssquared = demeanedstockreturns.^2
    lags = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]#collect[1:10]
    choice = rand(1:size(stockreturns,1))
    pyplot() # as GR wasn't working, maybe fix that
    plot_pricehistories(stocksval, marketval)
    plot_marketreturnhistogram(demeanedmarketreturns)
    plot_stockreturnhistogram(demeanedstockreturns, choice)
    plot_marketpacf(demeanedmarketreturns, lags)
    plot_stockpacf(demeanedstockreturns, choice, lags)
    plot_marketvolaclustering(demeanedmarketreturnssquared, lags)
    plot_stockvolaclustering(demeanedstockreturnssquared, choice, lags)
    plot_stockkurtoses(demeanedstockreturns)
    plot_lossgainratio(demeanedstockreturns, returnspctile)
    plot_volavolumecorr(stocksvolume, demeanedstockreturns)
end

# TODO: refactor plots - DRY
function plot_pricehistories(stocksval, marketval)
    periods = 1:length(marketval)
    nstocks = 1:100
    StatsPlots.plot(periods, transpose(stocksval)[:, nstocks], legend=:none,
    title="History of stock returns")
    plt = StatsPlots.plot!(
        periods, marketval, lw=3, lc=:black, label="Market index")
    png("returnshistory")
    display(plt)
end

function plot_marketreturnhistogram(asset)
    StatsPlots.histogram(asset, bins=100, title="Market returns vs normal distribution", label="Market returns", legend=:topright)
    x = findmin(asset)[1]:0.001:findmax(asset)[1]
    d = fit(Normal, asset)
    plt = StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
    png("marketreturnsvsnormal")
    display(plt)
end

function plot_stockreturnhistogram(demeanedstockreturns, choice)
    asset = demeanedstockreturns[choice, :]
    StatsPlots.histogram(asset, bins=100, title="Example of stock returns", label="Random stock's returns", legend=:topright)
    x = findmin(asset)[1]:0.001:findmax(asset)[1]
    d = fit(Normal, asset)
    plt = StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
    png("stockreturnsvsnormal")
    display(plt)
end

function plot_marketpacf(demeanedmarketreturns, lags)
    asset = demeanedmarketreturns
    pacfcoeffs = pacf(demeanedmarketreturns, lags)
    StatsPlots.bar(pacfcoeffs, title="Market autocorrelation", label="Partial correlation coefficients")
    plt = StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)),], lc=:red, label="critical values", legend=:bottomright)
    png("marketautocorrelation")
    display(plt)
end

function plot_stockpacf(demeanedstockreturns, choice, lags)
    asset = demeanedstockreturns[choice,: ]
    pacfcoeffs = pacf(asset, lags)
    StatsPlots.bar(pacfcoeffs, title="Stock autocorrelation (random stock example)", label="Partial autocorrelation coefficients, demeaned returns", legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)),], lc=:red, label="Critical values")
    png("stockautocorrelation")
    display(plt)
end

function plot_marketvolaclustering(demeanedmarketreturnssquared, lags)
    asset = demeanedmarketreturnssquared
    pacfcoeffs = pacf(asset, lags)
    StatsPlots.bar(pacfcoeffs, title="Volatility clustering, market", label="Partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)),], lc=:red, label="Critical values")
    png("marketvolaclustering")
    display(plt)
end

function plot_stockvolaclustering(demeanedstockreturnssquared, choice, lags)
    asset = demeanedstockreturnssquared[choice,:]
    pacfcoeffs = pacf(asset, lags)
    StatsPlots.bar(pacfcoeffs, title="Volatility clustering, random stock", label="partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10, [ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)), -ones(length(0:0.01:10)) .* 1.96 / sqrt(size(asset, 1)),], lc=:red, label="Critical values")
    png("stockvolaclustering")
    display(plt)
end

function plot_stockkurtoses(demeanedstockreturns)
    stockkurtoses = Func.calckurtoses(demeanedstockreturns)
    plt = StatsPlots.bar(sort!(stockkurtoses, rev=true), title="Kurtoses of stock returns", label="Kurtoses")
    png("kurtoses")
    display(plt)
end

function plot_lossgainratio(demeanedstockreturns, returnspercentile)
    ratios = Func.lossgainratio(demeanedstockreturns, returnspercentile)
    plot(ratios, title="Gain-loss asymmetry in stock returns", label="# large losses / # large gains")
    plt = plot!(0:0.01:200, ones(length(0:0.01:200)) * mean(ratios), label="Mean ratio")
    png("lossgainratio")
    display(plt)
end

function plot_volavolumecorr(volume, demeanedstockreturns)
    corrs = Func.volavolumecorr(volume, demeanedstockreturns)
    plt = bar(sort(corrs, rev=true), title="Volume-volatility correlation of stocks", label="Volume - absolute return correlation coefficients")
    png("volavolumecorr")
    display(plt)
end


end  # module Functions
