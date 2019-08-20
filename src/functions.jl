module Func
using Random, Test, StatsBase, ProgressMeter, StatsPlots, Distributions
using Base.Iterators, Parameters, HypothesisTests
using PyPlot

include("types.jl")
import .Types

# TODO: Write docstrings to be used with Documenter.jl.

# TODO: Change bracket, marketvol*randn() should be inside the bracket!

function marketmove!(marketvals, t::Int, params)
    @unpack drift, marketvol = params
    marketvals[t] = marketvals[t-1]*(1 + drift + marketvol * randn())
    return marketvals
end

function marketinit!(marketvals, params)
    @unpack marketstartval, perfwindow, drift, marketvol = params
    horizon = perfwindow[end]
    marketvals[1] = marketstartval
    for t in 2:horizon+1
        Func.marketmove!(marketvals, t, params)
    end
    return marketvals
end

function stockmove!(stockvals, t::Int, mktvals, betas, stockvolas)
    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    stockvals[:, t] = ((1 .+ mktreturn .* betas) .* stockvals[:, t-1]) +
        stockvolas .* randn(length(stockvolas)) .* stockvals[:, t-1]
    return stockvals
end

# Keeping old function around - it's 30% faster
# function betainit!(
#     betas,
#     bigm=Params.bigm,
#     betastd=Params.betastd,
#     betamean=Params.betamean)
#
#     betas .= randn(bigm) .* betastd .+ betamean
#
#     return betas
# end

function betainit!(betas, params)
    @unpack bigm, betastd, betamean  = params
    betas .= randn(bigm) .* betastd .+ betamean
    return betas
end

# Function for a discrete volatility range
function stockvolinit!(stockvol, params)
    @unpack stockvolrange, bigm = params
    stockvol .= rand(stockvolrange, bigm)
    return stockvol
end

function stockvalueinit!(stocks, marketvals, params)
    @unpack stockstartval, perfwindow = params
    horizon = perfwindow[end]
    stocks.value[:,1] .= stockstartval
    for t in 2:(horizon + 1)
        stocks.value .= Func.stockmove!(
        stocks.value, t, marketvals, stocks.beta, stocks.vol)
    end
    return stocks.value
end

function stockimpactinit!(impact, params)
    @unpack impactrange = params
    impact .= rand(impactrange, length(impact))
    return impact
end

function invhorizoninit!(horizons, params)
    @unpack perfwindow = params
    horizons .= rand(perfwindow, length(horizons))
    return horizons
end

function invthreshinit!(thresholds, params)
    @unpack thresholdmean, thresholdstd = params
    thresholds .= randn(length(thresholds)) .* thresholdstd .+ thresholdmean
    return thresholds
end

# TODO: Sort out the capital range to match other items, as in rand(range, n)
function invassetinit!(invassets, params)
    @unpack invcaprange, bigk = params
    ninvs = size(invassets, 1)
    capital = rand(invcaprange[1]:invcaprange[2], ninvs)
    for inv in 1:bigk
        invassets[inv, inv] = capital[inv]
    end
    for inv in (bigk+1):ninvs
        invassets[inv, rand(1:bigk)] = capital[inv]
    end
    return invassets
end

function fundcapitalinit!(fundvals, investorassets)
    fundvals[:,1] .= sum(investorassets, dims=1)[1:end-1]
    return fundvals
end

# TODO: Test of this function didn't throw error when ncols was misspecified
# as ncols=size(investorassets,1)- 1, i.e. wrong axis. Check why, write another
function fundstakeinit!(fundstakes, investorassets)
    ncols = size(investorassets, 2) - 1
    for col in 1:ncols
        fundstakes[col,:] = investorassets[:,col] ./ sum(investorassets[:,col])
    end
    return fundstakes
end

function fundholdinit!(holdings, initialcapital, stockvals, params)
    @unpack portfsizerange = params
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

function fundvalinit!(fundvals, holdings, stockvals, params)
    @unpack perfwindow = params
    horizon = perfwindow[end]
    kfunds = size(fundvals, 1)
    for t in 2:(horizon + 1) # t=1 is initialised by fundcapitalinit!
        for k in 1:kfunds
            fundvals[k, t] = sum(holdings[k, :] .* stockvals[:, t])
        end
    end
    return fundvals
end

function drawreviewers(params)
    @unpack bign, reviewprobability = params
    reviewers = rand(bign) .< reviewprobability
    return findall(reviewers)
end

# QUESTION: Is it sufficient to compare fund values or do we need to compare the
# value of investors' stakes over time?

# TODO: Test that only a single investors.asset is > 0
function perfreview(t, reviewers, investors, fundvals)
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

function liquidate!(holdings, stakes, stockvals, divestments)
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
function marketmake!(stocks, t, ordervals)
    netimpact = sum(ordervals, dims=1)' .* stocks.impact
    stocks.value[:, t] .= vec((1 .+ netimpact) .* stocks.value[:, t])
    nearzeroidx = findall(stocks.value[:, t] .< 0.000000001)
    stocks.value[nearzeroidx, t] .= 0.000000001
    return stocks.value
end

function executeorder!(stocks, t,  orders::Types.SellMarketOrder)
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
function executeorder!(stocks, t, orders::Types.BuyMarketOrder)
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
function disburse!(invassets, divestments, cashout)
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
function disburse!(funds::Types.EquityFund, sharesout, stockvals)
    for row in 1:size(sharesout, 1)
        fund = convert(Int64, sharesout[row, 1])
        funds.holdings[fund, :] += sharesout[row, 2:end]
        fundrevalue!(funds, stockvals, [fund])
    end
    return funds
end

# TODO: Didn't error when I changed last argument from range to variable
function fundrevalue!(funds, stockvals, targets)
    for k in targets
        funds.value[k, :] .= vec(funds.holdings[k, :]' * stockvals)
    end
    return funds.value
end

# TODO: respawn test should have thrown error when I ran it without updating
# inputs to fundholdinit!
function respawn!(funds, investors, t, stockvals, params)
    @unpack portfsizerange = params
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
        capital, stockvals[:, t], params)' .* stockvals[:, t])'
        buyorders.values = vcat(buyorders.values, buyorder)
        buyorders.funds = vcat(buyorders.funds, egg)
        # Set anchor investor's stake
        funds.stakes[egg, inv] = 1
    end
    return buyorders, investors.assets, funds.stakes
end

# NOTE: Bit of an issue: as many reinvestments as dead funds are randomly drawn,
# which reduces the flow-performance relationship, especially if many funds die

function reinvest!(investors, funds, stockvals, t, selection="probabilistic")
    kfunds = size(investors.assets, 2) - 1
    buyorders = Types.BuyMarketOrder(
    Array{Float64}(undef, 0, size(funds.holdings, 2)), Array{Float64}(undef, 0))
    reinvestors = findall(vec(investors.assets[:, end] .> 0))
    for reinv in reinvestors
        injection = investors.assets[reinv, end]
        horizonreturns = horizonreturncalc(funds.value,
            investors.horizon[reinv], t)
        if selection == "best"
            # Choose best-performing fund
            choice = bestperformer(horizonreturns)
        elseif selection == "goodenough"
            # Choose any fund above return threshold
            choice = goodenoughfund(horizonreturns, investors.threshold[reinv])
        elseif selection == "probabilistic"
            choice = probabilisticchoice(horizonreturns)
        elseif selection == "random"
            choice = randomchoice(horizonreturns)
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

function horizonreturncalc(fundvals, horizon, t)
    horizonreturns =
        (fundvals[:, t] - fundvals[:, t-horizon]) ./ fundvals[:, t-horizon]
    return horizonreturns
end

function bestperformer(horizonreturns)
    _, bestfund = findmax(horizonreturns)
    return bestfund
end

# TODO: Test for goodenoughfund function
function goodenoughfund(horizonreturns, threshold)
    candidates = findall(horizonreturns .> threshold)
    if size(candidates,1) == 0 # Choose best if none meet threshold
        _, goodenough = findmax(horizonreturns)
    else # Otherwise choose one at random from those that meet threshold
        goodenough = rand(candidates)
    end
    return goodenough
end

# NOTE: Could rationalise many in-loop operations by computation outside and
# passing as an argument
# TODO: Write a test for this function
function probabilisticchoice(horizonreturns)
    probability = similar(horizonreturns)
    indices = 1:length(horizonreturns)
    normaliser = sum(indices)
    returnrankorder = sortperm(horizonreturns) # Return indices, ascending order
    probability[returnrankorder] = collect(indices) / normaliser
    @assert sum(probability) ≈ 1
    distn = Categorical(probability)
    fundchoice = rand(distn)
    return fundchoice
end

function randomchoice(horizonreturns)
    fundchoice = rand(1:length(horizonreturns))
    return fundchoice
end

function initialise(market, stocks, investors, funds, params)
    Func.marketinit!(market.value, params)
    Func.betainit!(stocks.beta, params)
    Func.stockvolinit!(stocks.vol, params)
    Func.stockvalueinit!(stocks, market.value, params)
    Func.stockimpactinit!(stocks.impact, params)
    Func.invhorizoninit!(investors.horizon, params)
    Func.invthreshinit!(investors.threshold, params)
    Func.invassetinit!(investors.assets, params)
    Func.fundcapitalinit!(funds.value, investors.assets)
    Func.fundstakeinit!(funds.stakes, investors.assets)
    Func.fundholdinit!(funds.holdings, funds.value[:, 1], stocks.value, params)
    Func.fundvalinit!(funds.value, funds.holdings, stocks.value, params)
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

function modelrun(market, stocks, investors, funds, params,
        fundselector="bestperformer")
    @unpack perfwindow, bigt, bigk = params
    @showprogress for t in perfwindow[end]+2:bigt
        Func.marketmove!(market.value, t, params)
        Func.stockmove!(stocks.value, t, market.value, stocks.beta, stocks.vol)
        Func.fundrevalue!(funds, stocks.value, 1:bigk)
        # TODO: Test sellorders.investors == cashout[1,:], similar for buyorders
        reviewers = Func.drawreviewers(params)
        divestments = Func.perfreview(t, reviewers, investors, funds.value)
        if sum(divestments[1,:]) > 0
            sellorders, _, _ = Func.liquidate!(funds.holdings, funds.stakes,
            stocks.value[:, t], divestments)
            Func.trackvolume!(stocks.volume, sellorders, t)
            _, cashout = Func.executeorder!(stocks, t, sellorders)
            Func.disburse!(investors.assets, divestments, cashout)
            # TODO: Write test that check each step of modelrun()
            if any(sum(funds.stakes, dims=2) .== 0)
                buyorders, _, _ = Func.respawn!(funds, investors, t,
                    stocks.value, params)
                Func.trackvolume!(stocks.volume, sellorders, t)
                _, sharesout = Func.executeorder!(stocks, t, buyorders)
                Func.disburse!(funds, sharesout, stocks.value)
            end # respawn loop
            # Re-value funds prior to reinvestment decisions
            Func.fundrevalue!(funds, stocks.value, 1:bigk)
            # If there are still investors with spare cash, they reinvest it
            if any(investors.assets[:, end] .> 0)
                _, _, buyorders = Func.reinvest!(investors, funds, stocks.value,
                    t, fundselector)
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

function assetreturns(vector::Array{Float64,1}, perfwindow)
    startidx = perfwindow[end] + 1
    returnsvec = vecreturns(vector)[startidx:end]
    return returnsvec
end

function assetreturns(array::Array{Float64,2}, perfwindow)
    startidx = perfwindow[end] + 1
    nassets = size(array, 1)
    nperiods = size(array, 2) - perfwindow[end] -1
    returnsarray = zeros(nassets, nperiods)
    for asset in 1:nassets
        returnsarray[asset,:] .= vecreturns(array[asset, startidx:end])
    end
    return returnsarray
end

function subtractmean(returns::Array{Float64,1})
    demeaned = returns .- mean(returns)
    return demeaned
end

# TODO: Test demean functions
function demean(returns::Vector{Float64})
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

function calc_kurtoses(returns)
    bigm = size(returns, 1)
    kurtoses = zeros(bigm)
    for asset in 1:bigm
        kurtoses[asset] = kurtosis(returns[asset, :])
    end
    return kurtoses
end

# TODO: Test for trackvolume!
# NOTE: Stylised facts should only be calculated for period after initialisation
function trackvolume!(volume, order::Types.BuyMarketOrder, t)
    volume[:, t] += sum(order.values, dims=1)'
    return volume
end

# TODO: Test for trackvolume!
function trackvolume!(volume, order::Types.SellMarketOrder, t)
    volume[:, t] -= sum(order.values, dims=1)'
    return volume
end

# TODO: Test for volavolumecorr
function volavolumecorr(volumes, returns)
    bigm = size(returns, 1)
    correlation = zeros(bigm)
    for asset in 1:bigm
        correlation[asset] =
        StatsBase.cor(volumes[asset, :],
        broadcast(abs, returns[asset, :]))
    end
    return correlation
end

# TODO: Test for lossgainratio
function lossgainratio(returns, cutoff, params)
    @unpack perfwindow, bigm = params
    startidx = perfwindow[end] + 1
    ratios = zeros(bigm)
    for asset in 1:bigm
        stock = returns[asset, startidx:end]
        absolutevals = broadcast(abs, stock)
        largeamount = percentile(absolutevals, cutoff)
        upmoves = sum(stock .> largeamount)
        downmoves = sum(stock .< -largeamount)
        ratios[asset] = downmoves/upmoves
    end
    return ratios
end

function calc_pacf(demeanedreturns::Vector{Float64}, lags)
    pacfcoeffs = pacf(demeanedreturns, lags)
    confinterval = 1.96 / sqrt(length(demeanedreturns))
    return pacfcoeffs, confinterval
end

function calc_pacf(demeanedreturns::Matrix{Float64}, lags)
    bigm = size(demeanedreturns, 1)
    pacfcoeffs = similar(demeanedreturns, bigm, length(lags))
    for assetidx in 1:bigm
        pacfcoeffs[assetidx, :] .= pacf(demeanedreturns[assetidx, :], lags)
    end
    confinterval = 1.96 / sqrt(size(demeanedreturns, 1))
    return pacfcoeffs, confinterval
end

function calc_stylisedfacts(marketval, stocksval, stocksvolume, params;
        maxlag=10, returnspctile=5)
    @unpack perfwindow = params
    lags = collect(1:maxlag)
    # Returns from prices
    marketreturns = Func.assetreturns(marketval, perfwindow)
    stockreturns = Func.assetreturns(stocksval, perfwindow)
    demeanedmarketreturns = Func.demean(marketreturns)
    demeanedstockreturns = Func.demean(stockreturns)
    demeanedmarketreturnssquared = demeanedmarketreturns.^2
    demeanedstockreturnssquared = demeanedstockreturns.^2
    # Autocorrelation of returns
    stockpacf = calc_pacf(demeanedstockreturns, lags)
    marketpacf = calc_pacf(demeanedmarketreturns, lags)
    # Volatility clustering
    stockvolacluster = calc_pacf(demeanedstockreturnssquared, lags)
    marketvolacluster = calc_pacf(demeanedmarketreturnssquared, lags)
    # Kurtosis
    stockkurtoses = Func.calc_kurtoses(demeanedstockreturns)
    # Loss-gain ratio
    lossgainratios = Func.lossgainratio(demeanedstockreturns,
        returnspctile, params)
    # Volume-volatility correlation
    corrs = Func.volavolumecorr(stocksvolume[:, perfwindow[end]+2:end],
        demeanedstockreturns)

    # Z scores of mean being equal to hypotheses
    zScoreKurtoses = HypothesisTests.OneSampleZTest(stockkurtoses, 0.0)
    println(zScoreKurtoses)
    zScoreLossGain = HypothesisTests.OneSampleZTest(lossgainratios, 1.0)
    println(zScoreLossGain)
    zScoreVolumeVola = HypothesisTests.OneSampleZTest(corrs, 0.0)
    println(zScoreVolumeVola)
    flatreturns = collect(Base.Iterators.flatten(demeanedstockreturns))
    samplemean = StatsBase.mean(flatreturns)
    samplestd  = StatsBase.std(flatreturns)
    d = Distributions.Normal(samplemean, samplestd)
    resultKStest = HypothesisTests.ExactOneSampleKSTest(flatreturns,d)
    println(resultKStest)


    return demeanedstockreturns, demeanedmarketreturns, stockpacf, marketpacf,
        stockvolacluster, marketvolacluster, stockkurtoses,lossgainratios,
        corrs, zScoreKurtoses, zScoreLossGain, zScoreVolumeVola, resultKStest
end

function calc_zScore(sample, muZero)
    sampleMean = StatsBase.mean(sample)
    sampleSTD = StatsBase.std(sample)
    zScore = (sampleMean - muZero) / (sampleSTD / sqrt(length(sample)))
    return zScore
end

function plot_stylisedfacts(marketval, stocksval, stylefacts, params)
    @unpack plotpath = params
    demeanedstockreturns = stylefacts[1]
    demeanedmarketreturns = stylefacts[2]
    stockpacf = stylefacts[3]
    marketpacf = stylefacts[4]
    stockvolacluster = stylefacts[5]
    marketvolacluster = stylefacts[6]
    stockkurtoses = stylefacts[7]
    lossgainratios = stylefacts[8]
    corrs = stylefacts[9]
    choice = rand(1:size(stocksval,1)) # Could print random choice on plots
    # NOTE: confidence interval magic number
    confinterval = 1.96 / sqrt(size(demeanedstockreturns, 2))
    # TODO: Restrict the y-axis
    ## Comment out single examples and those not required for results ##
    # plot_stockpacf(demeanedstockreturns[choice, :], confinterval, plotpath)
    # plot_stockvolaclustering(stockvolacluster[1][choice,:], confinterval, plotpath)
    # plot_marketvolaclustering(demeanedmarketreturns.^2, confinterval, plotpath)
    # plot_marketreturnhistogram(demeanedmarketreturns, plotpath)
    # plot_stockreturnhistogram(demeanedstockreturns[choice, :], plotpath)
    # plot_marketpacf(demeanedmarketreturns, confinterval, plotpath)

    plot_pricehistories(stocksval, marketval, plotpath)

    plot_stockpacf(StatsBase.mean(stockpacf[1],dims=1)', confinterval, plotpath)

    plot_stockvolaclustering(StatsBase.mean(stockvolacluster[1],dims=1)',
        confinterval, plotpath)

    plot_stockkurtoses(stockkurtoses, plotpath)
    plot_lossgainratio(lossgainratios, plotpath)
    plot_volavolumecorr(corrs, plotpath)
    plot_logprobs(demeanedstockreturns, plotpath)

    # println(result)
end

# TODO: refactor plots - DRY
function plot_pricehistories(stocksval, marketval, plotpath)
    periods = 1:length(marketval)
    nstocks = 1:100
    StatsPlots.plot(periods, transpose(stocksval)[:, nstocks], legend=:none,
    title="Asset prices history (market in black)")
    plt = StatsPlots.plot!(
        periods, marketval, lw=3, lc=:black, label="Market index", ylims=(0,300))
    png(joinpath(plotpath, "returnshistory"))
    display(plt)
end

function plot_marketreturnhistogram(asset, plotpath)
    StatsPlots.histogram(asset, bins=100,
        title="Market returns vs normal distribution", label="Market returns", legend=:topright)
    x = findmin(asset)[1]:0.001:findmax(asset)[1]
    d = fit(Normal, asset)
    plt = StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
    png(joinpath(plotpath, "marketreturnsvsnormal"))
    display(plt)
end

function plot_stockreturnhistogram(demeanedstockreturns, plotpath)
    asset = demeanedstockreturns
    StatsPlots.histogram(asset, bins=100, title="Example of stock returns",
        label="Random stock's returns", legend=:topright)
    x = findmin(asset)[1]:0.001:findmax(asset)[1]
    d = fit(Normal, asset)
    plt = StatsPlots.plot!(x, pdf.(d, x), lc=:red, label="Normal distn")
    png(joinpath(plotpath, "stockreturnsvsnormal"))
    display(plt)
end


# TODO: X-axis should be Int 1:10, no duplication of legend, tell which stock#
function plot_marketpacf(pacfcoeffs, confinterval, plotpath)
    StatsPlots.bar(pacfcoeffs, title="Market autocorrelation",
        label="Partial correlation coefficients")
    plt = StatsPlots.plot!(0:0.01:10,
        [ones(length(0:0.01:10)) .* confinterval,
        -ones(length(0:0.01:10)) .* confinterval,],
        lc=:red, label="critical values", legend=:bottomright)
    png(joinpath(plotpath, "marketautocorrelation"))
    display(plt)
end

function plot_stockpacf(pacfcoeffs, confinterval, plotpath)
    StatsPlots.bar(pacfcoeffs,
        title="Stock autocorrelation (random stock example)",
        label="Partial autocorrelation coefficients, demeaned returns",
        legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10,
        [ones(length(0:0.01:10)) .* confinterval,
        -ones(length(0:0.01:10)) .* confinterval,],
        lc=:red, label="Critical values")
    png(joinpath(plotpath, "stockautocorrelation"))
    display(plt)
end

function plot_marketvolaclustering(pacfcoeffs, confinterval, plotpath)
    StatsPlots.bar(pacfcoeffs, title="Volatility clustering, market",
        label="Partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10,
        [ones(length(0:0.01:10)) .* confinterval,
        -ones(length(0:0.01:10)) .* confinterval,],
        lc=:red, label="Critical values")
    png(joinpath(plotpath, "marketvolaclustering"))
    display(plt)
end

function plot_stockvolaclustering(pacfcoeffs, confinterval, plotpath)
    StatsPlots.bar(pacfcoeffs, title="Volatility clustering, random stock",
        label="partial autocorrelation coefficients, squared demeaned returns", legend=:bottomright)
    plt = StatsPlots.plot!(0:0.01:10,
        [ones(length(0:0.01:10)) .* confinterval,
        -ones(length(0:0.01:10)) .* confinterval,],
        lc=:red, label="Critical values")
    png(joinpath(plotpath, "stockvolaclustering"))
    display(plt)
end

function plot_stockkurtoses(stockkurtoses, plotpath)
    plt = StatsPlots.bar(sort!(stockkurtoses, rev=true),
        title="Kurtoses of stock returns", label="Kurtoses")
    png(joinpath(plotpath, "kurtoses"))
    display(plt)
end

function plot_lossgainratio(ratios, plotpath)
    StatsPlots.plot(ratios, title="Gain-loss asymmetry in stock returns",
        label="# large losses / # large gains")
    plt = StatsPlots.plot!(0:0.01:200, ones(length(0:0.01:200)) * mean(ratios),
        label="Mean ratio")
    png(joinpath(plotpath, "lossgainratio"))
    display(plt)
end

function plot_volavolumecorr(correlations, plotpath)
    plt = StatsPlots.bar(sort(correlations, rev=true),
        title="Volume-volatility correlation of stocks",
        label="Volume - absolute return correlation coefficients")
    png(joinpath(plotpath, "volavolumecorr"))
    display(plt)
end

function plot_logprobs(returns, plotpath)
    stdevs = StatsBase.std(returns, dims=2)
    returns = returns ./ stdevs
    returns = collect(Base.Iterators.flatten(returns))
    returnshist = PyPlot.hist(returns, 100)
    returnbins = returnshist[1]
    steps = returnshist[2][2:end]
    returnprobs = returnbins ./ sum(returnbins)
    nmldistn = fit(Normal, returns)
    pdfnormal = pdf.(nmldistn, steps)
    normalprobs = pdfnormal ./ sum(pdfnormal)
    pyplot()
    StatsPlots.plot(steps, log10.(returnprobs), label=:"Actual returns",
        xlabel=:"Stocks' daily returns / st dev",
        ylabel=:"log10 of return probability")
    plt = StatsPlots.plot!(steps, log10.(normalprobs),
        label=:"Fitted normal distn")
    png(joinpath(plotpath, "logprobs"))
    display(plt)
end

end  # module Functions
