module Func
using Random

include("types.jl")
include("params.jl")
import .Types, .Params

# NOTE Used to have default parameters: drift=Params.drift, marketvol=Params.marketvol. Consider these for the functions in general, so that function arguments are minimal outside of testing.

# TODO: Add description to each function, according to some standard format
# that should probably inclucde input, output and purpose, plus comments on
# assumptions or non-obvious points.

function marketmove(
    t,
    marketvalue=Params.market.value,
    drift=Params.drift,
    marketvol=Params.marketvol)

    marketvalue[t] = marketvalue[t-1]*(1 + drift)
                   + marketvol * randn()

    return marketvalue
end

function marketinit!(
    marketval,
    marketstartval=Params.marketstartval,
    perfwindow=Params.perfwindow,
    drift=Params.drift,
    marketvol=Params.marketvol)

    horizon = perfwindow[end]
    marketval[1] = marketstartval
    for t in 2:horizon+1
        marketval = Func.marketmove(t, marketval, drift, marketvol)
    end

    return marketval
end

function stockmove(
    t::Int,
    mktvals,
    currentvals,
    betas,
    stockvolas)

    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    nextvals = ((1 .+ mktreturn .* betas) .* currentvals) + stockvolas .* randn(length(stockvolas)) .* currentvals

    return nextvals
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

    stockvol = rand(volrange, bigm)

    return stockvol
end

function stockvalueinit!(
    stocks,
    marketval,
    stockstartval=Params.stockstartval,
    perfwindow=Params.perfwindow)

    horizon=perfwindow[end]
    stocks.value[:,1] .= stockstartval
    for t in 2:(horizon + 1)
        stocks.value[:, t] = Func.stockmove(
        t, marketval, stocks.value[:, t-1], stocks.beta, stocks.vol)
    end

    return stocks.value
end

function stockimpactinit!(
    impacts,
    impactrange=Params.impactrange)

    impacts .= rand(impactrange, length(impacts))

    return impacts
end

function invhorizoninit!(
    emptyvals,
    hrange=Params.perfwindow)

    horizons = rand(hrange, length(emptyvals))

    return horizons
end

function invthreshinit!(
    emptyvals,
    mean=Params.thresholdmean,
    std=Params.thresholdstd)

    thresholds = randn(4) * std .+ mean

    return thresholds
end

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

    fundvals[:,1] = sum(investorassets, dims=1)[1:end-1]

    return fundvals
end

function fundstakeinit!(
    fundstakes,
    investorassets)

    ncols = size(investorassets, 1) - 1
    for col in 1:ncols
        fundstakes[col,:] = investorassets[:,col] ./ sum(investorassets[:,col])
    end

    return fundstakes
end

function fundholdinit!(
    holdings,
    capital,
    stockvals,
    portfsizerange=Params.portfsizerange)

    kfunds = size(holdings, 1)
    mstocks = size(holdings, 2)
    nofm = rand(portfsizerange, kfunds) # Number of stocks in each portfolio
    for k in 1:kfunds # Loop over funds
        selection = rand(1:mstocks, nofm[k]) # Randomly select w/ replacement
        for stock in selection # Loop over stocks
            holdings[k, stock] +=
            (capital[k] / length(selection)) / stockvals[stock]
        end
    end

    return holdings
end

function fundvalinit!(
    fundvals,
    holdings,
    stockvals,
    perfwindow=Params.perfwindow)

    horizon=perfwindow[end]
    kfunds = size(fundvals, 1)
    for t in 2:(horizon + 1) # t=1 is initialised by fundcapitalinit!
        for k in 1:kfunds
            fundvals[k, t] = sum(holdings[k, :] .* stockvals[:, t])
        end
    end

    return fundvals
end

function drawreviewers(
    bign=Params.bign)

    reviewers = rand(bign) .< 1/63 # On average, a fund reviews once a quarter

    return findall(reviewers)
end

# QUESTION: Is it sufficient to compare fund values or do we need to compare the
# value of investors' stakes over time?
function perfreview(
    t,
    reviewers,
    investors,
    fundvals)

    divestments = Array{Int64}(undef, 0, 2)
    for rev in reviewers
        horizon = investors.horizon[rev]
        fnds = findall(investors.assets[rev, :] .> 0)
        for fnd in fnds
            valchange =
            (fundvals[fnd,t] - fundvals[fnd, t-horizon]) / fundvals[fnd,horizon]
            if valchange < investors.threshold[rev]
                divestments = vcat(divestments, [rev fnd])
            end
        end
    end

    return divestments
end

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
        salevals = (holdings[fund, :] .* -stakes[fund, investor] .* stockvals)'
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

function marketmake!(
    stocks,
    t,
    ordervals)

    netimpact = sum(ordervals, dims=1)' .* stocks.impact
    stocks.value[:, t] .= vec((1 .+ netimpact) .* stocks.value[:, t])

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
function disburse!(
    funds::Types.EquityFund,
    sharesout,
    stockvals)

    for row in 1:size(sharesout, 1)
        fund = convert(Int64, sharesout[row, 1])
        funds.holdings[fund, :] += sharesout[row, 2:end]
        funds.value[fund, :] .= vec(funds.holdings[fund, :]' * stockvals)
    end

    return funds
end

# Note that targets can be both range or list
function fundrevalue!(
    funds,
    targets,
    stockvals)

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
        # Anchor investor notes her investment in the fund
        investors.assets[inv, egg] = capital
        investors.assets[inv, end] = 0
        # Re-use fund holdings init fctn to draw new stocks and generate order
        buyorder = (fundholdinit!(funds.holdings[egg, :]', capital,
        stockvals[:, t], portfsizerange)' .* stockvals[:, t])'
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
    t)

    kfunds = size(investors.assets,2) - 1

    buyorders = Types.BuyMarketOrder(
    Array{Float64}(undef, 0, size(funds.holdings, 2)), Array{Float64}(undef, 0))

    reinvestors = findall(vec(investors.assets[:, end] .> 0))
    for reinv in reinvestors
        injection = investors.assets[reinv, end]
        # Choose best-performing fund
        choice = bestperformer(funds.value, investors.horizon[reinv], t)
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
    t,
    perfwindow=Params.perfwindow)

    horizon=perfwindow[end]
    # Return between investor's horizon and now/t
    horizonreturns =
    (fundvals[:, t] - fundvals[:, t-horizon]) ./ fundvals[:, t-horizon]
    # Index of best-performing fund
    _, bestfund = findmax(horizonreturns)

    return bestfund
end


end  # module Functions
