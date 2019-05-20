
module Func
using Random

include("params.jl")
import .Params

function marketmove(currentval, drift=Params.drift, marketvol=Params.marketvol)
    nextval = currentval*(1 + drift) + marketvol * randn()
    return nextval
    # INCLUDE IN THIS A MARKET IMPACT FACTOR BASED ON PREVIOUS PERIOD TOTAL MARKET IMPACT
end

function marketinit!(marketval, marketstartval, horizon)
    marketval[1] = marketstartval
    for t in 2:horizon
        marketval[t] = Func.marketmove(marketval[t-1])
    end
    return marketval
end

function stockmove(t::Int, mktvals, currentvals, betas, stockvolas)
    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    nextvals = ((1 .+ mktreturn .* betas) .* currentvals) + stockvolas .* randn(length(stockvolas)) .* currentvals
    return nextvals
end

function betainit(bigm=Params.bigm, betastd=Params.betastd, betamean=Params.betamean)
    betas = randn(bigm) .* betastd .+ betamean
    return betas
end

function stockvolinit(bigm=Params.bigm, volrange=Params.stockvolrange)
    stockvol = volrange[1] .+ rand(bigm) .* (volrange[2] - volrange[1])
    return stockvol
end

function stockvalueinit!(
    stocks, stockstartval, horizon, marketval)
    stocks.value[:,1] .= stockstartval
    for t in 2:horizon
        stocks.value[:, t] = Func.stockmove(
        t, marketval, stocks.value[:, t-1], stocks.beta, stocks.vol)
    end
    return stocks.value
end

function invhorizoninit!(emptyvals, hrange)
    horizons = rand(hrange[1]:hrange[2], length(emptyvals))
    return horizons
end

function invthreshinit!(emptyvals, mean, std)
    thresholds = randn(4) * std .+ mean
    return thresholds
end

function invassetinit!(invassets, caprange, bigk)
    bign = size(invassets)[1]
    capital = rand(caprange[1]:caprange[2], bign)
    for inv in 1:bigk
        invassets[inv, inv] = capital[inv]
    end
    for inv in (bigk+1):bign
        invassets[inv, rand(1:bigk)] = capital[inv]
    end
    return invassets
end

function fundvalinit!(fundvals, investorassets)
    fundvals[:,1] = sum(investorassets, dims=1)[1:end-1]
    return fundvals
end

function fundstakeinit!(fundstakes, investorassets)
    ncols = size(investorassets)[1] - 1
    for col in 1:ncols
        fundstakes[col, :] = investorassets[:, col] ./ sum(investorassets[:, col])
    end
    return fundstakes
end

function fundholdinit!(holdings, portfsizerange, capital, stockvals)
    bigk = size(holdings)[1]
    bigm = size(holdings)[2]
    nstocks = rand(portfsizerange, bigk) # Number of stocks in each portfolio
    for k in 1:bigk # Loop over funds
        selection = rand(1:bigm, nstocks[k]) # Randomly select w/ replacement
        for stock in selection # Loop over
            holdings[k, stock] +=
            1 / length(selection) * (capital[k] / stockvals[stock])
        end
    end
    println("capital", capital)
    println("holdings", holdings)
    return holdings
end

end  # module Functions
