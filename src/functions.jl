
module Func
using Random


#include("params.jl")
include("types.jl")
#import .Params
import .Types

# NOTE Used to have default parameters: drift=Params.drift, marketvol=Params.marketvol. Consider these for the functions in general, so that function arguments are minimal outside of testing.
function marketmove(currentval, drift, marketvol)
    nextval = currentval*(1 + drift) + marketvol * randn()
    return nextval
    # INCLUDE IN THIS A MARKET IMPACT FACTOR BASED ON PREVIOUS PERIOD TOTAL MARKET IMPACT
end

function marketinit!(marketval, marketstartval, horizon, drift, marketvol)
    marketval[1] = marketstartval
    for t in 2:horizon
        marketval[t] = Func.marketmove(marketval[t-1], drift, marketvol)
    end
    return marketval
end

function stockmove(t::Int, mktvals, currentvals, betas, stockvolas)
    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    nextvals = ((1 .+ mktreturn .* betas) .* currentvals) + stockvolas .* randn(length(stockvolas)) .* currentvals
    return nextvals
end

# QUESTION: Should I provide default values from Params? This one used to be:
# betainit!(betas, bigm=Params.bigm, betastd=Params.betastd,
#    betamean=Params.betamean)
function betainit!(betas, bigm, betastd, betamean)
    betas .= randn(bigm) .* betastd .+ betamean
    return betas
end

# Function for a discrete volatility range
function stockvolinit!(stockvol, volrange, bigm)
    stockvol = rand(volrange, bigm)
    return stockvol
end

# Function for a continuous volatility range
#function stockvolinit(bigm=Params.bigm, volrange=Params.stockvolrange)
#    stockvol = volrange[1] .+ rand(bigm) .* (volrange[end] - volrange[1])
#    return stockvol
#end

function stockvalueinit!(
    stocks, stockstartval, horizon, marketval)
    stocks.value[:,1] .= stockstartval
    for t in 2:horizon
        stocks.value[:, t] = Func.stockmove(
        t, marketval, stocks.value[:, t-1], stocks.beta, stocks.vol)
    end
    return stocks.value
end

function stockimpactinit!(impact, impactrange, horizon)
    impact .= rand(impactrange, length(impact))
    return impact
end

function invhorizoninit!(emptyvals, hrange)
    horizons = rand(hrange, length(emptyvals))
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

function fundcapitalinit!(fundvals, investorassets)
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
    return holdings
end

function fundvalinit!(fundvals, holdings, stockvals, horizon)
    bigk = size(fundvals)[1]
    for t in 2:horizon # t=1 is initialised by fundcapitalinit!
        for k in 1:bigk
            fundvals[k, t] = sum(holdings[k, :] .* stockvals[:, t])
        end
    end
    return fundvals
end

function drawreviewers(bign)
    reviewers = rand(bign) .< 1/63 # On average, a fund reviews once a quarter
    return findall(reviewers)
end

function perfreview(t, reviewers, investors, fundvals)
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

function liquidate!(holdings, stakes, divestments)
    sellorders = Array{Float64}(undef, 0, size(holdings, 2) + 1)
    for row in 1:size(divestments, 1)
        investor = divestments[row, 1]
        fund = divestments[row, 2]

        # Note the minus
        sellorder = hcat((holdings[fund,:] .* -stakes[fund,investor])',investor)
        sellorders = vcat(sellorders, sellorder)

        remainder = (1 - stakes[fund, investor])
        #@assert cashout + sum(-sellorder[1:end-1]) == holdings[fund, :]
        # TODO: fix this assert
        holdings[fund, :] = holdings[fund, :] .* remainder

        stakes[fund, investor] = 0
        if sum(stakes[fund,: ]) > 0
            stakes[fund, :] = stakes[fund, :] ./ sum(stakes[fund,: ])
        end
    end
    return sellorders, holdings, stakes
end

# QUESTION: Is marketmake! a good candidate for multiple dispatch, where it uses
# one method for sales and one for purchases? Would I have to create an
# abstract 'Order' type that can be a sales or purchase order?a

# TODO: Add description to each function, according to some standard format
# that should probably inclucde input, output and purpose, plus comments on
# assumptions or non-obvious points.

# TODO: Market making for sell orders is now impact * number of shares,
# but it should be impact * value of shares.
# That means I need to recalculate the test! Bugger.


function marketmake!(stockvals, impacts, t, orders::Types.SellMarketOrder)
    cashout = Array{Float64}(undef, 0, 2)
    netimpact = sum(orders.values, dims=1)' .* impacts
    stockvals[:, t] = (1 .+ netimpact) .* stockvals[:, t]
    for order in 1:size(orders.values, 1)
        investor = orders.investors[order]
        amount = -sum(orders.values[order, :] .* stockvals[:, t])
        cashout = vcat(cashout, [investor amount])
    end
    return stockvals, cashout
end

function disburse!(invassets, divestments, cashout)
    for row in 1:size(divestments,1)
        inv = divestments[row, 1]
        fund = divestments[row, 2]
        @assert(inv == cashout[row, 1]) # Ensure divestment matches cashout
        invassets[inv, fund] = 0
        invassets[inv, end] = cashout[row, 2]
    end
    return invassets
end

# TODO: function marketmake!(stocks.value, orders::Purchases)
#     BODY PLEASE
# end

end  # module Functions
