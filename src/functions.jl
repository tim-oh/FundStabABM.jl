
module Func
using Random

include("params.jl")
import .Params

function marketmove(currentval, rng=MersenneTwister(), drift=Params.drift, marketvol=Params.marketvol)
    nextval = currentval*(1 + drift) + marketvol*randn(rng)
    return nextval
    # INCLUDE IN THIS A MARKET IMPACT FACTOR BASED ON PREVIOUS PERIOD TOTAL MARKET IMPACT
end

function marketinit!(marketvalue, pwc)
    for t in range(2, stop=pwc)
        marketvalue[t] = Func.marketmove(marketvalue[t-1])
    end
    return marketvalue
end

# FIXME: Replace magic randn(100) with Params.?
function stockmove(t::Int, mktvals, currentvals, betas,  stockvolas)
    mktreturn = (mktvals[t] - mktvals[t-1]) / mktvals[t-1]
    nextvals = ((1 .+ mktreturn .* betas) .* currentvals) + stockvolas .* randn(100) .* currentvals
    return nextvals
end

function betainit()
    betas = randn(rng, Params.big_m) .* Params.betastd .+ Params.betamean
    return betas
end

function stockvolinit()
    stockvol = Params.stockvolrange[1] .+ rand(rng, Params.big_m) .*
    (Params.stockvolrange[2] - Params.stockvolrange[1])
    return stockvol
end

function stockvalueinit!(stockvals, pwc, mktvals, betas, stockvolas)
    stockvals[:, 1] = ones(100) * 100
    for t in range(2, stop=pwc)
        stockvals[:, t] = Func.stockmove(
        t, mktvals, stockvals[:, t-1], betas, stockvolas, rng)
        #println()
    end
    return stockvals
end

function invfundinit!(assets, caprange=Params.capitalrange, rng=MersenneTwister())
    # Assign capital of first M investors to one specific fund
    for inv in range(1, stop=Params.big_k)
        assets[inv] = caprange[1] + rand(rng) * (caprange[2] - caprange[1])
    end
    if Params.big_n > Params.big_k
        # For the remaining investors assign capital to a random fund
        for inv in range((Params.big_k+1), stop=Params.big_n)
            pick = round(1 + rand(rng) * Params.big_k)
            pick = convert(Int16, pick)
            assets[pick] = caprange[1] + rand(rng) * (caprange[2] - caprange[1])
        end
    end
    return assets
end

end  # module Functions
