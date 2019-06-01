using Test, Random, LinearAlgebra

include("../src/types.jl")
include("../src/functions.jl")
include("../src/params.jl")
using .Types, .Func, .Params
#
# # Small parameter values for testing, same list as in /src/params.jl
# const bigk = 3 # Number of funds, 3
# const bign = 4 # Number of investors, 4
# const bigm = 5 # Number of stocks, 5
# const bigt = 6 # Number of time periods, 6
# const mktstartval = 100 # Market index starting value, 100
# const drift = 0.05 # Market index drift, 0.05
# const marketvol = 0.1 # Market volatility (std?), 0.1
# const perfwindow = 1:3 # Performance window range for investors (1:3)
# const betamean = 1 # Average stock beta, 1
# const betastd = 0.3 # Stock beta dispersion, 0.3
# const stockstartval = 100 # Starting value of stocks
# const stockvolrange = range(0.01, stop=0.1, step=0.01)
# # Range of stock volatility (0.1, 0.5). Consider if it should be continuous.
# const invcaprange = (10,1000) # Investors' range of initial capital, (10, 1000)
# const thresholdmean = 0 # Average investor return threshold for her fund
# const thresholdstd = 0.05 # Standard deviation of investor return thresholds
# const portfsizerange = 1:5 # Range of number of stocks in fund portfolio (1,5)
# const impactrange = 0.00001:0.00001:0.0001 # Stock price impact per currency unit

market = Types.MarketIndex(
    zeros(Params.bigt))
stocks = Types.Equity(
    zeros(Params.bigm, Params.bigt),
    zeros(Params.bigm),
    zeros(Params.bigm),
    zeros(Params.bigm))
investors = Types.RetailInvestor(
    zeros(Params.bign, Params.bigk + 1),
    zeros(Params.bign),
    zeros(Params.bign))
funds = Func.Types.EquityFund(
    zeros(Params.bigk, Params.bigm),
    zeros(Params.bigk, Params.bign),
    zeros(Params.bigk, Params.bigt))

Func.marketinit!(market.value)#, mktstartval, perfwindow, drift, marketvol)
println(market.value)

for t in (Params.perfwindow[end]+1):Params.bigt
    Func.marketmove(t, market.value)#, drift, marketvol)
end
println(market.value)

# marketinit!()
# betainit!()
# stockvolinit!()
# stockvalueinit!()
# stockimpactinit!()
# invhorizoninit!()
# invthreshinit!()
# invassetinit!()
# fundcapitalinit!()
# fundstakeinit!()
# fundholdinit!()
# fundvalinit!()