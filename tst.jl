using Test, Random
include("tmp.jl")
include("src/types.jl")
include("src/params.jl")

using .Types
using .Params
using .Tmp

# TODO: Create keybinding "switch to console, execute last command, switch back"
# TODO: Create keybinding "paste selection in console, execute, switch back"


stockreturns = Func.assetreturns(stocks.value)

nstocks = size(stockreturns,1)
kurtoses = zeros(nstocks)
for stock in 1:size(stockreturns,1)
        kurtoses[stock] = StatsBase.kurtosis(stockreturns[stock,:])
end
