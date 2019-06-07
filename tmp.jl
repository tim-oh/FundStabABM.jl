module Tmp
using Random

include("src/params.jl")
include("src/functions.jl")
import .Params
import .Func

# Calculate returns
marketreturns = Func.returns(market.value)
stockreturns = Func.returns(stocks.value)



end # module
