module Tmp
using Random, BenchmarkTools, Test, InteractiveUtils

include("src/params.jl")
include("src/functions.jl")
import .Params
import .Func

function flipcoin(randval::Float64)
    if randval<0.5
        return "H"
    else
        return "T"
    end
end
@show code_typed(flipcoin(rand()))

end # module
