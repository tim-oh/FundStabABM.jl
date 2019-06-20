#module Tmp
using Random, BenchmarkTools, Test, InteractiveUtils, JLD, PyPlot
using Base.Iterators, StatsPlots, StatsBase
#
# include("src/params.jl")
# include("src/functions.jl")
# import .Params
# import .Func

data = load("examplereturns.jld")
stocks = data["stockreturns"][:,100:end]
stdevs = StatsBase.std(stockreturns, dims=2)
stocks = stocks ./ stdevs

stocks = abs.(stocks)
stocks = collect(flatten(stocks))
stockshist = PyPlot.hist(stocks, 100)
returns = stockshist[1]
returnsprob = returns ./ sum(returns)
bins = stockshist[2][2:end]
logbins = log10.(bins)
logreturns = log10.(returnsprob)
StatsPlots.plot(logbins, logreturns, label=)
#end # module
