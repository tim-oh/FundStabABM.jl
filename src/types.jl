module Types
using Random

include("params.jl")
import .Params

abstract type Banker end

abstract type Asset end

abstract type AssetManager end

abstract type Investor end

struct Equity <: Asset
    value::Array{Float64}
    beta::Array{Float64}
    vol::Array{Float64}
end

struct MarketIndex <: Asset
    value::Array{Float64}
end

struct RetailInvestor <: Investor
    assets::Array{Float64} # N possible stakes in M funds and (1) cash
    horizon::Array{Int64} # Positive integer, draw from range, can use rand(S)
    threshold::Array{Float64} # ~N(0, 5)
end

struct EquityFund <: AssetManager
    holdings::Array{Float64}
    stakes::Array{Float64}
    value::Array{Float64}
end

struct MarketMaker <: Banker
    orderbook::Array{Float64}
end

end # module
