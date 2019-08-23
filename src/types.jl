module Types
using Random

include("params.jl")
import .Params

abstract type Asset end

abstract type AssetManager end

abstract type Investor end

abstract type Order end

struct Equity <: Asset
    value::Matrix{Float64} # Value history of each stock
    beta::Vector{Float64} # Stocks' betas
    vol::Vector{Float64} # Stocks' volatilities
    impact::Vector{Float64} # Symmetric price impact for 1 money unit bought/sold
    volume::Matrix{Float64} # Money value of the stock traded in a day
end

# NOTE: Mismatch between market value and what the index should be worth
# since the former determines the latter, while it should be vice versa.
# Implicit assumption is that there are other stocks not traded by our agents.
struct MarketIndex <: Asset
    value::Vector{Float64} # Index, used to drive stock values
end

struct RetailInvestor <: Investor
    assets::Matrix{Float64} # N possible stakes in M funds and (1) cash
    horizon::Vector{Int64} # Positive integer, draw from range, can use rand(S)
    threshold::Vector{Float64} # ~N(0, 5)
end

struct EquityFund <: AssetManager
    holdings::Matrix{Float64} # Units of each stock held in the fund's portfolio
    stakes::Matrix{Float64} # Share of assets each investor in the funds owns
    value::Matrix{Float64} # Fund's value history
end

mutable struct BuyMarketOrder <: Order
    values::Matrix{Float64} # Money amount spent on each stock for one of .funds
    funds::Vector{Int64} # Fund that a set of stocks (.values) gets bought for
end

mutable struct SellMarketOrder <: Order
    values::Matrix{Float64} # Money amount spent on each stock for one of .funds
    investors::Vector{Int64} # Investor that stocks (.values) are liquidated for
end

end # module
