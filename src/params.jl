module Params
using Random
const big_m = 100 # Number of stocks
const big_n = 500 # Number of investors
const big_t = 1000 # Number of time periods
const big_k = 50 # Number of funds
const marketstartval = 100 # Market index starting value
const drift = 0.05 # Market index drift
const marketvol = 0.1 # Market volatility (std?)
const perfwindow = (1, 3) # Performance window for investor
const betamean = 1 # Average stock beta
const betastd = 0.3 # Stock beta dispersion
const stockvolrange = (0, 0.2)
const capitalrange = (10, 10) # Investors' range of initial capital
rng = MersenneTwister(1999) # Default 'random' number
end #module
