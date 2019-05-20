module Params
using Random

const bigm = 100 # Number of stocks
const bign = 500 # Number of investors
const bigt = 1000 # Number of time periods
const bigk = 50 # Number of funds
const marketstartval = 100 # Market index starting value
const drift = 0.05 # Market index drift
const marketvol = 0.1 # Market volatility (std?)
const perfwindow = 1:3 # Performance window for investor
const betamean = 1 # Average stock beta
const betastd = 0.3 # Stock beta dispersion
const stockstartval = 100 # Price of stocks at the beginning
const stockvolrange = 0:0.2 # Range of stock volatilities
const invcaprange = (10, 10) # Investors' range of initial capital
const horizonrange = 1:100 # Range of investors' evaluation time horizon
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds
const portfsizerange = 10:100 # Range of number of stocks in funds' portfolio
rng = MersenneTwister(1999) # Default 'random' number
end #module

# TODO: turns ranges (a, b) into actual ranges, i.e. a:b
