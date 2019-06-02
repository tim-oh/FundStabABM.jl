module Params
using Random

const bigm = 100 # Number of stocks
const bign = 500 # Number of investors
const bigt = 1000 # Number of time periods
const bigk = 50 # Number of funds
const marketstartval = 100 # Market index starting value
const drift = 0.000307888 # Daily market index drift, 8%pa
# TODO: Verify that it should be sqrt(dt), and that this means the below
const marketvol = (1/250)^(1/2) # Market volatility
const marketvol = 0.1 # Market volatility (std?)
const perfwindow = 1:250 # Performance window for investor
const betamean = 1 # Average stock beta
const betastd = 0.3 # Stock beta dispersion
const stockstartval = 100 # Price of stocks at the beginning
const stockvolrange = 0.01:0.01:0.1 # Range of stock volatilitiess
const invcaprange = (50,150) # Investors' range of initial capital
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds
const portfsizerange = 10:100 # Range of number of stocks in funds' portfolio
const impactrange = 0.0001:0.0001:0.001 # Stock price impact per currency unit
end #module
