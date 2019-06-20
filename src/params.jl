module Params
using Parameters

default = @with_kw (
bigm = 200,
bign = 500,
bigt = 1100,
bigk = 200,

# Market parameters
marketstartval = 100,
drift = 0.000307888,
marketvol = 1,

# Stock parameters
betamean = 1,
betastd = 0.3,
stockstartval = 100,
stockvolrange = 0.005:0.005:0.05,

# Investor parameters
reviewprobability = 1/63,
perfwindow = 1:100,
invcaprange = (50,150),
thresholdmean = 0,
thresholdstd = 0.05,

# Fund parameters
portfsizerange = 10:100,
impactrange = 0.0001:0.0001:0.001,
plotpath = "./plots/")

# Agent and object numbers, length of time
const bigm = 200 # Number of stocks
const bign = 500 # Number of investors
const bigt = 1100 # Number of time periods
const bigk = 200 # Number of funds

# Market parameters
const marketstartval = 100 # Market index starting value
const drift = 0.000307888 # Daily market index drift, 8%pa
const marketvol = 1 # Market volatility

# Stock parameters
const betamean = 1 # Average stock beta
const betastd = 0.3 # Stock beta dispersion
const stockstartval = 100 # Price of stocks at the beginning
const stockvolrange = 0.001:0.001:0.01 # Range of stock volatilitiess

# Investor parameters
const reviewprobability = 1/63 # Daily likelihood of performance review
const perfwindow = 1:100 # Performance window for investor
const invcaprange = (50,150) # Investors' range of initial capital
const thresholdmean = 0 # Average investor return threshold for her fund
const thresholdstd = 0.05 # Standard deviation of investor return thresholds

# Fund parameters
const portfsizerange = 10:100 # Range of number of stocks in funds' portfolio
#const impactrange = 0.0001:0.0001:0.001 # Stock price impact per currency unit

const plotpath = "./plots/"

# Stylised fact-check parameters
# Number of lags?
# Size of confidence interval?


end #module
