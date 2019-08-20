module Params
using Parameters

default = @with_kw (
    bigm = 250, # Number of stocks
    bign = 1000, # Number of investors
    bigt = 1350, # Number of days
    bigk = 250, # Number of funds

    # Market parameters
    marketstartval = 100,
    drift = 0.000307888,
    marketvol = 0.01,

    # Stock parameters
    impactrange = 0.00075:0.0001:0.0075,
    betamean = 1,
    betastd = 0.3,
    stockstartval = 100,
    stockvolrange = 0.001:0.001:0.01,

    # Investor parameters
    reviewprobability = 1/63,
    perfwindow = 1:100,
    invcaprange = (50,150),
    thresholdmean = 0,
    thresholdstd = 0.05,

    # Fund parameters
    portfsizerange = 10:100,
    plotpath = "./plots/")


end #module
