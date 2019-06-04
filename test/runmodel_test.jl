

@test all(stocks.value .>= 0)
@test all(stocks.vol .>= 0)
@test all(stocks.impact .>= 0)

@test all(funds.holdings .>= 0)
@test all(funds.stakes .>= 0)
@test all(funds.value .>= 0)

@test all(investors.assets .>= 0)
@test all(investors.horizon .>= 0)
@test all(1 .>= investors.threshold .>= -1)

@test all(buyorder.values .>= 0)
@test all(buyorder.funds .>= 0)
@test all(sellorder.values .<= 0)
@test all(sellorder.investors .>= 0)
