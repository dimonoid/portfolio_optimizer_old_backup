from datetime import date, timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from yahoofinancials import YahooFinancials

# set stock ticker symbol
import weightFunction

stock_symbol = 'BRAG.TO'

# set date range for historical prices
end_time = date.today()
start_time = end_time - timedelta(days=365 * 10)

# reformat date range
end = end_time.strftime('%Y-%m-%d')
start = start_time.strftime('%Y-%m-%d')

# get daily stock prices over date range
json_prices = YahooFinancials(stock_symbol
                              ).get_historical_price_data(start, end, 'daily')

# transform json file to dataframe
prices = pd.DataFrame(json_prices[stock_symbol]
                      ['prices'])[['formatted_date', 'close']]

print(prices.describe())

# sort dates in descending order
prices.sort_index(ascending=False, inplace=True)

print(prices.head())
print()

# calculate daily logarithmic return
prices['returnsLog'] = (np.log(prices.close / prices.close.shift(-1)
                               )
                        )

prices['returns'] = (prices.close - prices.close.shift(-1)
                     )

# s=1
# for i in np.log(prices.close / prices.close.shift(-1)):
#    if str(i) == "nan":
#        continue
#    s*=(i+1)
#    print(s)

# calculate daily standard deviation of returns
daily_std = np.std(prices.returnsLog)

# annualized daily standard deviation
std = daily_std * 252.75 ** 0.5

# Plot histograms
fig, ax = plt.subplots(1, 1, figsize=(7, 5))
n, bins, patches = ax.hist(
    prices.returns.values,
    bins=100,
    alpha=0.7,
    color='blue'
)

# print(prices.returns.values)
sum = 0
for i in prices.returns.values:
    print(i)
    print(type(i))
    if np.isnan(i):
        continue
    # sum+=i
    sum += weightFunction.weight(i, 6)

print("Weight: " + str(sum))

ax.set_xlabel('log return of stock price')
ax.set_ylabel('frequency of log return')
ax.set_title('Historical Volatility for ' +
             stock_symbol)

# get x and y coordinate limits
x_corr = ax.get_xlim()
y_corr = ax.get_ylim()

# make room for text
header = y_corr[1] / 5
y_corr = (y_corr[0], y_corr[1] + header)
ax.set_ylim(y_corr[0], y_corr[1])

# print historical volatility on plot
x = x_corr[0] + (x_corr[1] - x_corr[0]) / 30
y = y_corr[1] - (y_corr[1] - y_corr[0]) / 15
ax.text(x, y, 'Annualized Volatility: ' + str(np.round(std * 100, 1)) + '%',
        fontsize=11, fontweight='bold')
x = x_corr[0] + (x_corr[1] - x_corr[0]) / 15
y -= (y_corr[1] - y_corr[0]) / 20

# save histogram plot of historical price volatility
fig.tight_layout()
fig.savefig('historical volatility.png')
fig.show()
