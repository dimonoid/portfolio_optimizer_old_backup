import csv
import pprint
import pandas as pd

# url = "https://www.tsx.com/json/company-directory/search/tsx/%5E*"

# import urllib.request, json
#
# with urllib.request.urlopen(url) as u:
#    data = json.loads(u.read().decode())
#    pprint.pprint(data)

data = pd.read_csv("WS-Trade-Info-main/etf_info/CA_ETF_WS.csv")

print(data['symbol'].head())

data['symbol'].to_csv('symbolsWealthsimple.csv', index=False, quoting=1, header=False)

# with open("symbolsTSX.csv", 'w') as myfile:
#    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
#    wr.writerow(tickers)
#
# print(len(tickers))
