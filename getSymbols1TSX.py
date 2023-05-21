# https://www.tsx.com/json/company-directory/search/tsx/%5E*

import csv
import pprint

url = "https://www.tsx.com/json/company-directory/search/tsx/%5E*"

import urllib.request, json

with urllib.request.urlopen(url) as u:
    data = json.loads(u.read().decode())
    pprint.pprint(data)

tickers = {}
for i in data['results']:
    for instr in i['instruments']:
        print(instr['symbol'].replace('.', '-') + ".TO")
        tickers[instr['symbol'].replace('.', '-') + ".TO"] = instr['name']

with open("symbolsTSX.csv", 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(tickers)

print(len(tickers))
