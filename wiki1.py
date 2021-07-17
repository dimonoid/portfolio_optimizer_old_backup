import urllib.request, json
import pprint
import csv

url = "https://sheet2api.com/wiki-api/List_of_Canadian_exchange-traded_funds/en/ETF+Table"


def main():
    with urllib.request.urlopen(url) as urlOpened:
        data = json.loads(urlOpened.read().decode())
    pprint.pprint(data)

    symbols = ["CLU-C.NE", "QIF.NE"]
    for s in data:
        aaa = s["Symbol"]
        aaa = aaa[5:]
        if "." not in aaa:
            aaa = aaa + ".TO"  # automatically assume it is TSX
        else:
            aaa = aaa.replace(".", "-") + ".TO"
        symbols.append(aaa)
    pprint.pprint(symbols)
    print(len(symbols))

    with open("symbols.csv", 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(symbols)


if __name__ == "__main__":
    main()
