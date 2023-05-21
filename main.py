# import matplotlib
import pandas

# import plotly
import plotly.graph_objects as go


# see https://plotly.com/python/time-series/

def main():
    snp = pandas.read_csv("^GSPC.csv")
    cadusd = pandas.read_csv("CAD=X.csv")
    df = pandas.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv')

    # print(snp)

    snpMax = 0
    snpMin = 100000
    cadusdMax = 0
    cadusdMin = 100000
    dfMax = 0
    dfMin = 100000
    for index, row in snp.iterrows():
        try:
            if row['Close'] > snpMax:
                snpMax = row['Close']
            if row['Close'] < snpMin:
                snpMin = row['Close']
            print(row['Date'], row['Close'], snpMax, snpMin)
        except Exception:
            pass
    for index, row in cadusd.iterrows():
        try:
            if row['Close'] > cadusdMax:
                cadusdMax = row['Close']
            if row['Close'] < cadusdMin:
                cadusdMin = row['Close']
            print(row['Date'], row['Close'], cadusdMax, cadusdMin)
        except Exception:
            pass
    for index, row in df.iterrows():
        try:
            if row['AAPL.Close'] > dfMax:
                dfMax = row['AAPL.Close']
            if row['AAPL.Close'] < dfMin:
                dfMin = row['AAPL.Close']
            print(row['Date'], row['AAPL.Close'], dfMax, dfMin)
        except Exception:
            pass
    fig = go.Figure(data=[go.Candlestick(x=snp['Date'],
                                         open=(snp['Open'] - snpMin) / (snpMax - snpMin),
                                         high=(snp['High'] - snpMin) / (snpMax - snpMin),
                                         low=(snp['Low'] - snpMin) / (snpMax - snpMin),
                                         close=(snp['Close'] - snpMin) / (snpMax - snpMin),
                                         name="SPN"
                                         ),
                          go.Candlestick(x=cadusd['Date'],
                                         open=(cadusd['Open'] - cadusdMin) / (cadusdMax - cadusdMin),
                                         high=(cadusd['High'] - cadusdMin) / (cadusdMax - cadusdMin),
                                         low=(cadusd['Low'] - cadusdMin) / (cadusdMax - cadusdMin),
                                         close=(cadusd['Close'] - cadusdMin) / (cadusdMax - cadusdMin),
                                         name="CADUSD"
                                         ),
                          go.Candlestick(x=snp['Date'],
                                         open=(snp['Open'] - snpMin) / (snpMax - snpMin) * cadusd['Open'],
                                         high=(snp['High'] - snpMin) / (snpMax - snpMin) * cadusd['High'],
                                         low=(snp['Low'] - snpMin) / (snpMax - snpMin) * cadusd['Low'],
                                         close=(snp['Close'] - snpMin) / (snpMax - snpMin) * cadusd['Close'],
                                         name="SPNCAD"
                                         ),
                          go.Candlestick(x=df['Date'],
                                         open=(df['AAPL.Open'] - dfMin) / (dfMax - dfMin),
                                         high=(df['AAPL.High'] - dfMin) / (dfMax - dfMin),
                                         low=(df['AAPL.Low'] - dfMin) / (dfMax - dfMin),
                                         close=(df['AAPL.Close'] - dfMin) / (dfMax - dfMin),
                                         name="AAPL"
                                         )
                          ]
                    )

    # fig = go.Figure(data=[go.Scatter(x=snp['Date'],
    #                                     open=(snp['Open'] - snpMin) / (snpMax - snpMin),
    #                                     high=(snp['High'] - snpMin) / (snpMax - snpMin),
    #                                     low=(snp['Low'] - snpMin) / (snpMax - snpMin),
    #                                     close=(snp['Close'] - snpMin) / (snpMax - snpMin),
    #                                     name="SPN"
    #                                     ),
    #                      go.Scatter(x=cadusd['Date'],
    #                                     open=(cadusd['Open'] - cadusdMin) / (cadusdMax - cadusdMin),
    #                                     high=(cadusd['High'] - cadusdMin) / (cadusdMax - cadusdMin),
    #                                     low=(cadusd['Low'] - cadusdMin) / (cadusdMax - cadusdMin),
    #                                     close=(cadusd['Close'] - cadusdMin) / (cadusdMax - cadusdMin),
    #                                     name="CADUSD"
    #                                     ),
    #                      go.Scatter(x=snp['Date'],
    #                                     open=(snp['Open'] - snpMin) / (snpMax - snpMin) * cadusd['Open'],
    #                                     high=(snp['High'] - snpMin) / (snpMax - snpMin) * cadusd['High'],
    #                                     low=(snp['Low'] - snpMin) / (snpMax - snpMin) * cadusd['Low'],
    #                                     close=(snp['Close'] - snpMin) / (snpMax - snpMin) * cadusd['Close'],
    #                                     name="SPNCAD"
    #                                     ),
    #                      go.Scatter(x=df['Date'],
    #                                     open=(df['AAPL.Open'] - dfMin) / (dfMax - dfMin),
    #                                     high=(df['AAPL.High'] - dfMin) / (dfMax - dfMin),
    #                                     low=(df['AAPL.Low'] - dfMin) / (dfMax - dfMin),
    #                                     close=(df['AAPL.Close'] - dfMin) / (dfMax - dfMin),
    #                                     name="AAPL"
    #                                     )
    #                      ]
    #                )

    fig.show()


if __name__ == "__main__":
    main()
