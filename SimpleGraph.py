import plotly.express as px
import plotly.graph_objects as go


def main():
    df = px.data.stocks()
    print(df.columns)
    AVG = []
    for index, row in df.iterrows():
        coeff = [1 / 6] * 6
        value = 0
        for n in range(len(row) - 1):
            value += row[n + 1] * coeff[n]
        print(value)
        AVG.append(value)

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=df["date"],
                             y=AVG,
                             mode='lines',
                             name='AVG',
                             text='AVG',
                             # line=dict(
                             #    color='firebrick',
                             #    width=4)
                             )
                  )
    for i in df.columns[1:]:
        fig.add_trace(go.Scatter(x=df["date"],
                                 y=df[i],
                                 mode='lines',
                                 name=i,
                                 text=i,
                                 )
                      )

    fig.show()
    print(weight(df['AAPL']))


def multiply(columns, coeff):
    return None


def weight(series):
    growth = series[-1] - series[0]

    # volatility =


if __name__ == "__main__":
    main()
