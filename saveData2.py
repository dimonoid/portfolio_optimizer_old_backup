# import adict as adict
import pandas as pd
import csv
from datetime import date, timedelta
import numpy as np
from yahoofinancials import YahooFinancials
from icecream import *

tickerList = 'CA_Stocks'
# tickerList = 'EmergingMarkets'
# tickerList = 'europe50'
# tickerList = 'tsx60'
# tickerList = 'QQQ'
# tickerList = 'USLargeCap'
# tickerList = 'Wealthsimple'
# tickerList = 'Wealthsimple'
tickerListTime = tickerList + '_10Years'

years = 10
days = 365 * years

# days = 90

repeatDebug = 10000


def main():
    global repeatDebug

    data = []
    with open('symbols' + tickerList + '.csv', newline='') as f:
        reader = csv.reader(f)
        for y in reader:
            for x in y:
                data.append(x)
        # data = list(reader)[0]

    print(data)

    pricesBuffered = []
    for d in data:
        print(d)
        with open('blacklist.csv', 'rt') as f:
            reader = csv.reader(f, delimiter=',')  # good point by @paco
            bl = False
            for row in reader:
                for field in row:
                    if field == d:
                        print("is in blacklist")
                        bl = True
                    if bl: break
                if bl: break
            if bl: continue
        try:

            # set stock ticker symbol
            stock_symbol = d

            # TODO DEBUG
            if repeatDebug > 0:
                repeatDebug -= 1
                print("repeatDebug=", repeatDebug)

            elif stock_symbol not in ['HXS.TO', 'ZQQ.TO']:
                # black list of errors, leveraged, or inverse etfs since it is cheating
                print("Skipping...")
                continue
            # TODO DEBUG

            # stock_symbol = 'DGRC.TO'  # temp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # if stock_symbol == 'CYBR-B.TO':
            #    print('Debug')
            # set date range for historical prices
            end_time = date.today()
            start_time = end_time - timedelta(days=days)

            # reformat date range
            end = end_time.strftime('%Y-%m-%d')
            start = start_time.strftime('%Y-%m-%d')

            # get daily stock prices over date range
            json_prices = YahooFinancials(stock_symbol
                                          ).get_historical_price_data(start, end, 'daily')

            dividends = {}

            try:  # in case there are no dividends
                for dividendDate in json_prices[stock_symbol]['eventsData']['dividends']:
                    # print("dividendDate=" + str(dividendDate))
                    # print("amount=" + str(json_prices[stock_symbol]['eventsData']['dividends'][dividendDate]['amount']))
                    dividends[dividendDate] = json_prices[stock_symbol]['eventsData']['dividends'][dividendDate][
                        'amount']
            except Exception as e:
                print("No dividends")
                print(e.args)

            # print(len(dividends))

            pricesTemp = json_prices[stock_symbol]['prices']
            pricesTempSorted = sorted(pricesTemp, key=lambda k: k['date'])  # should be from old to new

            numberOfShares = 1.0
            # try:
            for n, priceDict in enumerate(pricesTempSorted):
                if len(dividends) == 0:  # speed optimization
                    break

                priceOfOneShare = json_prices[stock_symbol]['prices'][n]['close']
                if np.equal(priceOfOneShare, None) or priceOfOneShare == 0:
                    json_prices[stock_symbol]['prices'][n]['close'] = 0  # will be processed later
                    continue

                # if price is ok, but dividend is None, then assume dividend=0
                if priceDict['formatted_date'] in dividends and \
                        not np.equal(dividends[priceDict['formatted_date']], None):
                    totalPrice = (priceOfOneShare + dividends[priceDict['formatted_date']]) * numberOfShares  # DRIP
                    numberOfShares = totalPrice / priceOfOneShare
                else:
                    totalPrice = priceOfOneShare * numberOfShares

                json_prices[stock_symbol]['prices'][n]['close'] = totalPrice

                # ic(numberOfShares, totalPrice)

                # print(cumulativeDividends, "=")
                # print(json_prices[stock_symbol]['prices'][n]['close'])
                # print('__________________')

            # transform json file to dataframe
            prices = pd.DataFrame(json_prices[stock_symbol]['prices']
                                  )[['formatted_date', 'close']]
            # except Exception as e:
            #    print(e)
            #    exit(0)
            # print(prices)
            # print(prices['close'].isna().sum())
            print(len(prices))
            # if prices['close'].isna().sum() > 0.05 * len(prices):
            # print("Too many NaN values")
            # continue

            if len(prices) < 400:
                print("Too small number of rows")
                # continue

            if np.isnan(prices.loc[prices['formatted_date'] == '2021-07-02'].iloc[0][
                            'close']):  # available for purchase TODO select last trading day
                print('No data for the last day')
                continue

            if json_prices[stock_symbol]['currency'] != 'CAD':  # in canadian dollars
                print("Not CAD, but " + json_prices[stock_symbol]['currency'])
                # continue

            if stock_symbol in {'FQC.TO', 'HQU.TO', 'HSU.TO', 'HFU.TO', 'HGU.TO', 'HMJU.TO', 'HXD.TO', 'HND.TO',
                                'FSF.TO', 'DXP.TO', 'DGRC.TO', 'HMJI.TO', 'HSD.TO', 'HUV.TO', 'HIU.TO'}:
                # black list of errors, leveraged, or inverse etfs since it is cheating
                continue

            prices.rename(columns={'close': stock_symbol}, inplace=True)  # rename in prparation to merging

            pricesBuffered.append(prices)
            # exit(0)  # debug
        except Exception as e:
            print(d + " Error " + str(e.args))
            # print(str(e.args[0]))
            if str(e.args[0]) in ["None of [Index(['formatted_date', 'close'], dtype='object')] are in the [columns]",
                                  "prices", "single positional indexer is out-of-bounds"]:
                with open("blacklist.csv", 'a') as myfile:
                    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                    wr.writerow([d])

    # print(pricesBuffered)

    mergedDf = pd.DataFrame()
    mergedDf['formatted_date'] = ""
    for d in pricesBuffered:
        mergedDf = pd.merge(mergedDf, d, how='outer', on='formatted_date')

    mergedDf = mergedDf.sort_values(by=['formatted_date'], ascending=True)

    mergedDf.to_csv('merged' + tickerListTime + '.csv', index=False)

    print(mergedDf)


if __name__ == "__main__":
    main()
