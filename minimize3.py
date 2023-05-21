import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as optimize

import weightFunction

# dataRaw = pd.read_csv('mergedTest.csv')
# dataRaw = pd.read_csv('mergedTest_1Year.csv')
# dataRaw = pd.read_csv('mergedTest_1Year.csv')

# dataRaw = pd.read_csv('mergedWealthsimple_1Year.csv')
# dataRaw = pd.read_csv('mergedWealthsimple_3Years.csv')
# dataRaw = pd.read_csv('mergedUSLargeCap_10Years.csv')
# dataRaw = pd.read_csv('mergedQQQ_10Years.csv')
# dataRaw = pd.read_csv('mergedtsx60_10Years.csv')
# dataRaw = pd.read_csv('mergedeurope50_10Years.csv')
# dataRaw = pd.read_csv('mergedEmergingMarkets_10Years.csv')
dataRaw = pd.read_csv('mergedCA_Stocks_10Years.csv')
# dataRaw = pd.read_csv('mergedWealthsimple_10Years.csv')
# dataRaw = pd.read_csv('mergedTest2_10YearsLite.csv')

# pd.set_option('max_columns', None)
# pd.set_option('max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.float_format', lambda x: '%.3f' % x)

# print(data)
dataRaw = dataRaw.sort_values(by=['formatted_date'], ascending=True)
dataRaw = dataRaw.reset_index(drop=True)
# print(data)

dataRaw = dataRaw.replace(0, np.nan)
dataRaw = dataRaw.fillna(method='ffill')
# dataRaw= dataRaw.replace(0, np.nan)
# dataRaw = dataRaw.fillna(method='bfill')  # TODO may cause problems with backtesting

# dataRaw.set_index('formatted_date')
# dates=data[]

# dataRaw.iloc[:,1:]
data = dataRaw.div(dataRaw.fillna(method='bfill').loc[0][
                       dataRaw.columns != 'formatted_date']).replace(np.nan,
                                                                     0)  # divide by first non np.nan or non 0 value

data.iloc[:, :] = np.log(data.iloc[:, :]).replace(-np.inf,
                                                  np.nan)  # TODO testing if it improves prediction accuracy

data['formatted_date'] = dataRaw['formatted_date']
data = data.sort_values(by=['formatted_date'], ascending=True)  # does nothing, I hope

# data = data.fillna(0)

print(data)

# global parameters
plotEvery = 0
plotEveryConst = 1000
recommendation = pd.DataFrame()
weightHistory = [0.] * 10
performanceHistory = [0.] * 10
performanceHistoryBefore = [0.] * 10
performanceHistoryAfter = [0.] * 10
override = False

# moneyToInvest = 34500
moneyToInvest = 5000
riskTolerance = 3.5

tickers = data.columns[:-1]
length = len(tickers)

dataDiff = data[tickers].diff(periods=1, axis=0).fillna(0)

valid = dataRaw.loc[:, tickers].notna()

fakePast = 0.3
fakePresent = 0.5
fakeFuture = 0.7

performanceTotal = -1
performanceBefore = -1
performanceAfter = -1

tradingDaysPerYearRatio = 1 / 252.75

pricesOnly = dataDiff.loc[:, tickers]


def visualize(modifiedCoeff, portfolioSeries, portfolioPunishedSeries, portfolioDiff, weight, HXS, ZQQ):
    global plotEvery
    global plotEveryConst
    global recommendation
    global weightHistory

    global performanceHistory
    global performanceHistoryBefore
    global performanceHistoryAfter

    global performanceTotal
    global performanceBefore
    global performanceAfter

    global fakePast
    global fakePresent
    global fakeFuture

    # s = portfolio[0]
    # portfolio -= 1
    # portfolio *= 10
    # portfolio += 1
    # portfolioPunishedSeries += portfolio[0]
    # HXS /= HXS[0]

    # HXS *= (portfolio[0] / HXS[0])
    # ZQQ *= (portfolio[0] / ZQQ[0])

    # portfolioPunishedSeries /= portfolioPunishedSeries[0]

    if plotEvery == 0 or override:  # display 1st time, then every loop, and then final solution at the end

        # data.iloc[0].

        # recommendation = pd.DataFrame(
        #    {'tickers': tickers.to_list(),
        #     'percentage': (params * 100).tolist(),
        #     'performance': data.iloc[-1][2:].sub(data.iloc[0][2:]).values
        #     })
        # recommendation = recommendation.sort_values(by=['percentage'], ascending=False)
        # print(recommendation.head(1000))
        table(modifiedCoeff)
        # calculate sum of all percentages, should be 80-100

        # s = 0
        # for i in range(len(recommendation)):
        #    s += recommendation.iloc[i][1]
        #    if s >= 0.8 * sum(params):
        #        print("stocks for 80% of recommendation=" + str(i + 1) + "/" + str(len(recommendation)), end=" ")
        #        break

        # print("sum=" + str(s))
        # print(portfolio.iloc[0])

        timeLength = len(portfolioSeries) - 1

        performanceTotal = (
                                   math.exp(
                                       portfolioSeries.iloc[round(fakeFuture * timeLength)] -
                                       portfolioSeries.iloc[round(fakePast * timeLength)]
                                   ) - 1
                           ) * 100 / (
                                   (fakeFuture - fakePast) * timeLength * tradingDaysPerYearRatio
                           )
        print(str(performanceTotal) + "%T")  # Total performance
        performanceBefore = (
                                    math.exp(
                                        portfolioSeries.iloc[round(fakePresent * timeLength)] -
                                        portfolioSeries.iloc[round(fakePast * timeLength)]
                                    ) - 1
                            ) * 100 / (
                                    (fakePresent - fakePast) * timeLength * tradingDaysPerYearRatio
                            )
        print(str(performanceBefore) + "%P")  # Past performance
        performanceAfter = (
                                   math.exp(
                                       portfolioSeries.iloc[round(fakeFuture * timeLength)] -
                                       portfolioSeries.iloc[round(fakePresent * timeLength)]
                                   ) - 1
                           ) * 100 / (
                                   (fakeFuture - fakePresent) * timeLength * tradingDaysPerYearRatio
                           )
        print(str(performanceAfter) + "%F")  # Future performance

        print(str(100 * performanceAfter / performanceBefore) + "% test/train ratio")

        # performanceTotal = ((portfolio.iloc[-1] / portfolio.iloc[0]) - 1) * 100
        # print(str(performanceTotal) + "%")  # performance over period

        # performanceBefore = ((portfolio.iloc[round(timeLength / 2)] / portfolio.iloc[0]) - 1) * 100
        # print(str(performanceBefore) + "%")  # 1st half
        # performanceAfter = ((portfolio.iloc[-1] / portfolio.iloc[round(timeLength / 2)]) - 1) * 100
        # print(str(performanceAfter) + "%")  # 2nd half

        # if plotEvery == 0:
        fig, axes = plt.subplots(2, 3, figsize=(20, 10))
        axes[0, 0].hist(
            portfolioDiff,
            bins=100,
            alpha=0.7,
            color='blue'
        )
        axes[0, 0].title.set_text('derivatives')

        axes[0, 1].plot(np.linspace(0, 1, len(portfolioSeries)), portfolioSeries, color='g', )
        axes[0, 1].plot(np.linspace(0, 1, len(portfolioPunishedSeries)), portfolioPunishedSeries, color='orange')
        axes[0, 1].plot(np.linspace(0, 1, len(HXS)), HXS, color='k')
        axes[0, 1].plot(np.linspace(0, 1, len(ZQQ)), ZQQ, color='r')

        axes[0, 1].axvline(x=fakePast, color='r')  # past limit
        axes[0, 1].axvline(x=fakePresent, color='k')  # fake present
        axes[0, 1].axvline(x=fakeFuture, color='g')  # future limit

        # axes[0, 1].plot(np.linspace(0, 1, len(ZQQ)), ZQQ*3, color='r')
        # axes[0, 1].plot(np.linspace(0, 1, len(portfolioPunishedDf)), portfolioPunishedDf)
        axes[0, 1].title.set_text('graph')
        axes[0, 1].autoscale()

        weightHistory.append(weight)
        weightHistory.pop(0)
        axes[0, 2].plot(np.linspace(0, 1, len(weightHistory)), weightHistory)
        axes[0, 2].title.set_text('weightHistory')

        # get x and y coordinate limits #######################################
        x_corr = axes[1, 1].get_xlim()
        y_corr = axes[1, 1].get_ylim()

        # make room for text
        header = y_corr[1] / 5
        y_corr = (y_corr[0], y_corr[1] + header)
        # axes[1, 1].set_ylim(y_corr[0], y_corr[1])

        # print historical volatility on plot
        x = x_corr[0] + (x_corr[1] - x_corr[0]) / 30
        y = y_corr[1] - (y_corr[1] - y_corr[0]) / 15
        #######################################################################
        performanceHistory.append(performanceTotal)
        performanceHistory.pop(0)
        axes[1, 1].plot(np.linspace(0, 1, len(performanceHistory)), performanceHistory)
        axes[1, 1].title.set_text('performanceHistoryTotal')
        # axes[1, 1].text(x, y, 'Performance total: ' + str(performanceTotal),
        #                fontsize=11, fontweight='bold')

        performanceHistoryBefore.append(performanceBefore)
        performanceHistoryBefore.pop(0)
        axes[1, 0].plot(np.linspace(0, 1, len(performanceHistoryBefore)), performanceHistoryBefore)
        axes[1, 0].title.set_text('performanceHistoryTrain')
        # axes[1, 0].text(x, y, 'Performance before: ' + str(performanceBefore), fontsize=11, fontweight='bold')

        performanceHistoryAfter.append(performanceAfter)
        performanceHistoryAfter.pop(0)
        axes[1, 2].plot(np.linspace(0, 1, len(performanceHistoryAfter)), performanceHistoryAfter)
        axes[1, 2].title.set_text('performanceHistoryTest')
        # axes[1, 2].text(x, y, 'Performance after: ' + str(performanceAfter), fontsize=11, fontweight='bold')

        fig.tight_layout()
        fig.show()

        plotEvery = plotEveryConst
    else:
        plotEvery -= 1


def table(modifiedCoeff):
    global recommendation
    global moneyToInvest
    global data

    # print(tickers.to_list())
    # print(data.iloc[-1][2:].to_list()) #price now
    # print(data.iloc[0][2:].to_list()) #first entry price
    # print(((data.iloc[-1][2:].divide(data.iloc[0][2:]) - 1) * 100).values.tolist())

    recommendation = pd.DataFrame(
        {'Ticker': tickers.to_list(),
         'Recommendation%': (modifiedCoeff.iloc[-1] * 100).tolist(),
         # 'Last year performance%': ((data.iloc[-1][:-1].divide(data.iloc[-250][:-1]) - 1) * 100).values.tolist(),
         'Last year performance%': (np.exp((data.iloc[-12][:-1] - (data.iloc[-250][:-1])).values.tolist()) - 1) * 100,
         # 'Price now': (dataRaw.iloc[-1][1:]).to_list(),
         # 'Number of shares to buy': [round(i) for i in (params * moneyToInvest / data.iloc[-1][:-1]).to_list()]
         }
    )

    recommendation = recommendation.sort_values(by=['Recommendation%'], ascending=False)
    # print(recommendation)
    recommendation['Recommendation%'] = \
        recommendation['Recommendation%'] * 100 / recommendation['Recommendation%'].sum(axis=0)
    print(recommendation)  # adjusted to 100%


def f(params):
    global plotEvery
    global weightHistory
    global performanceHistory
    global fakePast
    global fakePresent
    global fakeFuture
    global pricesOnly

    # print(params.tolist())

    params100Sum = sum(params)  # copied or referenced ?
    params = params / params100Sum  # to beat prices approaching 0, makes super low values useless #TODO just in case
    sumOfValidCoefficients = (valid * params).sum(axis=1)  # coefficient of parameters
    modifiedCoeff = pd.DataFrame(1 / sumOfValidCoefficients).dot(pd.DataFrame(params).T) * valid.values
    # params = modifiedCoeff
    # params = params * modifiedCoeff  # returns vector

    # pricesOnly = data.loc[:, tickers]
    # portfolio = (pricesOnly * params).sum(axis=1)  # calculate linear combination

    portfolioDiff = (pricesOnly * modifiedCoeff.values).sum(axis=1)

    # dataDiff

    # a = (pricesOnly.astype(bool).astype(int) * params)
    # b = a.divide(a.sum(axis=1), axis=0)

    # portfolio /= portfolio[0]

    # diff = (portfolio - portfolio.shift(1))
    # diffTest = (portfolio / portfolio.shift(1)) - 1
    # print(diff)
    # print(diffTest)
    # diff = diffTest  # TODO temp

    integralPunished = 0.0
    portfolioPunished = []
    # portfolioPunished.append(integralPunished)
    integral = 0.0
    portfolio = []
    # portfolio.append(integral)
    for i, d in enumerate(portfolioDiff):
        # print(i)
        # print(type(i))
        if np.isnan(d):
            continue

        if len(portfolioDiff) * fakePast <= i < len(portfolioDiff) * fakePresent:  # for back testing performance
            try:
                # integralPunished += d
                integralPunished += weightFunction.weight(d, riskTolerance)  # / sum(params)
                portfolioPunished.append(integralPunished)

                # integral += d
            except OverflowError:
                # integral += d  # TODO workaround, change later
                print("Overflow!!!")
        else:
            portfolioPunished.append(integralPunished)  # same value from last loop
            # integral += d

        integral += d
        portfolio.append(integral)

    # integral = 1
    # for i, d in enumerate(diff):
    #    # print(i)
    #    # print(type(i))
    #    if np.isnan(d):
    #        continue
    #
    #    if i <= len(diff) / 2:  # for back testing performance
    #        # try:
    #        integral *= weightFunction.weight(d, 4) + 1  # / sum(params)
    #        #integral *= (d + 1)
    #    portfolioPunished.append(integral)
    #        # print("d="+str(d))
    #        # print("integral="+str(integral))
    portfolioPunishedSeries = pd.Series(portfolioPunished)
    portfolioSeries = pd.Series(portfolio)
    # except OverflowError:
    #    integral += d #TODO workaround, change later

    # try:
    #    integral += weightFunction.weight(d, 6)  # / sum(params)
    # except OverflowError:
    #    integral += d
    # break
    # portfolio /= portfolio[0]
    ###################################
    # punishment = 0
    # rollingMax = portfolio[0]
    # for price in portfolio:
    #    if price < rollingMax:
    #        punishment += rollingMax - price
    #        # print(punishment)
    #    else:
    #        rollingMax = price
    ###################################

    # print("-integral=" + str(-integral))

    # print("sum=" + str(sum(params)))

    # print("test 100%=" + str((math.log((sum(params) - 1) ** 2 + 0.01) + 2) * 10000))
    # return -integral + (sum(params) - 1) ** 2 + punishment
    weight = -(integralPunished)  # + ((params100Sum - 1) ** 2) * 0.0001  # + punishment*0.001
    # print((-(integralPunished), ((params100Sum - 1) ** 2) * 0.0001), 'sum =', params100Sum)  # , punishment * 0.001))

    # weight = -integral + 1 / (sum(params100)+0.01)
    # return -integral + (math.log((sum(params) - 1) ** 2 + 0.01)+2)*10000
    # return -integral + punishment
    # weight = -integral

    # print("weight=" + str(weight))

    visualize(modifiedCoeff, portfolioSeries, portfolioPunishedSeries, portfolioDiff, weight,
              data['SHOP.TO'],
              # data['ZQQ.TO'])
              data['SHOP.TO'])

    return weight


def main(riskToleranceP, fakePastP, fakePresentP, fakeFutureP, initial_guess=np.array([1 / length] * length)):
    global override
    global riskTolerance
    global fakePast
    global fakePresent
    global fakeFuture

    riskTolerance = riskToleranceP  # override
    fakePast = fakePastP
    fakePresent = fakePresentP
    fakeFuture = fakeFutureP

    # initial_guess = np.array([0.96309659, 0.00242862, 0.52890531])

    # method='Powell', seems to work fine
    # L-BFGS-B

    # TNC
    # COBYLA slow
    # SLSQP
    # trust-constr

    # initial_guess = np.array(
    #    [7.76414045e-08, 1.00948821e-07, 5.62505260e-08, 1.13052379e-01,
    #     8.86762488e-08, 1.28133761e-07, 1.80109622e-08, 4.19877808e-08,
    #     4.37002995e-08, 2.34750763e-08, 5.98927898e-08, 1.03225570e-07,
    #     1.52960324e-07, 1.31610812e-08, 1.96278744e-08, 6.90804635e-09,
    #     3.33197141e-08, 4.05176795e-09, 1.96790610e-08, 5.74077348e-08,
    #     8.29009683e-09, 2.26744232e-07, 5.03068227e-08, 2.38633177e-07,
    #     8.90583854e-09, 1.15601554e-07, 5.36049772e-05, 2.52527894e-07,
    #     1.19025119e-09, 1.17640958e-08, 4.59398770e-08, 1.49768861e-05,
    #     1.43187628e-07, 1.91853206e-08, 5.08541606e-08, 7.61357121e-08,
    #     1.30894377e-07, 7.98925555e-08, 4.98812179e-08, 2.26784901e-07,
    #     4.33247081e-08, 5.01668686e-08, 3.64641269e-08, 2.30836708e-07,
    #     6.58426228e-05, 1.53799749e-08, 3.54729947e-08, 8.81403294e-08,
    #     2.80894471e-02, 4.11136424e-02, 2.57597926e-08, 1.08356884e-03,
    #     3.12869441e-08, 3.70112368e-08, 9.52222076e-08, 2.00861228e-08,
    #     3.62748684e-08, 5.74435318e-08, 2.99429075e-08, 7.74238707e-08,
    #     6.74669576e-08, 2.60755134e-04, 3.30907677e-08, 4.36335254e-03,
    #     1.55627127e-07, 7.71252599e-08, 9.92610118e-08, 9.69969145e-09,
    #     4.30585069e-02, 1.16891591e-07, 7.83095268e-08, 1.11607869e-07,
    #     1.86006344e-08, 5.23677235e-08, 1.54107888e-07, 2.72275829e-08,
    #     6.27633712e-05, 2.61778994e-03, 5.22897530e-08, 2.86106984e-08,
    #     1.10350211e-03, 1.25438848e-08, 1.97946189e-08, 1.76593207e-08,
    #     2.44728169e-08, 4.10779001e-09, 1.56155529e-08, 5.02768898e-08,
    #     3.49595768e-08, 8.39325461e-08, 5.52719783e-08, 9.30709708e-09,
    #     4.00631297e-08, 5.02411660e-07, 1.19476050e-07, 5.54895510e-08,
    #     4.39126706e-08, 8.83993316e-09, 8.02025680e-08, 7.27086078e-08,
    #     1.61810029e-07, 2.21100009e-08, 7.30631414e-09, 1.24702614e-02,
    #     1.15636316e-07, 4.48453449e-08, 4.28567724e-08, 7.27947342e-09,
    #     7.58444573e-08, 5.93438716e-09, 1.41740979e-07, 3.55536400e-07,
    #     7.65310265e-08, 6.32422360e-08, 1.93435706e-07, 1.51046282e-06,
    #     2.90923303e-08, 2.76466498e-04, 7.65915513e-09, 3.44272503e-08,
    #     4.06183757e-03, 1.43276271e-07, 7.54790768e-09, 6.11095954e-08,
    #     8.29918553e-08, 1.17152488e-07, 1.97344448e-03, 2.21711879e-08,
    #     8.67568930e-08, 3.70814020e-08, 1.67902207e-08, 6.93569859e-08,
    #     7.75757092e-08, 2.09500410e-07, 1.00918241e-07, 3.11376122e-08,
    #     3.75978494e-08, 1.06053996e-07, 4.73815352e-08, 1.91526524e-04,
    #     1.45498343e-07, 5.47065783e-08, 2.31699363e-08, 2.44084851e-08,
    #     6.83301255e-08, 4.74306364e-08, 1.06070330e-07, 8.84739475e-08,
    #     1.05493921e-07, 2.56623464e-07, 9.61735302e-08, 2.30086733e-07,
    #     4.00996591e-08, 8.34012601e-08, 1.14183544e-07, 1.05902015e-07,
    #     9.17049624e-06, 1.32004218e-07, 4.51461034e-07, 2.21075151e-08,
    #     1.92290274e-08, 1.49469340e-07, 3.46141709e-08, 1.47697916e-07,
    #     1.64823957e-07, 2.00777919e-07, 1.45134238e-07, 3.30900895e-07,
    #     2.61115988e-01, 4.60156863e-07, 5.05984842e-08, 2.68139831e-07,
    #     6.98627883e-08, 4.52623356e-08, 2.91639934e-07, 5.00496437e-08,
    #     5.14320470e-08, 7.60665027e-08, 2.07622161e-08, 3.57905142e-09,
    #     3.32722578e-08, 2.01554024e-07, 5.12342680e-08, 8.43573555e-08,
    #     1.21769263e-08, 2.69733254e-02, 9.87214444e-03, 3.94320425e-04,
    #     2.53841527e-07, 5.90791186e-05, 4.75513129e-08, 5.19953439e-07,
    #     5.59163859e-08, 9.89062501e-08, 3.08367359e-07, 1.16569945e-07,
    #     8.47688498e-08, 1.00724977e-07, 2.65302512e-08, 1.54266753e-04,
    #     1.89658355e-07, 1.27244829e-07, 5.92256454e-08, 1.94288688e-07,
    #     8.17584996e-08, 2.17514279e-07, 1.18164454e-01, 7.25950630e-08,
    #     5.58465937e-06, 1.56247865e-07, 3.40272085e-08, 2.24969472e-08,
    #     5.40486226e-08, 6.92881764e-08, 5.60098633e-08, 4.17107839e-03,
    #     2.75305567e-07, 2.77561085e-07, 5.89036341e-07, 1.00127391e-07,
    #     2.35686024e-08, 2.98057148e-06, 4.13933397e-07, 1.53571944e-07,
    #     8.93965029e-08, 5.44894161e-08, 3.71148634e-08, 2.54190596e-08,
    #     5.84699257e-08, 2.01145125e-08, 1.15250340e-07, 2.76689178e-05,
    #     5.23057215e-08, 1.12843424e-08, 5.82115448e-09, 1.80027967e-08,
    #     1.66331710e-08, 1.54199921e-08, 4.67819574e-08, 4.70165466e-08,
    #     1.30717395e-07, 1.06789423e-07, 8.29744461e-08, 1.06832991e-07,
    #     3.84121524e-08, 7.42540014e-09, 8.27694147e-08, 4.64858680e-08,
    #     1.17496906e-08, 1.84464470e-08, 1.97954706e-08, 8.93394184e-09,
    #     1.03383712e-08, 2.93884309e-03, 8.58230700e-09, 7.27386223e-08,
    #     6.92407939e-08, 1.46090437e-07, 9.23526345e-08, 1.13741714e-08,
    #     1.32832559e-07, 4.08772757e-08, 2.21819386e-07, 4.93478069e-08,
    #     1.41650015e-07, 1.43065693e-07, 3.51263174e-08, 2.34306682e-08,
    #     4.40445479e-08, 5.06615178e-08, 1.00754137e-07, 2.02619916e-08,
    #     6.94684427e-07, 9.86773310e-08, 6.23595064e-09, 4.81082206e-08,
    #     5.00336280e-08, 1.26207850e-08, 1.15108244e-07, 6.72843902e-08,
    #     1.50147659e-07, 2.82126478e-08, 1.18464832e-07, 2.88355413e-08,
    #     2.02393596e-08, 1.11809319e-07, 9.15594324e-08, 1.22781434e-07,
    #     8.34353100e-08, 4.61988538e-08, 4.05760999e-08, 8.69064215e-08,
    #     6.04826936e-09, 3.38407913e-08, 1.00522311e-08, 8.98532717e-09,
    #     4.61024747e-09, 8.36961073e-08, 7.32223511e-09, 1.30132792e-07,
    #     9.42254345e-08, 1.01782699e-08, 1.17758280e-07, 3.68942388e-07,
    #     7.57831765e-08, 1.91765531e-08, 6.32529730e-08, 3.01152404e-08,
    #     3.30360063e-07, 4.64075555e-08, 3.84637271e-08, 1.12784218e-07,
    #     2.19612048e-08, 1.10714670e-07, 8.80217470e-08, 2.60937700e-08,
    #     8.24751952e-08, 2.83062206e-08, 3.94703823e-08, 1.70962716e-07,
    #     3.20575242e-08, 5.57363901e-08, 5.02604720e-08, 5.42359856e-03,
    #     2.69248654e-08, 8.55448857e-08, 3.90419334e-08, 3.12200147e-08,
    #     2.70824634e-01, 7.25248755e-08, 7.59540970e-08, 1.84977072e-08,
    #     1.21676232e-07, 1.28238901e-07, 1.50988831e-07, 2.95790826e-08,
    #     2.78920125e-08, 7.76560388e-08, 9.20499011e-07, 3.72729119e-08,
    #     2.43663299e-08, 8.47435264e-08, 1.31827000e-07, 2.17395224e-08,
    #     8.79858006e-08, 9.61029947e-09, 1.12205741e-07, 6.96652060e-08,
    #     1.12661039e-07, 1.80563629e-08, 2.16892288e-08, 9.44250769e-08,
    #     1.85774114e-08, 7.37293186e-04, 4.63875920e-07, 5.89678750e-07,
    #     6.56326107e-08, 1.12039690e-07, 1.05781210e-07, 1.26062202e-03,
    #     1.27444366e-07, 2.86278287e-07, 1.14789509e-07, 1.30566498e-04,
    #     5.37400432e-08, 8.35925835e-08, 4.34987355e-08, 1.65067459e-07,
    #     1.28099163e-07, 9.23366128e-09, 2.25539491e-07, 1.76244101e-08,
    #     2.13456662e-08, 4.43106086e-09, 4.59527988e-08, 6.40848658e-08,
    #     7.66378083e-08, 1.32032369e-07, 5.63303385e-08, 1.34258935e-07,
    #     3.88725270e-05, 2.91723352e-08, 6.73387903e-08, 8.79350773e-08,
    #     4.18409025e-08, 9.58562359e-08, 1.08356290e-08, 5.37729973e-08,
    #     1.21420769e-07, 1.24123761e-08, 5.04086687e-08, 3.48407045e-08,
    #     9.80253787e-08, 1.86767318e-08, 1.08528879e-08, 1.35429753e-08,
    #     3.06015992e-08, 3.18211377e-08, 4.31345448e-09, 7.05777148e-09,
    #     2.49892863e-07, 8.20002034e-09, 8.19866895e-08, 1.73816452e-07,
    #     1.05196090e-07, 4.17829516e-08, 1.51567755e-08, 3.37903881e-07,
    #     6.28606129e-07, 3.73830916e-08, 4.71806721e-07, 1.26405099e-07,
    #     2.14376243e-08, 4.49112227e-08, 2.14877985e-08, 2.07015726e-08,
    #     4.27623401e-08, 1.02288406e-07, 4.02967702e-08, 4.64508319e-08,
    #     6.93052018e-08, 2.20686928e-08, 3.94208002e-08, 1.18238654e-07,
    #     4.56732209e-07, 1.14826929e-07, 1.05599250e-07, 1.65384621e-07,
    #     3.12362889e-08, 6.28562095e-08, 7.14491434e-08, 6.50424693e-08,
    #     8.21954715e-08, 7.28254168e-08, 1.42894661e-08, 8.33263115e-08,
    #     3.51152849e-08, 4.32216038e-08, 6.62808301e-07, 4.15677463e-08,
    #     1.18873287e-07, 1.65213015e-08, 1.21357280e-08, 4.33255197e-08,
    #     4.56027293e-08, 2.67651436e-08, 2.19608532e-07, 6.76627754e-07,
    #     9.37128615e-09, 8.77191054e-08, 1.76247693e-08, 5.03689139e-07,
    #     2.71566951e-07, 1.22894400e-07, 5.20962552e-07, 8.61454810e-08,
    #     3.65785947e-07, 1.05783468e-07, 8.94311703e-08, 5.13579789e-08,
    #     2.62946267e-08, 2.75686664e-08, 9.45016197e-08, 3.10886825e-05,
    #     5.25558125e-09, 4.55761476e-08, 5.68551701e-08, 4.82625977e-08,
    #     1.70978943e-07, 1.42599448e-07, 4.67629134e-07, 1.51429949e-07,
    #     3.78012616e-07, 3.88716428e-08, 1.13699736e-08, 8.83938315e-08,
    #     5.95617952e-08, 4.70183332e-08, 2.58072619e-07, 1.08181680e-07,
    #     9.34334979e-08, 1.22256104e-07, 9.04621898e-08, 7.39981247e-08,
    #     8.18112647e-04, 3.11516657e-07, 1.12443671e-07, 9.60706759e-08,
    #     4.87162149e-08, 1.19503291e-07, 5.39021644e-08, 1.17311226e-07,
    #     3.12422793e-07, 3.74178354e-08, 1.61714140e-08, 1.84874984e-08,
    #     7.47753438e-08, 6.35568551e-08, 6.75830051e-09, 4.08978681e-07,
    #     8.23794589e-08, 2.21550520e-07, 6.69681443e-08, 5.73285961e-09,
    #     1.03463460e-07, 1.15741131e-07, 7.98756585e-08, 4.12516554e-08,
    #     8.02104321e-08, 9.95070263e-08, 6.28878402e-08, 2.88241220e-08,
    #     1.10697055e-08, 2.06654683e-07, 4.39216047e-08, 3.91092572e-07,
    #     2.49464212e-07, 1.32748336e-08, 3.47546590e-08, 6.09939769e-08,
    #     1.30912682e-06, 6.27519539e-09, 2.82172780e-08, 4.72712969e-08,
    #     1.18192530e-07, 6.44227393e-08, 2.84880750e-07, 2.66313146e-08,
    #     6.50980965e-08, 5.26627823e-08, 1.04133454e-07, 1.20412942e-07,
    #     1.08362186e-07, 1.13714939e-05, 5.64577617e-07, 1.16010303e-07,
    #     1.02124518e-06, 5.14724052e-02, 7.49269996e-08, 1.28007418e-07,
    #     3.16399419e-08, 2.21079423e-08, 7.19147545e-09, 2.77564372e-08,
    #     1.01419331e-03, 1.75244512e-07, 9.12322758e-08, 1.12339325e-04,
    #     7.87169573e-07, 9.40099681e-08, 1.32190789e-07, 1.37245920e-07,
    #     1.09618487e-07, 8.72488012e-09, 3.98672461e-08, 5.47578728e-07,
    #     1.51294150e-07, 1.74158818e-07, 3.14957315e-08, 7.70781768e-08,
    #     5.51109565e-08, 1.11456402e-07, 4.94641597e-08, 4.76747165e-08,
    #     1.44557029e-07, 1.84034076e-04, 1.67599634e-08, 1.94804474e-05,
    #     3.50777535e-07, 1.09152787e-07, 1.60460145e-08, 8.31509021e-09,
    #     4.03739381e-08, 1.13609467e-07, 4.13824604e-08, 5.63343159e-08,
    #     1.64318745e-08, 1.55634541e-08, 4.89631168e-08, 1.74376525e-07,
    #     4.63399586e-08, 9.26477876e-08]
    # )

    def con(t):
        # print(t)
        return sum(t) - 1

    cons = {'type': 'eq', 'fun': con}

    result = optimize.minimize(f, initial_guess, method='SLSQP', bounds=[(0, 10) for i in range(length)],
                               constraints=cons)

    print(tickers.to_list())
    fitted_params = repr(result.x)
    print(fitted_params)

    if not result.success:
        print(result.message)

    print("Done! Printing best result...")
    override = True
    print(riskTolerance)
    f(result.x)
    override = False

    return performanceTotal, performanceBefore, performanceAfter, result.x


if __name__ == "__main__":
    main(riskToleranceP=4.5, fakePastP=0.7, fakePresentP=0.99, fakeFutureP=0.991)

# Start of Covid -> 0.855

# override = True
# f(np.array([5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 4.97518203e-18, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96092711e-05, 5.96086099e-05,
#       5.96088167e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 4.38531753e-02, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.71484086e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.67301085e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96087934e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96093112e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96090473e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.87512501e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 6.76801490e-01, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.90821170e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086143e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       9.57499951e-03, 8.35258287e-01, 5.96086099e-05, 5.96086099e-05,
#       3.25300059e-01, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 4.82852836e-02, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 3.57840205e-05,
#       4.70353618e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.95695070e-05, 5.96086099e-05, 5.95481636e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       3.00401290e-01, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05, 5.96086099e-05,
#       5.96086099e-05, 5.96086099e-05, 5.96086099e-05]))
