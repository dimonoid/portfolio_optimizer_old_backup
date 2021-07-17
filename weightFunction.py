import math
import sys

import matplotlib.pyplot as plt
import numpy as np


# 2.15<=risk
def weight(x, risk):
    r = 10 ** (risk - 5)
    if x >= 0:
        y = x
        return y
    else:
        # print("_")
        # print(x)
        # print(r)
        v = -(x / r)
        # v = min(v, 100+0.00001*v)
        # v = max(v, 1/100-0.00001*v)
        # print(v)
        y = (- math.exp(v) + 1) * r
        return y


def main():
    x = np.linspace(-1, 1, 101)
    y = [weight(i, 5) for i in x]

    plt.scatter(x, y)
    plt.show()


if __name__ == "__main__":
    main()
