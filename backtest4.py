import matplotlib.pyplot as plt
import numpy as np

import minimize3

saved = []
first = True


def main():
    global first
    partitions = [i for i in np.linspace(0, 1, 11, endpoint=True)]
    for n, i in enumerate(partitions):
        try:
            print(n, partitions[n], partitions[n + 2], partitions[n + 4])
            if first:  # speed boost
                saved.append(minimize3.main(riskToleranceP=4.5,
                                            fakePastP=partitions[n],
                                            fakePresentP=partitions[n + 2],
                                            fakeFutureP=partitions[n + 4]
                                            )
                             )
                first = False
            else:
                saved.append(minimize3.main(riskToleranceP=4.5,
                                            fakePastP=partitions[n],
                                            fakePresentP=partitions[n + 2],
                                            fakeFutureP=partitions[n + 4],
                                            initial_guess=saved[-1][3]
                                            )
                             )

        except IndexError as e:
            pass
            # print(e)
    # print(partitions)
    # minimize3.main(riskToleranceP=4.5, fakePastP=0.3, fakePresentP=0.5, fakeFutureP=0.7)
    data = []
    for i in saved:
        data.append(i[2])
        print(i[2])

    data.sort()

    xs = [x for x, y in enumerate(data)]
    ys = [y for x, y in enumerate(data)]

    print(sum(ys))  # integral of performance

    plt.plot(xs, ys)
    plt.xlabel("Cumulative probability")
    plt.ylabel("Total return, %")
    plt.show()


if __name__ == "__main__":
    main()
