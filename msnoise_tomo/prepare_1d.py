import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    #
    files = glob.glob("*_disp.txt")
    all = []
    Periods = np.arange(0, 10+1.1, 0.1)
    for file in files:
        D = np.loadtxt(file)
        isort = np.argsort(D[:,0])
        D = D[isort]
        per = D[:,0]
        disper = D[:,1]
        dci=np.interp(Periods, per, disper, left=np.nan, right=np.nan)
        all.append(dci)
    all = np.array(all).T
    all = pd.DataFrame(all)
    all.to_csv("all.csv")
    Xmean = all.copy().apply(func=np.nanmean, axis=1)
    Xstd = all.copy().apply(func=np.nanstd, axis=1)
    print Xstd.head()
    count = all.copy().count(axis=1)
    all.plot()
    plt.show()
    tmp = Xmean.to_frame()
    tmp.index = Periods
    tmp["std"] = Xstd
    #tmp["count"] = count
    print tmp.head()
    tmp.to_csv("mean.csv")


if __name__ == "__main__":
    main()