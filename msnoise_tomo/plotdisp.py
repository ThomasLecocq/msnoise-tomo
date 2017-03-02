import glob

import matplotlib.pyplot as plt
import pandas as pd


def main():
    alldf = []
    for file in glob.glob('TOMO_DISP/*'):
        alldf.append(pd.read_csv(file, index_col=0))

    alldf = pd.concat(alldf)


    alldf["mean"] = alldf.mean()
    alldf["median"] = alldf.median()

    alldf.plot(c='k')
    plt.show()