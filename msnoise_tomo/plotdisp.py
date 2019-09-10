import glob

import matplotlib.pyplot as plt
import pandas as pd


def main(filterid, comp):
    alldf = []
    for file in glob.glob('TOMO_DISP/%02i/%s/*' % (filterid, comp)):
        print(file)
        alldf.append(pd.read_csv(file, index_col=0))

    alldf = pd.concat(alldf)

    alldf["mean"] = alldf.mean()
    alldf["median"] = alldf.median()

    alldf.plot(c='k', lw="0.5", legend=False)
    plt.savefig("dispersions-f%02i-%s.pdf" % (filterid, comp))
    plt.show()