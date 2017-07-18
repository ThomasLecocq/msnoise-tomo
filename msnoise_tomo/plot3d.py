from msnoise.api import *


def main():
    db = connect()
    PER= get_config(db, "ftan_periods", plugin="Tomo")
    periods = np.array([float(pi) for pi in PER.split(',')])

    all = []
    for n, per in enumerate(periods):
        f = "tomo_%.1fs.txt"%per
        print(f)
        if not os.path.isfile(f):
            continue
        tmp = np.loadtxt(f)
        for i in tmp.shape[0]:
            for j in tmp.shape[1]:
                all.append([ i, j, n, tmp[i,j]] )
    print(all)






if __name__ == '__main__':
    main()