from msnoise.api import *


def main():
    db = connect()
    alldf = []
    for station1, station2 in get_station_pairs(db, used=True):
        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)
        pair = "%s_%s_MEAN.csv" % (sta1, sta2)
        fn = os.path.join("TOMO_DISP", pair)
        if not os.path.isfile(fn):
            continue
        tmp = pd.read_csv(fn, index_col=0,delimiter=',')

        alldf.append(tmp)
        del tmp
    alldf = pd.concat(alldf, axis=1)

    if not os.path.isdir("TOMO_FILES"):
        os.makedirs("TOMO_FILES")

    PER= get_config(db, "ftan_periods", plugin="Tomo")
    PER = np.array([float(pi) for pi in PER.split(',')])
    for per in PER:
        tmp = alldf.loc[per].to_frame()
        tmp["sta1"] = [t.split('_')[0] for t in tmp.index]
        tmp["sta2"] = [t.split('_')[1] for t in tmp.index]
        tmp["per"] = per
        tmp["error"] = 0
        dists = []
        for s1, s2 in zip(tmp["sta1"], tmp["sta2"]):
            n1,s1 = s1.split(".")
            s1 = get_station(db, n1, s1)
            n2,s2 = s2.split(".")
            s2 = get_station(db, n2, s2)
            d = get_interstation_distance(s1, s2)*1000
            dists.append(d)
        tmp["dist"] = dists
        print tmp.head()
        tmp = tmp[["sta1","sta2","per",per,"error","dist"]]
        print tmp.head()
        tmp = tmp.dropna()
        of = os.path.join("TOMO_FILES", "TestGroupVel_%.1fsGLISN.dat"%per)
        tmp.to_csv(of, index=False, header=False, sep=" ")
    df = []
    for s in get_stations(db):
        df.append(["%s.%s"%(s.net, s.sta), s.net, s.Y, s.X, 0])
    df = pd.DataFrame(df)
    of = os.path.join("TOMO_FILES", "GLISN_STACoord.dat")
    df.to_csv(of, index=False, header=False, sep=" ")

    of  = os.path.join("TOMO_FILES", "GLISNGrid.dat")
    f = open(of,'w')
    xstep = float(get_config(db, "xstep", plugin="Tomo"))
    ystep = float(get_config(db, "ystep", plugin="Tomo"))
    minlon = df[3].min() - xstep*2 - 0.005
    maxlon = df[3].max() + xstep*2 + 0.005
    minlat = df[2].min() - ystep*2 - 0.005
    maxlat = df[2].max() + ystep*2 + 0.005
    f.write("%f %f\n"%(minlon, maxlon))
    f.write("%f %f\n"%(minlat, maxlat))
    f.write("%f %f\n"%(xstep, ystep))
    f.close()


if __name__ == "__main__":
    main()
