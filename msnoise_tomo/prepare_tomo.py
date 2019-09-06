from msnoise.api import *


def main():
    db = connect()
    for filter in get_filters(db):
        filterid = filter.ref
        for comp in get_components_to_compute(db):
            alldf = []
            for station1, station2 in get_station_pairs(db, used=True):
                sta1 = "%s_%s" % (station1.net, station1.sta)
                sta2 = "%s_%s" % (station2.net, station2.sta)
                pair = "%s_%s_MEAN.csv" % (sta1, sta2)
                fn = os.path.join("TOMO_DISP", "%02i" % filterid, comp, pair)
                if not os.path.isfile(fn):
                    continue
                print("Reading", fn)                    
                tmp = pd.read_csv(fn, index_col=0, delimiter=',')
                tmp.columns = ["%s_%s"%(sta1, sta2)]
                alldf.append(tmp)
                del tmp
            alldf = pd.concat(alldf, axis=1)
            print(alldf.head())
            alldf = alldf.interpolate()
            print(alldf.head())
            if not os.path.isdir("TOMO_FILES"):
                os.makedirs("TOMO_FILES")
        
            PER= get_config(db, "ftan_periods", plugin="Tomo")
            PER = np.array([float(pi) for pi in PER.split(',')])
            
            for per in PER:
                
                tmp = alldf.copy().iloc[np.abs(alldf.index-per).argsort()[:1]].T
                print(tmp)
                print(tmp.columns)
                tmp.columns = [per,]
                
                # tmp = alldf.loc[per].to_frame()
                tmp["sta1"] = [".".join(t.split('_')[0:2]) for t in tmp.index]
                tmp["sta2"] = [".".join(t.split('_')[2:4]) for t in tmp.index]
                print(tmp.head())
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
                print("*"*100)
                print(tmp.head())
                # ,per, -> nan -> all get dropped.
                tmp = tmp.loc[:,["sta1","sta2","per",per, "error","dist"]]
                print("*" * 100)
                print(tmp.head())
                print("*" * 100)
                tmp = tmp.dropna()
                print("*" * 100)
                of = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "TestGroupVel_%.1fsGLISN.dat"%per)
                if not os.path.isdir(os.path.split(of)[0]):
                    os.makedirs(os.path.split(of)[0])
                tmp.to_csv(of, index=False, header=False, sep=" ")
            df = []
            for s in get_stations(db):
                df.append(["%s.%s"%(s.net, s.sta), s.net, s.Y, s.X, 0])
            df = pd.DataFrame(df)
            of = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "GLISN_STACoord.dat")
            df.to_csv(of, index=False, header=False, sep=" ")
        
            of = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "GLISNGrid.dat")
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
