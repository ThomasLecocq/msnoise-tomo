#from obspy.imaging.cm import viridis
import matplotlib.pyplot as plt
from .ftan_call import pickgroupdispcurv
from matplotlib.cm import hot_r as viridis
from msnoise.api import *


def main(pair, bmin, bmax, show):
    PLOTDIAGR=show
    PLOTRAWDISP=0
    PLOTDISPALL=0
    SAVEFILES=0
    if not os.path.isdir("TOMO_DISP"):
        os.makedirs("TOMO_DISP")
    db = connect()
    PER= get_config(db, "ftan_periods", plugin="Tomo")
    PER = np.array([float(pi) for pi in PER.split(',')])
    fmin = float(get_config(db, "ftan_fmin", plugin="Tomo"))
    fmax = float(get_config(db, "ftan_fmax", plugin="Tomo"))
    vgmin = float(get_config(db, "ftan_vgmin", plugin="Tomo"))
    vgmax = float(get_config(db, "ftan_vgmax", plugin="Tomo"))
    if bmin is None:
        bmin = float(get_config(db, "ftan_bmin", plugin="Tomo"))
    if bmax is None:
        bmax = float(get_config(db, "ftan_bmax", plugin="Tomo"))
    diagramtype = get_config(db, "ftan_diagramtype", plugin="Tomo")

    nfreq = int(get_config(db, "ftan_nfreq", plugin="Tomo"))
    ampmin = float(get_config(db, "ftan_ampmin", plugin="Tomo"))

    db = connect()

    while is_next_job(db, jobtype='TOMO_FTAN') or pair:
        SACfilelist = []
        if not pair:
            jobs = get_next_job(db, jobtype='TOMO_FTAN')

            for job in jobs:
                netsta1, netsta2 = job.pair.split(':')
                print(netsta1, netsta2)
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta1.replace('.','_'), netsta2.replace('.','_')))
                # SACfilelist.append(fn)
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta2.replace('.','_'), netsta1.replace('.','_')))
                # SACfilelist.append(fn)
                fn = os.path.join("TOMO_SAC", "%s_%s_MEAN.SAC"%(netsta1.replace('.','_'), netsta2.replace('.','_')))
                SACfilelist.append(fn)
        else:
            for pi in pair:
                netsta1, netsta2 = pi.split('_')
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta1.replace('.','_'), netsta2.replace('.','_')))
                # SACfilelist.append(fn)
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta2.replace('.','_'), netsta1.replace('.','_')))
                # SACfilelist.append(fn)
                fn = os.path.join("TOMO_SAC", "%s_%s_MEAN.SAC"%(netsta1.replace('.','_'), netsta2.replace('.','_')))
                SACfilelist.append(fn)
            pair = None
        print(SACfilelist)
        GVdisp = [{},]*len(SACfilelist)
        maxp = np.empty(len(SACfilelist))
        for i, filename in enumerate(SACfilelist):
            NET1, STA1, NET2, STA2, crap = os.path.split(filename)[1].split('_')

            st = read(filename)
            dist = st[0].stats.sac.dist
            dt = st[0].stats.delta

            per, disper = pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                                            diagramtype, nfreq, ampmin, dist)
            basename = "%s.%s_%s.%s_%s" % (NET1, STA1, NET2, STA2, crap)
            basename = basename.replace(".SAC","")
            for _ in ["write_amp.txt",
                      "write_disp.txt",
                      "write_FP.txt",
                      "write_ph.txt",
                      "write_TV.txt",
                      ]:
                os.system("rename %s %s"% (_, _.replace("write",basename)))
            if PLOTDIAGR:
                U = np.loadtxt('%s_TV.txt'%basename)
                P = np.loadtxt('%s_FP.txt'%basename)
                xmin = min(P)
                xmax = max(P)
                ymin = min(U)
                ymax = max(U)
                amp = np.loadtxt('%s_amp.txt'%basename).T
                iu = np.where( (disper>=ymin) & (disper<=ymax) )
                per = per[iu]
                disper = disper[iu]
                Per, Vitg = np.meshgrid(P,U)
                plt.figure()
                plt.contourf(Per, Vitg, amp, 35, cmap=viridis)
                plt.colorbar()
                #plt.contour(Per, Vitg, amp, 35, colors='k')

                plt.plot(per, disper,'-ok',lw=1.5)
                plt.xlim(xmin, xmax)
                plt.ylim(ymin, ymax)

                # TODO add axes labels depending on diagramtype
                plt.xlabel("Period (s)")

                plt.title("FTAN\n%s.%s - %s.%s"%(NET1, STA1, NET2, STA2))
                plt.show()

            GVdisp[i]["NET1"] = NET1
            GVdisp[i]["STA1"] = STA1
            GVdisp[i]["NET2"] = NET2
            GVdisp[i]["STA2"] = STA2
            GVdisp[i]["DISTKM"] = dist

            if diagramtype == 'PV':
                GVdisp[i]["PERIOD"] = per
                GVdisp[i]["GroupVel"] = disper
                GVdisp[i]["FREQ"] = 1./per
                GVdisp[i]["GroupTime"] = dist/disper

            elif diagramtype == 'FV':
                GVdisp[i]["PERIOD"] = 1./per
                GVdisp[i]["GroupVel"] = disper
                GVdisp[i]["FREQ"] = per
                GVdisp[i]["GroupTime"] = dist/disper

            elif diagramtype == 'FT':
                GVdisp[i]["PERIOD"] = 1./per
                GVdisp[i]["GroupVel"] = dist/disper
                GVdisp[i]["FREQ"] = per
                GVdisp[i]["GroupTime"] = disper

            elif diagramtype == 'PT':
                GVdisp[i]["PERIOD"] = per
                GVdisp[i]["GroupVel"] = dist/disper
                GVdisp[i]["FREQ"] = 1./per
                GVdisp[i]["GroupTime"] = disper
            # ICI
            Periods = np.arange(0, np.ceil(np.max(GVdisp[i]["PERIOD"]))+1.1, 0.1)
            dci=np.interp(Periods, GVdisp[i]["PERIOD"], GVdisp[i]["GroupVel"], left=np.nan, right=np.nan)
            dcii = np.interp(PER, Periods, dci, left=np.nan, right=np.nan)
            print(PER)
            print(dcii)
            df = pd.Series(dcii, index=PER, name="disp")
            df.plot()
            fn = os.path.join("TOMO_DISP", basename+".csv")
            df.to_csv(fn, header=[basename,])
            # LA
            maxp[i]=np.max(GVdisp[i]["PERIOD"])

        if PLOTRAWDISP:
            plt.figure()
            for i in range(len(GVdisp)):
                plt.plot(GVdisp[i]["PERIOD"], GVdisp[i]["GroupVel"], c='k', lw=0.5)
            plt.grid(True)
            plt.xlabel("Period (s)")
            plt.ylabel("Group Velocity (km/s)")
            plt.show()

        # Interpolation of the disp curv at specific periods

        Periods = np.arange(0, np.ceil(np.max(maxp))+1.1, 0.1)
        Disp = np.zeros((len(PER),len(GVdisp))) * np.nan

        for i in range(len(GVdisp)):
            P = GVdisp[i]["PERIOD"]
            DC = GVdisp[i]["GroupVel"]
            P, ip = np.unique(P, return_index=True)
            DC = DC[ip]
            if np.all(sorted(P) == P):
                dci=np.interp(Periods, P, DC, left=np.nan, right=np.nan)
            else:
                dci=np.zeros(len(Periods))*np.nan
            dcii = np.interp(PER, Periods, dci, left=np.nan, right=np.nan)

            Disp[:,i]=dcii

        if PLOTDISPALL:
            norm=np.sqrt(np.sum(np.isfinite(Disp),axis=1)) / np.max(np.sqrt(np.sum(np.isfinite(Disp),axis=1)))
            print(norm)
            plt.figure()
            plt.subplot(211)
            plt.plot(PER, Disp, lw=0.5)
            plt.errorbar(PER, np.nanmean(Disp, axis=1), np.nanstd(Disp, axis=1)/norm, lw=2)
            plt.ylim(Disp.min(), Disp.max())
            plt.xlim(0, np.max(maxp)+1)
            plt.grid(True)
            plt.ylabel('Group velocity (km/s)')

            plt.subplot(212)
            plt.plot(PER, np.sum(np.isfinite(Disp), axis=1), lw=2)
            plt.ylabel("# measurements")
            plt.xlabel('Period (s)')
            plt.xlim(0, np.max(maxp)+1)
            plt.grid(True)
            plt.show()

        if SAVEFILES:
            datapathout = "tmp/"
            prefix = "TestGroupVel_"
            sufix = "sGLISN.dat"

            for p in range(len(PER)):
                Dispfile = "%s%s%.1f%s"%(datapathout, prefix, PER[p], sufix)

                f = open(Dispfile,'w')

                for i in range(len(GVdisp)):
                    NET1 = GVdisp[i]["NET1"]
                    STA1 = GVdisp[i]["STA1"]
                    NET2 = GVdisp[i]["NET2"]
                    STA2 = GVdisp[i]["STA2"]
                    dist = GVdisp[i]["DISTKM"]

                    if not np.isnan(Disp[p,i]):
                        f.write( '%s %s %f %f %f %f\n'% (STA1,STA2,PER[p],Disp[p,i],0,dist))
                f.close()

if __name__ == "__main__":
    main()