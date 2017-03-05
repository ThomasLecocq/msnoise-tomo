import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import hot_r as viridis
from obspy.core import read

from .ftan_call import pickgroupdispcurv


def main():

    #SACfilelist = glob.glob(os.path.join(os.path.split(os.path.realpath(__file__))[0], r'*_Sym.SAC'))
    SACfilelist = glob.glob(r'TOMO_SAC/*.SAC')
    print(SACfilelist)
    PER=np.array([3.5, 5, 10, 20, 30, 50, 60]) # Interpolation periods
    PLOTDIAGR=1
    PLOTRAWDISP=1
    PLOTDISPALL=1
    SAVEFILES=1

    GVdisp = [{},]*len(SACfilelist)
    maxp = np.empty(len(SACfilelist))
    fmin = 0.1
    fmax = 1.0
    vgmin = 0.5
    vgmax = 3.0
    bmin = 0.00022
    bmax = 0.025
    diagramtype = 'PV'
    nfreq = 40
    ampmin = '0.05'

    for i, filename in enumerate(SACfilelist):
        NET1, STA1, NET2, STA2, crap = os.path.split(filename)[1].split('_')

        st = read(filename)
        dist = st[0].stats.sac.dist
        dt = st[0].stats.delta

        per, disper = pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                                        diagramtype, nfreq, ampmin, dist)
        if PLOTDIAGR:
            U = np.loadtxt('write_TV.txt')
            P = np.loadtxt('write_FP.txt')
            xmin = min(P)
            xmax = max(P)
            ymin = min(U)
            ymax = max(U)
            amp = np.loadtxt('write_amp.txt').T
            iu = np.where( (disper>=ymin) & (disper<=ymax) )
            per = per[iu]
            disper = disper[iu]
            Per, Vitg = np.meshgrid(P,U)
            plt.figure()
            plt.contourf(Per, Vitg, amp, 35, cmap=viridis)
            plt.colorbar()
            plt.contour(Per, Vitg, amp, 35, colors='k')

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
        datapathout = os.path.join(os.path.split(os.path.realpath(__file__))[0],"tmp/")
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