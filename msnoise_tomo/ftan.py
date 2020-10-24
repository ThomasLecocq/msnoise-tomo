# from obspy.imaging.cm import viridis
import matplotlib.pyplot as plt
from .ftan_call import pickgroupdispcurv
from matplotlib.cm import hot_r as viridis
from msnoise.api import *
import shutil
import os
import numpy as np
import pandas as pd


def main(pair, bmin, bmax, show):
    PLOTDIAGR = show
    PLOTRAWDISP = 0
    PLOTDISPALL = 0
    SAVEFILES = 0

    if not os.path.isdir("TOMO_DISP"):
        os.makedirs("TOMO_DISP")

    if not os.path.isdir("RAW_FTAN_FILES"):
        os.makedirs("RAW_FTAN_FILES")

    # Get parameters from the tomoconfig panel
    db = connect()
    PER = get_config(db, "ftan_periods", plugin="Tomo")
    # parse the USER's periods
    PER = np.array([float(pi) for pi in PER.split(',')])
    fmin = float(get_config(db, "ftan_fmin", plugin="Tomo"))
    fmax = float(get_config(db, "ftan_fmax", plugin="Tomo"))
    nfreq = int(get_config(db, "ftan_nfreq", plugin="Tomo"))
    vgmin = float(get_config(db, "ftan_vgmin", plugin="Tomo"))
    vgmax = float(get_config(db, "ftan_vgmax", plugin="Tomo"))

    # If bmin,bmax not passed to main() then get from db
    if bmin is None:
        bmin = float(get_config(db, "ftan_bmin", plugin="Tomo"))
    if bmax is None:
        bmax = float(get_config(db, "ftan_bmax", plugin="Tomo"))

    diagramtype = get_config(db, "ftan_diagramtype", plugin="Tomo")
    ampmin = float(get_config(db, "ftan_ampmin", plugin="Tomo"))
    params = get_params(db)

    db = connect()
    comps = params.components_to_compute

    filters = get_filters(db)

    freqs = 1/PER
    print("The requested periods in the TOMOCONFIG are: %s" % PER)
    print("This corresponds to frequencies: %s" % freqs)

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
                for filter in filters:
                    for comp in comps:
                        fn = os.path.join("TOMO_SAC", "%02i" % filter.ref, comp,
                        "%s_%s_MEAN.sac" % (netsta1.replace('.', '_'), netsta2.replace('.', '_')))
                        print(fn)
                        if os.path.isfile(fn):
                            SACfilelist.append(fn)
                        else:
                            print("no file named", fn)
        else:
            for pi in pair:
                netsta1, netsta2 = pi.split('_')
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta1.replace('.','_'), netsta2.replace('.','_')))
                # SACfilelist.append(fn)
                # fn = os.path.join("TOMO_SAC", "%s_%s_REAL.SAC"%(netsta2.replace('.','_'), netsta1.replace('.','_')))
                # SACfilelist.append(fn)
                fn = os.path.join("TOMO_SAC", "%s_%s_MEAN.sac" % (
                    netsta1.replace('.', '_'), netsta2.replace('.', '_')))
                SACfilelist.append(fn)
                break
            pair = None

        print("Will process the following SAC files")
        print(SACfilelist)
        GVdisp = [{}, ]*len(SACfilelist)
        # maxp = np.empty(len(SACfilelist))

        # Prepare for the interpolated dispersion curves
        iper = np.argsort(PER)
        PER = PER[iper]
        Disp = np.zeros((len(PER), len(GVdisp))) * np.nan

        # Run the FTAN computation and automated dispersion curve picking
        for i, filename in enumerate(SACfilelist):

            NET1, STA1, NET2, STA2, crap = os.path.split(filename)[
                                                         1].split('_')

            st = read(filename)
            dist = st[0].stats.sac.dist
            # dt = st[0].stats.delta

            # Setup the output dispersion curve file (i.e. the TOMO_DISP/ file)
            GVdisp[i]["NET1"] = NET1
            GVdisp[i]["STA1"] = STA1
            GVdisp[i]["NET2"] = NET2
            GVdisp[i]["STA2"] = STA2
            GVdisp[i]["DISTKM"] = dist

            # pick the dispserion curve from fmin to fmax with nfreq elements
            pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                              diagramtype, nfreq, ampmin, dist)

            # copy the files from C++ output to files with basename prepended
            basename = "%s.%s_%s.%s_%s" % (NET1, STA1, NET2, STA2, crap)
            basename = basename.replace(".SAC", "")
            basename = os.path.join("RAW_FTAN_FILES", basename)
            for _ in ["write_amp.txt",
                      "write_disp.txt",
                      "write_FP.txt",
                      "write_ph.txt",
                      "write_TV.txt",
                      ]:
                shutil.move(_, _.replace("write", basename))

            # Load and sort dispersion curve
            D = np.loadtxt('%s_disp.txt' % basename)
            if D.ndim == 2:
                # sort based on the first column (period/frequency)
                isort = np.argsort(D[:, 0])
                D = D[isort]
            else:
                print("Only one dispersion pick...check data!!!")

            # Process the dispersion curve
            U = np.loadtxt('%s_TV.txt' % basename)
            P = np.loadtxt('%s_FP.txt' % basename)
            # get axes limits
            xmin = min(P)
            xmax = max(P)
            ymin = min(U)
            ymax = max(U)

            # Eliminate any  picks outside the y-limits
            iu = np.where((D[:, 1] >= ymin) & (D[:, 1] <= ymax))

            # filter only the picks that are within the y-limits
            # iu = np.where( (disper>=ymin) & (disper<=ymax) )
            # print(iu[0])
            if len(iu[0]) != 0:  # need to check if iu is empty.
                # per    = per[iu]
                # disper = disper[iu]
                D = D[iu]
                # Make sure you have the correct vectors; convert as needed
                if diagramtype == 'PV':
                    GVdisp[i]["PERIOD"] = D[:, 0]
                    GVdisp[i]["FREQ"] = 1. / D[:, 0]
                    GVdisp[i]["GroupVel"] = D[:, 1]
                    GVdisp[i]["GroupTime"] = dist / D[:, 1]
                elif diagramtype == 'FV':
                    GVdisp[i]["PERIOD"] = 1. / D[:, 0]
                    GVdisp[i]["FREQ"] = D[:, 0]
                    GVdisp[i]["GroupVel"] = D[:, 1]
                    GVdisp[i]["GroupTime"] = dist / D[:, 1]
                elif diagramtype == 'FT':
                    GVdisp[i]["PERIOD"] = 1. / D[:, 0]
                    GVdisp[i]["FREQ"] = D[:, 0]
                    GVdisp[i]["GroupVel"] = dist / D[:, 1]
                    GVdisp[i]["GroupTime"] = D[:, 1]
                elif diagramtype == 'PT':
                    GVdisp[i]["PERIOD"] = D[:, 0]
                    GVdisp[i]["FREQ"] = 1. / D[:, 0]
                    GVdisp[i]["GroupVel"] = dist / D[:, 1]
                    GVdisp[i]["GroupTime"] = D[:, 1]

                # Interpolate the dispersion curve to the periods requested in tomoconfig
                # and store in big matrix for later plotting
                Disp[:, i] = interpolate_disp_curve(GVdisp[i]["PERIOD"],
                                                    GVdisp[i]["GroupVel"], PER)
                write_tomo_disp_file(filename, basename, Disp[:, i], PER)
                # write_tomo_disp_file(filename, basename,
                # GVdisp[i]["GroupVel"],
                # GVdisp[i]["PERIOD"])

            else:
                # per    = []
                # disper = []
                print("No dispersion picks within vmin-vmax range")

            # LA
            # maxp[i] = np.max(GVdisp[i]["PERIOD"])

            if PLOTDIAGR:
                plot_FTAN_result(filename, basename,
                                 D[:, 0], D[:, 1], diagramtype)

        # Done with loop over the SAC files. Things below here
        # plot or process all dispersion curves at once using the
        # GVdisp variable.

        # Plot all measured dispersion curves on a single plot
        if PLOTRAWDISP:
            plot_raw_dispersion_curves(GVdisp)

        # Interpolation of the disp curv at specific periods

        # Periods = np.arange(0, np.ceil(np.max(maxp))+1.1, 0.1)
        # Disp = np.zeros((len(PER),len(GVdisp))) * np.nan

        # for i in range(len(GVdisp)):
        #     P = GVdisp[i]["PERIOD"]
        #     DC = GVdisp[i]["GroupVel"]
        #     P, ip = np.unique(P, return_index=True)
        #     DC = DC[ip]
        #     if np.all(sorted(P) == P):
        #         dci=np.interp(Periods, P, DC, left=np.nan, right=np.nan)
        #     else:
        #         dci=np.zeros(len(Periods))*np.nan
        #     dcii = np.interp(PER, Periods, dci, left=np.nan, right=np.nan)

            # Disp[:,i]=dcii

        if PLOTDISPALL:
            plot_interp_dispersion_curves(PER, Disp)

        if SAVEFILES:
            write_interp_disp_curve(PER, Disp, GVdisp)

# ==============================================================================


def interpolate_disp_curve(x, y, PER):
    # This function will write the TOMO_DISP file from the dataframe

    # sort the periods prior to np.interp()
    ix = np.argsort(x)
    x = x[ix]
    y = y[ix]

    # Map automatically picked dispersion picks to user requested periods
    dcii = np.interp(PER, x, y, left=np.nan, right=np.nan)
    return dcii
# ==============================================================================


def write_tomo_disp_file(filename, basename, dcii, PER):
    # This function will write the TOMO_DISP file from the dataframe

    df = pd.Series(dcii, index=PER, name="disp")  # store as a dataframe

    # Write the dispersion curve to TOMO_DISP file
    fn = filename.replace("TOMO_SAC", "TOMO_DISP").replace(
        ".sac", ".csv").replace(".SAC", ".csv")
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0])
    df.to_csv(fn, header=[basename, ], float_format='%.4f')
# ==============================================================================


def plot_FTAN_result(filename, basename, per, disper, diagramtype):
    # This function will plot the FTAN matrix and overlay the dispersion curve

    # get FTAN matrix
    amp = np.loadtxt('%s_amp.txt' % basename).T

    # Process the dispersion curve
    U = np.loadtxt('%s_TV.txt' % basename)
    P = np.loadtxt('%s_FP.txt' % basename)
    # get axes limits
    xmin = min(P)
    xmax = max(P)
    ymin = min(U)
    ymax = max(U)

    # setup matrix for contour plot
    Per, Vitg = np.meshgrid(P, U)
    plt.figure()
    plt.contourf(Per, Vitg, amp, 35, cmap=viridis)
    plt.colorbar()
    # plt.contour(Per, Vitg, amp, 35, colors='k')

    plt.plot(per, disper, '-ok', lw=1.5)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    # Set axes labels depending on diagramtype
    if diagramtype == 'PV':
        plt.xlabel("Period (s)")
        plt.ylabel("Velocity (km/s)")
    elif diagramtype == 'FV':
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Velocity (km/s)")
    elif diagramtype == 'FT':
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Time (s)")
    elif diagramtype == 'PT':
        plt.xlabel("Period (s)")
        plt.ylabel("Time (s)")

    NET1, STA1, NET2, STA2, crap = os.path.split(filename)[1].split('_')

    plt.title("FTAN\n%s.%s - %s.%s" % (NET1, STA1, NET2, STA2))
    plt.show()
# ==============================================================================


def plot_raw_dispersion_curves(GVdisp):
    # This function will plot the raw dispersion curves that are come from write_disp.txt.
    # This file comes from the C++ code.
    # All dispersion curves and associated data have been put into GVdisp

    plt.figure()
    for i in range(len(GVdisp)):
        plt.plot(GVdisp[i]["PERIOD"], GVdisp[i]["GroupVel"], c='k', lw=0.5)
    plt.grid(True)
    plt.xlabel("Period (s)")
    plt.ylabel("Group Velocity (km/s)")
    plt.title("Raw dispersion curves: fmin:fmax, with nfreq")
    plt.show()
# ==============================================================================


def plot_interp_dispersion_curves(PER, Disp):
    # This functions makes two plots. The first is a plot of the average
    # dispersion curve with error bars related to the standard deviation of
    # the measurements at each period. The second plot is bar plot showing
    # the total number of measurements at each period for the interpolated
    # dispersion curve data.

    norm = np.sqrt(np.sum(np.isfinite(Disp), axis=1)) / \
                   np.max(np.sqrt(np.sum(np.isfinite(Disp), axis=1)))

    # print("Number of group velocity measurments: %f"%norm)
    print(norm)

    plt.figure()
    plt.subplot(211)
    plt.plot(PER, Disp, lw=0.5)
    plt.errorbar(PER, np.nanmean(Disp, axis=1),
                 np.nanstd(Disp, axis=1)/norm, lw=2)
    # plt.ylim(Disp.min(), Disp.max())
    # plt.xlim(0, np.max(PER))
    plt.grid(True)
    plt.ylabel('Group velocity (km/s)')

    plt.subplot(212)
    plt.plot(PER, np.sum(np.isfinite(Disp), axis=1), lw=2)
    plt.ylabel("# measurements")
    plt.xlabel('Period (s)')
    # plt.xlim(0, np.max(PER))
    plt.grid(True)
    plt.show()
# ==============================================================================


def write_interp_disp_curve(PER, Disp, GVdisp):
    # Write a file for each period requrest by the user.
    # Each row in the file is a dispersion measurement between two stations
    # PER = all of the periods requested by the user in tomoconfig
    # Disp = [np,nd] group velocities at each period (np) for each station pair (nd)

    datapathout = "GROUPVEL_FILES/"
    prefix = "GroupVel_"
    sufix = "s.dat"

    # make the tmp/ directory
    if not os.path.isdir(datapathout):
        os.makedirs(datapathout)

    # write a file for each period
    for p in range(len(PER)):
        Dispfile = "%s%.4f%s" % (prefix, PER[p], sufix)
        Dispfile = os.path.join(datapathout, Dispfile)
        f = open(Dispfile, 'w')  # create the new file
        f.write('%s %s %s %s %s %s\n' %
                ('STA1', 'STA2', 'PER', 'GRPV', 'STD', 'DISTKM'))
        for i in range(len(GVdisp)):
            NET1 = GVdisp[i]["NET1"]
            STA1 = GVdisp[i]["STA1"]
            NET2 = GVdisp[i]["NET2"]
            STA2 = GVdisp[i]["STA2"]
            dist = GVdisp[i]["DISTKM"]

            STA1 = "%s.%s" % (NET1, STA1)
            STA2 = "%s.%s" % (NET2, STA2)
            # Only write non-empty data
            if not np.isnan(Disp[p, i]):
                f.write('%s %s %f %f %f %f\n' %
                        (STA1, STA2, PER[p], Disp[p, i], 0.0, dist))
        f.close()  # close the file
        print("Finished writing %s" % Dispfile)


# ==============================================================================
if __name__ == "__main__":
    main()
