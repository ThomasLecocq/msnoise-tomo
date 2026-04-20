import os
# import sys
import traceback
import zipfile
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy.sparse
import scipy.io
import scipy.sparse.linalg
from matplotlib import cm
from msnoise.api import connect, get_config
from scipy import ndimage
from skimage import measure

from .fitellipse import fitellipse

import logging

from pyproj import Proj
_projections = {}  # make empy dictionary

kml = """\
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<GroundOverlay>
    <name>Tomo Result</name>
    <Icon>
        <href>{path}</href>
        <viewBoundScale>0.75</viewBoundScale>
    </Icon>
    <LatLonBox>
        <north>{lat_north}</north>
        <south>{lat_south}</south>
        <east>{lon_east}</east>
        <west>{lon_west}</west>
    </LatLonBox>
</GroundOverlay>
</kml>
"""


def loadH(G, lmd):

    print("lambda", lmd)
    H = np.exp(-lmd * np.sum(np.abs(np.sign(G.toarray().T)), axis=1))
    print(H.shape)
    H = np.diag(H)
    print(H.shape)
    return H


def loadG(nX, nY, file, gridfile):

    from .lib.libmk_MatPaths import path
    path(file, gridfile)
    DATA = np.fromfile('matG.bin', dtype=np.float)
    G = DATA.reshape(-1, int(nX * nY))
    return G


def loadF(nX, nY, smoothfactor, gridfile):

    from .lib.libmkMatSmoothing import smooth
    smooth("%.7f" % smoothfactor, gridfile)
    F = np.fromfile('matF.bin', dtype=np.float32)
    F = F.reshape(int(nX * nY), int(nX * nY))
    return F


def LoadSmoothParam(paramFile, tomofile):

    param = np.loadtxt(paramFile, skiprows=1, delimiter=" ")

    # For the first iteration
    a1 = param[0]
    b1 = param[1]
    l1 = param[2]  # lambda
    s1 = param[3]  # sigma --> Lcorr
    # For the second iteration
    a2 = param[4]
    b2 = param[5]
    l2 = param[6]
    s2 = param[7]

    logging.info("Smoothing parameter pairs")
    logging.info("Iteration 1 (alpha1,sigma1): (%f, %f)" % (a1, s1))
    logging.info("Iteration 2 (alpha2,sigma2): (%f, %f)" % (a2, s2))
    logging.info("Damping parameter pairs")
    logging.info("Iteration 1 (beta1,lambda1): (%f, %f)" % (b1, l1))
    logging.info("Iteration 2 (beta2,lambda2): (%f, %f)" % (b2, l2))

    return [a1, b1, l1, s1, a2, b2, l2, s2]


def initModel(gridfile, tomofile):

    grid = np.loadtxt(gridfile)
    print("Finished reading %s" % gridfile)
    print("Converting the grid corners to UTM.")
    print("Determing grid type to use.")

    # Notes: on the grid file format
    lonmin = grid[0,0]  # [deg] lower left corner
    latmin = grid[1,0]  # [deg] lower left corner
    lonmax = grid[0,1]  # [deg] upper right corner
    latmax = grid[1,1]  # [deg] upper right corner
    # dx = grid[2, 0]  # [deg] grid step sizes
    # dy = grid[2, 1]  # [deg] grid step sizes

    # Convert grid to UTM [m] and check if in single zone or not
    zmin, lmin, xmin, ymin = project((lonmin, latmin))  # lower left corner (lonmin,latmin)
    zmax, lmax, xmax, ymax = project((lonmax, latmax))  # upper right corner (lonmax,latmax)
    print("Lower left grid corner UTM zone: %d%s" % (zmin, lmin))
    print("Upper right grid corner UTM zone: %d%s" % (zmax, lmax))

    # mean Earth radius
    R = 6371
    # grid_type = "local"  # can force local computation by setting this variable
    grid_type = []

    # Choose which coordinate system to use (UTM, if single zone, otherwise local cartesian)
    if zmin != zmax or lmin != lmax or grid_type == "local":
        print("Tomography grid covers two UTM zones: (%d%s, %d%s)" % (zmin, lmin,zmax, lmax))
        print("Creating local cartesian grid centered in middle of grid.")
        grid_type = "local"

        lat0 = latmin + (latmax - latmin) / 2  # [deg] grid center
        lon0 = lonmin + (lonmax - lonmin) / 2  # [deg] grid center

        # Convert LL grid corner
        dlat = (latmin - lat0) * (np.pi/180)
        dlon = (lonmin - lon0) * (np.pi/180)
        # Convert to cartesian in [km]
        xmin = R * dlon * np.cos(dlat)
        ymin = R * dlat

        # Convert UR grid corner
        dlat = (latmax - lat0) * (np.pi/180)
        dlon = (lonmax - lon0) * (np.pi/180)
        # Convert to cartesian in [km]
        xmax = R * dlon * np.cos(dlat)
        ymax = R * dlat

    else:
        print("Tomography grid covers one UTM zone.")
        print("Using UTM coordinates for the cartesian grid.")
        grid_type = "utm"

        # Meters to km conversion for the UTM coordinates
        xmin = xmin / 1000
        xmax = xmax / 1000
        ymin = ymin / 1000
        ymax = ymax / 1000


    # Degrees to km for distance in degree with earth radius
    dx = R * grid[2, 0] * (np.pi / 180)  # [km] x-step size
    dy = R * grid[2, 1] * (np.pi / 180)  # [km] y-step size

    # Setup the grid in kilometers
    nX = np.floor(np.floor((xmax - xmin) / dx)) + 1
    nY = np.floor(np.floor((ymax - ymin) / dy)) + 1
    nX = int(nX)
    nY = int(nY)
    [X0, Y0] = np.meshgrid(xmin + np.arange(0, nX) * dx,
                           ymin + np.arange(0, nY) * dy)

    print("nX and nY from gridfile & initModel", nX, nY)

    # Write the new XYgrid file
    np.savetxt("XYgrid.dat", ((xmin,xmax),(ymin,ymax),(dx,dy)), delimiter=' ',fmt='%0.6f')
    # Note that the makeG and makeF functions have now hardcoded this filename, instead of 'gridfile' variable

    tomofile.write("Wrote XYgrid.dat\n")
    tomofile.write("x-min [km]: %0.6f\n" % xmin)
    tomofile.write("x-max [km]: %0.6f\n" % xmax)
    tomofile.write("dx [km]: %0.6f\n" % dx)
    tomofile.write("nx: %d\n" % nX)
    tomofile.write("y-min [km]: %0.6f\n" % ymin)
    tomofile.write("y-max [km]: %0.6f\n" % ymax)
    tomofile.write("dy [km]: %0.6f\n" % dy)
    tomofile.write("ny: %d\n" % nY)

    return [X0, Y0, nX, nY, dx, dy, grid_type]


def load_coorinate_file(stacoordfile, tomofile, grid_type, gridfile):
    # STA NET LAT LON EL
    STALOC = np.loadtxt(stacoordfile, dtype=str)
    nbsta = STALOC.shape[0]
    tomofile.write("Number of stations: %d\n" % nbsta)
    tomofile.write("NET.STA \t Network \t Northing \t Easting \t Elevation\n")

    # print("Station Information")
    # print("NET.STA \t NET \t LAT \t LON \t ELE")
    # print(STALOC)

    # grid_type = "local"  # can force local grid here

    if grid_type == "utm":
        for ii, station in enumerate(STALOC):
            # compute each station's UTM location
            z, l, x, y = project((float(STALOC[ii][3]), float(STALOC[ii][2])))
            STALOC[ii][3] = x/1000
            STALOC[ii][2] = y/1000
            tomofile.write("%s \t %s \t %s \t %s \t %s\n" % (STALOC[ii][0],
                                                             STALOC[ii][1],
                                                             y,
                                                             x,
                                                             STALOC[ii][4]))
    elif grid_type == "local":
        grid = np.loadtxt(gridfile)
        lonmin = grid[0,0]  # [deg] lower left corner
        latmin = grid[1,0]  # [deg] lower left corner
        lonmax = grid[0,1]  # [deg] upper right corner
        latmax = grid[1,1]  # [deg] upper right corner
        lat0 = latmin + (latmax - latmin) / 2  # [deg] grid center
        lon0 = lonmin + (lonmax - lonmin) / 2  # [deg] grid center

        lat = STALOC[:,2].astype(float)
        lon = STALOC[:,3].astype(float)

        # Convert station coordinates
        dlat = (lat - lat0) * (np.pi/180)
        dlon = (lon - lon0) * (np.pi/180)
        # Convert to cartesian in [km]
        R = 6371  # [km] mean Earth radius
        x = R * dlon * np.cos(dlat)
        y = R * dlat
        STALOC[:,3] = x
        STALOC[:,2] = y

        for ii, station in enumerate(STALOC):
            tomofile.write("%s \t %s \t %s \t %s \t %s\n" % (STALOC[ii][0],
                                                                   STALOC[ii][1],
                                                                   STALOC[ii][2],
                                                                   STALOC[ii][3],
                                                                   STALOC[ii][4]))

    # print("Station Information")
    # print("NET.STA \t NET \t Northing \t Easting \t ELE")
    # print(STALOC)

    return STALOC, nbsta


def load_dcfile(DCfile, STALOC, tomofile):
    nbsta = STALOC.shape[0]

    # STA1 STA2 PER VG eVG DIST
    DC = np.loadtxt(DCfile, dtype=str)
    # Data formatting

    nbCorr = DC.shape[0]
    Stations = np.zeros((nbCorr, 2), dtype=int)

    # Set the station indices for each data point
    for nbc in range(nbCorr):
        sta1 = DC[nbc][0]
        sta2 = DC[nbc][1]
        for nbs in range(nbsta):
            if sta1 == STALOC[nbs][0]:
                Stations[nbc, 0] = nbs

            if sta2 == STALOC[nbs][0]:
                Stations[nbc, 1] = nbs

    # p = np.array(DC[:, 2], dtype=float)  # [s] period
    Vg = np.array(DC[:, 3], dtype=float)  # [km/s] group velocity
    # error = np.array(DC[:, 4], dtype=float)  # [km/s] group velocity error
    dist = np.array(DC[:, 5], dtype=float)  # [km] interstation distance
    print('Dispersion data loaded')

    tomofile.write("Data loaded: %d\n" % nbCorr)

    return Stations, Vg, dist


def winnow_velocity_data(Vg, dist, var_lim, tomofile):
    MoyV_meas = np.median(Vg)
    VarV_meas = np.std(Vg)

    # Selection des tps d'arrivee compris entre +n et -n variance
    v_max = MoyV_meas + var_lim * VarV_meas
    v_min = MoyV_meas - var_lim * VarV_meas

    print('There are %i data in total.' % len(Vg))
    print("Winnowing based on group velocity with")
    print("V_min=%f [km/s] and V_max=%f [km/s]." % (v_min, v_max))

    sv = np.where((Vg <= v_max) & (Vg >= v_min))[0]
    Vg1 = Vg[sv]
    dist1 = dist[sv]

    print('There are %i data after winnowing.' % len(Vg1))

    tomofile.write("Removing data based on velocity outliers\n")
    tomofile.write("Num. data passed: %d\n" % len(Vg1))

    return Vg1, sv, dist1


def winnow_travel_time_data(dt3, nData1, Vg1, GGG, tomofile):
    MoyT_calc = np.mean(dt3)  # mean of the residuals
    VarT_calc = np.std(dt3)  # standard deviation

    # Selection of arrival times between +n and -n variance
    n = 2
    t_max = MoyT_calc + n * VarT_calc
    t_min = MoyT_calc - n * VarT_calc
    s = np.where((dt3 < t_max) & (dt3 > t_min))[0]
    nData2 = len(s)
    print('The selection procedure removed %d path(s) among %d possible data' % (nData1 - nData2, nData1))
    print('%d data now' % nData2)
    dt2 = dt3[s]  # remove these data points
    Vgmesur = Vg1[s]  # remove these velocities
    GG = GGG[s, :]  # remove these rows
    logging.info("Removed those %d rows of G-matrix" % (nData1 - nData2))

    tomofile.write("Removing outliers\n")
    tomofile.write("Num. data passed: %d\n" % nData2)

    return dt2, Vgmesur, GG, s


def make_smoothing_matrix(nX, nY, Lcorr1, gridfile):
    # Smoothing matrix
    F = loadF(nX, nY, Lcorr1, gridfile)
    F = scipy.sparse.lil_matrix(F)
    F1 = scipy.sparse.lil_matrix(F.T * F)
    print('F-Lcorr loaded, shape=', F1.shape)

    return F


def compute_raypaths(STALOC, Stations, sv):
    # Coordonnees des couples utilises et construction des rais
    x = np.array(STALOC[:, 3], dtype=float)
    y = np.array(STALOC[:, 2], dtype=float)
    # Station pair coordinates
    x1 = x[Stations[sv, 0]]
    y1 = y[Stations[sv, 0]]
    x2 = x[Stations[sv, 1]]
    y2 = y[Stations[sv, 1]]

    lpath = np.array([x1, y1, x2, y2]).T
    np.savetxt('lpath.txt', lpath, fmt="%.6f")
    print('lpath.txt written')
    nData1 = lpath.shape[0]
    print('There are %i data to invert.' % nData1)

    # compute distance between stations
    dist = (np.hypot(x1 - x2, y1 - y2))

    return x, y, dist


def compute_resolution_matrix(invG, invCd, G):
    # Computation of the resolution matrix (can be long)
    calcresol = True
    if calcresol:
        logging.info('Computing the resolution matrix')
        res_mat = invG * invCd * G
        np.savetxt('resolution.txt', res_mat.toarray(), fmt="%.12f")
        # plt.imshow(res_mat.toarray())
        # plt.colorbar()
        # plt.show()
        logging.info('Resolution matrix computed and saved: %s' % 'resolution.txt')
    else:
        res_mat = []

    return res_mat


def plot_ray_density(x, X, dx, nX, y, Y, dy, nY, G, PERIOD, d_cmap, plot_type):

    densitypath = np.sum(np.abs(np.sign(G.toarray())), axis=0).reshape(nY, nX)  # number of ray per cell
    logging.info("densitypath.shape: (%f, %f)" % densitypath.shape)

    # Smooth the raypath density image
    span = 1  # Size of the averaging window
    window = np.ones((span, span)) / (span * span)
    Dsity = ndimage.convolve(densitypath, window, mode='constant')
    # densitypath(Dsity==0)=NaN;

    # Make the path density figure
    if plot_type == "utm":
        # First need to convert things to units of meters
        X = X*1000; Y = Y*1000; dx = dx*1000; dy = dy*1000
        x = x*1000; y = y*1000

        plt.figure()
        plt.contourf(X + dx / 2, Y + dy / 2, densitypath, 30, origin='lower', cmap=d_cmap)
        cb = plt.colorbar()
        cb.set_label("Path density (#ray/pixel)")
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='w')
        plt.scatter(x, y, marker='^', c='k')
        plt.xlabel("Easting [m]")
        plt.ylabel("Northing [m]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))

    if plot_type == "local":

        plt.figure()
        plt.contourf(X + dx / 2, Y + dy / 2, densitypath, 30, origin='lower', cmap=d_cmap)
        cb = plt.colorbar()
        cb.set_label("Path density (#ray/pixel)")
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='w')
        plt.scatter(x, y, marker='^', c='k')
        plt.xlabel("Local x-coordinate [km]")
        plt.ylabel("Local y-coordinate [km]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))

    fig = plt.gcf()
    fig.set_size_inches(10,8)
    ax.axis('equal')
    plt.grid()
    plt.title("Period=%.4f s" % PERIOD)
    plt.savefig("result_density_%.4fs.png" % PERIOD, dpi=300)

    return Dsity

def plot_velocity_model(x, X, dx, y, Y, dy, M_vel, V0, vared2, Dsity, PERIOD, v_cmap, plot_type):

    # Set NaN grid cells with no rays passing through them
    M = set_vel_model_nans(Dsity, M_vel)

    if plot_type == "utm":
        # First need to convert things to units of meters
        X = X*1000; Y = Y*1000; dx = dx*1000; dy = dy*1000
        x = x*1000; y = y*1000

        plt.figure()
        cf = plt.contourf(X + dx/2, Y + dy/2, M, 30, origin='lower', cmap=v_cmap)
        cb = plt.colorbar(cf, format="%0.3f")
        cb.set_label("Group velocity (km/s)")
        plt.scatter(x, y, marker='^', c='k')  # plot stations (x,y)
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='k')
        plt.xlabel("Easting [m]")
        plt.ylabel("Northing [m]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))

    if plot_type == "local":

        plt.figure()
        cf = plt.contourf(X + dx/2, Y + dy/2, M, 30, origin='lower', cmap=v_cmap)
        cb = plt.colorbar(cf, format="%0.3f")
        cb.set_label("Group velocity (km/s)")
        plt.scatter(x, y, marker='^', c='k')  # plot stations (x,y)
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='k')
        plt.xlabel("Local x-coordinate [km]")
        plt.ylabel("Local y-coordinate [km]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))

    # Set title and save PNG file
    fig = plt.gcf()
    fig.set_size_inches(10,8)
    ax.axis('equal')
    plt.grid()
    plt.title("Period=%.4f s, V0=%.3f km/s, VarRed=%.2f%%" % (PERIOD, V0, vared2))
    plt.savefig("result_tomo_%.4fs.png" % PERIOD, dpi=300)


def write_tomo_kmz(x, X, dx, y, Y, dy, M_vel, Dsity, PERIOD, v_cmap):

    # TODO need to fix this function now that we do grid in kilometers

    M = set_vel_model_nans(Dsity, M_vel)

    lonlim = [np.amin(X) - dx, np.amax(X) + dx]
    latlim = [np.amin(Y) - dy, np.amax(Y) + dy]

    fig = plt.figure()
    cf = plt.contourf(X + dx/2, Y + dy/2, M, 30, origin='lower', cmap=v_cmap)
    plt.scatter(x, y, marker='^', c='k')
    plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='w')
    plt.xlim(lonlim[0], lonlim[1])
    plt.ylim(latlim[0], latlim[1])
    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.savefig("test.png", transparent=True, dpi=600)
    plt.close(fig)

    with zipfile.ZipFile('tomo-result_%.4fs.kmz' % PERIOD, 'w') as z:
        z.writestr('doc.kml', kml.format(
            path='files/test.png',
            lat_north=latlim[1],
            lat_south=latlim[0],
            lon_east=lonlim[1],
            lon_west=lonlim[0],
        ))
        z.write('test.png', 'files/test.png')

def plot_raypath_velocities(x, X, dx, y, Y, dy, Stations, sv, s, V0, Vgmesur, Dsity, PERIOD, v_cmap, plot_type):

    # The actual stations/data we used (first filter index, sv)
    x1 = x[Stations[sv, 0]]
    y1 = y[Stations[sv, 0]]
    x2 = x[Stations[sv, 1]]
    y2 = y[Stations[sv, 1]]
    # Paths from point 1 to point 2 (second filter index, s)
    x11 = x1[s]
    y11 = y1[s]
    x21 = x2[s]
    y21 = y2[s]

    # Setup the colorbar
    v = Vgmesur.T
    norm = mpl.colors.Normalize(vmin=np.min(v), vmax=np.max(v))
    m = cm.ScalarMappable(norm=norm, cmap=v_cmap)
    colors = m.to_rgba(v)
    m._A = []

    if plot_type == "utm":
        # First need to convert things to units of meters
        X = X*1000; Y = Y*1000; dx = dx*1000; dy = dy*1000
        x = x*1000; y = y*1000
        x11 = x11*1000; y11 = y11*1000; x21 = x21*1000; y21 = y21*1000

        plt.figure()
        for (a, b, c, d, C) in zip(x11, x21, y11, y21, colors):
            plt.plot([a, b], [c, d], color=C)
        cb = plt.colorbar(m, format="%0.3f")
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='k')
        plt.scatter(x, y, marker='^', c='k', zorder=5)
        plt.xlabel("Easting [m]")
        plt.ylabel("Northing [m]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%8.0f'))

    if plot_type == "local":

        plt.figure()
        for (a, b, c, d, C) in zip(x11, x21, y11, y21, colors):
            plt.plot([a, b], [c, d], color=C)
        cb = plt.colorbar(m, format="%0.3f")
        plt.contour(X + dx / 2, Y + dy / 2, Dsity, [1, ], colors='k')
        plt.scatter(x, y, marker='^', c='k', zorder=5)
        plt.xlabel("Local x-coordinate [km]")
        plt.ylabel("Local y-coordinate [km]")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))

    # Set title and save PNG file
    cb.set_label("Group velocity [km/s]")
    fig = plt.gcf()
    fig.set_size_inches(10,8)
    ax.axis('equal')
    plt.grid()
    plt.title("Period=%.4f s, V0=%.3f km/s" % (PERIOD, V0))
    plt.savefig("result_paths_%.4fs.png" % PERIOD, dpi=300)

def plot_residual_histogram(residual, dist, PERIOD):

    rmax = np.max(np.abs(residual))

    # 1D histogram bar plot
    plt.figure()
    plt.hist(residual, bins=20, range=(-rmax,rmax))
    plt.grid()
    plt.xlabel("Residual time [s]")
    plt.ylabel("Num. data")
    fig = plt.gcf()
    fig.set_size_inches(10,8)
    plt.title("Period=%.4f s" % PERIOD)
    plt.savefig("residuals_%.4fs.png" % PERIOD, dpi=300)

    # 2D heatmap
    plt.figure()
    plt.hist2d(residual, dist, bins=50, range=((-rmax,rmax),(0,np.max(dist))))
    plt.grid()
    plt.xlabel("Residual time [s]")
    plt.ylabel("Interstation distance [km]")
    cb = plt.colorbar(format="%0.3f")
    cb.set_label("Num. data")
    fig = plt.gcf()
    fig.set_size_inches(10,8)
    plt.title("Period=%.4f s" % PERIOD)
    plt.savefig("residuals_2D_%.4fs.png" % PERIOD, dpi=300)

def write_tomo_grid(M_vel, Dsity, PERIOD):

    M = set_vel_model_nans(Dsity, M_vel)

    # Tomo result
    filename = "tomo_%.4fs.txt" % PERIOD
    np.savetxt(filename, M)

def set_vel_model_nans(Dsity, M_vel):

    # find cell with no rays
    seuil = 0
    id = np.where(Dsity <= seuil)
    # print(id)
    # set the velocity model cells without rays to NaN
    M_vel[id] *= np.nan
    M = M_vel
    return M


# From https://gist.github.com/twpayne/4409500
def zone(coordinates):
    # coordinates = (lon, lat)
    if 56 <= coordinates[1] < 64 and 3 <= coordinates[0] < 12:
        return 32
    if 72 <= coordinates[1] < 84 and 0 <= coordinates[0] < 42:
        if coordinates[0] < 9:
            return 31
        elif coordinates[0] < 21:
            return 33
        elif coordinates[0] < 33:
            return 35
        return 37
    return int((coordinates[0] + 180) / 6) + 1


# From https://gist.github.com/twpayne/4409500
def letter(coordinates):
    # coordinates = (lon, lat)
    return 'CDEFGHJKLMNPQRSTUVWXX'[int((coordinates[1] + 80) / 8)]


# From https://gist.github.com/twpayne/4409500
def project(coordinates):
    # coordinates = (lon, lat)
    z = zone(coordinates)  # get utm zone for pyproj.Proj
    l = letter(coordinates)
    if z not in _projections:
        _projections[z] = Proj(proj='utm', zone=z, ellps='WGS84')
    x, y = _projections[z](coordinates[0], coordinates[1])
    if y < 0:
        y += 10000000
    return z, l, x, y

# From https://gist.github.com/twpayne/4409500
def unproject(z, l, x, y):
    if z not in _projections:
        _projections[z] = Proj(proj='utm', zone=z, ellps='WGS84')
    if l < 'N':
        y -= 10000000
    lon, lat = _projections[z](x, y, inverse=True)
    return (lon, lat)

def write_residuals(PERIOD, residual, dist):
    # write a CSV file of the residuals vs. interstation distance

    if len(residual) != len(dist):
        print("Residual vector not same length as distance vector.")

    file = "residuals_%.4fs.dat" % PERIOD
    file = open(file, "w")
    file.write("dist [km], residual [s]\n")
    for index in range(len(dist)):
        file.write("%0.3f, %0.5f\n" % (dist[index], residual[index]))
    file.close()

def ANSWT(gridfile, stacoordfile, DCfile, paramfile, PERIOD, show, v_cmap, d_cmap):
    # Create the tomography output file
    tomofile = "tomo_%.4fs.dat" % PERIOD
    tomofile = open(tomofile, "w")

    X, Y, nX, nY, dx, dy, grid_type = initModel(gridfile, tomofile)

    # Data loading
    STALOC, nbsta = load_coorinate_file(stacoordfile, tomofile, grid_type, gridfile)
    Stations, Vg, dist_dc = load_dcfile(DCfile, STALOC, tomofile)
    print("load_dcfile: dist")
    # print(dist_dc)

    # Winnow data outside of the var_lim standard deviations
    var_lim = 3
    Vg1, sv, dist_dc = winnow_velocity_data(Vg, dist_dc, var_lim, tomofile)
    V0 = np.median(Vg1)  # v0 = average velocity [km/s]
    tomofile.write("V0 [km/s]: %0.3f\n" % V0)  # write to the output file
    print("Median Group Velocity: V0=%0.3f [km/s]" % V0)
    m0 = V0 * np.ones((nX, nY))  # m0 = initial velocity model
    # print("m0.shape", m0.shape)
    logging.info("m0.shape: (%d, %d)" % m0.shape)

    # Get regularization terms
    alpha1, beta1, lambda1, Lcorr1, alpha2, beta2, lambda2, Lcorr2 = LoadSmoothParam(paramfile, tomofile)

    # Compute (x,y) pairs for data raypaths
    x, y, dist_hypot = compute_raypaths(STALOC, Stations, sv)

    # OLD WAY:
    logging.info("Loading G matrix from C-code")
    GGG = loadG(nX, nY, 'lpath.txt', "XYgrid.dat")
    logging.info("G.shape from old path code: (%d,%d) (matG.bin)" % GGG.shape)
    # print("G.shape from old path code: (%d,%d) (matG.bin)", GGG.shape())

    nData1 = GGG.shape[0]  # rumber of rows in G matrix

    # import pandas as pd
    # # NEW WAY
    # stations = pd.read_csv(stacoordfile, delimiter=' ', header=None)
    # stations = np.array([row for id, row in stations.iterrows()])

    # xmin = stations[:, 3].min() - dx * 2 - 0.005
    # xmax = stations[:, 3].max() + dx * 2 - 0.005
    #
    # ymin = stations[:, 2].min() - dy * 2 - 0.005
    # ymax = stations[:, 2].max() + dy * 2 + 0.005
    #
    # logging.info('Min lon:  %0.6f [deg]' % xmin)
    # logging.info('Max lon:  %0.6f [deg]' % xmax)
    # logging.info('Min lat:  %0.6f [deg]' % ymin)
    # logging.info('Max lat:  %0.6f [deg]' % ymax)
    # logging.info('dlon:     %0.6f [deg]' % dx)
    # logging.info('dlat:     %0.6f [deg]' % dy)

    # GGG = mkpath(xmin, xmax, dx, ymin, ymax, dy, stations, False)
    # GGG = GGG[sv]
    # print("G.shape from new path code", GGG.shape)
    # 
    # # tmp plot paths to be sure it's OK
    # G2 = GGG.sum(axis=0)
    # G2 = G2.reshape(nX, nY)
    # 
    # plt.figure()
    # plt.imshow(GGG1.sum(axis=0).reshape(nX, nY).T,
    #            extent=(xmin, xmax + dx, ymin, ymax + dy / 2.),
    #            origin='lower', aspect='auto', interpolation="none",
    #            cmap='inferno')
    # cb = plt.colorbar(orientation='vertical', shrink=0.7)
    # cb.set_label('Path length in cell')
    # 
    # plt.xlim(xmin, xmax + dx)
    # plt.ylim(ymin, ymax + dy / 2.)
    # plt.savefig("paths_old.png")
    # 
    # plt.figure()
    # plt.imshow(G2.T, extent=(xmin, xmax + dx, ymin, ymax + dy / 2.),
    #            origin='lower', aspect='auto', interpolation="none",
    #            cmap='inferno')
    # cb = plt.colorbar(orientation='vertical', shrink=0.7)
    # cb.set_label('Path length in cell')
    # 
    # plt.xlim(xmin, xmax + dx)
    # plt.ylim(ymin, ymax + dy / 2.)
    # plt.savefig("paths_new.png")
    # END NEW WAY


    # Remove any rows with zero length ray and the corresponding data (t_obs)
    # izero = (np.sum(GGG, axis=1) == 0)
    # t_obs[izero] = 0
    # print("Removed %d zero length rays." % np.sum(izero))

    GGG = scipy.sparse.lil_matrix(GGG)
    logging.info('GGG matrix loaded/made SPARSE')
    GG = GGG.copy()
    logging.info("GG.shape: (%d, %d)" % GG.shape)

    # sum along rows of G-matrix to get ray distance
    dist = GG.toarray().sum(axis=1) # distance in units of the grid
    # Compute the observed times in whatever units the grid is.
    t_obs = dist / Vg1  # [s] arrival time data
    # print("dist, dist_hypot, dist_dc")
    # print((dist, dist_hypot, dist_dc))

    print("\n----------------------------------------------------\n")
    print("Starting first iteration; velocity outliers removed!\n")
    print("----------------------------------------------------\n")

    # The first inversion is used to eliminate outliers
    T0 = GG * (1. / m0.flatten())  # compute predicted data for initial model
    dt = t_obs - T0  # compute travel-time residuals relative to reference model T0.
    G = GG * (1. / V0)  # normalize G by V0 (equation 14 in Mordret et al., 2013)
    # print("t_obs, T0, dt")
    # print(t_obs, T0, dt)

    # get the smoothing matrix
    F = make_smoothing_matrix(nX, nY, Lcorr1, "XYgrid.dat")

    # Compute the damping matrix
    H = loadH(G, float(lambda1))  # Damping matrix
    H = scipy.sparse.dia_matrix(H)  # convert to sparse representation
    logging.info('H loaded')
    H1 = H.T * H
    # print("h1 & F shape", H1.shape, F.shape)
    logging.info("H1.shape: (%d, %d)" % H1.shape)
    logging.info("F.shape: (%d, %d)" % F.shape)

    # compute the data covariance matrix
    invCd = scipy.sparse.dia_matrix(np.diag(np.ones(len(dt))))

    # Compute Q matrix as combination of damping and smoothing
    Q = (beta1 * (H1) + alpha1 * scipy.sparse.linalg.inv(F.astype(float)))
    logging.info('Q computed')

    # G1=G.T*G # possible to add data covariance matrix Cd (eq. 18 of Mordret et al. ,2013a)
    # print(type(G), type(invCd), type(G))
    logging.info("type(G): %s" % type(G))
    logging.info("type(G): %s" % type(invCd))

    G1 = np.multiply(G.T, invCd) * G
    # plt.imshow(G1.toarray())
    # plt.colorbar()
    # plt.show()
    # print(G1.shape)
    logging.info("G1.shape: (%d, %d)" % G1.shape)
    G1 = scipy.sparse.csc_matrix(G1)
    # print(Q.shape, G1.shape, G.shape, dt.shape, invCd.shape)
    logging.info("Q.shape: (%d, %d)" % Q.shape)
    logging.info("G1.shape: (%d, %d)" % G1.shape)
    logging.info("G.shape: (%d, %d)" % G.shape)
    logging.info("dt.shape: (%d)" % dt.shape)
    logging.info("invCd.shape: (%d, %d)" % invCd.shape)
    # print("*-" * 50)
    # print(type(G1 + Q), (G1 + Q).shape, G.T.shape)
    logging.info("type(G1 + Q): %s" % type(G1 + Q))
    logging.info("(G1 + Q).shape: (%d, %d)" % (G1 + Q).shape)
    logging.info("G.T.shape: (%d, %d)" % G.T.shape)

    invG = scipy.sparse.linalg.spsolve((G1 + Q), G.T)
    # plt.imshow(invG.toarray())
    # plt.colorbar()
    # plt.show()

    # Compute the linear inverse solution from the residuals
    m_inv = invG * invCd * dt
    logging.info('m_inv computed')

    m_vel1 = m0.flatten() / (m_inv + 1)  # we defined m=(u0-u)/u = u0/u - 1 so u=u0/(m+1) (eq. 5)
    t_calc1 = GG * (1. / m_vel1)  # compute new data from new velocity model

    dt3 = t_obs - t_calc1  # compute the new residuals from the smooth model
    RMS_l1 = np.std(dt3)  # compute the standard deviation of the residuals
    logging.info("Standard deviation of residuals: RMS_l1=%0.2e [s]" % RMS_l1)

    tomofile.write("Iteration 1:\n")
    tomofile.write("alpha1: %0.6f\n" % alpha1)
    tomofile.write("sigma1: %0.6f\n" % Lcorr1)
    tomofile.write("beta1: %0.6f\n" % beta1)
    tomofile.write("lambda1: %0.6f\n" % lambda1)
    tomofile.write("RMS1 [s]: %0.6f\n" % RMS_l1)

    print("\n----------------------------------------------------\n")
    print("Starting second iteration; outliers will be removed!\n")
    print("----------------------------------------------------\n")
    # 2nd iteration with the previous smooth model
    # m02 = model initial No 2
    m02 = m_vel1

    tomofile.write("Iteration 2:\n")

    # Remove outliers based on travel times
    dt2, Vgmesur, GG, s = winnow_travel_time_data(dt3, nData1, Vg1, GGG, tomofile)

    invCd = scipy.sparse.dia_matrix(np.diag(np.ones(len(dt2))))  # create
    m1 = 1. / m02  # create the slowness vector
    s0 = scipy.sparse.lil_matrix(np.repeat(m1.reshape(1, -1), len(dt2), axis=0))
    # Make the G matrix that is slowness times original G
    G = s0.multiply(GG)
    logging.info("G.toarray() dimension: (%d, %d)" % G.toarray().shape)
    logging.info("G.toarray().sum(axis=0): %d" % G.toarray().sum(axis=0).shape)
    logging.info("G.toarray().sum(axis=1): %d" % G.toarray().sum(axis=1).shape)

    # Get H matrix based on new G and lambda 1
    logging.info('Making H-matrix (damping) from C code for lambda: %0.2e' % lambda2)
    # H = loadH(G, float(lambda1))  # D.M. seems like it should be lambda2
    H = loadH(G, float(lambda2))
    H = scipy.sparse.dia_matrix(H)
    logging.info('H-matrix loaded and made SPARSE')
    H1 = H.T * H
    logging.info("Computed H'H")

    # Get F matrix based on Lcorr2 (sigma2)
    logging.info('Making F-matrix (smoothing) from C code for sigma: %0.2e' % Lcorr2)
    F = loadF(nX, nY, Lcorr2, "XYgrid.dat")  # Smoothing matrix
    F = scipy.sparse.lil_matrix(F)
    logging.info('F-sigma loaded')

    # Compute the Q matrix
    Q = (beta2 * (H1) + alpha2 * scipy.sparse.linalg.inv(F.astype(float)))

    # Compute inversion operator (equation 32 in Barmin et al., 2001)
    G1 = np.multiply(G.T, invCd) * G
    invG = scipy.sparse.linalg.spsolve((G1 + Q), G.T)  # (equation 18 in Mordret et al., 2013)

    m_inv2 = invG * invCd * dt2  # compute new model (equation 31 in Barmin et al., 2001)
    # (equation 17 in Mordret et al., 2013)

    # Model final
    m_vel2 = m02.flatten() / (1 + m_inv2)  # compute new model (equation 11 in Barmin et al., 2001)
    M_vel = m_vel2.reshape(nY, nX)  # final model
    logging.info("Final Mvel computed")

    t_calc2 = GG * (1. / m_vel2)
    residual = t_calc2 - t_obs[s]
    RMS_l2 = np.var(residual)
    RMS0_l = np.var(dt)

    logging.info("Variance of dt for background velocity model, RMS0_l = %0.2e " % RMS0_l)
    logging.info("Variance of dt for final velocity model, RMS_l2 = %0.2e " % RMS_l2)

    # compute the variance reduction
    vared2 = 100.0 * (1.0 - RMS_l2/RMS0_l)
    logging.info("Variance reduction from RMS0_l and RMS_l2: %0.2f%%" % vared2)

    tomofile.write("alpha2: %0.6f\n" % alpha2)
    tomofile.write("sigma2: %0.6f\n" % Lcorr2)
    tomofile.write("beta2: %0.6f\n" % beta2)
    tomofile.write("lambda2: %0.6f\n" % lambda2)
    tomofile.write("RMS2 [s]: %0.6f\n" % RMS_l2)
    tomofile.write("Var. Reduction [%%]: %0.6f\n" % vared2)

    # Write the residuals file
    dist = GG.toarray().sum(axis=1)  ## this is the final G-matrix after windowing
    write_residuals(PERIOD, residual, dist)

    print("Done with inversion...plotting results.")

    # Make figures even if the user did not ask for them to be shown
    doplot = 1
    plot_type = grid_type  # "utm" or "local"

    if doplot:

        # Make the raypath density plot
        Dsity = plot_ray_density(x, X, dx, nX, y, Y, dy, nY, G, PERIOD, d_cmap, plot_type)

        # Tomography model figure
        plot_velocity_model(x, X, dx, y, Y, dy, M_vel, V0, vared2, Dsity, PERIOD, v_cmap, plot_type)

        # Raypaths colored by velocities
        plot_raypath_velocities(x, X, dx, y, Y, dy, Stations, sv, s, V0, Vgmesur, Dsity, PERIOD, v_cmap, plot_type)

        # Residuals histogram
        plot_residual_histogram(residual, dist, PERIOD)

        # Write the KMZ file for GoogleEarth
        write_tomo_kmz(x, X, dx, y, Y, dy, M_vel, Dsity, PERIOD, v_cmap)

        # Write the velocity model to a text file
        write_tomo_grid(M_vel, Dsity, PERIOD)

        # vmin = np.mean(Vgmesur) - 1.5 * np.std(Vgmesur)
        # vmax = np.mean(Vgmesur) + 1.5 * np.std(Vgmesur)

    # show is the command line argument, default is show
    if show:
        plt.show()

    # Compute the resolution matrix (equation 34 in Barmin et al., 2001)
    # Each row in res_mat is a resolution map defining the resolution at one spatial node
    res_mat = compute_resolution_matrix(invG, invCd, G)  # (equation 19 in Mordret et al., 2013)

    # if False:
    #     plt.figure()
    #     ax = plt.subplot(111)
    # 
    #     Resmin = np.nan * np.ones(res_mat.shape[0])
    #     Resmax = Resmin.copy()
    #     Resdir = Resmin.copy()
    #     Areares = Resmin.copy()
    #     X0res = Resmin.copy()
    #     Y0res = Resmin.copy()
    #     Xp = X.T
    #     Yp = Y.T
    #     dx = (X[0, 1] - X[0, 0])
    #     dy = (Y[1, 1] - Y[0, 0])
    #     for ir in range(res_mat.shape[0]):
    #         resol = res_mat[ir, :].toarray()
    #         if resol.max() < 1e-2:
    #             continue
    #         seuilresol = 0.4 * resol.max()
    #         contours = measure.find_contours(resol.reshape(M_vel.shape), seuilresol)
    #         nbcontour = len(contours)
    #         ac = 0
    #         import shapely.geometry as sgeom
    #         largest = 0
    #         ilargest = 0
    #         for icont, contour in enumerate(contours):
    #             contour[:, 0] *= dx
    #             contour[:, 0] += X[0, 0]
    # 
    #             contour[:, 1] *= dy
    #             contour[:, 1] += Y[0, 0]
    # 
    #             c = sgeom.Polygon(contour)
    #             area = c.area * (111.12 ** 2)
    #             if area > largest:
    #                 largest = area
    #                 ilargest = icont
    #             ac += area
    #         Areares[ir] = ac
    #         if len(contours[ilargest][:, 0]) < 6:
    #             print("not enough point to fit an ellipse")
    #             Resmin[ir] = 111.11 * 0.02
    #             Resmax[ir] = 111.11 * 0.02
    #             Y0res[ir] = Yp[0][ir]
    #             X0res[ir] = Xp[0][ir]
    #             continue
    #         else:
    #             largest = contours[ilargest]
    #             try:
    #                 [zg, ag, bg, alphag] = fitellipse(contours[ilargest])
    #                 if ag:
    #                     Resmin[ir] = 111.11 * ag
    #                     Resmax[ir] = 111.11 * bg
    #                     Y0res[ir] = zg[0]
    #                     X0res[ir] = zg[1]
    #             except:
    #                 pass
    # 
    #     data1 = 2. * np.sqrt(Areares / np.pi)
    #     data1[data1 < 111.11 * pasgrille * 2] = 111.11 * pasgrille * 2
    # 
    #     bla1 = X.T
    #     bla2 = Y.T
    #     data2 = 111.11 * np.hypot(X0res - bla1.flatten().T, Y0res - bla2.flatten().T)
    #     data2 = data2.reshape(M_vel.shape)
    #     # data2[id]*=np.nan
    #     plt.contourf(X + dx / 2, Y + dy / 2, data2, 30, origin='lower', cmap='jet_r')
    #     plt.show()

    tomofile.close()  # close the output file


def main(per, a1, b1, l1, s1, a2, b2, l2, s2, filterid, comp, show):

    verbose = False
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # Smoothing and damping parameters
    db = connect()
    alpha1 = a1 if a1 is not None else float(get_config(db, "alpha1", plugin="Tomo"))
    beta1 = b1 if b1 is not None else float(get_config(db, "beta1", plugin="Tomo"))
    lambda1 = l1 if l1 is not None else float(get_config(db, "lambda1", plugin="Tomo"))
    sigma1 = s1 if s1 is not None else float(get_config(db, "sigma1", plugin="Tomo"))
    alpha2 = a2 if a2 is not None else float(get_config(db, "alpha2", plugin="Tomo"))
    beta2 = b2 if b2 is not None else float(get_config(db, "beta2", plugin="Tomo"))
    lambda2 = l2 if l2 is not None else float(get_config(db, "lambda2", plugin="Tomo"))
    sigma2 = s2 if s2 is not None else float(get_config(db, "sigma2", plugin="Tomo"))

    v_cmap = get_config(db, "v_cmap", plugin="Tomo")
    d_cmap = get_config(db, "d_cmap", plugin="Tomo")

    if per is None:
        PER = get_config(db, "ftan_periods", plugin="Tomo")
        periods = np.array([float(pi) for pi in PER.split(',')])
    else:
        periods = [float(per), ]

    # ANSWT inputs
    gridfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "Grid.dat")
    stacoordfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "STACoord.dat")

    for per in periods:
        DCfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "TestGroupVel_%.4fs.dat" % float(per))
        print("Processing %s" % DCfile)
        logging.info("Processing %s" % DCfile)
        PERIOD = per
        logging.info("Period from db_config: %.4f" % PERIOD)

        paramfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, 'ParamFile.txt')
        print("Writing parameters to %s" % paramfile)
        fid = open(paramfile, 'w')
        fid.write('%% alpha1 \t beta1 \t lambda1 \t Lcorr1 \t alpha2 \t beta2 \t lambda2 \t Lcorr2\n')
        fid.write('%f %f %f %f %f %f %f %f\n' % (alpha1, beta1, lambda1, sigma1, alpha2, beta2, lambda2, sigma2))
        fid.close()
        try:
            ANSWT(gridfile, stacoordfile, DCfile, paramfile, PERIOD, show, v_cmap, d_cmap)
        except:
            traceback.print_exc()
            print("!" * 80)
            print("Can't compute tomo for period=", per)
            print("!" * 80)


if __name__ == "__main__":
    main()
