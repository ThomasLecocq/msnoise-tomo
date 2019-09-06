# intersect.py
#
# Demonstrate how Shapely can be used to analyze and plot the intersection of
# a trajectory and regions in space.
import itertools
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
from shapely.geometry import LineString, Polygon, MultiPolygon


def mkpath(xmin, xmax, dx, ymin, ymax, dy, stations, greatcircle=False):
    if greatcircle:
        wgs84 = pyproj.Geod(ellps='WGS84')
    polygons = []
    # ids = {}
    X = np.arange(xmin, xmax, dx)
    Y = np.arange(ymin, ymax, dy)
    print("nx, ny", len(X),len(Y))
    # i = 0
    for ix, x in enumerate(X):
        for iy, y in enumerate(Y):
            polygons.append(Polygon([(x,y),(x+dx,y),(x+dx,y+dy),(x,y+dy)]))
            # ids[i] = [ix,iy]
            # i+=1
    patches = MultiPolygon(polygons=polygons)
    N = len(stations)
    G = np.zeros((N*(N-1)//2, len(X)*len(Y)))
    
    t = time.time()
    n = 0
    for sta1, sta2 in itertools.combinations(stations, 2):
        t2 = time.time()
        #~ print(sta1[0], sta2[0])
        x1 = sta1[3]
        x2 = sta2[3]
        y1 = sta1[2]
        y2 = sta2[2]
        if greatcircle:
            path = wgs84.npts(lon1=x1, lat1=y1,lon2=x2, lat2=y2, npts=25)
            path.insert(0, (x1,y1))
            path.append((x2,y2))
            vector = LineString(path)
        else:
            vector = LineString(((x1,y1), (x2, y2)))
        # intercepts = [[id, patch] for id, patch in enumerate(patches.geoms) if vector.intersects(patch)]
        # for id, patch in intercepts:
        #     G[n, id] += vector.intersection(patch).length
        for id, patch in enumerate(patches.geoms):
            if vector.intersects(patch):
                G[n, id] += vector.intersection(patch).length
        
        # x,y  = vector.xy
        # plt.plot(x,y,c='r',lw=1,zorder=10)
        # print("per pair: time taken", time.time() - t2)
        n+=1
    print("time taken", time.time() - t)
    return G





if __name__ == "__main__":
    gridfile=r"DATAFiles/GLISNGrid.txt"
    stations = pd.read_csv(r'DATAFiles/GLISN_STACoord.dat', delimiter=' ', header=None)
    print(stations[:1])
    stations = np.array([row for id,row in stations.iterrows()])
    
    xmin = stations[:,3].min()
    xmax = stations[:,3].max()
    
    ymin = stations[:,2].min()
    ymax = stations[:,2].max()
    
    grille=np.loadtxt(gridfile)
    dx = 1.0
    dy = 0.75
    print(xmin, xmax)
    nX=int(np.floor(np.floor((xmax-xmin)/dx))+1)
    nY=int(np.floor(np.floor((ymax-ymin)/dy))+1)
    
    G2 = mkpath(xmin, xmax, dx, ymin, ymax, dy, stations, True)
    
    plt.title("Sum of %i paths"%G2.shape[0])
    G2 = G2.sum(axis=0)
    G2 = G2.reshape(nX, nY)
    plt.imshow(G2.T, extent=(xmin, xmax+dx, ymin, ymax+dy/2.), 
               origin='lower', aspect='auto', interpolation="none",
               cmap='inferno')
    cb = plt.colorbar(orientation='vertical', shrink=0.7)
    cb.set_label('Path length in cell')
    
    stations = stations[:,2:4].astype(float)
    plt.scatter(stations[:,1], stations[:,0])
    plt.xlim(xmin, xmax+dx)
    plt.ylim(ymin, ymax+dy/2.)
    plt.show()