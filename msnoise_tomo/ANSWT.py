import os
import sys
import traceback
import zipfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from matplotlib import cm
from msnoise.api import connect, get_config
from scipy import ndimage
# from skimage import measure

from .fitellipse import fitellipse

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

def loadH(G,lmd):
    print("lambda", lmd)

    H=np.exp(-lmd*np.sum(np.abs(np.sign(G.toarray().T)),axis=1))
    print(H.shape)
    H = np.diag(H)
    print(H.shape)
    return H

def loadG(nX, nY, file, gridfile,ANSWT_toolbox_path):
    # if sys.platform[:3] == "win":
    #     command = os.path.join(ANSWT_toolbox_path, 'mk_MatPaths.exe') + " %s %s" %(file, gridfile)
    # else:
    #     command = os.path.join(ANSWT_toolbox_path, 'mk_MatPaths') + " %s %s" %(file, gridfile)
    # print(command)
    # os.system(command)
    from .intersect import mkpath
    from .lib.libmk_MatPaths import path
    path(file, gridfile)
    print('C code done')
    DATA = np.fromfile('matG.bin', dtype=np.float)
    G = DATA.reshape(-1,int(nX*nY))
    return G

def loadF(nX, nY, smoothfactor,gridfile,ANSWT_toolbox_path):
    # if sys.platform[:3] == "win":
    #     command = os.path.join(ANSWT_toolbox_path, 'mk_MatSmoothing.exe') + " %f %s " % (smooth, gridfile)
    # else:
    #     command = os.path.join(ANSWT_toolbox_path, 'mk_MatSmoothing') + " %f %s " % (smooth, gridfile)
    # print(command)
    # os.system(command)
    from .lib.libmkMatSmoothing import smooth
    smooth("%.7f"%smoothfactor, gridfile)
    F = np.fromfile('matF.bin', dtype=np.float32)
    F = F.reshape(int(nX*nY),int(nX*nY))
    return F

def LoadSmoothParam(paramFile):
    param=np.loadtxt(paramFile,skiprows=1, delimiter=" ")
    
    a1=param[0]
    b1=param[1]
    l1=param[2]
    s1=param[3]
    
    a2=param[4]
    b2=param[5]
    l2=param[6]
    s2=param[7]
    return [a1,b1,l1,s1,a2,b2,l2,s2]

def initModel(gridfile):
    grille=np.loadtxt(gridfile)
    nX=np.floor(np.floor((grille[0,1]-grille[0,0])/grille[2,0]))+1
    nY=np.floor(np.floor((grille[1,1]-grille[1,0])/grille[2,1]))+1
    [X0,Y0]=np.meshgrid(grille[0,0]+np.arange(0,nX)*grille[2,0], grille[1,0]+np.arange(0,nY)*grille[2,1])
    dx,dy=[grille[2,0], grille[2,1]]
    return [X0, Y0, nX,nY, dx, dy]

def ANSWT(gridfile,stacoordfile,DCfile,paramfile,PERIOD, show, v_cmap, d_cmap):
    path_tb=os.path.join(os.path.split(os.path.realpath(__file__))[0],'lib')
    print(path_tb)
    X,Y,nX,nY, dx, dy = initModel(gridfile)
    nX = int(nX)
    nY = int(nY)
    print(X)
    print("nX and Ny from gridfile & initModel", nX, nY)


    # Data loading


    #STA NET LAT LON EL
    STALOC = np.loadtxt(stacoordfile, dtype=str)
    print(STALOC)
    # STA1 STA2 PER VG eVG DIST
    DC = np.loadtxt(DCfile, dtype=str)
    # Data formatting

    nbsta=STALOC.shape[0]
    nbCorr=DC.shape[0]
    Stations = np.zeros((nbCorr, 2), dtype=int)

    for nbc in range(nbCorr):
        sta1 = DC[nbc][0]
        sta2 = DC[nbc][1]
        for nbs in range(nbsta):
            if sta1 == STALOC[nbs][0]:
                Stations[nbc,0]=nbs

            if sta2 == STALOC[nbs][0]:
                Stations[nbc,1]=nbs
    Vg = np.array(DC[:,3], dtype=float)
    dist = np.array(DC[:,5], dtype=float)

    print('Data loaded')
    alpha1,beta1,lambda1,Lcorr1,alpha2,beta2,lambda2,Lcorr2 = LoadSmoothParam(paramfile)
    F = loadF(nX,nY, Lcorr1,gridfile,path_tb) # Smoothing matrix
    F = scipy.sparse.lil_matrix(F)
    F1 = scipy.sparse.lil_matrix(F.T*F)
    print('F-Lcorr loaded, shape=', F1.shape)

    # Coordonnees des couples utilises et construction des rais
    x = np.array(STALOC[:,3], dtype=float)
    y = np.array(STALOC[:,2], dtype=float)

    MoyV_meas = np.median(Vg)
    VarV_meas = np.std(Vg)
    # Selection des tps d'arrivee compris entre +n et -n variance
    v_max=MoyV_meas+3*VarV_meas
    v_min=MoyV_meas-3*VarV_meas
    print(v_max, v_min)
    sv=np.where((Vg<=v_max) & (Vg>=v_min))[0]
    print('There are %i data in total.' % len(Vg))

    Vg1 = Vg[sv]
    dist1 = dist[sv]
    x1=x[Stations[sv,0]]
    x2=x[Stations[sv,1]]
    y1=y[Stations[sv,0]]
    y2=y[Stations[sv,1]]

    lpath = np.array([x1, y1, x2, y2]).T
    np.savetxt('lpath.txt', lpath, fmt="%.6f")
    print('lpath.txt written')

    dist=(np.hypot(x1-x2, y1-y2))

    t_obs = dist/Vg1
    nData1=lpath.shape[0]
    print('There are %i data to inverse.' % nData1)
    
    # OLD WAY:

    
    
    GGG=loadG(nX,nY,'lpath.txt',gridfile, path_tb)
    print("G.shape from old path code", GGG.shape)

    
    
    import pandas as pd
    # NEW WAY
    stations = pd.read_csv(stacoordfile, delimiter=' ', header=None)
    stations = np.array([row for id, row in stations.iterrows()])

    xmin = stations[:, 3].min() - dx * 2 - 0.005
    xmax = stations[:, 3].max() + dx * 2 - 0.005

    ymin = stations[:, 2].min() - dy * 2 - 0.005
    ymax = stations[:, 2].max() + dy * 2 + 0.005  

    from .intersect import mkpath
    # defined above
    # dx = 1.0
    # dy = 0.75
    print(xmin, xmax, dx, ymin, ymax, dy)
    # nX = int(np.floor(np.floor((xmax - xmin) / dx)) + 1)
    # nY = int(np.floor(np.floor((ymax - ymin) / dy)) + 1)

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
    
    izero = (np.sum(GGG, axis=1) == 0)
    t_obs[izero] = 0
    GGG=scipy.sparse.lil_matrix(GGG)
    print('GGG matrix loaded')
    GG = GGG.copy()

    # v0=vitesse moyenne (m/s)
    V0 = np.median(Vg1)
    print("V0",V0)
    # m0 = model initial, premiere inversion lissee pour eliminer les outliers
    m0 = V0*np.ones((nX,nY))
    print("GD.shape", GG.shape)
    print("m0.shape", m0.shape)
    T0 = GG*(1./m0.flatten())
    Dt = t_obs-T0
    G = 1./V0*GG
    # print G
    # print "*"*50
    H = loadH(G, float(lambda1))  #Damping matrix
    # print "*"*50
    # print H
    H = scipy.sparse.dia_matrix(H)
    print('H loaded')
    H1 = H.T*H
    invCd = scipy.sparse.dia_matrix(np.diag(np.ones(len(Dt))))
    print("h1 & F shape", H1.shape, F.shape)
    Q=(beta1*(H1)+alpha1*scipy.sparse.linalg.inv(F.astype(float)))
    # print 'Q', Q
    print('Q computed')

    #G1=G.T*G # possible to add data covariance matrix Cd (eq. 18 of Mordret et al. ,2013a)
    print(type(G), type(invCd), type(G))
    G1 = np.multiply(G.T, invCd) * G
    # plt.imshow(G1.toarray())
    # plt.colorbar()
    # plt.show()
    print(G1.shape)
    G1 = scipy.sparse.csc_matrix(G1)
    print(Q.shape, G1.shape, G.shape, Dt.shape, invCd.shape)
    print("*-"*50)
    print(type(G1+Q), (G1+Q).shape, G.T.shape)

    invG = scipy.sparse.linalg.spsolve((G1+Q), G.T)
    # plt.imshow(invG.toarray())
    # plt.colorbar()
    # plt.show()
    m_inv=invG*invCd*Dt
    print('m_inv computed')
    m_vel1 = m0.flatten()/(m_inv+1) # we defined m=(u0-u)/u = u0/u - 1 so u=u0/(m+1) (eq. 5)
    t_calc1=GG*(1./m_vel1)
    RMS_l1 = np.std(t_calc1-t_obs)
    print(RMS_l1)

    # 2eme iteration avec model lisse precedent
    # m02 = model initial No 2

    m02=m_vel1
    print(m02.shape)

    Dt3=t_obs-t_calc1

    MoyT_calc=np.mean(Dt3)
    VarT_calc=np.std(Dt3)
    # Selection des tps d'arrivee compris entre +n et -n variance
    t_max=MoyT_calc+2*VarT_calc
    t_min=MoyT_calc-2*VarT_calc
    s=np.where((Dt3<t_max) & (Dt3>t_min))[0]
    nData2=len(s)
    print('  The selection procedure removed ', nData1-nData2,' path(s) among ', nData1,' possible')
    Dt2=Dt3[s]
    Vgmesur=Vg1[s]
    GG=GGG[s,:]
    invCd = scipy.sparse.dia_matrix(np.diag(np.ones(len(Dt2))))
    m1 = 1./m02
    s0= scipy.sparse.lil_matrix(np.repeat(m1.reshape(1,-1), nData2, axis=0))
    G = s0.multiply(GG)

    H = loadH(G,  float(lambda1))
    H = scipy.sparse.dia_matrix(H)
    print('H loaded')
    H1 = H.T*H
    F = loadF(nX,nY, Lcorr2,gridfile,path_tb) # Smoothing matrix
    F = scipy.sparse.lil_matrix(F)
    F1 = scipy.sparse.lil_matrix(F.T*F)
    print('F-Lcorr loaded')
    # plt.figure()
    # plt.imshow(F.toarray(), vmin=0, vmax=1)
    # plt.colorbar()
    # plt.savefig("F.png")
    # plt.show()
    Q=(beta2*(H1)+alpha2*scipy.sparse.linalg.inv(F.astype(float)))
    G1 = np.multiply(G.T, invCd) * G
    invG = scipy.sparse.linalg.spsolve((G1+Q), G.T)
    m_inv2=invG*invCd*Dt2
    print('m_inv computed2')
    # Computation of the resolution matrix (can be long)
    calcresol = True
    if calcresol:
        print('Computation of the resolution matrix')
        MatRes=invG*invCd*G
        # plt.imshow(MatRes.toarray())
        # plt.colorbar()
        # plt.show()
        print('Resolution matrix computed')
    else:
        MatRes=[]

    # Model final
    m_vel2 = m02.flatten()/(m_inv2+1)
    t_calc2=GG*(1./m_vel2)
    M_vel= m_vel2.reshape(nY, nX) #  final model
    
    # plt.figure()
    # plt.imshow(M_vel)
    # plt.colorbar()
    # plt.savefig('mvel.png')
    # 
    # plt.figure()
    # plt.imshow(m_vel2.reshape(nX, nY))
    # plt.colorbar()
    # plt.savefig('mvel2.png')
    # 
    print("G.toarray()", G.toarray().shape)
    print("G.toarray().sum(axis=0)", G.toarray().sum(axis=0).shape)
    print("G.toarray().sum(axis=1)", G.toarray().sum(axis=1).shape)
    densitypath = np.sum(np.abs(np.sign(G.toarray())), axis=0).reshape(nY, nX) # number of ray per cell
    print(densitypath.shape)
    RMS_l2 = np.var((t_calc2-t_obs[s]))
    RMS0_l = np.var(Dt)
    print("Mvel computed")
    print(RMS0_l)
    # diminution de variance
    vared2 =100.*(1-RMS_l2/RMS0_l)
    print(vared2)

    doplot = 1
    if doplot:
        # Figures
        lonlim=[np.amin(X)-dx, np.amax(X)+dx]
        latlim=[np.amin(Y)-dy, np.amax(Y)+dy]
        span = 1 # Size of the averaging window
        window = np.ones((span,span))/(span*span)
        Dsity = ndimage.convolve(densitypath, window, mode='constant')
        # densitypath(Dsity==0)=NaN;
        seuil = 0
        id = np.where(Dsity <= seuil)
        print(id)
        M_vel[id] *= np.nan

        M = M_vel


        # Fig Path density
        plt.figure()
        plt.contourf(X+dx/2, Y+dy/2, densitypath, 30, origin='lower', cmap=d_cmap)
        cb = plt.colorbar()
        cb.set_label("Path density (#ray/pixel)")
        plt.contour(X+dx/2, Y+dy/2, Dsity, [1,], colors='k')
        plt.scatter(x,y,marker='^',c='k')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.title("Density Paths")
        plt.savefig("result_density_%.3fs.png"%PERIOD, dpi=300)

        # Tomo result

        np.savetxt('tomo_%.3fs.txt'%PERIOD, M)

        fig=plt.figure()
        vmin = np.mean(Vgmesur)-1.5*np.std(Vgmesur)
        vmax = np.mean(Vgmesur)+1.5*np.std(Vgmesur)

        cf = plt.contourf(X+dx/2, Y+dy/2, M, 30, origin='lower',
                     cmap=v_cmap)

        plt.scatter(x,y, marker='^',c='k')
        plt.contour(X+dx/2, Y+dy/2, Dsity, [1,], colors='w')
        plt.xlim(lonlim[0], lonlim[1])
        plt.ylim(latlim[0], latlim[1])
        plt.axis('off')
        plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
        plt.savefig("test.png", transparent=True, dpi=600)
        plt.close(fig)

        plt.figure()
        cf = plt.contourf(X+dx/2, Y+dy/2, M, 30, origin='lower',
                     cmap=v_cmap)
        plt.scatter(x,y, marker='^',c='k')
        plt.contour(X+dx/2, Y+dy/2, Dsity, [1,], colors='w')
        cb = plt.colorbar(cf)
        cb.set_label("Group velocity (km/s)")
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.title("Period = %.3f s, Vmean= %.3f km/s, VarRed= %.2f" %
                  (PERIOD, V0, vared2))
        plt.savefig("result_tomo_%.3fs.png"%PERIOD, dpi=300)

        with zipfile.ZipFile('tomo-result_%.3fs.kmz'%PERIOD, 'w') as z:
            z.writestr('doc.kml', kml.format(
                path='files/test.png',
                lat_north=latlim[1],
            lat_south=latlim[0],
                lon_east=lonlim[1],
            lon_west=lonlim[0],
                ))
            z.write('test.png', 'files/test.png')

        # Paths
        x11=x1[s]
        y11=y1[s]
        x21=x2[s]
        y21=y2[s]
        plt.figure()
        ax = plt.subplot(111)
        plt.scatter(x, y, marker='^',c='k')
        v = Vgmesur.T
        norm = mpl.colors.Normalize(vmin=np.min(v), vmax=np.max(v))
        # cmap = cm.jet_r
        m = cm.ScalarMappable(norm=norm, cmap=v_cmap)
        colors = m.to_rgba(v)
        for (a,b,c,d,C) in zip(x11, x21, y11, y21, colors):
            plt.plot([a,b],[c,d], color=C)
        m._A = []
        plt.colorbar(m)
        plt.contour(X+dx/2, Y+dy/2, Dsity, [1,], colors='k')
        plt.xlim(lonlim[0], lonlim[1])
        plt.ylim(latlim[0], latlim[1])
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.title("Period = %.3f s, Vmean= %.3f km/s" %
                  (PERIOD, V0))
        plt.savefig("result_paths_%.3fs.png"%PERIOD, dpi=300)
    if show:
        plt.show()

    if False:
        plt.figure()
        ax = plt.subplot(111)

        Resmin=np.nan*np.ones(MatRes.shape[0])
        Resmax=Resmin.copy()
        Resdir=Resmin.copy()
        Areares=Resmin.copy()
        X0res=Resmin.copy()
        Y0res=Resmin.copy()
        Xp=X.T
        Yp=Y.T
        dx = (X[0,1] - X[0,0])
        dy = (Y[1,1] - Y[0,0])
        for ir in range(MatRes.shape[0]):
            resol = MatRes[ir,:].toarray()
            if resol.max() < 1e-2:
                continue
            seuilresol = 0.4*resol.max()
            contours = measure.find_contours(resol.reshape(M_vel.shape), seuilresol)
            nbcontour = len(contours)
            ac = 0
            import shapely.geometry as sgeom
            largest = 0
            ilargest = 0
            for icont, contour in enumerate(contours):
                contour[:, 0] *= dx
                contour[:, 0] += X[0,0]

                contour[:, 1] *= dy
                contour[:, 1] += Y[0,0]

                c = sgeom.Polygon(contour)
                area = c.area*(111.12**2)
                if area > largest:
                    largest = area
                    ilargest = icont
                ac += area
            Areares[ir] = ac
            if len(contours[ilargest][:,0]) < 6:
                print("not enough point to fit an ellipse")
                Resmin[ir]=111.11*0.02
                Resmax[ir]=111.11*0.02
                Y0res[ir]=Yp[0][ir]
                X0res[ir]=Xp[0][ir]
                continue
            else:
                largest = contours[ilargest]
                try:
                    [zg, ag, bg, alphag] = fitellipse(contours[ilargest])
                    if ag:
                        Resmin[ir] = 111.11*ag
                        Resmax[ir] = 111.11*bg
                        Y0res[ir] = zg[0]
                        X0res[ir] = zg[1]
                except:
                    pass

        data1 = 2.*np.sqrt(Areares/np.pi)
        data1[data1<111.11*pasgrille*2]=111.11*pasgrille*2

        bla1 = X.T
        bla2 = Y.T
        data2=111.11*np.hypot(X0res-bla1.flatten().T, Y0res-bla2.flatten().T)
        data2 = data2.reshape(M_vel.shape)
        #data2[id]*=np.nan
        plt.contourf(X+dx/2, Y+dy/2, data2, 30, origin='lower',
                     cmap='jet_r')
        plt.show()


def main(per, a1, b1, l1, s1, a2, b2, l2, s2, filterid, comp, show):
    # Smoothing and damping parameters
    db = connect()
    alpha1 = a1 if a1 else float(get_config(db, "alpha1", plugin="Tomo"))
    beta1 = b1 if b1 else float(get_config(db, "beta1", plugin="Tomo"))
    lambda1 = l1 if l1 else float(get_config(db, "lambda1", plugin="Tomo"))
    sigma1 = s1 if s1 else float(get_config(db, "sigma1", plugin="Tomo"))

    alpha2 = a2 if a2 else float(get_config(db, "alpha2", plugin="Tomo"))
    beta2 = b2 if b2 else float(get_config(db, "beta2", plugin="Tomo"))
    lambda2 = l2 if l2 else float(get_config(db, "lambda2", plugin="Tomo"))
    sigma2 = s2 if s2 else float(get_config(db, "sigma2", plugin="Tomo"))
    
    v_cmap = get_config(db, "v_cmap", plugin="Tomo")
    d_cmap = get_config(db, "d_cmap", plugin="Tomo")
    
    
    if per is None:
        PER= get_config(db, "ftan_periods", plugin="Tomo")
        periods = np.array([float(pi) for pi in PER.split(',')])
    else:
        periods = [float(per),]

    # ANSWT inputs
    
    gridfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "Grid.dat")
    stacoordfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp, "STACoord.dat")
    
    for per in periods:
        DCfile=os.path.join("TOMO_FILES", "%02i" % filterid, comp, "TestGroupVel_%.3fs.dat"%float(per))
        PERIOD=per

        paramfile = os.path.join("TOMO_FILES", "%02i" % filterid, comp,'ParamFile.txt')
        fid = open(paramfile,'w');
        fid.write('%% alpha1 \t beta1 \t lambda1 \t Lcorr1 \t alpha2 \t beta2 \t lambda2 \t Lcorr2\n')
        fid.write('%f %f %f %f %f %f %f %f\n' %(alpha1,beta1,lambda1,sigma1,alpha2,beta2,lambda2,sigma2))
        fid.close()
        try:
            ANSWT(gridfile,stacoordfile,DCfile,paramfile,PERIOD, show, v_cmap, d_cmap)
        except:
            traceback.print_exc()
            print("!"*80)
            print("Can't compute tomo for period=", per)
            print("!"*80)

if __name__ == "__main__":
    main()
