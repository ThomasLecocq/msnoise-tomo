import os
import sys

import numpy as np


def pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, dist):
    if sys.platform[:3] == "win":
        ftan = os.path.join(os.path.split(os.path.realpath(__file__))[0],"lib", r"ftan.exe")
    else:
        ftan = os.path.join(os.path.split(os.path.realpath(__file__))[0],"lib", r"ftan")

    command = [ftan," ", filename, ' fmin=', fmin, ' fmax=', fmax,
        ' vgMin=', vgmin, ' vgMax=', vgmax, ' bmin=', bmin ,' bmax=', bmax,
        ' disp=none out=mat diag=', diagramtype, ' nfreq=', nfreq, ' ampMin=', ampmin]

    com = "".join([str(a) for a in command])
    print(com)
    os.system(com)

    V = np.loadtxt('write_TV.txt')
    P = np.loadtxt('write_FP.txt')
    amp=np.loadtxt('write_amp.txt')
    amp=amp.T

    index = np.unravel_index(np.argmax(amp), amp.shape)
    if len(index) == 2:
        iv, ip = index
    else:
        iv, ip = index[0]

    vinit = V[iv]
    pinit = P[ip]


    if diagramtype == 'PV':
        finit = "%.5f"% (1./pinit)
        vginit = vinit

    elif diagramtype == 'FV':
        finit = pinit
        vginit = vinit

    elif diagramtype == 'FT':
        finit = pinit
        vginit = dist/vinit

    elif diagramtype == 'PT':
        finit = 1. / pinit
        vginit = dist/vinit

    command = [ftan," ", filename, ' fmin=', fmin, ' fmax=', fmax,
        ' vgMin=', vgmin, ' vgMax=', vgmax, ' finit=', finit, ' vginit=',
               vginit, ' bmin=', bmin ,' bmax=', bmax,
        ' disp=cont out=mat diag=', diagramtype, ' nfreq=', nfreq,
               ' ampMin=', ampmin]

    com = "".join([str(a) for a in command])
    print(com)
    os.system(com)

    D = np.loadtxt('write_disp.txt')
    isort = np.argsort(D[:,0])
    D = D[isort]
    per = D[:,0]
    disper = D[:,1]

    command = [ftan, " ", filename, ' fmin=', fmin, ' fmax=', fmax,
               ' vgMin=', vgmin, ' vgMax=', vgmax, ' bmin=', bmin, ' bmax=',
               bmax,
               ' disp=all out=mat diag=', diagramtype, ' nfreq=', nfreq,
               ' ampMin=', ampmin]

    com = "".join([str(a) for a in command])
    print(com)
    os.system(com)
    seeds = np.loadtxt('write_disp.txt')
    isort = np.argsort(seeds[:, 0])
    seeds = seeds[isort]

    return per, disper, seeds


if __name__ == "__main__":
    filename = r"DK_NRS_DK_NUUG_Sym.SAC"
    fmin = 0.0066667
    fmax = 0.33333
    vgmin = 2.5
    vgmax = 5.0
    bmin = 0.0022
    bmax = 0.025
    diagramtype = 'PV'
    nfreq = 40
    ampmin = '0.05'
    dist = 1.2028e3

    pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, dist)