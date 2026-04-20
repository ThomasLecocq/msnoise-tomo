import os
import sys
import time

import numpy as np
from .lib.libvg_fta import ftan

def pickgroupdispcurv(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, dist, pinit=0, vinit=0):
    # if sys.platform[:3] == "win":
    #     ftan = os.path.join(os.path.split(os.path.realpath(__file__))[0],"lib", r"ftan.exe")
    # else:
    #     ftan = os.path.join(os.path.split(os.path.realpath(__file__))[0],"lib", r"ftan")
    
    # command = [ftan," ", filename, ' fmin=', fmin, ' fmax=', fmax,
    #     ' vgMin=', vgmin, ' vgMax=', vgmax, ' bmin=', bmin ,' bmax=', bmax,
    #     ' disp=none out=mat diag=', diagramtype, ' nfreq=', nfreq, ' ampMin=', ampmin]
    if pinit == 0 and vinit == 0:
        
        print("\nRunning FTAN the first time: %s\n"%filename)

		# run the C++ code on the SAC file 
        ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
             diagramtype, nfreq, ampmin, dist, disp="none")
        time.sleep(0.1)
    
    	# Grab the output from the C++ code: this is the FTAN amplitude matrix and axes
        amp = np.loadtxt('write_amp.txt') # this the FTAN matrix
        amp = amp.T # tranpose
    	
    	# Find the indices (velocity and period) of the maximum in the FTAN matrix
        index = np.unravel_index(np.argmax(amp), amp.shape)
        if len(index) == 2:
            iv, ip = index
        else:
            iv, ip = index[0]

        # Read the files properly depending on the diagramtype and get the SEED values
        if diagramtype == 'PV':
            P = np.loadtxt('write_FP.txt') # this is the period axis
            V = np.loadtxt('write_TV.txt') # this is the velocity axis
            print("Seed velocity: %f [km/s]"%V[iv])
            print("Seed period: %f [s]"%P[ip])  
            vinit = V[iv]
            pinit = P[ip]
        elif diagramtype == 'FV':
            F = np.loadtxt('write_FP.txt') # this is the frequency axis
            V = np.loadtxt('write_TV.txt') # this is the velocity axis
            print("Seed velocity: %f [km/s]"%V[iv])
            print("Seed frequency: %f [Hz]"%F[ip])    
            vinit = V[iv]
            pinit = 1/F[ip]
        elif diagramtype == 'FT':
            F = np.loadtxt('write_FP.txt') # this is the frequency axis
            T = np.loadtxt('write_TV.txt') # this is the time axis
            print("Seed time: %f"%T[iv])
            print("Seed frequency: %f"%F[ip])   
            vinit = dist/T[iv]
            pinit = 1/F[ip]
        elif diagramtype == 'PT':
            P = np.loadtxt('write_FP.txt') # this is the period axis
            T = np.loadtxt('write_TV.txt') # this is the timeaxis
            vinit = dist/T[iv]
            pinit = P[ip]
            print("Seed time: %f"%vinit)
            print("Seed period: %f"%pinit)	


    # # Set the starting points depending on which plot type
    # if diagramtype == 'PV':
    #     finit = "%.5f"% (1./pinit)
    # 	# finit = pinit
    #     vginit = vinit

    # elif diagramtype == 'FV':
    #     finit = pinit
    #     # finit = "%.5f"% (1./pinit)
    #     vginit = vinit

    # elif diagramtype == 'FT':
    #     finit = pinit
    #     vginit = dist/vinit

    # elif diagramtype == 'PT':
    #     finit = 1. / pinit
    #     vginit = dist/vinit

    # command = [ftan," ", filename, ' fmin=', fmin, ' fmax=', fmax,
    #     ' vgMin=', vgmin, ' vgMax=', vgmax, ' finit=', finit, ' vginit=',
    #            vginit, ' bmin=', bmin ,' bmax=', bmax,
    #     ' disp=cont out=mat diag=', diagramtype, ' nfreq=', nfreq,
    #            ' ampMin=', ampmin]
    # 
    # com = "".join([str(a) for a in command])
    # print(com)
    # os.system(com)
    # print("vinit: %f"%vinit)

    print("\nRunning FTAN a second time with the SEED\n")

    # Apply the FTAN C++ code now with the seeds pinit, vinit
    ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
         diagramtype, nfreq, ampmin, dist, disp="cont", tinit=pinit, vginit=vinit)
    time.sleep(0.1)
    # TODO still need to figure out how to use the tginit input in the C++ code. vinit doesn't seem to have an impact!

   	# # Process the automatically picked dispersion curve from the C++ code.
    # D = np.loadtxt('write_disp.txt')
    # if D.ndim == 2: # make sure that there is more than one pick
    #     isort  = np.argsort(D[:,0]) # sort based on the first column (period)
    #     D      = D[isort]
    #     per    = D[:,0]
    #     disper = D[:,1]
    # else:
    #     print("Only one dispersion pick...check data!!!")
    #     per = D[0]
    #     disper = D[1]

    # print(per)
    # print(disper)
    # print("Periods: %s"%per)
	# print("Velocties: %s"%disper)

    # command = [ftan, " ", filename, ' fmin=', fmin, ' fmax=', fmax,
    #            ' vgMin=', vgmin, ' vgMax=', vgmax, ' bmin=', bmin, ' bmax=',
    #            bmax,
    #            ' disp=all out=mat diag=', diagramtype, ' nfreq=', nfreq,
    #            ' ampMin=', ampmin]
    # 
    # com = "".join([str(a) for a in command])
    # print(com)
    # os.system(com)
    
    # ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
    #      diagramtype, nfreq, ampmin, dist, disp="cont", tinit=tinit, vginit=vginit)
    # time.sleep(0.1)
    
    # seeds = np.loadtxt('write_disp.txt')
    # isort = np.argsort(seeds[:, 0])
    # seeds = seeds[isort] # just sort the whole matrix based on period

    # return per, disper, seeds

    # print(per)
    # print(disper)

    # return per, disper, D
    return


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