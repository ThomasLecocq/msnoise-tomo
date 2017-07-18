import ctypes
import matplotlib.pyplot as plt
import numpy as np
from msnoise_tomo.lib.libmkMatSmoothing import libsmooth, LP_c_char
from msnoise_tomo.lib.libvg_fta import ftan

args = ["TEST", "1", r"D:\PythonForSource\MSNoise_Stack\msnoise-tomo\msnoise_tomo\DATAFiles\GLISNGrid.txt"]

argc = len(args)
argv = (LP_c_char * (argc + 1))()
for i, arg in enumerate(args):
    enc_arg = arg.encode('utf-8')
    argv[i] = ctypes.create_string_buffer(enc_arg)

libsmooth.main(argc, argv)



filename = r"D:\PythonForSource\MSNoise_Stack\msnoise-tomo\msnoise_tomo\DATAFiles\DK_NRS_DK_NUUG_Sym.SAC"
fmin = 0.0066667
fmax = 0.33333
vgmin = 2.5
vgmax = 5.0
bmin = 0.0022
bmax = 0.025
diagramtype = 'PV'
nfreq = 40
ampmin = 0.05
dist = 1.2028e3
ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, 0)

V = np.loadtxt('write_TV.txt')
P = np.loadtxt('write_FP.txt')
amp=np.loadtxt('write_amp.txt').T



plt.imshow(amp, aspect='auto', cmap="viridis")
plt.show()