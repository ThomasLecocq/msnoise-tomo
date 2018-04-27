import os
from future.utils import native_str
import ctypes
import numpy as np

libdir = os.path.join(os.path.dirname(__file__), "vg_fta.cpython-36m-x86_64-linux-gnu.so")
libfta = ctypes.CDLL(str(libdir))

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libfta.main.argtypes = [ctypes.c_int, # argc
                        LP_LP_c_char] # argv]
libfta.main.restype = ctypes.c_int


def ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, dist, disp):
    args = ["placeholder",
            filename,
            'fmin=%f' % float(fmin),
            'fmax=%f' % float(fmax),
            'vgMin=%f' % float(vgmin),
            'vgMax=%f' % float(vgmax),
            'bmin=%f' % float(bmin),
            'bmax=%f' % float(bmax),
            'disp=%s' % disp,
            'out=mat',
            'diag=%s' % diagramtype,
            'nfreq=%i' % int(nfreq),
            'ampMin=%f' % float(ampmin)
            ]

    argc = len(args)
    argv = (LP_c_char * (argc + 1))()
    for i, arg in enumerate(args):
        enc_arg = arg.encode('utf-8')
        argv[i] = ctypes.create_string_buffer(enc_arg)

    return libfta.main(argc, argv)