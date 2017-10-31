import os
from future.utils import native_str
import ctypes
import numpy as np

libdir = os.path.join(os.path.dirname(__file__), "vg_fta.cp36-win_amd64.pyd")
libfta = ctypes.CDLL(str(libdir))

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libfta.main.argtypes = [ctypes.c_int, # argc
                        LP_LP_c_char] # argv]
libfta.main.restype = ctypes.c_int


def ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                      diagramtype, nfreq, ampmin, dist):
    args = ["placeholder",
            filename,
            'fmin=%f' % fmin,
            'fmax=%f' % fmax,
            'vgMin=%f' % vgmin,
            'vgMax=%f' % vgmax,
            'bmin=%f' % bmin,
            'bmax=%f' % bmax,
            'disp=none',
            'out=mat',
            'diag=%s' % diagramtype,
            'nfreq=%i' % nfreq,
            'ampMin=%f' % ampmin
            ]

    argc = len(args)
    argv = (LP_c_char * (argc + 1))()
    for i, arg in enumerate(args):
        enc_arg = arg.encode('utf-8')
        argv[i] = ctypes.create_string_buffer(enc_arg)

    return libfta.main(argc, argv)