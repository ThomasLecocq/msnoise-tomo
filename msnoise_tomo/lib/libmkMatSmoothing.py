import os
from future.utils import native_str
import ctypes
import numpy as np

libdir = os.path.join(os.path.dirname(__file__), "mkMatSmoothing.cp36-win_amd64.pyd")
libsmooth = ctypes.CDLL(str(libdir))

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libsmooth.main.argtypes = [ctypes.c_int, # argc
                           LP_LP_c_char] # argv]
libsmooth.main.restype = ctypes.c_int
