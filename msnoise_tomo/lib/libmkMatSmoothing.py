import os
from future.utils import native_str
import ctypes
import numpy as np
from obspy.core.util.libnames import cleanse_pymodule_filename, _get_lib_name
import platform
lib = "mkMatSmoothing"

libname = _get_lib_name(lib, add_extension_suffix=True)
libname = os.path.join(os.path.dirname(__file__), libname)
print(libname)

libsmooth = ctypes.CDLL(str(libname))

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libsmooth.main.argtypes = [ctypes.c_int, # argc
                           LP_LP_c_char] # argv]
libsmooth.main.restype = ctypes.c_int
