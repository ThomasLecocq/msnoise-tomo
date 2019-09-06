import os
from future.utils import native_str
import ctypes
import numpy as np
from obspy.core.util.libnames import cleanse_pymodule_filename, _get_lib_name
import platform
lib = "mk_MatPaths"

libname = _get_lib_name(lib, add_extension_suffix=True)
libname = os.path.join(os.path.dirname(__file__), libname)
print("Path lib name:", libname)

libpath = ctypes.CDLL(str(libname))

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libpath.main.argtypes = [ctypes.c_int, # argc
                           LP_LP_c_char] # argv]
libpath.main.restype = ctypes.c_int

def path(file, gridfile):
    args = ["placeholder",
            file,
            gridfile]
    argc = len(args)
    argv = (LP_c_char * (argc + 1))()
    for i, arg in enumerate(args):
        enc_arg = arg.encode('utf-8')
        argv[i] = ctypes.create_string_buffer(enc_arg)
        
    return libpath.main(argc, argv)