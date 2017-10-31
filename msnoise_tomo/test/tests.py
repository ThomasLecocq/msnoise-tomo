import unittest
import traceback
import os
import ctypes
import datetime
import shutil
import glob
from obspy import read

class MSNoiseTomoTests(unittest.TestCase):

    def setUp(self):
        path = os.path.abspath(os.path.dirname(__file__))
        self.data_folder = os.path.join(path, 'data')

    def test_mksmooth(self):
        from msnoise_tomo.lib.libmkMatSmoothing import libsmooth, LP_c_char
        args = ["TEST", "1",
                os.path.join(self.data_folder, "GLISNGrid.txt")]

        argc = len(args)
        argv = (LP_c_char * (argc + 1))()
        for i, arg in enumerate(args):
            enc_arg = arg.encode('utf-8')
            argv[i] = ctypes.create_string_buffer(enc_arg)

        print("exit:", libsmooth.main(argc, argv))

    def test_ftan(self):
        from msnoise_tomo.lib.libvg_fta import ftan
        import numpy as np
        filename = os.path.join(self.data_folder, "DK_NRS_DK_NUUG_Sym.SAC")
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
        print("exit:", ftan(filename, fmin, fmax, vgmin, vgmax, bmin, bmax,
                            diagramtype, nfreq, ampmin, 0))

        V = np.loadtxt('write_TV.txt')
        P = np.loadtxt('write_FP.txt')
        amp = np.loadtxt('write_amp.txt').T
        print(amp.shape)


def main():
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")
    import os
    import sys


    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MSNoiseTomoTests)
    runner = unittest.TextTestRunner(verbosity=4)
    result = runner.run(suite)
    if not result.wasSuccessful():
        sys.exit(1)

if __name__ == '__main__':
    main()
