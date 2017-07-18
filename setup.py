
import sys
import os
from setuptools import setup, find_packages

def pre_install_build():
    print("ok, trying to build all that C/C++ codes")
    for module in ["configparser.cpp",
                   "fft_NR.cpp",
                   "fta_param.cpp",
                   "libfta.cpp",
                   "readsac.cpp"]:
        os.system("g++ -c src/%s" % module)

    if sys.platform[:3] == "win":
        os.system("gcc -m32 src/mk_MatSmoothing.c -lm -o msnoise_tomo/lib/mk_MatSmoothing.exe")
        os.system("gcc -m32 src/mk_MatPaths.c -lm -o msnoise_tomo/lib/mk_MatPaths.exe")
        os.system("g++ src/vg_fta.cpp libfta.o readsac.o configparser.o fta_param.o fft_NR.o -lm -o msnoise_tomo/lib/ftan.exe")
    else:
        os.system(
            "gcc -m32 src/mk_MatSmoothing.c -lm -o msnoise_tomo/lib/mk_MatSmoothing")
        os.system(
            "gcc -m32 src/mk_MatPaths.c -lm -o msnoise_tomo/lib/mk_MatPaths")
        os.system(
            "g++ src/vg_fta.cpp libfta.o readsac.o configparser.o fta_param.o fft_NR.o -lm -o msnoise_tomo/lib/ftan")
        # os.system("make -C ./src/ all")
        # os.system("mv ./src/mk_MatPaths ./src/mk_MatSmoothing ./msnoise_tomo/lib/")
        # os.system("g++ src/vg_fta.cpp libfta.o readsac.o configparser.o fta_param.o fft_NR.o -lm -o msnoise_tomo/lib/ftan")

pre_install_build()

setup(
    name='msnoise_tomo',
    version='0.1a',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['msnoise',],
    entry_points = {
        'msnoise.plugins.commands': [
            'tomo = msnoise_tomo.plugin_definition:tomo',
        ],
        'msnoise.plugins.jobtypes': [
            'register = msnoise_tomo.plugin_definition:register_job_types',
        ],
        'msnoise.plugins.table_def': [
            'TomoConfig = msnoise_tomo.tomo_table_def:TomoConfig',
        ],
        'msnoise.plugins.admin_view': [
            'TomoConfigView = msnoise_tomo.plugin_definition:TomoConfigView',
        ],

    },
    author = "Thomas Lecocq & MSNoise dev team",
    author_email = "Thomas.Lecocq@seismology.be",
    description = "A Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise",
    license = "SISPROBE LICENCE",
    url = "http://www.msnoise.org",
    keywords="noise monitoring seismic velocity change dvv dtt doublet stretching cross-correlation acoustics seismology"
)
