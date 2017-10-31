.. _installation:


Requirements
------------
MSNoise-TOMO is a plugin for MSNoise. You should therefore install MSNoise
beforehand: MSNoise_.

.. warning::
    Contrary to MSNoise, MSNoise-TOMO is ONLY compatible with Python 3.6.

MSNoise-TOMO uses C and C++ sources, which means you have to be able to compile
those source codes on your machine.

Linux
~~~~~

.. code-block:: sh

    sudo apt-get install build-essential

Windows
~~~~~~~

#. Install Microsoft Visual Studio Build Tools 2017 (VS_ bottom of
   the page).
#. Install TDM-GCC 64bits (TDM_) and make sure to add the install folder in
   the PATH.

To test the installations, in a new cmd.exe box, you should be able to run gcc:

.. code-block:: sh

    C:\Users\TLecocq>gcc

    gcc: fatal error: no input files
    compilation terminated.

Extra dependencies
~~~~~~~~~~~~~~~~~~

MSNoise-TOMO has some extra requirements:

- shapely
- pyproj

which can be installed using conda (recommended) or pip.



Installing MSNoise-TOMO
-----------------------

The package should as easy to install as:

.. code-block:: sh

    pip install msnoise-tomo.zip








.. _MSNoise: http://msnoise.org/doc/installation.html
.. _VS: https://www.visualstudio.com/downloads/
.. _TDM: http://tdm-gcc.tdragon.net/download
