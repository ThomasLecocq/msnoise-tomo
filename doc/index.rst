MSNoise-TOMO - Documentation
============================

This documentation covers the functionality of the MSNoise-TOMO package.

As for MSNoise, the goal of this extension "suite" is to provide researchers
with an efficient processing tool, while keeping the need for coding to a
minimum and avoiding being a black box. Moreover, as long as the in- and outputs
of each step are respected, they can easily be replaced with one's own codes !

The MSNoise-TOMO Plugin extends the MSNoise workflow right after the
`stack --ref` step.

If you use MSNoise for your research and prepare publications, please consider
citing MSNoise:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.

and

**Lecocq, T., Mordret, A. (in prep)**, MSNoise-TOMO.


Commercial Usage
================

If you plan to use MSNoise-TOMO for commercial purposes, please contact Thomas
Lecocq directly.

Installation / Adding the plugin to a project
=============================================

#. Install the package and requirements (see MSNoise_) 
#. In the current project folder, add msnoise_tomo to the plugins:
   e.g. ``msnoise config set plugins=msnoise_tomo``
#. run ``msnoise p tomo install``
#. if all goes well, the following command should work ``msnoise p tomo info``


Workflow
========

#. Reset the stack jobs and redo the REF stacks: ``msnoise reset STACK --all``
   and ``msnoise stack -r``
#. ``msnoise info -j`` should show ``TOMO_SAC`` jobs "T"o do.
#. Run ``msnoise p tomo prepare_ccf`` to create the SAC files necessary for the
   next step
#. ``msnoise info -j`` should show ``TOMO_FTAN`` jobs "T"o do.
#. ``msnoise p tomo iftan`` starts a GUI that allows you to check/save the 
   dispersion curves for individual sac files
#. Once satisfactory parameters are defined, add them to the database in the 
   ``tomo-config`` table
#. Run the ``msnoise p tomo ftan`` to compute the dispersion curves
   automatically for all your files (currently the same parameters are used for
   all components, filters, and distances).
#. ``msnoise p tomo prepare_tomo`` will create the input files for the period-
   map inversion procedure, for each of ``ftan_periods`` configured in the DB.
#. ``msnoise p tomo answt`` will compute a period map for each configured. 
   ``ftan_periods``. This step will output figures and a KMZ file to be opened
   in Google Earth. To compute only one period, pass the
   ``msnoise p tomo answt --per 5.0`` (provided 5.0 has been set in
   ``ftan_periods``).


Development & Miscellaneous
===========================

.. toctree::
    :maxdepth: 2

    contributors


.. toctree::
    :maxdepth: 2

    clickhelp/msnoise





Release Notes
=============



.. _PDF: http://msnoise.org/doc/MSNoise.pdf
.. _MSNoise: http://msnoise.org/doc/installation.html