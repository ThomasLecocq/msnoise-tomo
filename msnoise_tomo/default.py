from collections import OrderedDict
default = OrderedDict()

default['ftan_fmin'] = ["Minimum frequency for the FTAN matrix computation [Hz]",'0.1']
default['ftan_fmax'] = ["Maximum frequency for the FTAN matrix computation [Hz]",'1.0']
default['ftan_vgmin'] = ["Minimum group velocity for the FTAN matrix computation [km/s]",'0.5']
default['ftan_vgmax'] = ["Maximum group velocity for the FTAN matrix computation [km/s]",'5.0']
default['ftan_minWL'] = ["Minimum wavelength",'1.0']
default['ftan_bmin'] = ["",'0.0022']
default['ftan_bmax'] = ["",'0.025']
default['ftan_diagramtype'] = ["Choose the type of FTAN diagram to display -- options: PV, FV, PT, FT (F=frequency, P=period, V=velocity, T=time)",'PV']
default['ftan_periods'] = ["Choose the periods you want to store the dispersion picks for",'3.5, 5, 10, 20, 30, 50, 60']
default['ftan_nfreq'] = ["Number of frequencies beteween fmin,fmax for the FTAN matrix computation",'40']
default['ftan_ampmin'] = [" ",'0.05']

default["xstep"] = ["Step of the inversion grid in X direction [deg]","0.05"]
default["ystep"] = ["Step of the inversion grid in Y direction [deg]","0.05"]

default['alpha1']  = ["alpha1",'500']
default['beta1']   = ["beta1",'50']
default['lambda1'] = ["lambda1",'0.01']
default['sigma1']  = ["sigma1",'3']

default['alpha2']  = ["alpha2",'200']
default['beta2']   = ["beta2",'50']
default['lambda2'] = ["lambda2",'0.01']
default['sigma2']  = ["sigma2",'3']

default['v_cmap']  = ["Velocity Colormap, see <a href='https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html'>Matplotlib Colormaps</a> for choosing one", 'viridis']
default['d_cmap']  = ["Density Colormap, see <a href='https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html'>Matplotlib Colormaps</a> for choosing one", 'inferno']






