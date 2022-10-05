#!/usr/bin/env python
# coding: utf-8

#libraries




# os
import glob
import os

# Spectra
import specutils
from specutils.analysis import correlation
from specutils.fitting import continuum 
from specutils.analysis import template_comparison
from specutils.fitting.continuum import fit_generic_continuum
from specutils.manipulation import LinearInterpolatedResampler
from specutils.manipulation import FluxConservingResampler

# sncosmo
import sncosmo

# Astropy imports
import astropy
from astropy.modeling import models, fitting
from astropy.modeling.polynomial import Chebyshev1D
from astropy.table import QTable
from astropy.time import Time
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import quantity_support
import astropy.io.fits as pyfits
from astropy.nddata import StdDevUncertainty
from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity

# Plotting 
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns
from matplotlib.gridspec import GridSpec

#package imports
from analysis_packages import correlation
from analysis_packages.correlation import cross_correlation
from analysis_packages.distance_modulus import calculate_distance
from analysis_packages.explosion_date import get_explosion_date
from analysis_packages.k_corrections import kcorrections
import utils



# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.size"] = "12"
# quantity_support()  # for getting units on the axes below  
# get_ipython().run_line_magic('matplotlib', 'inline')

# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.size"] = "15"
# plt.rcParams['mathtext.fontset'] = 'stix'
# plt.rcParams['font.family'] = 'STIXGeneral'


# In[37]:


# set paths
ZTF_name ='ZTF19abqhobb'
z_input = 0.01815 #NEEDS TO BE SET EACH TIME


def main(ZTF_name, z_input):
    #class to load data 
    #get correlation coeff 
    return 



if __name__ == "__main__":
    ZTF_name = input('Please provide target name: ')
    z_input = input('Please provide target redshift: ')
    main(ZTF_name, z_input).run()