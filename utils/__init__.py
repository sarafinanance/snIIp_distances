from utils.format_data import make_Table_target
from utils.format_data import import_template
from utils.format_data import format_target
from utils.format_data import format_template
from utils.format_data import make_spec1d_target
from utils.format_data import make_spec1d_temps 

# Scipy
import scipy 
from scipy.optimize import curve_fit
from scipy import stats
from scipy.signal.windows import tukey
from scipy.interpolate import RectBivariateSpline, interp1d, LSQUnivariateSpline, UnivariateSpline, splint, BSpline, splrep, CubicSpline, splev
import scipy.integrate as integrate

# Numerical analysis
import numpy as np
import emcee
import random
import statistics
import itertools
import pandas as pd
import copy 
import statistics
from numpy import polyfit