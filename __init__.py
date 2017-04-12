# plotting
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

# numpy
import random as randommod
import numpy as np
from numpy import *
from numpy.random import *

# scipy
import scipy as sp
from scipy.interpolate import *
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
import scipy.ndimage
import scipy.ndimage.filters

# jit
from numba import jit
import threading

#import pyximport; pyximport.install()
import ctypes

import astropy as ap
from astropy.convolution import convolve_fft, AiryDisk2DKernel

# multiprocessing
import multiprocessing as mp

from copy import deepcopy

# FITS files
import pyfits as pf, h5py

# utilities
import os, time, sys, getpass, glob, fnmatch, cPickle, inspect, traceback, shelve
import functools

# tdpy
import tdpy.util
from tdpy.util import show, summgene, summ
import tdpy.mcmc

# HealPix
import healpy as hp
from healpy import ang2pix

# ignore warnings if not in diagnostic mode
import warnings
    
#from skimage.feature import blob_doh

# defualt options
#seterr(divide='raise', over='raise', invalid='raise')
#seterr(all='raise')
#seterr(under='ignore')
warnings.simplefilter('ignore')
np.set_printoptions(linewidth=180)
sns.set(context='poster', style='ticks', color_codes=True)
mpl.rc('image', interpolation='nearest', origin='lower')

# secondaries
## Symbolic Jacobian calculation
import sympy

## Probabilistic Graphical Model generation
import networkx as nx


