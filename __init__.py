# plotting
import matplotlib as mpl
mpl.use('Agg')
mpl.rc('image', interpolation='nearest', origin='lower')
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

# numpy
import random as randommod
import numpy as np
from numpy import *
from numpy.random import *
#from numpy.random import choice
# temp
#seterr(divide='raise', over='raise', invalid='raise')
#seterr(all='raise')
#seterr(under='ignore')

# scipy
import scipy as sp
from scipy.interpolate import *
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
import scipy.ndimage
import scipy.ndimage.filters

import astropy as ap
    
import shelve

# multiprocessing
import multiprocessing as mp

# pyfits
import pyfits as pf

# utilities
import os, time, sys, getpass, glob, fnmatch, cPickle, inspect
import functools

# Symbolic Jacobian calculation
import sympy

# tdpy
import tdpy.util
from tdpy.util import show, summgene, summ
import tdpy.mcmc

import networkx as nx

np.set_printoptions(linewidth=180)

# healpy
import healpy as hp
#from healpy.rotator import angdist
from healpy import ang2pix

# ignore warnings if not in diagnostic mode
import warnings
warnings.simplefilter('ignore')
    
import franlens

from astropy.convolution import convolve, AiryDisk2DKernel

from copy import deepcopy

