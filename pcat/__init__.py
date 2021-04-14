# plotting
import matplotlib as mpl

import matplotlib.pyplot as plt
import seaborn as sns

# numpy
import random as randommod
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons

# scipy
import scipy as sp
import scipy.interpolate
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
import scipy.ndimage

import scipy.fftpack

import scipy.ndimage.filters
import scipy.sparse

# jit
from numba import jit
import threading

#import pyximport; pyximport.install()
import ctypes

import subprocess as subp, psutil

import astropy
import astropy as ap
from astropy.convolution import convolve_fft, AiryDisk2DKernel

# multiprocessing
import multiprocessing as mp

from copy import deepcopy

# FITS files
import h5py

# utilities
import os, time, sys, getpass, glob, fnmatch, inspect, traceback, shelve
import pickle as cPickle
import functools

# tdpy
import tdpy.util
from tdpy.util import summgene
import tdpy.mcmc

# HealPix
#import healpy as hp
#from healpy import ang2pix

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

# secondaries
## Symbolic Jacobian calculation
#import sympy

## Probabilistic Graphical Model generation
import networkx as nx


