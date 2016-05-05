
# coding: utf-8

# In[ ]:

# plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
mpl.rc('image', interpolation='none', origin='lower')

import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

# numpy
import numpy as np
from numpy import *
from numpy.random import *
from numpy.random import choice

# scipy
import scipy as sp
from scipy import ndimage
from scipy.interpolate import *
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
from scipy import ndimage

# multiprocessing
import multiprocessing as mp

# healpy
import healpy as hp
from healpy.rotator import angdist
from healpy import ang2pix

# pyfits
import pyfits as pf

# utilities
import os, time, sys, datetime, warnings, getpass, glob, fnmatch, cPickle
import functools

# tdpy
import tdpy.util


# In[ ]:

def retr_psfn(globdata, psfipara, indxenertemp, thisangl, psfntype):

    if psfntype == 'singgaus':
        numbformpara = 1
    elif psfntype == 'singking':
        numbformpara = 2
    elif psfntype == 'doubgaus':
        numbformpara = 3
    elif psfntype == 'gausking':
        numbformpara = 4
    elif psfntype == 'doubking':
        numbformpara = 5
    
    thisangltemp = thisangl[None, :, None]

    indxpsfiparatemp = numbformpara * (indxenertemp[:, None] + globdata.numbener * globdata.indxevtt[None, :])
    
    if psfntype == 'singgaus':
        sigc = psfipara[indxpsfiparatemp]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(thisangltemp, sigc)

    elif psfntype == 'singking':
        sigc = psfipara[indxpsfiparatemp]
        gamc = psfipara[indxpsfiparatemp+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(thisangltemp, sigc, gamc)
        
    elif psfntype == 'doubgaus':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(thisangltemp, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        gamt = psfipara[indxpsfiparatemp+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(thisangltemp, frac, sigc, sigt, gamt)
        
    elif psfntype == 'doubking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        gamc = psfipara[indxpsfiparatemp+2]
        sigt = psfipara[indxpsfiparatemp+3]
        gamt = psfipara[indxpsfiparatemp+4]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_doubking(thisangltemp, frac, sigc, gamc, sigt, gamt)
            
    return psfn

