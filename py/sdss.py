
# coding: utf-8

# In[1]:

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

# pyfits
import pyfits as pf

# utilities
import os, time, sys, datetime, warnings, getpass, glob

# tdpy
import tdpy_util


# In[2]:

def writ_sdss():
    
    pixlsize = 0.396 # [arcsec]
    npixlside = 100

    maxmgang = npixlside * pixlsize / 3600. / 2.
    npixlheal = npixlside**2 * 12
    apix = deg2rad(pixlsize / 3600.)**2
    enerstrg = ['i', 'r', 'g']
    nener = len(enerstrg)
    iener = arange(nener)
    nevtt = 1

    from astropy.coordinates import ICRS, Galactic
    from astropy import units as u

    sdssdataflux = zeros((nener, npixlside, npixlside, nevtt))
    for i in iener:
        path = os.environ["PNTS_TRAN_DATA_PATH"] + '/frame-' + enerstrg[i] + '-001458-4-0700.fits'
        data, hdr = pf.getdata(path, 0, header=True)
        rasccntr = hdr['CRVAL1']
        declcntr = hdr['CRVAL2']

        rascbndr = zeros(4)
        declbndr = zeros(4)
        
        rascbndr[0] = rasccntr - hdr['CRPIX1'] * hdr['CD2_1'] - hdr['CRPIX2'] * hdr['CD2_2']
        declbndr[0] = declcntr - hdr['CRPIX1'] * hdr['CD1_1'] - hdr['CRPIX2'] * hdr['CD1_2']
        
        rascbndr[1] = rasccntr - hdr['CRPIX1'] * hdr['CD2_1'] - hdr['CRPIX2'] * hdr['CD2_2']
        declbndr[1] = declcntr + hdr['CRPIX1'] * hdr['CD1_1'] - hdr['CRPIX2'] * hdr['CD1_2']
        
        rascbndr[1] = rasccntr - hdr['CRPIX1'] * hdr['CD2_1'] - hdr['CRPIX2'] * hdr['CD2_2']
        declbndr[1] = declcntr - hdr['CRPIX1'] * hdr['CD1_1'] - hdr['CRPIX2'] * hdr['CD1_2']
        
        rascbndr[0] = rasccntr - hdr['CRPIX1'] * hdr['CD2_1'] - hdr['CRPIX2'] * hdr['CD2_2']
        declbndr[0] = declcntr - hdr['CRPIX1'] * hdr['CD1_1'] - hdr['CRPIX2'] * hdr['CD1_2']
        
    
        #CRPIX1  = 1.02450000000000E+03 / Column Pixel Coordinate of Reference
        #CRPIX2  = 7.44500000000000E+02 / Row Pixel Coordinate of Reference Pix
        #CRVAL1  = 6.42142612400000E+01 / DEC at Reference Pixel
        #CRVAL2  = 2.51207342810000E+02 / RA at Reference Pixel

        #CD1_1   = 4.75626416015645E-05 / DEC degrees per column pixel
        #CD1_2   = -9.9116868279565E-05 / DEC degrees per row pixel
        #CD2_1   = 9.91802898939385E-05 / RA  degrees per column pixel

        print 'rascbndr'
        print rascbndr
        print 'declbndr'
        print declbndr
        

        print 'data.shape'
        print data.shape
        
        calb, hdr = pf.getdata(path, 1, header=True)
        back, hdr = pf.getdata(path, 2, header=True)

        xaxi = arange(back['ALLSKY'].shape[1], dtype=float)
        yaxi = arange(back['ALLSKY'].shape[2], dtype=float)
        
        back = interp2d(xaxi, yaxi, back['ALLSKY'])(back['XINTERP'].flatten(), back['YINTERP'].flatten())

        data /= calb[None, :]
        data += back
        
        print 'amin(data)'
        print amin(data)
        print 'amax(data)'
        print amax(data)
        
        print 'amin(back)'
        print amin(back)
        print 'amax(back)'
        print amax(back)
        
        print
        
        data /= apix

        print 'Frame RA: ', rasccntr
        print 'Frame DEC: ', declcntr
        
        #rasccntr = rasccntr + rascpixlsize * (data.shape[0] - npixlside) / 2. / 3600.
        #declcntr = declcntr - declpixlsize * (data.shape[1] - npixlside) / 2. / 3600.
        
        objt = ICRS(ra=rasccntr, dec=declcntr, unit=(u.degree, u.degree))
        lgalcntr = objt.galactic.l.degree
        bgalcntr = objt.galactic.b.degree

        print 'Patch RA: ', rasccntr
        print 'Patch DEC: ', declcntr
        print 'Patch l: ', lgalcntr
        print 'Patch b: ', bgalcntr
        
        sdssdataflux[i, :, :, 0] = data[-npixlside:, -npixlside:]
        
    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/sdssdataflux.fits' 
    pf.writeto(path, sdssdataflux, clobber=True)

    sdssbackflux = ones((nener, npixlside, npixlside, nevtt)) / apix
    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/sdssbackflux.fits' 
    pf.writeto(path, sdssbackflux, clobber=True)
    
    sdssexpo = ones((nener, npixlside, npixlside, nevtt))
    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/sdssexpo.fits' 
    pf.writeto(path, sdssexpo, clobber=True)

    
#writ_sdss()


# In[ ]:



