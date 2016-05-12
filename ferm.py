# numpy
import numpy as np
from numpy import *
from numpy.random import *
from numpy.random import choice

import matplotlib.pyplot as plt

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
import os, time, sys, datetime, warnings, getpass, glob, inspect

# tdpy
import tdpy

import sympy

def retr_axes():

    reco = 8
    evtc = 128
    
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[0:-1])
    diffener = binsener[1:] - binsener[0:-1]
    numbener = meanener.size
    indxener = arange(numbener)
    minmener = amin(binsener)
    maxmener = amax(binsener)
    
    global indxevtt
    evtt = array([4, 16, 32, 64])
    numbevtt = evtt.size
    indxevtt = arange(numbevtt)
    
    nside = 256
    numbpixl = nside**2 * 12
    apix = 4. * pi / numbpixl
    
    return reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix


def make_maps():
    
    cmnd = 'mkdir -p $FERMI_DATA/exposure/pcat'
    os.system(cmnd)
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    global strgregi
    strgregi = ' ra=INDEF dec=INDEF rad=INDEF '
    
    
    # make file lists
    cmnd = 'rm $PCAT_DATA_PATH/phot_pass8.txt'
    os.system(cmnd)
    cmnd = 'ls -d -1 $FERMI_DATA/weekly/photon/*.fits >> $PCAT_DATA_PATH/phot_pass8.txt'
    os.system(cmnd)

    weekinit = 9
    weekfinl = 218
    listtimefrac = array([1., 0.75, 0.5, 0.25])
    numbtime = listtimefrac.size
    for t, timefrac in enumerate(listtimefrac):
        numbweek = (weekfinl - weekinit) * timefrac
        listweek = floor(linspace(weekinit, weekfinl - 1, numbweek)).astype(int)
        cmnd = 'rm $PCAT_DATA_PATH/phot_pass7_time%d.txt' % t
        os.system(cmnd)
        for week in listweek:
            cmnd = 'ls -d -1 $FERMI_DATA/weekly/p7v6c/*_w%03d_* >> $PCAT_DATA_PATH/phot_pass7_time%d.txt' % (week, t)
            os.system(cmnd)
   
    numbproc = numbtime + 1

    global rtag, reco, evtc, strgtime
    rtag = ['full'] + ['cmp%d' % t for t in range(numbtime)]
    reco = [8] + [7 for t in range(numbtime)] 
    evtc = [128] + [2 for t in range(numbtime)] 
    strgtime = ['tmin=INDEF tmax=INDEF'] + ['tmin=239155201 tmax=364953603' for t in range(numbtime)] 
        
    # process pool
    pool = mp.Pool(numbproc)

    # spawn the processes
    #pool.map(make_maps_sing, range(numbproc))
    make_maps_sing(0)
    pool.close()
    pool.join()


def make_maps_sing(indxprocwork):

    infl = '$PCAT_DATA_PATH/phot_pass%d_%s.txt' % (reco[indxprocwork], rtag[indxprocwork])
    spac = '$PCAT_DATA_PATH/spac_pass%d_%s.txt' % (reco[indxprocwork], rtag[indxprocwork])
        
    for m in indxevtt:
        
        if reco[indxprocwork] == 7:
            if m == 3:
                thisevtt = 1
            elif m == 2:
                thisevtt = 2
            else:
                continue
        if reco[indxprocwork] == 8:
            thisevtt = evtt[m]
                
        sele = '$PCAT_DATA_PATH/sele_pass%d_evtc%03d_evtt%03d_%s.fits' % (reco[indxprocwork], evtc[indxprocwork], thisevtt, rtag[indxprocwork])
        filt = '$PCAT_DATA_PATH/filt_pass%d_evtc%03d_evtt%03d_%s.fits' % (reco[indxprocwork], evtc[indxprocwork], thisevtt, rtag[indxprocwork])
        live = '$PCAT_DATA_PATH/live_pass%d_evtc%03d_evtt%03d_%s.fits' % (reco[indxprocwork], evtc[indxprocwork], thisevtt, rtag[indxprocwork])
        cnts = '$PCAT_DATA_PATH/cnts_pass%d_evtc%03d_evtt%03d_%s.fits' % (reco[indxprocwork], evtc[indxprocwork], thisevtt, rtag[indxprocwork])
        expo = '$PCAT_DATA_PATH/expo_pass%d_evtc%03d_evtt%03d_%s.fits' % (reco[indxprocwork], evtc[indxprocwork], thisevtt, rtag[indxprocwork])

        cmnd = 'gtselect infile=' + infl + ' outfile=' + sele + strgregi + \
            strgtime[indxprocwork] + ' emin=100 emax=100000 zmax=90 evclass=%d evtype=%d' % (evtc[indxprocwork], thisevtt)
        #os.system(cmnd)
        #print '%d, ' % indxprocwork + cmnd
        #print

        cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1"' + ' outfile=' + filt + ' roicut=no'
        os.system(cmnd)
        print '%d, ' % indxprocwork + cmnd
        print

        cmnd = 'gtbin evfile=' + filt + ' scfile=' + spac + ' outfile=' + cnts + \
            ' ebinalg=FILE ebinfile=/n/fink1/fermi/exposure/gcps_time/gtbndefn.fits algorithm=HEALPIX' + \
            ' hpx_ordering_scheme=RING coordsys=GAL hpx_order=8 hpx_ebin=yes'
        os.system(cmnd)
        print '%d, ' % indxprocwork + cmnd
        print

        cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live + ' dcostheta=0.025 binsz=1'
        os.system(cmnd)
        print '%d, ' % indxprocwork + cmnd
        print

        cmnd = 'gtexpcube2 infile=' + live + ' cmap=' + cnts + ' outfile=' + expo + ' irfs=CALDB evtype=%03d bincalc=CENTER' % thisevtt
        os.system(cmnd)
        print '%d, ' % indxprocwork + cmnd
        print
        print
        print

def writ_isot():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    # isotropic background
    path = os.environ["PCAT_DATA_PATH"] + '/iso_P8R2_ULTRACLEAN_V6_v06.txt'
    isotdata = loadtxt(path)
    enerisot = isotdata[:, 0] * 1e-3 # [GeV]
    isotflux = isotdata[:, 1] * 1e3 # [1/cm^2/s/sr/GeV]
    isotfluxheal = empty((numbener, numbpixl, numbevtt))
    nsampbins = 10
    enersamp = logspace(log10(amin(binsener)), log10(amax(binsener)), nsampbins * numbener)
    isotflux = interpolate.interp1d(enerisot, isotflux)(enersamp)
    for i in range(numbener):
        isotfluxheal[i, :, :] = trapz(isotflux[i*nsampbins:(i+1)*nsampbins],                              enersamp[i*nsampbins:(i+1)*nsampbins]) / diffener[i]
        
    path = os.environ["PCAT_DATA_PATH"] + '/fermisotflux.fits'
    pf.writeto(path, isotfluxheal, clobber=True)

    nfwpfluxtemp = tdpy.retr_nfwp(1., nside, norm=5.)
    nfwpspec = ones(numbener)
    nfwpfluxheal = zeros((numbener, numbpixl, numbevtt))
    for i in indxener:
        for m in indxevtt:
            nfwpfluxheal[i, :, m] = nfwpspec[i] * nfwpfluxtemp
    path4 = os.environ["PCAT_DATA_PATH"] + '/nfwpflux.fits'
    pf.writeto(path4, nfwpfluxheal, clobber=True)


def prep_maps():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    liststrgener = ['ENERGY1', 'ENERGY2', 'ENERGY3', 'ENERGY4', 'ENERGY5']
    liststrgchan = ['CHANNEL1', 'CHANNEL2', 'CHANNEL3', 'CHANNEL4', 'CHANNEL5']
    listdatatype = ['comp', 'full']
    listregitype = ['igal', 'ngal']
    
    cnts = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))
    flux = zeros((numbener, numbpixl, numbevtt))
    
    for k, datatype in enumerate(listdatatype):
    
        for m in indxevtt:
            if m < 2:
                continue
            path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_%s_evtt%d.fits' % (datatype, m)
            expoarry = pf.getdata(path, 1)
            for i in indxener:
                expo[i, :, m] = expoarry[liststrgener[i]]

            path = os.environ["PCAT_DATA_PATH"] + '/fermcnts_%s_evtt%d.fits' % (datatype, m)
            cntsarry = pf.getdata(path)
            for i in indxener:
                cnts[i, :, m] = cntsarry[liststrgchan[i]]

        indxexpo = where(expo > 0.) 
        flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
        flux /= diffener[:, None, None]

        for regitype in listregitype:
            if regitype == 'ngal':
                for i in indxener:
                    for m in indxevtt:
                        if m < 2:
                            continue
                        almc = hp.map2alm(flux[i, :, m])
                        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                        flux[i, :, m] = hp.alm2map(almc, nside)

                        almc = hp.map2alm(expo[i, :, m])
                        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                        expo[i, :, m] = hp.alm2map(almc, nside)

            path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, expo, clobber=True)

            path = os.environ["PCAT_DATA_PATH"] + '/fermflux_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, flux, clobber=True)


def writ_fdfm_doug():
    
    numbevtt = 4
    nside = 256
    numbener = 5
    numbpixl = nside**2 * 12
    fermfdfmfluxigal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxngal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxigal[0, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff0.fits')[:, None]
    fermfdfmfluxigal[1, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff1.fits')[:, None]
    fermfdfmfluxigal[2, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff2.fits')[:, None]
    fermfdfmfluxigal[3, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff3.fits')[:, None]
    fermfdfmfluxigal[4, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff4.fits')[:, None]
    
    path =os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    pf.writeto(path, fermfdfmfluxigal, clobber=True)
    
    for i in range(numbener):
        
        almc = hp.map2alm(fermfdfmfluxigal[i, :, 0])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        fermfdfmfluxngal[i, :, :] = hp.alm2map(almc, nside)[:, None]
        
    path =os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits'
    pf.writeto(path, fermfdfmfluxngal, clobber=True)
        

def writ_fdfm():
    
    nside = 256
    numbpixl = 12 * nside**2
    numbener = 3
    numbevtt = 4

    binsener = array([0.3, 1., 3., 10.])
    
    fermfdfmfluxigaltemp = tdpy.retr_fdfm(binsener, nside)

    fermfdfmfluxngaltemp = zeros((numbener, numbpixl))
    for i in range(numbener):
        almc = hp.map2alm(fermfdfmfluxigaltemp[i, :])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        fermfdfmfluxngaltemp[i, :] = hp.alm2map(almc, nside)

    fermfdfmfluxngal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxigal = zeros((numbener, numbpixl, numbevtt))
    for m in range(numbevtt):
        fermfdfmfluxigal[:, :, m] = fermfdfmfluxigaltemp
        fermfdfmfluxngal[:, :, m] = fermfdfmfluxngaltemp

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    pf.writeto(path, fermfdfmfluxigal, clobber=True)

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits'
    pf.writeto(path, fermfdfmfluxngal, clobber=True)


def retr_jcbn():
    
    lgl0, algl, bgl0, abgl, flx0, aflx = sympy.symbols('lgl0 algl bgl0 abgl flx0 aflx')
    matr = sympy.Matrix([[1, 1 - aflx / flx0, 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [1,    -aflx / flx0, 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [0,               0, 1, 1 - aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,               0, 1,    -aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,               0, 0,               0,                     1,            0], \
                         [0,               0, 0,               0,                    -1,            1]])

    jcbn = matr.det()
    return jcbn


def plot_maps():
    
    global numbpixl
    lgalheal, bgalheal, numbside, numbpixl, apix = tdpy.retr_heal(256)
    
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[:-1])
    global minmlgal, maxmlgal, minmbgal, maxmbgal, extt
    minmlgal = -10.
    maxmlgal = 10.
    minmbgal = -10.
    maxmbgal = 10.
    extt = [minmlgal, maxmlgal, minmbgal, maxmbgal]
    
    get_ipython().magic(u'matplotlib inline')

    path = os.environ["PCAT_DATA_PATH"] + '/fermflux_igal_full.fits'
    maps = sum(pf.getdata(path), 2)
    maps *= meanener[:, None]**2
    for i in range(5):
        maps[i, :] = hp.smoothing(maps[i, :], deg2rad(2.))

    plot_heal(maps[3, :] - maps[2, :])
    plot_heal(maps[2, :] - maps[1, :])


def plot_dust():
    
    redd = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/lambda_sfd_ebv.fits')['TEMPERATURE']
    
    numbside = 512
    numbpixl = numbside**2 * 12
    
    indxpixl = hp.ring2nest(numbside, arange(numbpixl))
    
    redd = redd[indxpixl]
    
    almc = hp.map2alm(redd)
    hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
    redd = hp.alm2map(almc, numbside)
    
    hist = plt.hist(redd, log=True)
    plt.show()
    
    redd[where(redd > 1.)] = 1.
    figr, axis = plt.subplots(figsize=(7, 7))
    cart = tdpy.retr_cart(redd, minmlgal=-minmgang, maxmlgal=minmgang, minmbgal=-minmgang, maxmbgal=minmgang)
    imag = axis.imshow(cart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=axis)
    plt.show()


#writ_isot()
#prep_maps()
#writ_fdfm()
#writ_fdfm_doug()
#plot_maps()
make_maps()

