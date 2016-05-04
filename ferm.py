
# coding: utf-8

# In[1]:

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


# In[2]:

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


# In[3]:

def make_maps():
    
    cmnd = 'mkdir -p $FERMI_DATA/exposure/ferm_line'
    os.system(cmnd)
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    numbproc = 40
    minmweek = 9
    maxmweek = 411
    weekpart = (maxmweek - minmweek) / numbproc
    
    global listweek
    listweek = []
    for k in range(numbproc):
        listweek.append(arange(minmweek + k * weekpart, minmweek + (k + 1) * weekpart))
    
    global strgregi
    strgregi = ' ra=INDEF dec=INDEF rad=INDEF '

    global rtag, strgtime
    
    for k in range(2):
        
        if k == 0:
            rtag = 'comp'
            reco = 7
            evtc = 2
            strgtime = 'tmin=239155201 tmax=364953603'
            
        if k == 1:
            rtag = 'full'
            reco = 8
            evtc = 128
            strgtime = 'tmin=INDEF tmax=INDEF'
            
        # process pool
        pool = mp.Pool(numbproc)

        # spawn the processes
        pool.map(make_maps_sing, range(numbproc))

        pool.close()
        pool.join()



# In[4]:

def make_maps_sing(indxproc):

    for t in listweek[indxproc]:
        
        print inspect.stack()[0][3] + ', proc%d is working on week %d...' % (indxproc, t)
        
        for m in indxevtt:
            
            if reco == 7:
                if m == 2:
                    thisevtt = 2
                elif m == 3:
                    thisevtt = 1
                else:
                    continue
            if reco == 8:
                thisevtt = evtt[m]
            
            if reco == 7:
                evnt = os.environ["FERMI_DATA"] +                     'weekly/p7v6c/lat_photon_weekly_w%03d_p202_v001.fits' % t
            if reco == 8:
                evnt = os.environ["FERMI_DATA"] +                     'weekly/photon/lat_photon_weekly_w%03d_p302_v001.fits' % t
                    
            spac = os.environ["FERMI_DATA"] +                 'weekly/spacecraft/lat_spacecraft_weekly_w%03d_p202_v001.fits' % t
                    
            sele = os.environ["FERMI_DATA"] + 'exposure/gcps_time/' + rtag +                 '/sele_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, t)
            filt = os.environ["FERMI_DATA"] + 'exposure/gcps_time/' + rtag +                 '/filt_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, t)
            live = os.environ["FERMI_DATA"] + 'exposure/gcps_time/' + rtag +                 '/live_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, t)
            cnts = os.environ["FERMI_DATA"] + 'exposure/gcps_time/' + rtag +                 '/cnts_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, t)
            expo = os.environ["FERMI_DATA"] + 'exposure/gcps_time/' + rtag +                 '/expo_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, t)

            cmnd = 'gtselect infile=' + evnt + ' outfile=' + sele + strgregi + strgtime +                 ' emin=100 emax=100000 zmax=90 evclass=%d evtype=%d' % (evtc, thisevtt)
            os.system(cmnd)

            cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1"' +                 ' outfile=' + filt + ' roicut=no'
            os.system(cmnd)

            cmnd = 'gtbin evfile=' + filt + ' scfile=' + spac + ' outfile=' + cnts +                 ' ebinalg=FILE ebinfile=/n/fink1/fermi/exposure/gcps_time/gtbndefn.fits algorithm=HEALPIX' +                 ' hpx_ordering_scheme=RING coordsys=GAL hpx_order=8 hpx_ebin=yes'
            os.system(cmnd)
            
            cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live +                 ' dcostheta=0.025 binsz=1'
            os.system(cmnd)
            
            cmnd = 'gtexpcube2 infile=' + live + ' cmap=' + cnts + ' outfile=' + expo +                 ' irfs=CALDB evtype=%03d bincalc=CENTER' % thisevtt
            os.system(cmnd)




# In[5]:

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




# In[6]:

def prep_maps_neww():
    
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




# In[7]:

def prep_maps():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener,         minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    liststrgener = ['ENERGY1', 'ENERGY2', 'ENERGY3', 'ENERGY4', 'ENERGY5']
    liststrgchan = ['CHANNEL1', 'CHANNEL2', 'CHANNEL3', 'CHANNEL4', 'CHANNEL5']

    listdatatype = ['comp', 'full']
    for datatype in listdatatype:

        if datatype == 'full':
            reco = 8
            evtc = 128
            weekinit = 11
            weekfinl = 409
            listtimefrac = array([1.])
        if datatype == 'comp':
            reco = 7
            evtc = 2
            weekinit = 9
            weekfinl = 217
            numbtime = 4
            listtimefrac = array([1., 0.75, 0.5, 0.25])

        numbtime = listtimefrac.size
                
        cntstemp = zeros((numbener, numbpixl, numbevtt))
        expotemp = zeros((numbener, numbpixl, numbevtt))
        
        cnts = zeros((numbener, numbpixl, numbevtt, numbtime))
        expo = zeros((numbener, numbpixl, numbevtt, numbtime))

        for t, timefrac in enumerate(listtimefrac):
            
            numbweek = (weekfinl - weekinit) * timefrac
            listweek = floor(linspace(weekinit, weekfinl - 1, numbweek)).astype(int)
            
            print 'numbweek'
            print numbweek
            print 'listweek'
            print listweek
            for w in listweek:

                print inspect.stack()[0][3] + ' is working on week %d...' % w 

                for m in indxevtt:
                    
                    if datatype == 'comp':
                        if m == 2:
                            thisevtt = 2
                        elif m == 3:
                            thisevtt = 1
                        else:
                            continue
                    if datatype == 'full':
                        thisevtt = evtt[m]

                    path = os.environ["FERMI_DATA"] + '/exposure/gcps_time/' + datatype +                         '/expo_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, w)

                    expoarry = pf.getdata(path, 1)
                    for i in indxener:
                        expotemp[i, :, m] = expoarry[liststrgener[i]]

                    path = os.environ["FERMI_DATA"] + '/exposure/gcps_time/' + datatype +                         '/cnts_pass%d_evtc%03d_evtt%03d_week%03d.fits' % (reco, evtc, thisevtt, w)

                    cntsarry = pf.getdata(path)
                    for i in indxener:
                        cntstemp[i, :, m] = cntsarry[liststrgchan[i]]

                expo[:, :, :, t] += expotemp
                cnts[:, :, :, t] += cntstemp
                
        flux = zeros((numbener, numbpixl, numbevtt, numbtime)) 
        indxexpo = where(expo > 0.) 
        flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
        flux /= diffener[:, None, None, None]


        listregitype = ['igal', 'ngal']
        for regitype in listregitype:

            if regitype == 'ngal':

                for i in indxener:
                    for m in indxevtt:
                        for t in range(numbtime):
                            almc = hp.map2alm(flux[i, :, m, t])
                            hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                            flux[i, :, m, t] = hp.alm2map(almc, nside)

                            almc = hp.map2alm(expo[i, :, m, t])
                            hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                            expo[i, :, m, t] = hp.alm2map(almc, nside)

            for t in range(numbtime):

                path = os.environ["PCAT_DATA_PATH"] + '/fermflux_' + regitype + '_' + datatype +                     '_time%d' % t + '.fits'
                pf.writeto(path, flux[:, :, :, t], clobber=True)
                path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_' + regitype + '_' + datatype +                     '_time%d' % t + '.fits'
                pf.writeto(path, expo[:, :, :, t], clobber=True)


# In[8]:

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
    
        
        


# In[9]:

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


# In[10]:

def retr_jcbn():
    
    lgl0, algl, bgl0, abgl, flx0, aflx = sympy.symbols('lgl0 algl bgl0 abgl flx0 aflx')
    matr = sympy.Matrix([[1, 1 - aflx / flx0, 0,               0, algl * aflx / flx0**2, -algl / flx0],                        [1,    -aflx / flx0, 0,               0, algl * aflx / flx0**2, -algl / flx0],                        [0,               0, 1, 1 - aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0],                        [0,               0, 1,    -aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0],                        [0,               0, 0,               0,                     1,            0],                        [0,               0, 0,               0,                    -1,            1]])

    jcbn = matr.det()
    return jcbn


# In[11]:

def plot_heal(maps):

    #satu = percentile(maps, 99)
    #satu = 1e9
    #maps[where(maps > satu)] = satu
    #maps[where(maps < -satu)] = -satu
    cart = tdpy.retr_cart(maps, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
    figr, axis = plt.subplots(figsize=(14, 14))
    axis.set_xlabel(r'$l$ [$^\circ$]')
    axis.set_ylabel(r'$b$ [$^\circ$]')
    imag = axis.imshow(cart, origin='lower', cmap='Reds', extent=extt)
    plt.show()
            


# In[12]:

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


# In[13]:

#writ_isot()
#make_maps()
#prep_maps()
#prep_maps_neww()
#writ_fdfm()
#writ_fdfm_doug()
#plot_maps()


# In[14]:

minmgang = 20.

data = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits')
get_ipython().magic(u'matplotlib inline')
print data.shape
for i in range(5):
    cart = tdpy.retr_cart(data[i, :, 0], minmlgal=-minmgang, maxmlgal=minmgang,                                minmbgal=-minmgang, maxmbgal=minmgang)
    figr, axis = plt.subplots(figsize=(7, 7))
    axis.set_xlabel(r'$l$ [$^\circ$]')
    axis.set_ylabel(r'$b$ [$^\circ$]')
    imag = axis.imshow(cart, origin='lower', cmap='Reds', extent=[-minmgang, minmgang, -minmgang, minmgang])
    plt.show()
    
    
redd = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/lambda_sfd_ebv.fits')['TEMPERATURE']

numbside = 512
numbpixl = numbside**2 * 12

indxpixl = hp.ring2nest(numbside, arange(numbpixl))

redd = redd[indxpixl]

almc = hp.map2alm(redd)
hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
redd = hp.alm2map(almc, numbside)

get_ipython().magic(u'matplotlib inline')
hist = plt.hist(redd, log=True)
plt.show()

redd[where(redd > 1.)] = 1.
figr, axis = plt.subplots(figsize=(7, 7))
cart = tdpy.retr_cart(redd, minmlgal=-minmgang, maxmlgal=minmgang, minmbgal=-minmgang, maxmbgal=minmgang)
imag = axis.imshow(cart, origin='lower', cmap='Reds')
plt.colorbar(imag, ax=axis)
plt.show()


# In[ ]:




# In[ ]:



