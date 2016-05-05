
# coding: utf-8

# In[1]:

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
import os, time, sys, datetime, warnings, getpass, glob, fnmatch

# tdpy
import tdpy.util

# pcat
from pcat.cnfg import *
from pcat.main import *
from pcat.samp import *
from pcat.util import *
from pcat.visu import *
from pcat.plot import *


# In[2]:

def retr_cnfg(               verbtype=1,               
              plotperd=50000, \
              
              makeplot=True, \
              
              diagsamp=False, \
              
              numbswep=1000000, \
              numbburn=None, \
              factthin=None, \
              
              regitype='ngal', \
              datatype='inpt', \
              randinit=True, \
              
              maxmgang=None, \
              
              minmspec=None, \
              maxmspec=None, \
              minmsind=None, \
              maxmsind=None, \
              
              minmfdfnnorm=None, \
              maxmfdfnnorm=None, \
              minmfdfnslop=None, \
              maxmfdfnslop=None, \

              fdfntype='singpowr', \
              
              psfntype='doubking', \
              
              colrprio=False, \
              
              numbpopl=1, \
              
              indxevttincl=arange(2, 4), \
              indxenerincl=arange(5), \
              
              maxmnumbpnts=array([200]), \

              initnumbpnts=None, \
              trueinfo=False, \
              pntscntr=False, \
              
              numbproc=None, \
              
              liketype='pois', \
              pixltype='heal', \
              exprtype='ferm', \
              
              lgalcntr=0., \
              bgalcntr=0., \
              

              
              margsize=None, \
              mocknumbpnts=None, \

              mockfdfnnorm=None, \
              mockfdfnslop=None, \

              maxmnormback=None, \
              minmnormback=None, \
              
              numbsidecart=None, \
              numbsideheal=None, \
              
              # temp
              maxmangleval=5., \
              
              spmrlbhl=2., \
            
              stdvlbhl=0.1, \
              stdvback=0.04, \
              stdvspec=0.15, \
              
              stdvfdfnnorm=0.05, \
              stdvfdfnslop=0.001, \
              stdvpsfipara=0.1, \
              
              fracrand=0.05, \
              
              mockpsfntype='doubking', \
              mocknormback=None, \
              
              exprfluxstrg=None, \
              listbackfluxstrg=None, \
              expostrg=None, \
              
              probprop=None, \
              
             ):
    
    
    # experiment-specific defaults
    if exprtype == 'ferm':
        
        # prior boundaries
        if maxmgang == None:
            maxmgang = 20.
            
        if minmsind == None:
            minmsind = array([1.])
        if maxmsind == None:
            maxmsind = array([3.])
            
            
        if minmfdfnnorm == None:
            minmfdfnnorm = 1e0
        if maxmfdfnnorm == None:
            maxmfdfnnorm = 1e2
        if minmfdfnslop == None:
            minmfdfnslop = 1.5
        if maxmfdfnslop == None:
            maxmfdfnslop = 2.5
            
        if margsize == None:
            margsize = 1.
        
    
    cnfg = dict()

    # verbosity level
    ## 0 - no output
    ## 1 - progress
    ## 2 - sampler diagnostics
    cnfg['verbtype'] = verbtype
    
    # plot settings
    ## MCMC time period over which a frame is produced
    cnfg['plotperd'] = plotperd
    ## flag to control generation of plots
    cnfg['makeplot'] = makeplot

    # flag to diagnose the sampler
    cnfg['diagsamp'] = diagsamp
    
    # MCMC setup
    ## number of sweeps, i.e., number of samples before thinning and including burn-in
    cnfg['numbswep'] = numbswep
    ## number of burn-in samples
    cnfg['numbburn'] = numbburn
    ## number of processes
    cnfg['numbproc'] = numbproc
    ## the factor by which to thin the chain
    cnfg['factthin'] = factthin
    
    # region type
    #- ngal - NGPR (North Galactic Polar Region)
    #- igal - IG (inner Galaxy)
    cnfg['regitype'] = regitype

    # data type
    #- mock - mock data
    #- inpt - input data
    cnfg['datatype'] = datatype
    
    # likelihood function type
    cnfg['liketype'] = liketype
    
    # experiment type
    cnfg['exprtype'] = exprtype
    
    # pixelization type
    cnfg['pixltype'] = pixltype
    
    # PSF model type
    cnfg['psfntype'] = psfntype
    
    if datatype == 'mock':
        cnfg['mockpsfntype'] = mockpsfntype
        
    # color priors
    cnfg['colrprio'] = colrprio
    
    # input data
    cnfg['exprfluxstrg'] = exprfluxstrg
    cnfg['listbackfluxstrg'] = listbackfluxstrg
    cnfg['expostrg'] = expostrg
    
    # flag to use truth information
    cnfg['trueinfo'] = trueinfo
    
    # Flux distribution function model
    cnfg['fdfntype'] = fdfntype
    
    # initial state setup
    ## number of point sources
    cnfg['initnumbpnts'] = initnumbpnts
    ## flag to draw the initial state from the prior
    cnfg['randinit'] = randinit

    # energy bins to be included
    cnfg['indxenerincl'] = indxenerincl
    
    # PSF bins to be included
    cnfg['indxevttincl'] = indxevttincl
    
    # number of Point Source (PS) populations
    cnfg['numbpopl'] = numbpopl
    
    # maximum angle from the PSs to evaluate the likelihood
    cnfg['maxmangleval'] = maxmangleval
    
    # parameter limits
    ## flux distribution function normalization
    cnfg['minmfdfnnorm'] = minmfdfnnorm
    cnfg['maxmfdfnnorm'] = maxmfdfnnorm
    ## flux distribution function power law index
    cnfg['minmfdfnslop'] = minmfdfnslop
    cnfg['maxmfdfnslop'] = maxmfdfnslop
    ## flux
    cnfg['minmsind'] = minmsind
    cnfg['maxmsind'] = maxmsind
    ## spectral power-law index
    cnfg['minmspec'] = minmspec
    cnfg['maxmspec'] = maxmspec
    ## background normalizations
    cnfg['minmnormback'] = minmnormback
    cnfg['maxmnormback'] = maxmnormback

    
    # model indicator limits
    cnfg['maxmnumbpnts'] = maxmnumbpnts
    
    # Region of interest
    ## image center
    cnfg['lgalcntr'] = lgalcntr
    cnfg['bgalcntr'] = bgalcntr
    ## half of the image size
    cnfg['maxmgang'] = maxmgang
    ## image margin size
    cnfg['margsize'] = margsize

    # proposals
    ## proposal scales
    cnfg['stdvfdfnnorm'] = stdvfdfnnorm
    cnfg['stdvfdfnslop'] = stdvfdfnslop
    cnfg['stdvpsfipara'] = stdvpsfipara
    cnfg['stdvback'] = stdvback
    cnfg['stdvlbhl'] = stdvlbhl
    cnfg['stdvspec'] = stdvspec

    ## fraction of heavy-tailed proposals
    cnfg['fracrand'] = fracrand
    
    ## maximum angle over which splits and merges are proposed
    cnfg['spmrlbhl'] = spmrlbhl
    
    # mock data setup
    ## mock parameters
    cnfg['mocknumbpnts'] = mocknumbpnts
    cnfg['mockfdfnnorm'] = mockfdfnnorm
    cnfg['mockfdfnslop'] = mockfdfnslop
    cnfg['mocknormback'] = mocknormback
    ## flag to position mock point sources at the image center
    cnfg['pntscntr'] = pntscntr
    ## mock image resolution
    cnfg['numbsidecart'] = numbsidecart
    cnfg['numbsideheal'] = numbsideheal

    # proposal frequencies
    cnfg['probprop'] = probprop

    return cnfg


# In[3]:

def cnfg_topo():
    
    cnfg = retr_cnfg(                      numbswep=10000,                      factthin=1,                      plotperd=20000,                      trueinfo=True,                      datatype='inpt',                      psfntype=psfntype,                      maxmgang=3.,                      minmspec=array([3e-10, 3e-11, 3e-12]),                      maxmspec=array([1e-6, 1e-7, 1e-8]),                      regitype='ngal',                      exprfluxstrg='fermflux_ngal.fits',                      listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'],                      expostrg='fermexpo_ngal.fits',                      maxmnormback=array([5., 5.]),                      minmnormback=array([0.2, 0.2]),                      stdvback=0.05,                      probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),                     )
                
    wrap(cnfg)
    
    labl = ['A', '$\psi_1$', '$\phi_1$', '$f_1$', '$\psi_2$', '$\phi_2$', '$f_2$', '$\psi_3$', '$\phi_3$', '$f_3$']
    labllist = list(itertools.combinations(labl, 2))
    indxlist = list(itertools.combinations(jpara, 2))
    ncomb = len(labllist)
    
    


# In[4]:

def cnfg_ferm_psfn_expr(psfntype):
     

    cnfg = retr_cnfg(                      numbswep=100000,                      factthin=1,                      plotperd=20000,                      trueinfo=True,                      datatype='inpt',                      psfntype=psfntype,                      maxmgang=10.,                      minmspec=array([3e-10, 3e-11, 3e-12]),                      maxmspec=array([1e-6, 1e-7, 1e-8]),                      regitype='ngal',                      exprfluxstrg='fermflux_ngal.fits',                      listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'],                      expostrg='fermexpo_ngal.fits',                      maxmnormback=array([5., 5.]),                      minmnormback=array([0.2, 0.2]),                      stdvback=0.05,                      probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),                     )
                
    wrap(cnfg)
    


# In[5]:

def cnfg_ferm_info():
    
    nruns = 2
    listlevi = zeros(nruns)
    listinfo = zeros(nruns)
    minmspec = logspace(-12., -7., nruns)

    for k in range(nruns):
        
        cnfg = retr_cnfg(                          psfntype='gausking',                          numbswep=50000,                          plotperd=50000,                          trueinfo=True,                          maxmgang=10.,                          maxmnumbpnts=array([3000]),                          colrprio=True,                          indxenerincl=arange(1),                          indxevttincl=arange(3, 4),                          minmspec=array([minmspec[k]]),                          maxmspec=array([3e-7]),                          regitype='ngal',                          maxmnormback=array([5., 5.]),                          minmnormback=array([0.2, 0.2]),                          listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'],                          expostrg='fermexpo_ngal.fits',                          stdvback=0.1,                          datatype='mock',                          mocknumbpnts=array([100]),                          numbsideheal=256,                          makeplot=False,                          mocknormback=ones((2, 3)),                         )
        
        gridchan = wrap(cnfg)
        numbproc = len(gridchan)
        for l in range(numbproc):
            listchan = gridchan[l]
            listlevi[k] = listchan[14]
            listinfo[k] = listchan[15]

    plot_minmspecinfo(minmspec, listinfo, listlevi)



# In[6]:

def cnfg_ferm_expr_igal(exprfluxstrg, expostrg):
      
    cnfg = retr_cnfg(                      psfntype='gausking',                      numbswep=3000000,                      numbburn=1500000,                      verbtype=1,                      makeplot=True,                      plotperd=50000,                      initnumbpnts=array([100]),                      maxmnumbpnts=array([600]),                      trueinfo=True,                      maxmgang=20.,                      colrprio=False,                      #indxenerincl=arange(1), \
                     #indxevttincl=arange(3, 4), \
                     #minmspec=array([1e-8]), \
                     #maxmspec=array([3e-6]), \
                     minmspec=array([1e-8, 1e-9, 1e-10, 1e-11, 1e-12]), \
                     maxmspec=array([3e-6, 3e-7, 3e-8, 3e-9, 3e-10]), \
                     regitype='igal', \
                     maxmnormback=array([2., 2.]), \
                     minmnormback=array([0.5, 0.5]), \
                     listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                     expostrg=expostrg, \
                     stdvback=0.1, \
                     datatype='inpt', \
                     exprfluxstrg=exprfluxstrg, \
                    )
        
    wrap(cnfg)


# In[7]:

def cnfg_ferm_mock_igal():
     
    cnfg = retr_cnfg(                      psfntype='singking',                      numbswep=1000000,                      plotperd=50000,                      numbsideheal=256,                      maxmgang=10.,                      minmspec=array([1e-9, 1e-10, 1e-11]),                      maxmspec=array([1e-6, 1e-7, 1e-8]),                      maxmnormback=array([5., 5.]),                      minmnormback=array([0.2, 0.2]),                      mocknormback=ones((2, 3)),                      regitype='igal',                      listbackfluxstrg=['fermisotflux.fits', 'fermfdfm.fits'],                      expostrg='fermexpo_igal.fits',                      stdvback=0.05,                      trueinfo=True,                      datatype='mock'                     )

    wrap(cnfg)
    


# In[8]:

def cnfg_ferm_expr_ngal(exprfluxstrg, expostrg):
     
    colrprio = False
    
    if colrprio:
        minmspec = array([1e-11])
        maxmspec = array([1e-7])
    else:
        minmspec = array([3e-9, 3e-10, 3e-11, 3e-12, 3e-13])
        maxmspec = array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
        
    cnfg = retr_cnfg(                      psfntype='gausking',                      numbswep=200000,                      numbburn=50000,                      verbtype=1,                      makeplot=True,                      plotperd=50000,                      initnumbpnts=array([100]),                      maxmnumbpnts=array([200]),                      trueinfo=True,                      maxmgang=20.,                      colrprio=colrprio,                      indxenerincl=arange(5),                      indxevttincl=arange(2, 4),                      minmspec=minmspec,                      maxmspec=maxmspec,                      regitype='ngal',                      maxmnormback=array([2., 2.]),                      minmnormback=array([0.5, 0.5]),                      listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'],                      expostrg=expostrg,                      stdvback=0.1,                      datatype='inpt',                      exprfluxstrg=exprfluxstrg,                     )
                
    wrap(cnfg)


# In[9]:

def cnfg_ferm_mock_ngal():
     
    colrprio = False
    
    if colrprio:
        minmspec = array([3e-11])
        maxmspec = array([1e-7])
        mockfdfnslop = array([[1.8]])
    else:
        minmspec = array([3e-9, 3e-10, 3e-11, 3e-12, 3e-13])
        maxmspec = array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
        mockfdfnslop = array([[1.8, 1.8, 1.8, 1.8, 1.8]])
      
    
    cnfg = retr_cnfg(                      psfntype='gausking',                      numbswep=200000,                      plotperd=50000,                      trueinfo=True,                      maxmgang=20.,                      colrprio=colrprio,                      verbtype=1,                      indxevttincl=arange(3, 4),                      indxenerincl=arange(5),                      maxmnumbpnts=array([200]),                      mocknumbpnts=array([100]),                      probprop=array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=float),                      minmspec=minmspec,                      maxmspec=maxmspec,                      regitype='ngal',                      maxmnormback=array([2., 2.]),                      minmnormback=array([0.5, 0.5]),                      listbackfluxstrg=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'],                      expostrg='fermexpo_ngal_comp.fits',                      stdvback=0.1,                      datatype='mock',                      numbsideheal=256,                      mockfdfnslop=mockfdfnslop,                      mockfdfnnorm=array([10.]),                      mocknormback=ones((2, 5)),                     )

    wrap(cnfg)
    


# In[10]:

def cnfg_sdss_mock():

    cnfg = retr_cnfg(psfntype='doubgaus',                      trueinfo=False,                      numbswep=100000,                      plotperd=20000,                      verbtype=1,                      minmspec=ones(3) * 1e3,                      maxmspec=ones(3) * 1e5,                      initnumbpnts=array([100]),                      exprtype='sdss',                      datatype='mock',                      pixltype='cart',                      regitype='mes5',                      stdvlbhl=2./3600.,                      lgalcntr=202.,                      bgalcntr=2.,                      mocknormback=ones((1, 3)),                      spmrlbhl=5./3600.,                      maxmnormback=array([1e3]),                      minmnormback=array([1e2]),                      maxmgang=30./3600.,                      numbsidecart=100,                      margsize=2./3600.,                      maxmangleval=10./3600.,                      listbackfluxstrg=['sdssisotflux.fits'],                      expostrg='sdssexpo.fits',                      stdvback=0.01,                      indxevttincl=arange(1),                      indxenerincl=arange(1)                     )

    wrap(cnfg)
    


# In[11]:

def cnfg_sdss_expr():

    cnfg = retr_cnfg(psfntype='doubgaus',                      trueinfo=False,                      numbswep=1000000,                      plotperd=20000,                      verbtype=1,                      minmspec=ones(3) * 1e3,                      maxmspec=ones(3) * 1e5,                      initnumbpnts=array([10]),                      maxmnumbpnts=20,                      exprtype='sdss',                      datatype='inpt',                      pixltype='cart',                      regitype='mes5',                      stdvlbhl=2./3600.,                      lgalcntr=202.,                      bgalcntr=2.,                      spmrlbhl=0.5/3600.,                      stdvspec=0.05,                      maxmnormback=array([1e3]),                      minmnormback=array([1e2]),                      margsize=2./3600.,                      maxmgang=30./3600.,                      maxmangleval=10./3600.,                      exprfluxstrg='sdssflux.fits',                      listbackfluxstrg=['sdssisotflux.fits'],                      expostrg='sdssexpo.fits',                      stdvback=1e-4,                      indxevttincl=arange(1),                      indxenerincl=arange(1)                     )

    wrap(cnfg)


# In[12]:

if __name__ == '__main__':
    
    pass

    #cnfg_ferm_info()
    
    #cnfg_ferm_psfn_mock('gausking')
    #cnfg_ferm_psfn_mock('doubking')

    #cnfg_ferm_psfn_expr('gausking')
    #cnfg_ferm_psfn_expr('doubking')
    
    #cnfg_ferm_expr_igal('fermflux_igal_comp_time0.fits', 'fermexpo_igal_comp_time0.fits')
    #cnfg_ferm_mock_igal()
    
    #cnfg_ferm_expr_ngal('fermflux_comp_ngal.fits', 'fermexpo_comp_ngal.fits')
    #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time0.fits', 'fermexpo_ngal_comp_time0.fits')
    #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time1.fits', 'fermexpo_ngal_comp_time1.fits')
    #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time2.fits', 'fermexpo_ngal_comp_time2.fits')
    #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time3.fits', 'fermexpo_ngal_comp_time3.fits')
    #cnfg_ferm_expr_ngal('fermflux_ngal_full.fits', 'fermexpo_ngal_full.fits')
    
    cnfg_ferm_mock_ngal()
    
    #cnfg_sdss_mock()
    #cnfg_sdss_expr()


# In[ ]:




# In[ ]:



