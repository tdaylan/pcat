
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
import os, time, sys, datetime, warnings, getpass, glob, fnmatch, cPickle
import functools

# tdpy
import tdpy.util
import tdpy.mcmc

# pcat
from util import *
from visu import *


def work(globdata, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for the process
    seed()
    
    # construct the run tag
    globdata.rtag = retr_rtag(globdata, indxprocwork)
    
    # initialize the sample vector 
    if globdata.randinit or not globdata.trueinfo:
        if globdata.initnumbpnts != None:
            thisnumbpnts = globdata.initnumbpnts
        else:
            thisnumbpnts = empty(globdata.numbpopl, dtype=int)
            for l in globdata.indxpopl:
                thisnumbpnts[l] = choice(arange(globdata.minmnumbpnts, globdata.maxmnumbpnts[l] + 1))
    else:
        thisnumbpnts = globdata.truenumbpnts
        
    globdata.thisindxpntsfull = []
    globdata.thisindxpntsempt = []
    for l in globdata.indxpopl:
        globdata.thisindxpntsfull.append(range(thisnumbpnts[l]))
        globdata.thisindxpntsempt.append(range(thisnumbpnts[l], globdata.maxmnumbpnts[l]))
    globdata.thisindxsamplgal, globdata.thisindxsampbgal, globdata.thisindxsampspec, \
        globdata.thisindxsampsind, globdata.thisindxsampcomp = retr_indx(globdata, globdata.thisindxpntsfull)
    
    globdata.drmcsamp = zeros((globdata.maxmsampsize, 2))
    
    globdata.drmcsamp[globdata.indxsampnumbpnts, 0] = thisnumbpnts
    globdata.drmcsamp[globdata.indxsampfdfnnorm, 0] = rand(globdata.numbpopl)
    if globdata.trueinfo and globdata.datatype == 'mock':
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = cdfn_atan(globdata.truefdfnslop, globdata.minmfdfnslop, globdata.factfdfnslop)
    else:
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = rand(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener))
    globdata.drmcsamp[globdata.indxsampnormback, 0] = rand(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener))
    if globdata.randinit or not globdata.trueinfo or globdata.truepsfipara == None:
        globdata.drmcsamp[globdata.indxsamppsfipara, 0] = rand(globdata.numbpsfipara)
    else:
        for k in globdata.indxmodlpsfipara:
            globdata.drmcsamp[globdata.indxsamppsfipara[k], 0] = cdfn_psfipara(globdata, globdata.truepsfipara[k], k)
        
    for l in globdata.indxpopl:
        if globdata.randinit or not globdata.trueinfo:
            globdata.drmcsamp[globdata.thisindxsampcomp[l], 0] = rand(globdata.thisindxsampcomp[l].size)
        else:
            globdata.drmcsamp[globdata.thisindxsamplgal[l], 0] = copy(cdfn_self(globdata.truelgal[l], 
            -globdata.maxmgangmarg, 
            2. * globdata.maxmgangmarg))
            globdata.drmcsamp[globdata.thisindxsampbgal[l], 0] = copy(cdfn_self(globdata.truebgal[l], 
            -globdata.maxmgangmarg, 
            2. * globdata.maxmgangmarg))
            for i in globdata.indxenerfdfn:
                fdfnslop = icdf_atan(globdata.drmcsamp[globdata.indxsampfdfnslop[l, i], 0], globdata.minmfdfnslop, globdata.factfdfnslop)
                spec = cdfn_spec(globdata, globdata.truespec[l][0, i, :], fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])
                globdata.drmcsamp[globdata.thisindxsampspec[l][i, :], 0] = copy(spec)
            if globdata.colrprio:
                globdata.drmcsamp[globdata.thisindxsampsind[l], 0] = copy(cdfn_atan(globdata.truesind[l],                 globdata.minmsind,                 globdata.factsind))
                flux = globdata.drmcsamp[globdata.thisindxsampspec[l][globdata.indxenerfdfn, :], 0]
                sind = globdata.drmcsamp[globdata.thisindxsampsind[l], 0]
                globdata.drmcsamp[globdata.thisindxsampspec[l], 0] = retr_spec(flux, sind)

    
    globdata.thissampvarb, thisindxpixlpnts, thiscnts, globdata.thispntsflux,         thismodlflux, globdata.thismodlcnts = pars_samp(globdata, globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])

        
    if globdata.verbtype > 2:
        print 'thisindxpntsfull: ', globdata.thisindxpntsfull
        print 'thisindxpntsempt: ', globdata.thisindxpntsempt  
        print 'thisindxsamplgal: ', globdata.thisindxsamplgal
        print 'thisindxsampbgal: ', globdata.thisindxsampbgal
        print 'thisindxsampspec: '
        print globdata.thisindxsampspec
        if globdata.colrprio:
            print 'thisindxsampsind: ', globdata.thisindxsampsind
        print 'thisindxsampcomp: ', globdata.thisindxsampcomp


        


    globdata.nextpntsflux = zeros_like(globdata.thispntsflux)
    globdata.nextmodlflux = zeros_like(globdata.thispntsflux)
    globdata.nextmodlcnts = zeros_like(globdata.thispntsflux)
    globdata.nextllik = zeros_like(globdata.thispntsflux)

    globdata.nextsampvarb = copy(globdata.thissampvarb)
    
    if globdata.verbtype > 1:
        print 'thissampvarb: ', globdata.thissampvarb
        
    # sampler setup
    # auxiliary variable standard globdata.deviation for merge/split
    globdata.maxmdistpnts = 2. # [deg]
 
    listchan = rjmc(globdata, indxprocwork)
    
    timereal = time.time() - timereal
    timeproc = time.clock() - timeproc
    
    listchan.append(timereal)
    listchan.append(timeproc)
    
    return listchan


def wrap(cnfg):
    
    timetotlreal = time.time()
    timetotlproc = time.clock()
    
    globdata = globdatastrt()
    
    globdata.verbtype = cnfg['verbtype']
    globdata.plotperd = cnfg['plotperd']
    globdata.makeplot = cnfg['makeplot']
    globdata.diagsamp = cnfg['diagsamp']
    
    globdata.numbproc = cnfg['numbproc']
    globdata.numbswep = cnfg['numbswep']
    globdata.numbburn = cnfg['numbburn']
    globdata.factthin = cnfg['factthin']
    
    globdata.stdvfdfnnorm = cnfg['stdvfdfnnorm']
    globdata.stdvfdfnslop = cnfg['stdvfdfnslop']
    globdata.stdvpsfipara = cnfg['stdvpsfipara']
    globdata.stdvback = cnfg['stdvback']
    globdata.stdvlbhl = cnfg['stdvlbhl']
    globdata.stdvspec = cnfg['stdvspec']
    
    globdata.spmrlbhl = cnfg['spmrlbhl']
    
    globdata.fracrand = cnfg['fracrand']
    
    globdata.datatype = cnfg['datatype']
    
    if globdata.datatype == 'mock':
        mockpsfntype = cnfg['mockpsfntype']
    
    globdata.psfntype = cnfg['psfntype']
    
    globdata.liketype = cnfg['liketype']
    globdata.exprtype = cnfg['exprtype']
    globdata.pixltype = cnfg['pixltype']
    
    globdata.regitype = cnfg['regitype']
    globdata.randinit = cnfg['randinit']
    
    globdata.fdfntype = cnfg['fdfntype']

    globdata.colrprio = cnfg['colrprio']
    globdata.numbpopl = cnfg['numbpopl']
    globdata.indxenerincl = cnfg['indxenerincl']
    globdata.indxevttincl = cnfg['indxevttincl']
    
    globdata.maxmangleval = cnfg['maxmangleval']

    globdata.minmspec = cnfg['minmspec']
    globdata.maxmspec = cnfg['maxmspec']
    if globdata.colrprio:
        globdata.minmsind = cnfg['minmsind']
        globdata.maxmsind = cnfg['maxmsind']
    
    globdata.minmfdfnnorm = cnfg['minmfdfnnorm']
    globdata.maxmfdfnnorm = cnfg['maxmfdfnnorm']
    globdata.minmfdfnslop = cnfg['minmfdfnslop']
    globdata.maxmfdfnslop = cnfg['maxmfdfnslop']
    

    globdata.maxmnormback = cnfg['maxmnormback']
    globdata.minmnormback = cnfg['minmnormback']
    
    if globdata.datatype == 'mock':
        
        mockfdfnslop = cnfg['mockfdfnslop']
        mocknumbpnts = cnfg['mocknumbpnts']
        mockfdfnnorm = cnfg['mockfdfnnorm']
        mocknormback = cnfg['mocknormback']

    globdata.maxmnumbpnts = cnfg['maxmnumbpnts']
    
    globdata.probprop = cnfg['probprop']
    
    exprfluxstrg = cnfg['exprfluxstrg']
    listbackfluxstrg = cnfg['listbackfluxstrg']
    globdata.expostrg = cnfg['expostrg']
    
    globdata.initnumbpnts = cnfg['initnumbpnts']
    globdata.trueinfo = cnfg['trueinfo']
    
    globdata.margsize = cnfg['margsize']
    globdata.maxmgang = cnfg['maxmgang']
    globdata.lgalcntr = cnfg['lgalcntr']
    globdata.bgalcntr = cnfg['bgalcntr']
    
    globdata.numbsideheal = cnfg['numbsideheal']
    globdata.numbsidecart = cnfg['numbsidecart']
    
    if globdata.colrprio:
        if globdata.minmspec.size != 1 or globdata.minmspec.size != 1:
            print 'Spectral limits must be numpy scalars when color priors are used!'
        
    # the number of processes (each with globdata.numbswep samples)
    if globdata.numbproc == None:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu':
            globdata.numbproc = 10
        else:
            globdata.numbproc = 1
        
    # date and time
    globdata.strgtime = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
    if globdata.verbtype > 0:
        print 'PCAT started at ', globdata.strgtime
        print 'Initializing...'
    
    globdata.strgfluxunit = retr_strgfluxunit(globdata)
        
    # number of bins
    globdata.numbspec = 10
    
    globdata.numbbins = 10
    
    if globdata.datatype == 'mock':
        mocknormback = mocknormback[:, globdata.indxenerincl]
        if not globdata.colrprio:
            mockfdfnslop = mockfdfnslop[:, globdata.indxenerincl]

    if not globdata.colrprio:
        globdata.minmspec = globdata.minmspec[globdata.indxenerincl]
        globdata.maxmspec = globdata.maxmspec[globdata.indxenerincl]
        

    globdata.minmnumbpnts = 1

    globdata.numbback = len(listbackfluxstrg)
    globdata.indxback = arange(globdata.numbback)
    
    if globdata.numbback == 1:
        globdata.probprop[5] = 0.
        
    # PSF class axis
    globdata.numbevtt = globdata.indxevttincl.size
    globdata.indxevtt = arange(globdata.numbevtt)
    
    # energy axis
    globdata.numbener = globdata.indxenerincl.size
    globdata.indxenerinclbins = empty(globdata.numbener+1, dtype=int)
    globdata.indxenerinclbins[0:-1] = globdata.indxenerincl
    globdata.indxenerinclbins[-1] = globdata.indxenerincl[-1] + 1
    globdata.binsener = array([0.1, 0.3, 1., 3., 10., 100.])[globdata.indxenerinclbins]
    globdata.diffener = (roll(globdata.binsener, -1) - globdata.binsener)[0:-1]

    globdata.meanener = sqrt(roll(globdata.binsener, -1) * globdata.binsener)[0:-1]
    globdata.indxener = arange(globdata.numbener, dtype=int)
    
    if globdata.colrprio: 
        # temp
        globdata.indxenerfdfn = array([globdata.numbener / 2])
        f
        globdata.minmspec = globdata.minmspec * (globdata.meanener[globdata.indxenerfdfn] / globdata.meanener)**2
        globdata.maxmspec = globdata.maxmspec * (globdata.meanener[globdata.indxenerfdfn] / globdata.meanener)**2
        
        if globdata.datatype == 'mock':
            
            mockfdfnslop = tile(mockfdfnslop, (1, globdata.numbener))
        
    else:
        globdata.indxenerfdfn = globdata.indxener

    globdata.enerfdfn = globdata.meanener[globdata.indxenerfdfn]
        
    if globdata.colrprio:
        globdata.indxenerprio = globdata.indxenerfdfn
    else:
        globdata.indxenerprio = globdata.indxener

    globdata.indxspecpivt = globdata.numbspec / 2
    
    # temp
    if globdata.exprtype == 'sdss':
        globdata.diffener = ones(globdata.numbener)
        
    # angular globdata.deviation
    globdata.numbangl = 100
    if globdata.exprtype == 'sdss':
        globdata.maxmangldisp = deg2rad(5. / 3600.) # [rad]
    if globdata.exprtype == 'ferm':
        globdata.maxmangldisp = deg2rad(5.) # [rad]
    globdata.angldisp = linspace(0., globdata.maxmangldisp, globdata.numbangl) # [rad]
    maxmangl = deg2rad(3.5 * globdata.maxmgang) # [rad]
    angl = linspace(0., maxmangl, globdata.numbangl) # [rad]

    if globdata.exprtype == 'ferm':
        globdata.strganglunit = '[deg]'
    if globdata.exprtype == 'sdss':
        globdata.strganglunit = '[arcsec]'
        
    if globdata.exprtype == 'ferm':
        globdata.fermpsfn, fermpsfipara = retr_fermpsfn(globdata)

    # energy bin string
    globdata.enerstrg, globdata.binsenerstrg = retr_enerstrg(globdata)
    
    # PSF class string
    globdata.evttstrg = []
    for m in globdata.indxevtt:
        globdata.evttstrg.append('PSF%d' % globdata.indxevttincl[m])
        
    globdata.mrkralph = 0.8

    globdata.spltjcbnfact = log(2.**(2 - globdata.numbener))
    
    # number of components
    globdata.numbcomp = 2 + globdata.numbener
    if globdata.colrprio:
        globdata.numbcomp += 1
    globdata.numbcompcolr = 4
    globdata.jcbnsplt = 2.**(2 - globdata.numbener)
    
    # limits on the number of point sources
    globdata.factfdfnnorm = log(globdata.maxmfdfnnorm / globdata.minmfdfnnorm)
    

    if globdata.datatype == 'mock':

        # mock PSF parameters
        if mockpsfntype == 'singgaus':
            globdata.numbmockformpara = 1
        elif mockpsfntype == 'singking':
            globdata.numbmockformpara = 2 
        elif mockpsfntype == 'doubgaus':
            globdata.numbmockformpara = 3
        elif mockpsfntype == 'gausking':
            globdata.numbmockformpara = 4
        elif mockpsfntype == 'doubking':
            globdata.numbmockformpara = 5

        globdata.numbmockpsfiparaevtt = globdata.numbener * globdata.numbmockformpara
        globdata.numbmockpsfipara = globdata.numbmockpsfiparaevtt * globdata.numbevtt
        globdata.indxmockpsfipara = arange(globdata.numbmockpsfipara)   


    if globdata.regitype == 'igal':
        globdata.longlabl = '$l$'
        globdata.latilabl = '$b$'
    else:
        globdata.longlabl = r'$\nu$'
        globdata.latilabl = r'$\mu$'
        
    if globdata.exprtype == 'ferm':
        globdata.longlabl += r' [$^\circ$]'
        globdata.latilabl += r' [$^\circ$]'
    if globdata.exprtype == 'sdss':
        globdata.longlabl += ' [arcsec]'
        globdata.latilabl += ' [arcsec]'
    
    # construct the PSF model
    retr_psfimodl(globdata)

    # proposals
    retr_propmodl(globdata)
    
    # factors in the prior expression
    globdata.priofactlgalbgal = 2. * log(1. / 2. / globdata.maxmgang)
    globdata.priofactfdfnslop = globdata.numbener * log(1. / (arctan(globdata.maxmfdfnslop) - arctan(globdata.minmfdfnslop)))
    globdata.priofactfdfnnorm = log(1. / (log(globdata.maxmfdfnnorm) - log(globdata.minmfdfnnorm)))

    # sample vector indices  
    globdata.indxsampnumbpnts = arange(globdata.numbpopl)
    globdata.indxsampfdfnnorm = arange(globdata.numbpopl) + amax(globdata.indxsampnumbpnts) + 1
    globdata.indxsampfdfnslop = arange(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener)) + amax(globdata.indxsampfdfnnorm) + 1
    globdata.indxsamppsfipara = arange(globdata.numbpsfipara) + amax(globdata.indxsampfdfnslop) + 1
    globdata.indxsampnormback = arange(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener)) + amax(globdata.indxsamppsfipara) + 1

    globdata.fluxpivt = sqrt(globdata.minmspec * globdata.maxmspec)
    
    globdata.maxmnumbcomp = globdata.maxmnumbpnts * globdata.numbcomp
    globdata.indxsampcompinit = amax(globdata.indxsampnormback) + 1
    
    # maximum number of parameters
    globdata.maxmsampsize = int(globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.numbpopl)
    
    if globdata.numbburn == None:
        globdata.numbburn = globdata.numbswep / 5
    if globdata.factthin == None:
        globdata.factthin = min(globdata.maxmsampsize * 5, globdata.numbswep / 2)



    
    if globdata.verbtype > 1:
        print 'probprop: '
        print vstack((arange(globdata.numbprop), globdata.strgprop, globdata.probprop)).T


    # run tag
    globdata.rtag = retr_rtag(globdata, None)
    
    # plots
    if globdata.makeplot:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
            plotfold = '/n/pan/www/tansu/png/pcat/'
        else:
            plotfold = os.environ["PCAT_DATA_PATH"] + '/png/'
        globdata.plotpath = plotfold + globdata.strgtime + '_' + globdata.rtag + '/'
        cmnd = 'mkdir -p ' + globdata.plotpath
        os.system(cmnd)

    globdata.factfdfnslop = arctan(globdata.maxmfdfnslop) - arctan(globdata.minmfdfnslop)
    if globdata.colrprio:
        globdata.factsind = arctan(globdata.maxmsind) - arctan(globdata.minmsind)

    # number of samples to be saved
    globdata.numbsamp = (globdata.numbswep - globdata.numbburn) / globdata.factthin

    # rescale the positional update scale
    globdata.stdvlbhl /= 2. * globdata.maxmgang

    # determine proposal probabilities
    globdata.probpropminm = copy(globdata.probprop)
    globdata.probpropmaxm = copy(globdata.probprop)
    
    ## do not allow death or merge when the number of PS is at its minimum or 
    ## births or splits when the number of PS is at its maximum
    globdata.probpropminm[[globdata.indxpropdeth, globdata.indxpropmerg]] = 0.
    globdata.probpropmaxm[[globdata.indxpropbrth, globdata.indxpropsplt]] = 0.
    globdata.probprop /= sum(globdata.probprop)
    globdata.probpropmaxm /= sum(globdata.probpropmaxm)
    globdata.probpropminm /= sum(globdata.probpropminm)
    
    # population indices
    globdata.indxpopl = arange(globdata.numbpopl)

    globdata.maxmgangmarg = globdata.maxmgang + globdata.margsize
    
    globdata.minmlgalmarg = -globdata.maxmgangmarg
    globdata.maxmlgalmarg = globdata.maxmgangmarg
    globdata.minmbgalmarg = -globdata.maxmgangmarg
    globdata.maxmbgalmarg = globdata.maxmgangmarg
    globdata.minmlgal = -globdata.maxmgang
    globdata.maxmlgal = globdata.maxmgang
    globdata.minmbgal = -globdata.maxmgang
    globdata.maxmbgal = globdata.maxmgang
    

    # input data
    if globdata.datatype == 'inpt':
        
        path = os.environ["PCAT_DATA_PATH"] + '/' + exprfluxstrg
        exprflux = pf.getdata(path)

        globdata.indxenerinclfull = arange(exprflux.shape[0])
        globdata.indxevttinclfull = arange(exprflux.shape[2])

        if globdata.pixltype == 'heal':
            globdata.numbpixlheal = exprflux.shape[1]
            globdata.numbsideheal = int(sqrt(globdata.numbpixlheal / 12))
        else:
            globdata.numbsidecart = exprflux.shape[1]
            exprflux = exprflux.reshape((exprflux.shape[0], globdata.numbsidecart**2, exprflux.shape[3]))
            
    else:
        
        if globdata.exprtype == 'ferm':
            globdata.indxenerinclfull = arange(5)
            globdata.indxevttinclfull = arange(4)
         
    if globdata.pixltype == 'heal':
        globdata.numbpixlheal = globdata.numbsideheal**2 * 12
        globdata.apix = 4. * pi / globdata.numbpixlheal
    if globdata.pixltype == 'cart':
        globdata.binslgalcart = linspace(globdata.minmlgal, globdata.maxmlgal, globdata.numbsidecart + 1)
        globdata.binsbgalcart = linspace(globdata.minmbgal, globdata.maxmbgal, globdata.numbsidecart + 1)
        globdata.lgalcart = (globdata.binslgalcart[0:-1] + globdata.binslgalcart[1:]) / 2.
        globdata.bgalcart = (globdata.binsbgalcart[0:-1] + globdata.binsbgalcart[1:]) / 2.
        globdata.apix = deg2rad(2. * globdata.maxmgang / globdata.numbsidecart)**2
        
    # temp
    globdata.tracsamp = False
    
    # center of the ROI
    if globdata.regitype == 'igal':
        globdata.cntrlghp, globdata.cntrbghp = 0., 0.
    else:
        globdata.cntrlghp, globdata.cntrbghp = 0., 90.
    
    
    # plot settings
    
    ## marker size
    globdata.minmmrkrsize = 50
    globdata.maxmmrkrsize = 500
    
    ## roi
    globdata.exttrofi = [globdata.minmlgal, globdata.maxmlgal, globdata.minmbgal, globdata.maxmbgal]
    if globdata.exprtype == 'sdss':
        globdata.exttrofi *= 3600.
        globdata.frambndr = globdata.maxmgang * 3600.
        globdata.frambndrmarg = globdata.maxmgangmarg * 3600.
    else:
        globdata.frambndr = globdata.maxmgang
        globdata.frambndrmarg = globdata.maxmgangmarg
      
    
    
    ## FDM normalization prior limits
    globdata.factnormback = log(globdata.maxmnormback / globdata.minmnormback)

    # sky coordinates
    globdata.binslbhl = linspace(-globdata.maxmgang, globdata.maxmgang, globdata.numbbins + 1)

    globdata.binsspec = zeros((globdata.numbener, globdata.numbspec + 1))
    for i in globdata.indxener:
        globdata.binsspec[i, :] = logspace(log10(globdata.minmspec[i]), log10(globdata.maxmspec[i]), globdata.numbspec + 1)
    globdata.meanspec = sqrt(globdata.binsspec[:, 1:] * globdata.binsspec[:, 0:-1])
    globdata.diffspec = globdata.binsspec[:, 1:] - globdata.binsspec[:, 0:-1]

    if globdata.colrprio:
        globdata.binssindunit = arange(globdata.numbbins + 1.) / globdata.numbbins
        globdata.meansindunit = (globdata.binssindunit[:-1] + globdata.binssindunit[1:]) / 2.
        globdata.binssind = icdf_atan(globdata.binssindunit, globdata.minmsind, globdata.factsind)
        globdata.meansind = icdf_atan(globdata.meansindunit, globdata.minmsind, globdata.factsind)
        globdata.diffsind = globdata.binssind[1:] - globdata.binssind[:-1]

    # temp
    globdata.numbspecprox = 1
    globdata.meanspecprox = globdata.binsspec[globdata.indxenerfdfn, globdata.numbspec / 2]
    
    if globdata.verbtype > 1:
        print 'indxsampnumbpnts: ', globdata.indxsampnumbpnts
        print 'indxsampfdfnnorm: ', globdata.indxsampfdfnnorm
        print 'indxsampfdfnslop: ', globdata.indxsampfdfnslop
        print 'indxsamppsfipara: ', globdata.indxsamppsfipara
        print 'indxsampnormback: '
        print globdata.indxsampnormback
        print 'indxsampcompinit: ', globdata.indxsampcompinit

        
    # region of interest
    if globdata.pixltype == 'heal':
        
        lgalheal, bgalheal, globdata.numbsideheal, globdata.numbpixlheal, globdata.apix = tdpy.util.retr_healgrid(globdata.numbsideheal)

        globdata.indxpixlrofi = where((abs(lgalheal) < globdata.maxmgang) & (abs(bgalheal) < globdata.maxmgang))[0]
        globdata.indxpixlrofimarg = where((abs(lgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal) &                               (abs(bgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal))[0]

        
        globdata.lgalgrid = lgalheal[globdata.indxpixlrofi]
        globdata.bgalgrid = bgalheal[globdata.indxpixlrofi]
        
        path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%03d.p' % globdata.maxmgang
        if os.path.isfile(path):
            fobj = open(path, 'rb')
            globdata.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            globdata.pixlcnvt = zeros(globdata.numbpixlheal, dtype=int)
            for k in range(globdata.indxpixlrofimarg.size):
                dist = retr_dist(globdata, lgalheal[globdata.indxpixlrofimarg[k]],                                  bgalheal[globdata.indxpixlrofimarg[k]], globdata.lgalgrid, globdata.bgalgrid)
                globdata.pixlcnvt[globdata.indxpixlrofimarg[k]] = argmin(dist)

            fobj = open(path, 'wb')
            cPickle.dump(globdata.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()


            
    else:
        isidecart = arange(globdata.numbsidecart)
        temp = meshgrid(isidecart, isidecart, indexing='ij')
        globdata.bgalgrid = globdata.bgalcart[temp[1].flatten()]
        globdata.lgalgrid = globdata.lgalcart[temp[0].flatten()]
        
    globdata.numbpixl = globdata.lgalgrid.size
    globdata.indxpixl = arange(globdata.numbpixl)
    globdata.fullindx = meshgrid(globdata.indxener, globdata.indxpixl, globdata.indxevtt, indexing='ij')
    
    globdata.filtindx = meshgrid(globdata.indxenerincl, globdata.indxpixl, globdata.indxevttincl, indexing='ij')
    
    globdata.numbpixlsave = 1000
    globdata.indxpixlsave = choice(arange(globdata.numbpixlsave), size=globdata.numbpixlsave)


    if globdata.pixltype == 'heal':
        healindx = meshgrid(globdata.indxenerinclfull, globdata.indxpixlrofi, globdata.indxevttinclfull, indexing='ij')
        

    if globdata.datatype == 'inpt' and globdata.pixltype == 'heal':
        exprflux = exprflux[healindx]


    if globdata.datatype == 'inpt':
        exprflux = exprflux[globdata.filtindx]
 
    if globdata.datatype == 'mock' or globdata.exprtype == 'ferm':
        mockpsfipara = fermpsfipara

    # globdata.exposure
    if globdata.expostrg == 'unit':
        globdata.expo = ones((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/' + globdata.expostrg
        globdata.expo = pf.getdata(path)

        if globdata.pixltype == 'heal':
            globdata.expo = globdata.expo[healindx]
        else:
            globdata.expo = globdata.expo.reshape((globdata.expo.shape[0], -1, globdata.expo.shape[-1]))
            
        globdata.expo = globdata.expo[globdata.filtindx]
    

    # backgrounds
    globdata.backflux = []
    globdata.backfluxmean = []
    for c, backfluxstrg in enumerate(listbackfluxstrg):
        path = os.environ["PCAT_DATA_PATH"] + '/' + backfluxstrg
        backfluxtemp = pf.getdata(path)
        if globdata.pixltype == 'heal':
            backfluxtemp = backfluxtemp[healindx]
        else:
            backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        backfluxtemp = backfluxtemp[globdata.filtindx]
        globdata.backflux.append(backfluxtemp)
        globdata.backfluxmean.append(mean(sum(backfluxtemp * globdata.expo, 2) / sum(globdata.expo, 2), 1))
        

        
    # test plot
    # temp
    if globdata.datatype == 'inpt' and globdata.makeplot and False:
        for i in globdata.indxener:
            for m in globdata.indxevtt:
                globdata.backfluxtemp = zeros(globdata.numbpixl)
                for c in globdata.indxback:
                    globdata.backfluxtemp[:] += globdata.backflux[c][i, :, m]
                resicnts = (exprflux[i, :, m] - globdata.backfluxtemp) * globdata.expo[i, :, m] * globdata.apix * globdata.diffener[i]

                
                resicntstemp = tdpy.util.retr_cart(resicnts, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(resicntstemp, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testresiflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                cart = tdpy.util.retr_cart(exprflux[i, :, m], indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testexprflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                cart = tdpy.util.retr_cart(globdata.backfluxtemp, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testglobdata.backflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)

                cart = tdpy.util.retr_cart(globdata.expo[i, :, m], indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testglobdata.expo%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                
                
    # get 3FGL catalog
    if globdata.exprtype == 'ferm':
        globdata.fgl3lgal, globdata.fgl3bgal, globdata.fgl3spec, fgl3gang,             globdata.fgl3cnts, globdata.fgl3timevari, globdata.fgl3sind, globdata.fgl3spectype, globdata.fgl3scur, globdata.fgl3scut = retr_fgl3(globdata)
        
        if globdata.regitype == 'ngal':
            rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
            globdata.fgl3bgal, globdata.fgl3lgal = rad2deg(rttr(deg2rad(90. - globdata.fgl3bgal), deg2rad(globdata.fgl3lgal)))
            globdata.fgl3bgal = 90. - globdata.fgl3bgal

        globdata.indxfgl3rofi = arange(globdata.fgl3lgal.size, dtype=int)
        for i in globdata.indxener:
            globdata.indxfgl3rofi = intersect1d(where((globdata.fgl3spec[0, i, :] > globdata.minmspec[i]) & (globdata.fgl3spec[0, i, :] < globdata.maxmspec[i]))[0], globdata.indxfgl3rofi)
        globdata.indxfgl3rofi = intersect1d(where((abs(globdata.fgl3lgal) < globdata.maxmgangmarg) & (abs(globdata.fgl3bgal) < globdata.maxmgangmarg))[0], globdata.indxfgl3rofi)

        globdata.indxfgl3timevari = where(globdata.fgl3timevari > 72.44)[0]
        

        if globdata.makeplot:
            plot_fgl3(globdata)
        
        globdata.fgl3lgal = globdata.fgl3lgal[globdata.indxfgl3rofi]
        globdata.fgl3bgal = globdata.fgl3bgal[globdata.indxfgl3rofi]
        globdata.fgl3spec = globdata.fgl3spec[:, :, globdata.indxfgl3rofi]
        globdata.fgl3cnts = globdata.fgl3cnts[:, globdata.indxfgl3rofi, :]
        globdata.fgl3spectype = globdata.fgl3spectype[globdata.indxfgl3rofi]
        globdata.fgl3scur = globdata.fgl3scur[globdata.indxfgl3rofi]
        globdata.fgl3scut = globdata.fgl3scut[globdata.indxfgl3rofi]
        globdata.fgl3timevari = globdata.fgl3timevari[globdata.indxfgl3rofi]

        globdata.indxfgl3timevari = where(globdata.fgl3timevari > 72.44)[0]

        globdata.fgl3numbpnts = globdata.fgl3lgal.size
        
        for i in globdata.indxener:
            if (globdata.fgl3spec[0, i, :] > globdata.maxmspec[i]).any():
                print 'maxmspec %d is bad!' % i


    # true data
    if globdata.datatype == 'inpt':
        globdata.truenumbpnts = None
        globdata.truefdfnnorm = None
        globdata.truefdfnslop = None
        globdata.truenormback = None
    
    # get count data

    ## input data
    if globdata.datatype == 'inpt':
        globdata.datacnts = exprflux * globdata.expo * globdata.apix * globdata.diffener[:, None, None] # [1]

    ## mock data
    if globdata.datatype == 'mock':

        if mocknumbpnts == None:
            mocknumbpnts = empty(globdata.numbpopl)
            for l in globdata.indxpopl:
                mocknumbpnts[l] = random_integers(globdata.minmnumbpnts, globdata.maxmnumbpnts[l])
        
        if mockfdfnnorm == None:
            mockfdfnnorm = icdfn_logt(rand(globdata.numbpopl), globdata.minmfdfnnorm, globdata.factfdfnnorm)
            
        if mockfdfnslop == None:
            temp = rand(globdata.numbener * globdata.numbpopl)
            mockfdfnslop = icdf_atan(temp, globdata.minmfdfnslop, globdata.factfdfnslop)
            
        if mockpsfipara == None: 
            mockpsfntype = psfntpye
            numbmockpsfipara = globdata.numbpsfipara
            mockpsfipara = empty(numbmockpsfipara)
            for k in arange(numbmockpsfipara):
                mockpsfipara[k] = icdf_psfipara(globdata, rand(), k)
   
        if mocknormback == None:
            for c in globdata.indxback:
                mocknormback[c, :] = icdf_logt(rand(globdata.numbener),                                                     globdata.minmnormback[c], globdata.factnormback[c])

        mockcnts = [[] for l in globdata.indxpopl]
        mocklgal = [[] for l in globdata.indxpopl]
        mockbgal = [[] for l in globdata.indxpopl]
        mockspec = [[] for l in globdata.indxpopl]
        if globdata.colrprio:
            mocksind = [[] for l in globdata.indxpopl]
        for l in globdata.indxpopl:
            mocklgal[l] = icdf_self(rand(mocknumbpnts[l]), -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            mockbgal[l] = icdf_self(rand(mocknumbpnts[l]), -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg) 
            mockspec[l] = empty((globdata.numbener, mocknumbpnts[l]))
            for i in globdata.indxenerfdfn:
                mockspec[l][i, :] = icdf_spec(globdata, rand(mocknumbpnts[l]),                                               mockfdfnslop[l, i], globdata.minmspec[i], globdata.maxmspec[i])
            if globdata.colrprio:
                mocksind[l] = icdf_atan(rand(mocknumbpnts[l]), globdata.minmsind, globdata.factsind)
                mockspec[l] = retr_spec(mockspec[l][globdata.indxenerfdfn, :], mocksind[l])
            indxpixltemp = retr_indxpixl(globdata, mockbgal[l], mocklgal[l])
            mockcnts[l] = mockspec[l][:, :, None] * globdata.expo[:, indxpixltemp, :] * globdata.diffener[:, None, None]
        mockpntsflux = retr_pntsflux(globdata, concatenate(mocklgal), concatenate(mockbgal),                                                concatenate(mockspec, axis=1), mockpsfipara, mockpsfntype)
        mocktotlflux = retr_rofi_flux(globdata, mocknormback, mockpntsflux, globdata.fullindx)
        mocktotlcnts = mocktotlflux * globdata.expo * globdata.apix * globdata.diffener[:, None, None] # [1]

        if globdata.verbtype > 1:
            print 'mocknumbpnts: ', mocknumbpnts
            print 'mockfdfnnorm: ', mockfdfnnorm
            print 'mockfdfnslop: ', mockfdfnslop
            print 'mockpsfipara: '
            print mockpsfipara
            print 'mocknormback: '
            print mocknormback
            
            
        globdata.datacnts = zeros((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
        for i in globdata.indxener:
            for k in range(globdata.numbpixl):
                for m in globdata.indxevtt:
                    globdata.datacnts[i, k, m] = poisson(mocktotlcnts[i, k, m])
              
        if globdata.trueinfo:
            
            globdata.truelgal = []
            globdata.truebgal = []
            globdata.truespec = []
            if globdata.colrprio:
                globdata.truesind = []
            for i in globdata.indxpopl:
                globdata.truelgal.append(mocklgal[l])
                globdata.truebgal.append(mockbgal[l])
                globdata.truespec.append(mockspec[l])
                if globdata.colrprio:
                    globdata.truesind.append(mocksind[l])
                    
            globdata.indxtruepntstimevari = [array([])] * globdata.numbpopl
                    
            globdata.truenumbpnts = mocknumbpnts
            globdata.truefdfnnorm = mockfdfnnorm
            globdata.truefdfnslop = mockfdfnslop
            globdata.truenormback = mocknormback
            globdata.truecnts = mockcnts
            
            globdata.truespec = []
            for l in globdata.indxpopl:
                globdata.truespectemp = empty((3, globdata.numbener, mocknumbpnts[l]))
                globdata.truespectemp[:] = mockspec[l][None, :, :]
                globdata.truespec.append(globdata.truespectemp)
                
            globdata.truepsfipara = mockpsfipara
                

        #plot_pntsdiff()
    #plot_datacntshist()

    if amax(abs(globdata.datacnts - globdata.datacnts.astype(int)) / globdata.datacnts) > 1e-3:
        print 'Fractional counts!'
        
    if amin(globdata.datacnts) < 0.:
        print 'Negative counts!'

    globdata.expotemp = mean(globdata.expo, 1)
    globdata.minmcnts = globdata.minmspec * amin(globdata.expotemp, 1) * globdata.diffener
    globdata.maxmcnts = globdata.maxmspec * amax(globdata.expotemp, 1) * globdata.diffener
    globdata.binscnts = zeros((globdata.numbener, globdata.numbspec + 1))
    for i in globdata.indxener:
        globdata.binscnts[i, :] = logspace(log10(globdata.minmcnts[i]), log10(globdata.maxmcnts[i]), globdata.numbspec + 1) # [1]
        
    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            globdata.truelabl = 'Mock data'
        else:
            globdata.truelabl = '3FGL'
            
    ## Real data
    if globdata.datatype == 'inpt':
        if globdata.trueinfo:
            globdata.truenumbpnts = array([globdata.fgl3numbpnts], dtype=int)
            globdata.truelgal = [globdata.fgl3lgal]
            globdata.truebgal = [globdata.fgl3bgal]
            globdata.truespec = [globdata.fgl3spec]
            if globdata.colrprio:
                globdata.truesind = [globdata.fgl3sind]
            indxpixltemp = retr_indxpixl(globdata, globdata.truebgal[0], globdata.truelgal[0])
            globdata.truecnts = [globdata.fgl3spec[0, :, :, None] * globdata.expo[:, indxpixltemp, :] *                             globdata.diffener[:, None, None]]
            globdata.indxtruepntstimevari = [globdata.indxfgl3timevari]
            if globdata.exprtype == 'ferm':
                globdata.truespec = [globdata.fgl3spec]
                
    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            globdata.truepsfn = retr_psfn(globdata, globdata.truepsfipara, globdata.indxener,                                           globdata.angldisp, mockpsfntype)
        else:
            if globdata.exprtype == 'ferm':
                globdata.truepsfn = globdata.fermpsfn
            if globdata.exprtype == 'sdss':
                globdata.truepsfn = sdsspsfn
                
        truefwhm = retr_fwhm(globdata, globdata.truepsfn)
        
        truebackcnts = []
        globdata.truesigm = []
        for l in globdata.indxpopl:
            indxpixltemp = retr_indxpixl(globdata, globdata.truebgal[l], globdata.truelgal[l])
            truebackcntstemp = zeros((globdata.numbener, globdata.truenumbpnts[l], globdata.numbevtt))
            for c in globdata.indxback:
                truebackcntstemp += globdata.backflux[c][:, indxpixltemp, :] * globdata.expo[:, indxpixltemp, :] *                                 globdata.diffener[:, None, None] * pi * truefwhm[:, None, :]**2 / 4.
            truebackcnts.append(truebackcntstemp)
            globdata.truesigm.append(globdata.truecnts[l] / sqrt(truebackcntstemp))
        
        if globdata.verbtype > 1:
            print 'truelgal: ', globdata.truelgal
            print 'truebgal: ', globdata.truebgal
            print 'truespec: '
            print globdata.truespec
            print 'truenumbpnts: ', globdata.truenumbpnts
            print 'truefdfnnorm: ', globdata.truefdfnnorm
            print 'truefdfnslop: ', globdata.truefdfnslop
            print 'truenormback: '
            print globdata.truenormback
            print
        
    globdata.datafluxmean = mean(sum(globdata.datacnts, 2) / sum(globdata.expo, 2), 1) / globdata.apix / globdata.diffener
    globdata.datacntsmean = mean(sum(globdata.datacnts, 2), 1)
    globdata.datacntssatu = ceil((amax(sum(globdata.datacnts, 2), 1) - globdata.datacntsmean) * 0.05 + globdata.datacntsmean)
    globdata.resicntssatu = ceil(globdata.datacntssatu * 0.5)
    
    # auxiliary variables for plots
    if globdata.pixltype == 'heal':
        for i in globdata.indxener:
            for m in globdata.indxevtt:
                globdata.datacntscarttemp = tdpy.util.retr_cart(globdata.datacnts[i, :, m],                                                                 globdata.indxpixlrofi,                                                                 numbsideinpt=globdata.numbsideheal,                                                                 minmlgal=globdata.minmlgal,                                                                 maxmlgal=globdata.maxmlgal,                                                                 minmbgal=globdata.minmbgal,                                                                 maxmbgal=globdata.maxmbgal)
                if i == 0 and m == 0:
                    globdata.datacntscart = zeros((globdata.datacntscarttemp.shape[0], globdata.datacntscarttemp.shape[1], globdata.numbener, globdata.numbevtt))
                globdata.datacntscart[:, :, i, m] = globdata.datacntscarttemp
        
    else:
        globdata.datacntscart = globdata.datacnts.reshape((globdata.numbener, globdata.numbsidecart, globdata.numbsidecart, globdata.numbevtt))
        globdata.datacntscart = swapaxes(swapaxes(globdata.datacntscart, 0, 2), 0, 1)

    for i in globdata.indxener:
        indxdatacntscartsatu = where(globdata.datacntscart[:, :, i, :] > globdata.datacntssatu[i])
        globdata.datacntscart[indxdatacntscartsatu[0],                               indxdatacntscartsatu[1], i, indxdatacntscartsatu[2]] = globdata.datacntssatu[i]

    # make a look-up table of nearby pixels for each pixel
    path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%03d.p' % globdata.maxmgang
    if os.path.isfile(path):
        fobj = open(path, 'rb')
        globdata.indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        globdata.indxpixlprox = [[] for h in range(globdata.numbspecprox)]
        for j in globdata.indxpixl:
            dist = retr_dist(globdata, globdata.lgalgrid[j], globdata.bgalgrid[j], globdata.lgalgrid, globdata.bgalgrid)
            for h in range(globdata.numbspecprox):
                # temp
                globdata.indxpixlproxtemp = where(dist < deg2rad(globdata.maxmangleval))[0]
                globdata.indxpixlprox[h].append(globdata.indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(globdata.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()


    # temp
    #plot_3fgl_thrs()
    #plot_intr()
    #plot_king()
    #plot_look()
    
    if globdata.verbtype > 0:
        print 'Sampling...'
    
    timereal = zeros(globdata.numbproc)
    timeproc = zeros(globdata.numbproc)
    if globdata.numbproc == 1:
        gridchan = [work(globdata, 0)]
    else:
        if globdata.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        globdata.lock = mp.Manager().Lock()

        # process pool
        pool = mp.Pool(globdata.numbproc)
        
        # spawn the processes
        workpart = functools.partial(work, globdata)
        indxproc = arange(globdata.numbproc)
        gridchan = pool.map(workpart, indxproc)

        pool.close()
        pool.join()

    for k in range(globdata.numbproc):
        timereal[k] = gridchan[k][17]
        timeproc[k] = gridchan[k][18]


    if globdata.verbtype > 0:
        print 'Accumulating samples from all processes...'
        tim0 = time.time()

    # parse the sample bundle
    listsampvarb = zeros((globdata.numbsamp, globdata.numbproc, globdata.maxmsampsize))
    listindxprop = zeros((globdata.numbswep, globdata.numbproc))
    listchro = zeros((globdata.numbswep, globdata.numbproc, 4))
    listllik = zeros((globdata.numbsamp, globdata.numbproc))
    listlpri = zeros((globdata.numbsamp, globdata.numbproc))
    listaccp = zeros((globdata.numbswep, globdata.numbproc))
    listmodlcnts = zeros((globdata.numbsamp, globdata.numbproc, globdata.numbpixlsave))
    listpntsfluxmean = zeros((globdata.numbsamp, globdata.numbproc, globdata.numbener))
    listindxpntsfull = []
    listindxsampmodi = zeros((globdata.numbswep, globdata.numbproc), dtype=int)
    
    globdata.listauxipara = empty((globdata.numbswep, globdata.numbproc, globdata.numbcomp))
    globdata.listlaccfrac = empty((globdata.numbswep, globdata.numbproc))
    globdata.listnumbpair = empty((globdata.numbswep, globdata.numbproc))
    globdata.listjcbnfact = empty((globdata.numbswep, globdata.numbproc))
    globdata.listcombfact = empty((globdata.numbswep, globdata.numbproc))

    levi = 0.
    info = 0.
    
    for k in range(globdata.numbproc):
        globdata.rtag = retr_rtag(globdata, k)
        listchan = gridchan[k]
        listsampvarb[:, k, :] = listchan[0]
        listindxprop[:, k] = listchan[1]
        listchro[:, k, :] = listchan[2]
        listllik[:, k] = listchan[3]
        listlpri[:, k] = listchan[4]
        listaccp[:, k] = listchan[5]
        listmodlcnts[:, k, :] = listchan[6]    
        listindxpntsfull.append(listchan[7])
        listindxsampmodi[:, k] = listchan[8]
        globdata.listauxipara[:, k, :] = listchan[9]
        globdata.listlaccfrac[:, k] = listchan[10]
        globdata.listnumbpair[:, k] = listchan[11]
        globdata.listjcbnfact[:, k] = listchan[12]
        globdata.listcombfact[:, k] = listchan[13]
        levi += listchan[14]
        info += listchan[15]
        listpntsfluxmean[:, k, :] = listchan[16]

    listindxprop = listindxprop.flatten()
    globdata.listauxipara = globdata.listauxipara.reshape((globdata.numbswep * globdata.numbproc, globdata.numbcomp))
    globdata.listlaccfrac = globdata.listlaccfrac.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listnumbpair = globdata.listnumbpair.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listjcbnfact = globdata.listjcbnfact.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listcombfact = globdata.listcombfact.reshape(globdata.numbswep * globdata.numbproc)
    
    listchro = listchro.reshape((globdata.numbproc * globdata.numbswep, 4)) 
    listaccp = listaccp.flatten()
    listindxsampmodi = listindxsampmodi.flatten()
    
    globdata.rtag = retr_rtag(globdata, None)

    listnumbpnts = listsampvarb[:, :, globdata.indxsampnumbpnts].astype(int).reshape(globdata.numbsamp * globdata.numbproc, -1)
    listfdfnnorm = listsampvarb[:, :, globdata.indxsampfdfnnorm].reshape(globdata.numbsamp * globdata.numbproc, -1)
    listfdfnslop = listsampvarb[:, :, globdata.indxsampfdfnslop].reshape(globdata.numbsamp * globdata.numbproc, globdata.numbpopl, globdata.numbener)
    listpsfipara = listsampvarb[:, :, globdata.indxsamppsfipara].reshape(globdata.numbsamp * globdata.numbproc, -1)
    listnormback = listsampvarb[:, :, globdata.indxsampnormback].reshape(globdata.numbsamp * globdata.numbproc, globdata.numbback, globdata.numbener)
    
    listpntsfluxmean = listpntsfluxmean.reshape(globdata.numbsamp * globdata.numbproc, globdata.numbener)
    
    listlgal = [[] for l in globdata.indxpopl]
    listbgal = [[] for l in globdata.indxpopl]
    listspec = [[] for l in globdata.indxpopl]
    if globdata.colrprio:
        listsind = [[] for l in globdata.indxpopl]
    listspechist = empty((globdata.numbsamp * globdata.numbproc, globdata.numbpopl, globdata.numbener, globdata.numbspec))
    for k in range(globdata.numbproc):
        for j in range(globdata.numbsamp):            
            n = k * globdata.numbsamp + j
            indxsamplgal, indxsampbgal, indxsampspec,                 indxsampsind, indxsampcomp = retr_indx(globdata, listindxpntsfull[k][j])
            for l in globdata.indxpopl:
                listlgal[l].append(listsampvarb[j, k, indxsamplgal[l]])
                listbgal[l].append(listsampvarb[j, k, indxsampbgal[l]])
                listspec[l].append(listsampvarb[j, k, indxsampspec[l]])
                if globdata.colrprio:
                    listsind[l].append(listsampvarb[j, k, indxsampsind[l]])
                for i in globdata.indxener:
                    listspechist[n, l, i, :] = histogram(listspec[l][n][i, :], globdata.binsspec[i, :])[0]
        
    
    for l in globdata.indxpopl:
        
        listlgaltemp = zeros((globdata.numbsamp, globdata.maxmnumbpnts[l])) - 1.
        listbgaltemp = zeros((globdata.numbsamp, globdata.maxmnumbpnts[l])) - 1.
        listspectemp = zeros((globdata.numbsamp, globdata.numbener, globdata.maxmnumbpnts[l])) - 1.
        for k in range(globdata.numbsamp):
            listlgaltemp[k, 0:listlgal[l][k].size] = listlgal[l][k]
            listbgaltemp[k, 0:listbgal[l][k].size] = listbgal[l][k]
            listspectemp[k, :, 0:listspec[l][k].shape[1]] = listspec[l][k]

        listlgal[l] = listlgaltemp
        listbgal[l] = listbgaltemp 
        listspec[l] = listspectemp    

        
        
    if globdata.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        tim0 = time.time()

    # posterior maps
    pntsprob = zeros((globdata.numbpopl, globdata.numbener, globdata.numbpixl, globdata.numbspec))
    for k in range(globdata.numbsamp):
        for l in globdata.indxpopl:
            for i in globdata.indxener:
                for h in range(globdata.numbspec):
                    indxpnts = where((globdata.binsspec[i, h] < listspec[l][k, i, :]) & (listspec[l][k, i, :] < globdata.binsspec[i, h+1]))[0]
                    hpixl = retr_indxpixl(globdata, listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    
        
    if globdata.verbtype > 0:
        print 'Performing Gelman-Rubin convergence test...'
        tim0 = time.time()

    gmrbstat = zeros(globdata.numbpixlsave)
    for n in range(globdata.numbpixlsave):
        gmrbstat[n] = tdpy.mcmc.gmrb_test(listmodlcnts[:, :, n])


            
    pathprobcatl = os.environ["PCAT_DATA_PATH"] + '/probcatl_' + globdata.strgtime + '_' + globdata.rtag + '.fits'  
    
    head = pf.Header()
    head['numbener'] = (globdata.numbener, 'Number of energy bins')
    head['numbevtt'] = (globdata.numbevtt, 'Number of PSF class bins')
    head['numbpopl'] = (globdata.numbpopl, 'Number of PS population')
    #head['numbpsfipara'] = globdata.numbpsfipara
    #head['numbformpara'] = globdata.numbformpara
    
    head['numbsamp'] = globdata.numbsamp
    head['numbburn'] = globdata.numbburn
    head['numbswep'] = globdata.numbswep
    head['factthin'] = globdata.factthin
    head['numbpopl'] = globdata.numbpopl
    head['numbproc'] = globdata.numbproc
    
    head['maxmgang'] = globdata.maxmgang
    head['lgalcntr'] = globdata.lgalcntr
    head['bgalcntr'] = globdata.bgalcntr
    head['numbspec'] = globdata.numbspec
    #head['numbsideheal'] = globdata.numbsideheal
    #head['numbsidecart'] = globdata.numbsidecart
    
    head['minmlgal'] = globdata.minmlgal
    head['maxmlgal'] = globdata.maxmlgal
    head['minmbgal'] = globdata.minmbgal
    head['maxmbgal'] = globdata.maxmbgal
    
    head['datatype'] = globdata.datatype
    head['regitype'] = globdata.regitype
    head['psfntype'] = globdata.psfntype
    #if globdata.datatype == 'mock':
    #    head['mockpsfntype'] = mockpsfntype
    head['exprtype'] = globdata.exprtype
    head['pixltype'] = globdata.pixltype
    
    head['colrprio'] = globdata.colrprio
    head['trueinfo'] = globdata.trueinfo
    head['margsize'] = globdata.margsize
    head['strgtime'] = globdata.strgtime
    
    head['levi'] = levi
    head['info'] = info
    
    listhdun = []
    listhdun.append(pf.PrimaryHDU(header=head))
    for l in globdata.indxpopl:
        listhdun.append(pf.ImageHDU(listlgal[l]))
        listhdun[-1].header['EXTNAME'] = 'lgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listbgal[l]))
        listhdun[-1].header['EXTNAME'] = 'bgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listspec[l]))
        listhdun[-1].header['EXTNAME'] = 'specpopl%d' % l
        if globdata.colrprio:
            listhdun.append(pf.ImageHDU(listsind[l]))
            listhdun[-1].header['EXTNAME'] = 'sindpopl%d' % l
        listhdun.append(pf.ImageHDU(globdata.maxmnumbpnts[l]))
        listhdun[-1].header['EXTNAME'] = 'maxmnumbpntspopl%d' % l

    
    listhdun.append(pf.ImageHDU(listnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'numbpnts'

    listhdun.append(pf.ImageHDU(listfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'fdfnnorm'
    
    listhdun.append(pf.ImageHDU(listfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'fdfnslop'
    
    listhdun.append(pf.ImageHDU(listpsfipara))
    listhdun[-1].header['EXTNAME'] = 'psfipara'
    
    listhdun.append(pf.ImageHDU(listnormback))
    listhdun[-1].header['EXTNAME'] = 'normback'
    
    listhdun.append(pf.ImageHDU(listspechist))
    listhdun[-1].header['EXTNAME'] = 'spechist'
    
    
    listhdun.append(pf.ImageHDU(listllik))
    listhdun[-1].header['EXTNAME'] = 'llik'
    listhdun.append(pf.ImageHDU(listlpri))
    listhdun[-1].header['EXTNAME'] = 'lpri'
    
    # convergence diagnostics
    listhdun.append(pf.ImageHDU(gmrbstat))
    listhdun[-1].header['EXTNAME'] = 'gmrbstat'
    listhdun.append(pf.ImageHDU(listmodlcnts))
    listhdun[-1].header['EXTNAME'] = 'modlcnts'
    
    
    listhdun.append(pf.ImageHDU(pntsprob))
    listhdun[-1].header['EXTNAME'] = 'pntsprob'
    
    
    # truth information
    if globdata.trueinfo:
        listhdun.append(pf.ImageHDU(globdata.truenumbpnts))
        listhdun[-1].header['EXTNAME'] = 'truenumbpnts'

        for l in globdata.indxpopl:
            listhdun.append(pf.ImageHDU(globdata.truelgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truelgalpopl%d' % l

            listhdun.append(pf.ImageHDU(globdata.truebgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truebgalpopl%d' % l

            listhdun.append(pf.ImageHDU(globdata.truespec[l]))
            listhdun[-1].header['EXTNAME'] = 'truespecpopl%d' % l

        listhdun.append(pf.ImageHDU(globdata.truefdfnnorm))
        listhdun[-1].header['EXTNAME'] = 'truefdfnnorm'

        listhdun.append(pf.ImageHDU(globdata.truefdfnslop))
        listhdun[-1].header['EXTNAME'] = 'truefdfnslop'

        listhdun.append(pf.ImageHDU(globdata.truenormback))
        listhdun[-1].header['EXTNAME'] = 'truenormback'

        listhdun.append(pf.ImageHDU(globdata.truepsfipara))
        listhdun[-1].header['EXTNAME'] = 'truepsfipara'

        
    # boundaries
    listhdun.append(pf.ImageHDU(globdata.minmspec))
    listhdun[-1].header['EXTNAME'] = 'minmspec'
    listhdun.append(pf.ImageHDU(globdata.maxmspec))
    listhdun[-1].header['EXTNAME'] = 'maxmspec'
    listhdun.append(pf.ImageHDU(globdata.binsspec))
    listhdun[-1].header['EXTNAME'] = 'binsspec'
    listhdun.append(pf.ImageHDU(globdata.meanspec))
    listhdun[-1].header['EXTNAME'] = 'meanspec'

    if globdata.colrprio:
        listhdun.append(pf.ImageHDU(globdata.minmsind))
        listhdun[-1].header['EXTNAME'] = 'minmsind'
        listhdun.append(pf.ImageHDU(globdata.maxmsind))
        listhdun[-1].header['EXTNAME'] = 'maxmsind'
        listhdun.append(pf.ImageHDU(globdata.binssind))
        listhdun[-1].header['EXTNAME'] = 'binssind'
        listhdun.append(pf.ImageHDU(globdata.meansind))
        listhdun[-1].header['EXTNAME'] = 'meansind'
    
    listhdun.append(pf.ImageHDU(globdata.binsener))
    listhdun[-1].header['EXTNAME'] = 'binsener'
    listhdun.append(pf.ImageHDU(globdata.meanener))
    listhdun[-1].header['EXTNAME'] = 'meanener'
    
    
    listhdun.append(pf.ImageHDU(globdata.indxenerincl))
    listhdun[-1].header['EXTNAME'] = 'indxenerincl'
    
    listhdun.append(pf.ImageHDU(globdata.indxevttincl))
    listhdun[-1].header['EXTNAME'] = 'indxevttincl'
    
   
 
    # utilities
    listhdun.append(pf.ImageHDU(listpntsfluxmean))
    listhdun[-1].header['EXTNAME'] = 'listpntsfluxmean'
    
    listhdun.append(pf.ImageHDU(globdata.probprop))
    listhdun[-1].header['EXTNAME'] = 'probprop'

    listhdun.append(pf.ImageHDU(listindxprop))
    listhdun[-1].header['EXTNAME'] = 'indxprop'
    
    listhdun.append(pf.ImageHDU(listchro))
    listhdun[-1].header['EXTNAME'] = 'chro'
    
    listhdun.append(pf.ImageHDU(listaccp))
    listhdun[-1].header['EXTNAME'] = 'accp'
    
    listhdun.append(pf.ImageHDU(listindxsampmodi))
    listhdun[-1].header['EXTNAME'] = 'sampmodi'
    
    listhdun.append(pf.ImageHDU(globdata.listauxipara))
    listhdun[-1].header['EXTNAME'] = 'auxipara'
    
    listhdun.append(pf.ImageHDU(globdata.listlaccfrac))
    listhdun[-1].header['EXTNAME'] = 'laccfrac'
    
    listhdun.append(pf.ImageHDU(globdata.listnumbpair))
    listhdun[-1].header['EXTNAME'] = 'numbpair'
    
    listhdun.append(pf.ImageHDU(globdata.listjcbnfact))
    listhdun[-1].header['EXTNAME'] = 'jcbnfact'
    
    listhdun.append(pf.ImageHDU(globdata.listcombfact))
    listhdun[-1].header['EXTNAME'] = 'combfact'
    
    pf.HDUList(listhdun).writeto(pathprobcatl, clobber=True, output_verify='ignore')

    # temp
    print pf.info(pathprobcatl)
    
    if globdata.makeplot:
        plot_post(pathprobcatl)

    timetotlreal = time.time() - timetotlreal
    timetotlproc = time.clock() - timetotlproc
     
    if globdata.verbtype > 0:
        for k in range(globdata.numbproc):
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, timereal[k], timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (timetotlreal, sum(timeproc))
        print 'The ensemble of catalogs is at ' + pathprobcatl
        if globdata.makeplot:
            print 'The plots are in ' + globdata.plotpath
        
    return gridchan
    
    
def plot_samp(globdata):

    globdata.thisresicnts = globdata.datacnts - globdata.thismodlcnts

    thispsfn = retr_psfn(globdata, globdata.thissampvarb[globdata.indxsamppsfipara],                          globdata.indxener, globdata.angldisp, globdata.psfntype)
    if globdata.pixltype == 'cart':
        thispsfn = thispsfn.reshape((globdata.numbener, -1, globdata.numbevtt))
    thisfwhm = retr_fwhm(globdata, thispsfn)

    plot_psfn(globdata, thispsfn)

    globdata.thisbackcntsmean = empty((globdata.numbener, globdata.numbevtt))
    for c in globdata.indxback:
        globdata.thisbackcntsmean += mean(globdata.thissampvarb[globdata.indxsampnormback[c, :, None, None]] *                                           globdata.backflux[c] * globdata.expo *                                           globdata.diffener[:, None, None] * pi * thisfwhm[:, None, :]**2 / 4., 1)

    thiscnts = []
    for l in globdata.indxpopl:
        indxpixltemp = retr_indxpixl(globdata, globdata.thissampvarb[globdata.thisindxsampbgal[l]], globdata.thissampvarb[globdata.thisindxsamplgal[l]])
        cntstemp = globdata.thissampvarb[globdata.thisindxsampspec[l]][:, :, None] *             globdata.expo[:, indxpixltemp, :] * globdata.diffener[:, None, None]
        thiscnts.append(cntstemp)

        #if globdata.thissampvarb[globdata.indxsampnumbpnts[l]] > 1:
        if globdata.colrprio:
            plot_histsind(globdata, l)
        plot_scatpixl(globdata, l)
       
        if globdata.trueinfo:
            indxmodl, globdata.trueindxpntsbias, globdata.trueindxpntsmiss = pair_catl(globdata, l, globdata.thissampvarb[globdata.thisindxsamplgal[l]], \
                globdata.thissampvarb[globdata.thisindxsampbgal[l]], globdata.thissampvarb[globdata.thisindxsampspec[l]])

            thisspecmtch = globdata.thissampvarb[globdata.thisindxsampspec[l]][:, indxmodl]
            thisspecmtch[:, globdata.trueindxpntsmiss] = 0.
            if globdata.verbtype > 1:
                print 'thisspecmtch: (popl%d) ' % l
                print thisspecmtch
            plot_scatspec(globdata, l, thisspecmtch=thisspecmtch)
        plot_histspec(globdata, l)
        plot_histcnts(globdata, l, thiscnts)
        plot_compfrac(globdata)

    for i in globdata.indxener:
        
        plot_datacnts(globdata, i, None)
        #plot_catl(globdata, i, None, thiscnts)
        plot_modlcnts(globdata, i, None)
        plot_resicnts(globdata, i, None, globdata.thisresicnts)

        #for m in indxevtt:
        #    plot_datacnts(globdata, i, m)
        #    plot_catl(globdata, i, m, thiscnts)
        #    plot_modlcnts(globdata, i, m)
        #    plot_resicnts(globdata, i, m, globdata.thisresicnts)
        
    #if globdata.numbener == 3:
    #    plot_datacnts(globdata, None, None)
        
    #plot_fwhm(globdata, thisfwhm)
    #plot_backcntsmean(globdata, globdata.thisbackcntsmean)
    
    tempsampvarb, tempppixl, tempcnts, temppntsflux,         tempmodlflux, tempmodlcnts = pars_samp(globdata, globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])
    globdata.errrmodlcnts = globdata.thismodlcnts - tempmodlcnts
    
    for i in globdata.indxener:
        plot_errrcnts(globdata, i, None, globdata.errrmodlcnts)

    if amax(abs(globdata.errrmodlcnts)) > 0.1 and False:
        print 'Approximation error went above the limit!'
    
    
def rjmc(globdata, indxprocwork):

    # sweeps to be saved
    boolsave = zeros(globdata.numbswep, dtype=bool)
    indxswepsave = arange(globdata.numbburn, globdata.numbswep, globdata.factthin)
    boolsave[indxswepsave] = True
    
    sampindx = zeros(globdata.numbswep, dtype=int)
    sampindx[indxswepsave] = arange(globdata.numbsamp)

    listsampvarb = zeros((globdata.numbsamp, globdata.maxmsampsize)) + -1.
    listindxprop = zeros(globdata.numbswep)
    listchro = zeros((globdata.numbswep, 4))
    listllik = zeros(globdata.numbsamp)
    listlprising = zeros(globdata.numbsamp)
    listlpri = zeros(globdata.numbsamp)
    listaccp = zeros(globdata.numbswep, dtype=bool)
    listaccpspec = []
    listindxsampmodi = zeros(globdata.numbswep, dtype=int)
    listmodlcnts = zeros((globdata.numbsamp, globdata.numbpixlsave))
    listpntsfluxmean = zeros((globdata.numbsamp, globdata.numbener))
    listindxpntsfull = []
    
    globdata.listauxipara = zeros((globdata.numbswep, globdata.numbcomp))
    globdata.listlaccfrac = zeros(globdata.numbswep)
    globdata.listnumbpair = zeros(globdata.numbswep)
    globdata.listjcbnfact = zeros(globdata.numbswep)
    globdata.listcombfact = zeros(globdata.numbswep)

    # initialize the chain
    retr_llik(globdata, init=True)
    retr_lpri(globdata, init=True)

    # current sample index
    thiscntr = -1
    
    globdata.cntrswep = 0
    while globdata.cntrswep < globdata.numbswep:
        
        timeinit = time.time()
        
        if globdata.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % globdata.cntrswep

        thismakefram = (globdata.cntrswep % globdata.plotperd == 0) and indxprocwork == int(float(globdata.cntrswep) / globdata.numbswep * globdata.numbproc) and globdata.makeplot
        globdata.reje = False
    
        # choose a proposal type
        retr_thisindxprop(globdata, globdata.drmcsamp[:, 0])
            
        # save the proposal type
        listindxprop[globdata.cntrswep] = globdata.thisindxprop
        if globdata.verbtype > 1:
            print 'thisindxprop: ', globdata.strgprop[globdata.thisindxprop]
        
        if globdata.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop(globdata)
        timefinl = time.time()
        listchro[globdata.cntrswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if globdata.numbproc > 1:
                globdata.lock.acquire()
            print 'Process %d started making a frame' % indxprocwork
            plot_samp(globdata)
            print 'Process %d finished making a frame' % indxprocwork
            if globdata.numbproc > 1:
                globdata.lock.release()
            
        # reject the sample if proposal is outside the prior
        if globdata.thisindxprop != globdata.indxpropbrth and             globdata.thisindxprop != globdata.indxpropdeth and not globdata.reje:
            if where((globdata.drmcsamp[globdata.indxsampmodi, 1] < 0.) |                      (globdata.drmcsamp[globdata.indxsampmodi, 1] > 1.))[0].size > 0:
                globdata.reje = True
        if globdata.thisindxprop == globdata.indxproppsfipara:
            if globdata.psfntype == 'doubking':
                if globdata.nextpsfipara[1] > globdata.nextpsfipara[3]:
                    globdata.reje = True
            elif globdata.psfntype == 'doubgaus':
                if globdata.nextpsfipara[1] > globdata.nextpsfipara[2]:
                    globdata.reje = True
                
  
            
        if not globdata.reje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri(globdata)
            timefinl = time.time()
            listchro[globdata.cntrswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(globdata)          
            timefinl = time.time()
            listchro[globdata.cntrswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(globdata.deltllik + globdata.deltlpri + globdata.laccfrac)

            if globdata.verbtype > 1:
                print 'deltlpri'
                print globdata.deltlpri
                print 'deltllik'
                print globdata.deltllik
                print 'laccfrac'
                print globdata.laccfrac
                print
                
        else:
            accpprob = 0.
    
    
        # accept the sample
        if accpprob >= rand():

            if globdata.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(globdata)

            listaccp[globdata.cntrswep] = True

        # reject the sample
        else:

            if globdata.verbtype > 1:
                print 'Rejected.'

            listaccp[globdata.cntrswep] = False
             
        # sanity checks
        if where((globdata.drmcsamp[1:, 0] > 1.) | (globdata.drmcsamp[1:, 0] < 0.))[0].size > 0:
            print 'Unit sample vector went outside [0,1]!'
        for l in globdata.indxpopl:
            for i in globdata.indxenerfdfn:
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] < globdata.minmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went below the prior range!'
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] > globdata.maxmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went above the prior range!'          

        # save the sample
        if boolsave[globdata.cntrswep]:
            listsampvarb[sampindx[globdata.cntrswep], :] = globdata.thissampvarb
            listmodlcnts[sampindx[globdata.cntrswep], :] = globdata.thismodlcnts[0, globdata.indxpixlsave, 0]
            listpntsfluxmean[sampindx[globdata.cntrswep], :] = mean(sum(globdata.thispntsflux * globdata.expo, 2) / sum(globdata.expo, 2), 1)
            listindxpntsfull.append(globdata.thisindxpntsfull)
            listllik[sampindx[globdata.cntrswep]] = sum(globdata.thisllik)
            
            lpri = 0.
            for l in globdata.indxpopl:
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[l]]
                fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[l]]
                lpri += numbpnts * globdata.priofactlgalbgal + globdata.priofactfdfnslop + globdata.priofactfdfnnorm - log(fdfnnorm)
                for i in globdata.indxenerprio:
                    flux = globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]]
                    fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[l, i]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_spec(globdata, flux, fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])))
            listlpri[sampindx[globdata.cntrswep]] = lpri
            
            
            if globdata.tracsamp:
                
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[globdata.cntrswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            tim1 = time.time()
            listchro[globdata.cntrswep, 0] = tim1 - timeinit

        # log the progress
        if globdata.verbtype > 0:
            thiscntr = tdpy.util.show_prog(globdata.cntrswep, globdata.numbswep, thiscntr, indxprocwork=indxprocwork)
            
            
        if globdata.diagsamp:
            plot_datacnts(0, 0, nextstat=True)
            plot_resicnts(0, 0, thisresicnts, nextstat=True)
        
        if globdata.verbtype > 1:
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
            
        
        
        # update the sweep counter
        globdata.cntrswep += 1

    
    if globdata.verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    
    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(globdata.datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    minmlistllik = amin(listllik)
    levi = -log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    info = mean(listllik) - levi

    listchan = [listsampvarb, listindxprop, listchro, listllik, listlpri, listaccp, listmodlcnts, listindxpntsfull, listindxsampmodi, \
        globdata.listauxipara, globdata.listlaccfrac, globdata.listnumbpair, globdata.listjcbnfact, globdata.listcombfact, \
        levi, info, listpntsfluxmean]
    
    return listchan

