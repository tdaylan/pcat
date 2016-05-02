
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
import functools

# tdpy
import tdpy_util

# pnts_tran
from pnts_tran.cnfg import *
from pnts_tran.main import *
from pnts_tran.samp import *
from pnts_tran.util import *
from pnts_tran.visu import *
from pnts_tran.plot import *



# In[2]:

def work(globdata, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for the process
    seed()
    
    # construct the run tag
    globdata.rtag = retr_globdata.rtag(indxprocwork)
    
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
    globdata.thisindxsamplgal, globdata.thisindxsampbgal, globdata.thisindxsampspec,         globdata.thisindxsampsind, globdata.thisindxsampcomp = retr_indx(globdata.thisindxpntsfull)
      
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

    globdata.drmcsamp = zeros((globdata.maxmsampsize, 2))
    
    globdata.drmcsamp[globdata.indxsampnumbpnts, 0] = thisnumbpnts
    globdata.drmcsamp[globdata.indxsampfdfnnorm, 0] = rand(globdata.numbpopl)
    if globdata.trueinfo and globdata.datatype == 'mock':
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = cdfn_atan(mocksampvarb[globdata.indxsampfdfnslop], globdata.minmfdfnslop, globdata.factfdfnslop)
    else:
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = rand(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener))
    globdata.drmcsamp[globdata.indxsampnormback, 0] = rand(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener))
    if globdata.randinit or not globdata.trueinfo or globdata.truepsfipara == None:
        globdata.drmcsamp[globdata.indxsamppsfipara, 0] = rand(globdata.numbmodlpsfipara)
    else:
        for k in globdata.indxpsfipara:
            globdata.drmcsamp[globdata.indxsamppsfipara[k], 0] = cdfn_psfipara(globdata.truepsfipara[k], k)
        
    for l in globdata.indxpopl:
        if globdata.pntscntr:
            globdata.drmcsamp[globdata.thisindxsampcomp[l], 0] = 0.5
        else:
            if globdata.randinit or not globdata.trueinfo:
                globdata.drmcsamp[globdata.thisindxsampcomp[l], 0] = rand(globdata.thisindxsampcomp[l].size)
            else:
                globdata.drmcsamp[globdata.thisindxsamplgal[l], 0] = copy(cdfn_self(globdata.truelgal[l], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg))
                globdata.drmcsamp[globdata.thisindxsampbgal[l], 0] = copy(cdfn_self(globdata.truebgal[l], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg))  
                if globdata.datatype == 'mock':
                    globdata.drmcsamp[globdata.thisindxsampspec[l], 0] = copy(mocksamp[mockindxsampspec[l]])
                    if globdata.colrprio:
                        globdata.drmcsamp[globdata.thisindxsampsind[l], 0] = copy(mocksamp[mockindxsampsind[l]])
                else:
                    for i in globdata.indxenerfdfn:
                        fdfnsloptemp = icdf_atan(globdata.drmcsamp[globdata.indxsampfdfnslop[l, i], 0], globdata.minmfdfnslop, globdata.factfdfnslop)
                        globdata.drmcsamp[globdata.thisindxsampspec[l][i, :], 0] = copy(cdfn_spec(globdata.truespec[l][0, i, :],                                                                                 fdfnsloptemp, globdata.minmspec[i], globdata.maxmspec[i]))


    
    globdata.thissampvarb, thisppixl, thiscnts, globdata.thispntsflux,         thismodlflux, globdata.thismodlcnts = pars_samp(globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])



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
 
    listchan = rjmc(indxprocwork)
    
    timereal = time.time() - timereal
    timeproc = time.clock() - timeproc
    
    listchan.append(timereal)
    listchan.append(timeproc)
    
    return listchan


class globdatastrt(object):
    
    def __init__(self):
        pass


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
        globdata.mockpsfntype = cnfg['mockpsfntype']
    
    globdata.modlpsfntype = cnfg['modlpsfntype']
    
    globdata.liketype = cnfg['liketype']
    globdata.exprtype = cnfg['exprtype']
    globdata.pixltype = cnfg['pixltype']
    
    globdata.regitype = cnfg['regitype']
    globdata.randinit = cnfg['randinit']
    globdata.pntscntr = cnfg['pntscntr']

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
        
        globdata.mockfdfnslop = cnfg['mockfdfnslop']
        globdata.mocknumbpnts = cnfg['mocknumbpnts']
        globdata.mockfdfnnorm = cnfg['mockfdfnnorm']
        globdata.mocknormback = cnfg['mocknormback']

    globdata.maxmnumbpnts = cnfg['maxmnumbpnts']
    
    globdata.probprop = cnfg['probprop']
    
    exprfluxstrg = cnfg['exprfluxstrg']
    listbackfluxstrg = cnfg['listbackfluxstrg']
    expostrg = cnfg['expostrg']
    
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
    globdata.datetimestrg = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
    if globdata.verbtype > 0:
        print 'PNTS_TRAN started at ', globdata.datetimestrg
        print 'Initializing...'
    
    globdata.strgfluxunit = retr_strgfluxunit(globdata)
        
    # number of bins
    globdata.numbspec = 10
    
    globdata.numbbins = 10
    
    if globdata.datatype == 'mock':
        globdata.mocknormback = globdata.mocknormback[:, globdata.indxenerincl]
        if not globdata.colrprio:
            globdata.mockfdfnslop = globdata.mockfdfnslop[:, globdata.indxenerincl]

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
            
            globdata.mockfdfnslop = tile(globdata.mockfdfnslop, (1, globdata.numbener))
        
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
        globdata.fermpsfn = retr_fermpsfn(globdata)
        
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
    
    # PSF parameters
    if globdata.modlpsfntype == 'singgaus':
        globdata.numbmodlformpara = 1
    elif globdata.modlpsfntype == 'singking':
        globdata.numbmodlformpara = 2 
    elif globdata.modlpsfntype == 'doubgaus':
        globdata.numbmodlformpara = 3
    elif globdata.modlpsfntype == 'gausking':
        globdata.numbmodlformpara = 4
    elif globdata.modlpsfntype == 'doubking':
        globdata.numbmodlformpara = 5
        
    globdata.numbmodlpsfiparaevtt = globdata.numbener * globdata.numbmodlformpara
    globdata.numbmodlpsfipara = globdata.numbmodlpsfiparaevtt * globdata.numbevtt
    globdata.indxmodlpsfipara = arange(globdata.numbmodlpsfipara)   

    if globdata.datatype == 'mock':

        # mock PSF parameters
        if globdata.mockpsfntype == 'singgaus':
            globdata.numbmockformpara = 1
        elif globdata.mockpsfntype == 'singking':
            globdata.numbmockformpara = 2 
        elif globdata.mockpsfntype == 'doubgaus':
            globdata.numbmockformpara = 3
        elif globdata.mockpsfntype == 'gausking':
            globdata.numbmockformpara = 4
        elif globdata.mockpsfntype == 'doubking':
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
        
    # factors in the prior expression
    globdata.priofactlgalbgal = 2. * log(1. / 2. / globdata.maxmgang)
    globdata.priofactfdfnslop = globdata.numbener * log(1. / (arctan(globdata.maxmfdfnslop) - arctan(globdata.minmfdfnslop)))
    globdata.fdfnnormfact = log(1. / (log(globdata.maxmfdfnnorm) - log(globdata.minmfdfnnorm)))

    # sample vector indices  
    globdata.indxsampnumbpnts = arange(globdata.numbpopl)
    globdata.indxsampfdfnnorm = arange(globdata.numbpopl) + amax(globdata.indxsampnumbpnts) + 1
    globdata.indxsampfdfnslop = arange(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener)) + amax(globdata.indxsampfdfnnorm) + 1
    globdata.indxsamppsfipara = arange(globdata.numbmodlpsfipara) + amax(globdata.indxsampfdfnslop) + 1
    globdata.indxsampnormback = arange(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener)) + amax(globdata.indxsamppsfipara) + 1

    globdata.fluxpivt = sqrt(globdata.minmspec * globdata.maxmspec)
    
    globdata.maxmnumbcomp = globdata.maxmnumbpnts * globdata.numbcomp
    globdata.indxcompinit = amax(globdata.indxsampnormback) + 1
    
    # maximum number of parameters
    globdata.maxmsampsize = globdata.indxcompinit + globdata.maxmnumbcomp * globdata.numbpopl
    
    if globdata.numbburn == None:
        globdata.numbburn = globdata.numbswep / 5
    if globdata.factthin == None:
        globdata.factthin = min(int(globdata.maxmsampsize) * 5, globdata.numbswep / 2)

    globdata.strgprop = ['fdfnnorm', 'fdfnslop', 'psfipara', 'normback', 'brth',                 'deth', 'splt', 'merg', 'lgal', 'bgal', 'spec', 'sind']

    globdata.indxpropfdfnnorm = 0
    globdata.indxpropfdfnslop = 1
    globdata.indxproppsfipara = 2
    globdata.indxpropnormback = 3
    globdata.indxpropbrth = 4
    globdata.indxpropdeth = 5
    globdata.indxpropsplt = 6
    globdata.indxpropmerg = 7
    globdata.indxproplgal = 8
    globdata.indxpropbgal = 9
    globdata.indxpropspec = 10
    globdata.indxpropsind = 11

    if globdata.probprop == None:

        probfdfnnorm = array([1.])
        if globdata.colrprio:
            probfdfnslop = array([1.])
        else:
            #probfdfnslop = array([1.] * globdata.numbener)
            probfdfnslop = array([1.])
            
        #probpsfipara = array([1.] * globdata.numbmodlpsfipara)
        #probnormback = array([1.] * globdata.numbback * globdata.numbener)
        
        probpsfipara = array([1.])
        probnormback = array([1.])
        
        probbrth = array([0.1 * sum(globdata.maxmnumbpnts)])
        probdeth = array([0.1 * sum(globdata.maxmnumbpnts)])
        probsplt = array([0. * sum(globdata.maxmnumbpnts)])
        probmerg = array([0. * sum(globdata.maxmnumbpnts)])
        
        problgal = array([sum(globdata.maxmnumbpnts) / 2.])
        probbgal = array([sum(globdata.maxmnumbpnts) / 2.])
        if globdata.colrprio:
            probspec = array([sum(globdata.maxmnumbpnts) / 2.])
            probsind = array([sum(globdata.maxmnumbpnts) / 2.])
        else:
            probspec = array([sum(globdata.maxmnumbpnts) / 2.])
            #probspec = array([sum(globdata.maxmnumbpnts) / 2.] * globdata.numbener)
            probsind = array([0.])
            
                            
        globdata.probprop = concatenate((probfdfnnorm, probfdfnslop,                                 probpsfipara, probnormback,                                 probbrth, probdeth,                                 probsplt, probmerg,                                 problgal, probbgal,                                 probspec, probsind))
        

        globdata.probprop /= sum(globdata.probprop)
        
    # number of proposal types
    globdata.numbprop = globdata.probprop.size
   
    if globdata.verbtype > 1:
        print 'probprop: '
        print vstack((arange(globdata.numbprop), globdata.strgprop, globdata.probprop)).T


    # run tag
    globdata.rtag = retr_rtag(globdata, None)
    
    # plots
    if globdata.makeplot:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
            plotfold = '/n/pan/www/tansu/png/pnts_tran/'
        else:
            plotfold = os.environ["PNTS_TRAN_DATA_PATH"] + '/png/'
        globdata.plotpath = plotfold + globdata.datetimestrg + '_' + globdata.rtag + '/'
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
        
        path = os.environ["PNTS_TRAN_DATA_PATH"] + '/' + exprfluxstrg
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
    
    globdata.minmmodlpsfipara, globdata.maxmmodlpsfipara,         globdata.factmodlpsfipara, globdata.strgmodlpsfipara,         globdata.scalmodlpsfipara, globdata.indxmodlpsfipara = retr_psfimodl(globdata, globdata.modlpsfntype, 'modl')


 
    if globdata.verbtype > 1:
        print 'indxsampnumbpnts: ', globdata.indxsampnumbpnts
        print 'indxsampfdfnnorm: ', globdata.indxsampfdfnnorm
        print 'indxsampfdfnslop: ', globdata.indxsampfdfnslop
        print 'indxsamppsfipara: ', globdata.indxsamppsfipara
        print 'indxsampnormback: '
        print globdata.indxsampnormback
        print 'indxcompinit: ', globdata.indxcompinit

        
    # region of interest
    if globdata.pixltype == 'heal':
        
        lgalheal, bgalheal, globdata.numbsideheal, globdata.numbpixlheal, globdata.apix = tdpy_util.retr_healgrid(globdata.numbsideheal)

        globdata.indxpixlrofi = where((abs(lgalheal) < globdata.maxmgang) & (abs(bgalheal) < globdata.maxmgang))[0]
        globdata.indxpixlrofimarg = where((abs(lgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal) &                               (abs(bgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal))[0]

        
        globdata.lgalgrid = lgalheal[globdata.indxpixlrofi]
        globdata.bgalgrid = bgalheal[globdata.indxpixlrofi]
        
        path = os.environ["PNTS_TRAN_DATA_PATH"] + '/pixlcnvt_%03d.p' % globdata.maxmgang
        if os.path.isfile(path):
            fobj = open(path, 'rb')
            globdata.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            globdata.pixlcnvt = zeros(globdata.numbpixlheal, dtype=int)
            for k in range(globdata.indxpixlrofimarg.size):
                dist = retr_dist(lgalheal[globdata.indxpixlrofimarg[k]], bgalheal[globdata.indxpixlrofimarg[k]], globdata.lgalgrid, globdata.bgalgrid)
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
        
    # temp
    if globdata.exprtype == 'ferm':
        
        if psfntype == 'singgaus':
            jfermformpara = array([1])
        elif psfntype == 'singking':
            jfermformpara = array([1, 2])
        elif psfntype == 'doubgaus':
            jfermformpara = array([0, 1, 3])
        elif psfntype == 'gausking':
            jfermformpara = array([0, 1, 3, 4])
        elif psfntype == 'doubking':
            jfermformpara = arange(5)

        jfermpsfipara = tile(jfermformpara, globdata.numbener) + repeat(globdata.indxener, jfermformpara.size) * nfermformpara
        jfermpsfipara = tile(jfermpsfipara, globdata.numbevtt) + repeat(globdata.indxevtt, jfermpsfipara.size) * globdata.numbener * nfermformpara

        globdata.truepsfipara = fermpsfipara[jfermpsfipara]
        
    else:   
        
        globdata.truepsfipara = None

    # temp
    globdata.truepsfipara = array([1., deg2rad(0.1), deg2rad(2.), 10.,                           1., deg2rad(0.1), deg2rad(2.), 10.,                           1., deg2rad(0.1), deg2rad(2.), 10.,                           1., deg2rad(0.1), deg2rad(2.), 10.,                           1., deg2rad(0.1), deg2rad(2.), 10.])

    print 'truepsfipara'
    print globdata.truepsfipara
    
    # exposure
    if globdata.expostrg == 'unit':
        globdata.expo = ones((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    else:
        path = os.environ["PNTS_TRAN_DATA_PATH"] + '/' + globdata.expostrg
        globdata.expo = pf.getdata(path)

        if globdata.pixltype == 'heal':
            globdata.expo = globdata.expo[healindx]
        else:
            globdata.expo = globdata.expo.reshape((globdata.expo.shape[0], -1, globdata.expo.shape[-1]))
            
        globdata.expo = globdata.expo[globdata.filtindx]
    

    # backgrounds
    globdata.backflux = []
    globdata.backfluxmean = []
    for c, globdata.backfluxstrg in enumerate(listglobdata.backfluxstrg):
        path = os.environ["PNTS_TRAN_DATA_PATH"] + '/' + globdata.backfluxstrg
        globdata.backfluxtemp = pf.getdata(path)
        if globdata.pixltype == 'heal':
            globdata.backfluxtemp = globdata.backfluxtemp[healindx]
        else:
            globdata.backfluxtemp = globdata.backfluxtemp.reshape((globdata.backfluxtemp.shape[0], -1, globdata.backfluxtemp.shape[-1]))
        globdata.backfluxtemp = globdata.backfluxtemp[globdata.filtindx]
        globdata.backflux.append(globdata.backfluxtemp)
        globdata.backfluxmean.append(mean(sum(globdata.backfluxtemp * globdata.expo, 2) / sum(globdata.expo, 2), 1))
        

        
    # test plot
    # temp
    if globdata.datatype == 'inpt' and globdata.makeplot and False:
        for i in globdata.indxener:
            for m in globdata.indxevtt:
                globdata.backfluxtemp = zeros(globdata.numbpixl)
                for c in globdata.indxback:
                    globdata.backfluxtemp[:] += globdata.backflux[c][i, :, m]
                resicnts = (exprflux[i, :, m] - globdata.backfluxtemp) * globdata.expo[i, :, m] * globdata.apix * globdata.diffener[i]

                
                resicntstemp = tdpy_util.retr_cart(resicnts, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(resicntstemp, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testresiflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                cart = tdpy_util.retr_cart(exprflux[i, :, m], indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testexprflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                cart = tdpy_util.retr_cart(globdata.backfluxtemp, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testglobdata.backflux%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)

                cart = tdpy_util.retr_cart(globdata.expo[i, :, m], indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal,                                                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                                    minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                figr, axis = plt.subplots()
                imag = axis.imshow(cart, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
                cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
                plt.savefig(globdata.plotpath + 'testglobdata.expo%d%d_' % (i, m) + globdata.rtag + '.png')
                plt.close(figr)
                
                
                
    # get 3FGL catalog
    if globdata.exprtype == 'ferm':
        globdata.fgl3lgal, globdata.fgl3bgal, globdata.fgl3spec, fgl3gang,             globdata.fgl3cnts, globdata.fgl3timevari, globdata.fgl3sind, globdata.fgl3spectype, globdata.fgl3scur, globdata.fgl3scut = retr_fgl3()
        
        if globdata.regitype == 'ngal':
            rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
            globdata.fgl3bgal, globdata.fgl3lgal = rad2deg(rttr(deg2rad(90. - globdata.fgl3bgal), deg2rad(globdata.fgl3lgal)))
            globdata.fgl3bgal = 90. - globdata.fgl3bgal

        globdata.indxfgl3rofi = arange(globdata.fgl3lgal.size, dtype=int)
        for i in globdata.indxener:
            globdata.indxfgl3rofi = intersect1d(where((globdata.fgl3spec[0, i, :] > globdata.minmspec[i]) & (globdata.fgl3spec[0, i, :] < globdata.maxmspec[i]))[0], globdata.indxfgl3rofi)
        globdata.indxfgl3rofi = intersect1d(where((abs(globdata.fgl3lgal) < globdata.maxmgangmarg) & (abs(globdata.fgl3bgal) < globdata.maxmgangmarg))[0], globdata.indxfgl3rofi)

        globdata.indxglobdata.fgl3timevari = where(globdata.fgl3timevari > 72.44)[0]
        

        plot_fgl3()
        
        globdata.fgl3lgal = globdata.fgl3lgal[globdata.indxfgl3rofi]
        globdata.fgl3bgal = globdata.fgl3bgal[globdata.indxfgl3rofi]
        globdata.fgl3spec = globdata.fgl3spec[:, :, globdata.indxfgl3rofi]
        globdata.fgl3cnts = globdata.fgl3cnts[:, globdata.indxfgl3rofi, :]
        globdata.fgl3spectype = globdata.fgl3spectype[globdata.indxfgl3rofi]
        globdata.fgl3scur = globdata.fgl3scur[globdata.indxfgl3rofi]
        globdata.fgl3scut = globdata.fgl3scut[globdata.indxfgl3rofi]
        globdata.fgl3timevari = globdata.fgl3timevari[globdata.indxfgl3rofi]

        globdata.indxglobdata.fgl3timevari = where(globdata.fgl3timevari > 72.44)[0]

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

        if globdata.mocknumbpnts == None:
            globdata.mocknumbpnts = empty(globdata.numbpopl)
            for l in globdata.indxpopl:
                globdata.mocknumbpnts[l] = random_integers(globdata.minmnumbpnts, globdata.maxmnumbpnts[l])
        
        mockindxpntsfull = []    
        for l in globdata.indxpopl:
            mockindxpntsfull.append(range(globdata.mocknumbpnts[l]))
          
        mockindxsamplgal, mockindxsampbgal, mockindxsampspec, mockindxsampsind, mockindxsampcomp = retr_indx(mockindxpntsfull)

        mocksamp = zeros(globdata.maxmsampsize)
        mocksamp[globdata.indxsampnumbpnts] = globdata.mocknumbpnts

        if globdata.mockfdfnnorm != None:
            mocksamp[globdata.indxsampfdfnnorm] = cdfn_logt(globdata.mockfdfnnorm, globdata.minmfdfnnorm, globdata.factfdfnnorm)
        else:
            mocksamp[globdata.indxsampfdfnnorm] = rand(globdata.numbener)
        
            
        if globdata.mockfdfnslop != None:
            mocksamp[globdata.indxsampfdfnslop] = cdfn_atan(globdata.mockfdfnslop, globdata.minmfdfnslop, globdata.factfdfnslop)
        else:
            mocksamp[globdata.indxsampfdfnslop] = rand(globdata.numbener)
            
            
        if globdata.truepsfipara != None: 
            for k in globdata.indxpsfipara:
                if globdata.exprtype == 'ferm':
                    mocksamp[globdata.indxsamppsfipara[k]] = cdfn_psfipara(globdata.truepsfipara[k], k)
        else:
            mocksamp[globdata.indxsamppsfipara] = rand(globdata.numbmockpsfipara)
            
        for c in globdata.indxback:
            mocksamp[globdata.indxsampnormback[c, :]] = cdfn_logt(globdata.mocknormback[c, :], globdata.minmnormback[c], globdata.factnormback[c])

        if globdata.pntscntr:
            for l in globdata.indxpopl:
                mocksamp[mockindxsampcomp[l]] = 0.5
        else:
            for l in globdata.indxpopl:
                mocksamp[mockindxsampcomp[l]] = rand(mockindxsampcomp[l].size)

        if globdata.verbtype > 1:
            print 'mocksamp: '
            for k in range(mocksamp.size):
                print mocksamp[k]
            print

        mocksampvarb, mockppixl, mockcnts, mockpntsflux, mocktotlflux, mocktotlcnts = pars_samp(mockindxpntsfull, mocksamp)

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
                globdata.truelgal.append(mocksampvarb[mockindxsamplgal[l]])
                globdata.truebgal.append(mocksampvarb[mockindxsampbgal[l]])
                if globdata.colrprio:
                    globdata.truesind.append(mocksampvarb[mockindxsampsind[l]])
                    
            globdata.jtruepntstimevari = [array([])] * globdata.numbpopl
                    
            globdata.truenumbpnts = globdata.mocknumbpnts
            globdata.truefdfnnorm = globdata.mockfdfnnorm
            globdata.truefdfnslop = globdata.mockfdfnslop
            globdata.truenormback = globdata.mocknormback
            globdata.truecnts = mockcnts
            
            globdata.truespec = []
            for l in globdata.indxpopl:
                globdata.truespectemp = empty((3, globdata.numbener, globdata.mocknumbpnts[l]))
                globdata.truespectemp[:] = mocksampvarb[mockindxsampspec[l][None, :, :]]
                globdata.truespec.append(globdata.truespectemp)
                

        #plot_pntsdiff()
         
    #plot_globdata.datacntshist()
    

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
            globdata.truecnts = [globdata.fgl3spec[0, :, :, None] * globdata.expo[:, retr_pixl(globdata.truebgal[0], globdata.truelgal[0]), :] *                             globdata.diffener[:, None, None]]
            globdata.jtruepntstimevari = [globdata.indxglobdata.fgl3timevari]
            if globdata.exprtype == 'ferm':
                globdata.truespec = [globdata.fgl3spec]
                
    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            globdata.truepsfn = retr_psfn(globdata, globdata.truepsfipara, globdata.indxener,                                           globdata.angldisp, globdata.mockpsfntype, 'mock')
        else:
            if globdata.exprtype == 'ferm':
                globdata.truepsfn = globdata.fermpsfn
            if globdata.exprtype == 'sdss':
                globdata.truepsfn = sdsspsfn
                
        truefwhm = retr_fwhm(globdata.truepsfn)
        
        truebackcnts = []
        globdata.truesigm = []
        for l in globdata.indxpopl:
            ppixl = retr_pixl(globdata.truebgal[l], globdata.truelgal[l])
            truebackcntstemp = zeros((globdata.numbener, globdata.truenumbpnts[l], globdata.numbevtt))
            for c in globdata.indxback:
                truebackcntstemp += globdata.backflux[c][:, ppixl, :] * globdata.expo[:, ppixl, :] *                                 globdata.diffener[:, None, None] * pi * truefwhm[:, None, :]**2 / 4.
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
                globdata.datacntscarttemp = tdpy_util.retr_cart(globdata.datacnts[i, :, m],                                                                 globdata.indxpixlrofi,                                                                 numbsideinpt=globdata.numbsideheal,                                                                 minmlgal=globdata.minmlgal,                                                                 maxmlgal=globdata.maxmlgal,                                                                 minmbgal=globdata.minmbgal,                                                                 maxmbgal=globdata.maxmbgal)
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
    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/indxpixlprox_%03d.p' % globdata.maxmgang
    if os.path.isfile(path):
        fobj = open(path, 'rb')
        globdata.indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        globdata.indxpixlprox = [[] for h in range(globdata.numbspecprox)]
        for j in globdata.indxpixl:
            dist = retr_dist(globdata.lgalgrid[j], globdata.bgalgrid[j], globdata.lgalgrid, globdata.bgalgrid)
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
        gridchan = [work(0)]
        
    else:

        if globdata.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        globdata.lock = mp.Lock()

        # process pool
        pool = mp.Pool(globdata.numbproc)
        
        # spawn the processes
        workpart = functools.partial(work, globdata)
        gridchan = pool.map(workpart, range(globdata.numbproc))
        #gridchan = pool.map(work, range(globdata.numbproc))
        
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
    listmodlcnts = zeros((globdata.numbsamp, globdata.numbproc, ngpixl))
    listpntsfluxmean = zeros((globdata.numbsamp, globdata.numbproc, globdata.numbener))
    listindxpntsfull = []
    listindxsampmodi = zeros((globdata.numbswep, globdata.numbproc), dtype=int)
    
    listauxipara = empty((globdata.numbswep, globdata.numbproc, globdata.numbcomp))
    listlaccfrac = empty((globdata.numbswep, globdata.numbproc))
    listnumbpair = empty((globdata.numbswep, globdata.numbproc))
    listjcbnfact = empty((globdata.numbswep, globdata.numbproc))
    listcombfact = empty((globdata.numbswep, globdata.numbproc))

    levi = 0.
    info = 0.
    
    for k in range(globdata.numbproc):
        globdata.rtag = retr_globdata.rtag(k)
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
        listauxipara[:, k, :] = listchan[9]
        listlaccfrac[:, k] = listchan[10]
        listnumbpair[:, k] = listchan[11]
        listjcbnfact[:, k] = listchan[12]
        listcombfact[:, k] = listchan[13]
        levi += listchan[14]
        info += listchan[15]
        listpntsfluxmean[:, k, :] = listchan[16]

    listindxprop = listindxprop.flatten()
    listauxipara = listauxipara.reshape((globdata.numbswep * globdata.numbproc, globdata.numbcomp))
    listlaccfrac = listlaccfrac.reshape(globdata.numbswep * globdata.numbproc)
    listnumbpair = listnumbpair.reshape(globdata.numbswep * globdata.numbproc)
    listjcbnfact = listjcbnfact.reshape(globdata.numbswep * globdata.numbproc)
    listcombfact = listcombfact.reshape(globdata.numbswep * globdata.numbproc)
    
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
            indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(listindxpntsfull[k][j])
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
                    hpixl = retr_pixl(listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    
        
    if globdata.verbtype > 0:
        print 'Performing Gelman-Rubin convergence test...'
        tim0 = time.time()

    gmrbstat = zeros(ngpixl)
    for n in range(ngpixl):
        gmrbstat[n] = gmrb_test(listmodlcnts[:, :, n])


            
    pathprobcatl = os.environ["PNTS_TRAN_DATA_PATH"] + '/probcatl_' + globdata.datetimestrg + '_' + globdata.rtag + '.fits'  
    
    head = pf.Header()
    head['numbener'] = (globdata.numbener, 'Number of energy bins')
    head['numbevtt'] = (globdata.numbevtt, 'Number of PSF class bins')
    head['numbpopl'] = (globdata.numbpopl, 'Number of PS population')
    head['numbmodlpsfipara'] = globdata.numbmodlpsfipara
    head['numbmodlformpara'] = globdata.numbmodlformpara
    
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
    head['numbsideheal'] = globdata.numbsideheal
    head['numbsidecart'] = globdata.numbsidecart
    
    head['minmlgal'] = globdata.minmlgal
    head['maxmlgal'] = globdata.maxmlgal
    head['minmbgal'] = globdata.minmbgal
    head['maxmbgal'] = globdata.maxmbgal
    
    head['datatype'] = globdata.datatype
    head['regitype'] = globdata.regitype
    head['modlpsfntype'] = globdata.modlpsfntype
    if globdata.datatype == 'mock':
        head['mockpsfntype'] = globdata.mockpsfntype
    head['exprtype'] = globdata.exprtype
    head['pixltype'] = globdata.pixltype
    
    head['colrprio'] = globdata.colrprio
    head['trueinfo'] = globdata.trueinfo
    head['margsize'] = globdata.margsize
    head['datetimestrg'] = globdata.datetimestrg
    
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
        # temp
        if globdata.colrprio:
            listhdun.append(pf.ImageHDU(listsind[l]))
            listhdun[-1].header['EXTNAME'] = 'sindpopl%d' % l

    
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
    
    listhdun.append(pf.ImageHDU(listauxipara))
    listhdun[-1].header['EXTNAME'] = 'auxipara'
    
    listhdun.append(pf.ImageHDU(listlaccfrac))
    listhdun[-1].header['EXTNAME'] = 'laccfrac'
    
    listhdun.append(pf.ImageHDU(listnumbpair))
    listhdun[-1].header['EXTNAME'] = 'numbpair'
    
    listhdun.append(pf.ImageHDU(listjcbnfact))
    listhdun[-1].header['EXTNAME'] = 'jcbnfact'
    
    listhdun.append(pf.ImageHDU(listcombfact))
    listhdun[-1].header['EXTNAME'] = 'combfact'
    
    pf.HDUList(listhdun).writeto(pathprobcatl, clobber=True)
    
    if globdata.makeplot:
        plot_post(pathprobcatl)

    timetotlreal = time.time() - timetotlreal
    timetotlproc = time.clock() - timetotlproc
     
    if globdata.verbtype > 0:
        for k in range(globdata.numbproc):
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, timereal[k], timeproc[k])
        print 'PNTS_TRAN has run in %d real seconds, %d CPU seconds.' % (timetotlreal, sum(timeproc))
        print 'The ensemble of catalogs is at ' + pathprobcatl
        if globdata.makeplot:
            print 'The plots are in ' + globdata.plotpath
        
    return gridchan
    


# In[ ]:



