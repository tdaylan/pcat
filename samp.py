
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



# In[ ]:

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
    
    globdata.thisindxswep = 0
    while globdata.thisindxswep < globdata.numbswep:
        
        timeinit = time.time()
        
        if globdata.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % globdata.thisindxswep

        thismakefram = (globdata.thisindxswep % globdata.plotperd == 0) and             indxprocwork == int(float(globdata.thisindxswep) / globdata.numbswep * globdata.numbproc)             and globdata.makeplot
        globdata.reje = False
    
        # choose a proposal type
        retr_indxprop(globdata, globdata.drmcsamp[:, 0])
            
        # save the proposal type
        listindxprop[globdata.thisindxswep] = globdata.thisindxprop
        if globdata.verbtype > 1:
            print 'indxprop: ', strgprop[indxprop]
        
        
        if globdata.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop(globdata)
        timefinl = time.time()
        listchro[globdata.thisindxswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if globdata.numbproc > 1:
                lock.acquire()
            print 'Process %d started making a frame' % indxprocwork
            plot_samp(globdata)
            print 'Process %d finished making a frame' % indxprocwork
            if globdata.numbproc > 1:
                lock.release()
            
        # reject the sample if proposal is outside the prior
        if globdata.thisindxprop != globdata.indxpropbrth and             globdata.thisindxprop != globdata.indxpropdeth and not globdata.reje:
            if where((globdata.drmcsamp[globdata.indxsampmodi, 1] < 0.) |                      (globdata.drmcsamp[globdata.indxsampmodi, 1] > 1.))[0].size > 0:
                globdata.reje = True
        if globdata.thisindxprop == globdata.indxproppsfipara:
            if modlpsfntype == 'doubking':
                if globdata.nextpsfipara[1] > globdata.nextpsfipara[3]:
                    globdata.reje = True
            elif modlpsfntype == 'doubgaus':
                if globdata.nextpsfipara[1] > globdata.nextpsfipara[2]:
                    globdata.reje = True
                
  
            
        if not globdata.reje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri(globdata)
            timefinl = time.time()
            listchro[globdata.thisindxswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(globdata)          
            timefinl = time.time()
            listchro[globdata.thisindxswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(deltllik + deltlpri + laccfrac)

            if globdata.verbtype > 1:
                print 'deltlpri'
                print deltlpri
                print 'deltllik'
                print deltllik
                print 'laccfrac'
                print laccfrac
                print
                
        else:
            accpprob = 0.
    
    
        # accept the sample
        if accpprob >= rand():

            if globdata.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp()

            listaccp[globdata.thisindxswep] = True

        # reject the sample
        else:

            if globdata.verbtype > 1:
                print 'Rejected.'

            listaccp[globdata.thisindxswep] = False
             
        # sanity checks
        if where((globdata.drmcsamp[1:, 0] > 1.) | (globdata.drmcsamp[1:, 0] < 0.))[0].size > 0:
            print 'Unit sample vector went outside [0,1]!'
        for l in indxpopl:
            for i in globdata.indxenerfdfn:
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] < globdata.minmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went below the prior range!'
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] > globdata.maxmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went above the prior range!'          
        if amax(abs(errrmodlcnts)) > 0.1 and False:
            print 'Approximation error went above the limit!'

        # save the sample
        if boolsave[globdata.thisindxswep]:
            listsampvarb[sampindx[globdata.thisindxswep], :] = globdata.thissampvarb
            listmodlcnts[sampindx[globdata.thisindxswep], :] = thismodlcnts[0, gpixl, 0]
            listpntsfluxmean[sampindx[globdata.thisindxswep], :] = mean(sum(thispntsflux * expo, 2) / sum(expo, 2), 1)
            listindxpntsfull.append(globdata.thisindxpntsfull)
            listllik[sampindx[globdata.thisindxswep]] = sum(thisllik)
            
            lpri = 0.
            for l in indxpopl:
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[l]]
                fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[l]]
                lpri += numbpnts * priofactlgalbgal + priofactfdfnslop + fdfnnormfact - log(fdfnnorm)
                for i in globdata.indxenerprio:
                    flux = globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]]
                    fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[l, i]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_spec(flux, fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])))
            listlpri[sampindx[globdata.thisindxswep]] = lpri
            
            
            if tracsamp:
                
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[globdata.thisindxswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            tim1 = time.time()
            listchro[globdata.thisindxswep, 0] = tim1 - timeinit

        # log the progress
        if globdata.verbtype > 0:
            thiscntr = tdpy.util.show_prog(globdata.thisindxswep, globdata.numbswep,                                            thiscntr, indxprocwork=indxprocwork)
            
            
        if diagsamp:
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
        j += 1

    
    if globdata.verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    
    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    minmlistllik = amin(listllik)
    levi = -log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    info = mean(listllik) - levi

    listchan = [listsampvarb, listindxprop, listchro, listllik, listlpri, listaccp,                 listmodlcnts, listindxpntsfull, listindxsampmodi,                 globdata.listauxipara, globdata.listlaccfrac, globdata.listnumbpair,                 globdata.listjcbnfact, globdata.listcombfact, levi, info, listpntsfluxmean]
    
    return listchan



def retr_prop(globdata):
  
    globdata.thisindxsamplgal, globdata.thisindxsampbgal,         globdata.thisindxsampspec, globdata.thisindxsampsind,         globdata.thisindxsampcomp = retr_indx(globdata, globdata.thisindxpntsfull)
    
    if globdata.verbtype > 2:
        print 'retr_prop(): '

        print 'drmcsamp'
        print globdata.drmcsamp
        
        print 'globdata.thissampvarb: '
        for k in range(globdata.thissampvarb.size):
            if k == globdata.indxsampcompinit:
                print
            if k > globdata.indxsampcompinit and (k - globdata.indxsampcompinit) % globdata.numbcomp == 0:
                print
            print globdata.thissampvarb[k]
        print
            
        print 'thisindxpntsfull: ', globdata.thisindxpntsfull
        print 'thisindxpntsempt: ', thisindxpntsempt  
        print 'thisindxsamplgal: ', globdata.thisindxsamplgal
        print 'thisindxsampbgal: ', globdata.thisindxsampbgal
        print 'thisindxsampspec: '
        print globdata.thisindxsampspec
        if globdata.colrprio:
            print 'thisindxsampsind: ', globdata.thisindxsampsind
        print 'thisindxsampcomp: ', globdata.thisindxsampcomp
        print
        
    # hyper-parameter changes
    # mean number of point sources
    if globdata.thisindxprop == globdata.indxpropfdfnnorm:
        globdata.indxsampmodi = globdata.indxsampfdfnnorm[globdata.indxpoplmodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvfdfnnorm)
        globdata.nextsampvarb[globdata.indxsampfdfnnorm] = icdf_logt(globdata.drmcsamp[globdata.indxsampmodi, -1], minmfdfnnorm, factfdfnnorm)
        if globdata.colrprio:
            globdata.indxenermodi = globdata.indxenerfdfn
        else:
            globdata.indxenermodi = globdata.indxener
        
    # flux distribution function slope
    if globdata.thisindxprop == globdata.indxpropfdfnslop:
        if globdata.colrprio:
            globdata.indxenermodi = globdata.indxenerfdfn
        else:
            globdata.indxenermodi = choice(globdata.indxener)
        globdata.indxsampmodi = globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvfdfnslop)
        globdata.nextsampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]] =             icdf_atan(globdata.drmcsamp[globdata.indxsampmodi, -1], minmfdfnslop, factfdfnslop)
        if globdata.colrprio:
            globdata.indxsampmodi = concatenate((globdata.indxsampmodi,                                                  globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxener, :].flatten()))
        else:
            globdata.indxsampmodi = concatenate((array([globdata.indxsampmodi]),                                                  globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :]))
            
            
        if globdata.verbtype > 2:
            print 'globdata.indxpoplmodi'
            print globdata.indxpoplmodi
            print 'globdata.indxenermodi'
            print globdata.indxenermodi
            print 'nextsampvarb[globdata.indxsampfdfnslop]'
            print globdata.nextsampvarb[globdata.indxsampfdfnslop]
            print 'indxsampmodi'
            print globdata.indxsampmodi
        
        
            
    # PSF parameter change 
    if globdata.thisindxprop == globdata.indxproppsfipara:
        
        # index of the PSF parameter to change
        indxpsfiparamodi = choice(ipsfipara)

        # the energy bin of the PS flux map to be modified
        globdata.indxenermodi = array([(indxpsfiparamodi % numbpsfiparaevtt) // nformpara])
        globdata.indxsampmodi = globdata.indxsamppsfipara[indxpsfiparamodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvpsfipara)
        globdata.nextpsfipara = copy(globdata.thissampvarb[globdata.indxsamppsfipara])
        globdata.nextpsfipara[indxpsfiparamodi] = icdf_psfipara(globdata.drmcsamp[globdata.indxsampmodi, -1], indxpsfiparamodi)

        globdata.modilgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]]
        globdata.modibgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]
        globdata.modispec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi]]
        
        if globdata.verbtype > 1:
            
            print 'globdata.thissampvarb[globdata.indxsamppsfipara]: ', globdata.thissampvarb[globdata.indxsamppsfipara]
            print 'nextpsfipara: ', globdata.nextpsfipara
            print 'indxpsfiparamodi: ', indxpsfiparamodi
            print 'globdata.thissampvarb[globdata.indxsampmodi]: ', globdata.thissampvarb[globdata.indxsampmodi]
            print 'nextpsfipara: ', globdata.nextpsfipara[indxpsfiparamodi]
            print 

        
    # background changes
    
    # diffuse model
    if globdata.thisindxprop == globdata.indxpropnormback:

        # determine the sample index to be changed
        globdata.indxenermodi = choice(globdata.indxener)
        globdata.indxbackmodi = choice(globdata.indxback)
        globdata.indxsampmodi = globdata.indxsampnormback[globdata.indxbackmodi, globdata.indxenermodi]
        
        # propose
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvback)

        # transform back from the unit space
        globdata.nextsampvarb[globdata.indxsampmodi] = icdf_logt(globdata.drmcsamp[globdata.indxsampmodi, -1],                                                minmnormback[globdata.indxbackmodi], factnormback[globdata.indxbackmodi])

        if globdata.verbtype > 1:
            print 'indxsampmodi: ', globdata.indxsampmodi
            print 'nextsampvarb[globdata.indxsampmodi]: ', globdata.nextsampvarb[globdata.indxsampmodi]

    
    # birth
    if globdata.thisindxprop == globdata.indxpropbrth:

        # change the number of PS
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxbrth = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + thisindxpntsempt[globdata.indxpoplmodi][0] * globdata.numbcomp
        
        # sample auxiliary variables
        if globdata.colrprio:
            numbauxipara = globdata.numbcompcolr
        else:
            numbauxipara = globdata.numbcomp
        auxipara = rand(numbauxipara)

        if globdata.colrprio:
            globdata.drmcsamp[indxbrth:indxbrth+2, -1] = auxipara[0:2]
            globdata.drmcsamp[indxbrth+2+globdata.indxenerfdfn, -1] = auxipara[-2]
            globdata.drmcsamp[indxbrth+globdata.numbcomp-1, -1] = auxipara[-1]
        else:
            globdata.drmcsamp[indxbrth:indxbrth+globdata.numbcomp, -1] = auxipara

        # sample indices to be modified
        globdata.indxsampmodi = arange(indxbrth, indxbrth + globdata.numbcomp, dtype=int)

        # modification catalog
        globdata.modilgal = empty(1)
        globdata.modibgal = empty(1)
        if globdata.colrprio:
            globdata.modiflux = empty(1)
            globdata.modisind = empty(1)
        globdata.modispec = zeros((globdata.numbener, 1))
        
        globdata.modilgal[0] = icdf_self(globdata.drmcsamp[indxbrth, -1], -maxmgangmarg, 2. * maxmgangmarg)
        globdata.modibgal[0] = icdf_self(globdata.drmcsamp[indxbrth+1, -1], -maxmgangmarg, 2. * maxmgangmarg)

        if globdata.colrprio:
            globdata.modiflux[0] = icdf_spec(globdata, globdata.drmcsamp[indxbrth+2+globdata.indxenerfdfn, -1],                                     globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenerfdfn]],                                     globdata.minmspec[globdata.indxenerfdfn], globdata.maxmspec[globdata.indxenerfdfn])
            globdata.modisind[0] = icdf_atan(globdata.drmcsamp[indxbrth+globdata.numbcomp-1, -1], minmsind, factsind)
            globdata.modispec[:, 0] = retr_spec(globdata.modiflux[0], globdata.modisind[0]).flatten()
        else:
            for i in globdata.indxener:
                globdata.modispec[i, 0] = icdf_spec(globdata, globdata.drmcsamp[indxbrth+2+i, -1],                                            globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]],                                            globdata.minmspec[i], globdata.maxmspec[i])
    
        if globdata.verbtype > 1:
            print 'auxipara: ', auxipara
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print
            
    # kill
    if globdata.thisindxprop == globdata.indxpropdeth:
        
        # change the number of PS
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        global killindxpnts
        killindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        globdata.indxsampmodi = array([])
            
        # modification catalog
        globdata.modilgal = empty(1)
        globdata.modibgal = empty(1)
        globdata.modispec = zeros((globdata.numbener, 1))
        globdata.modilgal[0] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][killindxindxpnts]]
        globdata.modibgal[0] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][killindxindxpnts]]
        globdata.modispec[:, 0] = -globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, killindxindxpnts]]

        if globdata.verbtype > 1:
            print 'killindxpnts: ', killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print
            
  
    # split
    if globdata.thisindxprop == globdata.indxpropsplt:
        
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] + 1
        
        # determine which point source to split
        #global spltindxindxpnts        
        spltindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        spltindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        globdata.indxsampchd0 = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + globdata.thisindxpntsfull[globdata.indxpoplmodi][spltindxindxpnts] * globdata.numbcomp
        indxfinl0 = globdata.indxsampchd0 + globdata.numbcomp
        globdata.indxsampchd1 = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + thisindxpntsempt[globdata.indxpoplmodi][0] * globdata.numbcomp
        indxfinl1 = globdata.indxsampchd1 + globdata.numbcomp
        
        # determine the modified sample vector indices
        globdata.indxsampmodi = concatenate((arange(globdata.indxsampchd0, indxfinl0, dtype=int), arange(globdata.indxsampchd1, indxfinl1, dtype=int)))
        
        thislgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][spltindxindxpnts]]
        thisbgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][spltindxindxpnts]]
        thisspec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, spltindxindxpnts]]
        
        if globdata.verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxsampchd0: ', globdata.indxsampchd0
            print 'indxfinl0: ', indxfinl0
            print 'indxsampchd1: ', globdata.indxsampchd1
            print 'indxfinl1: ', indxfinl1
            if pixltype == 'heal':
                print 'thislgal: ', thislgal
                print 'thisbgal: ', thisbgal
            else:
                print 'thislgal: ', 3600. * thislgal
                print 'thisbgal: ', 3600. * thisbgal
            print 'thisspec: ', thisspec
            
            
        # determine the new components
        auxipara = empty(globdata.numbcomp)
        auxipara[0:2] = rand(2) * spmrlbhl
        auxipara[2:] = (exp(rand(globdata.numbener)) - 1.) / (exp(1.) - 1.) * (globdata.maxmspec - globdata.minmspec) + globdata.minmspec
        
        if globdata.verbtype > 1:
            if pixltype == 'heal':
                print 'auxipara[0]: ', auxipara[0]
                print 'auxipara[1]: ', auxipara[1]
            else:
                print 'auxipara[0]: ', 3600. * auxipara[0]
                print 'auxipara[1]: ', 3600. * auxipara[1]
            print 'auxipara[2:]: ', auxipara[2:]
            print
            
        nextlgal0 = thislgal + auxipara[0]
        nextlgal1 = thislgal - auxipara[0]
        nextbgal0 = thisbgal + auxipara[1]
        nextbgal1 = thisbgal - auxipara[1]
        nextspec0 = (thisspec + auxipara[2:]) / 2.
        nextspec1 = (thisspec - auxipara[2:]) / 2.
        
        if globdata.verbtype > 1:
            if pixltype == 'heal':
                print 'nextlgal0: ', nextlgal0
                print 'nextlgal1: ', nextlgal1
                print 'nextbgal0: ', nextbgal0
                print 'nextbgal1: ', nextbgal1
            else:
                print 'nextlgal0: ', 3600. * nextlgal0
                print 'nextlgal1: ', 3600. * nextlgal1
                print 'nextbgal0: ', 3600. * nextbgal0
                print 'nextbgal1: ', 3600. * nextbgal1
            print 'nextspec0: ', nextspec0
            print 'nextspec1: ', nextspec1

            



        if abs(nextlgal0) > maxmgangmarg or abs(nextlgal1) > maxmgangmarg or         abs(nextbgal0) > maxmgangmarg or abs(nextbgal1) > maxmgangmarg or         where((nextspec0 > globdata.maxmspec) | (nextspec0 < globdata.minmspec))[0].size > 0 or         where((nextspec1 > globdata.maxmspec) | (nextspec1 < globdata.minmspec))[0].size > 0:
            globdata.reje = True
                
        if not globdata.reje:

            
            lgal = concatenate((array([nextlgal0, nextlgal1]), setdiff1d(globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgal0, nextbgal1]), setdiff1d(globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]], thisbgal)))
            pairlist = retr_pairlist(lgal, bgal)


            ## first new component
            globdata.drmcsamp[globdata.indxsampchd0, -1] = cdfn_self(nextlgal0, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd0+1, -1] = cdfn_self(nextbgal0, -maxmgangmarg, 2. * maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd0+2+i, -1] = cdfn_spec(nextspec0[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])

            ## second new component
            globdata.drmcsamp[globdata.indxsampchd1, -1] = cdfn_self(nextlgal1, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd1+1, -1] = cdfn_self(nextbgal1, -maxmgangmarg, 2. * maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd1+2+i, -1] = cdfn_spec(nextspec1[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])


            # construct the modification catalog
            globdata.modilgal = empty(3)
            globdata.modibgal = empty(3)
            globdata.modispec = empty((globdata.numbener, 3))

            ## component to be removed
            globdata.modilgal[0] = thislgal
            globdata.modibgal[0] = thisbgal
            globdata.modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            globdata.modilgal[1] = nextlgal0
            globdata.modibgal[1] = nextbgal0
            globdata.modispec[:, 1] = nextspec0.flatten()

            # second component to be added
            globdata.modilgal[2] = nextlgal1
            globdata.modibgal[2] = nextbgal1
            globdata.modispec[:, 2] = nextspec1.flatten()

        
    # merge
    if globdata.thisindxprop == globdata.indxpropmerg:
        
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]], globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]])
            
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]
        pairlist = retr_pairlist(lgal, bgal)
        
        if globdata.verbtype > 1:
            print 'lgal'
            print lgal
            print 'bgal'
            print bgal
            print 'pairlist'
            print pairlist
            
            
        if len(pairlist) == 0:
            globdata.reje = True
        else:
            globdata.reje = False
            jpair = choice(arange(len(pairlist)))
            mergindxindxpnts0 = pairlist[jpair][0]
            mergindxindxpnts1 = pairlist[jpair][1]
  
        if not globdata.reje:

            # fisrt PS index to be merged
            mergindxchd0 = globdata.thisindxpntsfull[globdata.indxpoplmodi][mergindxindxpnts0]
            mergindxsampinit0 = globdata.indxsampcompinit + mergindxchd0 * globdata.numbcomp

            # second PS index to be merged
            globdata.mergindxchd1 = globdata.thisindxpntsfull[globdata.indxpoplmodi][mergindxindxpnts1]
            mergindxsampinit1 = globdata.indxsampcompinit + globdata.mergindxchd1 * globdata.numbcomp

            # determine the modified sample vector indices
            globdata.indxsampchd0 = globdata.indxsampcompinit + globdata.numbcomp * mergindxchd0
            indxfinl0 = globdata.indxsampchd0 + globdata.numbcomp
            globdata.indxsampchd1 = globdata.indxsampcompinit + globdata.numbcomp * globdata.mergindxchd1
            indxfinl1 = globdata.indxsampchd1 + globdata.numbcomp

            globdata.indxsampmodi = arange(globdata.indxsampchd0, indxfinl0)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxchd0, globdata.mergindxchd1], dtype=int))

            thislgal0 = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][mergindxindxpnts0]]
            thisbgal0 = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][mergindxindxpnts0]]
            thisspec0 = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, mergindxindxpnts0]]

            thislgal1 = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][mergindxindxpnts1]]
            thisbgal1 = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][mergindxindxpnts1]]
            thisspec1 = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, mergindxindxpnts1]]

            # auxiliary component
            auxipara = zeros(globdata.numbcomp)
            auxipara[0] = (thislgal0 - thislgal1) / 2.
            auxipara[1] = (thisbgal0 - thisbgal1) / 2.
            auxipara[2:] = thisspec0 - thisspec1

            # merged PS
            nextlgal = (thislgal0 + thislgal1) / 2.
            nextbgal = (thisbgal0 + thisbgal1) / 2.
            nextspec = thisspec0 + thisspec1
            
            globdata.drmcsamp[globdata.indxsampchd0, -1] = cdfn_self(nextlgal, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd0+1, -1] = cdfn_self(nextbgal, -maxmgangmarg, 2. * maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd0+2+i, -1] = cdfn_spec(nextspec[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])

            # construct the modification catalog
            globdata.modilgal = empty(3)
            globdata.modibgal = empty(3)
            globdata.modispec = empty((globdata.numbener, 3))

            ## first component to be merged
            globdata.modilgal[0] = thislgal0
            globdata.modibgal[0] = thisbgal0
            globdata.modispec[:, 0] = -thisspec0.flatten()

            ## first component to be merged
            globdata.modilgal[1] = thislgal1
            globdata.modibgal[1] = thisbgal1
            globdata.modispec[:, 1] = -thisspec1.flatten()

            ## component to be added
            globdata.modilgal[2] = nextlgal
            globdata.modibgal[2] = nextbgal
            globdata.modispec[:, 2] = nextspec.flatten()

            if globdata.verbtype > 1:
                print 'mergindxchd0: ', mergindxchd0
                print 'mergindxindxpnts0: ', mergindxindxpnts0
                print 'mergindxchd1: ', globdata.mergindxchd1
                print 'mergindxindxpnts1: ', mergindxindxpnts1
                print 'indxsampchd0: ', globdata.indxsampchd0
                print 'indxfinl0: ', indxfinl0
                print 'indxsampchd1: ', globdata.indxsampchd1
                print 'indxfinl1: ', indxfinl1
                if pixltype == 'heal':
                    print 'thislgal0: ', thislgal0
                    print 'thisbgal0: ', thisbgal0
                    print 'thislgal1: ', thislgal1
                    print 'thisbgal1: ', thisbgal1
                else:
                    print 'thislgal0: ', 3600. * thislgal0
                    print 'thisbgal0: ', 3600. * thisbgal0
                    print 'thislgal1: ', 3600. * thislgal1
                    print 'thisbgal1: ', 3600. * thisbgal1 
                print 'thisspec0: ', thisspec0
                print 'thisspec1: ', thisspec1

                if pixltype == 'heal':
                    print 'nextlgal: ', nextlgal
                    print 'nextbgal: ', nextbgal
                    print 'auxipara[0]: ', auxipara[0]
                    print 'auxipara[1]: ', auxipara[1]
                else:
                    print 'nextlgal: ', 3600. * nextlgal
                    print 'nextbgal: ', 3600. * nextbgal
                    print 'auxipara[0]: ', 3600. * auxipara[0]
                    print 'auxipara[1]: ', 3600. * auxipara[1]
                print 'nextspec: ', nextspec
                print 'auxipara[2:]: ', auxipara[2:]
                print

    # component change
    if globdata.thisindxprop >= globdata.indxproplgal:     
        
        if globdata.thisindxprop == globdata.indxproplgal or globdata.thisindxprop == globdata.indxpropbgal:
            if globdata.thisindxprop == globdata.indxproplgal:
                globdata.indxcompmodi = 0
            else:
                globdata.indxcompmodi = 1
            globdata.indxenermodi = globdata.indxener
        else:
            if globdata.colrprio:
                globdata.indxenermodi = globdata.indxener
                if globdata.thisindxprop == globdata.indxpropspec:
                    globdata.indxcompmodi = 2 + globdata.indxenerfdfn
                elif globdata.thisindxprop == globdata.indxpropsind:
                    globdata.indxcompmodi = array([2 + globdata.numbener])
            else:
                globdata.indxenermodi = choice(globdata.indxener)
                globdata.indxcompmodi = globdata.indxenermodi + 2
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        globdata.indxsampmodiinit = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + modiindxpnts * globdata.numbcomp
        
        # sample index to be modified
        globdata.indxsampmodi = globdata.indxsampmodiinit + globdata.indxcompmodi
        if globdata.colrprio:
            globdata.indxsampmodispec = globdata.indxsampmodiinit + 2 + globdata.indxener
        
        # propose
        if globdata.thisindxprop == globdata.indxpropspec:
            retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvspec)
        else:
            retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvlbhl) 

        # modification catalog
        globdata.modilgal = empty(2)
        globdata.modibgal = empty(2)
        globdata.modispec = empty((globdata.indxenermodi.size, 2))
  
        if globdata.colrprio:
            thisflux = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenerfdfn, modiindxindxpnts]]
            thissind = globdata.thissampvarb[globdata.thisindxsampsind[globdata.indxpoplmodi][modiindxindxpnts]]
            thisspec = retr_spec(thisflux, thissind)
        else:
            thisspec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, modiindxindxpnts]]
            
        globdata.modispec[:, 0] = -thisspec.flatten()
        if globdata.thisindxprop == globdata.indxproplgal or globdata.thisindxprop == globdata.indxpropbgal:
            if globdata.indxcompmodi == 0:
                globdata.modilgal[0] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modilgal[1] = icdf_self(globdata.drmcsamp[globdata.indxsampmodi, -1], -maxmgangmarg, 2. * maxmgangmarg)
                globdata.modibgal[:] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
            else:
                globdata.modilgal[:] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modibgal[0] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modibgal[1] = icdf_self(globdata.drmcsamp[globdata.indxsampmodi, -1], -maxmgangmarg, 2. * maxmgangmarg)
            globdata.modispec[:, 1] = thisspec.flatten()
        else:
            globdata.modilgal[:] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
            globdata.modibgal[:] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
            if globdata.colrprio:
                if globdata.thisindxprop == globdata.indxpropspec:
                    globdata.modiflux = icdf_spec(globdata, globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenerfdfn]], 
                                         globdata.minmspec[globdata.indxenerfdfn], globdata.maxmspec[globdata.indxenerfdfn])
                    globdata.modisind = globdata.thissampvarb[globdata.thisindxsampsind[globdata.indxpoplmodi][modiindxindxpnts]]
                else:
                    globdata.modiflux = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenerfdfn, modiindxindxpnts]]
                    globdata.modisind = icdf_atan(globdata.drmcsamp[globdata.indxsampmodi, -1], minmsind, factsind)

                globdata.modispec[:, 1] = retr_spec(globdata.modiflux, globdata.modisind).flatten()
            else:
                specunit = globdata.drmcsamp[globdata.indxsampmodi, -1]
                fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]]
                globdata.modispec[:, 1] = icdf_spec(globdata, specunit, fdfnslop,                                                     globdata.minmspec[globdata.indxenermodi],                                                     globdata.maxmspec[globdata.indxenermodi])

        # log
        if globdata.verbtype > 1:
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print 'indxcompmodi: ', globdata.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts


    # energy bin in which to evaluate the log-likelihood
    if globdata.indxpropbrth <= globdata.thisindxprop <= globdata.indxpropmerg:
        globdata.indxenermodi = arange(globdata.numbener)

    if type(globdata.indxenermodi) == int64:
        globdata.indxenermodi = array([globdata.indxenermodi])

    if globdata.verbtype > 1:
        print 'indxsampmodi: ', globdata.indxsampmodi
        print 'globdata.indxenermodi: ', globdata.indxenermodi

    # auxiliary variable density fraction and jacobian
    global laccfrac
    if (globdata.thisindxprop == globdata.indxpropsplt or globdata.thisindxprop == globdata.indxpropmerg) and not globdata.reje:

        spltcombfact = log(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]]**2 / len(pairlist))
        
        if globdata.thisindxprop == globdata.indxpropsplt:
            thiscombfact = spltcombfact 
            thisjcbnfact = spltjcbnfact
        else:
            thiscombfact = -spltcombfact 
            thisjcbnfact = -spltjcbnfact


        laccfrac = thisjcbnfact + thiscombfact

        globdata.listnumbpair[globdata.thisindxswep] = len(pairlist)
        globdata.listjcbnfact[globdata.thisindxswep] = thisjcbnfact
        globdata.listcombfact[globdata.thisindxswep] = thiscombfact
        globdata.listauxipara[globdata.thisindxswep, :] = auxipara
        globdata.listlaccfrac[globdata.thisindxswep] = laccfrac

    else:
        laccfrac = 0.  
        

