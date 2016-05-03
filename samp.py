
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
import tdpy_util.util

# pnts_tran
from pnts_tran.cnfg import *
from pnts_tran.main import *
from pnts_tran.samp import *
from pnts_tran.util import *
from pnts_tran.visu import *
from pnts_tran.plot import *



# In[ ]:

def rjmc(indxprocwork):

    # sweeps to be saved
    boolsave = zeros(globdata.numbswep, dtype=bool)
    indxswepsave = arange(globdata.numbburn, globdata.numbswep, globdata.factthin)
    boolsave[indxswepsave] = True
    
    sampindx = zeros(globdata.numbswep, dtype=int)
    sampindx[indxswepsave] = arange(globdata.numbsamp)

    listsampvarb = zeros((globdata.numbsamp, maxmsampsize)) + -1.
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
    
    global listauxipara, listlaccfrac, listcombfact, listjcbnfact, listnumbpair
    listauxipara = zeros((globdata.numbswep, globdata.numbcomp))
    listlaccfrac = zeros(globdata.numbswep)
    listnumbpair = zeros(globdata.numbswep)
    listjcbnfact = zeros(globdata.numbswep)
    listcombfact = zeros(globdata.numbswep)

    # initialize the chain
    retr_llik(init=True)
    retr_lpri(init=True)

    # current sample index
    thiscntr = -1
    
    globdata.thisindxswep = 0
    while globdata.thisindxswep < globdata.numbswep:
        
        timeinit = time.time()
        
        if verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % globdata.thisindxswep

        thismakefram = (globdata.thisindxswep % plotperd == 0) and             indxprocwork == int(float(globdata.thisindxswep) / globdata.numbswep * globdata.numbproc)             and globdata.makeplot
        globdata.reje = False
    
        # choose a proposal type
        retr_indxprop(globdata.drmcsamp[:, 0])
            
        # save the proposal type
        listindxprop[globdata.thisindxswep] = thisindxprop
        if verbtype > 1:
            print 'indxprop: ', strgprop[indxprop]
        
        
        if verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop()
        timefinl = time.time()
        listchro[globdata.thisindxswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if globdata.numbproc > 1:
                lock.acquire()
            print 'Process %d started making a frame' % indxprocwork
            plot_samp()
            print 'Process %d finished making a frame' % indxprocwork
            if globdata.numbproc > 1:
                lock.release()
            
        # reject the sample if proposal is outside the prior
        if thisindxprop != thisindxpropbrth and thisindxprop != thisindxpropdeth and not globdata.reje:
            if where((globdata.drmcsamp[indxsampmodi, 1] < 0.) | (globdata.drmcsamp[indxsampmodi, 1] > 1.))[0].size > 0:
                globdata.reje = True
        if thisindxprop == thisindxproppsfipara:
            if modlpsfntype == 'doubking':
                if nextpsfipara[1] > nextpsfipara[3]:
                    globdata.reje = True
            elif modlpsfntype == 'doubgaus':
                if nextpsfipara[1] > nextpsfipara[2]:
                    globdata.reje = True
                
  
            
        if not globdata.reje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri()
            timefinl = time.time()
            listchro[globdata.thisindxswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik()          
            timefinl = time.time()
            listchro[globdata.thisindxswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(deltllik + deltlpri + laccfrac)

            if verbtype > 1:
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

            if verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp()

            listaccp[globdata.thisindxswep] = True

        # reject the sample
        else:

            if verbtype > 1:
                print 'Rejected.'

            listaccp[globdata.thisindxswep] = False
             
        # sanity checks
        if where((globdata.drmcsamp[1:, 0] > 1.) | (globdata.drmcsamp[1:, 0] < 0.))[0].size > 0:
            print 'Unit sample vector went outside [0,1]!'
        for l in indxpopl:
            for i in indxenerfdfn:
                if where(thissampvarb[thisindxsampspec[l][i, :]] < minmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went below the prior range!'
                if where(thissampvarb[thisindxsampspec[l][i, :]] > maxmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went above the prior range!'          
        if amax(abs(errrmodlcnts)) > 0.1 and False:
            print 'Approximation error went above the limit!'

        # save the sample
        if boolsave[globdata.thisindxswep]:
            listsampvarb[sampindx[globdata.thisindxswep], :] = thissampvarb
            listmodlcnts[sampindx[globdata.thisindxswep], :] = thismodlcnts[0, gpixl, 0]
            listpntsfluxmean[sampindx[globdata.thisindxswep], :] = mean(sum(thispntsflux * expo, 2) / sum(expo, 2), 1)
            listindxpntsfull.append(thisindxpntsfull)
            listllik[sampindx[globdata.thisindxswep]] = sum(thisllik)
            
            lpri = 0.
            for l in indxpopl:
                numbpnts = thissampvarb[indxsampnumbpnts[l]]
                fdfnnorm = thissampvarb[indxsampfdfnnorm[l]]
                lpri += numbpnts * priofactlgalbgal + priofactfdfnslop + fdfnnormfact - log(fdfnnorm)
                for i in indxenerprio:
                    flux = thissampvarb[thisindxsampspec[l][i, :]]
                    fdfnslop = thissampvarb[indxsampfdfnslop[l, i]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_spec(flux, fdfnslop, minmspec[i], maxmspec[i])))
            listlpri[sampindx[globdata.thisindxswep]] = lpri
            
            
            if tracsamp:
                
                numbpnts = thissampvarb[indxsampnumbpnts[0]]
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
        if verbtype > 0:
            thiscntr = tdpy_util.util.show_prog(globdata.thisindxswep, globdata.numbswep,                                            thiscntr, indxprocwork=indxprocwork)
            
            
        if diagsamp:
            plot_datacnts(0, 0, nextstat=True)
            plot_resicnts(0, 0, thisresicnts, nextstat=True)
        
        if verbtype > 1:
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

    
    if verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    
    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    minmlistllik = amin(listllik)
    levi = -log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    info = mean(listllik) - levi

    listchan = [listsampvarb, listindxprop, listchro, listllik, listlpri, listaccp,                 listmodlcnts, listindxpntsfull, listindxsampmodi,                 listauxipara, listlaccfrac, listnumbpair, listjcbnfact, listcombfact, levi, info, listpntsfluxmean]
    
    
    return listchan



def retr_prop():

    # temp
    global thislgal, thisbgal, thisspec, nextlgal0, nextlgal1, nextbgal0, nextbgal1, nextspec0, nextspec1
    global thislgal0, thisbgal0, thisspec0, thislgal1, thisbgal1, thisspec1, nextlgal, nextbgal, nextspec
    global laccfrac, thisjcbnfact, thiscombfact
        
    global nextsampvarb, nextpsfipara, mpixl, indxpsfiparamodi,         modilgal, modibgal, modispec, modiflux, modisind,         indxinit0, indxinit1, nextflux0, nextflux1,         indxenermodi, indxbackmodi,         nextflux, nextsind, nextpntsflux, modiindxindxpnts, nextindxsampslot,         indxsampmodispec, indxsampmodi, indxcompmodi,         reje, mergindxpnts1
    
    global thisindxsamplgal, thisindxsampbgal, thisindxsampspec, thisindxsampcomp 
    thisindxsamplgal, thisindxsampbgal, thisindxsampspec, thisindxsampsind, thisindxsampcomp = retr_indx(thisindxpntsfull)
    
    if verbtype > 2:
        print 'retr_prop(): '

        print 'globdata.drmcsamp'
        print globdata.drmcsamp
        
        print 'thissampvarb: '
        for k in range(thissampvarb.size):
            if k == indxcompinit:
                print
            if k > indxcompinit and (k - indxcompinit) % globdata.numbcomp == 0:
                print
            print thissampvarb[k]
        print
            
        print 'thisindxpntsfull: ', thisindxpntsfull
        print 'thisindxpntsempt: ', thisindxpntsempt  
        print 'thisindxsamplgal: ', thisindxsamplgal
        print 'thisindxsampbgal: ', thisindxsampbgal
        print 'thisindxsampspec: '
        print thisindxsampspec
        if globdata.colrprio:
            print 'thisindxsampsind: ', thisindxsampsind
        print 'thisindxsampcomp: ', thisindxsampcomp
        print
        
    # hyper-parameter changes
    # mean number of point sources
    if thisindxprop == thisindxpropfdfnnorm:
        indxsampmodi = indxsampfdfnnorm[indxpoplmodi]
        retr_gaus(indxsampmodi, stdvfdfnnorm)
        nextsampvarb[indxsampfdfnnorm] = icdf_logt(globdata.drmcsamp[indxsampmodi, -1], minmfdfnnorm, factfdfnnorm)
        if globdata.colrprio:
            indxenermodi = indxenerfdfn
        else:
            indxenermodi = indxener
        
    # flux distribution function slope
    if thisindxprop == thisindxpropfdfnslop:
        if globdata.colrprio:
            indxenermodi = indxenerfdfn
        else:
            indxenermodi = choice(indxener)
        indxsampmodi = indxsampfdfnslop[indxpoplmodi, indxenermodi]
        retr_gaus(indxsampmodi, stdvfdfnslop)
        nextsampvarb[indxsampfdfnslop[indxpoplmodi, indxenermodi]] = icdf_atan(globdata.drmcsamp[indxsampmodi, -1], minmfdfnslop, factfdfnslop)
        if globdata.colrprio:
            indxsampmodi = concatenate((indxsampmodi, thisindxsampspec[indxpoplmodi][indxener, :].flatten()))
        else:
            indxsampmodi = concatenate((array([indxsampmodi]), thisindxsampspec[indxpoplmodi][indxenermodi, :]))
            
            
        if verbtype > 2:
            print 'indxpoplmodi'
            print indxpoplmodi
            print 'indxenermodi'
            print indxenermodi
            print 'nextsampvarb[indxsampfdfnslop]'
            print nextsampvarb[indxsampfdfnslop]
            print 'indxsampmodi'
            print indxsampmodi
        
        
            
    # PSF parameter change 
    if thisindxprop == thisindxproppsfipara:
        
        # index of the PSF parameter to change
        indxpsfiparamodi = choice(ipsfipara)

        # the energy bin of the PS flux map to be modified
        indxenermodi = array([(indxpsfiparamodi % numbpsfiparaevtt) // nformpara])
        indxsampmodi = indxsamppsfipara[indxpsfiparamodi]
        retr_gaus(indxsampmodi, stdvpsfipara)
        nextpsfipara = copy(thissampvarb[indxsamppsfipara])
        nextpsfipara[indxpsfiparamodi] = icdf_psfipara(globdata.drmcsamp[indxsampmodi, -1], indxpsfiparamodi)

        modilgal = thissampvarb[thisindxsamplgal[indxpoplmodi]]
        modibgal = thissampvarb[thisindxsampbgal[indxpoplmodi]]
        modispec = thissampvarb[thisindxsampspec[indxpoplmodi]]
        
        if verbtype > 1:
            
            print 'thissampvarb[indxsamppsfipara]: ', thissampvarb[indxsamppsfipara]
            print 'nextpsfipara: ', nextpsfipara
            print 'indxpsfiparamodi: ', indxpsfiparamodi
            print 'thissampvarb[indxsampmodi]: ', thissampvarb[indxsampmodi]
            print 'nextpsfipara: ', nextpsfipara[indxpsfiparamodi]
            print 

        
    # background changes
    
    # diffuse model
    if thisindxprop == thisindxpropnormback:

        # determine the sample index to be changed
        indxenermodi = choice(indxener)
        indxbackmodi = choice(indxback)
        indxsampmodi = indxsampnormback[indxbackmodi, indxenermodi]
        
        # propose
        retr_gaus(indxsampmodi, stdvback)

        # transform back from the unit space
        nextsampvarb[indxsampmodi] = icdf_logt(globdata.drmcsamp[indxsampmodi, -1],                                                minmnormback[indxbackmodi], factnormback[indxbackmodi])

        if verbtype > 1:
            print 'indxsampmodi: ', indxsampmodi
            print 'nextsampvarb[indxsampmodi]: ', nextsampvarb[indxsampmodi]

       
    
        
    # birth
    if thisindxprop == thisindxpropbrth:

        # change the number of PS
        nextsampvarb[indxsampnumbpnts[indxpoplmodi]] = thissampvarb[indxsampnumbpnts[indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxbrth = indxcompinit + maxmglobdata.numbcomp * indxpoplmodi + thisindxpntsempt[indxpoplmodi][0] * globdata.numbcomp
        
        # sample auxiliary variables
        if globdata.colrprio:
            numbauxipara = globdata.numbcompcolr
        else:
            numbauxipara = globdata.numbcomp
        auxipara = rand(numbauxipara)

        if globdata.colrprio:
            globdata.drmcsamp[indxbrth:indxbrth+2, -1] = auxipara[0:2]
            globdata.drmcsamp[indxbrth+2+indxenerfdfn, -1] = auxipara[-2]
            globdata.drmcsamp[indxbrth+globdata.numbcomp-1, -1] = auxipara[-1]
        else:
            globdata.drmcsamp[indxbrth:indxbrth+globdata.numbcomp, -1] = auxipara

        # sample indices to be modified
        indxsampmodi = arange(indxbrth, indxbrth + globdata.numbcomp, dtype=int)

        # modification catalog
        modilgal = empty(1)
        modibgal = empty(1)
        if globdata.colrprio:
            modiflux = empty(1)
            modisind = empty(1)
        modispec = zeros((globdata.numbener, 1))
        
        modilgal[0] = icdf_self(globdata.drmcsamp[indxbrth, -1], -maxmgangmarg, 2. * maxmgangmarg)
        modibgal[0] = icdf_self(globdata.drmcsamp[indxbrth+1, -1], -maxmgangmarg, 2. * maxmgangmarg)

        if globdata.colrprio:
            modiflux[0] = icdf_spec(globdata.drmcsamp[indxbrth+2+indxenerfdfn, -1],                                     thissampvarb[indxsampfdfnslop[indxpoplmodi, indxenerfdfn]],                                     minmspec[indxenerfdfn], maxmspec[indxenerfdfn])
            modisind[0] = icdf_atan(globdata.drmcsamp[indxbrth+globdata.numbcomp-1, -1], minmsind, factsind)
            modispec[:, 0] = retr_spec(modiflux[0], modisind[0]).flatten()
        else:
            for i in indxener:
                modispec[i, 0] = icdf_spec(globdata.drmcsamp[indxbrth+2+i, -1],                                            thissampvarb[indxsampfdfnslop[indxpoplmodi, i]],                                            minmspec[i], maxmspec[i])
    
        if verbtype > 1:
            print 'auxipara: ', auxipara
            print 'modilgal: ', modilgal
            print 'modibgal: ', modibgal
            print 'modispec: '
            print modispec
            print
            
    # kill
    if thisindxprop == thisindxpropdeth:
        
        # change the number of PS
        nextsampvarb[indxsampnumbpnts[indxpoplmodi]] = thissampvarb[indxsampnumbpnts[indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(thissampvarb[indxsampnumbpnts[indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        global killindxpnts
        killindxpnts = thisindxpntsfull[indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        indxsampmodi = array([])
            
        # modification catalog
        modilgal = empty(1)
        modibgal = empty(1)
        modispec = zeros((globdata.numbener, 1))
        modilgal[0] = thissampvarb[thisindxsamplgal[indxpoplmodi][killindxindxpnts]]
        modibgal[0] = thissampvarb[thisindxsampbgal[indxpoplmodi][killindxindxpnts]]
        modispec[:, 0] = -thissampvarb[thisindxsampspec[indxpoplmodi][:, killindxindxpnts]]

        if verbtype > 1:
            print 'killindxpnts: ', killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', modilgal
            print 'modibgal: ', modibgal
            print 'modispec: '
            print modispec
            print
            
  
    # split
    if thisindxprop == thisindxpropsplt:
        
        nextsampvarb[indxsampnumbpnts[indxpoplmodi]] = thissampvarb[indxsampnumbpnts[indxpoplmodi]] + 1
        
        # determine which point source to split
        #global spltindxindxpnts        
        spltindxindxpnts = choice(arange(thissampvarb[indxsampnumbpnts[indxpoplmodi]], dtype=int))
        spltindxpnts = thisindxpntsfull[indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        indxinit0 = indxcompinit + maxmglobdata.numbcomp * indxpoplmodi + thisindxpntsfull[indxpoplmodi][spltindxindxpnts] * globdata.numbcomp
        indxfinl0 = indxinit0 + globdata.numbcomp
        indxinit1 = indxcompinit + maxmglobdata.numbcomp * indxpoplmodi + thisindxpntsempt[indxpoplmodi][0] * globdata.numbcomp
        indxfinl1 = indxinit1 + globdata.numbcomp
        
        # determine the modified sample vector indices
        indxsampmodi = concatenate((arange(indxinit0, indxfinl0, dtype=int), arange(indxinit1, indxfinl1, dtype=int)))
        
        thislgal = thissampvarb[thisindxsamplgal[indxpoplmodi][spltindxindxpnts]]
        thisbgal = thissampvarb[thisindxsampbgal[indxpoplmodi][spltindxindxpnts]]
        thisspec = thissampvarb[thisindxsampspec[indxpoplmodi][:, spltindxindxpnts]]
        
        if verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxinit0: ', indxinit0
            print 'indxfinl0: ', indxfinl0
            print 'indxinit1: ', indxinit1
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
        auxipara[2:] = (exp(rand(globdata.numbener)) - 1.) / (exp(1.) - 1.) * (maxmspec - minmspec) + minmspec
        
        if verbtype > 1:
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
        
        if verbtype > 1:
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

            



        if abs(nextlgal0) > maxmgangmarg or abs(nextlgal1) > maxmgangmarg or         abs(nextbgal0) > maxmgangmarg or abs(nextbgal1) > maxmgangmarg or         where((nextspec0 > maxmspec) | (nextspec0 < minmspec))[0].size > 0 or         where((nextspec1 > maxmspec) | (nextspec1 < minmspec))[0].size > 0:
            globdata.reje = True
                
        if not globdata.reje:

            
            lgal = concatenate((array([nextlgal0, nextlgal1]), setdiff1d(thissampvarb[thisindxsamplgal[indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgal0, nextbgal1]), setdiff1d(thissampvarb[thisindxsampbgal[indxpoplmodi]], thisbgal)))
            pairlist = retr_pairlist(lgal, bgal)


            ## first new component
            globdata.drmcsamp[indxinit0, -1] = cdfn_self(nextlgal0, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[indxinit0+1, -1] = cdfn_self(nextbgal0, -maxmgangmarg, 2. * maxmgangmarg)
            for i in indxener:
                globdata.drmcsamp[indxinit0+2+i, -1] = cdfn_spec(nextspec0[i], thissampvarb[indxsampfdfnslop[indxpoplmodi, i]], minmspec[i], maxmspec[i])

            ## second new component
            globdata.drmcsamp[indxinit1, -1] = cdfn_self(nextlgal1, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[indxinit1+1, -1] = cdfn_self(nextbgal1, -maxmgangmarg, 2. * maxmgangmarg)
            for i in indxener:
                globdata.drmcsamp[indxinit1+2+i, -1] = cdfn_spec(nextspec1[i], thissampvarb[indxsampfdfnslop[indxpoplmodi, i]], minmspec[i], maxmspec[i])


            # construct the modification catalog
            modilgal = empty(3)
            modibgal = empty(3)
            modispec = empty((globdata.numbener, 3))

            ## component to be removed
            modilgal[0] = thislgal
            modibgal[0] = thisbgal
            modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            modilgal[1] = nextlgal0
            modibgal[1] = nextbgal0
            modispec[:, 1] = nextspec0.flatten()

            # second component to be added
            modilgal[2] = nextlgal1
            modibgal[2] = nextbgal1
            modispec[:, 2] = nextspec1.flatten()

        
    # merge
    if thisindxprop == thisindxpropmerg:
        
        nextsampvarb[indxsampnumbpnts[indxpoplmodi]] = thissampvarb[indxsampnumbpnts[indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([thissampvarb[thisindxsamplgal[indxpoplmodi]], thissampvarb[thisindxsampbgal[indxpoplmodi]]])
            
        lgal = thissampvarb[thisindxsamplgal[indxpoplmodi]]
        bgal = thissampvarb[thisindxsampbgal[indxpoplmodi]]
        pairlist = retr_pairlist(lgal, bgal)
        
        if verbtype > 1:
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
            mergindxpnts0 = thisindxpntsfull[indxpoplmodi][mergindxindxpnts0]
            mergindxsampinit0 = indxcompinit + mergindxpnts0 * globdata.numbcomp

            # second PS index to be merged
            mergindxpnts1 = thisindxpntsfull[indxpoplmodi][mergindxindxpnts1]
            mergindxsampinit1 = indxcompinit + mergindxpnts1 * globdata.numbcomp

            # determine the modified sample vector indices
            indxinit0 = indxcompinit + globdata.numbcomp * mergindxpnts0
            indxfinl0 = indxinit0 + globdata.numbcomp
            indxinit1 = indxcompinit + globdata.numbcomp * mergindxpnts1
            indxfinl1 = indxinit1 + globdata.numbcomp

            indxsampmodi = arange(indxinit0, indxfinl0)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxpnts0, mergindxpnts1], dtype=int))

            thislgal0 = thissampvarb[thisindxsamplgal[indxpoplmodi][mergindxindxpnts0]]
            thisbgal0 = thissampvarb[thisindxsampbgal[indxpoplmodi][mergindxindxpnts0]]
            thisspec0 = thissampvarb[thisindxsampspec[indxpoplmodi][:, mergindxindxpnts0]]

            thislgal1 = thissampvarb[thisindxsamplgal[indxpoplmodi][mergindxindxpnts1]]
            thisbgal1 = thissampvarb[thisindxsampbgal[indxpoplmodi][mergindxindxpnts1]]
            thisspec1 = thissampvarb[thisindxsampspec[indxpoplmodi][:, mergindxindxpnts1]]

            # auxiliary component
            auxipara = zeros(globdata.numbcomp)
            auxipara[0] = (thislgal0 - thislgal1) / 2.
            auxipara[1] = (thisbgal0 - thisbgal1) / 2.
            auxipara[2:] = thisspec0 - thisspec1

            # merged PS
            nextlgal = (thislgal0 + thislgal1) / 2.
            nextbgal = (thisbgal0 + thisbgal1) / 2.
            nextspec = thisspec0 + thisspec1
            
            globdata.drmcsamp[indxinit0, -1] = cdfn_self(nextlgal, -maxmgangmarg, 2. * maxmgangmarg)
            globdata.drmcsamp[indxinit0+1, -1] = cdfn_self(nextbgal, -maxmgangmarg, 2. * maxmgangmarg)
            for i in indxener:
                globdata.drmcsamp[indxinit0+2+i, -1] = cdfn_spec(nextspec[i], thissampvarb[indxsampfdfnslop[indxpoplmodi, i]], minmspec[i], maxmspec[i])

            # construct the modification catalog
            modilgal = empty(3)
            modibgal = empty(3)
            modispec = empty((globdata.numbener, 3))

            ## first component to be merged
            modilgal[0] = thislgal0
            modibgal[0] = thisbgal0
            modispec[:, 0] = -thisspec0.flatten()

            ## first component to be merged
            modilgal[1] = thislgal1
            modibgal[1] = thisbgal1
            modispec[:, 1] = -thisspec1.flatten()

            ## component to be added
            modilgal[2] = nextlgal
            modibgal[2] = nextbgal
            modispec[:, 2] = nextspec.flatten()

            if verbtype > 1:
                print 'mergindxpnts0: ', mergindxpnts0
                print 'mergindxindxpnts0: ', mergindxindxpnts0
                print 'mergindxpnts1: ', mergindxpnts1
                print 'mergindxindxpnts1: ', mergindxindxpnts1
                print 'indxinit0: ', indxinit0
                print 'indxfinl0: ', indxfinl0
                print 'indxinit1: ', indxinit1
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
    if thisindxprop >= thisindxproplgal:     
        
        if thisindxprop == thisindxproplgal or thisindxprop == thisindxpropbgal:
            if thisindxprop == thisindxproplgal:
                indxcompmodi = 0
            else:
                indxcompmodi = 1
            indxenermodi = indxener
        else:
            if globdata.colrprio:
                indxenermodi = indxener
                if thisindxprop == thisindxpropspec:
                    indxcompmodi = 2 + indxenerfdfn
                elif thisindxprop == thisindxpropsind:
                    indxcompmodi = array([2 + globdata.numbener])
            else:
                indxenermodi = choice(indxener)
                indxcompmodi = indxenermodi + 2
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(thissampvarb[indxsampnumbpnts[indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = thisindxpntsfull[indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        indxsampmodiinit = indxcompinit + maxmglobdata.numbcomp * indxpoplmodi + modiindxpnts * globdata.numbcomp
        
        # sample index to be modified
        indxsampmodi = indxsampmodiinit + indxcompmodi
        if globdata.colrprio:
            indxsampmodispec = indxsampmodiinit + 2 + indxener
        
        # propose
        if thisindxprop == thisindxpropspec:
            retr_gaus(indxsampmodi, stdvspec)
        else:
            retr_gaus(indxsampmodi, stdvlbhl) 

        # modification catalog
        modilgal = empty(2)
        modibgal = empty(2)
        modispec = empty((indxenermodi.size, 2))
  
        if globdata.colrprio:
            thisflux = thissampvarb[thisindxsampspec[indxpoplmodi][indxenerfdfn, modiindxindxpnts]]
            thissind = thissampvarb[thisindxsampsind[indxpoplmodi][modiindxindxpnts]]
            thisspec = retr_spec(thisflux, thissind)
        else:
            thisspec = thissampvarb[thisindxsampspec[indxpoplmodi][indxenermodi, modiindxindxpnts]]
            
        modispec[:, 0] = -thisspec.flatten()
        if thisindxprop == thisindxproplgal or thisindxprop == thisindxpropbgal:
            if indxcompmodi == 0:
                modilgal[0] = thissampvarb[thisindxsamplgal[indxpoplmodi][modiindxindxpnts]]
                modilgal[1] = icdf_self(globdata.drmcsamp[indxsampmodi, -1], -maxmgangmarg, 2. * maxmgangmarg)
                modibgal[:] = thissampvarb[thisindxsampbgal[indxpoplmodi][modiindxindxpnts]]
            else:
                modilgal[:] = thissampvarb[thisindxsamplgal[indxpoplmodi][modiindxindxpnts]]
                modibgal[0] = thissampvarb[thisindxsampbgal[indxpoplmodi][modiindxindxpnts]]
                modibgal[1] = icdf_self(globdata.drmcsamp[indxsampmodi, -1], -maxmgangmarg, 2. * maxmgangmarg)
            modispec[:, 1] = thisspec.flatten()
        else:
            modilgal[:] = thissampvarb[thisindxsamplgal[indxpoplmodi][modiindxindxpnts]]
            modibgal[:] = thissampvarb[thisindxsampbgal[indxpoplmodi][modiindxindxpnts]]
            if globdata.colrprio:
                if thisindxprop == thisindxpropspec:
                    modiflux = icdf_spec(globdata.drmcsamp[indxsampmodi, -1], thissampvarb[indxsampfdfnslop[indxpoplmodi, indxenerfdfn]], 
                                         minmspec[indxenerfdfn], maxmspec[indxenerfdfn])
                    modisind = thissampvarb[thisindxsampsind[indxpoplmodi][modiindxindxpnts]]
                else:
                    modiflux = thissampvarb[thisindxsampspec[indxpoplmodi][indxenerfdfn, modiindxindxpnts]]
                    modisind = icdf_atan(globdata.drmcsamp[indxsampmodi, -1], minmsind, factsind)

                modispec[:, 1] = retr_spec(modiflux, modisind).flatten()
            else:
                modispec[:, 1] = icdf_spec(globdata.drmcsamp[indxsampmodi, -1],                                            thissampvarb[indxsampfdfnslop[indxpoplmodi, indxenermodi]],                                            minmspec[indxenermodi], maxmspec[indxenermodi])

                
        # log
        if verbtype > 1:
            print 'modilgal: ', modilgal
            print 'modibgal: ', modibgal
            print 'modispec: '
            print modispec
            print 'indxcompmodi: ', indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts


    # energy bin in which to evaluate the log-likelihood
    if thisindxpropbrth <= thisindxprop <= thisindxpropmerg:
        indxenermodi = arange(globdata.numbener)

    if type(indxenermodi) == int64:
        indxenermodi = array([indxenermodi])

    if verbtype > 1:
        print 'indxsampmodi: ', indxsampmodi
        print 'indxenermodi: ', indxenermodi

    # auxiliary variable density fraction and jacobian
    global laccfrac
    if (indxprop == thisindxpropsplt or thisindxprop == thisindxpropmerg) and not globdata.reje:

        spltcombfact = log(thissampvarb[indxsampnumbpnts[indxpoplmodi]]**2 / len(pairlist))
        
        if thisindxprop == thisindxpropsplt:
            thiscombfact = spltcombfact 
            thisjcbnfact = spltjcbnfact
        else:
            thiscombfact = -spltcombfact 
            thisjcbnfact = -spltjcbnfact


        laccfrac = thisjcbnfact + thiscombfact

        listnumbpair[globdata.thisindxswep] = len(pairlist)
        listjcbnfact[globdata.thisindxswep] = thisjcbnfact
        listcombfact[globdata.thisindxswep] = thiscombfact
        listauxipara[globdata.thisindxswep, :] = auxipara
        listlaccfrac[globdata.thisindxswep] = laccfrac

    else:
        laccfrac = 0.  
        

