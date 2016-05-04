
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
import tdpy_util.util

# pnts_tran
from pnts_tran.cnfg import *
from pnts_tran.main import *
from pnts_tran.samp import *
from pnts_tran.util import *
from pnts_tran.visu import *
from pnts_tran.plot import *



# In[ ]:

def retr_fwhm(globdata, psfn):
    
    fwhm = zeros((globdata.numbener, globdata.numbevtt))
    for i in globdata.indxener:
        for m in globdata.indxevtt:
            indxangldispgood = argsort(psfn[i, :, m])
            intpfwhm = max(amax(psfn[i, :, m]) / 2., amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, indxangldispgood, m]) and intpfwhm < amax(psfn[i, indxangldispgood, m]):
                fwhm[i, m] = 2. * interp1d(psfn[i, indxangldispgood, m], globdata.angldisp[indxangldispgood])(intpfwhm) # [rad]
    return fwhm


def retr_spec(globdata, flux, sind):

    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])
        
    spec = flux[None, :] * (globdata.meanener[:, None] / globdata.enerfdfn)**(-sind[None, :])
    
    return spec


def retr_indx(globdata, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampsind = []
    indxsampcomp = []
    for l in globdata.indxpopl:
        indxsamplgaltemp = globdata.indxcompinit + globdata.maxmnumbcomp * l +             array(indxpntsfull[l], dtype=int) * globdata.numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], globdata.numbener, 0) +                          repeat(arange(globdata.numbener), len(indxpntsfull[l])).reshape(globdata.numbener, -1))
        if globdata.colrprio:
            indxsampsind.append(indxsamplgaltemp + 2 + globdata.numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, globdata.numbcomp) + tile(arange(globdata.numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(globdata, lgal, bgal, spec, psfipara):
    
    numbpnts = lgal.size
    
    dist = empty((globdata.numbpixl, numbpnts))
    for k in range(numbpnts):
        dist[:, k] = retr_dist(globdata, lgal[k], bgal[k], globdata.lgalgrid, globdata.bgalgrid)

    # convolve with the PSF
    pntsflux = empty((numbpnts, globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    for k in range(numbpnts):
        psfn = retr_psfn(globdata, psfipara, globdata.indxener, dist[:, k], globdata.psfntype, 'modl')
        pntsflux[k, :, :, :] = spec[:, k, None, None] * psfn

    # sum contributions from all PS
    pntsfluxtemp = sum(pntsflux, 0) 

    return pntsfluxtemp


def retr_rofi_flux(globdata, normback, pntsflux, tempindx):

    modlflux = pntsflux[tempindx]
    for c in indxback:
        modlflux += normback[c, :, None, None] * globdata.backflux[c][tempindx]        
    
    return modlflux


def cdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))

    if flux <= fluxbrek:
        
        fluxunit = norm / (1. - fdfnsloplowr) *             ((flux / fluxbrek)**(1. - fdfnsloplowr) - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
        
    else:
        
        fluxunit = norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -             norm / (1. - fdfnslopuppr) * (1. - (flux / fluxbrek)**(1. - fdfnslopuppr))
       
    return fluxunit


def pdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))

    if flux <= fluxbrek:
        
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnsloplowr)
        
    else:
        
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnslopuppr)
        
    return pdfnflux


def icdf_spec_brok(globdata, fluxunit, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):
    
    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) -                  (1. - (globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    
    fluxunitbrek = norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
    
    if fluxunit < fluxunitbrek:
        
        flux = fluxbrek * (fluxunit * (1. - fdfnsloplowr) / norm +                        (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))**(1. / (1. - fdfnsloplowr))
        
    else:
        
        flux = fluxbrek * (1. - (norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -                                  fluxunit) * (1. - fdfnslopuppr) / norm)**(1. / (1. - fdfnslopuppr))

    return flux


def cdfn_spec(globdata, flux, fdfnslop, minmspectemp, maxmspectemp):
        
    fluxunit = (flux**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) / (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop))
        
    return fluxunit


def icdf_spec(globdata, fluxunit, fdfnslop, minmspectemp, maxmspectemp):
    
    flux = (fluxunit * (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) + minmspectemp**(1. - fdfnslop))**(1. / (1. - fdfnslop))
    
    return flux


def pdfn_spec(globdata, flux, fdfnslop, minmspectemp, maxmspectemp):
  
    pdfnflux = (1. - fdfnslop) / (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) * flux**(-fdfnslop)
          
    return pdfnflux


def icdf_self(paraunit, minmpara, factpara):
    para = factpara * paraunit + minmpara
    return para


def cdfn_self(para, minmpara, factpara):
    paraunit = (para - minmpara) / factpara
    return paraunit


def icdf_logt(paraunit, minmpara, factpara):
    para = minmpara * exp(paraunit * factpara)
    return para


def cdfn_logt(para, minmpara, factpara):
    paraunit = log(para / minmpara) / factpara
    return paraunit


def icdf_atan(paraunit, minmpara, factpara):
    para = tan(factpara * paraunit + arctan(minmpara))
    return para


def cdfn_atan(para, minmpara, factpara):
    paraunit = (arctan(para) - arctan(minmpara)) / factpara
    return paraunit


def icdf_psfipara(globdata, psfiparaunit, thisindxpsfipara):

    minmpsfipara = globdata.minmpsfipara[thisindxpsfipara]
    factpsfipara = globdata.factpsfipara[thisindxpsfipara]
    scalpsfipara = globdata.scalpsfipara[thisindxpsfipara]
        
    if scalpsfipara == 'self':
        psfipara = icdf_self(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfipara = icdf_logt(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfipara = icdf_atan(psfiparaunit, minmpsfipara, factpsfipara)

    return psfipara


def cdfn_psfipara(globdata, psfipara, thisindxpsfipara):
    
    minmpsfipara = globdata.minmpsfipara[thisindxpsfipara]
    factpsfipara = globdata.factpsfipara[thisindxpsfipara]
    scalpsfipara = globdata.scalpsfipara[thisindxpsfipara]

    if scalpsfipara == 'self':
        psfiparaunit = cdfn_self(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfiparaunit = cdfn_logt(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfiparaunit = cdfn_atan(psfipara, minmpsfipara, factpsfipara)
    
    return psfiparaunit


def retr_indxprop(globdata, samp):
    
    # choose the population to be modified
    globdata.indxpoplmodi = choice(globdata.indxpopl)
    
    if thissampvarb[indxsampnumbpnts[globdata.indxpoplmodi]] == maxmnumbpnts[globdata.indxpoplmodi]:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probpropmaxm)
    elif thissampvarb[indxsampnumbpnts[globdata.indxpoplmodi]] == minmnumbpnts:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probpropminm)
    else:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probprop)
        

def retr_indxpixl(globdata, bgal, lgal):

    if globdata.pixltype == 'heal':
        indxpixl = globdata.pixlcnvt[ang2pix(globdata.numbsideheal, deg2rad(90. - bgal), deg2rad(lgal))]
    else:
        
        indxlgcr = floor(globdata.numbsidecart * (lgal - globdata.minmlgal) / 2. / globdata.maxmgang).astype(int)
        indxbgcr = floor(globdata.numbsidecart * (bgal - globdata.minmbgal) / 2. / globdata.maxmgang).astype(int)

        if isscalar(indxlgcr):
            if indxlgcr < 0:
                indxlgcr = 0
            if indxlgcr >= globdata.numbsidecart:
                indxlgcr = globdata.numbsidecart - 1
        else:
            indxlgcr[where(indxlgcr < 0)] = 0
            indxlgcr[where(indxlgcr >= globdata.numbsidecart)] = globdata.numbsidecart - 1
            
        if isscalar(indxbgcr):
            if indxbgcr < 0:
                indxbgcr = 0
            if indxbgcr >= globdata.numbsidecart:
                indxbgcr = globdata.numbsidecart - 1
        else:
            indxbgcr[where(indxbgcr < 0)] = 0
            indxbgcr[where(indxbgcr >= globdata.numbsidecart)] = globdata.numbsidecart - 1
            
        indxpixl = indxbgcr * globdata.numbsidecart + indxlgcr

    return indxpixl


def retr_llik(globdata, init=False):

    global thisllik, nextmodlflux, nextmodlcnts, deltllik, nextllik, thismodlcnts, modiindx
    
    if init:

        thisllik = datacnts * log(thismodlcnts) - thismodlcnts
        
    elif globdata.thisindxprop >= indxproppsfipara:

        # determine pixels over which to evaluate the log-likelihood
        if globdata.thisindxprop == indxpropnormback:
            indxpixlmodi = ipixl
            
        if globdata.thisindxprop == indxproppsfipara or globdata.thisindxprop >= indxpropbrth:
            
            if globdata.thisindxprop == indxproppsfipara:
                numbpnts = int(sum(thissampvarb[indxsampnumbpnts]))
                lgal = thissampvarb[concatenate(thisindxsamplgal)]
                bgal = thissampvarb[concatenate(thisindxsampbgal)]
                spec = thissampvarb[concatenate(thisindxsampspec)[globdata.indxenermodi, :]]
            else:
                numbpnts = modilgal.shape[0]
                lgal = modilgal
                bgal = modibgal
                spec = modispec
                
            thisindxpixlprox = []
            for k in range(numbpnts):
                indxspecproxtemp = argmin(spec[0, k] - meanspecprox)
                indxpixltemp = retr_indxpixl(globdata, bgal[k], lgal[k])
                thisindxpixlprox.append(globdata.indxpixlprox[indxspecproxtemp][indxpixltemp])
            indxpixlmodi = unique(concatenate(thisindxpixlprox))

            #print 'thisindxpixlprox[0].size: '
            #print thisindxpixlprox[0].size
            #print 'indxpixlmodi.size'
            #print indxpixlmodi.size
            #print
            
        # construct the mesh grid for likelihood evaluation
        if globdata.thisindxprop >= indxproppsfipara:
            modiindx = meshgrid(globdata.indxenermodi, indxpixlmodi, globdata.indxevtt, indexing='ij')

        # update the point source flux map
        if globdata.thisindxprop == indxproppsfipara or globdata.thisindxprop >= indxpropbrth:

            if globdata.thisindxprop == indxproppsfipara:
                nextpntsflux[modiindx] = 0.
            else:
                nextpntsflux[modiindx] = thispntsflux[modiindx]
                
            for k in range(numbpnts):
                
                # calculate the distance to the pixels to be updated
                dist = retr_dist(globdata, lgal[k], bgal[k], globdata.lgalgrid[thisindxpixlprox[k]],                                      globdata.bgalgrid[thisindxpixlprox[k]])

                # evaluate the PSF over the set of data cubes to be updated
                if globdata.thisindxprop == indxproppsfipara:
                    temppsfipara = nextpsfipara
                else:
                    temppsfipara = thissampvarb[indxsamppsfipara]
                psfn = retr_psfn(globdata, temppsfipara, globdata.indxenermodi, dist, globdata.psfntype, 'modl')
                                
                # update the data cubes
                for i in range(globdata.indxenermodi.size):
                    nextpntsflux[globdata.indxenermodi[i], thisindxpixlprox[k], :] += spec[i, k] * psfn[i, :, :]
            
        if globdata.thisindxprop == indxpropnormback:
            
            normbacktemp = empty((numbback, 1))
            for c in indxback:
                if c == indxbackmodi:
                    normbacktemp[c, 0] = nextsampvarb[indxsampnormback[c, globdata.indxenermodi]]
                else:
                    normbacktemp[c, 0] = thissampvarb[indxsampnormback[c, globdata.indxenermodi]]
                
            nextmodlflux[modiindx] = retr_rofi_flux(globdata, normbacktemp, thispntsflux, modiindx)

        if globdata.thisindxprop == indxproppsfipara or globdata.thisindxprop >= indxpropbrth:
            normback = thissampvarb[indxsampnormback[meshgrid(indxback, globdata.indxenermodi, indexing='ij')]]
            nextmodlflux[modiindx] = retr_rofi_flux(globdata, normback, nextpntsflux, modiindx)

        nextmodlcnts[modiindx] = nextmodlflux[modiindx] * globdata.expo[modiindx] *             apix * globdata.diffener[globdata.indxenermodi, None, None] # [1]
        nextllik[modiindx] = datacnts[modiindx] * log(nextmodlcnts[modiindx]) - nextmodlcnts[modiindx]
            
        if not isfinite(nextllik[modiindx]).any():
            warnings.warn('Log-likelihood went NAN!')
            
        deltllik = sum(nextllik[modiindx] - thisllik[modiindx])
        
    else:
        
        deltllik = 0.
        
    
def retr_fdfn(globdata, fdfnnorm, fdfnslop, i):
               
    fluxhistmodl = fdfnnorm / diffspec[i, indxspecpivt] * diffspec[i, :] * (meanspec[i, :] / fluxpivt[i])**(-fdfnslop)
          
    return fluxhistmodl


def retr_lpri(globdata, init=False):
        
    global nextlpri, thislpri, deltlpri

    if init:
        thislpri = zeros((numbpopl, globdata.numbener))
        
        for i in globdata.indxenerfdfn:
            for l in globdata.indxpopl:
                fluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[l]], thissampvarb[indxsampfdfnslop[l, i]], i)
                fluxhist = histogram(thissampvarb[thisindxsampspec[l][i, :]], binsspec[i, :])[0]
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                thislpri[l, i] = sum(lprbpois)            
      
        nextlpri = copy(thislpri)
                
    else:
        nextlpri = copy(thislpri)
        if globdata.thisindxprop == indxpropfdfnnorm or globdata.thisindxprop == indxpropfdfnslop             or globdata.thisindxprop >= indxpropbrth and globdata.thisindxprop <= indxpropmerg:
              
            if globdata.colrprio:
                indxenertemp = globdata.indxenerfdfn
            else:
                indxenertemp = globdata.indxenermodi
            for i in indxenertemp:
                if globdata.thisindxprop == indxpropfdfnnorm:
                    fdfnnorm = nextsampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = thissampvarb[indxsampfdfnslop[globdata.indxpoplmodi, i]]
                elif globdata.thisindxprop == indxpropfdfnslop:
                    fdfnnorm = thissampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = nextsampvarb[indxsampfdfnslop[globdata.indxpoplmodi, i]]
                else:
                    fdfnnorm = thissampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = thissampvarb[indxsampfdfnslop[globdata.indxpoplmodi, i]]
                fluxhistmodl = retr_fdfn(fdfnnorm, fdfnslop, i)
                              
                fluxhist = histogram(thissampvarb[thisindxsampspec[globdata.indxpoplmodi][i, :]], binsspec[i, :])[0] 
                if globdata.thisindxprop == indxpropbrth:
                    fluxhist += histogram(modispec[i, 0], binsspec[i, :])[0]
                elif globdata.thisindxprop == indxpropdeth:
                    fluxhist -= histogram(-modispec[i, 0], binsspec[i, :])[0]
                elif globdata.thisindxprop == indxpropsplt:
                    fluxhist -= histogram(-modispec[i, 0], binsspec[i, :])[0]
                    fluxhist += histogram(modispec[i, 1:3], binsspec[i, :])[0]
                elif globdata.thisindxprop == indxpropmerg:
                    fluxhist -= histogram(-modispec[i, 0:2], binsspec[i, :])[0]
                    fluxhist += histogram(modispec[i, 2], binsspec[i, :])[0]
                
                if False:
                    
                    prevfluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]], thissampvarb[indxsampfdfnslop[globdata.indxpoplmodi, i]], i)
                    get_ipython().magic(u'matplotlib inline')
                    plt.loglog(meanspec[i, :], fluxhistmodl)
                    plt.loglog(meanspec[i, :], fluxhist)
                    plt.loglog(meanspec[i, :], prevfluxhistmodl, ls='--')
                    plt.show()


                
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                nextlpri[globdata.indxpoplmodi, i] = sum(lprbpois)


                
            
            deltlpri = sum(nextlpri[globdata.indxpoplmodi, globdata.indxenermodi] - thislpri[globdata.indxpoplmodi, globdata.indxenermodi])

                      
        else:
            deltlpri = 0.

        
def pars_samp(globdata, indxpntsfull, samp):
    

    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[indxsampnumbpnts] = samp[indxsampnumbpnts]
    sampvarb[indxsampfdfnnorm] = icdf_logt(samp[indxsampfdfnnorm], minmfdfnnorm, factfdfnnorm)
    sampvarb[indxsampfdfnslop] = icdf_atan(samp[indxsampfdfnslop], minmfdfnslop, factfdfnslop)
         
    for c in indxback:
        sampvarb[indxsampnormback[c, :]] = icdf_logt(samp[indxsampnormback[c, :]], minmnormback[c], factnormback[c])
    for k in ipsfipara:
        sampvarb[indxsamppsfipara[k]] = icdf_psfipara(globdata, samp[indxsamppsfipara[k]], k, 'modl')

    cnts = []
    listspectemp = []
    indxpixlpnts = []
    for l in globdata.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg) 
        for i in globdata.indxenerfdfn:
            sampvarb[indxsampspec[l][i, :]] = icdf_spec(samp[indxsampspec[l][i, :]], sampvarb[indxsampfdfnslop[l, i]], globdata.minmspec[i], globdata.maxmspec[i])
            
        if globdata.colrprio:
            sampvarb[indxsampsind[l]] = icdf_atan(samp[indxsampsind[l]], globdata.minmsind, globdata.factsind)
            sampvarb[indxsampspec[l]] = retr_spec(sampvarb[indxsampspec[l][globdata.indxenerfdfn, :]], sampvarb[indxsampsind[l]])
            
        indxpixlpntstemp = retr_indxpixl(globdata, sampvarb[indxsampbgal[l]], sampvarb[indxsamplgal[l]])
    
        cntstemp = sampvarb[indxsampspec[l]][:, :, None] *             globdata.expo[:, indxpixlpntstemp, :] *             globdata.diffener[:, None, None]
            
        cnts.append(cntstemp)
        indxpixlpnts.append(indxpixlpntstemp)
        listspectemp.append(sampvarb[indxsampspec[l]])
        

    pntsflux = retr_pntsflux(globdata, sampvarb[concatenate(indxsamplgal)],                                            sampvarb[concatenate(indxsampbgal)],                                            concatenate(listspectemp, axis=1), sampvarb[indxsamppsfipara])
    
    totlflux = retr_rofi_flux(sampvarb[indxsampnormback], pntsflux, fullindx)
    totlcnts = totlflux * globdata.expo * apix * globdata.diffener[:, None, None] # [1]
    
    return sampvarb, indxpixlpnts, cnts, pntsflux, totlflux, totlcnts
    
    
def retr_mrkrsize(globdata, spec, indxenertemp):

    mrkrsize = (spec - globdata.minmspec[indxenertemp]) / (globdata.maxmspec[indxenertemp] - globdata.minmspec[indxenertemp]) *                     (maxmmrkrsize - minmmrkrsize) + minmmrkrsize
        
    return mrkrsize


def retr_scalfromangl(globdata, thisangl, i, m):
    scalangl = 2. * arcsin(.5 * sqrt(2. - 2 * cos(thisangl))) / globdata.fermscalfact[i, m]
    return scalangl

def retr_anglfromscal(globdata, scalangl, i, m):
    thisangl = arccos(1. - 0.5 * (2. * sin(scalangl * globdata.fermscalfact[i, m] / 2.))**2)
    return thisangl


def retr_fermpsfn(globdata):
   
    name = os.environ["PNTS_TRAN_DATA_PATH"] + '/irf/psf_P8R2_SOURCE_V6_PSF.fits'
    irfn = pf.getdata(name, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']

    globdata.numbfermformpara = 5
    globdata.numbfermscalpara = 3
    
    fermscal = zeros((globdata.numbevtt, globdata.numbfermscalpara))
    fermform = zeros((globdata.numbener, globdata.numbevtt, globdata.numbfermformpara))
    fermpsfipara = zeros((globdata.numbener * globdata.numbfermformpara * globdata.numbevtt))
    for m in globdata.indxevtt:
        fermscal[m, :] = pf.getdata(name, 2 + 3 * globdata.indxevttincl[m])['PSFSCALE']
        irfn = pf.getdata(name, 1 + 3 * globdata.indxevttincl[m])
        for k in range(5):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(globdata.meanener)

    globdata.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * globdata.meanener[:, None])**fermscal[None, :, 2])**2 +                         fermscal[None, :, 1]**2)
    
    # convert N_tail to f_core
    for m in globdata.indxevtt:
        for i in globdata.indxener:
            fermform[i, m, 1] = retr_anglfromscal(globdata, fermform[i, m, 1], i, m) # [rad]
            fermform[i, m, 3] = retr_anglfromscal(globdata, fermform[i, m, 3], i, m) # [rad]
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)
    
    # store the fermi PSF parameters
    for m in globdata.indxevtt:
        for k in range(globdata.numbfermformpara):
            indxfermpsfiparatemp = m * globdata.numbfermformpara * globdata.numbener + globdata.indxener * globdata.numbfermformpara + k
            fermpsfipara[indxfermpsfiparatemp] = fermform[:, m, k]
        
    frac = fermform[:, :, 0]
    sigc = fermform[:, :, 1]
    gamc = fermform[:, :, 2]
    sigt = fermform[:, :, 3]
    gamt = fermform[:, :, 4]
    
    fermpsfn = retr_doubking(globdata.angldisp[None, :, None], frac[:, None, :], sigc[:, None, :], gamc[:, None, :],                          sigt[:, None, :], gamt[:, None, :])
            
    return fermpsfn, fermpsfipara


def updt_samp(globdata):
    
    global thissampvarb, thispntsflux, thismodlcnts, thisindxpntsfull,         thisindxpntsempt, thisllik, thislpri

    if globdata.thisindxprop == indxpropfdfnnorm:
        thissampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]] = nextsampvarb[indxsampfdfnnorm[globdata.indxpoplmodi]]
        thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]

    if globdata.thisindxprop == indxpropfdfnslop:
 
        # update the unit sample vector
        drmcsamp[thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :], -1] =             cdfn_spec(thissampvarb[thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :]],                       nextsampvarb[indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]],                       globdata.minmspec[globdata.indxenermodi], globdata.maxmspec[globdata.indxenermodi])
        
        # update the sample vector
        thissampvarb[indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]] =             nextsampvarb[indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]]
            
        # update the prior register
        thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]

    # likelihood updates
    if globdata.thisindxprop >= indxproppsfipara:
        thisllik[modiindx] = nextllik[modiindx]
        thismodlcnts[modiindx] = nextmodlcnts[modiindx]
        
    if globdata.thisindxprop == indxproppsfipara:
        thissampvarb[indxsamppsfipara[indxpsfiparamodi]] = nextpsfipara[indxpsfiparamodi]
        
    if globdata.thisindxprop == indxpropnormback:
        thissampvarb[indxsampnormback[indxbackmodi, globdata.indxenermodi]] =             nextsampvarb[indxsampnormback[indxbackmodi, globdata.indxenermodi]]
        
    if globdata.thisindxprop >= indxpropbrth or globdata.thisindxprop == indxproppsfipara:
        thispntsflux[modiindx] = nextpntsflux[modiindx]
        
    # transdimensinal updates
    if globdata.thisindxprop >= indxpropbrth and globdata.thisindxprop <= indxpropmerg:
        thissampvarb[indxsampnumbpnts[globdata.indxpoplmodi]] = nextsampvarb[indxsampnumbpnts[globdata.indxpoplmodi]]
        thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]
        
    # birth
    if globdata.thisindxprop == indxpropbrth:
        
        # update the PS index lists
        thisindxpntsfull[globdata.indxpoplmodi].append(thisindxpntsempt[globdata.indxpoplmodi][0])
        del thisindxpntsempt[globdata.indxpoplmodi][0]

        # update the components
        thissampvarb[indxsampmodi[0]] = modilgal
        thissampvarb[indxsampmodi[1]] = modibgal
        if globdata.colrprio:
            thissampvarb[indxsampmodi[2+globdata.indxener]] = modispec
            thissampvarb[indxsampmodi[2+globdata.numbener]] = modisind
        else:
            thissampvarb[indxsampmodi[2:]] = modispec
            
        
    # death
    if globdata.thisindxprop == indxpropdeth:
        
        # update the PS index lists
        thisindxpntsempt[globdata.indxpoplmodi].append(killindxpnts)
        thisindxpntsfull[globdata.indxpoplmodi].remove(killindxpnts)


    # split
    if globdata.thisindxprop == indxpropsplt:

        # update the PS index lists
        thisindxpntsfull[globdata.indxpoplmodi].append(thisindxpntsempt[globdata.indxpoplmodi][0])
        del thisindxpntsempt[globdata.indxpoplmodi][0]
        
        # update the components
        # first component
        thissampvarb[indxinit0] = modilgal[1]
        thissampvarb[indxinit0+1] = modibgal[1]
        thissampvarb[indxinit0+2:indxinit0+2+globdata.numbener] = modispec[:, 1]
  
        # second component
        thissampvarb[indxinit1] = modilgal[2]
        thissampvarb[indxinit1+1] = modibgal[2]
        thissampvarb[indxinit1+2:indxinit1+2+globdata.numbener] = modispec[:, 2]
        
    # merge
    if globdata.thisindxprop == indxpropmerg:
        
        # update the PS index lists
        thisindxpntsfull[globdata.indxpoplmodi].remove(mergindxpnts1)
        thisindxpntsempt[globdata.indxpoplmodi].append(mergindxpnts1)

        # update the component
        thissampvarb[indxsampmodi[0]] = modilgal[2]
        thissampvarb[indxsampmodi[1]] = modibgal[2]
        thissampvarb[indxsampmodi[2:]] = modispec[:, 2]
        
        
    # component change
    if globdata.thisindxprop >= indxproplgal:  
        if globdata.thisindxprop == indxproplgal:
            thissampvarb[indxsampmodi] = modilgal[1]
        elif globdata.thisindxprop == indxpropbgal:
            thissampvarb[indxsampmodi] = modibgal[1]
        else:
            if globdata.colrprio:
                if False:
                    print 'hey'
                    print 'modispec'
                    print modispec
                    print 'modisind'
                    print modisind
                thissampvarb[indxsampmodispec] = modispec[:, 1]
                if globdata.thisindxprop == indxpropsind:
                    thissampvarb[indxsampmodi] = modisind
            else:
                thissampvarb[indxsampmodi] = modispec[0, 1]

    # update the unit sample vector
    if indxsampmodi.size > 0:
        drmcsamp[indxsampmodi, 0] = drmcsamp[indxsampmodi, 1]


def retr_postvarb(listvarb):

    shap = zeros(len(listvarb.shape), dtype=int)
    shap[0] = 3
    shap[1:] = listvarb.shape[1:]
    shap = list(shap)
    postvarb = zeros(shap)
    
    postvarb[0, :] = percentile(listvarb, 50., axis=0)
    postvarb[1, :] = percentile(listvarb, 16., axis=0)
    postvarb[2, :] = percentile(listvarb, 84., axis=0)

    return postvarb


def retr_errrvarb(postvarb):

    errr = abs(postvarb[0, :] - postvarb[1:3, :])

    return errr


def retr_pairlist(lgal, bgal):
    
    pairlist = []
    for k in range(lgal.size):
        indxpnts = k + 1 + where((lgal[k+1:] < lgal[k] + spmrlbhl) &                               (lgal[k+1:] > lgal[k] - spmrlbhl) &                               (bgal[k+1:] < bgal[k] + spmrlbhl) &                               (bgal[k+1:] > bgal[k] - spmrlbhl))[0]
        for l in range(indxpnts.size):
            pairlist.append([k, indxpnts[l]])
            
    return pairlist


def pair_catl(globdata, modllgal, modlbgal, modlspec):

    indxmodl = zeros_like(truelgal, dtype=int) - 1
    dir2 = array([modllgal, modlbgal])
    for k in range(truelgal.size):
        dir1 = array([truelgal[k], truebgal[k]])
        dist = angdist(dir1, dir2, lonlat=True)
        jdist = argmin(dist) 
        if dist[jdist] < deg2rad(0.5):
            indxmodl[k] = jdist

    jtruepntsbias = where(amax(abs(modlspec[:, indxmodl] - truespec) / truespec, axis=0) > 1.2)[0]
    jtruepntsmiss = where(indxmodl == -1)[0]
    
    return indxmodl, jtruepntsbias, jtruepntsmiss


def retr_fgl3(globdata):
        
    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/catl/gll_psc_v16.fit'

    fgl3 = pf.getdata(path)
    
    fgl3numbpnts = fgl3['glon'].size
    
    fgl3lgal = fgl3['glon']
    fgl3lgal = ((fgl3lgal - 180.) % 360.) - 180.

    fgl3bgal = fgl3['glat']

    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3sind = fgl3['Spectral_Index']
    
    fgl3spectype = fgl3['SpectrumType']
    fgl3scur = fgl3['beta']
    fgl3scut = fgl3['Cutoff'] * 1e-3
    
    fgl3timevari = fgl3['Variability_Index']
    
    fgl3spectemp = stack((fgl3['Flux100_300'],                           fgl3['Flux300_1000'],                           fgl3['Flux1000_3000'],                           fgl3['Flux3000_10000'],                           fgl3['Flux10000_100000']))[globdata.indxenerincl, :] / globdata.diffener[:, None]
    fgl3specstdv = stack((fgl3['Unc_Flux100_300'],                           fgl3['Unc_Flux300_1000'],                           fgl3['Unc_Flux1000_3000'],                           fgl3['Unc_Flux3000_10000'],                           fgl3['Unc_Flux10000_100000']))[globdata.indxenerincl, :, :] / globdata.diffener[:, None, None]
    
    fgl3spec = zeros((3, globdata.numbener, fgl3numbpnts))
    fgl3spec[0, :, :] = fgl3spectemp
    fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdv[:, :, 0]
    fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdv[:, :, 1]
    
    # get PS counts
    indxpixlfgl3 = retr_indxpixl(globdata, fgl3bgal, fgl3lgal)
    fgl3cnts = fgl3spec[0, :, :, None] * globdata.expo[:, indxpixlfgl3, :] * globdata.diffener[:, None, None]
    fgl3gang = rad2deg(arccos(cos(deg2rad(fgl3lgal)) * cos(deg2rad(fgl3bgal))))
    
    return fgl3lgal, fgl3bgal, fgl3spec, fgl3gang, fgl3cnts,         fgl3timevari, fgl3sind, fgl3spectype, fgl3scur, fgl3scut
    
    
def retr_rtag(globdata, indxprocwork):
    
    if indxprocwork == None:
        rtag = 'AA_%d_%d_%d_%d_%s_%s_%s' % (globdata.numbproc, globdata.numbswep, globdata.numbburn,                                             globdata.factthin, globdata.datatype, globdata.regitype,                                             globdata.psfntype)
    else:
        rtag = '%02d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, globdata.numbproc, globdata.numbswep,                                               globdata.numbburn, globdata.factthin, globdata.datatype,                                               globdata.regitype, globdata.psfntype)
        
    return rtag


def retr_gaus(indxsamp, stdv):
    
    if fracrand > 0.:
        if rand() < fracrand:
            drmcsamp[indxsamp, 1] = rand()
        else:
            drmcsamp[indxsamp, 1] = drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        drmcsamp[indxsamp, 1] = drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_dist(globdata, lgal0, bgal0, lgal1, bgal1):
    
    if globdata.pixltype == 'heal':
        dir1 = array([lgal0, bgal0])
        dir2 = array([lgal1, bgal1])
        dist = angdist(dir1, dir2, lonlat=True) # [rad]
    else:
        dist = deg2rad(sqrt((lgal0 - lgal1)**2 + (bgal0 - bgal1)**2))

    return dist
       
    
def retr_propstrg():
    
    propstrg = [r'$\mu$', r'$\alpha$', r'$\beta$',                  'PSF', 'Diffuse Norm', 'Isotropic Norm', 'Birth', 'Death',                  'Split', 'Merge', 'l-update', 'b-update', 'f-update']
    
    return propstrg
    
    
def gmrb_test(griddata):
    
    withvari = mean(var(griddata, 0))
    btwnvari = griddata.shape[0] * var(mean(griddata, 0))
    wgthvari = (1. - 1. / griddata.shape[0]) * withvari + btwnvari / griddata.shape[0]
    psrf = sqrt(wgthvari / withvari)

    return psrf


def retr_psfn(globdata, psfipara, indxenertemp, thisangl):

    thisangltemp = thisangl[None, :, None]

    indxpsfiparatemp = indxenertemp[:, None] * globdata.numbformpara +         globdata.indxevtt[None, :] * globdata.numbpsfiparaevtt
    
    if globdata.psfntype == 'singgaus':
        sigc = psfipara[indxpsfiparatemp]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(thisangltemp, sigc)

    elif globdata.psfntype == 'singking':
        sigc = psfipara[indxpsfiparatemp]
        gamc = psfipara[indxpsfiparatemp+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(thisangltemp, sigc, gamc)
        
    elif globdata.psfntype == 'doubgaus':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(thisangltemp, frac, sigc, sigt)

    elif globdata.psfntype == 'gausking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        gamt = psfipara[indxpsfiparatemp+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(thisangltemp, frac, sigc, sigt, gamt)
        
    elif globdata.psfntype == 'doubking':
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


def retr_singgaus(scaldevi, sigc):
    
    psfn = 1. / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_singking(scaldevi, sigc, gamc):
    
    psfn = 1. / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc)

    return psfn


def retr_doubgaus(scaldevi, frac, sigc, sigt):
    
    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) +                 (1. - frac) / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_gausking(scaldevi, frac, sigc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) +                 (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) +                 (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def chsq_fdfnslop(globdata, para, i):

    fdfnslop = para[0]
    fdfnnormnorm = para[1]
    
    fluxhistmodl = fdfnnormnorm * diffspec[i, :] * pdfn_spec(meanspec[i, :], fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])

    chsq = sum(((fluxhistmodl.flatten()[jspecfgl3] - fgl3spechist[i, jspecfgl3]) / fgl3spechist[i, jspecfgl3])**2)
    
    return chsq


def retr_psfimodl(globdata):
    
    numbformpara = globdata.numbformpara
    numbpsfiparaevtt = globdata.numbpsfiparaevtt

    minmformpara = zeros(numbformpara)
    maxmformpara = zeros(numbformpara)
    factformpara = zeros(numbformpara)
    scalformpara = zeros(numbformpara, dtype=object)
    if globdata.exprtype == 'ferm':
        minmanglpsfn = deg2rad(0.01)
        maxmanglpsfn = deg2rad(3.)
        minmgamm = 2.
        maxmgamm = 20.
        minmpsfnfrac = 0.
        maxmpsfnfrac = 1.
    if globdata.exprtype == 'sdss':
        minmanglpsfn = deg2rad(0.01 / 3600.)
        maxmanglpsfn = deg2rad(2. / 3600.)
    if globdata.psfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif globdata.psfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif globdata.psfntype == 'doubgaus':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmanglpsfn
        maxmformpara[2] = maxmanglpsfn
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'logt'
        strgformpara = ['$f_c', r'$\sigma_c', r'$\sigma_t']
    elif globdata.psfntype == 'gausking':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmanglpsfn
        maxmformpara[2] = maxmanglpsfn
        minmformpara[3] = minmgamm
        maxmformpara[3] = maxmgamm
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'logt'
        scalformpara[3] = 'atan'
        strgformpara = ['$f_g', r'$\sigma_g', r'$\sigma_k', r'$\gamma']
    elif globdata.psfntype == 'doubking':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmgamm
        maxmformpara[2] = maxmgamm
        minmformpara[3] = minmanglpsfn
        maxmformpara[3] = maxmanglpsfn
        minmformpara[4] = minmgamm
        maxmformpara[4] = maxmgamm
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'atan'
        scalformpara[3] = 'logt'
        scalformpara[4] = 'atan'
        strgformpara = ['$f_c', r'$\sigma_c', r'$\gamma_c', r'$\sigma_t', r'$\gamma_t']
    
    for k in range(numbformpara):
        if scalformpara[k] == 'self':
            factformpara[k] = maxmformpara[k] - minmformpara[k]
        if scalformpara[k] == 'logt':
            factformpara[k] = log(maxmformpara[k] / minmformpara[k])
        if scalformpara[k] == 'atan':
            factformpara[k] = arctan(maxmformpara[k]) - arctan(minmformpara[k])
            
    minmpsfipara = tile(tile(minmformpara, globdata.numbener), globdata.numbevtt)
    maxmpsfipara = tile(tile(maxmformpara, globdata.numbener), globdata.numbevtt)
    scalpsfipara = tile(tile(scalformpara, globdata.numbener), globdata.numbevtt)
    factpsfipara = tile(tile(factformpara, globdata.numbener), globdata.numbevtt)
    
    strgpsfipara = [strgformpara[k] + '^{%d%d}$' % (globdata.indxenerincl[i], globdata.indxevttincl[m])                         for m in globdata.indxevtt for i in globdata.indxener for k in range(numbformpara)]

    indxpsfiparainit = (arange(globdata.numbevtt)[:, None] * numbpsfiparaevtt + arange(globdata.numbener)[None, :] * numbformpara).flatten()

    for k in arange(numbformpara):
        if globdata.psfntype == 'singgaus' or globdata.psfntype == 'singking':
            strgpsfipara[indxpsfiparainit[k]] += ' ' + globdata.strganglunit
        elif globdata.psfntype == 'doubgaus' or globdata.psfntype == 'gausking':
            strgpsfipara[indxpsfiparainit[k]+1] += ' ' + globdata.strganglunit
            strgpsfipara[indxpsfiparainit[k]+2] += ' ' + globdata.strganglunit
        elif globdata.psfntype == 'doubking':
            strgpsfipara[indxpsfiparainit[k]+1] += ' ' + globdata.strganglunit
            strgpsfipara[indxpsfiparainit[k]+3] += ' ' + globdata.strganglunit

    return minmpsfipara, maxmpsfipara, factpsfipara, strgpsfipara, scalpsfipara, indxpsfiparainit


def retr_strgfluxunit(globdata):
    
    if globdata.exprtype == 'ferm':
        strgfluxunit = r'[1/cm$^2$/s/GeV]'
    if globdata.exprtype == 'sdss':
        strgfluxunit = '[nMgy]'
        
    return strgfluxunit
     
    
def retr_enerstrg(globdata):
    
    if globdata.exprtype == 'ferm':
        binsenerstrg = []
        enerstrg = []
        for i in globdata.indxener:
            binsenerstrg.append('%.3g GeV - %.3g GeV' % (globdata.binsener[i], globdata.binsener[i+1]))
            enerstrg.append('%.3g GeV' % globdata.meanener[i])
            
    if globdata.exprtype == 'sdss':
        globdata.binsenerstrg = ['i-band', 'r-band', 'g-band']
        enerstrg = ['i-band', 'r-band', 'g-band']
        
    return enerstrg, binsenerstrg



# In[ ]:




# In[ ]:



