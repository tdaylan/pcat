
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
import tdpy_util

# pnts_tran
from pnts_tran.cnfg import *
from pnts_tran.main import *
from pnts_tran.samp import *
from pnts_tran.util import *
from pnts_tran.visu import *
from pnts_tran.plot import *



# In[ ]:

def retr_fwhm(globdata, psfn):
    
    fwhm = zeros((numbener, numbevtt))
    for i in indxener:
        for m in indxevtt:
            jangldisp = argsort(psfn[i, :, m])
            intpfwhm = max(amax(psfn[i, :, m]) / 2., amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, jangldisp, m]) and intpfwhm < amax(psfn[i, jangldisp, m]):
                fwhm[i, m] = 2. * interp1d(psfn[i, jangldisp, m], angldisp[jangldisp])(intpfwhm) # [rad]
    return fwhm


def retr_spec(globdata, flux, sind):

    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])
        
    spec = flux[None, :] * (meanener[:, None] / enerfdfn)**(-sind[None, :])
    
    return spec


def retr_indx(globdata, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampsind = []
    indxsampcomp = []
    for l in indxpopl:
        indxsamplgaltemp = indxcompinit + maxmnumbcomp * l + array(indxpntsfull[l], dtype=int) * numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], numbener, 0) +                          repeat(arange(numbener), len(indxpntsfull[l])).reshape(numbener, -1))
        if colrprio:
            indxsampsind.append(indxsamplgaltemp + 2 + numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, numbcomp) + tile(arange(numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(globdata, lgal, bgal, spec, psfipara):
    
    numbpnts = lgal.size
    
    dist = empty((numbpixl, numbpnts))
    for k in range(numbpnts):
        dist[:, k] = retr_dist(lgal[k], bgal[k], lgalgrid, bgalgrid)

    # convolve with the PSF
    pntsflux = empty((numbpnts, numbener, numbpixl, numbevtt))
    for k in range(numbpnts):
        psfn = retr_psfn(psfipara, indxener, dist[:, k], psfntype=modlpsfntype)
        pntsflux[k, :, :, :] = spec[:, k, None, None] * psfn

    # sum contributions from all PS
    pntsfluxtemp = sum(pntsflux, 0) 

    return pntsfluxtemp


def retr_rofi_flux(globdata, normback, pntsflux, tempindx):

    modlflux = pntsflux[tempindx]
    for c in indxback:
        modlflux += normback[c, :, None, None] * backflux[c][tempindx]        
    
    return modlflux


def cdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))

    if flux <= fluxbrek:
        
        fluxunit = norm / (1. - fdfnsloplowr) *             ((flux / fluxbrek)**(1. - fdfnsloplowr) - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
        
    else:
        
        fluxunit = norm / (1. - fdfnsloplowr) * (1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -             norm / (1. - fdfnslopuppr) * (1. - (flux / fluxbrek)**(1. - fdfnslopuppr))
       
    return fluxunit


def pdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))

    if flux <= fluxbrek:
        
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnsloplowr)
        
    else:
        
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnslopuppr)
        
    return pdfnflux


def icdf_spec_brok(globdata, fluxunit, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):
    
    norm = 1. / ((1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) -                  (1. - (maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    
    fluxunitbrek = norm / (1. - fdfnsloplowr) * (1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
    
    if fluxunit < fluxunitbrek:
        
        flux = fluxbrek * (fluxunit * (1. - fdfnsloplowr) / norm +                        (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))**(1. / (1. - fdfnsloplowr))
        
    else:
        
        flux = fluxbrek * (1. - (norm / (1. - fdfnsloplowr) * (1. - (minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -                                  fluxunit) * (1. - fdfnslopuppr) / norm)**(1. / (1. - fdfnslopuppr))

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


def icdf_psfipara(globdata, psfiparaunit, indxpsfipara):
    
    if scalpsfipara[indxpsfipara] == 'self':
        psfipara = icdf_self(psfiparaunit, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])
    if scalpsfipara[indxpsfipara] == 'logt':
        psfipara = icdf_logt(psfiparaunit, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])
    if scalpsfipara[indxpsfipara] == 'atan':
        psfipara = icdf_atan(psfiparaunit, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])

    return psfipara


def cdfn_psfipara(globdata, psfipara, indxpsfipara):
    
    if scalpsfipara[indxpsfipara] == 'self':
        psfiparaunit = cdfn_self(psfipara, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])
    if scalpsfipara[indxpsfipara] == 'logt':
        psfiparaunit = cdfn_logt(psfipara, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])
    if scalpsfipara[indxpsfipara] == 'atan':
        psfiparaunit = cdfn_atan(psfipara, minmpsfipara[indxpsfipara], factpsfipara[indxpsfipara])
    
    return psfiparaunit


def retr_indxprop(globdata, samp):
    
    indxpoplmodi = choice(indxpopl)
    
    numbpnts = int(thissampvarb[indxsampnumbpnts[indxpoplmodi]])
        
    if numbpnts == maxmnumbpnts[indxpoplmodi]:
        indxprop = choice(prop, p=probpropmaxm)
    elif numbpnts == minmnumbpnts:
        indxprop = choice(prop, p=probpropminm)
    else:
        indxprop = choice(prop, p=probprop)
        
    return indxprop
        

def retr_pixl(globdata, bgal, lgal):

    if pixltype == 'heal':
        pixl = pixlcnvt[ang2pix(numbsideheal, deg2rad(90. - bgal), deg2rad(lgal))]
    else:
        
        jlgcr = floor(nsidecart * (lgal - minmlgal) / 2. / maxmgang).astype(int)
        jbgcr = floor(nsidecart * (bgal - minmbgal) / 2. / maxmgang).astype(int)

        if isscalar(jlgcr):
            if jlgcr < 0:
                jlgcr = 0
            if jlgcr >= nsidecart:
                jlgcr = nsidecart - 1
        else:
            jlgcr[where(jlgcr < 0)] = 0
            jlgcr[where(jlgcr >= nsidecart)] = nsidecart - 1
            
        if isscalar(jbgcr):
            if jbgcr < 0:
                jbgcr = 0
            if jbgcr >= nsidecart:
                jbgcr = nsidecart - 1
        else:
            jbgcr[where(jbgcr < 0)] = 0
            jbgcr[where(jbgcr >= nsidecart)] = nsidecart - 1
            
        pixl = jbgcr * nsidecart + jlgcr

    return pixl


def retr_llik(globdata, init=False):

    global thisllik, nextmodlflux, nextmodlcnts, deltllik, nextllik, thismodlcnts, modiindx
    
    if init:

        thisllik = datacnts * log(thismodlcnts) - thismodlcnts
        
    elif indxprop >= indxproppsfipara:

        # determine pixels over which to evaluate the log-likelihood
        if indxprop == indxpropnormback:
            mpixl = ipixl
            
        if indxprop == indxproppsfipara or indxprop >= indxpropbrth:
            
            if indxprop == indxproppsfipara:
                numbpnts = int(sum(thissampvarb[indxsampnumbpnts]))
                lgal = thissampvarb[concatenate(thisindxsamplgal)]
                bgal = thissampvarb[concatenate(thisindxsampbgal)]
                spec = thissampvarb[concatenate(thisindxsampspec)[indxenermodi, :]]
            else:
                numbpnts = modilgal.shape[0]
                lgal = modilgal
                bgal = modibgal
                spec = modispec
                
            xpixlclos = []
            for k in range(numbpnts):
                jspecclos = argmin(spec[0, k] - meanspecprox)
                jpixl = retr_pixl(bgal[k], lgal[k])
                xpixlclos.append(ypixl[jspecclos][jpixl])
                
            mpixl = unique(concatenate(xpixlclos))


            #print 'xpixlclos[0].size: '
            #print xpixlclos[0].size
            #print 'mpixl.size'
            #print mpixl.size
            #print
            
        # construct the mesh grid for likelihood evaluation
        if indxprop >= indxproppsfipara:
            modiindx = meshgrid(indxenermodi, mpixl, indxevtt, indexing='ij')

        # update the point source flux map
        if indxprop == indxproppsfipara or indxprop >= indxpropbrth:

            if indxprop == indxproppsfipara:
                nextpntsflux[modiindx] = 0.
            else:
                nextpntsflux[modiindx] = thispntsflux[modiindx]
                
            for k in range(numbpnts):
                
                # calculate the distance to the pixels to be updated
                dist = retr_dist(lgal[k], bgal[k], lgalgrid[xpixlclos[k]], bgalgrid[xpixlclos[k]])

                # evaluate the PSF over the set of data cubes to be updated
                if indxprop == indxproppsfipara:
                    temppsfipara = nextpsfipara
                else:
                    temppsfipara = thissampvarb[indxsamppsfipara]
                psfn = retr_psfn(temppsfipara, indxenermodi, dist, psfntype=modlpsfntype)
                                
                # update the data cubes
                for i in range(indxenermodi.size):
                    nextpntsflux[indxenermodi[i], xpixlclos[k], :] += spec[i, k] * psfn[i, :, :]
            
        if indxprop == indxpropnormback:
            
            normbacktemp = empty((numbback, 1))
            for c in indxback:
                if c == indxbackmodi:
                    normbacktemp[c, 0] = nextsampvarb[indxsampnormback[c, indxenermodi]]
                else:
                    normbacktemp[c, 0] = thissampvarb[indxsampnormback[c, indxenermodi]]
                
            nextmodlflux[modiindx] = retr_rofi_flux(normbacktemp, thispntsflux, modiindx)

        if indxprop == indxproppsfipara or indxprop >= indxpropbrth:
            nextmodlflux[modiindx] = retr_rofi_flux(thissampvarb[indxsampnormback[meshgrid(indxback, indxenermodi, indexing='ij')]],                                                     nextpntsflux, modiindx)

        nextmodlcnts[modiindx] = nextmodlflux[modiindx] * expo[modiindx] * apix * diffener[indxenermodi, None, None] # [1]
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
        thislpri = zeros((numbpopl, numbener))
        
        for i in indxenerfdfn:
            for l in indxpopl:
                fluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[l]], thissampvarb[indxsampfdfnslop[l, i]], i)
                fluxhist = histogram(thissampvarb[thisindxsampspec[l][i, :]], binsspec[i, :])[0]
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                thislpri[l, i] = sum(lprbpois)            
      
        nextlpri = copy(thislpri)
                
    else:
        nextlpri = copy(thislpri)
        if indxprop == indxpropfdfnnorm or indxprop == indxpropfdfnslop             or indxprop >= indxpropbrth and indxprop <= indxpropmerg:
              
            if colrprio:
                indxenertemp = indxenerfdfn
            else:
                indxenertemp = indxenermodi
            for i in indxenertemp:
                if indxprop == indxpropfdfnnorm:
                    fdfnnorm = nextsampvarb[indxsampfdfnnorm[indxpoplmodi]]
                    fdfnslop = thissampvarb[indxsampfdfnslop[indxpoplmodi, i]]
                elif indxprop == indxpropfdfnslop:
                    fdfnnorm = thissampvarb[indxsampfdfnnorm[indxpoplmodi]]
                    fdfnslop = nextsampvarb[indxsampfdfnslop[indxpoplmodi, i]]
                else:
                    fdfnnorm = thissampvarb[indxsampfdfnnorm[indxpoplmodi]]
                    fdfnslop = thissampvarb[indxsampfdfnslop[indxpoplmodi, i]]
                fluxhistmodl = retr_fdfn(fdfnnorm, fdfnslop, i)
                              
                fluxhist = histogram(thissampvarb[thisindxsampspec[indxpoplmodi][i, :]], binsspec[i, :])[0] 
                if indxprop == indxpropbrth:
                    fluxhist += histogram(modispec[i, 0], binsspec[i, :])[0]
                elif indxprop == indxpropdeth:
                    fluxhist -= histogram(-modispec[i, 0], binsspec[i, :])[0]
                elif indxprop == indxpropsplt:
                    fluxhist -= histogram(-modispec[i, 0], binsspec[i, :])[0]
                    fluxhist += histogram(modispec[i, 1:3], binsspec[i, :])[0]
                elif indxprop == indxpropmerg:
                    fluxhist -= histogram(-modispec[i, 0:2], binsspec[i, :])[0]
                    fluxhist += histogram(modispec[i, 2], binsspec[i, :])[0]
                
                if False:
                    
                    prevfluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[indxpoplmodi]], thissampvarb[indxsampfdfnslop[indxpoplmodi, i]], i)
                    get_ipython().magic(u'matplotlib inline')
                    plt.loglog(meanspec[i, :], fluxhistmodl)
                    plt.loglog(meanspec[i, :], fluxhist)
                    plt.loglog(meanspec[i, :], prevfluxhistmodl, ls='--')
                    plt.show()


                
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                nextlpri[indxpoplmodi, i] = sum(lprbpois)


                
            
            deltlpri = sum(nextlpri[indxpoplmodi, indxenermodi] - thislpri[indxpoplmodi, indxenermodi])

                      
        else:
            deltlpri = 0.

        
def pars_samp(globdata, indxpntsfull, samp):
    
    cnts = []
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[indxsampnumbpnts] = samp[indxsampnumbpnts]
    sampvarb[indxsampfdfnnorm] = icdf_logt(samp[indxsampfdfnnorm], minmfdfnnorm, factfdfnnorm)
    sampvarb[indxsampfdfnslop] = icdf_atan(samp[indxsampfdfnslop], minmfdfnslop, factfdfnslop)
         
    for c in indxback:
        sampvarb[indxsampnormback[c, :]] = icdf_logt(samp[indxsampnormback[c, :]], minmnormback[c], factnormback[c])
    for k in ipsfipara:
        sampvarb[indxsamppsfipara[k]] = icdf_psfipara(samp[indxsamppsfipara[k]], k)

    listspectemp = []
    for l in indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -maxmgangmarg, 2. * maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -maxmgangmarg, 2. * maxmgangmarg) 
        for i in indxenerfdfn:
            sampvarb[indxsampspec[l][i, :]] = icdf_spec(samp[indxsampspec[l][i, :]], sampvarb[indxsampfdfnslop[l, i]], minmspec[i], maxmspec[i])
            
        if colrprio:
            sampvarb[indxsampsind[l]] = icdf_atan(samp[indxsampsind[l]], minmsind, factsind)
            sampvarb[indxsampspec[l]] = retr_spec(sampvarb[indxsampspec[l][indxenerfdfn, :]], sampvarb[indxsampsind[l]])
            
        listspectemp.append(sampvarb[indxsampspec[l]])
        
        ppixl = retr_pixl(sampvarb[indxsampbgal[l]], sampvarb[indxsamplgal[l]])
    
        cntstemp = sampvarb[indxsampspec[l]][:, :, None] * expo[:, ppixl, :] * diffener[:, None, None]
        cnts.append(cntstemp)

    pntsflux = retr_pntsflux(sampvarb[concatenate(indxsamplgal)],                                            sampvarb[concatenate(indxsampbgal)],                                            concatenate(listspectemp, axis=1), sampvarb[indxsamppsfipara])
    
    totlflux = retr_rofi_flux(sampvarb[indxsampnormback], pntsflux, fullindx)
    totlcnts = totlflux * expo * apix * diffener[:, None, None] # [1]
    
    return sampvarb, ppixl, cnts, pntsflux, totlflux, totlcnts
    
    
def retr_mrkrsize(globdata, spec, indxenertemp):

    mrkrsize = (spec - minmspec[indxenertemp]) / (maxmspec[indxenertemp] - minmspec[indxenertemp]) *                     (maxmmrkrsize - minmmrkrsize) + minmmrkrsize
        
    return mrkrsize


def retr_scalfromangl(globdata, thisangl, i, m):
    scalangl = 2. * arcsin(.5 * sqrt(2. - 2 * cos(thisangl))) / fermscalfact[i, m]
    return scalangl

def retr_anglfromscal(globdata, scalangl, i, m):
    thisangl = arccos(1. - 0.5 * (2. * sin(scalangl * fermscalfact[i, m] / 2.))**2)
    return thisangl


def retr_fermpsfn(globdata):
   
    name = os.environ["PNTS_TRAN_DATA_PATH"] + '/irf/psf_P8R2_SOURCE_V6_PSF.fits'
    irfn = pf.getdata(name, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']

    global nfermformpara
    nfermformpara = 5
    nfermscalpara = 3
    
    global fermpsfipara
    fermscal = zeros((numbevtt, nfermscalpara))
    fermform = zeros((numbener, numbevtt, nfermformpara))
    fermpsfipara = zeros((numbener * nfermformpara * numbevtt))
    for m in indxevtt:
        fermscal[m, :] = pf.getdata(name, 2 + 3 * indxevttincl[m])['PSFSCALE']
        irfn = pf.getdata(name, 1 + 3 * indxevttincl[m])
        for k in range(5):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(meanener)

            
    global fermscalfact
    fermscalfact = sqrt((fermscal[None, :, 0] * (10. * meanener[:, None])**fermscal[None, :, 2])**2 +                         fermscal[None, :, 1]**2)
    
    # convert N_tail to f_core
    for m in indxevtt:
        for i in indxener:
            fermform[i, m, 1] = retr_anglfromscal(fermform[i, m, 1], i, m) # [rad]
            fermform[i, m, 3] = retr_anglfromscal(fermform[i, m, 3], i, m) # [rad]
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)
    
    # store the fermi PSF parameters
    for m in indxevtt:
        for k in range(nfermformpara):
            fermpsfipara[m*nfermformpara*numbener+indxener*nfermformpara+k] = fermform[:, m, k]
        
    frac = fermform[:, :, 0]
    sigc = fermform[:, :, 1]
    gamc = fermform[:, :, 2]
    sigt = fermform[:, :, 3]
    gamt = fermform[:, :, 4]
    
    psfn = retr_doubking(angldisp[None, :, None], frac[:, None, :], sigc[:, None, :], gamc[:, None, :],                          sigt[:, None, :], gamt[:, None, :])
            
    return psfn


def updt_samp(globdata):
    
    global thissampvarb, thispntsflux, thismodlcnts, thisindxpntsfull,         thisindxpntsempt, thisllik, thislpri

    if indxprop == indxpropfdfnnorm:
        thissampvarb[indxsampfdfnnorm[indxpoplmodi]] = nextsampvarb[indxsampfdfnnorm[indxpoplmodi]]
        thislpri[indxpoplmodi, indxenermodi] = nextlpri[indxpoplmodi, indxenermodi]

    if indxprop == indxpropfdfnslop:
 
        # update the unit sample vector
        drmcsamp[thisindxsampspec[indxpoplmodi][indxenermodi, :], -1] =             cdfn_spec(thissampvarb[thisindxsampspec[indxpoplmodi][indxenermodi, :]],                       nextsampvarb[indxsampfdfnslop[indxpoplmodi, indxenermodi]],                       minmspec[indxenermodi], maxmspec[indxenermodi])
        
        # update the sample vector
        thissampvarb[indxsampfdfnslop[indxpoplmodi, indxenermodi]] =             nextsampvarb[indxsampfdfnslop[indxpoplmodi, indxenermodi]]
            
        # update the prior register
        thislpri[indxpoplmodi, indxenermodi] = nextlpri[indxpoplmodi, indxenermodi]

    # likelihood updates
    if indxprop >= indxproppsfipara:
        thisllik[modiindx] = nextllik[modiindx]
        thismodlcnts[modiindx] = nextmodlcnts[modiindx]
        
    if indxprop == indxproppsfipara:
        thissampvarb[indxsamppsfipara[indxpsfiparamodi]] = nextpsfipara[indxpsfiparamodi]
        
    if indxprop == indxpropnormback:
        thissampvarb[indxsampnormback[indxbackmodi, indxenermodi]] =             nextsampvarb[indxsampnormback[indxbackmodi, indxenermodi]]
        
    if indxprop >= indxpropbrth or indxprop == indxproppsfipara:
        thispntsflux[modiindx] = nextpntsflux[modiindx]
        
    # transdimensinal updates
    if indxprop >= indxpropbrth and indxprop <= indxpropmerg:
        thissampvarb[indxsampnumbpnts[indxpoplmodi]] = nextsampvarb[indxsampnumbpnts[indxpoplmodi]]
        thislpri[indxpoplmodi, indxenermodi] = nextlpri[indxpoplmodi, indxenermodi]
        
    # birth
    if indxprop == indxpropbrth:
        
        # update the PS index lists
        thisindxpntsfull[indxpoplmodi].append(thisindxpntsempt[indxpoplmodi][0])
        del thisindxpntsempt[indxpoplmodi][0]

        # update the components
        thissampvarb[indxsampmodi[0]] = modilgal
        thissampvarb[indxsampmodi[1]] = modibgal
        if colrprio:
            thissampvarb[indxsampmodi[2+indxener]] = modispec
            thissampvarb[indxsampmodi[2+numbener]] = modisind
        else:
            thissampvarb[indxsampmodi[2:]] = modispec
            
        
    # death
    if indxprop == indxpropdeth:
        
        # update the PS index lists
        thisindxpntsempt[indxpoplmodi].append(killindxpnts)
        thisindxpntsfull[indxpoplmodi].remove(killindxpnts)


    # split
    if indxprop == indxpropsplt:

        # update the PS index lists
        thisindxpntsfull[indxpoplmodi].append(thisindxpntsempt[indxpoplmodi][0])
        del thisindxpntsempt[indxpoplmodi][0]
        
        # update the components
        # first component
        thissampvarb[indxinit0] = modilgal[1]
        thissampvarb[indxinit0+1] = modibgal[1]
        thissampvarb[indxinit0+2:indxinit0+2+numbener] = modispec[:, 1]
  
        # second component
        thissampvarb[indxinit1] = modilgal[2]
        thissampvarb[indxinit1+1] = modibgal[2]
        thissampvarb[indxinit1+2:indxinit1+2+numbener] = modispec[:, 2]
        
    # merge
    if indxprop == indxpropmerg:
        
        # update the PS index lists
        thisindxpntsfull[indxpoplmodi].remove(mergindxpnts1)
        thisindxpntsempt[indxpoplmodi].append(mergindxpnts1)

        # update the component
        thissampvarb[indxsampmodi[0]] = modilgal[2]
        thissampvarb[indxsampmodi[1]] = modibgal[2]
        thissampvarb[indxsampmodi[2:]] = modispec[:, 2]
        
        
    # component change
    if indxprop >= indxproplgal:  
        if indxprop == indxproplgal:
            thissampvarb[indxsampmodi] = modilgal[1]
        elif indxprop == indxpropbgal:
            thissampvarb[indxsampmodi] = modibgal[1]
        else:
            if colrprio:
                if False:
                    print 'hey'
                    print 'modispec'
                    print modispec
                    print 'modisind'
                    print modisind
                thissampvarb[indxsampmodispec] = modispec[:, 1]
                if indxprop == indxpropsind:
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
    
    fgl3spectemp = stack((fgl3['Flux100_300'],                           fgl3['Flux300_1000'],                           fgl3['Flux1000_3000'],                           fgl3['Flux3000_10000'],                           fgl3['Flux10000_100000']))[indxenerincl, :] / diffener[:, None]
    fgl3specstdv = stack((fgl3['Unc_Flux100_300'],                           fgl3['Unc_Flux300_1000'],                           fgl3['Unc_Flux1000_3000'],                           fgl3['Unc_Flux3000_10000'],                           fgl3['Unc_Flux10000_100000']))[indxenerincl, :, :] / diffener[:, None, None]
    
    fgl3spec = zeros((3, numbener, fgl3numbpnts))
    fgl3spec[0, :, :] = fgl3spectemp
    fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdv[:, :, 0]
    fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdv[:, :, 1]
    
    # get PS counts
    ppixl = retr_pixl(fgl3bgal, fgl3lgal)
    fgl3cnts = fgl3spec[0, :, :, None] * expo[:, ppixl, :] * diffener[:, None, None]
    fgl3gang = rad2deg(arccos(cos(deg2rad(fgl3lgal)) * cos(deg2rad(fgl3bgal))))
    
    return fgl3lgal, fgl3bgal, fgl3spec, fgl3gang, fgl3cnts,         fgl3timevari, fgl3sind, fgl3spectype, fgl3scur, fgl3scut
    
    
def retr_rtag(globdata, indxprocwork):
    
    if indxprocwork == None:
        rtag = 'A_%d_%d_%d_%d_%s_%s_%s' % (numbproc, numbswep, numbburn, factthin, datatype, regitype, modlpsfntype)
    else:
        rtag = '%d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, numbproc, numbswep, numbburn, factthin, datatype, regitype, modlpsfntype)
        
    return rtag


def retr_gaus(indxsamp, stdv):
    
    if fracrand > 0.:
        if rand() < fracrand:
            drmcsamp[indxsamp, 1] = rand()
        else:
            drmcsamp[indxsamp, 1] = drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        drmcsamp[indxsamp, 1] = drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_dist(lgal0, bgal0, lgal1, bgal1):
    
    if pixltype == 'heal':
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


def retr_psfn(globdata, psfipara, indxenertemp, thisangl, psfntype='doubking'):

    thisangltemp = thisangl[None, :, None]

    indxpsfipara = indxenertemp[:, None] * numbformpara + indxevtt[None, :] * numbpsfiparaevtt
    
    if modlpsfntype == 'singgaus':
        sigc = psfipara[indxpsfipara]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(thisangltemp, sigc)

    elif psfntype == 'singking':
        sigc = psfipara[indxpsfipara]
        gamc = psfipara[indxpsfipara+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(thisangltemp, sigc, gamc)
        
    elif psfntype == 'doubgaus':
        frac = psfipara[indxpsfipara]
        sigc = psfipara[indxpsfipara+1]
        sigt = psfipara[indxpsfipara+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(thisangltemp, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfipara[indxpsfipara]
        sigc = psfipara[indxpsfipara+1]
        sigt = psfipara[indxpsfipara+2]
        gamt = psfipara[indxpsfipara+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(thisangltemp, frac, sigc, sigt, gamt)
        
    elif psfntype == 'doubking':
        frac = psfipara[indxpsfipara]
        sigc = psfipara[indxpsfipara+1]
        gamc = psfipara[indxpsfipara+2]
        sigt = psfipara[indxpsfipara+3]
        gamt = psfipara[indxpsfipara+4]
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
    
    fluxhistmodl = fdfnnormnorm * diffspec[i, :] * pdfn_spec(meanspec[i, :], fdfnslop, minmspec[i], maxmspec[i])

    chsq = sum(((fluxhistmodl.flatten()[jspecfgl3] - fgl3spechist[i, jspecfgl3]) / fgl3spechist[i, jspecfgl3])**2)
    
    return chsq


def retr_psfimodl(globdata, numbformpara, exprtype, psfntype, numbener, numbevtt, indxenerincl, indxevttincl):
    
    if exprtype == 'ferm':
        strganglunit = '[deg]'
    if exprtype == 'sdss':
        strganglunit = '[arcsec]'

    numbpsfiparaevtt = numbformpara * numbener
    
    minmformpara = zeros(numbformpara)
    maxmformpara = zeros(numbformpara)
    factformpara = zeros(numbformpara)
    scalformpara = zeros(numbformpara, dtype=object)
    if exprtype == 'ferm':
        minmanglpsfn = deg2rad(0.01)
        maxmanglpsfn = deg2rad(3.)
        minmgamm = 2.
        maxmgamm = 20.
        minmpsfnfrac = 0.
        maxmpsfnfrac = 1.
    if exprtype == 'sdss':
        minmanglpsfn = deg2rad(0.01 / 3600.)
        maxmanglpsfn = deg2rad(2. / 3600.)
    if psfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif psfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif psfntype == 'doubgaus':
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
    elif psfntype == 'gausking':
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
    elif psfntype == 'doubking':
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
            
    minmpsfipara = tile(tile(minmformpara, numbener), numbevtt)
    maxmpsfipara = tile(tile(maxmformpara, numbener), numbevtt)
    scalpsfipara = tile(tile(scalformpara, numbener), numbevtt)
    factpsfipara = tile(tile(factformpara, numbener), numbevtt)
    
    indxener = arange(numbevtt)
    indxener = arange(numbener)
    strgpsfipara = [strgformpara[k] + '^{%d%d}$' % (indxenerincl[i], indxevttincl[m])                         for m in indxevtt for i in indxener for k in range(numbformpara)]

    indxpsfipara = (arange(numbevtt)[:, None] * numbpsfiparaevtt + arange(numbener)[None, :] * numbformpara).flatten()

    for k in range(indxpsfipara.size):
        if psfntype == 'singgaus' or psfntype == 'singking':
            strgpsfipara[indxpsfipara[k]] += ' ' + strganglunit
        elif psfntype == 'doubgaus' or psfntype == 'gausking':
            strgpsfipara[indxpsfipara[k]+1] += ' ' + strganglunit
            strgpsfipara[indxpsfipara[k]+2] += ' ' + strganglunit
        elif psfntype == 'doubking':
            strgpsfipara[indxpsfipara[k]+1] += ' ' + strganglunit
            strgpsfipara[indxpsfipara[k]+3] += ' ' + strganglunit

    
    return minmpsfipara, maxmpsfipara, factpsfipara, strgpsfipara, scalpsfipara, indxpsfipara


def retr_strgfluxunit(globdata):
    
    if exprtype == 'ferm':
        strgfluxunit = r'[1/cm$^2$/s/GeV]'
    if exprtype == 'sdss':
        strgfluxunit = '[nMgy]'
        
    return strgfluxunit
     
    
def retr_enerstrg(globdata):
    
    if exprtype == 'ferm':
        binsenerstrg = []
        enerstrg = []
        for i in indxener:
            binsenerstrg.append('%.3g GeV - %.3g GeV' % (binsener[i], binsener[i+1]))
            enerstrg.append('%.3g GeV' % meanener[i])
    if exprtype == 'sdss':
        binsenerstrg = ['i-band', 'r-band', 'g-band']
        enerstrg = ['i-band', 'r-band', 'g-band']
        
    return enerstrg, binsenerstrg



# In[ ]:




# In[ ]:



