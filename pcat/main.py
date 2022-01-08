# plotting
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

# numpy
import numpy as np

# scipy
import scipy as sp
import scipy.interpolate
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
import scipy.fftpack
import scipy.sparse

# jit
from numba import jit

import ctypes

import astropy
import astropy as ap
from astropy.convolution import convolve_fft, AiryDisk2DKernel

import pickle

# multiprocessing
import multiprocessing as mp

from copy import deepcopy

# utilities
import os, time, sys, glob, fnmatch, inspect, traceback, functools

# HealPix
import healpy as hp

# ignore warnings if not in diagnostic mode
import warnings
    
#seterr(divide='raise', over='raise', invalid='raise')
#seterr(all='raise')
#seterr(under='ignore')
#warnings.simplefilter('ignore')
#np.set_printoptions(linewidth=180)
#sns.set(context='poster', style='ticks', color_codes=True)

import h5py

# utilities

# secondaries
## Symbolic Jacobian calculation
#import sympy

# tdpy
import tdpy
from tdpy.util import summgene



# photometry related

### find the spectra of sources
def retr_spec(gdat, flux, sind=None, curv=None, expc=None, sindcolr=None, elin=None, edisintp=None, sigm=None, gamm=None, spectype='powr', plot=False):
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if plot:
            meanener = gdat.meanpara.enerplot
        else:
            meanener = gdat.meanpara.ener

        if gmod.spectype == 'gaus':
            spec = 1. / edis[None, :] / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanpara.ener[:, None] - elin[None, :]) / edis[None, :])**2)
        if gmod.spectype == 'voig':
            args = (gdat.meanpara.ener[:, None] + 1j * gamm[None, :]) / np.sqrt(2.) / sigm[None, :]
            spec = 1. / sigm[None, :] / np.sqrt(2. * pi) * flux[None, :] * real(scipy.special.wofz(args))
        if gmod.spectype == 'edis':
            edis = edisintp(elin)[None, :]
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'pvoi':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'lore':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :])
        if gmod.spectype == 'colr':
            if plot:
                spec = np.zeros((gdat.numbenerplot, flux.size))
            else:
                spec = np.empty((gdat.numbener, flux.size))
                for i in gdat.indxener:
                    if i < gdat.indxenerpivt:
                        spec[i, :] = flux * (gdat.meanpara.ener[i] / gdat.enerpivt)**(-sindcolr[i])
                    elif i == gdat.indxenerpivt:
                        spec[i, :] =  flux
                    else:
                        spec[i, :] = flux * (gdat.meanpara.ener[i] / gdat.enerpivt)**(-sindcolr[i-1])
        if gmod.spectype == 'curv':
            spec = flux[None, :] * meanener[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if gmod.spectype == 'expc':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :]) * np.exp(-(meanener - gdat.enerpivt)[:, None] / expc[None, :])
    
    return spec


### find the surface brightness due to one point source
def retr_sbrtpnts(gdat, lgal, bgal, spec, psfnintp, indxpixlelem):
    
    # calculate the distance to all pixels from each point source
    dist = retr_angldistunit(gdat, lgal, bgal, indxpixlelem)
    
    # interpolate the PSF onto the pixels
    if gdat.kernevaltype == 'ulip':
        psfntemp = psfnintp(dist)
    if gdat.kernevaltype == 'bspx':
        pass

    # scale by the PS spectrum
    sbrtpnts = spec[:, None, None] * psfntemp
    
    return sbrtpnts


def retr_psfnwdth(gdat, psfn, frac):
    '''
    Return the PSF width
    '''

    wdth = np.zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            psfntemp = psfn[i, :, m]
            indxanglgood = np.argsort(psfntemp)
            intpwdth = max(frac * np.amax(psfntemp), np.amin(psfntemp))
            if intpwdth >= np.amin(psfntemp[indxanglgood]) and intpwdth <= np.amax(psfntemp[indxanglgood]):
                wdthtemp = sp.interpolate.interp1d(psfntemp[indxanglgood], gdat.binspara.angl[indxanglgood], fill_value='extrapolate')(intpwdth)
            else:
                wdthtemp = 0.
            wdth[i, m] = wdthtemp
                        
    return wdth


# lensing-related
def samp_lgalbgalfromtmpl(gdat, probtmpl):
    
    indxpixldraw = np.random.choice(gdat.indxpixl, p=probtmpl)
    lgal = gdat.lgalgrid[indxpixldraw] + randn(gdat.sizepixl)
    bgal = gdat.bgalgrid[indxpixldraw] + randn(gdat.sizepixl)
    
    return lgal, bgal


## custom random variables, pdfs, cdfs and icdfs
### probability distribution functions
def retr_lprbpois(data, modl):
    
    lprb = data * np.log(modl) - modl - sp.special.gammaln(data + 1)
    
    return lprb
    
        
### probability density functions
def pdfn_self(xdat, minm, maxm):
    
    pdfn = 1. / (maxm - minm)
    
    return pdfn


def pdfn_expo(xdat, maxm, scal):

    if (xdat > maxm).any():
        pdfn = 0.
    else:
        pdfn = 1. / scal / (1. - np.exp(-maxm / scal)) * np.exp(-xdat / scal)

    return pdfn


def pdfn_dexp(xdat, maxm, scal):
    
    pdfn = 0.5 * pdfn_expo(np.fabs(xdat), maxm, scal)

    return pdfn


def pdfn_dpow(xdat, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(xdat):
        xdat = np.array([xdat])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / \
                                            (1. - sloplowr) + (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)
    
    pdfn = np.empty_like(xdat)
    indxlowr = np.where(xdat <= brek)[0]
    indxuppr = np.where(xdat > brek)[0]
    if indxlowr.size > 0:
        pdfn[indxlowr] = faca * brek**(sloplowr - slopuppr) * xdat[indxlowr]**(-sloplowr)
    if indxuppr.size > 0:
        pdfn[indxuppr] = faca * xdat[indxuppr]**(-slopuppr)
    
    return pdfn


def pdfn_powr(xdat, minm, maxm, slop):
  
    norm = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop))
    
    pdfn = norm * xdat**(-slop)
    
    return pdfn


def pdfn_logt(xdat, minm, maxm):
    
    pdfn =  1. / (np.log(maxm) - np.log(minm)) / xdat
    
    return pdfn


def pdfn_igam(xdat, slop, cutf):
    
    pdfn = sp.stats.invgamma.pdf(xdat, slop - 1., scale=cutf)
    
    return pdfn


def pdfn_lnor(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(np.log(xdat), np.log(mean), stdv)

    return pdfn


def pdfn_gaus(xdat, mean, stdv):
    
    pdfn = 1. / np.sqrt(2. * pi) / stdv * np.exp(-0.5 * ((xdat - mean) / stdv)**2)

    return pdfn


def pdfn_lgau(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(np.log(xdat), np.log(mean), stdv)

    return pdfn


def pdfn_atan(para, minmpara, maxmpara):

    pdfn = 1. / (para**2 + 1.) / (np.arctan(maxmpara) - np.arctan(minmpara))
    
    return pdfn


def cdfn_paragenrscalbase(gdat, strgmodl, paragenrscalbase, thisindxparagenrbase):
    
    gmod = getattr(gdat, strgmodl)

    scalparagenrbase = gmod.scalpara.genrbase[thisindxparagenrbase]
    
    if scalparagenrbase == 'self' or scalparagenrbase == 'logt' or scalparagenrbase == 'atan':
        
        listminmparagenrscalbase = gmod.minmpara.genrbase[thisindxparagenrbase]
        factparagenrscalbase = gmod.factparagenrscalbase[thisindxparagenrbase]

        if scalparagenrbase == 'self':
            paragenrscalbaseunit = cdfn_self(paragenrscalbase, listminmparagenrscalbase, factparagenrscalbase)
        elif scalparagenrbase == 'logt':
            paragenrscalbaseunit = cdfn_logt(paragenrscalbase, listminmparagenrscalbase, factparagenrscalbase)

        elif scalparagenrbase == 'atan':
            gmod.listmaxmparagenrscalbase = gmod.listmaxmparagenrscalbase[thisindxparagenrbase]
            paragenrscalbaseunit = cdfn_atan(paragenrscalbase, listminmparagenrscalbase, gmod.listmaxmparagenrscalbase)
    
    elif scalparagenrbase == 'gaus' or scalparagenrbase == 'eerr':
        gmod.listmeanparagenrscalbase = gmod.listmeanparagenrscalbase[thisindxparagenrbase]
        gmod.liststdvparagenrscalbase = gmod.liststdvparagenrscalbase[thisindxparagenrbase]
        if scalparagenrbase == 'eerr':
            gmod.cdfnlistminmparagenrscalbaseunit = gmod.cdfnlistminmparagenrscalbaseunit[thisindxparagenrbase]
            gmod.listparagenrscalbaseunitdiff = gmod.listparagenrscalbaseunitdiff[thisindxparagenrbase]
            paragenrscalbaseunit = cdfn_eerr(paragenrscalbase, gmod.listmeanparagenrscalbase, gmod.liststdvparagenrscalbase, \
                                                                            gmod.cdfnlistminmparagenrscalbaseunit, gmod.listparagenrscalbaseunitdiff)
        else:
            paragenrscalbaseunit = cdfn_gaus(paragenrscalbase, gmod.listmeanparagenrscalbase, gmod.liststdvparagenrscalbase)

    elif scalparagenrbase == 'pois':
        paragenrscalbaseunit = paragenrscalbase
    
    if gdat.booldiagmode:
        if paragenrscalbaseunit == 0:
            print('Warning. CDF is zero.')

    return paragenrscalbaseunit


def icdf_paragenrscalfull(gdat, strgmodl, paragenrunitfull, indxparagenrfullelem):
    
    gmod = getattr(gdat, strgmodl)

    # tobechanged
    # temp -- change zeros to empty
    paragenrscalfull = np.zeros_like(paragenrunitfull)
    for scaltype in gdat.listscaltype:
        listindxparagenrbasescal = gmod.listindxparagenrbasescal[scaltype]
        if len(listindxparagenrbasescal) == 0:
            continue
        paragenrscalfull[listindxparagenrbasescal] = icdf_paragenrscalbase(gdat, strgmodl, paragenrunitfull[listindxparagenrbasescal], scaltype, listindxparagenrbasescal)
        
    if not np.isfinite(paragenrscalfull).all():
        raise Exception('')
        
    if indxparagenrfullelem is not None:
        for l in gmod.indxpopl:
            for g in gmod.indxparagenrelemsing[l]:
                indxparagenrfulltemp = indxparagenrfullelem[l][gmod.namepara.genrelem[l][g]]
                if indxparagenrfulltemp.size == 0:
                    continue
                paragenrscalfull[indxparagenrfulltemp] = icdf_trap(gdat, strgmodl, paragenrunitfull[indxparagenrfulltemp], paragenrscalfull, \
                                                                    gmod.listscalparagenrelem[l][g], gmod.namepara.genrelem[l][g], l)
    
                if gdat.booldiagmode:
                    if not np.isfinite(paragenrscalfull[indxparagenrfulltemp]).all():
                        raise Exception('')

    if not np.isfinite(paragenrscalfull).all():
        raise Exception('')
    
    return paragenrscalfull

    
def icdf_paragenrscalbase(gdat, strgmodl, paragenrunitbase, scaltype, indxparagenrbasescal):
    
    gmod = getattr(gdat, strgmodl)
    
    if scaltype == 'self' or scaltype == 'logt' or scaltype == 'atan':
        minmparagenrscalbase = gmod.minmpara.genrbase[indxparagenrbasescal]
        factparagenrscalbase = gmod.factpara.genrbase[indxparagenrbasescal]

    if scaltype == 'self':
        paragenrscalbase = tdpy.icdf_self(paragenrunitbase, minmparagenrscalbase, factparagenrscalbase)
    elif scaltype == 'logt':
        paragenrscalbase = tdpy.icdf_logt(paragenrunitbase, minmparagenrscalbase, factparagenrscalbase)
    elif scaltype == 'atan':
        listmaxmparagenrscalbase = gmod.listmaxmparagenrscalbase[indxparagenrbasescal]
        paragenrscalbase = tdpy.icdf_atan(paragenrunitbase, minmparagenrscalbase, listmaxmparagenrscalbase)
    elif scaltype == 'gaus' or scaltype == 'eerr':
        listmeanparagenrscalbase = gmod.listmeanparagenrscalbase[indxparagenrbasescal]
        liststdvparagenrscalbase = gmod.liststdvparagenrscalbase[indxparagenrbasescal]
        if scaltype == 'eerr':
            cdfnminmparagenrscalbaseunit = gmod.cdfnminmparagenrscalbaseunit[indxparagenrbasescal]
            listparagenrscalbaseunitdiff = gmod.listparagenrscalbaseunitdiff[indxparagenrbasescal]
            paragenrscalbase = tdpy.icdf_eerr(paragenrunitbase, listmeanparagenrscalbase, liststdvparagenrscalbase, cdfnminmparagenrscalbaseunit, listparagenrscalbaseunitdiff)
        else:
            paragenrscalbase = tdpy.icdf_gaus(paragenrunitbase, listmeanparagenrscalbase, liststdvparagenrscalbase)
    elif scaltype == 'pois':
        paragenrscalbase = paragenrunitbase
    
    if gdat.booldiagmode:
        if not np.isfinite(paragenrscalbase).all():
            print('scaltype')
            print(scaltype)
            print('paragenrscalbase')
            print(paragenrscalbase)
            print('type(paragenrscalbase)')
            print(type(paragenrscalbase))
            print('paragenrscalbase.dtype')
            print(paragenrscalbase.dtype)
            raise Exception('')

    return paragenrscalbase


def icdf_trap(gdat, strgmodl, cdfn, paragenrscalfull, scalcomp, nameparagenrelem, l):
    
    gmod = getattr(gdat, strgmodl)
    
    if scalcomp == 'self' or scalcomp == 'powr' or scalcomp == 'dpowslopbrek' or scalcomp == 'logt':
        minm = getattr(gmod.minmpara, nameparagenrelem)
    
    if scalcomp != 'self':
        maxm = getattr(gmod.maxmpara, nameparagenrelem)
    
    if scalcomp == 'powr':
        slop = paragenrscalfull[getattr(gmod.indxpara, 'slopprio%spop%d' % (nameparagenrelem, l))]
        if gdat.booldiagmode:
            if not np.isfinite(slop):
                raise Exception('')
            if maxm < minm:
                raise Exception('')
        icdf = tdpy.icdf_powr(cdfn, minm, maxm, slop)

    if scalcomp == 'dpowslopbrek':
        distbrek = paragenrscalfull[getattr(gmod.indxpara, 'brekprio' + nameparagenrelem)[l]]
        sloplowr = paragenrscalfull[getattr(gmod.indxpara, 'sloplowrprio' + nameparagenrelem)[l]]
        slopuppr = paragenrscalfull[getattr(gmod.indxpara, 'slopupprprio' + nameparagenrelem)[l]]
        icdf = tdpy.icdf_dpow(cdfn, minm, maxm, distbrek, sloplowr, slopuppr)
    
    if scalcomp == 'expo':
        sexp = getattr(gmod, nameparagenrelem + 'distsexppop%d' % l)
        icdf = tdpy.icdf_expo(cdfn, maxm, sexp)
    
    if scalcomp == 'self':
        fact = getattr(gmod.factpara, nameparagenrelem)
        icdf = tdpy.icdf_self_fact(cdfn, minm, fact)
    
    if scalcomp == 'logt':
        icdf = tdpy.icdf_logt(cdfn, minm, fact)
    
    if scalcomp == 'dexp':
        scal = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distscal')[l]]
        icdf = tdpy.icdf_dexp(cdfn, maxm, scal)
    
    if scalcomp == 'lnormeanstdv':
        distmean = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distmean')[l]]
        diststdv = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'diststdv')[l]]
        icdf = tdpy.icdf_lnor(cdfn, distmean, diststdv)
    
    if scalcomp == 'igam':
        slop = paragenrscalfull[getattr(gmod.indxpara, 'slopprio' + nameparagenrelem)[l]]
        cutf = getattr(gdat, 'cutf' + nameparagenrelem)
        icdf = tdpy.icdf_igam(cdfn, slop, cutf)
    
    if scalcomp == 'gaus':
        distmean = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distmean')[l]]
        diststdv = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'diststdv')[l]]
        icdf = tdpy.icdf_gaus(cdfn, distmean, diststdv)
    
    if gdat.booldiagmode:
        if not np.isfinite(icdf).all():
            print('icdf')
            print(icdf)
            raise Exception('')

    return icdf


def cdfn_trap(gdat, gdatmodi, strgmodl, icdf, indxpoplthis):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    
    gmod.listscalparagenrelem = gmod.listscalparagenrelem[indxpoplthis]
    cdfn = np.empty_like(icdf)
    for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[indxpoplthis]):
        
        if gmod.listscalparagenrelem[k] == 'self' or gmod.listscalparagenrelem[k] == 'dexp' or gmod.listscalparagenrelem[k] == 'expo' \
                                                                        or gmod.listscalparagenrelem[k] == 'powr' or gmod.listscalparagenrelem[k] == 'dpowslopbrek':
            minm = getattr(gdat.fitt.minm, nameparagenrelem)
            if gmod.listscalparagenrelem[k] == 'powr':
                maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                slop = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'slop')[indxpoplthis]]
                cdfn[k] = cdfn_powr(icdf[k], minm, maxm, slop)
            elif gmod.listscalparagenrelem[k] == 'dpowslopbrek':
                maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                brek = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'distbrek')[indxpoplthis]]
                sloplowr = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'sloplowr')[indxpoplthis]]
                slopuppr = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'slopuppr')[indxpoplthis]]
                cdfn[k] = cdfn_dpow(icdf[k], minm, maxm, brek, sloplowr, slopuppr)
            else:
                fact = getattr(gdat.fitt, 'fact' + nameparagenrelem)
                cdfn[k] = cdfn_self(icdf[k], minm, fact)
        if gmod.listscalparagenrelem[k] == 'lnormeanstdv':
            distmean = gdatmodi.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'diststdv')[indxpoplthis]]
            cdfn[k] = cdfn_lnor(icdf[k], distmean, slop)
        if gmod.listscalparagenrelem[k] == 'igam':
            slop = gdatmodi.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'slop')[indxpoplthis]]
            cutf = getattr(gdat, 'cutf' + nameparagenrelem)
            cdfn[k] = cdfn_igam(icdf[k], slop, cutf)
        if gmod.listscalparagenrelem[k] == 'gaus':
            distmean = gdatmodi.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'diststdv')[indxpoplthis]]
            cdfn[k] = cdfn_gaus(icdf[k], distmean, diststdv)
    
    return cdfn


### update sampler state
def updt_stat(gdat, gdatmodi):
   
    if gdat.typeverb > 1:
        print('updt_stat()')
    
    # update the sample and the unit sample vectors
    gdatmodi.this.lpritotl = gdatmodi.next.lpritotl
    gdatmodi.this.lliktotl = gdatmodi.next.lliktotl
    gdatmodi.this.lpostotl = gdatmodi.next.lpostotl
    gdatmodi.this.paragenrscalfull[gdatmodi.indxsampmodi] = np.copy(gdatmodi.next.paragenrscalfull[gdatmodi.indxsampmodi])
    gdatmodi.this.paragenrunitfull[gdatmodi.indxsampmodi] = np.copy(gdatmodi.next.paragenrunitfull[gdatmodi.indxsampmodi])
    if gdatmodi.this.indxproptype > 0:
        gdatmodi.this.indxelemfull = deepcopy(gdatmodi.next.indxelemfull)
        gdatmodi.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gdatmodi.this.indxelemfull, 'fitt')


def initcompfromstat(gdat, gdatmodi, namerefr):
    
    for l in gmod.indxpopl:
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
            minm = getattr(gdat.fitt.minmpara, nameparagenrelem)
            maxm = getattr(gdat.fitt.maxmpara, nameparagenrelem)
            try:
                comp = getattr(gdat, namerefr + nameparagenrelem)[l][0, :]
                if gmod.listscalparagenrelem[l][g] == 'self' or gmod.listscalparagenrelem[l][g] == 'logt':
                    fact = getattr(gdat.fitt, 'fact' + nameparagenrelem)
                    if gmod.listscalparagenrelem[l][g] == 'self':
                        compunit = cdfn_self(comp, minm, fact)
                    if gmod.listscalparagenrelem[l][g] == 'logt':
                        compunit = cdfn_logt(comp, minm, fact)
                if gmod.listscalparagenrelem[l][g] == 'expo':
                    scal = getattr(gdat.fitt, 'gangdistsexp')
                    maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                    compunit = cdfn_expo(icdf, maxm, scal)
                if gmod.listscalparagenrelem[l][g] == 'powr' or gmod.listscalparagenrelem[l][g] == 'igam':
                    slop = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'slop')[l]]
                    if gmod.listscalparagenrelem[l][g] == 'powr':
                        compunit = cdfn_powr(comp, minm, maxm, slop)
                    if gmod.listscalparagenrelem[l][g] == 'igam':
                        cutf = getattr(gdat, 'cutf' + nameparagenrelem)
                        compunit = cdfn_igam(comp, slop, cutf)
                if gmod.listscalparagenrelem[l][g] == 'dpowslopbrek':
                    brek = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'distbrek')[l]]
                    sloplowr = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'sloplowr')[l]]
                    slopuppr = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'slopuppr')[l]]
                    compunit = cdfn_powr(comp, minm, maxm, brek, sloplowr, slopuppr)
                if gmod.listscalparagenrelem[l][g] == 'gaus':
                    distmean = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'distmean')[l]]
                    diststdv = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt, 'indxparagenrbase' + nameparagenrelem + 'diststdv')[l]]
                    compunit = cdfn_gaus(comp, distmean, diststdv)
            except:
                if gdat.typeverb > 0:
                    print('Initialization from the reference catalog failed for %s. Sampling randomly...' % nameparagenrelem)
                compunit = np.random.rand(gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int))
            gdatmodi.this.paragenrunitfull[gdatmodi.this.indxparagenrfullelem[l][nameparagenrelem]] = compunit


### find the set of pixels in proximity to a position on the map
def retr_indxpixlelemconc(gdat, strgmodl, dictelem, l):
    
    gmod = getattr(gdat, strgmodl)
    
    lgal = dictelem[l]['lgal']
    bgal = dictelem[l]['bgal']
    varbampl = dictelem[l][gmod.nameparagenrelemampl[l]]
    
    if gmod.typeelemspateval[l] == 'locl':
        listindxpixlelem = [[] for k in range(lgal.size)]
        for k in range(lgal.size):
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
            
            indxfluxproxtemp = np.digitize(varbampl[k], gdat.binspara.prox)
            if indxfluxproxtemp > 0:
                indxfluxproxtemp -= 1
            if indxfluxproxtemp == gdat.binspara.prox.size - 1:
                print('Warning! Index of the proximity pixel list overflew. Taking the largest list...')
                indxfluxproxtemp -= 1
            indxpixlelem = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
            if isinstance(indxpixlelem, int):
                indxpixlelem = gdat.indxpixl
            listindxpixlelem[k] = indxpixlelem

        listindxpixlelemconc = np.unique(np.concatenate(listindxpixlelem))
    else:
        listindxpixlelemconc = gdat.indxpixl
        listindxpixlelem = gdat.indxpixl
    
    return listindxpixlelem, listindxpixlelemconc


### find the distance between two points on the map
def retr_angldistunit(gdat, lgal, bgal, indxpixlelem, retranglcosi=False):
   
    if gdat.typepixl == 'heal':
        xdat, ydat, zaxi = retr_unit(lgal, bgal)
        anglcosi = gdat.xdatgrid[indxpixlelem] * xdat + gdat.ydatgrid[indxpixlelem] * ydat + gdat.zaxigrid[indxpixlelem] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = np.arccos(anglcosi)
            return angldist
    
    else:
        angldist = np.sqrt((lgal - gdat.lgalgrid[indxpixlelem])**2 + (bgal - gdat.bgalgrid[indxpixlelem])**2)
        
        return angldist
    

### find the pixel index of a point on the map
def retr_indxpixl(gdat, bgal, lgal):

    if gdat.typepixl == 'heal':
        indxpixl = gdat.pixlcnvt[hp.ang2pix(gdat.numbsideheal, np.pi / 2. - bgal, lgal)]
        if gdat.booldiagmode:
            if (indxpixl == -1).any():  
                raise Exception('pixlcnvt went negative!')

    if gdat.typepixl == 'cart':
        indxlgcr = np.floor(gdat.numbsidecart * (lgal - gdat.minmlgaldata) / 2. / gdat.maxmgangdata).astype(int)
        indxbgcr = np.floor(gdat.numbsidecart * (bgal - gdat.minmbgaldata) / 2. / gdat.maxmgangdata).astype(int)

        if np.isscalar(indxlgcr):
            if indxlgcr < 0:
                indxlgcr = 0
            if indxlgcr >= gdat.numbsidecart:
                indxlgcr = gdat.numbsidecart - 1
        else:
            indxlgcr[np.where(indxlgcr < 0)] = 0
            indxlgcr[np.where(indxlgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        if np.isscalar(indxbgcr):
            if indxbgcr < 0:
                indxbgcr = 0
            if indxbgcr >= gdat.numbsidecart:
                indxbgcr = gdat.numbsidecart - 1
        else:
            indxbgcr[np.where(indxbgcr < 0)] = 0
            indxbgcr[np.where(indxbgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        indxpixl = indxlgcr * gdat.numbsidecart + indxbgcr
    
    # convert to an index of non-zero exposure pixels
    #indxpixl = gdat.indxpixlroficnvt[indxpixl]

    return indxpixl


## obtain count maps
def retr_cntp(gdat, sbrt):
   
    cntp = sbrt * gdat.expo * gdat.apix
    if gdat.enerdiff:
        cntp *= gdat.deltener[:, None, None] 
        
    return cntp


## plotting
### construct path for plots
def retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, strgplot, nameinte=''):
    
    if strgmodl == 'true' or strgstat == '':
        path = gdat.pathinit + nameinte + strgplot + '.pdf'
    elif strgstat == 'pdfn' or strgstat == 'mlik':
        path = gdat.pathplotrtag + strgpdfn + '/finl/' + nameinte + strgstat + strgplot + '.pdf'
    elif strgstat == 'this':
        path = gdat.pathplotrtag + strgpdfn + '/fram/' + nameinte + strgstat + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


### determine the marker size
def retr_mrkrsize(gdat, strgmodl, compampl, nameparagenrelemampl):
    
    gmod = getattr(gdat, strgmodl)
    minm = getattr(gdat.minmpara, nameparagenrelemampl) 
    maxm = getattr(gdat.maxmpara, nameparagenrelemampl)
    mrkrsize = (np.sqrt(compampl) - np.sqrt(minm)) / (np.sqrt(maxm) - np.sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


## experiment specific
def retr_psfphubb(gmod):

    # temp
    gmod.psfpexpr = np.array([0.080, 0.087]) / gdat.anglfact


def retr_psfpchan(gmod):

    # temp
    #gmod.psfpexpr = np.array([0.25, 0.3, 0.4, 0.6, 0.7]) / gdat.anglfact
    if gdat.numbenerfull == 5:
        gmod.psfpexpr = np.array([0.424 / gdat.anglfact, 2.75, 0.424 / gdat.anglfact, 2.59, 0.440 / gdat.anglfact, 2.47, 0.457 / gdat.anglfact, 2.45, 0.529 / gdat.anglfact, 3.72])
    if gdat.numbenerfull == 2:
        gmod.psfpexpr = np.array([0.427 / gdat.anglfact, 2.57, 0.449 / gdat.anglfact, 2.49])
    #gdat.psfpchan = gmod.psfpexpr[(2 * gdat.indxenerincl[:, None] + np.arange(2)[None, :]).flatten()] 
    #gmod.psfpexpr = np.array([0.25 / gdat.anglfact, 
    #                       0.30 / gdat.anglfacti\
    #                       0.40 / gdat.anglfacti\
    #                       0.60 / gdat.anglfacti\
    #                       0.70 / gdat.anglfacti
    #gmod.psfpexpr = np.array([0.35 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.e-1, 2.])
    #gmod.psfpexpr = np.array([0.25 / gdat.anglfact, 2.0e-1, 1.9, \
    #                       0.30 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.40 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.60 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.70 / gdat.anglfact, 1.0e-1, 2.0])
   

def retr_psfpsdyn(gmod):

    gmod.psfpexpr = np.array([0.05])
   

def retr_psfpferm(gmod):
   
    if gdat.anlytype.startswith('rec8'):
        path = gdat.pathdata + 'expr/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
    irfn = astropy.io.fits.getdata(path, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = np.sqrt(minmener * maxmener)

    numbpsfpscal = 3
    numbpsfpform = 5
    
    fermscal = np.zeros((gdat.numbevtt, numbpsfpscal))
    fermform = np.zeros((gdat.numbener, gdat.numbevtt, numbpsfpform))
    
    strgpara = ['score', 'gcore', 'stail', 'gtail', 'ntail']
    for m in gdat.indxevtt:
        if gdat.anlytype.startswith('rec8'):
            irfn = astropy.io.fits.getdata(path, 1 + 3 * gdat.indxevttincl[m])
            fermscal[m, :] = astropy.io.fits.getdata(path, 2 + 3 * gdat.indxevttincl[m])['PSFSCALE']
        else:
            if m == 1:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 0:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
            else:
                continue
            irfn = astropy.io.fits.getdata(path, 1)
            fermscal[m, :] = astropy.io.fits.getdata(path, 2)['PSFSCALE']
        for k in range(numbpsfpform):
            fermform[:, m, k] = sp.interpolate.interp1d(enerirfn, np.mean(irfn[strgpara[k]].squeeze(), axis=0), fill_value='extrapolate')(gdat.meanpara.ener)
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = np.sqrt((fermscal[None, :, 0] * (10. * gdat.meanpara.ener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gmod.psfpexpr = np.zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            gmod.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
    

def retr_refrchaninit(gdat):
    
    gdat.indxrefr = np.arange(gdat.numbrefr)
    
    gdat.dictrefr = []
    for q in gdat.indxrefr:
        gdat.dictrefr.append(dict())
    
    gdat.refr.namepara.elemsign = ['flux', 'magt']
    
    gdat.refr.lablelem = ['Xue+2011', 'Wolf+2008']
    
    gdat.listnamerefr += ['xu11', 'wo08']
    
    setattr(gdat, 'plotminmotyp', 0.)
    setattr(gdat, 'plottmaxmotyp', 1.)
    setattr(gmod.lablrootpara, 'otyp', 'O')
    setattr(gdat, 'scalotypplot', 'self')
    
    setattr(gmod.lablrootpara, 'otypxu11', 'O')
    for name in gdat.listnamerefr:
        setattr(gdat, 'plotminmotyp' + name, 0.)
        setattr(gdat, 'plotmaxmotyp' + name, 1.)
    
    if gdat.strgcnfg == 'pcat_chan_inpt_home4msc':
        with open(gdat.pathinpt + 'ECDFS_Cross_ID_Hsu2014.txt', 'r') as thisfile:
            for k, line in enumerate(thisfile):
                if k < 18:
                    continue
                rasccand =line[2]
                declcand =line[2]
       
    gdat.refr.namepara.elem[0] += ['lgal', 'bgal', 'flux', 'sind', 'otyp', 'lumi']
    gdat.refr.namepara.elem[1] += ['lgal', 'bgal', 'magt', 'reds', 'otyp']


def retr_refrchanfinl(gdat):
    
    booltemp = False
    if gdat.anlytype.startswith('extr'):
        if gdat.numbsidecart == 300:
            gdat.numbpixllgalshft[0] = 1490
            gdat.numbpixlbgalshft[0] = 1430
        else:
            booltemp = True
    elif gdat.anlytype.startswith('home'):
        gdat.numbpixllgalshft[0] = 0
        gdat.numbpixlbgalshft[0] = 0
    
        if gdat.numbsidecart == 600:
            pass
        elif gdat.numbsidecart == 100:
            indxtile = int(gdat.anlytype[-4:])
            numbsidecntr = int(gdat.anlytype[8:12])
            numbtileside = numbsidecntr / gdat.numbsidecart
            indxtilexaxi = indxtile // numbtileside
            indxtileyaxi = indxtile % numbtileside
            gdat.numbpixllgalshft[0] += indxtilexaxi * gdat.numbsidecart
            gdat.numbpixlbgalshft[0] += indxtileyaxi * gdat.numbsidecart
        elif gdat.numbsidecart == 300:
            gdat.numbpixllgalshft[0] += 150
            gdat.numbpixlbgalshft[0] += 150
        else:
            booltemp = True
    else:
        booltemp = True

    if booltemp:
        raise Exception('Reference elements cannot be aligned with the spatial axes!')
    
    ## WCS object for rotating reference elements into the ROI
    if gdat.numbener == 2:
        gdat.listpathwcss[0] = gdat.pathinpt + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    else:
        gdat.listpathwcss[0] = gdat.pathinpt + '0.5-0.91028_flux_%sMs.img' % gdat.anlytype[4]
    
    # Xue et al. (2011)
    #with open(gdat.pathinpt + 'chancatl.txt', 'r') as thisfile:
    pathfile = gdat.pathinpt + 'Xue2011.fits'
    hdun = pf.open(pathfile)
    hdun.info()
    lgalchan = hdun[1].data['_Glon'] / 180. * pi
    bgalchan = hdun[1].data['_Glat'] / 180. * pi
    fluxchansoft = hdun[1].data['SFlux']
    fluxchanhard = hdun[1].data['HFlux']
    objttypechan = hdun[1].data['Otype']
    gdat.refrlumi[0][0] = hdun[1].data['Lx']
    
    # position
    gdat.refr.dictelem[0]['lgal'] = lgalchan
    gdat.refr.dictelem[0]['bgal'] = bgalchan

    # spectra
    gdat.refrspec = [[np.zeros((3, gdat.numbener, lgalchan.size))]]
    if gdat.numbener == 2:
        gdat.refrspec[0][0, 0, :] = fluxchansoft * 0.624e9
        gdat.refrspec[0][0, 1, :] = fluxchanhard * 0.624e9 / 16.
    else:
        gdat.refrspec[0][0, :, :] = 2. * fluxchansoft[None, :] * 0.624e9
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :]
   
    # fluxes
    gdat.refrflux[0] = gdat.refrspec[0][:, gdat.indxenerpivt, :]

    # spectral indices
    if gdat.numbener > 1:
        gdat.refrsind[0] = -np.log(gdat.refrspec[0][0, 1, :] / gdat.refrspec[0][0, 0, :]) / np.log(np.sqrt(7. / 2.) / np.sqrt(0.5 * 2.))

    ## object type
    objttypechantemp = np.zeros(lgalchan.size) - 1.
    indx = np.where(objttypechan == 'AGN')[0]
    objttypechantemp[indx] = 0.165
    indx = np.where(objttypechan == 'Galaxy')[0]
    objttypechantemp[indx] = 0.495
    indx = np.where(objttypechan == 'Star')[0]
    objttypechantemp[indx] = 0.835
    gdat.refrotyp[0][0] = objttypechantemp

    # Wolf et al. (2011)
    path = gdat.pathdata + 'inpt/Wolf2008.fits'
    data = astropy.io.fits.getdata(path)
    gdat.refrlgal[1] = np.deg2rad(data['_Glon'])
    gdat.refrlgal[1] = ((gdat.refrlgal[1] - pi) % (2. * pi)) - pi
    gdat.refrbgal[1] = np.deg2rad(data['_Glat'])
    gdat.refrmagt[1][0] = data['Rmag']
    gdat.refrreds[1][0] = data['MCz']
  
    #listname = []
    #for k in range(data['MCclass'].size):
    #    if not data['MCclass'][k] in listname:
    #        listname.append(data['MCclass'][k])
    listname = ['Galaxy', 'Galaxy  (Uncl!)', 'QSO     (Gal?)', 'Galaxy  (Star?)', 'Star', 'Strange Object', 'QSO', 'WDwarf']
    gdat.refrotyp[1][0] = np.zeros_like(gdat.refrreds[1][0]) - 1. 
    for k, name in enumerate(listname):
        indx = np.where(data['MCclass'] == name)[0]
        gdat.refrotyp[1][0][indx] = k / 10.
    
    # error budget
    for name in ['lgal', 'bgal', 'sind', 'otyp', 'lumi', 'magt', 'reds']:
        refrtile = [[] for q in gdat.indxrefr]
        refrfeat = getattr(gdat.refr, name)
        for q in gdat.indxrefr:
            if len(refrfeat[q]) > 0:
                refrtile[q] = np.tile(refrfeat[q], (3, 1))
        setattr(gdat.refr, name, refrtile)
        

def retr_refrferminit(gdat):
    
    gdat.listnamerefr += ['ac15', 'ma05']
    gdat.indxrefr = np.arange(gdat.numbrefr)
    
    gdat.refr.lablelem = ['Acero+2015', 'Manchester+2005']

    gdat.refr.namepara.elemsign = ['flux', 'flux0400']
    
    setattr(gmod.lablrootpara, 'curvac15', '%s_{3FGL}' % gdat.lablcurv)
    setattr(gmod.lablrootpara, 'expcac15', 'E_{c,3FGL}')
    
    for name in gdat.listnamerefr:
        setattr(gdat.minmpara, 'curv' + name, -1.)
        setattr(gdat.maxmpara, 'curv' + name, 1.)
        setattr(gdat.minmpara, 'expc' + name, 0.1)
        setattr(gdat.maxmpara, 'expc' + name, 10.)
   
    gdat.refr.namepara.elem[0] += ['lgal', 'bgal', 'flux', 'sind', 'curv', 'expc', 'tvar', 'etag', 'styp', 'sindcolr0001', 'sindcolr0002']
    gdat.refr.namepara.elem[1] += ['lgal', 'bgal', 'flux0400', 'per0', 'per1']


def retr_refrfermfinl(gdat):

    gdat.minmstyp = -0.5
    gdat.maxmstyp = 3.5
    gdat.lablstyp = 'S'
    gmod.scalstypplot = 'self'
    
    gdat.minmtvar = 0.
    gdat.maxmtvar = 400.
    gdat.labltvar = 'T'
    gmod.scaltvarplot = 'logt'
    
    # Acero+2015
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = astropy.io.fits.getdata(path)
    
    gdat.refr.dictelem[0]['lgal'] = np.deg2rad(fgl3['glon'])
    gdat.refr.dictelem[0]['lgal'] = np.pi - ((gdat.refr.dictelem[0]['lgal'] - np.pi) % (2. * np.pi))
    gdat.refr.dictelem[0]['bgal'] = np.deg2rad(fgl3['glat'])
    
    gdat.refr.numbelemfull = gdat.refr.dictelem[0]['lgal'].size

    gdat.refrspec = [np.empty((3, gdat.numbener, gdat.refr.dictelem[0]['lgal'].size))]
    gdat.refrspec[0][0, :, :] = np.stack((fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    
    fgl3specstdvtemp = np.stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.refrspec[0][np.where(np.isfinite(gdat.refrspec[0]) == False)] = 0.
    
    gdat.refrflux[0] = gdat.refrspec[0][:, gdat.indxenerpivt, :]
    gdat.refrsindcolr0001[0] = -np.log(gdat.refrspec[0][:, 1, :] / gdat.refrflux[0]) / np.log(gdat.meanpara.ener[1] / gdat.enerpivt)
    gdat.refrsindcolr0002[0] = -np.log(gdat.refrspec[0][:, 2, :] / gdat.refrflux[0]) / np.log(gdat.meanpara.ener[2] / gdat.enerpivt)
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = np.deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(np.cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(np.sin(fgl3anglstdv))

    gdat.refretag[0] = np.zeros(gdat.refr.dictelem[0]['lgal'].size, dtype=object)
    for k in range(gdat.refr.dictelem[0]['lgal'].size):
        gdat.refretag[0][k] = '%s, %s, %s' % (fgl3['Source_Name'][k], fgl3['CLASS1'][k], fgl3['ASSOC1'][k])
    gdat.refrtvar[0] = fgl3['Variability_Index']
    
    gdat.refrstyp[0] = np.zeros_like(gdat.refr.dictelem[0]['lgal']) - 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PowerLaw        ')] = 0
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'LogParabola     ')] = 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLExpCutoff     ')] = 2
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLSuperExpCutoff')] = 3
    indx = np.where(gdat.refrstyp[0] == -1)[0]
    if indx.size > 0:
        raise Exception('')
    gdat.refrsind[0] = fgl3['Spectral_Index']
    gdat.refrcurv[0] = fgl3['beta']
    gdat.refrexpc[0] = fgl3['Cutoff'] * 1e-3
    
    gdat.refrcurv[0][np.where(np.logical_not(np.isfinite(gdat.refrcurv[0])))] = -10.
    gdat.refrexpc[0][np.where(np.logical_not(np.isfinite(gdat.refrexpc[0])))] = 0.
    
    gdat.refrsind[0] = np.tile(gdat.refrsind[0], (3, 1)) 
    gdat.refrcurv[0] = np.tile(gdat.refrcurv[0], (3, 1)) 
    gdat.refrexpc[0] = np.tile(gdat.refrexpc[0], (3, 1)) 

    # Manchester+2005
    path = gdat.pathdata + 'inpt/Manchester2005.fits'
    data = astropy.io.fits.getdata(path)
   
    gdat.refrlgal[1] = np.deg2rad(data['glon'])
    gdat.refrlgal[1] = ((gdat.refrlgal[1] - np.pi) % (2. * np.pi)) - np.pi
    gdat.refrbgal[1] = np.deg2rad(data['glat'])
    
    gdat.refrper0[1] = data['P0']
    gdat.refrper1[1] = data['P1']
    gdat.refrflux0400[1] = data['S400']
    #gdat.refrdism[1] = data['DM']
    #gdat.refrdlos[1] = data['Dist']

    # error budget
    for name in ['lgal', 'bgal', 'per0', 'per1', 'flux0400', 'tvar', 'styp']:
        refrtile = [[] for q in gdat.indxrefr]
        refrfeat = getattr(gdat.refr, name)
        for q in gdat.indxrefr:
            if len(refrfeat[q]) > 0:
                refrtile[q] = np.tile(refrfeat[q], (3, 1))
        setattr(gdat.refr, name, refrtile)


def retr_singgaus(scaldevi, sigc):
    
    psfn = 1. / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_singking(scaldevi, sigc, gamc):
    
    psfn = 1. / 2. / np.pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc)

    return psfn


def retr_doubgaus(scaldevi, frac, sigc, sigt):
    
    psfn = frac / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_gausking(scaldevi, frac, sigc, sigt, gamt):

    psfn = frac / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / np.pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / np.pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) + \
                            (1. - frac) / 2. / np.pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_lgalbgal(gang, aang):
    
    lgal = gang * np.cos(aang)
    bgal = gang * np.sin(aang)

    return lgal, bgal


def retr_gang(lgal, bgal):
    
    gang = np.arccos(np.cos(lgal) * np.cos(bgal))

    return gang


def retr_aang(lgal, bgal):

    aang = np.arctan2(bgal, lgal)

    return aang


def show_paragenrscalfull(gdat, gdatmodi, strgstat='this', strgmodl='fitt', indxsampshow=None):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    print('strgmodl: ' + strgmodl)
    print('strgstat: ' + strgstat)
    print('%5s %20s %30s %30s %15s' % ('index', 'namepara', 'paragenrunitfull', 'paragenrscalfull', 'scalpara'))
    for k in gmod.indxparagenrfull:
        
        if indxsampshow is not None and not k in indxsampshow:
            continue
        
        if gmod.numbparaelem > 0:
            
            booltemp = False
            for l in gmod.indxpopl:
                if k == gmod.indxparagenrelemsing[l][0]:
                    booltemp = True
            if booltemp:
                print('')
        print('%5d %20s %30g %30g %15s' % (k, gmod.namepara.genrfull[k], gmodstat.paragenrunitfull[k], gmodstat.paragenrscalfull[k], gmod.scalpara.genrfull[k]))
    

def prop_stat(gdat, gdatmodi, strgmodl, thisindxelem=None, thisindxpopl=None, brth=False, deth=False):
 
    if gdat.typeverb > 1:
        print('prop_stat()')
    
    #indxproptype
    # within, birth, death, split, merge
    # 0, 1, 2, 3, 4
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodthis = getattr(gdatobjt, 'this')
    gmodnext = getattr(gdatobjt, 'next')
    
    if gmod.numbparaelem > 0:
        if gdat.booldiagmode:
            for l in gmod.indxpopl:
                if len(gmodthis.indxelemfull[l]) > len(set(gmodthis.indxelemfull[l])):
                    raise Exception('Repeating entry in the element index list!')

        thisindxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmodthis.indxelemfull, strgmodl)
        setattr(gmodthis, 'indxparagenrfullelem', thisindxparagenrfullelem)
    else:
        thisindxparagenrfullelem = None
    
    gdatmodi.this.boolpropfilt = True 

    # index of the population in which a transdimensional proposal will be attempted
    if gmod.numbparaelem > 0:
        if thisindxpopl is None:
            gdatmodi.indxpopltran = np.random.choice(gmod.indxpopl)
        else:
            gdatmodi.indxpopltran = thisindxpopl
        numbelemtemp = gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]]
    
    # forced death or birth does not check for the prior on the dimensionality on purpose!
    if gmod.numbparaelem > 0 and (deth or brth or np.random.rand() < gdat.probtran) and \
                        not (numbelemtemp == gmod.minmpara.numbelem[gdatmodi.indxpopltran] and numbelemtemp == gmod.maxmpara.numbelem[gdatmodi.indxpopltran]):

        if brth or deth or np.random.rand() < gdat.probbrde or \
                            numbelemtemp == gmod.maxmpara.numbelem[gdatmodi.indxpopltran] and numbelemtemp == 1 or numbelemtemp == 0:
            
            ## births and deaths
            if numbelemtemp == gmod.maxmpara.numbelem[gdatmodi.indxpopltran] or deth:
                gdatmodi.this.indxproptype = 2
            elif numbelemtemp == gmod.minmpara.numbelem[gdatmodi.indxpopltran] or brth:
                gdatmodi.this.indxproptype = 1
            else:
                if np.random.rand() < 0.5:
                    gdatmodi.this.indxproptype = 1
                else:
                    gdatmodi.this.indxproptype = 2

        else:
            ## splits and merges
            if numbelemtemp == gmod.minmpara.numbelem[gdatmodi.indxpopltran] or numbelemtemp < 2:
                gdatmodi.this.indxproptype = 3
            elif numbelemtemp == gmod.maxmpara.numbelem[gdatmodi.indxpopltran]:
                gdatmodi.this.indxproptype = 4
            else:
                if np.random.rand() < 0.5:
                    gdatmodi.this.indxproptype = 3
                else:
                    gdatmodi.this.indxproptype = 4
    else:
        
        if gdat.booldiagmode and (gdatmodi.stdp > 1e2).any():
            raise Exception('')

        thisindxparagenrfullelemconc = []
        for l in gmod.indxpopl:
            thisindxparagenrfullelemconc.append(thisindxparagenrfullelem[l]['full'])

        # get the indices of the current parameter vector
        if gmod.numbparaelem > 0:
            thisindxsampfull = np.concatenate([gmod.indxparagenrbasestdv] + thisindxparagenrfullelemconc)
        else:
            thisindxsampfull = gmod.indxparagenrbasestdv
        
        thisstdp = gdatmodi.stdp[gdat.indxstdppara[thisindxsampfull]]
        if not np.isfinite(thisstdp).all():
            raise Exception('')
        gdatmodi.this.indxproptype = 0
    
    if gdat.booldiagmode and gdat.probspmr == 0 and gdatmodi.this.indxproptype > 2:
        raise Exception('')

    if gdat.typeverb > 1:
        print('gdatmodi.this.indxproptype')
        print(gdatmodi.this.indxproptype)

    if gdatmodi.this.indxproptype == 0:
        gmodnext.paragenrunitfull = np.copy(gmodthis.paragenrunitfull)
        if gmod.numbparaelem > 0:
            gmodnext.indxelemfull = gmodthis.indxelemfull
    if gdatmodi.this.indxproptype > 0:
        gmodnext.paragenrunitfull = np.copy(gmodthis.paragenrunitfull)
        gmodnext.paragenrscalfull = np.copy(gmodthis.paragenrscalfull)
        if gmod.numbparaelem > 0:
            gmodnext.indxelemfull = deepcopy(gmodthis.indxelemfull)
    
    if gdatmodi.this.indxproptype == 0:
        
        ## proposal scale
        if False:
            # amplitude-dependent proposal scale
            for l in gmod.indxpopl:
                thiscompampl = gmodthis.paragenrscalfull[thisindxparagenrfullelem[indxelemfull][gmod.nameparagenrelemampl[l]][l]]
                compampl = gmodnext.paragenrscalfull[thisindxparagenrfullelem[gmod.nameparagenrelemampl[l]][l][indxelemfull]]
                minmcompampl = getattr(gmod.minmpara, gmod.nameparagenrelemampl[l])
                thiscompunit = gmodthis.paragenrscalfull[thisindxparagenrfullelem[gmod.nameparagenrelemampl[l]][l][indxelemfull]]
                compunit = gmodnext.paragenrscalfull[thisindxparagenrfullelem[gmod.nameparagenrelemampl[l]][l][indxelemfull]]
                if nameparagenrelem == gmod.nameparagenrelemampl[l]:
                    # temp -- this only works if compampl is powr distributed
                    gdatmodi.this.stdp = stdpcomp / (thiscompampl / minmcompampl)**2.
                    gdatmodi.this.stdv = stdpcomp / (compampl / minmcompampl)**2.
                    gdatmodi.this.ltrp += np.sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.this.stdv**2 - 1. / gdatmodi.this.stdv**2))
                else:
                    gdatmodi.this.stdp = stdpcomp / (np.minimum(thiscompampl, compampl) / minmcompampl)**0.5
        
        ## propose a step
        diffparagenrunitfull = np.random.normal(size=thisindxsampfull.size) * thisstdp
        gmodnext.paragenrunitfull[thisindxsampfull] = gmodthis.paragenrunitfull[thisindxsampfull] + diffparagenrunitfull
        
        if gdat.booldiagmode:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(gmodnext.paragenrunitfull).all():
                raise Exception('')

        indxsamplowr = np.where(gmodnext.paragenrunitfull[gmod.numbpopl:] < 0.)[0]
        if indxsamplowr.size > 0:
            gmodnext.paragenrunitfull[gmod.numbpopl+indxsamplowr] = abs(gmodnext.paragenrunitfull[gmod.numbpopl+indxsamplowr]) % 1.
        
        if gdat.booldiagmode:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

        indxsampuppr = np.where(gmodnext.paragenrunitfull[gmod.numbpopl:] > 1.)[0]
        if indxsampuppr.size > 0:
            gmodnext.paragenrunitfull[gmod.numbpopl+indxsampuppr] = (gmodnext.paragenrunitfull[gmod.numbpopl+indxsampuppr] - 1.) % 1.
        
        if gdat.booldiagmode:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(gmodnext.paragenrunitfull).all():
                raise Exception('')

        gmodnext.paragenrscalfull = icdf_paragenrscalfull(gdat, strgmodl, gmodnext.paragenrunitfull, thisindxparagenrfullelem)

        if gdat.booldiagmode:
            if not np.isfinite(gmodnext.paragenrunitfull).all():
                raise Exception('')
        
            if np.amin(gmodnext.paragenrunitfull[gmod.numbpopl:]) < 0.:
                raise Exception('')
        
            if np.amax(gmodnext.paragenrunitfull[gmod.numbpopl:]) > 1.:
                raise Exception('')
        
            if not np.isfinite(gmodnext.paragenrscalfull).all():
                raise Exception('')
        
    if gdatmodi.this.indxproptype > 0:
        gdatmodi.indxsamptran = []
        if gdatmodi.this.indxproptype == 1:
            gdatmodi.this.auxipara = np.random.rand(gmod.numbparagenrelemsing[gdatmodi.indxpopltran])
        elif gdatmodi.this.indxproptype != 2:
            gdatmodi.this.auxipara = np.empty(gmod.numbparagenrelemsing[gdatmodi.indxpopltran])
    
    if gdatmodi.this.indxproptype == 1 or gdatmodi.this.indxproptype == 3:
       
        # find an empty slot in the element list
        for u in range(gmod.maxmpara.numbelem[gdatmodi.indxpopltran]):
            if not u in gdatmodi.this.indxelemfull[gdatmodi.indxpopltran]:
                break
        gdatmodi.indxelemmodi = [u]
        gdatmodi.indxelemfullmodi = [gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]].astype(int)]
       
        # sample indices to add the new element
        gdatmodi.indxparagenrfullelemaddd = retr_indxparaelem(gmod, gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.append(gdatmodi.indxparagenrfullelemaddd)
        gmodnext.indxelemfull[gdatmodi.indxpopltran].append(gdatmodi.indxelemmodi[0])
    if gdatmodi.this.indxproptype == 1:
        
        # sample auxiliary variables
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0]] = gdatmodi.this.auxipara
    
    # death
    if gdatmodi.this.indxproptype == 2:
        
        # occupied element index to be killed
        if thisindxelem is None:
            dethindxindxelem = np.random.choice(np.arange(gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]], dtype=int))
        else:
            dethindxindxelem = thisindxelem

        # element index to be killed
        gdatmodi.indxelemmodi = []
        gdatmodi.indxelemfullmodi = []
        if gdat.typeverb > 1:
            print('dethindxindxelem')
            print(dethindxindxelem)

        gdatmodi.indxelemmodi.append(gmodthis.indxelemfull[gdatmodi.indxpopltran][dethindxindxelem])
        gdatmodi.indxelemfullmodi.append(dethindxindxelem)
        # parameter indices to be killed
        indxparagenrfullelemdeth = retr_indxparaelem(gmod, gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.append(indxparagenrfullelemdeth)
        
        gdatmodi.this.auxipara = gmodthis.paragenrscalfull[indxparagenrfullelemdeth]

    if gdatmodi.this.indxproptype > 2:
        gdatmodi.comppare = np.empty(gmod.numbparagenrelemsing[gdatmodi.indxpopltran])
        gdatmodi.compfrst = np.empty(gmod.numbparagenrelemsing[gdatmodi.indxpopltran])
        gdatmodi.compseco = np.empty(gmod.numbparagenrelemsing[gdatmodi.indxpopltran])
    
    # split
    if gdatmodi.this.indxproptype == 3:
        
        # find the probability of splitting elements
        gdatmodi.indxelemfullsplt = np.random.choice(np.arange(gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]], dtype=int))
        gdatmodi.indxelemsplt = gmodthis.indxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullsplt]
        gdatmodi.indxelemfullmodi.insert(0, gdatmodi.indxelemfullsplt)
        gdatmodi.indxelemmodi.insert(0, gdatmodi.indxelemsplt)

        # sample indices for the first element
        gdatmodi.indxparagenrfullelemfrst = retr_indxparaelem(gmod, l, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.insert(0, gdatmodi.indxparagenrfullelemfrst)
        
        # sample indices for the second element
        gdatmodi.indxsampseco = gdatmodi.indxparagenrfullelemaddd
        
        # take the parent element parameters
        for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            gdatmodi.comppare[k] = np.copy(gmodthis.paragenrscalfull[thisindxparagenrfullelem[gdatmodi.indxpopltran][nameparagenrelem][gdatmodi.indxelemfullmodi[0]]])
        
        # draw the auxiliary parameters
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if gmod.boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.this.auxipara[g] = np.random.randn() * gdat.radispmr
            elif g == gmod.indxparagenrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.this.auxipara[g] = np.random.rand()
            else:
                gdatmodi.this.auxipara[g] = icdf_trap(gdat, strgmodl, np.random.rand(), gmodthis.paragenrscalfull, gmod.listscalparagenrelem[gdatmodi.indxpopltran][g], \
                                                                                                     gmod.namepara.genrelem[gdatmodi.indxpopltran][g], l)

        # determine the new parameters
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.this.auxipara[1]) * gdatmodi.this.auxipara[0]
        else:
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.this.auxipara[2]) * gdatmodi.this.auxipara[0]
            gdatmodi.compfrst[1] = gdatmodi.comppare[1] + (1. - gdatmodi.this.auxipara[2]) * gdatmodi.this.auxipara[1]
        gdatmodi.compfrst[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] = gdatmodi.this.auxipara[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] * \
                                                                                                        gdatmodi.comppare[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.this.auxipara[1] * gdatmodi.this.auxipara[0]
        else:
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.this.auxipara[2] * gdatmodi.this.auxipara[0]
            gdatmodi.compseco[1] = gdatmodi.comppare[1] - gdatmodi.this.auxipara[2] * gdatmodi.this.auxipara[1]
        gdatmodi.compseco[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] = (1. - gdatmodi.this.auxipara[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]) * \
                                                                                                        gdatmodi.comppare[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]
        for g in range(gmod.numbparagenrelemsing[gdatmodi.indxpopltran]):
            if not gmod.boolcompposi[gdatmodi.indxpopltran][g] and g != gmod.indxparagenrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.compfrst[g] = gdatmodi.comppare[g]
                gdatmodi.compseco[g] = gdatmodi.this.auxipara[g]
       
        # place the new parameters into the sample vector
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0]] = cdfn_trap(gdat, gdatmodi, strgmodl, gdatmodi.compfrst, gdatmodi.indxpopltran)
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0]] = gdatmodi.compfrst
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[1]] = cdfn_trap(gdat, gdatmodi, strgmodl, gdatmodi.compseco, gdatmodi.indxpopltran)
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[1]] = gdatmodi.compseco
        
        # check for prior boundaries
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            if np.fabs(gdatmodi.compfrst[0]) > gdat.maxmelin or np.fabs(gdatmodi.compseco[0]) > gdat.maxmelin:
                gdatmodi.this.boolpropfilt = False
        else:
            if np.fabs(gdatmodi.compfrst[0]) > maxmlgal or np.fabs(gdatmodi.compseco[0]) > maxmlgal or \
                                                                    np.fabs(gdatmodi.compfrst[1]) > maxmbgal or np.fabs(gdatmodi.compseco[1]) > maxmbgal:
                gdatmodi.this.boolpropfilt = False
        if gdatmodi.compfrst[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] < getattr(gmod.minmpara, gmod.nameparagenrelemampl[gdatmodi.indxpopltran]) or \
           gdatmodi.compseco[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] < getattr(gmod.minmpara, gmod.nameparagenrelemampl[gdatmodi.indxpopltran]):
            gdatmodi.this.boolpropfilt = False
        if gdat.typeverb > 1:
            if not gdatmodi.this.boolpropfilt:
                print('Rejecting the proposal due to a split that falls out of the prior...')
    
    if gdatmodi.this.indxproptype == 4:
        
        # determine the index of the primary element to be merged (in the full element list)
        gdatmodi.indxelemfullmergfrst = np.random.choice(np.arange(len(gmodthis.indxelemfull[gdatmodi.indxpopltran])))

        ## first element index to be merged
        gdatmodi.mergindxelemfrst = gmodthis.indxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullmergfrst]
         
        # find the probability of merging this element with the others 
        probmerg = retr_probmerg(gdat, gdatmodi, gmodthis.paragenrscalfull, thisindxparagenrfullelem, gdatmodi.indxpopltran, 'seco', typeelem=gmod.typeelem)
        
        indxelemfulltemp = np.arange(len(gmodthis.indxelemfull[gdatmodi.indxpopltran]))
        if gdat.booldiagmode:
            if indxelemfulltemp.size < 2:
                raise Exception('')
        gdatmodi.indxelemfullmergseco = np.random.choice(np.setdiff1d(indxelemfulltemp, np.array([gdatmodi.indxelemfullmergfrst])), p=probmerg)
        gdatmodi.indxelemfullmodi = np.sort(np.array([gdatmodi.indxelemfullmergfrst, gdatmodi.indxelemfullmergseco]))
        
        # parameters of the first element to be merged
        for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            ## first
            gdatmodi.compfrst[k] = gmodthis.paragenrscalfull[thisindxparagenrfullelem[gdatmodi.indxpopltran][nameparagenrelem][gdatmodi.indxelemfullmodi[0]]]
        
        # determine indices of the modified elements in the sample vector
        ## first element
        # temp -- this would not work for multiple populations !
        gdatmodi.indxparagenrfullelemfrst = retr_indxparaelem(gmod, l, gdatmodi.mergindxelemfrst)
        gdatmodi.indxsamptran.append(gdatmodi.indxparagenrfullelemfrst)

        ## second element index to be merged
        gdatmodi.mergindxelemseco = gmodthis.indxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullmergseco]
       
        ## second element
        gdatmodi.indxparagenrfullelemseco = retr_indxparaelem(gmod, l, gdatmodi.mergindxelemseco)
        gdatmodi.indxsamptran.append(gdatmodi.indxparagenrfullelemseco)
        
        # parameters of the elements to be merged
        for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            ## second
            gdatmodi.compseco[k] = gmodthis.paragenrscalfull[thisindxparagenrfullelem[gdatmodi.indxpopltran][nameparagenrelem][gdatmodi.indxelemfullmodi[1]]]

        # indices of the element to be merged
        gdatmodi.indxelemmodi = [gdatmodi.mergindxelemfrst, gdatmodi.mergindxelemseco]

        # auxiliary parameters
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.this.auxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
        else:
            gdatmodi.this.auxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
            gdatmodi.this.auxipara[1] = gdatmodi.compseco[1] - gdatmodi.compfrst[1]
        gdatmodi.this.auxipara[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] / \
                                        (gdatmodi.compfrst[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] + gdatmodi.compseco[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]) 
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if not gmod.boolcompposi[gdatmodi.indxpopltran][g] and g != gmod.indxparagenrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.this.auxipara[g] = gdatmodi.compseco[g]

        # merged element
        gdatmodi.comppare[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] + \
                                                                                                gdatmodi.compseco[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]
        if gdatmodi.comppare[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]] > getattr(gdat, 'maxm' + gmod.nameparagenrelemampl[gdatmodi.indxpopltran]):
            gdatmodi.this.boolpropfilt = False
            if gdat.typeverb > 1:
                print('Proposal rejected due to falling outside the prior.')
            return

        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.this.auxipara[1]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
        else:
            gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.this.auxipara[2]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
            gdatmodi.comppare[1] = gdatmodi.compfrst[1] + (1. - gdatmodi.this.auxipara[2]) * (gdatmodi.compseco[1] - gdatmodi.compfrst[1])
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if gmod.boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g] + (1. - gdatmodi.this.auxipara[gmod.indxparagenrelemampl[gdatmodi.indxpopltran]]) * \
                                                                                            (gdatmodi.compseco[g] - gdatmodi.compfrst[g])
            elif g == gmod.indxparagenrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g] + gdatmodi.compseco[g]
            else:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g]

        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0]] = cdfn_trap(gdat, gdatmodi, strgmodl, gdatmodi.comppare, gdatmodi.indxpopltran)
        gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0]] = gdatmodi.comppare

        # calculate the proposed list of pairs
        if gdat.typeverb > 1:
            print('mergindxfrst: ', gdatmodi.mergindxelemfrst)
            print('gdatmodi.indxelemfullmergfrst: ', gdatmodi.indxelemfullmergfrst)
            print('mergindxseco: ', gdatmodi.mergindxelemseco)
            print('gdatmodi.indxelemfullmergseco: ', gdatmodi.indxelemfullmergseco)
            print('indxparagenrfullelemfrst: ', gdatmodi.indxparagenrfullelemfrst)
            print('indxparagenrfullelemseco: ', gdatmodi.indxparagenrfullelemseco)

    if gdat.typeverb > 1 and (gdatmodi.this.indxproptype == 3 or gdatmodi.this.boolpropfilt and gdatmodi.this.indxproptype == 4):
        
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            print('elinfrst: ', gdatmodi.compfrst[0])
            print('amplfrst: ', gdatmodi.compfrst[1])
            print('elinseco: ', gdatmodi.compseco[0])
            print('amplseco: ', gdatmodi.compseco[1])
            print('elinpare: ', gdatmodi.comppare[0])
            print('fluxpare: ', gdatmodi.comppare[1])
            print('auxipara[0][0]: ', gdatmodi.this.auxipara[0])
            print('auxipara[0][1]: ', gdatmodi.this.auxipara[1])
        else:
            print('lgalfrst: ', gdat.anglfact * gdatmodi.compfrst[0])
            print('bgalfrst: ', gdat.anglfact * gdatmodi.compfrst[1])
            print('amplfrst: ', gdatmodi.compfrst[2])
            print('lgalseco: ', gdat.anglfact * gdatmodi.compseco[0])
            print('bgalseco: ', gdat.anglfact * gdatmodi.compseco[1])
            print('amplseco: ', gdatmodi.compseco[2])
            print('lgalpare: ', gdat.anglfact * gdatmodi.comppare[0])
            print('bgalpare: ', gdat.anglfact * gdatmodi.comppare[1])
            print('fluxpare: ', gdatmodi.comppare[2])
            print('auxipara[0][0]: ', gdat.anglfact * gdatmodi.this.auxipara[0])
            print('auxipara[0][1]: ', gdat.anglfact * gdatmodi.this.auxipara[1])
            print('auxipara[0][2]: ', gdatmodi.this.auxipara[2])
                
    if gmod.numbparaelem > 0 and gdatmodi.this.indxproptype > 0 and gdatmodi.this.boolpropfilt:
        # change the number of elements
        if gdatmodi.this.indxproptype == 1 or gdatmodi.this.indxproptype == 3:
            gmodnext.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]] = gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]] + 1
        if gdatmodi.this.indxproptype == 2 or gdatmodi.this.indxproptype == 4:
            gmodnext.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]] = gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]] - 1
        gmodnext.paragenrunitfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]] = gmodnext.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]]
        
        # remove the element from the occupied element list
        if (gdatmodi.this.indxproptype == 2 or gdatmodi.this.indxproptype == 4):
            for a, indxelem in enumerate(gdatmodi.indxelemmodi):
                if a == 0 and gdatmodi.this.indxproptype == 2 or a == 1 and gdatmodi.this.indxproptype == 4:
                    gmodnext.indxelemfull[gdatmodi.indxpopltran].remove(indxelem)
    
    if gdatmodi.this.indxproptype == 0:
        gdatmodi.indxsampmodi = thisindxsampfull
    else:
        if gdatmodi.this.indxproptype == 1:
            gdatmodi.indxsampmodi = np.concatenate((np.array([gmod.indxpara.numbelem[gdatmodi.indxpopltran]]), gdatmodi.indxsamptran[0]))
        if gdatmodi.this.indxproptype == 2:
            gdatmodi.indxsampmodi = [gmod.indxpara.numbelem[gdatmodi.indxpopltran]]
        if gdatmodi.this.indxproptype == 3:
            gdatmodi.indxsampmodi = np.concatenate((np.array([gmod.indxpara.numbelem[gdatmodi.indxpopltran]]), \
                                                                            gdatmodi.indxsamptran[0], gdatmodi.indxsamptran[1]))
        if gdatmodi.this.indxproptype == 4:
            gdatmodi.indxsampmodi = np.concatenate((np.array([gmod.indxpara.numbelem[gdatmodi.indxpopltran]]), gdatmodi.indxsamptran[0]))
    
    if gmod.numbparaelem > 0:
        if gdatmodi.this.indxproptype == 0:
            indxparagenrfullelem = thisindxparagenrfullelem
        else:
            indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmodnext.indxelemfull, strgmodl)
    if gdat.typeverb > 1:
        print('gdatmodi.indxsampmodi')
        print(gdatmodi.indxsampmodi)
        if gmod.numbparaelem > 0:
            print('gmodthis.indxelemfull')
            print(gmodthis.indxelemfull)
            print('gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]].astype(int)')
            print(gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]].astype(int))
            if gdatmodi.this.indxproptype > 0:
                print('gdatmodi.indxelemmodi')
                print(gdatmodi.indxelemmodi)
                print('gdatmodi.indxelemfullmodi')
                print(gdatmodi.indxelemfullmodi)
                print('gdatmodi.this.boolpropfilt')
                print(gdatmodi.this.boolpropfilt)
            print('indxparagenrfullelem')
            print(indxparagenrfullelem)
    
    if gdatmodi.this.indxproptype == 1:
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            gmodnext.paragenrscalfull[gdatmodi.indxsamptran[0][g]] = icdf_trap(gdat, strgmodl, gdatmodi.this.auxipara[g], gmodthis.paragenrscalfull, \
                                                                            gmod.listscalparagenrelem[gdatmodi.indxpopltran][g], \
                                                                            gmod.namepara.genrelem[gdatmodi.indxpopltran][g], gdatmodi.indxpopltran)

    if gdat.booldiagmode:
        if gmod.numbparaelem > 0:
            for l in gmod.indxpopl:
                if gmodthis.paragenrunitfull[gmod.indxpara.numbelem[l]] != round(gmodthis.paragenrunitfull[gmod.indxpara.numbelem[l]]):
                    print('l')
                    print(l)
                    print('gmod.indxpara.numbelem')
                    print(gmod.indxpara.numbelem)
                    print('gmodthis.paragenrunitfull')
                    print(gmodthis.paragenrunitfull)
                    raise Exception('')
                if gmodthis.paragenrscalfull[gmod.indxpara.numbelem[l]] != round(gmodthis.paragenrscalfull[gmod.indxpara.numbelem[l]]):
                    raise Exception('')
                if gmodnext.paragenrunitfull[gmod.indxpara.numbelem[l]] != round(gmodnext.paragenrunitfull[gmod.indxpara.numbelem[l]]):
                    raise Exception('')
                if gmodnext.paragenrscalfull[gmod.indxpara.numbelem[l]] != round(gmodnext.paragenrscalfull[gmod.indxpara.numbelem[l]]):
                    raise Exception('')

        if strgmodl == 'fitt':
            diffparagenrscalfull = abs(gmodnext.paragenrscalfull - gmodthis.paragenrscalfull)
            #size = np.where(((gmodthis.paragenrscalfull == 0.) & (diffparagenrscalfull > 0.)) | ((gmodthis.paragenrscalfull != 0.) & (diffparagenrscalfull / gmodthis.paragenrscalfull > 0)))[0].size
            size = np.where(diffparagenrscalfull != 0.)[0].size
            if gdatmodi.this.indxproptype == 1:
                if size - 1 != gmod.numbparagenrelemsing[gdatmodi.indxpopltran]:
                    raise Exception('')
    

def calc_probprop(gdat, gdatmodi):
    
    gmod = gdat.fitt

    # calculate the factor to multiply the acceptance rate, i.e., 
    ## probability of the auxiliary parameters,
    if gdatmodi.this.indxproptype == 0:
        gdatmodi.this.lpau = 0.
    elif gdatmodi.this.indxproptype == 1 or gdatmodi.this.indxproptype == 2:
        gdatmodi.this.lpau = gdatmodi.next.lpritotl - gdatmodi.this.lpritotl
        lpautemp = 0.5 * gdat.priofactdoff * gmod.numbparagenrelemsing[gdatmodi.indxpopltran]
        if gdatmodi.this.indxproptype == 1:
            gdatmodi.this.lpau += lpautemp
        if gdatmodi.this.indxproptype == 2:
            gdatmodi.this.lpau -= lpautemp
    elif gdatmodi.this.indxproptype == 3 or gdatmodi.this.indxproptype == 4:
        gdatmodi.this.lpau = 0.
        dictelemtemp = [dict()]
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if gmod.gmod.boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.this.lpau += -0.5 * np.log(2. * np.pi * gdat.radispmr**2) - 0.5 * (gdatmodi.this.auxipara[g] / gdat.radispmr)**2
            elif g != gmod.indxparagenrelemampl[gdatmodi.indxpopltran]:
                dictelemtemp[0][nameparagenrelem] = gdatmodi.this.auxipara[g]
                gdatmodi.this.lpau += retr_lprielem(gdat, 'fitt', gdatmodi.indxpopltran, g, \
                                            gmod.namepara.genrelem[gdatmodi.indxpopltran][g], gmod.listscalparagenrelem[gdatmodi.indxpopltran][g], \
                                            gdatmodi.this.paragenrscalfull, dictelemtemp, [1])
        if gdatmodi.this.indxproptype == 4:
            gdatmodi.this.lpau *= -1.

    if gdatmodi.this.indxproptype > 2 and gdatmodi.this.boolpropfilt:
        ## the ratio of the probability of the reverse and forward proposals, and
        if gdatmodi.this.indxproptype == 3:
            gdatmodi.this.probmergtotl = retr_probmerg(gdat, gdatmodi, gdatmodi.next.paragenrscalfull, gdatmodi.next.indxparagenrfullelem, gdatmodi.indxpopltran, 'pair', \
                                                                               typeelem=gmod.typeelem)
            gdatmodi.this.ltrp = np.log(gdatmodi.this.numbelem[gdatmodi.indxpopltran] + 1) + np.log(gdatmodi.this.probmergtotl)

        else:
            gdatmodi.this.probmergtotl = retr_probmerg(gdat, gdatmodi, gdatmodi.this.paragenrscalfull, gdatmodi.this.indxparagenrfullelem, gdatmodi.indxpopltran, 'pair', \
                                                                               typeelem=gmod.typeelem)
            
            gdatmodi.this.ltrp = -np.log(gdatmodi.this.numbelem[gdatmodi.indxpopltran]) - np.log(gdatmodi.this.probmergtotl)
        
        ## Jacobian
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.this.ljcb = np.log(gdatmodi.comppare[1])
        else:
            gdatmodi.this.ljcb = np.log(gdatmodi.comppare[2])
        if gdatmodi.this.indxproptype == 4:
            gdatmodi.this.ljcb *= -1.
        
    else:
        gdatmodi.this.ljcb = 0.
        gdatmodi.this.ltrp = 0.
    
    for l in gmod.indxpopl:
        if gdatmodi.this.indxproptype > 0:
            setattr(gdatmodi, 'auxiparapop%d' % l, gdatmodi.this.auxipara)


def retr_indxparagenrfullelem(gdat, indxelemfull, strgmodl):

    gmod = getattr(gdat, strgmodl)
    
    ## element parameters
    if gmod.numbparaelem > 0:
        indxparagenrfullelem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            indxparagenrfulltemp = gmod.indxparagenrfulleleminit + gmod.numbparagenrelemcuml[l] + np.array(indxelemfull[l], dtype=int) * gmod.numbparagenrelemsing[l]
            
            cntr = tdpy.cntr()
            
            indxparagenrfullelem[l] = dict()
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                indxparagenrfullelem[l][nameparagenrelem] = indxparagenrfulltemp + cntr.incr()
            indxparagenrfullelem[l]['full'] = np.repeat(indxparagenrfulltemp, gmod.numbparagenrelemsing[l]) + np.tile(gmod.indxparagenrelemsing[l], len(indxelemfull[l]))
        
        if gdat.booldiagmode:
            for l in gmod.indxpopl:
                if len(indxparagenrfullelem[l]['full']) > 0:
                    if np.amax(indxparagenrfullelem[l]['full']) > gmod.numbparagenrelem[l] + gmod.numbparagenrbase:
                        print('strgmodl')
                        print(strgmodl)
                        print('strgstat')
                        print(strgstat)
                        print('gmod.numbparagenrbase')
                        print(gmod.numbparagenrbase)
                        print('gmod.numbparagenrelem[l]')
                        print(gmod.numbparagenrelem[l])
                        print('indxparagenrfullelem[l][full]')
                        summgene(indxparagenrfullelem[l]['full'])
                        print('gdat.fitt.minmpara.numbelempop0')
                        print(gdat.fitt.minmpara.numbelempop0)
                        print('gdat.fitt.maxmpara.numbelempop0')
                        print(gdat.fitt.maxmpara.numbelempop0)
                        raise Exception('Element parameter indices are bad.') 
        
    else:
        indxparagenrfullelem = None
    
    return indxparagenrfullelem
    

def retr_weigmergodim(gdat, elin, elinothr):
    
    weigmerg = np.exp(-0.5 * ((elin - elinothr) / gdat.radispmr)**2)
    
    return weigmerg


def retr_weigmergtdim(gdat, lgal, lgalothr, bgal, bgalothr):
    
    weigmerg = np.exp(-0.5 * (((lgal - lgalothr) / gdat.radispmr)**2 + ((bgal - bgalothr) / gdat.radispmr)**2))
    
    return weigmerg


def retr_probmerg(gdat, gdatmodi, paragenrscalfull, indxparagenrfullelem, indxpopltran, strgtype, typeelem=None):
    
    # calculate the weights
    if strgtype == 'seco':
        numb = 1
    if strgtype == 'pair':
        numb = 2
    listweigmerg = []
    for a in range(numb):
        if gmod.typeelem[indxpopltran].startswith('lghtline'):
            elintotl = paragenrscalfull[indxparagenrfullelem['elin'][indxpopltran]]
            elin = elintotl[gdatmodi.indxelemfullmodi[0]]
            elinothr = np.concatenate((elintotl[:gdatmodi.indxelemfullmodi[0]], elintotl[gdatmodi.indxelemfullmodi[0]+1:]))
            weigmerg = retr_weigmergodim(gdat, elin, elinothr)
        else:
            lgaltotl = paragenrscalfull[indxparagenrfullelem['lgal'][indxpopltran]]
            bgaltotl = paragenrscalfull[indxparagenrfullelem['bgal'][indxpopltran]]
            lgal = lgaltotl[gdatmodi.indxelemfullmodi[0]]
            bgal = bgaltotl[gdatmodi.indxelemfullmodi[0]]
            lgalothr = np.concatenate((lgaltotl[:gdatmodi.indxelemfullmodi[0]], lgaltotl[gdatmodi.indxelemfullmodi[0]+1:]))
            bgalothr = np.concatenate((bgaltotl[:gdatmodi.indxelemfullmodi[0]], bgaltotl[gdatmodi.indxelemfullmodi[0]+1:]))
            weigmerg = retr_weigmergtdim(gdat, lgal, lgalothr, bgal, bgalothr)
        listweigmerg.append(weigmerg) 

    # determine the probability of merging the second element given the first element
    if strgtype == 'seco':
        probmerg = listweigmerg[0] / np.sum(listweigmerg[0])
    
    # determine the probability of merging the pair
    if strgtype == 'pair':
        if gmod.typeelem[indxpopltran].startswith('lghtline'):
            weigpair = retr_weigmergtdim(gdat, elin, elintotl[gdatmodi.indxelemfullmodi[1]])
        else:
            weigpair = retr_weigmergtdim(gdat, lgal, lgaltotl[gdatmodi.indxelemfullmodi[1]], bgal, bgaltotl[gdatmodi.indxelemfullmodi[1]])
        probmerg = weigpair / np.sum(listweigmerg[0]) + weigpair / np.sum(listweigmerg[1])
        
    if gdat.booldiagmode:
        if not np.isfinite(probmerg).all():
            raise Exception('Merge probability is infinite.')

    return probmerg

    
def retr_indxparaelem(gmod, l, u):

    indxsamppnts = gmod.indxparagenrfulleleminit + gmod.numbparagenrelemcuml[l] + u * gmod.numbparagenrelemsing[l] + gmod.indxparagenrelemsing[l]

    return indxsamppnts


def gang_detr():

    gang, aang, lgal, bgal = sympy.symbols('gang aang lgal bgal')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])


def retr_psfn(gdat, psfp, indxenertemp, thisangl, typemodlpsfn, strgmodl):

    gmod = getattr(gdat, strgmodl)
    
    indxpsfpinit = gmod.numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    
    if gdat.typeexpr == 'ferm':
        scalangl = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(gdat.binspara.angl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        scalangl = thisangl[None, :, None]
    
    if typemodlpsfn == 'singgaus':
        sigc = psfp[indxpsfpinit]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    
    elif typemodlpsfn == 'singking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(scalangl, sigc, gamc)
    
    elif typemodlpsfn == 'doubking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        gamt = psfp[indxpsfpinit+3]
        frac = psfp[indxpsfpinit+4]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        frac = frac[:, None, :]
        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
        if gdat.typeexpr == 'ferm':
            psfnnorm = retr_doubking(scalanglnorm, frac, sigc, gamc, sigt, gamt)
        
    # normalize the PSF
    if gdat.typeexpr == 'ferm':
        fact = 2. * np.pi * np.trapz(psfnnorm * np.sin(gdat.binspara.angl[None, :, None]), gdat.binspara.angl, axis=1)[:, None, :]
        psfn /= fact

    return psfn


def retr_unit(lgal, bgal):

    xdat = np.cos(bgal) * np.cos(lgal)
    ydat = -np.cos(bgal) * np.sin(lgal)
    zaxi = np.sin(bgal)

    return xdat, ydat, zaxi


def retr_psec(gdat, conv):

    # temp
    conv = conv.reshape((gdat.numbsidecart, gdat.numbsidecart))
    psec = (abs(scipy.fftpack.fft2(conv))**2)[:gdat.numbsidecarthalf, :gdat.numbsidecarthalf] * 1e-3
    psec = psec.flatten()

    return psec
   

def retr_psecodim(gdat, psec):
    
    psec = psec.reshape((gdat.numbsidecarthalf, gdat.numbsidecarthalf))
    psecodim = np.zeros(gdat.numbsidecarthalf)
    for k in gdat.indxmpolodim:
        indxmpol = np.where((gdat.meanpara.mpol > gdat.binspara.mpolodim[k]) & (gdat.meanpara.mpol < gdat.binspara.mpolodim[k+1]))
        psecodim[k] = np.mean(psec[indxmpol])
    psecodim *= gdat.meanpara.mpolodim**2
    
    return psecodim


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / np.sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / np.sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff
    

def retr_condcatl(gdat):
  
    # setup
    ## number of stacked samples
    numbstks = 0
    indxtupl = []
    indxstks = []
    indxstksparagenrscalfull = []
    for n in gdat.indxsamptotl:
        indxstks.append([])
        indxstkssamptemp = []
        for l in gmod.indxpopl:
            indxstks[n].append([])
            for k in range(len(gdat.listpostindxelemfull[n][l])):
                indxstks[n][l].append(numbstks)
                indxstkssamptemp.append(numbstks)
                indxtupl.append([n, l, k])
                numbstks += 1
        indxstkssamp.append(np.array(indxstkssamptemp))
    
    if gdat.typeverb > 1:
        print('indxstks')
        print(indxstks)
        print('indxtupl')
        print(indxtupl)
        print('indxstkssamp')
        print(indxstksparagenrscalfull)
        print('numbstks')
        print(numbstks)

    cntr = 0 
    arrystks = np.zeros((numbstks, gmod.numbparagenrelemtotl))
    for n in gdat.indxsamptotl:
        indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gdat.listpostindxelemfull[n], 'fitt') 
        for l in gmod.indxpopl:
            for k in np.arange(len(gdat.listpostindxelemfull[n][l])):
                for m, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    arrystks[indxstks[n][l][k], m] = gdat.listpostparagenrscalfull[n, gmodstat.indxparagenrfullelem[l][nameparagenrelem][k]]

    if gdat.typeverb > 0:
        print('Constructing the distance matrix for %d stacked samples...' % arrystks.shape[0])
        timeinit = gdat.functime()
    
    gdat.distthrs = np.empty(gmod.numbparagenrelemtotl)
    for k, nameparagenrelem in enumerate(gmod.namepara.elem):
       # temp
       l = 0
       gdat.distthrs[k] = gdat.stdp[getattr(gdat, 'indxstdppop%d' % l + nameparagenrelem)]
    
    # construct lists of samples for each proposal type
    listdisttemp = [[] for k in range(gmod.numbparagenrelemtotl)]
    indxstksrows = [[] for k in range(gmod.numbparagenrelemtotl)]
    indxstkscols = [[] for k in range(gmod.numbparagenrelemtotl)]
    thisperc = 0
    cntr = 0
    for k in gmod.indxparagenrelemtotl:
        for n in range(numbstks):
            dist = np.fabs(arrystks[n, k] - arrystks[:, k])
            indxstks = np.where(dist < gdat.distthrs[k])[0]
            if indxstks.size > 0:
                for j in indxstks:
                    cntr += 1
                    listdisttemp[k].append(dist[j])
                    indxstksrows[k].append(n)
                    indxstkscols[k].append(j)
            
            nextperc = np.floor(100. * float(k * numbstks + n) / numbstks / gmod.numbparagenrelemtotl)
            if nextperc > thisperc:
                thisperc = nextperc
            if cntr > 1e6:
                break
        
        listdisttemp[k] = np.array(listdisttemp[k])
        indxstksrows[k] = np.array(indxstksrows[k])
        indxstkscols[k] = np.array(indxstkscols[k])

        if cntr > 1e6:
            break
    
    listdist = [[] for k in range(gmod.numbparagenrelemtotl)]
    for k, nameparagenrelem in enumerate(gmod.namepara.elem):
        listdist[k] = scipy.sparse.csr_matrix((listdisttemp[k], (indxstksrows[k], indxstkscols[k])), shape=(numbstks, numbstks))
    
    listindxstkspair = []
    indxstksleft = []

    if gdat.typeverb > 0:
        timefinl = gdat.functime()
    
    indxstksleft = range(numbstks)

    # list of sample lists of the labeled element
    indxstksassc = []
    cntr = 0
    
    gdat.prvlthrs = 0.05

    while len(indxstksleft) > 0:
        
        # count number of associations
        numbdist = np.zeros(numbstks, dtype=int) - 1
        for p in range(len(indxstksleft)):
            indxindx = np.where((listdist[0][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmlgal < gdat.anglassc) & \
                             (listdist[1][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmbgal < gdat.anglassc))[0]
            numbdist[indxstksleft[p]] = indxindx.size
            
        prvlmaxmesti = np.amax(numbdist) / float(gdat.numbsamptotl)
        
        if prvlmaxmesti < gdat.prvlthrs:
            break

        # determine the element with the highest number of neighbors
        indxstkscntr = np.argmax(numbdist)
        indxsamptotlcntr = indxtupl[indxstkscntr][0]
        indxpoplcntr = indxtupl[indxstkscntr][1]
        indxelemcntr = indxtupl[indxstkscntr][2]

        # add the central element sample
        indxstksassc.append([])
        indxstksassc[cntr].append(indxstkscntr)
        indxstksleft.remove(indxstkscntr)

        if gdat.typeverb > 1:
            print('Match step %d' % cntr)
            print('numbdist')
            print(numbdist)
            print('indxstkscntr')
            print(indxstkscntr)
            print('indxstksleft')
            print(indxstksleft)
        
        # add the associated element samples
        if len(indxstksleft) > 0:
            for n in gdat.indxsamptotl:
                
                indxstkstemp = np.intersect1d(np.array(indxstksleft), indxstksparagenrscalfull[n])
                
                if n == indxsamptotlcntr:
                    continue
                
                if indxstkstemp.size > 0:
                    totl = np.zeros_like(indxstkstemp)
                    for k in gmod.indxparagenrelemtotl:
                        temp = listdist[k][indxstkscntr, indxstkstemp].tonp.array()[0]
                        totl = totl + temp**2

                    indxleft = np.argsort(totl)[0]
                    
                    indxstksthis = indxstkstemp[indxleft]
                
                    thisbool = True
                    for k in gmod.indxparagenrelemtotl:
                        if listdist[k][indxstkscntr, indxstksthis] > gdat.distthrs[k]:
                            thisbool = False

                    if thisbool:
                        indxstksassc[cntr].append(indxstksthis)
                        indxstksleft.remove(indxstksthis)
            
                # temp
                #if gdat.makeplot:
                #    gdatmodi = tdpy.gdatstrt()
                #    gdatmodi.this.indxelemfull = deepcopy(listindxelemfull[n])
                #    for r in range(len(indxstksassc)): 
                #        calc_poststkscond(gdat, indxstksassc)
                #    gdatmodi.this.indxelemfull = [[] for l in gmod.indxpopl]
                #    for indxstkstemp in indxstksleft:
                #        indxsamptotlcntr = indxtupl[indxstkstemp][0]
                #        indxpoplcntr = indxtupl[indxstkstemp][1]
                #        indxelemcntr = indxtupl[indxstkstemp][2]
                #        gdatmodi.this.paragenrscalfull = gdat.listparagenrscalfull[indxsamptotlcntr, :]
                #        gdatmodi.this.indxelemfull[].append()

                #    plot_genemaps(gdat, gdatmodi, 'this', 'cntpdata', strgpdfn, indxenerplot=0, indxevttplot=0, cond=True)
                
            cntr += 1
        
    gdat.dictglob['poststkscond'] = []
    gdat.dictglob['liststkscond'] = []
    # for each condensed element
    for r in range(len(indxstksassc)): 
        gdat.dictglob['liststkscond'].append([])
        gdat.dictglob['liststkscond'][r] = {}
        gdat.dictglob['poststkscond'].append([])
        gdat.dictglob['poststkscond'][r] = {}
        for strgfeat in gmod.namepara.genrelem:
            gdat.dictglob['liststkscond'][r][strgfeat] = []

        # for each associated sample associated with the central stacked sample 
        for k in range(len(indxstksassc[r])):
            indxsamptotlcntr = indxtupl[indxstksassc[r][k]][0]
            indxpoplcntr = indxtupl[indxstksassc[r][k]][1]
            indxelemcntr = indxtupl[indxstksassc[r][k]][2]
            
            for strgfeat in gmod.namepara.genrelem:
                temp = getattr(gdat, 'list' + strgfeat)
                if temp[indxsamptotlcntr][indxpoplcntr].size > 0:
                    temp = temp[indxsamptotlcntr][indxpoplcntr][..., indxelemcntr]
                    gdat.dictglob['liststkscond'][r][strgfeat].append(temp)

    for r in range(len(gdat.dictglob['liststkscond'])):
        for strgfeat in gmod.namepara.genrelem:
            arry = np.stack(gdat.dictglob['liststkscond'][r][strgfeat], axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat] = np.zeros(([3] + list(arry.shape[1:])))
            gdat.dictglob['poststkscond'][r][strgfeat][0, ...] = median(arry, axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][1, ...] = percennp.tile(arry, 16., axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][2, ...] = percennp.tile(arry, 84., axis=0)
            
    gdat.numbstkscond = len(gdat.dictglob['liststkscond'])

    gdat.indxstkscond = np.arange(gdat.numbstkscond)
    gdat.prvl = np.empty(gdat.numbstkscond)
    for r in gdat.indxstkscond:
        gdat.prvl[r] = len(gdat.dictglob['liststkscond'][r]['deltllik'])
    gdat.prvl /= gdat.numbsamptotl
    gdat.minmprvl = 0.
    gdat.maxmprvl = 1.
    retr_axis(gdat, 'prvl')
    gdat.histprvl = np.histogram(gdat.prvl, bins=gdat.binspara.prvl)[0]
    if gdat.makeplot:
        pathcond = getattr(gdat, 'path' + strgpdfn + 'finlcond')
        for k, nameparagenrelem in enumerate(gmod.namepara.elem):
            path = pathcond + 'histdist' + nameparagenrelem 
            listtemp = np.copy(listdist[k].tonp.array()).flatten()
            listtemp = listtemp[np.where(listtemp != 1e20)[0]]
            tdpy.mcmc.plot_hist(path, listtemp, r'$\Delta \tilde{' + getattr(gmod.lablrootpara, nameparagenrelem) + '}$')
            path = pathcond + 'histprvl'
            tdpy.mcmc.plot_hist(path, gdat.prvl, r'$p$')
    gdat.prvlthrs = 0.1 
    gdat.indxprvlhigh = np.where(gdat.prvl > gdat.prvlthrs)[0]
    gdat.numbprvlhigh = gdat.indxprvlhigh.size


def retr_conv(gdat, defl):
    
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    # temp
    conv = abs(np.gradient(defl[:, :, 0], gdat.sizepixl, axis=0) + np.gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) / 2.
    conv = conv.flatten()
    
    return conv


def retr_invm(gdat, defl):
    
    # temp
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    invm = (1. - np.gradient(defl[:, :, 0], gdat.sizepixl, axis=0)) * (1. - np.gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) - \
                                                np.gradient(defl[:, :, 0], gdat.sizepixl, axis=1) * np.gradient(defl[:, :, 1], gdat.sizepixl, axis=0)
    invm = invm.flatten()
    return invm


def setp_indxswepsave(gdat):

    gdat.indxswep = np.arange(gdat.numbswep)
    gdat.boolsave = np.zeros(gdat.numbswep, dtype=bool)
    gdat.indxswepsave = np.arange(gdat.numbburn, gdat.numbburn + gdat.numbsamp * gdat.factthin, gdat.factthin)
    gdat.boolsave[gdat.indxswepsave] = True
    gdat.indxsampsave = np.zeros(gdat.numbswep, dtype=int) - 1
    gdat.indxsampsave[gdat.indxswepsave] = np.arange(gdat.numbsamp)
    

def retr_cntspnts(gdat, listposi, spec):
    
    cnts = np.zeros((gdat.numbener, spec.shape[1]))
    
    if gdat.boolbinsspat:
        lgal = listposi[0]
        bgal = listposi[1]
        indxpixlpnts = retr_indxpixl(gdat, bgal, lgal)
    else:
        elin = listposi[0]
        indxpixlpnts = np.zeros_like(elin, dtype=int)
    for k in range(spec.shape[1]):
        cnts[:, k] += spec[:, k] * gdat.expototl[:, indxpixlpnts[k]]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None]
    cnts = np.sum(cnts, axis=0)

    return cnts


def retr_mdencrit(gdat, adissour, adishost, adishostsour):
    
    mdencrit = gdat.factnewtlght / 4. / np.pi * adissour / adishostsour / adishost
        
    return mdencrit


def retr_massfrombein(gdat, adissour, adishost, adishostsour):

    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    massfrombein = np.pi * adishost**2 * mdencrit

    return massfrombein


def retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut):
    
    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    
    fracacutasca = acut / asca
    
    factmcutfromdefs = np.pi * adishost**2 * mdencrit * asca * retr_mcutfrommscl(fracacutasca)

    return factmcutfromdefs


def retr_mcut(gdat, defs, asca, acut, adishost, mdencrit):
    
    mscl = defs * np.pi * adishost**2 * mdencrit * asca
    fracacutasca = acut / asca
    mcut = mscl * retr_mcutfrommscl(fracacutasca)
    
    return mcut


def retr_mcutfrommscl(fracacutasca):
    
    mcut = fracacutasca**2 / (fracacutasca**2 + 1.)**2 * ((fracacutasca**2 - 1.) * np.log(fracacutasca) + fracacutasca * np.pi - (fracacutasca**2 + 1.))

    return mcut


def retr_negalogt(varb):
    
    negalogt = sign(varb) * np.log10(np.fabs(varb))
    
    return negalogt


def retr_gradmaps(gdat, maps):
    
    # temp -- this does not work with vanishing exposure
    maps = maps.reshape((gdat.numbsidecart, gdat.numbsidecart))
    grad = np.dstack((np.gradient(maps, gdat.sizepixl, axis=0), np.gradient(maps, gdat.sizepixl, axis=1))).reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    grad = grad.reshape((gdat.numbpixlcart, 2))

    return grad


def retr_spatmean(gdat, inpt, boolcntp=False):
    
    listspatmean = [[] for b in gdat.indxspatmean]
    listspatstdv = [[] for b in gdat.indxspatmean]
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        if boolcntp:
            cntp = inpt[gdat.listindxcubespatmean[b]]
        else:
            cntp = inpt[gdat.listindxcubespatmean[b]] * gdat.expo[gdat.listindxcubespatmean[b]] * gdat.apix
            if gdat.enerdiff:
                cntp *= gdat.deltener[:, None, None]
        spatmean = np.mean(np.sum(cntp, 2), axis=1) / gdat.apix
        spatstdv = np.sqrt(np.sum(cntp, axis=(1, 2))) / gdat.numbdata / gdat.apix
        if gdat.boolcorrexpo:
            spatmean /= gdat.expototlmean
            spatstdv /= gdat.expototlmean
        if gdat.enerdiff:
            spatmean /= gdat.deltener
            spatstdv /= gdat.deltener
        listspatmean[b] = spatmean
        listspatstdv[b] = spatstdv

    return listspatmean, listspatstdv


def retr_rele(gdat, maps, lgal, bgal, defs, asca, acut, indxpixlelem, absv=True, cntpmodl=None):
    
    grad = retr_gradmaps(gdat, maps)
        
    defl = retr_defl(gdat, indxpixlelem, lgal, bgal, defs, asca=asca, acut=acut)

    prod = grad * defl
    if cntpmodl is not None:
        prod /= cntpmodl[:, None]
    dotstemp = np.sum(prod, 1)
    if absv:
        dotstemp = np.fabs(dotstemp)
    else:
        dotstemp = dotstemp
    
    dots = np.mean(dotstemp)
    
    return dots


def retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn, strgmome='pmea', indxvarb=None, indxlist=None):
    
    if strgvarb.startswith('cntpdata'):
        varb = getattr(gdat, strgvarb)
    elif strgvarb.startswith('histcntpdata'):
        varb = getattr(gdat, strgvarb)
    else:
        if strgmodl == 'true':
            gmod = getattr(gdat, strgmodl)
            gmodstat = getattr(gmod, strgstat)
            varb = getattr(gmodstat, strgvarb)
        if strgmodl == 'fitt':
            if strgstat == 'this':
                if strgmome == 'errr':
                    varb = getattr(gdatmodi, strgstat + 'errr' + strgvarb)
                else:
                    varb = getattr(gdatmodi, strgstat + strgvarb)
            if strgstat == 'pdfn':
                varb = getattr(gdat, strgmome + strgpdfn + strgvarb)

    if indxlist is not None:
        varb = varb[indxlist]

    if indxvarb is not None:
        if strgmome == 'errr':
            varb = varb[[slice(None)] + indxvarb]
        else:
            varb = varb[indxvarb]

    return np.copy(varb)


def setp_indxpara(gdat, typesetp, strgmodl='fitt'):
    
    print('setp_indxpara(): Building parameter indices for model %s with type %s...' % (strgmodl, typesetp))

    gmod = getattr(gdat, strgmodl)
    
    if typesetp == 'init':
        
        if strgmodl == 'fitt':
            gmod.lablmodl = 'Model'
        if strgmodl == 'true':
            gmod.lablmodl = 'True'

        # transdimensional element populations
        
        gmod.numbpopl = len(gmod.typeelem)
        gmod.indxpopl = np.arange(gmod.numbpopl)
        
        if gdat.typeexpr != 'user':
            # background component
            gmod.numbback = 0
            gmod.indxback = []
            for c in range(len(gmod.typeback)):
                if isinstance(gmod.typeback[c], str):
                    if gmod.typeback[c].startswith('bfunfour') or gmod.typeback[c].startswith('bfunwfou'):
                        namebfun = gmod.typeback[c][:8]
                        ordrexpa = int(gmod.typeback[c][8:])
                        numbexpa = 4 * ordrexpa**2
                        indxexpa = np.arange(numbexpa)
                        del gmod.typeback[c]
                        for k in indxexpa:
                            gmod.typeback.insert(c+k, namebfun + '%04d' % k)
            gmod.numbback = len(gmod.typeback)
            gmod.indxback = np.arange(gmod.numbback)
            gmod.numbbacktotl = np.sum(gmod.numbback)
            gmod.indxbacktotl = np.arange(gmod.numbbacktotl)
            
            # galaxy components
            gmod.indxsersfgrd = np.arange(gmod.numbsersfgrd)

            # name of the generative element parameter used for the amplitude
            gmod.nameparagenrelemampl = [[] for l in gmod.indxpopl]
            gmod.indxparagenrelemampl = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lghtpntspuls':
                    gmod.nameparagenrelemampl[l] = 'per0'
                    gmod.indxparagenrelemampl[l] = 2
                elif gmod.typeelem[l] == 'lghtpntsagnntrue':
                    gmod.nameparagenrelemampl[l] = 'lum0'
                    gmod.indxparagenrelemampl[l] = 2
                elif gmod.typeelem[l].startswith('lghtline'):
                    gmod.nameparagenrelemampl[l] = 'flux'
                    gmod.indxparagenrelemampl[l] = 1
                elif gmod.typeelem[l].startswith('lghtpnts'):
                    gmod.nameparagenrelemampl[l] = 'flux'
                    gmod.indxparagenrelemampl[l] = 2
                elif gmod.typeelem[l].startswith('lghtgausbgrd'):
                    gmod.nameparagenrelemampl[l] = 'flux'
                    gmod.indxparagenrelemampl[l] = 2
                if gmod.typeelem[l] == 'lens':
                    gmod.nameparagenrelemampl[l] = 'defs'
                    gmod.indxparagenrelemampl[l] = 2
                if gmod.typeelem[l].startswith('clus'):
                    gmod.nameparagenrelemampl[l] = 'nobj'
                    gmod.indxparagenrelemampl[l] = 2
                if gmod.typeelem[l] == 'lens':
                    gmod.nameparagenrelemampl[l] = 'defs'
                if gmod.typeelem[l] == 'clus':
                    gmod.nameparagenrelemampl[l] = 'nobj'
                if len(gmod.nameparagenrelemampl[l]) == 0:
                    raise Exception('Amplitude feature undefined.')
    
        for featpara in gdat.listfeatpara:
            for strggrop in gdat.liststrggroppara:
                setattr(gmod, 'list' + featpara + 'para' + strggrop, [])
    
    if typesetp == 'finl':
        
        # number of elements in the current state of the true model
        if strgmodl == 'true':
            gmod.numbelem = np.zeros(gmod.numbpopl)
            for l in gmod.indxpopl:
                gmod.numbelem[l] += getattr(gmod.maxmpara, 'numbelempop%d' % l)
            gmod.numbelemtotl = np.sum(gmod.numbelem) 
         
        # element setup
        ## flag to calculate the kernel approximation errors
        boolcalcerrr = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelemspateval[l] == 'locl' and gdat.numbpixlfull < 1e5:
                # temp
                boolcalcerrr[l] = False
            else:
                boolcalcerrr[l] = False
        setp_varb(gdat, 'boolcalcerrr', valu=boolcalcerrr, strgmodl=strgmodl)
        
        # maximum number of elements for each population
        gmod.maxmpara.numbelem = np.zeros(gmod.numbpopl, dtype=int)
        for l in gmod.indxpopl:
            gmod.maxmpara.numbelem[l] = getattr(gmod.maxmpara, 'numbelempop%d' % l)
        
        # maximum number of elements summed over all populations
        gmod.maxmpara.numbelemtotl = np.sum(gmod.maxmpara.numbelem) 

        ## sorting feature
        nameparaelemsort = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            # feature to be used to sort elements
            if gmod.typeelem[l].startswith('lght'):
                nameparaelemsort[l] = 'flux'
            if gmod.typeelem[l] == 'lens':
                nameparaelemsort[l] = 'defs'
            if gmod.typeelem[l].startswith('clus'):
                nameparaelemsort[l] = 'nobj'
    
        ## label extensions
        gmod.lablelemextn = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gdat.numbgrid > 1:
                if gmod.typeelem[l] == 'lghtpnts':
                    gmod.lablelemextn[l] = r'\rm{fps}'
                if gmod.typeelem[l] == 'lghtgausbgrd':
                    gmod.lablelemextn[l] = r'\rm{bgs}'
            else:
                if gmod.typeelem[l].startswith('lghtpntspuls'):
                    gmod.lablelemextn[l] = r'\rm{pul}'
                if gmod.typeelem[l].startswith('lghtpntsagnn'):
                    gmod.lablelemextn[l] = r'\rm{agn}'
                elif gmod.typeelem[l] == 'lghtpnts':
                    gmod.lablelemextn[l] = r'\rm{pts}'
            if gmod.typeelem[l] == 'lens':
                gmod.lablelemextn[l] = r'\rm{sub}'
            if gmod.typeelem[l].startswith('clus'):
                gmod.lablelemextn[l] = r'\rm{cls}'
            if gmod.typeelem[l].startswith('lghtline'):
                gmod.lablelemextn[l] = r'\rm{lin}'
    
        gmod.indxpoplgrid = [[] for y in gdat.indxgrid]
        for y in gdat.indxgrid: 
            for indx, typeelemtemp in enumerate(gmod.typeelem):
                # foreground grid (image plane) -- the one np.where the data is measured
                if y == 0:
                    if typeelemtemp.startswith('lght') and not typeelemtemp.endswith('bgrd') or typeelemtemp.startswith('clus'):
                        gmod.indxpoplgrid[y].append(indx)
                # foreground mass grid
                if y == 1:
                    if typeelemtemp.startswith('lens'):
                        gmod.indxpoplgrid[y].append(indx)
                # background grid (source plane)
                if y == 2:
                    if typeelemtemp.endswith('bgrd'):
                        gmod.indxpoplgrid[y].append(indx)
        
        indxgridpopl = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            for y in gdat.indxgrid:
                if l in gmod.indxpoplgrid[y]:
                    indxgridpopl[l] = y
    
        calcelemsbrt = False
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lghtpnts'):
                calcelemsbrt = True
    
        if 'lghtgausbgrd' in gmod.typeelem:
            calcelemsbrtbgrd = True
        else:
            calcelemsbrtbgrd = False

        if gmod.boollenssubh:
            calcelemdefl = True
        else:
            calcelemdefl = False

        ## element Boolean flags
        gmod.boolelemlght = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght'):
                gmod.boolelemlght[l] = True
            else:
                gmod.boolelemlght[l] = False
        gmod.boolelemlghtanyy = True in gmod.boolelemlght
        
        gmod.boolelemlens = False
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lens'):
                gmod.boolelemlens = True
        
        gmod.boolelemsbrtdfnc = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.maxmpara.numbelem[l] > 0 and (gmod.typeelem[l].startswith('lght') and not gmod.typeelem[l].endswith('bgrd') or gmod.typeelem[l].startswith('clus')):
                gmod.boolelemsbrtdfnc[l] = True
            else:
                gmod.boolelemsbrtdfnc[l] = False
        gmod.boolelemsbrtdfncanyy = True in gmod.boolelemsbrtdfnc

        gmod.boolelemdeflsubh = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lens':
                gmod.boolelemdeflsubh[l] = True
            else:
                gmod.boolelemdeflsubh[l] = False
        gmod.boolelemdeflsubhanyy = True in gmod.boolelemdeflsubh

        gmod.boolelemsbrtextsbgrd = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght') and gmod.typeelem[l].endswith('bgrd'):
                gmod.boolelemsbrtextsbgrd[l] = True
            else:
                gmod.boolelemsbrtextsbgrd[l] = False
        gmod.boolelemsbrtextsbgrdanyy = True in gmod.boolelemsbrtextsbgrd
        
        if gmod.boolelemsbrtextsbgrdanyy:
            gmod.indxpopllens = 1
        else:
            gmod.indxpopllens = 0

        gmod.boolelemsbrtpnts = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght') and gmod.typeelem[l] != 'lghtline' or gmod.typeelem[l] == 'clus':
                gmod.boolelemsbrtpnts[l] = True
            else:
                gmod.boolelemsbrtpnts[l] = False
        gmod.boolelemsbrtpntsanyy = True in gmod.boolelemsbrtpnts

        # temp -- because there is currently no extended source
        gmod.boolelemsbrt = gmod.boolelemsbrtdfnc
    
        gmod.boolelempsfn = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lghtpnts') or gmod.typeelem[l] == 'clus':
                gmod.boolelempsfn[l] = True
            else:
                gmod.boolelempsfn[l] = False
        gmod.boolelempsfnanyy = True in gmod.boolelempsfn
        
        spectype = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.boolelemlght[l]:
                spectype[l] = 'powr'
            else:
                spectype[l] = 'none'
        setp_varb(gdat, 'spectype', valu=spectype, strgmodl=strgmodl)
    
        minmgwdt = 2. * gdat.sizepixl
        maxmgwdt = gdat.maxmgangdata / 4.
        setp_varb(gdat, 'gwdt', minm=minmgwdt, maxm=maxmgwdt, strgmodl=strgmodl)
        setp_varb(gdat, 'aerr', minm=-100, maxm=100, strgmodl=strgmodl, popl='full')
    
        if gmod.boolelemlghtanyy:
            # flux
            if gdat.typeexpr == 'ferm':
                minmflux = 1e-9
                maxmflux = 1e-6
            if gdat.typeexpr == 'tess':
                minmflux = 1.
                maxmflux = 1e3
            if gdat.typeexpr == 'chan':
                if gdat.anlytype == 'spec':
                    minmflux = 1e4
                    maxmflux = 1e7
                else:
                    minmflux = 3e-9
                    maxmflux = 1e-6
            if gdat.typeexpr == 'gene':
                minmflux = 0.1
                maxmflux = 100.
            if gdat.typeexpr == 'hubb':
                minmflux = 1e-20
                maxmflux = 1e-17
            if gdat.typeexpr == 'fire':
                minmflux = 1e-20
                maxmflux = 1e-17
            setp_varb(gdat, 'flux', limt=[minmflux, maxmflux], strgmodl=strgmodl)
            
            if gdat.typeexpr == 'ferm':
                setp_varb(gdat, 'brekprioflux', limt=[3e-9, 1e-6], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'sloplowrprioflux', limt=[0.5, 3.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'slopupprprioflux', limt=[0.5, 3.], popl=l, strgmodl=strgmodl)
            
            if gdat.boolbinsener:
                ### spectral parameters
                if gdat.typeexpr == 'ferm':
                    sind = [1., 3.]
                    minmsind = 1.
                    maxmsind = 3.
                if gdat.typeexpr == 'chan':
                    minmsind = 0.4
                    maxmsind = 2.4
                    sind = [0.4, 2.4]
                if gdat.typeexpr == 'hubb':
                    minmsind = 0.5
                    maxmsind = 2.5
                    sind = [0.4, 2.4]
                if gdat.typeexpr != 'fire':
                    setp_varb(gdat, 'sind', limt=[minmsind, maxmsind], strgmodl=strgmodl)
                    setp_varb(gdat, 'curv', limt=[-1., 1.], strgmodl=strgmodl)
                    setp_varb(gdat, 'expc', limt=[0.1, 10.], strgmodl=strgmodl)
                    setp_varb(gdat, 'sinddistmean', limt=sind, popl='full', strgmodl=strgmodl)
                    #### standard deviations should not be too small
                    setp_varb(gdat, 'sinddiststdv', limt=[0.3, 2.], popl='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'curvdistmean', limt=[-1., 1.], popl='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'curvdiststdv', limt=[0.1, 1.], popl='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'expcdistmean', limt=[1., 8.], popl='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'expcdiststdv', limt=[0.01 * gdat.maxmener, gdat.maxmener], popl='full', strgmodl=strgmodl)
                    for i in gdat.indxenerinde:
                        setp_varb(gdat, 'sindcolr0001', limt=[-2., 6.], strgmodl=strgmodl)
                        setp_varb(gdat, 'sindcolr0002', limt=[0., 8.], strgmodl=strgmodl)
                        setp_varb(gdat, 'sindcolr%04d' % i, limt=[-5., 10.], strgmodl=strgmodl)
        
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lghtpntspuls':
                setp_varb(gdat, 'gang', limt=[1e-1 * gdat.sizepixl, gdat.maxmgangdata], strgmodl=strgmodl)
                setp_varb(gdat, 'geff', limt=[0., 0.4], strgmodl=strgmodl)
                setp_varb(gdat, 'dglc', limt=[10., 3e3], strgmodl=strgmodl)
                setp_varb(gdat, 'phii', limt=[0., 2. * np.pi], strgmodl=strgmodl)
                setp_varb(gdat, 'thet', limt=[0., np.pi], strgmodl=strgmodl)
                setp_varb(gdat, 'per0distmean', limt=[5e-4, 1e1], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'magfdistmean', limt=[1e7, 1e16], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'per0diststdv', limt=[1e-2, 1.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'magfdiststdv', limt=[1e-2, 1.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'gangslop', limt=[0.5, 4.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'dglcslop', limt=[0.5, 2.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'spatdistcons', limt=[1e-4, 1e-2], popl='full')
                setp_varb(gdat, 'bgaldistscal', limt=[0.5 / gdat.anglfact, 5. / gdat.anglfact], popl='full', strgmodl=strgmodl)
            if gmod.typeelem[l] == 'lghtpntsagnntrue':
                setp_varb(gdat, 'dlos', limt=[1e7, 1e9], strgmodl=strgmodl)
                setp_varb(gdat, 'dlosslop', limt=[-0.5, -3.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'lum0', limt=[1e43, 1e46], strgmodl=strgmodl)
                setp_varb(gdat, 'lum0distbrek', limt=[1e42, 1e46], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'lum0sloplowr', limt=[0.5, 3.], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'lum0slopuppr', limt=[0.5, 3.], popl=l, strgmodl=strgmodl)
        
        # construct background surface brightness templates from the user input
        gmod.sbrtbacknorm = [[] for c in gmod.indxback]
        gmod.boolunifback = np.ones(gmod.numbback, dtype=bool)
        for c in gmod.indxback:
            gmod.sbrtbacknorm[c] = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            if gmod.typeback[c] == 'data':
                gmod.sbrtbacknorm[c] = np.copy(gdat.sbrtdata)
                gmod.sbrtbacknorm[c][np.where(gmod.sbrtbacknorm[c] == 0.)] = 1e-100
            elif isinstance(gmod.typeback[c], float):
                gmod.sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + gmod.typeback[c]
            elif isinstance(gmod.typeback[c], list) and isinstance(gmod.typeback[c], float):
                gmod.sbrtbacknorm[c] = retr_spec(gdat, np.array([gmod.typeback[c]]), sind=np.array([gmod.typeback[c]]))[:, 0, None, None]
            elif isinstance(gmod.typeback[c], np.ndarray) and gmod.typeback[c].ndim == 1:
                gmod.sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + gmod.typeback[c][:, None, None]
            elif gmod.typeback[c].startswith('bfunfour') or gmod.typeback[c].startswith('bfunwfou'):
                indxexpatemp = int(gmod.typeback[c][8:]) 
                indxterm = indxexpatemp // ordrexpa**2
                indxexpaxdat = (indxexpatemp % ordrexpa**2) // ordrexpa + 1
                indxexpaydat = (indxexpatemp % ordrexpa**2) % ordrexpa + 1
                if namebfun == 'bfunfour':
                    ampl = 1.
                    func = gdat.meanpara.bgalcart 
                if namebfun == 'bfunwfou':
                    functemp = np.exp(-0.5 * (gdat.meanpara.bgalcart / (1. / gdat.anglfact))**2)
                    ampl = np.sqrt(functemp)
                    func = functemp
                argslgal = 2. * np.pi * indxexpaxdat * gdat.meanpara.lgalcart / gdat.maxmgangdata
                argsbgal = 2. * np.pi * indxexpaydat * func / gdat.maxmgangdata
                if indxterm == 0:
                    termfrst = np.sin(argslgal)
                    termseco = ampl * np.sin(argsbgal)
                if indxterm == 1:
                    termfrst = np.sin(argslgal)
                    termseco = ampl * np.cos(argsbgal)
                if indxterm == 2:
                    termfrst = np.cos(argslgal)
                    termseco = ampl * np.sin(argsbgal)
                if indxterm == 3:
                    termfrst = np.cos(argslgal)
                    termseco = ampl * np.cos(argsbgal)
                gmod.sbrtbacknorm[c] = (termfrst[None, :] * termseco[:, None]).flatten()[None, :, None] * \
                                                            np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                
            else:
                path = gdat.pathinpt + gmod.typeback[c]
                gmod.sbrtbacknorm[c] = astropy.io.fits.getdata(path)
                
                if gdat.typepixl == 'cart':
                    if not gdat.boolforccart:
                        if gmod.sbrtbacknorm[c].shape[2] != gdat.numbsidecart:
                            raise Exception('Provided background template must have the chosen image dimensions.')
                    
                    gmod.sbrtbacknorm[c] = gmod.sbrtbacknorm[c].reshape((gmod.sbrtbacknorm[c].shape[0], -1, gmod.sbrtbacknorm[c].shape[-1]))
        
                if gdat.typepixl == 'cart' and gdat.boolforccart:
                    sbrtbacknormtemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                    for i in gdat.indxenerfull:
                        for m in gdat.indxevttfull:
                            sbrtbacknormtemp[i, :, m] = tdpy.retr_cart(gmod.sbrtbacknorm[c][i, :, m], \
                                                    numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).flatten()
                    gmod.sbrtbacknorm[c] = sbrtbacknormtemp

            # determine spatially uniform background templates
            for i in gdat.indxenerfull:
                for m in gdat.indxevttfull:
                    if np.std(gmod.sbrtbacknorm[c][i, :, m]) > 1e-6:
                        gmod.boolunifback[c] = False

        boolzero = True
        gmod.boolbfun = False
        for c in gmod.indxback:
            if np.amin(gmod.sbrtbacknorm[c]) < 0. and isinstance(gmod.typeback[c], str) and not gmod.typeback[c].startswith('bfun'):
                booltemp = False
                raise Exception('Background templates must be positive-definite every where.')
        
            if not np.isfinite(gmod.sbrtbacknorm[c]).all():
                raise Exception('Background template is not finite.')

            if np.amin(gmod.sbrtbacknorm[c]) > 0. or gmod.typeback[c] == 'data':
                boolzero = False
            
            if isinstance(gmod.typeback[c], str) and gmod.typeback[c].startswith('bfun'):
                gmod.boolbfun = True
        
        if boolzero and not gmod.boolbfun:
            raise Exception('At least one background template must be positive everynp.where.')
        
        # temp -- does not take into account dark hosts
        gmod.boolhost = gmod.typeemishost != 'none'
    
        # type of PSF evaluation
        if gmod.maxmpara.numbelemtotl > 0 and gmod.boolelempsfnanyy:
            if gmod.typeemishost != 'none' or not gmod.boolunifback.all():
                # the background is not convolved by a kernel and point sources exist
                typeevalpsfn = 'full'
            else:
                # the background is not convolved by a kernel and point sources exist
                typeevalpsfn = 'kern'
        else:
            if gmod.typeemishost != 'none' or not gmod.boolunifback.all():
                # the background is convolved by a kernel, no point source exists
                typeevalpsfn = 'conv'
            else:
                # the background is not convolved by a kernel, no point source exists
                typeevalpsfn = 'none'
        setp_varb(gdat, 'typeevalpsfn', valu=typeevalpsfn, strgmodl=strgmodl)
        
        if gdat.typeverb > 1:
            print('gmod.typeevalpsfn')
            print(gmod.typeevalpsfn)
        
        gmod.boolapplpsfn = gmod.typeevalpsfn != 'none'

        ### PSF model
        if gmod.typeevalpsfn != 'none':
            
            if gmod.typemodlpsfn == 'singgaus':
                numbpsfpform = 1
            elif gmod.typemodlpsfn == 'singking':
                numbpsfpform = 2
            elif gmod.typemodlpsfn == 'doubgaus':
                numbpsfpform = 3
            elif gmod.typemodlpsfn == 'gausking':
                numbpsfpform = 4
            elif gmod.typemodlpsfn == 'doubking':
                numbpsfpform = 5
            
            gmod.numbpsfptotl = numbpsfpform
            
            if gdat.boolpriopsfninfo:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        meansigc = gmod.psfpexpr[i * gmod.numbpsfptotl + m * gmod.numbpsfptotl * gdat.numbener]
                        stdvsigc = meansigc * 0.1
                        setp_varb(gdat, 'sigcen%02devt%d' % (i, m), mean=meansigc, stdv=stdvsigc, lablroot='$\sigma$', scal='gaus', \
                                                                                                            strgmodl=strgmodl)
                        
                        if gmod.typemodlpsfn == 'doubking' or gmod.typemodlpsfn == 'singking':
                            meangamc = gmod.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvgamc = meangamc * 0.1
                            setp_varb(gdat, 'gamcen%02devt%d' % (i, m), mean=meangamc, stdv=stdvgamc, strgmodl=strgmodl)
                            if gmod.typemodlpsfn == 'doubking':
                                meansigt = gmod.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                                stdvsigt = meansigt * 0.1
                                setp_varb(gdat, 'sigten%02devt%d' % (i, m), mean=meansigt, stdv=stdvsigt, strgmodl=strgmodl)
                                meangamt = gmod.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 3]
                                stdvgamt = meangamt * 0.1
                                setp_varb(gdat, 'gamten%02devt%d' % (i, m), mean=meangamt, stdv=stdvgamt, strgmodl=strgmodl)
                                meanpsff = gmod.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 4]
                                stdvpsff = meanpsff * 0.1
                                setp_varb(gdat, 'psffen%02devt%d' % (i, m), mean=meanpsff, stdv=stdvpsff, strgmodl=strgmodl)
            else:
                if gdat.typeexpr == 'gene':
                    minmsigm = 0.01 / gdat.anglfact
                    maxmsigm = 0.1 / gdat.anglfact
                if gdat.typeexpr == 'ferm':
                    minmsigm = 0.1
                    maxmsigm = 10.
                if gdat.typeexpr == 'hubb':
                    minmsigm = 0.01 / gdat.anglfact
                    maxmsigm = 0.1 / gdat.anglfact
                if gdat.typeexpr == 'chan':
                    minmsigm = 0.1 / gdat.anglfact
                    maxmsigm = 2. / gdat.anglfact
                minmgamm = 1.5
                maxmgamm = 20.
                setp_varb(gdat, 'sigc', minm=minmsigm, maxm=maxmsigm, lablroot='$\sigma_c$', ener='full', evtt='full', strgmodl=strgmodl)

                setp_varb(gdat, 'sigt', minm=minmsigm, maxm=maxmsigm, ener='full', evtt='full', strgmodl=strgmodl)
                setp_varb(gdat, 'gamc', minm=minmgamm, maxm=maxmgamm, ener='full', evtt='full', strgmodl=strgmodl)
                setp_varb(gdat, 'gamt', minm=minmgamm, maxm=maxmgamm, ener='full', evtt='full', strgmodl=strgmodl)
                
            setp_varb(gdat, 'psff', minm=0., maxm=1., ener='full', evtt='full', strgmodl=strgmodl)
 
        # background
        ## number of background parameters
        numbbacp = 0
        for c in gmod.indxback:
            if gmod.boolspecback[c]:
                numbbacp += 1
            else:
                numbbacp += gdat.numbener
   
        ## background parameter indices
        gmod.indxbackbacp = np.zeros(numbbacp, dtype=int)
        indxenerbacp = np.zeros(numbbacp, dtype=int)
        cntr = 0
        for c in gmod.indxback:
            if gmod.boolspecback[c]:
                gmod.indxbackbacp[cntr] = c
                cntr += 1
            else:
                for i in gdat.indxener:
                    indxenerbacp[cntr] = i
                    gmod.indxbackbacp[cntr] = c
                    cntr += 1
        
        # indices of background parameters for each background component
        gmod.indxbacpback = [[] for c in gmod.indxback]
        for c in gmod.indxback:
            gmod.indxbacpback[c] = np.where((gmod.indxbackbacp == c))[0]
                
        # list of names of diffuse components
        gmod.listnamediff = []
        for c in gmod.indxback:
            gmod.listnamediff += ['back%04d' % c]
        if gmod.typeemishost != 'none':
            for e in gmod.indxsersfgrd:
                gmod.listnamediff += ['hostisf%d' % e]
        if gmod.boollens:
            gmod.listnamediff += ['lens']
        
        # list of names of emission components
        listnameecom = deepcopy(gmod.listnamediff)
        for l in gmod.indxpopl:
            if gmod.boolelemsbrt[l]:
                if strgmodl == 'true' and gmod.numbelem[l] > 0 or strgmodl == 'fitt' and gmod.maxmpara.numbelem[l] > 0:
                    if not 'dfnc' in listnameecom:
                        listnameecom += ['dfnc']
                    if not 'dfncsubt' in listnameecom:
                        listnameecom += ['dfncsubt']
        gmod.listnameecomtotl = listnameecom + ['modl']
        
        for c in gmod.indxback:
            setp_varb(gdat, 'cntpback%04d' % c, lablroot='$C_{%d}$' % c, minm=1., maxm=100., scal='logt', strgmodl=strgmodl)
        
        gmod.listnamegcom = deepcopy(gmod.listnameecomtotl)
        if gmod.boollens:
            gmod.listnamegcom += ['bgrd']
            if gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
                gmod.listnamegcom += ['bgrdgalx', 'bgrdexts']
        
        numbdiff = len(gmod.listnamediff)
        convdiff = np.zeros(numbdiff, dtype=bool)
        for k, namediff in enumerate(gmod.listnamediff):
            if not (gdat.boolthindata or gmod.typeevalpsfn == 'none' or gmod.typeevalpsfn == 'kern'):
                if namediff.startswith('back'):
                    indx = int(namediff[-4:])
                    convdiff[k] = not gmod.boolunifback[indx] 
                else:
                    convdiff[k] = True
        
        # element parameters that correlate with the statistical significance of the element
        gmod.namepara.elemsign = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght'):
                gmod.namepara.elemsign[l] = 'flux'
            if gmod.typeelem[l] == 'lens':
                gmod.namepara.elemsign[l] = 'defs'
            if gmod.typeelem[l].startswith('clus'):
                gmod.namepara.elemsign[l] = 'nobj'
    
        if gdat.typeverb > 0:
            if strgmodl == 'true':
                strgtemp = 'true'
            if strgmodl == 'fitt':
                strgtemp = 'fitting'
            print('Building elements for the %s model...' % strgtemp)
        
        # define the names and scalings of element parameters
        gmod.namepara.genrelem = [[] for l in gmod.indxpopl]
        gmod.listscalparagenrelem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            
            if gmod.typeelem[l].startswith('lghtline'):
                gmod.namepara.genrelem[l] = ['elin']
                gmod.listscalparagenrelem[l] = ['logt']
            elif gmod.typespatdist[l] == 'diskscal':
                gmod.namepara.genrelem[l] = ['lgal', 'bgal']
                gmod.listscalparagenrelem[l] = ['self', 'dexp']
            elif gmod.typespatdist[l] == 'gangexpo':
                gmod.namepara.genrelem[l] = ['gang', 'aang']
                gmod.listscalparagenrelem[l] = ['expo', 'self']
            elif gmod.typespatdist[l] == 'glc3':
                gmod.namepara.genrelem[l] = ['dglc', 'thet', 'phii']
                gmod.listscalparagenrelem[l] = ['powr', 'self', 'self']
            else:
                gmod.namepara.genrelem[l] = ['lgal', 'bgal']
                gmod.listscalparagenrelem[l] = ['self', 'self']
            
            # amplitude
            if gmod.typeelem[l] == 'lghtpntsagnntrue':
                gmod.namepara.genrelem[l] += ['lum0']
                gmod.listscalparagenrelem[l] += ['dpowslopbrek']
            elif gmod.typeelem[l] == 'lghtpntspuls':
                gmod.namepara.genrelem[l] += ['per0']
                gmod.listscalparagenrelem[l] += ['lnormeanstdv']
            elif gmod.typeelem[l].startswith('lght'):
                gmod.namepara.genrelem[l] += ['flux']
                gmod.listscalparagenrelem[l] += [gmod.typeprioflux[l]]
            elif gmod.typeelem[l] == 'lens':
                gmod.namepara.genrelem[l] += ['defs']
                gmod.listscalparagenrelem[l] += ['powr']
            elif gmod.typeelem[l].startswith('clus'):
                gmod.namepara.genrelem[l] += ['nobj']
                gmod.listscalparagenrelem[l] += ['powr']
           
            # shape
            if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
                gmod.namepara.genrelem[l] += ['gwdt']
                gmod.listscalparagenrelem[l] += ['powr']
            if gmod.typeelem[l] == 'lghtlinevoig':
                gmod.namepara.genrelem[l] += ['sigm']
                gmod.listscalparagenrelem[l] += ['logt']
                gmod.namepara.genrelem[l] += ['gamm']
                gmod.listscalparagenrelem[l] += ['logt']
            
            # others
            if gmod.typeelem[l] == 'lghtpntspuls':
                gmod.namepara.genrelem[l] += ['magf']
                gmod.listscalparagenrelem[l] += ['lnormeanstdv']
                gmod.namepara.genrelem[l] += ['geff']
                gmod.listscalparagenrelem[l] += ['self']
            elif gmod.typeelem[l] == 'lghtpntsagnntrue':
                gmod.namepara.genrelem[l] += ['dlos']
                gmod.listscalparagenrelem[l] += ['powr']

            if gdat.numbener > 1 and gmod.typeelem[l].startswith('lghtpnts'):
                if gmod.spectype[l] == 'colr':
                    for i in gdat.indxener:
                        if i == 0:
                            continue
                        gmod.namepara.genrelem[l] += ['sindcolr%04d' % i]
                        gmod.listscalparagenrelem[l] += ['self']
                else:
                    gmod.namepara.genrelem[l] += ['sind']
                    gmod.listscalparagenrelem[l] += ['self']
                    if gmod.spectype[l] == 'curv':
                        gmod.namepara.genrelem[l] += ['curv']
                        gmod.listscalparagenrelem[l] += ['self']
                    if gmod.spectype[l] == 'expc':
                        gmod.namepara.genrelem[l] += ['expc']
                        gmod.listscalparagenrelem[l] += ['self']
            if gmod.typeelem[l] == 'lens':
                if gdat.variasca:
                    gmod.namepara.genrelem[l] += ['asca']
                    gmod.listscalparagenrelem[l] += ['self']
                if gdat.variacut:
                    gmod.namepara.genrelem[l] += ['acut']
                    gmod.listscalparagenrelem[l] += ['self']
        
        # names of element parameters for each scaling
        gmod.namepara.genrelemscal = [{} for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            for scaltype in gdat.listscaltype:
                gmod.namepara.genrelemscal[l][scaltype] = []
                for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    if scaltype == gmod.listscalparagenrelem[l][k]:
                        gmod.namepara.genrelemscal[l][scaltype].append(nameparagenrelem)

        # variables for which whose marginal distribution and pair-correlations will be plotted
        gmod.namepara.derielemodim = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            gmod.namepara.derielemodim[l] = deepcopy(gmod.namepara.genrelem[l])
            gmod.namepara.derielemodim[l] += ['deltllik']
            if gdat.boolbinsspat:
                if not 'lgal' in gmod.namepara.derielemodim[l]:
                    gmod.namepara.derielemodim[l] += ['lgal']
                if not 'bgal' in gmod.namepara.derielemodim[l]:
                    gmod.namepara.derielemodim[l] += ['bgal']
                if not 'gang' in gmod.namepara.derielemodim[l]:
                    gmod.namepara.derielemodim[l] += ['gang']
                if not 'aang' in gmod.namepara.derielemodim[l]:
                    gmod.namepara.derielemodim[l] += ['aang']
            if gmod.typeelem[l].startswith('lght'):
                gmod.namepara.derielemodim[l] += ['cnts']
                if gdat.typeexpr == 'ferm':
                    gmod.namepara.derielemodim[l] + ['sbrt0018']
                
            if gmod.typeelem[l] == 'lghtpntsagnntrue':
                gmod.namepara.derielemodim[l] += ['reds']
                gmod.namepara.derielemodim[l] += ['lumi']
                gmod.namepara.derielemodim[l] += ['flux']
            if gmod.typeelem[l] == 'lghtpntspuls':
                gmod.namepara.derielemodim[l] += ['lumi']
                gmod.namepara.derielemodim[l] += ['flux']
                gmod.namepara.derielemodim[l] += ['mass']
                gmod.namepara.derielemodim[l] += ['dlos']
            if gmod.typeelem[l] == 'lens':
                gmod.namepara.derielemodim[l] += ['mcut', 'diss', 'rele', 'reln', 'relk', 'relf', 'relm', 'reld', 'relc']
        
            #for k in range(len(gmod.namepara.derielemodim[l])):
            #    gmod.namepara.derielemodim[l][k] += 'pop%d' % l
            
            # check later
            # temp
            #if strgmodl == 'fitt':
            #    for q in gdat.indxrefr: 
            #        if gmod.nameparagenrelemampl[l] in gdat.refr.namepara.elem[q]:
            #            gmod.namepara.derielemodim[l].append('aerr' + gdat.listnamerefr[q])
        
        if gdat.typeverb > 1:
            print('gmod.namepara.derielemodim')
            print(gmod.namepara.derielemodim)
        
        # derived element parameters
        gmod.namepara.derielem = gmod.namepara.derielemodim[:]
        
        if gdat.typeverb > 1:
            print('gmod.namepara.derielem')
            print(gmod.namepara.derielem)
        
        # derived parameters
        gmod.listnameparaderitotl = [temptemp for temp in gmod.namepara.derielem for temptemp in temp]
        #gmod.listnameparaderitotl += gmod.namepara.scal
        
        for namediff in gmod.listnamediff:
            gmod.listnameparaderitotl += ['cntp' + namediff]
        
        if gdat.typeverb > 1:
            print('gmod.listnameparaderitotl')
            print(gmod.listnameparaderitotl)

        if strgmodl == 'fitt':
            # add reference element parameters that are not available in the fitting model
            gdat.refr.namepara.elemonly = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
            gmod.namepara.extrelem = [[] for l in gmod.indxpopl]
            for q in gdat.indxrefr: 
                if gdat.refr.numbelem[q] == 0:
                    continue
                for name in gdat.refr.namepara.elem[q]:
                    for l in gmod.indxpopl:
                        if gmod.typeelem[l].startswith('lght') and (name == 'defs' or name == 'acut' or name == 'asca' or name == 'mass'):
                            continue
                        if gmod.typeelem[l] == ('lens') and (name == 'cnts' or name == 'flux' or name == 'spec' or name == 'sind'):
                            continue
                        if not name in gmod.namepara.derielemodim[l]:
                            nametotl = name + gdat.listnamerefr[q]
                            if name == 'etag':
                                continue
                            gmod.namepara.derielemodim[l].append(nametotl)
                            
                            if gdat.refr.numbelem[q] == 0:
                                continue

                            gdat.refr.namepara.elemonly[q][l].append(name)
                            if not nametotl in gmod.namepara.extrelem[l]:
                                gmod.namepara.extrelem[l].append(nametotl) 
                            #if name == 'reds':
                            #    for nametemp in ['lumi', 'dlos']:
                            #        nametemptemp = nametemp + gdat.listnamerefr[q]
                            #        if not nametemptemp in gmod.namepara.extrelem[l]:
                            #            gmod.namepara.derielemodim[l].append(nametemp + gdat.listnamerefr[q])
                            #            gmod.namepara.extrelem[l].append(nametemptemp)
            
            if gdat.typeverb > 1:
                print('gdat.refr.namepara.elemonly')
                print(gdat.refr.namepara.elemonly)
        
            if gdat.typeexpr == 'chan' and gdat.typedata == 'inpt':
                for l in gmod.indxpopl:
                    if gmod.typeelem[l] == 'lghtpnts':
                        gmod.namepara.extrelem[l].append('lumiwo08')
                        gmod.namepara.derielemodim[l].append('lumiwo08')
            
            if gdat.typeverb > 1:
                print('gmod.namepara.extrelem')
                print(gmod.namepara.extrelem)

        # defaults
        gmod.liststrgpdfnmodu = [[] for l in gmod.indxpopl]
        gmod.namepara.genrelemmodu = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght'): 
                if gdat.typeexpr == 'ferm' and gdat.lgalcntr == 0.:
                    if l == 1:
                        gmod.liststrgpdfnmodu[l] += ['tmplnfwp']
                        gmod.namepara.genrelemmodu[l] += ['lgalbgal']
                    if l == 2:
                        gmod.liststrgpdfnmodu[l] += ['tmplnfwp']
                        gmod.namepara.genrelemmodu[l] += ['lgalbgal']
        
        gmod.namepara.elem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            for liststrg in [gmod.namepara.genrelem[l], gmod.namepara.derielemodim[l]]:
                for strgthis in liststrg:
                    if not strgthis in gmod.namepara.elem[l]:
                        gmod.namepara.elem[l].append(strgthis)
        
        # temp
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lghtline'):
                gmod.namepara.genrelem[l] += ['spec']
            if gmod.typeelem[l].startswith('lght'):
                gmod.namepara.genrelem[l] += ['spec', 'specplot']
            if gmod.typeelem[l] == 'lens':
                gmod.namepara.genrelem[l] += ['deflprof']
        
        #gmod.namepara.genrelemeval = [[] for l in gmod.indxpopl]
        #for l in gmod.indxpopl:
        #    if gmod.typeelem[l].startswith('clus'):
        #        gmod.namepara.genrelemeval[l] = ['lgal', 'bgal', 'nobj']
        #    if gmod.typeelem[l] == 'clusvari':
        #        gmod.namepara.genrelemeval[l] += ['gwdt']
        #    if gmod.typeelem[l] == 'lens':
        #        gmod.namepara.genrelemeval[l] = ['lgal', 'bgal', 'defs', 'asca', 'acut']
        #    if gmod.typeelem[l].startswith('lghtline'):
        #        gmod.namepara.genrelemeval[l] = ['elin', 'spec']
        #    elif gmod.typeelem[l] == 'lghtgausbgrd':
        #        gmod.namepara.genrelemeval[l] = ['lgal', 'bgal', 'gwdt', 'spec']
        #    elif gmod.typeelem[l].startswith('lght'):
        #        gmod.namepara.genrelemeval[l] = ['lgal', 'bgal', 'spec']
        
        ## element legends
        lablpopl = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gdat.numbgrid > 1:
                if gmod.typeelem[l] == 'lghtpnts':
                    lablpopl[l] = 'FPS'
                if gmod.typeelem[l] == 'lghtgausbgrd':
                    lablpopl[l] = 'BGS'
            else:
                if gmod.typeelem[l] == 'lghtpntspuls':
                    lablpopl[l] = 'Pulsar'
                elif gmod.typeelem[l].startswith('lghtpntsagnn'):
                    lablpopl[l] = 'AGN'
                elif gmod.typeelem[l].startswith('lghtpnts'):
                    lablpopl[l] = 'PS'
            if gmod.typeelem[l] == 'lens':
                lablpopl[l] = 'Subhalo'
            if gmod.typeelem[l].startswith('clus'):
                lablpopl[l] = 'Cluster'
            if gmod.typeelem[l].startswith('lghtline'):
                lablpopl[l]= 'Line'
        setp_varb(gdat, 'lablpopl', valu=lablpopl, strgmodl=strgmodl)

        if strgmodl == 'true':
            gmod.indxpoplassc = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                if gmod.numbpopl == 3 and gmod.typeelem[1] == 'lens':
                    gmod.indxpoplassc[l] = [l]
                else:
                    gmod.indxpoplassc[l] = gmod.indxpopl

        # variables for which two dimensional histograms will be calculated
        gmod.namepara.genrelemcorr = [[] for l in gmod.indxpopl]
        if gdat.boolplotelemcorr:
            for l in gmod.indxpopl:
                for strgfeat in gmod.namepara.derielemodim[l]:
                    gmod.namepara.genrelemcorr[l].append(strgfeat)
        
        # number of element parameters
        if gmod.numbpopl > 0:
            gmod.numbparagenrelemsing = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparaderielemsing = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparaelemsing = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparagenrelem = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparagenrelemcuml = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparagenrelemcumr = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparaderielem = np.zeros(gmod.numbpopl, dtype=int)
            gmod.numbparaelem = np.zeros(gmod.numbpopl, dtype=int)
            for l in gmod.indxpopl:
                # number of generative element parameters for a single element of a specific population
                gmod.numbparagenrelemsing[l] = len(gmod.namepara.genrelem[l])
                # number of derived element parameters for a single element of a specific population
                gmod.numbparaderielemsing[l] = len(gmod.namepara.derielem[l])
                # number of element parameters for a single element of a specific population
                gmod.numbparaelemsing[l] = len(gmod.namepara.elem[l])
                # number of generative element parameters for all elements of a specific population
                gmod.numbparagenrelem[l] = gmod.numbparagenrelemsing[l] * gmod.maxmpara.numbelem[l]
                # number of generative element parameters up to the beginning of a population
                gmod.numbparagenrelemcuml[l] = np.sum(gmod.numbparagenrelem[:l])
                # number of generative element parameters up to the end of a population
                gmod.numbparagenrelemcumr[l] = np.sum(gmod.numbparagenrelem[:l+1])
                # number of derived element parameters for all elements of a specific population
                gmod.numbparaderielem[l] = gmod.numbparaderielemsing[l] * gmod.numbelem[l]
                # number of element parameters for all elements of a specific population
                gmod.numbparaelem[l] = gmod.numbparaelemsing[l] * gmod.numbelem[l]
            # number of generative element parameters summed over all populations
            gmod.numbparagenrelemtotl = np.sum(gmod.numbparagenrelem)
            # number of derived element parameters summed over all populations
            gmod.numbparaderielemtotl = np.sum(gmod.numbparaderielem)
            # number of element parameters summed over all populations
            gmod.numbparaelemtotl = np.sum(gmod.numbparaderielem)
        
        gmod.indxparagenrelemsing = []
        for l in gmod.indxpopl:
            gmod.indxparagenrelemsing.append(np.arange(gmod.numbparagenrelemsing[l]))
        
        gmod.indxparaderielemsing = []
        for l in gmod.indxpopl:
            gmod.indxparaderielemsing.append(np.arange(gmod.numbparaderielemsing[l]))
        
        gmod.indxparaelemsing = []
        for l in gmod.indxpopl:
            gmod.indxparaelemsing.append(np.arange(gmod.numbparaelemsing[l]))

        # size of the auxiliary variable propobability density vector
        if gmod.maxmpara.numbelemtotl > 0:
            gmod.numblpri = 3 + gmod.numbparagenrelem * gmod.numbpopl
        else:
            gmod.numblpri = 0
        if gdat.penalpridiff:
            gmod.numblpri += 1
        indxlpri = np.arange(gmod.numblpri)

        # append the population tags to element parameter names
        #for l in gmod.indxpopl:
        #    gmod.namepara.genrelem[l] = [gmod.namepara.genrelem[l][g] + 'pop%d' % l for g in gmod.indxparagenrelemsing[l]]
        
        gmod.boolcompposi = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            gmod.boolcompposi[l] = np.zeros(gmod.numbparagenrelemsing[l], dtype=bool)
            if gmod.typeelem[l].startswith('lghtline'):
                gmod.boolcompposi[l][0] = True
            else:
                gmod.boolcompposi[l][0] = True
                gmod.boolcompposi[l][1] = True
        
        # list of strings across all populations
        ## all (generative and derived) element parameters
        gmod.numbparaelem = len(gmod.namepara.elem)
        gmod.indxparaelem = np.arange(gmod.numbparaelem)
        
        # flattened list of generative element parameters
        gmod.listnameparagenfelem = []
        for l in gmod.indxpopl:
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                gmod.listnameparagenfelem.append(nameparagenrelem + 'pop%d' % l)
        
        # concatenated list of flattened generative and derived element parameters
        gmod.listnameparatotlelem = gmod.listnameparagenfelem + gmod.namepara.derielem

        gmod.numbparaelem = np.empty(gmod.numbpopl, dtype=int)
        for l in gmod.indxpopl:
            gmod.numbparaelem[l] = len(gmod.namepara.elem[l])
        
        numbdeflsubhplot = 2
        numbdeflsingplot = numbdeflsubhplot
        if gmod.numbparaelem > 0:
            numbdeflsingplot += 3

        gmod.convdiffanyy = True in convdiff

        cntr = tdpy.cntr()
        
        if gmod.boollens:
            adishost = gdat.adisobjt(redshost)
            adissour = gdat.adisobjt(redssour)
            adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
            massfrombein = retr_massfrombein(gdat, adissour, adishost, adishostsour)
            mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
        
        # object of parameter indices
        gmod.indxpara = tdpy.gdatstrt()
        
        # define parameter indices
        if gmod.numbparaelem > 0:

            # number of elements
            #gmod.indxpara.numbelem = np.empty(gmod.numbpopl, dtype=int)
            for l in gmod.indxpopl:
                indx = cntr.incr()
                setattr(gmod.indxpara, 'numbelempop%d' % l, indx)
                #gmod.indxpara.numbelem[l] = indx
            
            # hyperparameters
            ## mean number of elements
            if gmod.typemodltran == 'pois':
                #gmod.indxpara.meanelem = np.empty(gmod.numbpopl, dtype=int)
                for l in gmod.indxpopl:
                    if gmod.maxmpara.numbelem[l] > 0:
                        indx = cntr.incr()
                        setattr(gmod.indxpara, 'meanelempop%d' % l, indx)
                        #gmod.indxpara.meanelem[l] = indx

            ## parameters parametrizing priors on element parameters
            liststrgvarb = []
            for l in gmod.indxpopl:
                if gmod.maxmpara.numbelem[l] > 0:
                    for strgpdfnelemgenr, strgfeat in zip(gmod.listscalparagenrelem[l], gmod.namepara.genrelem[l]):
                        if strgpdfnelemgenr == 'expo' or strgpdfnelemgenr == 'dexp':
                            liststrgvarb += [strgfeat + 'distscal']
                        if strgpdfnelemgenr == 'powr':
                            liststrgvarb += ['slopprio' + strgfeat + 'pop%d' % l]
                        if strgpdfnelemgenr == 'dpow':
                            liststrgvarb += [strgfeat + 'distbrek']
                            liststrgvarb += [strgfeat + 'sloplowr']
                            liststrgvarb += [strgfeat + 'slopuppr']
                        if strgpdfnelemgenr == 'gausmean' or strgpdfnelemgenr == 'lnormean':
                            liststrgvarb += [strgfeat + 'distmean']
                        if strgpdfnelemgenr == 'gausstdv' or strgpdfnelemgenr == 'lnorstdv':
                            liststrgvarb += [strgfeat + 'diststdv']
                        if strgpdfnelemgenr == 'gausmeanstdv' or strgpdfnelemgenr == 'lnormeanstdv':
                            liststrgvarb += [nameparagenrelem + 'distmean', nameparagenrelem + 'diststdv']
            for strgvarb in liststrgvarb:
                setattr(gmod.indxpara, strgvarb, np.zeros(gmod.numbpopl, dtype=int) - 1)

            for l in gmod.indxpopl:
                strgpopl = 'pop%d' % l
                if gmod.maxmpara.numbelem[l] > 0:
                    for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                        
                        if gmod.listscalparagenrelem[l][k] == 'self':
                            continue
                        indx = cntr.incr()

                        if gmod.listscalparagenrelem[l][k] == 'dpow':
                            for nametemp in ['brek', 'sloplowr', 'slopuppr']:
                                strg = '%s' % nametemp + nameparagenrelem
                                setattr(gmod.indxpara, strg, indx)
                                setattr(gmod.indxpara, strg, indx)
                        else:
                            if gmod.listscalparagenrelem[l][k] == 'expo' or gmod.listscalparagenrelem[l][k] == 'dexp':
                                strghypr = 'scal'
                            if gmod.listscalparagenrelem[l][k] == 'powr':
                                strghypr = 'slop'
                            if gmod.listscalparagenrelem[l][k] == 'gausmean' or gmod.listscalparagenrelem[l][k] == 'gausmeanstdv' or \
                                            gmod.listscalparagenrelem[l][k] == 'lnormean' or gmod.listscalparagenrelem[l][k] == 'lnormeanstdv':
                                strghypr = 'mean'
                            if gmod.listscalparagenrelem[l][k] == 'gausstdv' or gmod.listscalparagenrelem[l][k] == 'gausmeanstdv' or \
                                            gmod.listscalparagenrelem[l][k] == 'lnorstdv' or gmod.listscalparagenrelem[l][k] == 'lnormeanstdv':
                                strghypr = 'stdv'
                            strg = strghypr + 'prio' + nameparagenrelem + 'pop%d' % l
                            setattr(gmod.indxpara, strg, indx)
        
        # group PSF parameters
        if gmod.typeevalpsfn == 'kern' or gmod.typeevalpsfn == 'full':
            for m in gdat.indxevtt:
                for i in gdat.indxener:
                    setattr(gmod.indxpara, 'sigcen%02devt%d' % (i, m), cntr.incr())
                    if gmod.typemodlpsfn == 'doubking' or gmod.typemodlpsfn == 'singking':
                        setattr(gmod.indxpara, 'gamcen%02devt%d' % (i, m), cntr.incr())
                        if gmod.typemodlpsfn == 'doubking':
                            setattr(gmod.indxpara, 'sigten%02devt%d' % (i, m), cntr.incr())
                            setattr(gmod.indxpara, 'gamten%02devt%d' % (i, m), cntr.incr())
                            setattr(gmod.indxpara, 'ffenen%02devt%d' % (i, m), cntr.incr())
            
            gmod.indxpara.psfp = []
            for strg, valu in gmod.indxpara.__dict__.items():
                if strg.startswith('sigce') or strg.startswith('sigte') or strg.startswith('gamce') or strg.startswith('gamte') or strg.startswith('psffe'):
                    gmod.indxpara.psfp.append(valu)
            gmod.indxpara.psfp = np.array(gmod.indxpara.psfp) 

            gmod.numbpsfptotlevtt = gdat.numbevtt * gmod.numbpsfptotl
            gmod.numbpsfptotlener = gdat.numbener * gmod.numbpsfptotl
            numbpsfp = gmod.numbpsfptotl * gdat.numbener * gdat.numbevtt
            indxpsfpform = np.arange(numbpsfpform)
            indxpsfptotl = np.arange(gmod.numbpsfptotl)
   
            indxpsfp = np.arange(numbpsfp)
            gmod.indxpara.psfp = np.sort(gmod.indxpara.psfp)
            gmod.indxparapsfpinit = gmod.indxpara.psfp[0]
        
        # group background parameters
        gmod.indxpara.bacp = []
        for c in gmod.indxback:
            if gmod.boolspecback[c]:
                indx = cntr.incr()
                setattr(gmod.indxpara, 'bacpback%04d' % c, indx)
                gmod.indxpara.bacp.append(indx)
            else:
                for i in gdat.indxener:
                    indx = cntr.incr()
                    setattr(gmod.indxpara, 'bacpback%04den%02d' % (c, i), indx)
                    gmod.indxpara.bacp.append(indx)
        gmod.indxpara.bacp = np.array(gmod.indxpara.bacp)

        # temp
        #gmod.indxpara.anglsour = []
        #gmod.indxpara.anglhost = []
        #gmod.indxpara.angllens = []
        
        if gmod.typeemishost != 'none':
            gmod.indxpara.specsour = []
            gmod.indxpara.spechost = []

        if gmod.boollens:
            gmod.indxpara.lgalsour = cntr.incr()
            gmod.indxpara.bgalsour = cntr.incr()
            gmod.indxpara.fluxsour = cntr.incr()
            if gdat.numbener > 1:
                gmod.indxpara.sindsour = cntr.incr()
            gmod.indxpara.sizesour = cntr.incr()
            gmod.indxpara.ellpsour = cntr.incr()
            gmod.indxpara.anglsour = cntr.incr()
        if gmod.typeemishost != 'none' or gmod.boollens:
            for e in gmod.indxsersfgrd: 
                if gmod.typeemishost != 'none':
                    setattr(gmod.indxpara, 'lgalhostisf%d' % e, cntr.incr())
                    setattr(gmod.indxpara, 'bgalhostisf%d' % e, cntr.incr())
                    setattr(gmod.indxpara, 'fluxhostisf%d' % e, cntr.incr())
                    if gdat.numbener > 1:
                        setattr(gmod.indxpara, 'sindhostisf%d' % e, cntr.incr())
                    setattr(gmod.indxpara, 'sizehostisf%d' % e, cntr.incr())
                if gmod.boollens:
                    setattr(gmod.indxpara, 'beinhostisf%d' % e, cntr.incr())
                if gmod.typeemishost != 'none':
                    setattr(gmod.indxpara, 'ellphostisf%d' % e, cntr.incr())
                    setattr(gmod.indxpara, 'anglhostisf%d' % e, cntr.incr())
                    setattr(gmod.indxpara, 'serihostisf%d' % e, cntr.incr())
        if gmod.boollens:
            gmod.indxpara.sherextr = cntr.incr()
            gmod.indxpara.sangextr = cntr.incr()
            gmod.indxpara.sour = []
        
        if gmod.boollens and gmod.typeemishost == 'none':
            raise Exception('Lensing cannot be modeled without host galaxy emission.')
    
        # collect groups of parameters
        if gdat.typeexpr == 'hubb':
            gmod.listnamecomplens = ['hostlght', 'hostlens', 'sour', 'extr']
            for namecomplens in gmod.listnamecomplens:
                setattr(gmod, 'liststrg' + namecomplens, [])
                setattr(gmod.indxpara, namecomplens, [])
        if gmod.boollens or gmod.typeemishost != 'none':
            gmod.liststrghostlght += ['lgalhost', 'bgalhost', 'ellphost', 'anglhost']
            gmod.liststrghostlens += ['lgalhost', 'bgalhost', 'ellphost', 'anglhost']
        if gmod.typeemishost != 'none':
            gmod.liststrghostlght += ['fluxhost', 'sizehost', 'serihost']
            if gdat.numbener > 1:
                gmod.liststrghostlght += ['sindhost']
        if gmod.boollens:
            gmod.liststrghostlens += ['beinhost']
            gmod.liststrgextr += ['sherextr', 'sangextr']
            gmod.liststrgsour += ['lgalsour', 'bgalsour', 'fluxsour', 'sizesour', 'ellpsour', 'anglsour']
            if gdat.numbener > 1:
                gmod.liststrgsour += ['sindsour']
        
        for strg, valu in gmod.__dict__.items():
            
            if isinstance(valu, list) or isinstance(valu, np.ndarray):
                continue
            
            if gdat.typeexpr == 'hubb':
                for namecomplens in gmod.listnamecomplens:
                    for strgtemp in getattr(gmod, 'liststrg' + namecomplens):
                        if strg[12:].startswith(strgtemp):
                            
                            if isinstance(valu, list):
                                for valutemp in valu:
                                    gmod['indxparagenr' + namecomplens].append(valutemp)
                            else:
                                gmod['indxparagenr' + namecomplens].append(valu)
            
            # remove indxpara. from strg
            strg = strg[12:]
            
            if strg.startswith('fluxsour') or strg.startswith('sindsour'):
                gmod.indxpara.specsour.append(valu)

            if strg.startswith('fluxhost') or strg.startswith('sindhost'):
                gmod.indxpara.spechost.append(valu)
        
        if gmod.boollens or gmod.boolhost:
            gmod.indxpara.host = gmod.indxparahostlght + gmod.indxparahostlens
            gmod.indxpara.lens = gmod.indxpara.host + gmod.indxpara.sour + gmod.indxpara.extr

        ## number of model spectral parameters for each population
        #numbspep = np.empty(gmod.numbpopl, dtype=int)
        #liststrgspep = [[] for l in range(gmod.numbpopl)]
        #for l in gmod.indxpopl:
        #    if gdat.numbener > 1:
        #        liststrgspep[l] += ['sind']
        #        if gmod.spectype[l] == 'expc':
        #            liststrgspep[l] += ['expc']
        #        if gmod.spectype[l] == 'curv':
        #            liststrgspep[l] = ['curv']
        #    numbspep[l] = len(liststrgspep[l]) 
        

def setp_paragenrscalbase(gdat, strgmodl='fitt'):
    '''
    Setup labels and scales for base parameters
    '''
    
    print('setp_paragenrscalbase(): Building the %s model base paremeter names and scales...' % strgmodl)
    gmod = getattr(gdat, strgmodl)
    
    listlablback = []
    listlablback = []
    for nameback in gmod.listnameback:
        if nameback == 'isot':
            listlablback.append('Isotropic')
            listlablback.append(r'$\mathcal{I}$')
        if nameback == 'fdfm':
            listlablback.append('FDM')
            listlablback.append(r'$\mathcal{D}$')
        if nameback == 'dark':
            listlablback.append('NFW')
            listlablback.append(r'$\mathcal{D}_{dark}$')
        if nameback == 'part':
            listlablback.append('Particle Back.')
            listlablback.append(r'$\mathcal{I}_p$')

    # background templates
    listlablsbrt = deepcopy(listlablback)
    numblablsbrt = 0
    for l in gmod.indxpopl:
        if gmod.boolelemsbrt[l]:
            listlablsbrt.append(gmod.lablpopl[l])
            listlablsbrt.append(gmod.lablpopl[l] + ' subt')
            numblablsbrt += 2
    if gmod.boollens:
        listlablsbrt.append('Source')
        numblablsbrt += 1
    if gmod.typeemishost != 'none':
        for e in gmod.indxsersfgrd:
            listlablsbrt.append('Host %d' % e)
            numblablsbrt += 1
    if gmod.numbpopl > 0:
        if 'clus' in gmod.typeelem or 'clusvari' in gmod.typeelem:
            listlablsbrt.append('Uniform')
            numblablsbrt += 1
    
    listlablsbrtspec = ['Data']
    listlablsbrtspec += deepcopy(listlablsbrt)
    if len(listlablsbrt) > 1:
        listlablsbrtspec.append('Total Model')
    
    numblablsbrtspec = len(listlablsbrtspec)
    
    # number of generative parameters per element, depends on population
    #numbparaelem = gmod.numbparagenrelem + numbparaelemderi

    # maximum total number of parameters
    #numbparagenrfull = gmod.numbparagenrbase + gmod.numbparaelem
    
    #numbparaelemkind = gmod.numbparagenrbase
    #for l in gmod.indxpopl:
    #    numbparaelemkind += gmod.numbparagenrelemsing[l]
    
    #nameparagenrbase
    #gmod.namepara.genrelem
    
    #listnameparaderifixd
    #listnameparaderielem
    
    #gmod.namepara.genrelemextd = gmod.namepara.genrelem * maxm.numbelem
    #listnameparaderielemextd = gmod.namepara.genrelem * maxm.numbelem
    
    gmod.listindxparakindscal = {}
    for scaltype in gdat.listscaltype:
        gmod.listindxparakindscal[scaltype] = np.where(scaltype == gmod.listscalparakind)[0]

    #
    ## stack
    ## gmod.listnameparastck
    #gmod.listnameparastck = np.zeros(gmod.maxmnumbpara, dtype=object)
    #gmod.listscalparastck = np.zeros(gmod.maxmnumbpara, dtype=object)
    #
    #gmod.listnameparastck[gmod.indxparagenrbase] = gmod.nameparagenrbase
    #gmod.listscalparastck[gmod.indxparagenrbase] = gmod.listscalparagenrbase
    #for k in range(gmod.numbparaelem):
    #    for l in gmod.indxpopl:  
    #        if k >= gmod.numbparagenrelemcuml[l]:
    #            indxpopltemp = l
    #            indxelemtemp = (k - gmod.numbparagenrelemcuml[indxpopltemp]) // gmod.numbparagenrelemsing[indxpopltemp]
    #            gmod.indxparagenrelemtemp = (k - gmod.numbparagenrelemcuml[indxpopltemp]) % gmod.numbparagenrelemsing[indxpopltemp]
    #            break
    #    gmod.listnameparastck[gmod.numbparagenrbase+k] = '%spop%d%04d' % (gmod.namepara.genrelem[indxpopltemp][gmod.indxparagenrelemtemp], indxpopltemp, indxelemtemp)
    #    gmod.listscalparastck[gmod.numbparagenrbase+k] = gmod.listscalparagenrelem[indxpopltemp][gmod.indxparagenrelemtemp]
    #
    #
    #if np.where(gmod.listscalpara == 0)[0].size > 0:
    #    print('gmod.listscalpara[gmod.indxparagenrbase]')
    #    print(gmod.listscalpara[gmod.indxparagenrbase])
    #    raise Exception('')
    #
    ## labels and scales for variables
    if gmod.boollens:
        setattr(gmod.lablrootpara, 'masssubhintg', r'$M_{\rm{sub}}$')
        setattr(gmod.lablrootpara, 'masssubhdelt', r'$\rho_{\rm{sub}}$')
        setattr(gmod.lablrootpara, 'masssubhintgbein', r'$M_{\rm{sub,E}}$')
        setattr(gmod.lablrootpara, 'masssubhdeltbein', r'$\rho_{\rm{sub,E}}$')
        setattr(gmod.lablrootpara, 'masssubhintgunit', '$10^9 M_{\odot}$')
        setattr(gmod.lablrootpara, 'masssubhdeltunit', '$M_{\odot}$/kpc')
        setattr(gmod.lablrootpara, 'masssubhintgbeinunit', '$10^9 M_{\odot}$')
        setattr(gmod.lablrootpara, 'masssubhdeltbeinunit', '$M_{\odot}$/kpc')
        setattr(gmod.lablrootpara, 'fracsubhintg', r'f_{\rm{sub}}')
        setattr(gmod.lablrootpara, 'fracsubhdelt', r'f_{\rho,\rm{sub}}')
        setattr(gmod.lablrootpara, 'fracsubhintgbein', r'$f_{\rm{sub,E}}$')
        setattr(gmod.lablrootpara, 'fracsubhdeltbein', r'$f_{\rho,\rm{sub,E}}$')
        for e in gmod.indxsersfgrd:
            setattr(gmod.lablrootpara, 'masshostisf%dbein' % e, r'$M_{\rm{hst,%d,C}}$' % e)
            setattr(gmod.lablrootpara, 'masshostisf%dintg' % e, r'$M_{\rm{hst,%d<}}$' % e)
            setattr(gmod.lablrootpara, 'masshostisf%ddelt' % e, r'$M_{\rm{hst,%d}}$' % e)
            setattr(gmod.lablrootpara, 'masshostisf%dintgbein' % e, r'$M_{\rm{hst,E,%d<}}$' % e)
            setattr(gmod.lablrootpara, 'masshostisf%ddeltbein' % e, r'$M_{\rm{hst,E,%d}}$' % e)
        for namevarb in ['fracsubh', 'masssubh']:
            for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                for nameeval in ['', 'bein']:
                    setattr(gdat, 'scal' + namevarb + strgcalcmasssubh + nameeval, 'logt')
        for e in gmod.indxsersfgrd:
            setattr(gdat, 'scalmasshostisf%d' % e + 'bein', 'logt')
            for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                for nameeval in ['', 'bein']:
                    setattr(gdat, 'scalmasshostisf%d' % e + strgcalcmasssubh + nameeval, 'logt')
    
    # scalar variable setup
    gdat.lablhistcntplowrdfncsubten00evt0 = 'N_{pix,l}'
    gdat.lablhistcntphigrdfncsubten00evt0 = 'N_{pix,h}'
    gdat.lablhistcntplowrdfncen00evt0 = 'N_{pix,l}'
    gdat.lablhistcntphigrdfncen00evt0 = 'N_{pix,h}'
    
    gdat.lablbooldfncsubt = 'H'
    
    gdat.lablpriofactdoff = r'$\alpha_{p}$'
    gmod.scalpriofactdoff = 'self'

    gdat.minmreds = 0.
    gdat.maxmreds = 1.5
    
    gdat.minmmagt = 19.
    gdat.maxmmagt = 28.

    gmod.scalpara.numbelem = 'logt'
    gmod.scalpara.lliktotl = 'logt'

    gdat.lablener = 'E'
    #gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    
    # width of the Gaussian clusters
    gdat.lablgwdt = r'\sigma_G'
    
    gdat.lablgang = r'\theta'
    gdat.lablaang = r'\phi'
    gdat.labllgalunit = gdat.lablgangunit
    gdat.lablbgalunit = gdat.lablgangunit
   
    gdat.lablanglfromhost = r'\theta_{\rm{0,hst}}'
    gdat.lablanglfromhostunit = gdat.lablgangunit

    gdat.labldefs = r'\alpha_s'
    gdat.lablflux = 'f'
    gdat.lablnobj = 'p'
    
    gdat.lablelin = r'\mathcal{E}'
    
    gdat.lablsbrt = r'\Sigma'
    
    gdat.labldeflprof = r'\alpha_a'
    gdat.labldeflprofunit = u'$^{\prime\prime}$'
    
    gdat.strgenerkevv = 'keV'
    gdat.strgenergevv = 'GeV'
    gdat.strgenerergs = 'erg'
    gdat.strgenerimum = '\mu m^{-1}'

    gdat.labldefsunit = u'$^{\prime\prime}$'
    gdat.lablprat = 'cm$^{-2}$ s$^{-1}$'
    

    ### labels for derived fixed dimensional parameters
    if gdat.boolbinsener:
        for i in gdat.indxener:
            setattr(gmod.lablrootpara, 'fracsdenmeandarkdfncsubten%02d' % i, 'f_{D/ST,%d}' % i)
    else:
        gmod.lablrootpara.fracsdenmeandarkdfncsubt = 'f_{D/ST}'
        setattr(gmod.lablrootpara, 'fracsdenmeandarkdfncsubt', 'f_{D/ST}')
    
    ### labels for background units
    if gdat.typeexpr == 'ferm':
        for nameenerscaltype in ['en00', 'en01', 'en02', 'en03']:
            
            for labltemptemp in ['flux', 'sbrt']:

                # define the label
                if nameenerscaltype == 'en00':
                    strgenerscal = '%s' % labltemp
                if nameenerscaltype == 'en01':
                    strgenerscal = 'E%s' % labltemp
                if nameenerscaltype == 'en02':
                    strgenerscal = 'E^2%s' % labltemp
                if nameenerscaltype == 'en03':
                    strgenerscal = '%s' % labltemp
                labl = '%s' % strgenerscal

                for nameenerunit in ['gevv', 'ergs', 'kevv', 'imum']:
                    
                    strgenerunit = getattr(gdat, 'strgener' + nameenerunit)

                    if nameenerscaltype == 'en00':
                        strgenerscalunit = '%s$^{-1}$' % strgenerunit
                    if nameenerscaltype == 'en01':
                        strgenerscalunit = '' 
                    if nameenerscaltype == 'en02':
                        strgenerscalunit = '%s' % strgenerunit
                    if nameenerscaltype == 'en03':
                        strgenerscalunit = '%s' % strgenerunit
                    
                    # define the label unit
                    for namesoldunit in ['ster', 'degr']:
                        if labltemptemp == 'flux':
                            lablunit = '%s %s' % (strgenerscalunit, gdat.lablprat)
                            setattr(gmod.lablunitpara, 'lablflux' + nameenerscaltype + nameenerunit + 'unit', lablunit)
                        else:
                            if namesoldunit == 'ster':
                                lablunit = '%s %s sr$^{-1}$' % (strgenerscalunit, gdat.lablprat)
                            if namesoldunit == 'degr':
                                lablunit = '%s %s deg$^{-2}$' % (strgenerscalunit, gdat.lablprat)
                            setattr(gmod.lablunitpara, 'sbrt' + nameenerscaltype + nameenerunit + namesoldunit + 'unit', lablunit)

        if gdat.boolbinsener:
            gdat.lablfluxunit = getattr(gmod.lablunitpara, 'fluxen00' + gdat.nameenerunit + 'unit')
            gdat.lablsbrtunit = getattr(gmod.lablunitpara, 'sbrten00' + gdat.nameenerunit + 'sterunit')

    gdat.lablexpo = r'$\epsilon$'
    gdat.lablexpounit = 'cm$^2$ s'
    
    gdat.lablprvl = '$p$'
    
    gdat.lablreds = 'z'
    gdat.lablmagt = 'm_R'
    
    gdat.lablper0 = 'P_0'
    gmod.scalper0plot = 'logt'
  
    gdat.labldglc = 'd_{gc}'
    gmod.scaldglcplot = 'logt'
    
    gdat.labldlos = 'd_{los}'
    gmod.scaldlosplot = 'logt'
    if gdat.typeexpr == 'ferm':
        gdat.labldlosunit = 'kpc'
        gdat.labllumi = r'L_{\gamma}'
    if gdat.typeexpr == 'chan':
        gdat.labldlosunit = 'Mpc'
        gdat.labllumi = r'L_{X}'
        gdat.labllum0 = r'L_{X, 0}'
    
    gdat.lablgeff = r'\eta_{\gamma}'
    gmod.scalgeffplot = 'logt'
    
    gmod.scallumiplot = 'logt'
    gdat.labllumiunit = 'erg s$^{-1}$'
    gdat.labllum0unit = 'erg s$^{-1}$'
    
    gdat.lablthet = r'\theta_{gc}'
    gmod.scalthetplot = 'self'
    
    gdat.lablphii = r'\phi_{gc}'
    gmod.scalphiiplot = 'self'
    
    setattr(gmod.lablrootpara, 'magf', 'B')
    setattr(gdat, 'scalmagfplot', 'logt')
    
    setattr(gmod.lablrootpara, 'per1', 'P_1')
    if gdat.typedata == 'inpt':
        gdat.minmpara.per0 = 1e-3
        gdat.maxmpara.per0 = 1e1
        gdat.minmpara.per1 = 1e-20
        gdat.maxmpara.per1 = 1e-10
        gdat.minmpara.per1 = 1e-20
        gdat.maxmpara.per1 = 1e-10
        gdat.minmpara.flux0400 = 1e-1
        gdat.maxmpara.flux0400 = 1e4
    setattr(gdat, 'scalper1plot', 'logt')
    setattr(gmod.lablrootpara, 'flux0400', 'S_{400}')
    setattr(gdat, 'scalflux0400plot', 'logt')
    
    for q in gdat.indxrefr:
        setattr(gmod.lablrootpara, 'aerr' + gdat.listnamerefr[q], '\Delta_{%d}' % q)
    gdat.lablsigm = '\sigma_l'
    gdat.lablgamm = '\gamma_l'

    gdat.lablbcom = '\eta'
    
    gdat.lablinfopost = 'D_{KL}'
    gdat.lablinfopostunit = 'nat'
    gdat.lablinfoprio = 'D_{KL,pr}'
    gdat.lablinfopriounit = 'nat'
    
    gdat.labllevipost = '\ln P(D)'
    gdat.labllevipostunit = 'nat'
    gdat.lablleviprio = '\ln P_{pr}(D)'
    gdat.labllevipriounit = 'nat'
    
    gdat.lablsind = 's'
    if gdat.boolbinsener:
        for i in gdat.indxenerinde:
            setattr(gmod.lablrootpara, 'sindcolr%04d' % i, 's_%d' % i)

    gdat.lablexpcunit = gdat.strgenerunit
    
    gdat.labllliktotl = r'\ln P(D|M)'
    
    gdat.labllpripena = r'\ln P(N)'
    
    gdat.lablasca = r'\theta_s'
    gdat.lablascaunit = gdat.lablgangunit
    gdat.lablacut = r'\theta_c'
    gdat.lablacutunit = gdat.lablgangunit
    
    gdat.lablmcut = r'M_{c,n}'
    gdat.lablmcutunit = r'$M_{\odot}$'
    
    gdat.lablmcutcorr = r'\bar{M}_{c,n}'
    gdat.lablmcutcorrunit = r'$M_{\odot}$'
    
    gdat.lablspec = gdat.lablflux
    gdat.lablspecunit = gdat.lablfluxunit
    gdat.lablspecplot = gdat.lablflux
    gdat.lablspecplotunit = gdat.lablfluxunit
    gdat.lablcnts = 'C'
    gdat.labldeltllik = r'\Delta_n \ln P(D|M)'
    gdat.labldiss = r'\theta_{sa}'
    gdat.labldissunit = gdat.lablgangunit
    
    gdat.lablrele = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_l| \rangle'
    
    gdat.lablrelc = r'\langle\vec{\alpha}_n \cdot \vec{\nabla} k_l \rangle'
    
    gdat.lablreld = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_d| \rangle'
    
    gdat.lablreln = r'\langle \Delta \theta_{pix} |\hat{\alpha}_n \cdot \vec{\nabla} k_l| / \alpha_{s,n} \rangle'
    
    gdat.lablrelm = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle'
    gdat.lablrelk = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle'
    gdat.lablrelf = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle / k_m'
    
    for q in gdat.indxrefr:
        for l in gmod.indxpopl:
            setp_varb(gdat, 'fdispop%dpop%d' % (l, q), minm=0., maxm=1., lablroot='$F_{%d%d}$' % (l, q))
            setp_varb(gdat, 'cmplpop%dpop%d' % (l, q), minm=0., maxm=1., lablroot='$C_{%d%d}$' % (l, q))
                    
    if gdat.typeexpr == 'chan':
        if gdat.anlytype == 'spec':
            gdat.minmspec = 1e-2
            gdat.maxmspec = 1e1
        else:
            gdat.minmspec = 1e-11
            gdat.maxmspec = 1e-7
    else:
        gdat.minmspec = 1e-11
        gdat.maxmspec = 1e-7
    
    if gdat.typeexpr == 'ferm':
        gdat.minmlumi = 1e32
        gdat.maxmlumi = 1e36
    elif gdat.typeexpr == 'chan':
        if gdat.typedata == 'inpt':
            gdat.minmlum0 = 1e42
            gdat.maxmlum0 = 1e46
        gdat.minmlumi = 1e41
        gdat.maxmlumi = 1e45
    
    try:
        gdat.minmdlos
    except:
        if gdat.typeexpr == 'chan':
            gdat.minmdlos = 1e7
            gdat.maxmdlos = 1e9
        else:
            gdat.minmdlos = 6e3
            gdat.maxmdlos = 1.1e4
    
    if gdat.typeexpr == 'ferm':
        gdat.minmcnts = 1e1
        gdat.maxmcnts = 1e5
    if gdat.typeexpr == 'chan':
        if gdat.numbpixlfull == 1:
            gdat.minmcnts = 1e4
            gdat.maxmcnts = 1e8
        else:
            gdat.minmcnts = 1.
            gdat.maxmcnts = 1e3
    if gdat.typeexpr == 'hubb':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3
    if gdat.typeexpr == 'fire':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3

    gdat.minmspecplot = gdat.minmspec
    gdat.maxmspecplot = gdat.maxmspec
    
    gdat.minmdeltllik = 1.
    gdat.maxmdeltllik = 1e3
    gdat.minmdiss = 0.
    gdat.maxmdiss = gdat.maxmgangdata * np.sqrt(2.)
    
    gdat.minmrele = 1e-3
    gdat.maxmrele = 1e1

    gdat.minmreln = 1e-3
    gdat.maxmreln = 1.

    gdat.minmrelk = 1e-3
    gdat.maxmrelk = 1.

    gdat.minmrelf = 1e-5
    gdat.maxmrelf = 1e-1

    gdat.minmrelm = 1e-3
    gdat.maxmrelm = 1e1

    gdat.minmreld = 1e-3
    gdat.maxmreld = 1e1

    gdat.minmrelc = 1e-3
    gdat.maxmrelc = 1.

    gdat.minmmcut = 3e7
    gdat.maxmmcut = 2e9
    gdat.minmmcutcorr = gdat.minmmcut
    gdat.maxmmcutcorr = gdat.maxmmcut

    if gdat.boolbinsspat:
        gdat.minmbein = 0.
        gdat.maxmbein = 1. / gdat.anglfact
    
    # scalar variables
    if gdat.boolbinsspat:
        gdat.minmdeflprof = 1e-3 / gdat.anglfact
        gdat.maxmdeflprof = 0.1 / gdat.anglfact
    
    #gdat.minmfracsubh = 0.
    #gdat.maxmfracsubh = 0.3
    #gmod.scalfracsubh = 'self'

    #gdat.minmmasshost = 1e10
    #gdat.maxmmasshost = 1e13
    #gmod.scalmasshost = 'self'
    #
    #gdat.minmmasssubh = 1e8
    #gdat.maxmmasssubh = 1e10
    #gmod.scalmasssubh = 'self'

    # collect groups of parameter indices into lists
    ## labels and scales for base parameters
    gmod.nameparagenrbase = []
    for name, k in gmod.indxpara.__dict__.items():
        if not np.isscalar(k):
            print('name')
            print(name)
            print('temp: no nonscalar should be here!')
            continue
        gmod.nameparagenrbase.append(name)
    gmod.numbparagenrbase = len(gmod.nameparagenrbase)
    gmod.indxparagenrbase = np.arange(gmod.numbparagenrbase)
    gmod.indxparagenrbasestdv = gmod.indxparagenrbase[gmod.numbpopl:]
    ## list of scalar variable names
    gmod.namepara.scal = list(gmod.nameparagenrbase) 
    gmod.namepara.scal += ['lliktotl']

    # derived parameters
    print('Determining the list of derived, fixed-dimensional parameter names...')
    gmod.namepara.genrelemextd = [[[] for g in gmod.indxparagenrelemsing[l]] for l in gmod.indxpopl]
    gmod.namepara.derielemextd = [[[] for k in gmod.indxparaderielemsing[l]] for l in gmod.indxpopl]
    gmod.namepara.genrelemflat = []
    gmod.namepara.derielemflat = []
    gmod.namepara.genrelemextdflat = []
    gmod.namepara.derielemextdflat = []
    for l in gmod.indxpopl:
        for g in gmod.indxparagenrelemsing[l]:
            gmod.namepara.genrelemflat.append(gmod.namepara.genrelem[l][g] + 'pop%d' % l)
            for d in range(gmod.maxmpara.numbelem[l]):
                gmod.namepara.genrelemextd[l][g].append(gmod.namepara.genrelem[l][g] + 'pop%d' % l + '%04d' % d)
                gmod.namepara.genrelemextdflat.append(gmod.namepara.genrelemextd[l][g][d])
        for k in gmod.indxparaderielemsing[l]:  
            gmod.namepara.derielemflat.append(gmod.namepara.derielem[l][k] + 'pop%d' % l)
            for d in range(gmod.maxmpara.numbelem[l]):
                gmod.namepara.derielemextd[l][k].append(gmod.namepara.derielem[l][k] + 'pop%d' % l + '%04d' % d)
                gmod.namepara.derielemextdflat.append(gmod.namepara.derielemextd[l][k][d])

    # list of element parameter names (derived and generative), counting label-degenerate element parameters only once 
    gmod.namepara.elem = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        gmod.namepara.elem[l].extend(gmod.namepara.genrelem[l])
        gmod.namepara.elem[l].extend(gmod.namepara.derielem[l])
    
    gmod.namepara.elemflat = []
    for l in gmod.indxpopl:
        gmod.namepara.elemflat.extend(gmod.namepara.elem[l])

    gmod.namepara.genrelemdefa = deepcopy(gmod.namepara.elemflat)
    if gmod.boolelemlghtanyy:
        for strgfeat in ['sind', 'curv', 'expc'] + ['sindcolr%04d' % i for i in gdat.indxenerinde]:
            if not strgfeat in gmod.namepara.genrelemdefa:
                gmod.namepara.genrelemdefa.append(strgfeat)

    # list of flattened generative element parameter names, counting label-degenerate element parameters only once
    gmod.namepara.genrelemkind = gmod.namepara.genrelemflat + gmod.namepara.derielemflat
    gmod.numbparagenrelemkind = len(gmod.namepara.genrelemkind)
    #gmod.inxparagenrscalelemkind = np.arange(gmod.numbparagenrelemkind)
    gmod.inxparagenrscalelemkind = tdpy.gdatstrt()
    
    gmod.numbparagenrelemextdflat = len(gmod.namepara.genrelemextdflat)
    gmod.indxparagenrelemextdflat = np.arange(gmod.numbparagenrelemextdflat)

    # list of parameter names (derived and generative), counting label-degenerate element parameters only once, element lists flattened
    gmod.namepara.kind = gmod.nameparagenrbase + gmod.listnameparaderitotl + gmod.namepara.genrelemflat + gmod.namepara.derielemflat
    
    gmod.numbparakind = len(gmod.namepara.kind)
    gmod.indxparakind = np.arange(gmod.numbparakind)

    # list of generative parameter names, separately including all label-degenerate element parameters, element lists flattened
    gmod.namepara.genrscalfull = gmod.nameparagenrbase + gmod.namepara.genrelemextdflat
    gmod.namepara.genrscalfull = np.array(gmod.namepara.genrscalfull)
    gmod.numbparagenrfull = len(gmod.namepara.genrscalfull)
    gmod.indxparagenrfull = np.arange(gmod.numbparagenrfull)

    # list of generative parameter names, counting label-degenerate element parameters only once, element lists flattened
    gmod.listnameparagenrscal = gmod.nameparagenrbase + gmod.namepara.genrelemflat
    gmod.numbparagenr = len(gmod.listnameparagenrscal)
    gmod.indxparagenr = np.arange(gmod.numbparagenr)

    # list of parameter names (derived and generative), element lists flattened
    gmod.listnameparatotl = gmod.nameparagenrbase + gmod.listnameparaderitotl + \
                                                                        gmod.namepara.genrelemextdflat + gmod.namepara.derielemextdflat
    
    gmod.nameparagenrbase = np.array(gmod.nameparagenrbase)

    for e in gmod.indxsersfgrd:
        gmod.namepara.scal += ['masshost' + strgsersfgrd + 'bein']
        for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
            gmod.namepara.scal += ['masshost' + strgsersfgrd + strgcalcmasssubh + 'bein']
    if gmod.numbparaelem > 0:
        if gmod.boollenssubh:
            for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                gmod.namepara.scal += ['masssubh' + strgcalcmasssubh + 'bein', 'fracsubh' + strgcalcmasssubh + 'bein'] 
    if gmod.numbparaelem > 0:
        gmod.namepara.scal += ['lpripena']
    if False and gmod.boolelemsbrtdfncanyy:
        for strgbins in ['lowr', 'higr']:
            gmod.namepara.scal += ['histcntp%sdfncen00evt0' % strgbins]
            gmod.namepara.scal += ['histcntp%sdfncsubten00evt0' % strgbins]
        for i in gdat.indxener:
            gmod.namepara.scal += ['fracsdenmeandarkdfncsubten%02d' % i]
        gmod.namepara.scal += ['booldfncsubt']
    if gmod.numbparaelem > 0:
        for q in gdat.indxrefr:
            if gdat.boolasscrefr[q]:
                for l in gmod.indxpopl:
                    gmod.namepara.scal += ['cmplpop%dpop%d' % (l, q)]
                    gmod.namepara.scal += ['fdispop%dpop%d' % (q, l)]
    
    gmod.numbvarbscal = len(gmod.namepara.scal)
    gmod.indxvarbscal = np.arange(gmod.numbvarbscal)
    
    # determine total label
    gmod.listnameparaglob = gmod.namepara.kind + gmod.namepara.genrelemextdflat + gmod.namepara.derielemextdflat
    gmod.listnameparaglob += ['cntpmodl']
    for l in gmod.indxpopl:
        for g in gmod.indxparagenrelemsing[l]:
            if not gmod.namepara.genrelem[l][g] in gmod.listnameparaglob:
                gmod.listnameparaglob.append(gmod.namepara.genrelem[l][g])
                gmod.listnameparaglob.append(gmod.namepara.derielem[l][g])

    for name in gmod.listnameparaglob:
        lablroot = getattr(gmod.lablrootpara, name)
        lablunit = getattr(gmod.lablunitpara, name)
        labltotl = tdpy.retr_labltotlsing(lablroot, lablunit)
        setattr(gmod.labltotlpara, name, labltotl)
    
    # define fact
    for l in gmod.indxpopl:
        for k in gmod.indxparakind:
            name = gmod.namepara.kind[k]
            scal = getattr(gmod.scalpara, name)
            if scal == 'self' or scal == 'logt':
                minm = getattr(gmod.minmpara, name)
                maxm = getattr(gmod.maxmpara, name)
                if scal == 'self':
                    fact = maxm - minm
                if scal == 'logt':
                    fact = np.log(maxm / minm)
                
                if fact == 0:
                    print('name')
                    print(name)
                    raise Exception('')
                setattr(gmod.factpara, name, fact)

    if gmod.numbparaelem > 0:
        gmod.indxparagenrfulleleminit = gmod.indxparagenrbase[-1] + 1
    else:
        gmod.indxparagenrfulleleminit = -1
    

    ## arrays of parameter features (e.g., minm, maxm, labl, scal, etc.)
    for featpara in gdat.listfeatparalist:
        
        gmodfeat = getattr(gmod, featpara + 'para')
        
        ### elements
        #for strgtypepara in gdat.liststrgtypepara:
        #    listname = getattr(gmod.namepara, strgtypepara + 'elem')
        #    listfeat = [[] for l in gmod.indxpopl]
        #    listfeatflat = []

        #    for l in gmod.indxpopl:
        #        
        #        numb = getattr(gmod, 'numbpara' + strgtypepara + 'elemsing')[l]
        #        listfeat[l] = [[] for k in range(numb)]
        #        for k in range(numb):
        #            scal = getattr(gmod.scalpara, listname[l][k])
        #            if featpara == 'fact' and not (scal == 'self' or scal == 'logt'):
        #                continue
        #            if featpara == 'mean' and (scal != 'gaus' and scal != 'lnor'):
        #                continue
        #            if featpara == 'stdv' and (scal != 'gaus' and scal != 'lnor'):
        #                continue
        #            
        #            if strgtypepara == 'genr':
        #                strgextn = 'pop%d' % l
        #            else:
        #                strgextn = ''
        #            print('featpara')
        #            print(featpara)
        #            print('listname')
        #            print(listname)
        #            listfeat[l][k] = getattr(gmodfeat, listname[l][k] + strgextn)
        #            listfeatflat.append(listfeat[l][k])
        #    setattr(gmodfeat, strgtypepara + 'elem', listfeat)
        #    setattr(gmodfeat, strgtypepara + 'elemflat', listfeatflat)
        
        ### groups of parameters inside the parameter vector
        ### 'base': all fixed-dimensional generative parameters
        ### 'full': all generative parameters
        for strggroppara in ['base', 'full']:
            indx = getattr(gmod, 'indxparagenr' + strggroppara)
            feat = [0. for k in indx]
        
            for attr, valu in gmod.indxpara.__dict__.items():
                
                if not np.isscalar(valu):
                    continue

                scal = getattr(gmod.scalpara, attr)
                if not (scal == 'self' or scal == 'logt') and featpara == 'fact':
                    continue

                if scal != 'gaus' and (featpara == 'mean' or featpara == 'stdv'):
                    print('Mean or Std for non-Gaussian')
                    continue
                
                if featpara == 'name':
                    feat[valu] = attr
                else:
                    feat[valu] = getattr(gmodfeat, attr)
            
            feat = np.array(feat)
            setattr(gmodfeat, 'genr' + strggroppara, feat)
    
    
    #print('gmod.minmpara')
    #for attr, varb in gmod.minmpara.__dict__.items():
    #    print(attr, varb)
    #print('gmod.maxmpara')
    #for attr, varb in gmod.maxmpara.__dict__.items():
    #    print(attr, varb)
    #print('gmod.scalpara')
    #for attr, varb in gmod.scalpara.__dict__.items():
    #    print(attr, varb)
    #raise Exception('')

    ## population groups
    ### number of elements
    for strgvarb in ['numbelem', 'meanelem']:
        listindxpara = []
        if strgmodl == 'true':
            listpara = []
        for strg, valu in gmod.indxpara.__dict__.items():
            if strg.startswith(strgvarb + 'p'):
                listindxpara.append(valu)
                if strgmodl == 'true':
                    listpara.append(getattr(gmod.this, strg))
        listindxpara = np.array(listindxpara)
        setattr(gmod.indxpara, strgvarb, listindxpara)
        if strgmodl == 'true':
            listpara = np.array(listpara)
            setattr(gmod, strgvarb, listpara)
        
    ### parameters of priors for element parameters
    gmod.indxpara.prioelem = []
    for strg, valu in gmod.indxpara.__dict__.items():
        if strg == 'dist' and np.isscalar(valu):
            gmod.indxpara.prioelem.append(valu)
    gmod.indxpara.prioelem = np.array(gmod.indxpara.prioelem) 
    
    ### hyperparameters
    if gmod.typemodltran == 'pois':
        gmod.indxpara.hypr = np.array(list(gmod.indxpara.prioelem) + list(gmod.indxpara.meanelem))
    else:
        gmod.indxpara.hypr = gmod.indxpara.prioelem
        
    ## generative base parameter indices for each scaling
    gmod.listindxparagenrbasescal = dict()
    for scaltype in gdat.listscaltype:
        gmod.listindxparagenrbasescal[scaltype] = np.where(np.array(gmod.scalpara.genrbase) == scaltype)[0]

    if gdat.booldiagmode:
        if np.where(gmod.scalpara.genrfull == 0)[0].size > 0:
            raise Exception('')


def plot_lens(gdat):
    
    if gmod.boolelemdeflsubh:
        xdat = gdat.binspara.angl[1:] * gdat.anglfact
        lablxdat = gdat.labltotlpara.gang
        
        listdeflscal = np.array([4e-2, 4e-2, 4e-2]) / gdat.anglfact
        listanglscal = np.array([0.05, 0.1, 0.05]) / gdat.anglfact
        listanglcutf = np.array([1.,    1.,  10.]) / gdat.anglfact
        listasym = [False, False, False]
        listydat = []
        for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
            listydat.append(retr_deflcutf(gdat.binspara.angl[1:], deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
        
        for scalxdat in ['self', 'logt']:
            path = gdat.pathinitintr + 'deflcutf' + scalxdat + '.pdf'
            tdpy.plot_gene(path, xdat, listydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, \
                                                             lablydat=r'$\alpha_n$ [$^{\prime\prime}$]', limtydat=[1e-3, 1.5e-2], limtxdat=[None, 2.])
       
        # pixel-convoltuion of the Sersic profile
        # temp -- y axis labels are wrong, should be per solid angle
        xdat = gdat.binspara.lgalsers * gdat.anglfact
        for n in range(gdat.numbindxsers + 1):
            for k in range(gdat.numbhalfsers + 1):
                if k != 5:
                    continue
                path = gdat.pathinitintr + 'sersprofconv%04d%04d.pdf' % (n, k)
                tdpy.plot_gene(path, xdat, gdat.sersprof[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                #path = gdat.pathinitintr + 'sersprofcntr%04d%04d.pdf' % (n, k)
                #tdpy.plot_gene(path, xdat, gdat.sersprofcntr[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.pdf' % (n, k)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.pdf' % (n, k)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], scalxdat='logt', \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
       
        xdat = gdat.binspara.angl * gdat.anglfact
        listspec = np.array([1e-19, 1e-18, 1e-18, 1e-18]) / gdat.anglfact
        listsize = np.array([0.3, 1., 1., 1.]) / gdat.anglfact
        listindx = np.array([4., 2., 4., 10.])
        listydat = []
        listlabl = []
        for spec, size, indx in zip(listspec, listsize, listindx):
            listydat.append(spec * retr_sbrtsersnorm(gdat.binspara.angl, size, indxsers=indx))
            listlabl.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
        path = gdat.pathinitintr + 'sersprof.pdf'
        tdpy.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, \
                                                                                                   listlegd=listlegd, listhlin=1e-7, limtydat=[1e-8, 1e0])
    
        minmredshost = 0.01
        maxmredshost = 0.4
        minmredssour = 0.01
        maxmredssour = 2.
        numbreds = 200
        retr_axis(gdat, 'redshost')
        retr_axis(gdat, 'redssour')
        
        gdat.meanpara.adishost = np.empty(numbreds)
        for k in range(numbreds):
            gdat.meanpara.adishost[k] = gdat.adisobjt(gdat.meanpara.redshost[k])
        
        asca = 0.1 / gdat.anglfact
        acut = 1. / gdat.anglfact
    
        minmmass = np.zeros((numbreds + 1, numbreds + 1))
        maxmmass = np.zeros((numbreds + 1, numbreds + 1))
        for k, redshost in enumerate(gdat.binspara.redshost):
            for n, redssour in enumerate(gdat.binspara.redssour):
                if redssour > redshost:
                    adishost = gdat.adisobjt(redshost)
                    adissour = gdat.adisobjt(redssour)
                    adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
                    factmcutfromdefs = retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut)
                    minmmass[n, k] = np.log10(factmcutfromdefs * gdat.minmdefs)
                    maxmmass[n, k] = np.log10(factmcutfromdefs * gdat.maxmdefs)
       
        #valulevl = np.linspace(7.5, 9., 5)
        valulevl = [7.0, 7.3, 7.7, 8., 8.6]
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        cont = axis.contour(gdat.binspara.redshost, gdat.binspara.redssour, minmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=20, fmt='%.3g')
        axis.set_xlabel(r'$z_{\rm{hst}}$')
        axis.set_ylabel(r'$z_{\rm{src}}$')
        axis.set_title(r'$M_{c,min}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsminm.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        valulevl = np.linspace(9., 11., 20)
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=9., vmax=11.)
        cont = axis.contour(gdat.binspara.redshost, gdat.binspara.redssour, maxmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=15, fmt='%.3g')
        axis.set_xlabel('$z_{hst}$')
        axis.set_ylabel('$z_{src}$')
        axis.set_title(r'$M_{c,max}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsmaxm.pdf'
        plt.colorbar(imag) 
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adishost * gdat.sizepixl * 1e-3)
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adishost * 2. * gdat.maxmgangdata * 1e-3)
        axis.set_xlabel('$z_h$')
        axis.set_yscale('log')
        axis.set_ylabel(r'$\lambda$ [kpc]')
        path = gdat.pathinitintr + 'wlenreds.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        fracacutasca = np.logspace(-1., 2., 20)
        mcut = retr_mcutfrommscl(fracacutasca)
        axis.lognp.log(fracacutasca, mcut)
        axis.set_xlabel(r'$\tau_n$')
        axis.set_ylabel(r'$M_{c,n} / M_{0,n}$')
        axis.axhline(1., ls='--')
        path = gdat.pathinitintr + 'mcut.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
       

def retr_listrtagprev(strgcnfg, pathpcat):
    
    # list of PCAT run plot outputs
    pathimag = pathpcat + '/imag/'
    listrtag = fnmatch.filter(os.listdir(pathimag), '2*')
    
    listrtagprev = []
    for rtag in listrtag:
        strgstat = pathpcat + '/data/outp/' + rtag
        
        if chec_statfile(pathpcat, rtag, 'gdatmodipost', typeverb=0) and strgcnfg + '_' + rtag[16:].split('_')[-1] == rtag[16:]:
            listrtagprev.append(rtag) 
    
    listrtagprev.sort()

    return listrtagprev


def make_legd(axis, offs=None, loca=1, numbcols=1, ptch=None, line=None):
   
    hand, labl = axis.get_legend_handles_labels()
    legd = axis.legend(hand, labl, fancybox=True, frameon=True, bbox_to_anchor=offs, bbox_transform=axis.transAxes, ncol=numbcols, loc=loca, labelspacing=1, handlelength=2)
    legd.get_frame().set_fill(True)
    legd.get_frame().set_facecolor('white')


def setp_namevarbsing(gdat, gmod, strgmodl, strgvarb, popl, ener, evtt, back, isfr, iele):
    
    if popl == 'full':
        indxpopltemp = gmod.indxpopl
    elif popl != 'none':
        indxpopltemp = [popl]
    
    if ener == 'full':
        indxenertemp = gdat.indxener
    elif ener != 'none':
        indxenertemp = [ener]
    
    if evtt == 'full':
        indxevtttemp = gdat.indxevtt
    elif evtt != 'none':
        indxevtttemp = [evtt]
    
    if back == 'full':
        gmod.indxbacktemp = gmod.indxback
    elif isinstance(back, int):
        gmod.indxbacktemp = np.array([back])
    
    liststrgvarb = []
    if iele != 'none':
        for l in gmod.indxpopl:
            if iele == 'full':
                listiele = np.arange(gmod.maxmpara.numbelem)
            else:
                listiele = [iele]
            for k in listiele:
                liststrgvarb.append(strgvarb + 'pop%d%04d' % (l, k))
    
    if popl != 'none' and ener == 'none' and evtt == 'none' and back == 'none' and iele == 'none':
        for l in indxpopltemp:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none' and isfr != 'none':
        for e in indxisfrtemp:
            liststrgvarb.append(strgvarb + 'isf%d' % e)
    
    if popl == 'none' and ener != 'none' and evtt != 'none' and back == 'none':
        for i in indxenertemp:
            for m in indxevtttemp:
                liststrgvarb.append(strgvarb + 'en%02devt%d' % (i, m))
    
    if popl == 'none' and ener != 'none' and evtt == 'none' and back != 'none':
        for c in gmod.indxbacktemp:
            for i in indxenertemp:
                liststrgvarb.append(strgvarb + 'back%04den%02d' % (c, i))
    
    if popl == 'none' and ener == 'none' and evtt == 'none' and back != 'none':
        for c in gmod.indxbacktemp:
            liststrgvarb.append(strgvarb + 'back%04d' % c)
    
    if popl == 'none' and ener != 'none' and evtt == 'none' and back == 'none':
        for i in indxenertemp:
            liststrgvarb.append(strgvarb + 'en%02d' % i)
    
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none' and isfr == 'none':
        liststrgvarb.append(strgvarb)
    
    if gdat.booldiagmode:
        for strgvarb in liststrgvarb:
            if liststrgvarb.count(strgvarb) != 1:
                print('liststrgvarb')
                print(liststrgvarb)
                print('popl')
                print(popl)
                print('ener')
                print(ener)
                print('evtt')
                print(evtt)
                print('back')
                print(back)
                print('isfr')
                print(isfr)
                print('iele')
                print(iele)
                raise Exception('')
    
    return liststrgvarb


def setp_varb(gdat, strgvarbbase, valu=None, minm=None, maxm=None, scal='self', lablroot=None, lablunit='', mean=None, stdv=None, cmap=None, numbbins=10, \
                        popl='none', ener='none', evtt='none', back='none', isfr='none', iele='none', \
                        boolinvr=False, \
                        strgmodl=None, strgstat=None, \
                        ):
    '''
    Set up variable values across all models (true and fitting) as well as all populations, energy bins, 
    event bins, background components, and Sersic components 
    '''
    
    # determine the list of models
    if strgmodl is None:
        if gdat.typedata == 'mock':
            liststrgmodl = ['true', 'fitt', 'plot']
        else:
            liststrgmodl = ['fitt', 'plot']
    else:
        if strgmodl == 'true' or strgmodl == 'plot' or strgmodl == 'refr':
            liststrgmodl = [strgmodl]
        else:
            liststrgmodl = ['fitt', 'plot']
    print('liststrgmodl')
    print(liststrgmodl)
    for strgmodl in liststrgmodl:
        
        if strgmodl == 'plot':
            gmod = gdat.fitt
            gmodoutp = gdat
        else:
            gmod = getattr(gdat, strgmodl)
            gmodoutp = gmod

        # get the list of names of the variable
        liststrgvarbnone = setp_namevarbsing(gdat, gmod, strgmodl, strgvarbbase, popl, ener, evtt, back, isfr, 'none')
        
        if iele != 'none':
            liststrgvarb = setp_namevarbsing(gdat, gmod, strgmodl, strgvarbbase, popl, ener, evtt, back, isfr, iele)
        else:
            liststrgvarb = liststrgvarbnone

        # set the values of each variable in the list
        for strgvarb in liststrgvarb:
            if minm is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.minmpara, strgvarb, minm)
            
            if maxm is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.maxmpara, strgvarb, maxm)
            
            if mean is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.meanpara, strgvarb, mean)
            
            if stdv is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.meanpara, strgvarb, stdv)
            
            if valu is not None:
                if strgstat is None:
                    print('strgvarb')
                    print(strgvarb)
                    print('strgmodl')
                    print(strgmodl)
                    print('valu')
                    print(valu)
                    print('')
                    setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, valu)
                elif strgstat == 'this':
                    setp_varbcore(gdat, strgmodl, gmodoutp.this, strgvarb, valu)
                    
            if scal is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.scalpara, strgvarb, scal)

            if lablroot is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.lablrootpara, strgvarb, lablroot)

            if lablunit is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.lablunitpara, strgvarb, lablunit)

            if cmap is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp.cmappara, strgvarb, cmap)

            setp_varbcore(gdat, strgmodl, gmodoutp.numbbinspara, strgvarb, numbbins)
            
            # create limt, bins, mean, and delt
            if minm is not None and maxm is not None or mean is not None and stdv is not None:
                # determine minima and maxima for Gaussian or log-Gaussian distributed parameters
                if mean is not None:
                    minm = mean - gdat.numbstdvgaus * stdv
                    maxm = mean + gdat.numbstdvgaus * stdv
                
                # uniformly-distributed
                if scal == 'self' or scal == 'pois' or scal == 'gaus':
                    binsunif = np.linspace(minm, maxm, numbbins + 1)
                if scal == 'logt' or scal == 'powr':
                    binsunif = np.linspace(np.log10(minm), np.log10(maxm), numbbins + 1)
                    if gdat.booldiagmode:
                        if minm <= 0.:
                            raise Exception('')
                if scal == 'asnh':
                    binsunif = np.linspace(np.arcsinh(minm), np.arcsinh(maxm), numbbins + 1)
                
                if boolinvr:
                    binsunif = binsunif[::-1]
                
                meanparaunif = (binsunif[1:] + binsunif[:-1]) / 2.
                
                if scal == 'self' or scal == 'pois' or scal == 'gaus':
                    meanpara = meanparaunif
                    bins = binsunif
                    minmunif = minm
                    maxmunif = maxm
                if scal == 'logt' or scal == 'powr':
                    meanpara = 10**meanparaunif
                    bins = 10**binsunif
                    minmunif = np.log10(minm)
                    maxmunif = np.log10(maxm)
                if scal == 'asnh':
                    meanpara = np.sinh(meanparaunif)
                    bins = np.sinh(binsunif)
                    minmunif = np.arcsinh(minm)
                    maxmunif = np.arcsinh(maxm)

                delt = np.diff(bins) 
                limt = np.array([minm, maxm]) 
                
                # 'self' is not yet defined
                if scal == 'asnh' or scal == 'logt' or scal == 'powr':
                    listvalutickmajr, listlabltickmajr, listvalutickminr, listlabltickminr = tdpy.retr_valulabltick(minm, maxm, scal)
                    setattr(gmodoutp.labltickmajrpara, strgvarb, listlabltickmajr)
                    setattr(gmodoutp.valutickmajrpara, strgvarb, listvalutickmajr)
                    setattr(gmodoutp.labltickminrpara, strgvarb, listlabltickminr)
                    setattr(gmodoutp.valutickminrpara, strgvarb, listvalutickminr)
                
                #labltick = np.empty(gdat.numbtickcbar, dtype=object)
                #for k in range(gdat.numbtickcbar):
                #    if scal == 'asnh':
                #        valutick[k] = np.sinh(tickunif[k])
                #    if scal == 'logt' or scal == 'powr':
                #        valutick[k] = 10**(tickunif[k])

                #    # avoid very small, but nonzero central values in the residual count color maps
                #    if strgcbar == 'cntpresi' and np.fabs(valutick[k]) < 1e-5:
                #        valutick[k] = 0.

                #    if strgcbar == 'cntpdata' and np.amax(valutick) > 1e3:
                #        labltick[k] = '%d' % valutick[k]
                #    else:
                #        labltick[k] = '%.3g' % valutick[k]

                setattr(gmodoutp.limtpara, strgvarb, limt)
                setattr(gmodoutp.binspara, strgvarb, bins)
                setattr(gmodoutp.meanpara, strgvarb, meanpara)
                setattr(gmodoutp.deltpara, strgvarb, delt)
       

def retr_ticklabltemp(gdat, strgcbar):
    
    minm = getattr(gdat.minmpara, strgcbar)
    maxm = getattr(gdat.maxmpara, strgcbar)
    scal = getattr(gdat.scalpara, strgcbar)
    numb = gdat.numbtickcbar - 1
    retr_axis(gdat, strgcbar, numb=numb)

    minmscal = minm
    if scal == 'asnh':
        minmscal = np.arcsinh(minmscal)
    if scal == 'logt':
        minmscal = np.log10(minmscal)
    maxmscal = maxm
    if scal == 'asnh':
        maxmscal = np.arcsinh(maxmscal)
    if scal == 'logt':
        maxmscal = np.log10(maxmscal)

    tickscal = np.linspace(minmscal, maxmscal, gdat.numbtickcbar)
    labl = np.empty(gdat.numbtickcbar, dtype=object)
    tick = np.copy(tickscal)
    for k in range(gdat.numbtickcbar):
        if scal == 'asnh':
            tick[k] = np.sinh(tickscal[k])
        elif scal == 'logt':
            tick[k] = 10**(tickscal[k])

        # avoid very small, but nonzero central values in the residual count color maps
        if strgcbar == 'cntpresi' and np.fabs(tick[k]) < 1e-5:
            tick[k] = 0.

        if strgcbar == 'cntpdata' and np.amax(tick) > 1e3:
            labl[k] = '%d' % tick[k]
        else:
            labl[k] = '%.3g' % tick[k]
    setattr(gdat.tickpara, strgcbar, tick)


def retr_axistemp(gdat, strgvarb, strgmodl=None, boolinvr=False):
    
    if strgmodl is None:
        listgdattemp = [gdat]
        for strgmodl in gdat.liststrgmodl:
            listgdattemp.append(getattr(gdat, strgmodl))
    elif strgmodl == 'fitt' or strgmodl == 'true':
        listgdattemp = [getattr(gdat, strgmodl)]
    elif strgmodl == 'allm':
        listgdattemp = []
        for strgmodl in gdat.liststrgmodl:
            listgdattemp = getattr(gdat, strgmodl)
    
    for gdattemp in listgdattemp:
        minm = getattr(gdattemp.minmpara, strgvarb)
        maxm = getattr(gdattemp.maxmpara, strgvarb)
        numb = getattr(gdattemp.numbbinspara, strgvarb)
        scal = getattr(gdattemp.scalpara, strgvarb)

        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            binsscal = np.linspace(minm, maxm, numb + 1)
        if scal == 'logt':
            print('minm')
            print(minm)
            print('maxm')
            print(maxm)
            print('strgvarb')
            print(strgvarb)
            binsscal = np.linspace(np.log10(minm), np.log10(maxm), numb + 1)
            print('')
            if gdat.booldiagmode:
                if minm <= 0.:
                    raise Exception('')

        if scal == 'asnh':
            binsscal = np.linspace(np.arcsinh(minm), np.arcsinh(maxm), numb + 1)
        
        if boolinvr:
            binsscal = binsscal[::-1]
        
        meanvarbscal = (binsscal[1:] + binsscal[:-1]) / 2.
        
        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            meanvarb = meanvarbscal
            bins = binsscal
        if scal == 'logt':
            meanvarb = 10**meanvarbscal
            bins = 10**binsscal
        if scal == 'asnh':
            meanvarb = np.sinh(meanvarbscal)
            bins = np.sinh(binsscal)

        delt = np.diff(bins) 
        limt = np.array([np.amin(bins), np.amax(bins)]) 
        
        setattr(gdattemp.limtpara, strgvarb, limt)
        setattr(gdattemp.binspara, strgvarb, bins)
        setattr(gdattemp.meanpara, strgvarb, meanvarb)
        setattr(gdattemp.deltpara, strgvarb, delt)


def setp_varbcore(gdat, strgmodl, gdattemp, strgvarbtemp, valu):
    
    # check if the variable is defined by the user
    try:
        valutemp = getattr(gdattemp, strgvarbtemp)
        if valutemp is None:
            raise
        if gdat.typeverb > 0:
            print('Received custom value for %s, %s: %s' % (strgvarbtemp, strgmodl, valutemp))
    # if not defined or defined as None, define it
    except:
        
        setattr(gdattemp, strgvarbtemp, valu)


def intp_sinc(gdat, lgal, bgal):

    intpsinc = 4. * gdat.numbsidepsfn**2 * np.sum(gdat.temppsfn * sinc(gdat.numbsidepsfn * (gdat.gridpsfnlgal + lgal) - gdat.gridpsfnlgal) * \
                                                                                                      sinc(gdat.numbsidepsfn * (gdat.gridpsfnbgal + bgal) - gdat.gridpsfnbgal))

    return intpsinc


def retr_fluxbrgt(gdat, lgal, bgal, flux):

    if lgal.size == 0:
        fluxbrgt = np.array([0.])
        fluxbrgtassc = np.array([0.])
    else:
        indxbrgt = np.argmax(flux)
        fluxbrgt = flux[indxbrgt]

    return fluxbrgt, fluxbrgtassc


def init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot):
    
    figrsize = (gdat.sizeimag, gdat.sizeimag)
    figr, axis = plt.subplots(figsize=figrsize)
    
    nameplot = strgplot

    if gdat.numbener > 1:
        nameplot += 'en%02d' % gdat.indxenerincl[indxenerplot]
    
    if gdat.numbener > 1:
        if indxevttplot == -1:
            nameplot += 'evtA'
        else:
            nameplot += 'evt%d' % gdat.indxevttincl[indxevttplot]
    
    if gdat.fitt.numbpopl > 1:
        if indxpoplplot == -1:
            nameplot += 'popA'
        else:
            nameplot += 'pop%d' % indxpoplplot

    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, nameplot)
    
    print('gdat.fitt.labltotlpara.lgalpop0')
    print(gdat.fitt.labltotlpara.lgalpop0)
    print('gdat.fitt.labltotlpara.bgalpop0')
    print(gdat.fitt.labltotlpara.bgalpop0)
    axis.set_xlabel(gdat.fitt.labltotlpara.lgalpop0)
    axis.set_ylabel(gdat.fitt.labltotlpara.bgalpop0)
    titl = ''
    if indxenerplot is not None and gdat.numbener > 1 and strgplot.endswith('cnts'):
        titl = gdat.strgener[indxenerplot]
    if indxevttplot is not None and gdat.numbevtt > 1 and strgplot.endswith('cnts'):
        titl += ' ' + gdat.strgevtt[indxevttplot]
    axis.set_title(titl)

    return figr, axis, path


def draw_frambndr(gdat, axis):
    
    outr = max(gdat.frambndrmodl, gdat.frambndrdata)
    axis.set_xlim([-outr, outr])
    axis.set_ylim([-outr, outr])
    innr = min(gdat.frambndrmodl, gdat.frambndrdata)
    axis.axvline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axvline(-innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(-innr, ls='--', alpha=gdat.alphbndr, color='black')


def retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot=None, indxevttplot=-1, booltdim=False, imag=None):
    
    draw_frambndr(gdat, axis)
    
    # take the relevant energy and PSF bins
    if indxenerplot is not None:
        if indxevttplot == -1:
            maps = np.sum(maps[indxenerplot, ...], axis=1)
        else:
            maps = maps[indxenerplot, :, indxevttplot]
    
    # project the map to 2D
    if gdat.typepixl == 'heal':
        maps = tdpy.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
    
    if gdat.typepixl == 'cart':
        shap = [gdat.numbsidecart] + list(maps.shape)
        shap[1] = gdat.numbsidecart
        shapflat = list(maps.shape)
        shapflat[0] = gdat.numbpixlfull
        mapstemp = np.zeros(shapflat)
        if maps.size == gdat.indxpixlrofi.size:
            mapstemp[gdat.indxpixlrofi, ...] = maps
        else:
            mapstemp[:, ...] = maps
        maps = mapstemp.reshape(shap).swapaxes(0, 1)

    # temp -- this is needed to bring the Fermi-LAT map to the right direction
    #maps = fliplr(maps)

    # rescale the map
    if strgmodl is not None:
        gmod = getattr(gdat, strgmodl)
    else:
        gmod = gdat

    scal = getattr(gdat.scalpara, strgcbar)
    cmap = getattr(gdat.cmappara, strgcbar)
    vmin = getattr(gdat.minmpara, strgcbar)
    vmax = getattr(gdat.maxmpara, strgcbar)
    if scal == 'asnh':
        maps = np.arcsinh(maps)
    if scal == 'logt':
        maps = np.log10(maps)
    if imag is None:
        imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', vmin=vmin, vmax=vmax, alpha=gdat.alphmaps)
        return imag
    else:
        imag.set_data(maps)
    

def make_cbar(gdat, axis, imag, strgvarb):

    # make a color bar
    valutickmajr = getattr(gdat.valutickmajrpara, strgvarb)
    labltickmajr = getattr(gdat.labltickmajrpara, strgvarb)
    
    print('valutickmajr')
    print(valutickmajr)
    print('labltickmajr')
    print(labltickmajr)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05, aspect=15)
    cbar.set_ticks(valutickmajr)
    cbar.set_ticklabels(labltickmajr)
    
    return cbar


def make_legdmaps(gdat, strgstat, strgmodl, axis, mosa=False, assc=False):
    
    gmod = getattr(gdat, strgmodl)
    
    # transdimensional elements
    if strgmodl == 'fitt' and (strgstat == 'pdfn' and gdat.boolcondcatl or strgstat == 'this') and gmod.numbparaelem > 0:
        for l in gmod.indxpopl:
            colr = retr_colr(gdat, strgstat, strgmodl, l)
            if strgstat == 'pdfn':
                labl = 'Condensed %s %s' % (gmod.legd, gmod.legdpopl[l])
            else:
                labl = 'Sample %s %s' % (gmod.legd, gmod.legdpopl[l])
            if not gmod.maxmpara.numbelem[l] == 0:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                        label=labl, marker=gmod.listelemmrkr[l], lw=gdat.mrkrlinewdth, color=colr)
    
    for q in gdat.indxrefr:
        if not np.amax(gdat.refr.numbelem[q]) == 0:
            if assc:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                       label=gdat.refr.lablhits[q], marker=gdat.refr.listmrkrhits[q], lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                       label=gdat.refr.lablmiss[q], marker=gdat.refr.listmrkrmiss[q], lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
            else:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                       label=gdat.refr.lablelem[q], marker=gdat.refr.listmrkrmiss[q], lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
    
    # fixed-dimensional objects
    if strgmodl == 'fitt':
        if gmod.boollens:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                           label='%s Source' % gmod.lablmodl, marker='<', lw=gdat.mrkrlinewdth, color=gmod.colr)
        
        if gmod.typeemishost != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                           label='%s Host' % gmod.lablmodl, marker='s', lw=gdat.mrkrlinewdth, color=gmod.colr)
    
    if gdat.typedata == 'mock':
        if gmod.boollens:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Source' % gdat.refr.labl, marker='>', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
        
        if gmod.typeemishost != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Host' % gdat.refr.labl, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
    
    temphand, temp = axis.get_legend_handles_labels()
    numblabl = len(temp)
    
    if numblabl == 4:
        numbcols = 2
    else:
        numbcols = 3
    if mosa:
        axis.legend(bbox_to_anchor=[1., 1.15], loc='center', ncol=numbcols)
    else:
        axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=numbcols)
        

def supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot=-1, assc=False):
    
    gmod = getattr(gdat, strgmodl)
    gmodstat = getattr(gmod, strgstat)
    
    # associations with the reference elements
    for q in gdat.indxrefr:
        if gdat.refr.numbelem[q] > 0:
            if indxpoplplot == -1:
                listindxpoplplot = gmod.indxpopl
            else:
                listindxpoplplot = [indxpoplplot]
            for l in listindxpoplplot:
                reframpl = gdat.refr.dictelem[q][gdat.refr.nameparagenrelemampl[q]][0, :]
                mrkrsize = retr_mrkrsize(gdat, strgmodl, reframpl, gdat.refr.nameparagenrelemampl[q])
                lgal = np.copy(gdat.refr.dictelem[q]['lgal'][0, :])
                bgal = np.copy(gdat.refr.dictelem[q]['bgal'][0, :])
                numbelem = int(gdat.refr.numbelem[q])
                
                if gdatmodi is not None and gmod.numbparaelem > 0 and assc:   
                    ### hit
                    indx = gdatmodi.this.indxelemrefrasschits[q][l]
                    if indx.size > 0:
                        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, label=gdat.refr.lablhits, \
                                                                      marker=gdat.refrlistmrkrhits[q], lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
                    ### missed
                    indx = gdatmodi.this.indxelemrefrasscmiss[q][l]
                else:
                    indx = np.arange(lgal.size)
                
                if indx.size > 0: 
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, facecolor='none', \
                                                             label=gdat.refr.listlablmiss, marker=gdat.refr.listmrkrmiss[q], \
                                                             lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
        
            sizexoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            sizeyoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            if 'etag' in gdat.refr.namepara.elem[q]:
                for k in range(indx.size):
                    axis.text(gdat.anglfact * lgal[indx[k]] + sizexoff, gdat.anglfact * bgal[indx[k]] + sizeyoff, gdat.refretag[q][indx[k]], \
                                                                                            verticalalignment='center', horizontalalignment='center', \
                                                                                                             color='red', fontsize=1)

    # temp -- generalize this to input refrlgalhost vs.
    if gdat.typedata == 'mock':
        ## host galaxy position
        if gmod.typeemishost != 'none':
            for e in gmod.indxsersfgrd:
                lgalhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'lgalhostisf%d' % (e))]
                bgalhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'bgalhostisf%d' % (e))]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, facecolor='none', alpha=0.7, \
                                             label='%s Host %d' % (gdat.refr.labl, e), s=300, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
        if gmod.boollens:
            ## host galaxy Einstein radius
            for e in gmod.indxsersfgrd:
                truelgalhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'lgalhostisf%d' % (e))]
                truebgalhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'bgalhostisf%d' % (e))]
                truebeinhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % (e))]
                axis.add_patch(plt.Circle((gdat.anglfact * truelgalhost, \
                                           gdat.anglfact * truebgalhost), \
                                           gdat.anglfact * truebeinhost, \
                                           edgecolor=gdat.refr.colr, facecolor='none', lw=gdat.mrkrlinewdth))
            
        if gmod.boollens:
            ## source galaxy position
            axis.scatter(gdat.anglfact * gmodstat.paragenrscalfull[gmod.indxpara.lgalsour], \
                                                    gdat.anglfact * gmodstat.paragenrscalfull[gmod.indxpara.bgalsour], \
                                                    facecolor='none', \
                                                    alpha=0.7, \
                                                    #alpha=gdat.alphelem, \
                                                    label='%s Source' % gdat.refr.labl, s=300, marker='>', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
        
    # model catalog
    if indxpoplplot == -1:
        listindxpoplplot = gmod.indxpopl
    else:
        listindxpoplplot = [indxpoplplot]
    for l in listindxpoplplot:
        if gdatmodi is not None:
            if gmod.numbparaelem > 0:
                colr = retr_colr(gdat, strgstat, strgmodl, l)
                mrkrsize = retr_mrkrsize(gdat, strgmodl, gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[gmod.nameparagenrelemampl[l]][l]], gmod.nameparagenrelemampl[l])
                if 'lgal' in gdatmodi.this.indxparagenrfullelem:
                    lgal = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l]['lgal']]
                    bgal = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l]['bgal']]
                else:
                    gang = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l]['gang']]
                    aang = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l]['aang']]
                    lgal, bgal = retr_lgalbgal(gang, aang)
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphelem, label='Sample', marker=gmod.listelemmrkr[l], \
                                                                                                                           lw=gdat.mrkrlinewdth, color=colr)

            ## source
            if gmod.boollens:
                lgalsour = gdatmodi.this.paragenrscalfull[gmod.indxpara.lgalsour]
                bgalsour = gdatmodi.this.paragenrscalfull[gmod.indxpara.bgalsour]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, facecolor='none', \
                                                  alpha=gdat.alphelem, \
                                                  label='%s Source' % gmod.lablpara, s=300, marker='<', lw=gdat.mrkrlinewdth, color=gmod.colr)
    
            if gmod.typeemishost != 'none':
                ## host
                lgalhost = [[] for e in gmod.indxsersfgrd]
                bgalhost = [[] for e in gmod.indxsersfgrd]
                for e in gmod.indxsersfgrd:
                    lgalhost[e] = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'lgalhostisf%d' % (e))]
                    bgalhost[e] = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'bgalhostisf%d' % (e))]
                    axis.scatter(gdat.anglfact * lgalhost[e], gdat.anglfact * bgalhost[e], facecolor='none', \
                                                     alpha=gdat.alphelem, \
                                                     label='%s Host' % gmod.lablpara, s=300, marker='s', lw=gdat.mrkrlinewdth, color=gmod.colr)
                    if gmod.boollens:
                        beinhost = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % (e))]
                        axis.add_patch(plt.Circle((gdat.anglfact * lgalhost[e], gdat.anglfact * bgalhost[e]), \
                                                       gdat.anglfact * beinhost, edgecolor=gmod.colr, facecolor='none', \
                                                       lw=gdat.mrkrlinewdth, ls='--'))
                
    # temp
    if strgstat == 'pdfn' and gdat.boolcondcatl and gmod.numbparaelem > 0:
        lgal = np.zeros(gdat.numbprvlhigh)
        bgal = np.zeros(gdat.numbprvlhigh)
        ampl = np.zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                lgal[cntr] = gdat.dictglob['poststkscond'][r]['lgal'][0]
                bgal[cntr] = gdat.dictglob['poststkscond'][r]['bgal'][0]
                # temp -- this does not allow sources with different spectra to be assigned to the same stacked sample
                ampl[cntr] = gdat.dictglob['poststkscond'][r][gmod.nameparagenrelemampl[l]][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, strgmodl, ampl, gmod.nameparagenrelemampl[l])
        
        colr = retr_colr(gdat, strgstat, strgmodl, l)
        axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, \
                                    label='Condensed', marker=gmod.listelemmrkr[l], color='black', lw=gdat.mrkrlinewdth)
        for r in gdat.indxstkscond:
            lgal = np.array([gdat.dictglob['liststkscond'][r]['lgal']])
            bgal = np.array([gdat.dictglob['liststkscond'][r]['bgal']])
            axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, \
                                                marker=gmod.listelemmrkr[l], color='black', alpha=0.1, lw=gdat.mrkrlinewdth)


def retr_colr(gdat, strgstat, strgmodl, indxpopl=None):
    
    if strgmodl == 'true':
        if indxpopl is None:
            colr = gdat.refr.colr
        else:
            colr = gdat.refr.colrelem[indxpopl]
    if strgmodl == 'fitt':
        if strgstat == 'this' or strgstat == 'pdfn':
            if indxpopl is None:
                colr = gmod.colr
            else:
                colr = gmod.colrelem[indxpopl]
        if strgstat == 'mlik':
            colr = 'r'
    
    return colr


def retr_levipost(listllik):
    
    minmlistllik = np.amin(listllik)
    levipost = np.log(np.mean(1. / np.exp(listllik - minmlistllik))) + minmlistllik
    
    return levipost


def retr_infofromlevi(pmeallik, levi):
    
    info = pmeallik - levi

    return info


def retr_jcbn():
    
    fluxpare, lgalpare, bgalpare, fluxauxi, lgalauxi, bgalauxi = sympy.symbols('fluxpare lgalpare bgalpare fluxauxi lgalauxi bgalauxi')
    
    matr = sympy.Matrix([[ fluxpare,      fluxauxi, 0,            0, 0,            0], \
                         [-fluxpare, 1  - fluxauxi, 0,            0, 0,            0], \
                         [-lgalauxi,             0, 1, 1 - fluxauxi, 0,            0], \
                         [-lgalauxi,             0, 1,    -fluxauxi, 0,            0], \
                         [-bgalauxi,             0, 0,            0, 1, 1 - fluxauxi], \
                         [-bgalauxi,             0, 0,            0, 1,    -fluxauxi]])

    jcbn = matr.det()

    return jcbn

# f1 = uf f0
# f2 = (1 - uf) f0
# x1 = x0 + (1 - uf) ux
# x2 = x0 - uf ux
# y1 = y0 + (1 - uf) uy
# y2 = y0 - uf uy

# f1/uf f1/f0 f1/x0 f1/ux f1/y0 f1/uy
# f2/uf f2/f0 f2/x0 f2/ux f2/y0 f2/uy
# x1/uf x1/f0 x1/x0 x1/ux x1/y0 x1/uy
# x2/uf x2/f0 x2/x0 x2/ux x2/y0 x2/uy
# y1/uf y1/f0 y1/x0 y1/ux y1/y0 y1/uy
# y2/uf y2/f0 y2/x0 y2/ux y2/y0 y2/uy

#  f0     uf 0      0 0      0
# -f0 1 - uf 0      0 0      0
# -ux      0 1 1 - uf 0      0
# -ux      0 1    -uf 0      0
# -uy      0 0      0 1 1 - uf
# -uy      0 0      0 1    -uf

# f0
#retr_jcbn()

def retr_angldist(gdat, lgalfrst, bgalfrst, lgalseco, bgalseco):
    
    # temp -- heal does not work when the dimension of lgalfrst is 1
    if gdat.typepixl == 'heal':
        dir1 = np.array([lgalfrst, bgalfrst])
        dir2 = np.array([lgalseco, bgalseco])
        angldist = hp.rotator.angdist(dir1, dir2)
    else:
        angldist = np.sqrt((lgalfrst - lgalseco)**2 + (bgalfrst - bgalseco)**2)

    return angldist


def retr_deflextr(gdat, indxpixlelem, sher, sang):
    
    factcosi = sher * np.cos(2. * sang)
    factsine = sher * np.cos(2. * sang)
    defllgal = factcosi * gdat.lgalgrid[indxpixlelem] + factsine * gdat.bgalgrid[indxpixlelem]
    deflbgal = factsine * gdat.lgalgrid[indxpixlelem] - factcosi * gdat.bgalgrid[indxpixlelem]
    
    return np.vstack((defllgal, deflbgal)).T 


def readfile(path):

    print('Reading %s...' % path)

    filepick = open(path + '.p', 'rb')
    filearry = h5py.File(path + '.h5', 'r')
    gdattemptemp = pickle.load(filepick)
    
    for attr in filearry:
        setattr(gdattemptemp, attr, filearry[attr][()])

    filepick.close()
    filearry.close()
    
    if 'gdatfinl' in path or 'gdatinit' in path:
        if hasattr(gdattemptemp, 'edis') and gdattemptemp.edis is not None and hasattr(gdattemptemp, 'binsener'):
            gdattemptemp.edisintp = sp.interpolate.interp1d(gdattemptemp.binsener, gdattemptemp.edis, fill_value='extrapolate')
        gdattemptemp.adisobjt = sp.interpolate.interp1d(gdattemptemp.redsintp, gdattemptemp.adisintp, fill_value='extrapolate')
        gdattemptemp.redsfromdlosobjt = sp.interpolate.interp1d(gdattemptemp.adisintp * gdattemptemp.redsintp, \
                                                                                        gdattemptemp.redsintp, fill_value='extrapolate')
    
    return gdattemptemp


def init_stat(gdat):
    
    # construct the initial state
    if gdat.typeverb > 0:
        print('Initializing the sampler state...')
        print('inittype')
        print(gdat.inittype)
    
    gmod = gdat.fitt
    
    ## initialization
    ### initialize the unit sample vector randomly
    gmod.this.paragenrunitfull = np.random.rand(gmod.numbparagenrfull)
    gmod.this.paragenrscalfull = np.empty(gmod.numbparagenrfull)

    ## impose user-specified initial state
    ### number of elements
    ## create dummy indxparagenrfullelem 
    gmod.this.indxparagenrfullelem = None
    if gmod.numbparaelem > 0:
        if gdat.inittype == 'refr':
            for l in gmod.indxpopl:
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = gmod.paragenrunitfull[gmod.indxpara.numbelem[l]]
        else:
            for l in gmod.indxpopl:
                if gmod.typemodltran == 'pois':
                    meanelemtemp = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, \
                                                gmod.this.indxparagenrfullelem)[gmod.indxpara.meanelem[l]]
                
                print('temp -- user input is not working for numbelem')
                #namevarb = 'numbelempop%d' % l
                #initvalu = getattr(gmod.init, namevarb)
                #if initvalu > gmod.maxmpara.numbelem[l] or initvalu < gmod.minmpara.numbelem[l]:
                #    raise Exception('Bad initial number of elements...')
                #gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = initvalu
                
                if gmod.typemodltran == 'pois':
                    gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = np.random.poisson(meanelemtemp)
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = round(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]])
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = \
                                        min(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]], gmod.maxmpara.numbelem[l])
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = \
                                        max(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]], gmod.minmpara.numbelem[l])
                gmod.this.paragenrscalfull[gmod.indxpara.numbelem[l]] = gmod.this.paragenrscalfull[gmod.indxpara.numbelem[l]]
    
    if gdat.booldiagmode:
        if gdat.typedata == 'mock' and gdat.inittype == 'refr':
            for l in gmod.indxpopl:
                if gmod.paragenrunitfull[gmod.indxpara.numbelem[l]] > gmod.maxmpara.numbelem[l]:
                    raise Exception('')

    if gmod.numbparaelem > 0:
        gmod.this.indxelemfull = []
        for l in gmod.indxpopl:
            gmod.this.indxelemfull.append(list(range(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]].astype(int))))
        gmod.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmod.this.indxelemfull, 'fitt')

    if gdat.inittype == 'reco':
        if gdat.namerecostat is not None:
            strgcnfg = gdat.namerecostat
        else:
            strgcnfg = gdat.strgcnfg
        path = gdat.pathoutp + 'stat_' + strgcnfg + '.h5'
        if os.path.exists(path):
            boolinitreco = True
            thisfile = h5py.File(path, 'r')
            if gdat.typeverb > 0:
                print('Initializing from the state %s...' % path)
                print('Likelihood:')
                print(thisfile['lliktotl'][...])
                
                # find the number of populations provided
                maxmindxpopl = 0
                for l in range(10):
                    for attr in thisfile:
                        if attr.startswith('lgalpop'):
                            gmod.indxpopl = int(attr[7])
                            if gmod.indxpopl > maxmindxpopl:
                                maxmindxpopl = gmod.indxpopl
                numbpoplinpt = maxmindxpopl + 1
                
                if numbpoplinpt != gmod.numbpopl:
                    print('State file and fitting metamodel have different number of populations.')
                
                # find the number of elements provided
                cntr = np.zeros(gmod.numbpoplinpt, dtype=int)
                for attr in thisfile:
                    if attr.startswith('lgalpop'):
                        gmod.indxpopl = int(attr[7])
                        cntr[indxpopl] += 1
                if gdat.typeverb > 0:
                    print('Number of elements found:')
                    print(cntr)

            for attr in thisfile:
                for k, gmod.nameparagenrbase in enumerate(gmod.nameparagenrbase):
                    if gmod.nameparagenrbase == attr:
                        if gmod.nameparagenrbase.startswith('numbelem'):
                            try:
                                indxpopltemp = int(gmod.nameparagenrbase[-1])
                                initnumbelem = getattr(gdat, 'initnumbelempop%d' % indxpopltemp)
                                print('Initial condition for the number of elements conflicts with the state file. Defaulting to the argument...')
                            except:
                                initnumbelem = thisfile[attr][()]
                            gmod.this.paragenrunitfull[k] = initnumbelem
                        else:
                            gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', thisfile[attr][()], k)
                        if gmod.this.paragenrunitfull[k] == 0.:
                            print('Warning CDF is zero.')
                        if not np.isfinite(thisfile[attr][()]):
                            raise Exception('Retreived state parameter is not finite.')
                        if (gmod.numbparaelem == 0 or gmod.numbparaelem > 0 and not k in gmod.indxpara.numbelem) and \
                                        (not np.isfinite(gmod.this.paragenrunitfull[k]) or gmod.this.paragenrunitfull[k] < 0. or \
                                        gmod.this.paragenrunitfull[k] > 1.):
                            raise Exception('CDF of the retreived state parameter is bad.')
            if gmod.numbparaelem > 0:
                for l in gmod.indxpopl:
                    maxm.numbelem = getattr(gdat.fitt.maxm, 'numbelempop%d' % l)
                    if gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] > maxm.numbelem:
                        gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = maxm.numbelem
                        if gdat.typeverb > 0:
                            print('Tapering off the element list...')

                gmod.this.indxelemfull = []
                for l in gmod.indxpopl:
                    gmod.this.indxelemfull.append(list(range(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]].astype(int))))
                if gdat.typeverb > 0:
                    print('gmod.this.paragenrunitfull[gmod.indxpara.numbelem]')
                    print(gmod.this.paragenrunitfull[gmod.indxpara.numbelem])
            
            gmod.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmod.this.indxelemfull, 'fitt')
            gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrfullelem)
            
            if (gmod.this.paragenrunitfull == 0).all():
                raise Exception('Bad initialization.')
    
            if gmod.numbparaelem > 0 and gmod.this.indxparagenrfullelem is not None:
                for nameparagenrelem in gmod.namepara.elem:
                    initcomp = [[] for l in gmod.indxpopl]
                    for l in gmod.indxpopl:
                        initcomp[l] = np.empty(len(gmod.this.indxelemfull[l]))
                        for k in range(len(gmod.this.indxelemfull[l])):
                            namefiel = '%spop%d%04d' % (nameparagenrelem, l, k)
                            for attr in thisfile:
                                if namefiel == attr:
                                    initcomp[l][k] = thisfile[namefiel][()]
                    setattr(gdat, 'init' + nameparagenrelem, initcomp)
                initcompfromstat(gdat, gdatmodi, 'init')
            thisfile.close()
        else:
            boolinitreco = False
            if gdat.typeverb > 0:
                print('Could not find the state file, %s, to initialize the sampler.' % path)
    
    if gdat.inittype == 'refr':
        if gdat.typedata == 'inpt':
            for l in gmod.indxpopl:
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = gdat.refr.numbelem[l]
        if gdat.typedata == 'mock':
            for k, gmod.nameparagenrbase in enumerate(gmod.nameparagenrbase):
                if not (gdat.inittype == 'pert' and gmod.nameparagenrbase.startswith('numbelem')) and \
                                                                gmod.nameparagenrbase in gmod.nameparagenrbase:
                    gmod.indxpara.true = np.where(gmod.nameparagenrbase == gmod.nameparagenrbase)[0]
                    gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gmodstat.paragenrscalfull[gmod.indxpara.true], k)
        if gmod.numbparaelem > 0:
            gmod.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmod.this.indxelemfull, 'fitt')
        if gdat.typeverb > 1:
            show_paragenrscalfull(gdat, gdatmodi)
        if gmod.this.indxparagenrfullelem is not None:
            print('Initializing elements from the reference element parameters...')
            show_paragenrscalfull(gdat, gdatmodi)
            gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrfullelem)
            show_paragenrscalfull(gdat, gdatmodi)
            initcompfromstat(gdat, gdatmodi, 'refr')
        gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrfullelem)
    
    ## impose user-specified individual initial values
    for k, gmod.nameparagenrbase in enumerate(gmod.nameparagenrbase):
        if gmod.nameparagenrbase.startswith('numbelem'):
            continue
        if gdat.inittype == 'reco' or  gdat.inittype == 'refr' or gdat.inittype == 'pert':
            try:
                getattr(gdat, 'init' + gmod.nameparagenrbase)
                print('Conflicting initial state arguments detected, init keyword takes precedence.')
            except:
                pass
        try:
            raise Exception('')
            initvalu = getattr(gdat, 'init' + gmod.nameparagenrbase)
            gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', initvalu, k)
            if gdat.typeverb > 0:
                print('Received initial condition for %s: %.3g' % (gmod.nameparagenrbase, initvalu))
        except:
            pass
    
    ## PSF
    if gdat.initpsfp is not None:
        print('Initializing the metamodel PSF from the provided initial state...')
        if gdat.initpsfp.size != gmod.indxpara.psfp.size:
            raise Exception('')
        for k, gmod.nameparagenrbase in enumerate(gmod.nameparagenrbase):
            if k in gmod.indxpara.psfp:
                gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gdat.initpsfp[k-gmod.indxpara.psfp[0]], k)
    if gdat.initpsfprefr:
        print('Initializing the metamodel PSF from the reference state...')
        for k, gmod.nameparagenrbase in enumerate(gmod.nameparagenrbase):
            if k in gmod.indxpara.psfp:
                gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gmod.psfpexpr[k-gmod.indxpara.psfp[0]], k)

    if gdat.inittype == 'rand' or gdat.inittype == 'reco' and not boolinitreco:
        if gdat.typeverb > 0:
            print('Initializing from a random state...')
        gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrfullelem)
    
    if gmod.numbparaelem > 0:
        gmod.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmod.this.indxelemfull, 'fitt')

    # check the initial unit sample vector for bad entries
    if gmod.numbparaelem > 0:
        indxsampdiff = np.setdiff1d(gmod.indxparagenrfull, gmod.indxpara.numbelem)
        
        if np.logical_not(np.isfinite(gmod.this.paragenrunitfull[indxsampdiff])).any():
            raise Exception('')
        indxsampbaddlowr = np.where((gmod.this.paragenrunitfull[indxsampdiff] <= 0.) | np.logical_not(np.isfinite(gmod.this.paragenrunitfull[indxsampdiff])))[0]
        indxsampbadduppr = np.where(gmod.this.paragenrunitfull[indxsampdiff] >= 1.)[0]
        indxsampbaddlowr = indxsampdiff[indxsampbaddlowr]
        indxsampbadduppr = indxsampdiff[indxsampbadduppr]
    else:
        indxsampbaddlowr = np.where(gmod.this.paragenrunitfull <= 0.)[0]
        indxsampbadduppr = np.where(gmod.this.paragenrunitfull >= 1.)[0]
    
    indxsampbadd = np.concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print('Initial value caused unit sample vector to go outside the unit interval...')
        show_paragenrscalfull(gdat, gdatmodi, indxsampshow=indxsampbadd)
        gmod.this.paragenrunitfull[indxsampbadd] = np.random.rand(indxsampbadd.size)
        raise Exception('')
    
    gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrfullelem)
    indxbadd = np.where(np.logical_not(np.isfinite(gmod.this.paragenrscalfull)))[0]
    if indxbadd.size > 0:
        raise Exception('')


def writfile(gdattemp, path):
    
    filepick = open(path + '.p', 'wb')
    filearry = h5py.File(path + '.h5', 'w')
    
    gdattemptemp = tdpy.gdatstrt()
    for attr, valu in gdattemp.__dict__.items():
        if attr.endswith('psfnintp'):
            continue
        
        if isinstance(valu, np.ndarray) and valu.dtype != np.dtype('O') and valu.dtype != np.dtype('<U4'):# or isinstance(valu, str) or \
                                       #isinstance(valu, float) or isinstance(valu, bool) or isinstance(valu, int) or isinstance(valu, np.float):
            
            filearry.create_dataset(attr, data=valu)
        else:
            # temp -- make sure interpolation objects are not written.
            if attr != 'adisobjt' and attr != 'redsfromdlosobjt' and attr != 'edisintp':
                setattr(gdattemptemp, attr, valu)
    
    print('Writing to %s...' % path)

    pickle.dump(gdattemptemp, filepick, protocol=pickle.HIGHEST_PROTOCOL)
    filepick.close()
    filearry.close()
   

def retr_deflcutf(angl, defs, asca, acut, asym=False):

    fracanglasca = angl / asca
    
    deflcutf = defs / fracanglasca
    
    # second term in the NFW deflection profile
    fact = np.ones_like(fracanglasca)
    indxlowr = np.where(fracanglasca < 1.)[0]
    indxuppr = np.where(fracanglasca > 1.)[0]
    fact[indxlowr] = np.arccosh(1. / fracanglasca[indxlowr]) / np.sqrt(1. - fracanglasca[indxlowr]**2)
    fact[indxuppr] = np.arccos(1. / fracanglasca[indxuppr]) / np.sqrt(fracanglasca[indxuppr]**2 - 1.)
    
    if asym:
        deflcutf *= np.log(fracanglasca / 2.) + fact
    else:
        fracacutasca = acut / asca
        factcutf = fracacutasca**2 / (fracacutasca**2 + 1)**2 * ((fracacutasca**2 + 1. + 2. * (fracanglasca**2 - 1.)) * fact + \
                np.pi * fracacutasca + (fracacutasca**2 - 1.) * np.log(fracacutasca) + np.sqrt(fracanglasca**2 + fracacutasca**2) * (-np.pi + (fracacutasca**2 - 1.) / fracacutasca * \
                np.log(fracanglasca / (np.sqrt(fracanglasca**2 + fracacutasca**2) + fracacutasca))))
        deflcutf *= factcutf
       
    return deflcutf


def initchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi.this, 'chro' + name, gdat.functime())
    

def stopchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi.this, 'chro' + name, gdat.functime() - getattr(gdatmodi.this, 'chro' + name))


def retr_defl(gdat, indxpixlelem, lgal, bgal, angllens, ellp=None, angl=None, rcor=None, asca=None, acut=None):
    
    # translate the grid
    lgaltran = gdat.lgalgrid[indxpixlelem] - lgal
    bgaltran = gdat.bgalgrid[indxpixlelem] - bgal
    
    if acut is not None:
        defs = angllens
        angl = np.sqrt(lgaltran**2 + bgaltran**2)
        defl = retr_deflcutf(angl, defs, asca, acut)
        defllgal = lgaltran / angl * defl
        deflbgal = bgaltran / angl * defl

    else:
        bein = angllens

        # rotate the grid
        lgalrttr = np.cos(angl) * lgaltran - np.sin(angl) * bgaltran
        bgalrttr = np.sin(angl) * lgaltran + np.cos(angl) * bgaltran
        
        axisrati = 1. - ellp
        facteccc = np.sqrt(1. - axisrati**2)
        factrcor = np.sqrt(axisrati**2 * lgalrttr**2 + bgalrttr**2)
        defllgalrttr = bein * axisrati / facteccc *  np.arctan(facteccc * lgalrttr / factrcor)
        deflbgalrttr = bein * axisrati / facteccc * np.arctanh(facteccc * bgalrttr / factrcor)
        
        # totate back vector to original basis
        defllgal = np.cos(angl) * defllgalrttr + np.sin(angl) * deflbgalrttr
        deflbgal = -np.sin(angl) * defllgalrttr + np.cos(angl) * deflbgalrttr
   
    defl = np.vstack((defllgal, deflbgal)).T
    
    return defl


def retr_lpriselfdist(gdat, strgmodl, feat, strgfeat):
    
    minm = getattr(gmod.minmpara, strgfeat)
    maxm = getattr(gmod.maxmpara, strgfeat)
    
    lpri = np.sum(np.log(pdfn_self(feat, minm, maxm)))
    
    return lpri


def retr_lprilogtdist(gdat, strgmodl, feat, strgfeat):
    
    minm = getattr(gmod.minmpara, strgfeat)
    maxm = getattr(gmod.maxmpara, strgfeat)
    
    lpri = np.sum(np.log(pdfn_logt(feat, minm, maxm)))
    
    return lpri


def retr_lpripowrdist(gdat, strgmodl, feat, strgfeat, paragenrscalfull, l):
    
    gmod = getattr(gdat, strgmodl)
    
    minm = getattr(gmod.minmpara, strgfeat)
    maxm = getattr(gmod.maxmpara, strgfeat)
    
    slop = paragenrscalfull[getattr(gmod.indxpara, 'slopprio' + strgfeat + 'pop%d' % l)]
    
    lpri = np.sum(np.log(pdfn_powr(feat, minm, maxm, slop)))
    
    return lpri


def retr_lpridpowdist(gdat, strgmodl, feat, strgfeat, paragenrscalfull, l):
    
    minm = getattr(gmod.minmpara, strgfeat)
    maxm = getattr(gmod.maxmpara, strgfeat)
    
    brek = paragenrscalfull[getattr(gmod.indxpara, strgfeat + 'distbrek')[l]]
    sloplowr = paragenrscalfull[getattr(gmod.indxpara, 'sloplowrprio' + strgfeat)[l]]
    slopuppr = paragenrscalfull[getattr(gmod.indxpara, 'slopupprprio' + strgfeat)[l]]
    
    lpri = np.sum(np.log(pdfn_dpow(feat, minm, maxm, brek, sloplowr, slopuppr)))
    
    return lpri


def retr_lprigausdist(gdat, strgmodl, feat, strgfeat, paragenrscalfull, l):
    
    distmean = paragenrscalfull[getattr(gmod.indxpara, strgfeat + 'distmean')[l]]
    diststdv = paragenrscalfull[getattr(gmod.indxpara, strgfeat + 'diststdv')[l]]
    
    lpri = np.sum(np.log(pdfn_gaus(feat, distmean, diststdv)))
    
    return lpri


def retr_lpriigamdist(gdat, strgmodl, feat, strgfeat, paragenrscalfull, l):
    
    slop = paragenrscalfull[getattr(gmod.indxpara, strgfeat + 'slop')[l]]
    cutf = getattr(gmod, 'cutf' + strgfeat)
    
    lpri = np.sum(np.log(pdfn_igam(feat, slop, cutf)))

    return lpri


def traptdim(gdat, arry):
    
    s1 = arry[0, 0] + arry[-1, 0] + arry[0, -1] + arry[-1, -1]
    s2 = np.sum(arry[1:-1, 0]) + np.sum(arry[1:-1, -1]) + np.sum(arry[0, 1:-1]) + np.sum(arry[-1, 1:-1])
    s3 = np.sum(arry[1:-1, 1:-1])
    summ = (s1 + 2*s2 + 4*s3) * gdat.apix
    
    return summ


def retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons=None):
    
    pdfnspatprio = pdfnspatpriotemp
    if spatdistcons is not None:
        pdfnspatprio += spatdistcons

    summ = traptdim(gdat, pdfnspatprio)
    pdfnspatprio /= summ
    lpdfspatprio = np.log(pdfnspatprio)
    lpdfspatprioobjt = sp.interpolate.RectBivariateSpline(gdat.binspara.bgalcart, gdat.binspara.lgalcart, lpdfspatprio)
    
    return lpdfspatprio, lpdfspatprioobjt


def retr_gdatobjt(gdat, gdatmodi, strgmodl, boolinit=False):
    
    if strgmodl == 'true':
        gdatobjt = gdat.true
    elif strgmodl == 'fitt' and boolinit:
        gdatobjt = gdat.fitt
    else:
        gdatobjt = gdatmodi

    return gdatobjt


def proc_samp(gdat, gdatmodi, strgstat, strgmodl, fast=False, boolinit=False):
   
    gmod = getattr(gdat, strgmodl)
    
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl, boolinit=boolinit)
    gmodstat = getattr(gdatobjt, strgstat)
    
    initchro(gdat, gdatmodi, 'pars')

    # grab the sample vector
    indxpara = np.arange(gmodstat.paragenrscalfull.size) 

    if gdat.booldiagmode:
        if not np.isfinite(gmodstat.paragenrscalfull).all():
            raise Exception('')

    if gmod.typeevalpsfn != 'none' and (strgmodl == 'true' or boolinit or gdat.boolmodipsfn):
        psfp = gmodstat.paragenrscalfull[gmod.indxpara.psfp]
        if gdat.booldiagmode:
            if np.where(psfp == 0)[0].size == psfp.size:
                raise Exception('')
        setattr(gmodstat, 'psfp', psfp)
    bacp = gmodstat.paragenrscalfull[gmod.indxpara.bacp]
   
    if gmod.numbparaelem > 0:
        
        # temp -- this may slow down execution
        gmodstat.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gmodstat.indxelemfull, strgmodl)

        gmodstat.numbelem = np.empty(gmod.numbpopl, dtype=int)
        indxelem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            gmodstat.numbelem[l] = gmodstat.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int)
            indxelem[l] = np.arange(gmodstat.numbelem[l])
            gmodstat.numbelem[l] = np.sum(gmodstat.numbelem[l])
        gmodstat.numbelemtotl = np.sum(gmodstat.numbelem) 

        gmodstat.dictelem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            gmodstat.dictelem[l] = dict()
            for strgfeat in gmod.namepara.genrelemdefa:
                gmodstat.dictelem[l][strgfeat] = []
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                gmodstat.dictelem[l][nameparagenrelem] = gmodstat.paragenrscalfull[gmodstat.indxparagenrfullelem[l][nameparagenrelem]]
                if gdat.booldiagmode:
                    if ((abs(gmodstat.paragenrscalfull[gmodstat.indxparagenrfullelem[l][nameparagenrelem]]) < 1e-100 ) & (abs(gmodstat.paragenrscalfull[gmodstat.indxparagenrfullelem[l][nameparagenrelem]]) > 0.)).any():
                        raise Exception('')

                    if gmodstat.numbelem[l] != len(gmodstat.dictelem[l][nameparagenrelem]):
                        print('l')
                        print(l)
                        print('numbelem')
                        print(numbelem)
                        print('gmodstat.dictelem')
                        print(gmodstat.dictelem)
                        print('nameparagenrelem')
                        print(nameparagenrelem)
                        raise Exception('')
    
        if gdat.boolbinsener:
            if gdat.typeverb > 2:
                print('Calculating element spectra...')
            initchro(gdat, gdatmodi, 'spec')
            for l in gmod.indxpopl:
                for strgfeat in gmod.namepara.genrelem[l]:
                    sindcolr = [gmodstat.dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], sind=gmodstat.dictelem[l]['sind'], curv=gmodstat.dictelem[l]['curv'], \
                              expc=gmodstat.dictelem[l]['expc'], sindcolr=sindcolr, spectype=gmod.spectype[l])
                    if gmod.typeelem[l].startswith('lghtline'):
                        if gmod.typeelem[l] == 'lghtlinevoig':
                            gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], elin=gmodstat.dictelem[l]['elin'], sigm=gmodstat.dictelem[l]['sigm'], \
                                                                                                 gamm=gmodstat.dictelem[l]['gamm'], spectype=gmod.spectype[l])
                        else:
                            gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], elin=gmodstat.dictelem[l]['elin'], \
                                                                                                edisintp=gdat.edisintp, spectype=gmod.spectype[l])

            stopchro(gdat, gdatmodi, 'spec')
        
        if gdat.typeverb > 2:
            print('Element features:')
            for l in gmod.indxpopl:
                print('l')
                print(l)
                for strgfeat in gmod.namepara.genrelem[l]:
                    print(strgfeat)
                    print(gmodstat.dictelem[l][strgfeat])
    
        if gdat.booldiagmode:
            for l in gmod.indxpopl:
                for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    if (gmod.listscalparagenrelem[l][g] != 'gaus' and not gmod.listscalparagenrelem[l][g].startswith('lnor')) and  \
                       (gmod.listscalparagenrelem[l][g] != 'expo' and (gmodstat.dictelem[l][nameparagenrelem] < getattr(gmod.minmpara, nameparagenrelem)).any()) or \
                                        (gmodstat.dictelem[l][nameparagenrelem] > getattr(gmod.maxmpara, nameparagenrelem)).any():
                        
                        print('l, g')
                        print(l, g)
                        print('nameparagenrelem')
                        print(nameparagenrelem)
                        print('gmodstat.dictelem[l][nameparagenrelem]')
                        summgene(gmodstat.dictelem[l][nameparagenrelem])
                        print('getattr(gmod, minm + nameparagenrelem)')
                        print(getattr(gmod.minmpara, nameparagenrelem))
                        print('getattr(gmod, maxm + nameparagenrelem)')
                        print(getattr(gmod.maxmpara, nameparagenrelem))
                        print('gmod.listscalparagenrelem[l][g]')
                        print(gmod.listscalparagenrelem[l][g])
                        raise Exception('')
           
        # calculate element spectra
        # temp
        if gdat.booldiagmode:
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lens':
                    if gdat.variasca:
                        indx = np.where(gmodstat.paragenrscalfull[gmodstat.indxparagenrfullelem[l]['acut']] < 0.)[0]
                        if indx.size > 0:
                            raise Exception('')
                    if gdat.variacut:
                        indx = np.where(gmodstat.paragenrscalfull[gmodstat.indxparagenrfullelem[l]['asca']] < 0.)[0]
                        if indx.size > 0:
                            raise Exception('')
    
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght'):
                    
                # evaluate horizontal and vertical position for elements whose position is a power law in image-centric radius
                if gmod.typespatdist[l] == 'glc3':
                    gmodstat.dictelem[l]['dlos'], gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal'] = retr_glc3(gmodstat.dictelem[l]['dglc'], \
                                                                                                        gmodstat.dictelem[l]['thet'], gmodstat.dictelem[l]['phii'])
                
                if gmod.typespatdist[l] == 'gangexpo':
                    gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal'], = retr_lgalbgal(gmodstat.dictelem[l]['gang'], \
                                                                                                        gmodstat.dictelem[l]['aang'])
                    
                    if gdat.booldiagmode:
                        if gmodstat.numbelem[l] > 0:
                            if np.amin(gmodstat.dictelem[l]['lgal']) < gmod.minmlgal or \
                               np.amax(gmodstat.dictelem[l]['lgal']) > gmod.maxmlgal or \
                               np.amin(gmodstat.dictelem[l]['bgal']) < gmod.minmbgal or \
                               np.amax(gmodstat.dictelem[l]['bgal']) > gmod.maxmbgal:
                                raise Exception('Bad coordinates!')

                if gmod.typespatdist[l] == 'los3':
                    gmodstat.dictelem[l]['dglc'], gmodstat.dictelem[l]['thet'], gmodstat.dictelem[l]['phii'] = retr_los3(gmodstat.dictelem[l]['dlos'], \
                                                                                                        gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal'])

                # evaluate flux for pulsars
                if gmod.typeelem[l] == 'lghtpntspuls':
                    gmodstat.dictelem[l]['lumi'] = retr_lumipuls(gmodstat.dictelem[l]['geff'], gmodstat.dictelem[l]['magf'], gmodstat.dictelem[l]['per0'])
                if gmod.typeelem[l] == 'lghtpntsagnntrue':
                    gmodstat.dictelem[l]['reds'] = gdat.redsfromdlosobjt(gmodstat.dictelem[l]['dlos'])
                    gmodstat.dictelem[l]['lumi'] = gmodstat.dictelem[l]['lum0'] * (1. + gmodstat.dictelem[l]['reds'])**4
                if gmod.typeelem[l] == 'lghtpntspuls' or gmod.typeelem[l] == 'lghtpntsagnntrue':
                    gmodstat.dictelem[l]['flux'] = retr_flux(gdat, gmodstat.dictelem[l]['lumi'], gmodstat.dictelem[l]['dlos'])
                # evaluate spectra
                if gmod.typeelem[l].startswith('lghtline'):
                    if gmod.typeelem[l] == 'lghtlinevoig':
                        gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], elin=gmodstat.dictelem[l]['elin'], sigm=gmodstat.dictelem[l]['sigm'], \
                                                                                                          gamm=gmodstat.dictelem[l]['gamm'], spectype=gmod.spectype[l])
                    else:
                        gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], elin=gmodstat.dictelem[l]['elin'], edisintp=gdat.edisintp, spectype=gmod.spectype[l])
                else:
                    sindcolr = [gmodstat.dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], sind=gmodstat.dictelem[l]['sind'], curv=gmodstat.dictelem[l]['curv'], \
                                       expc=gmodstat.dictelem[l]['expc'], sindcolr=sindcolr, spectype=gmod.spectype[l])
        

    stopchro(gdat, gdatmodi, 'pars')
    
    ### loglikelihood
    initchro(gdat, gdatmodi, 'modl')
    
    if gmod.boollens:
        lgalsour = gmodstat.paragenrscalfull[gmod.indxpara.lgalsour]
        bgalsour = gmodstat.paragenrscalfull[gmod.indxpara.bgalsour]
    
    if gdat.typeverb > 2:
        print('Evaluating the likelihood...')
    
    # process a sample vector and the occupancy list to calculate secondary variables
    if gmod.boollens:
        fluxsour = gmodstat.paragenrscalfull[gmod.indxpara.fluxsour]
        if gdat.numbener > 1:
            sindsour = gmodstat.paragenrscalfull[gmod.indxpara.sindsour]
        sizesour = gmodstat.paragenrscalfull[gmod.indxpara.sizesour]
        ellpsour = gmodstat.paragenrscalfull[gmod.indxpara.ellpsour]
        anglsour = gmodstat.paragenrscalfull[gmod.indxpara.anglsour]
    if gmod.typeemishost != 'none':
        lgalhost = [[] for e in gmod.indxsersfgrd]
        bgalhost = [[] for e in gmod.indxsersfgrd]
        fluxhost = [[] for e in gmod.indxsersfgrd]
        if gdat.numbener > 1:
            sindhost = [[] for e in gmod.indxsersfgrd]
        sizehost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            lgalhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'lgalhostisf%d' % e)]
            bgalhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'bgalhostisf%d' % e)]
            fluxhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'fluxhostisf%d' % e)]
            if gdat.numbener > 1:
                sindhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sindhostisf%d' % e)]
            sizehost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sizehostisf%d' % e)]
    if gmod.boollens:
        beinhost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            beinhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % e)]
    if gmod.typeemishost != 'none':
        ellphost = [[] for e in gmod.indxsersfgrd]
        anglhost = [[] for e in gmod.indxsersfgrd]
        serihost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            ellphost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'ellphostisf%d' % e)]
            anglhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'anglhostisf%d' % e)]
            serihost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'serihostisf%d' % e)]
    if gmod.boollens:
        numbpixltemp = gdat.numbpixlcart
        defl = np.zeros((numbpixltemp, 2))
        
    # determine the indices of the pixels over which element kernels will be evaluated
    if gdat.boolbinsspat:
        if gmod.numbparaelem > 0:
            listindxpixlelem = [[] for l in gmod.indxpopl]
            listindxpixlelemconc = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                if gmodstat.numbelem[l] > 0:
                    listindxpixlelem[l], listindxpixlelemconc[l] = retr_indxpixlelemconc(gdat, strgmodl, gmodstat.dictelem, l)
                    
    if gmod.boollens:
        sherextr = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sherextr')]
        sangextr = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sangextr')]
       
        ## host halo deflection
        initchro(gdat, gdatmodi, 'deflhost')
        deflhost = [[] for e in gmod.indxsersfgrd]
            
        indxpixlmiss = gdat.indxpixlcart

        for e in gmod.indxsersfgrd:
            if gdat.typeverb > 2:
                print('Evaluating the deflection field due to host galaxy %d' % e)
                print('lgalhost[e]')
                print(lgalhost[e])
                print('bgalhost[e]')
                print(bgalhost[e])
                print('beinhost[e]')
                print(beinhost[e])
                print('ellphost[e]')
                print(ellphost[e])
                print('anglhost[e]')
                print(anglhost[e])

            deflhost[e] = retr_defl(gdat, indxpixlmiss, lgalhost[e], bgalhost[e], beinhost[e], ellp=ellphost[e], angl=anglhost[e])
             
            if gdat.booldiagmode:
                indxpixltemp = slice(None)
            
            setattr(gmodstat, 'deflhostisf%d' % e, deflhost[e])
       
            if gdat.typeverb > 2:
                print('deflhost[e]')
                summgene(deflhost[e])
                
            defl += deflhost[e]
            if gdat.typeverb > 2:
                print('After adding the host deflection...')
                print('defl')
                summgene(defl)
        if gdat.booldiagmode:
            if not np.isfinite(deflhost).all():
                raise Exception('')
        
        stopchro(gdat, gdatmodi, 'deflhost')

        ## external shear
        initchro(gdat, gdatmodi, 'deflextr')
        deflextr = []
        indxpixltemp = gdat.indxpixlcart
        deflextr = retr_deflextr(gdat, indxpixltemp, sherextr, sangextr)
        defl += deflextr
        if gdat.typeverb > 2:
            print('After adding the external deflection...')
            print('defl')
            summgene(defl)
        stopchro(gdat, gdatmodi, 'deflextr')
    
    # Boolean flag to indicate that the object to convolve the image will be needed
    boolneedpsfnconv = gdat.typepixl == 'cart' and (gmod.typeevalpsfn == 'conv' or gmod.typeevalpsfn == 'full')
    
    ## Boolean flag to indicate that the object to convolve the image will be constructed
    boolcalcpsfnconv = strgmodl == 'true' or boolinit or gdat.boolmodipsfn
    
    # get the convolution object
    if boolneedpsfnconv and boolcalcpsfnconv:
        initchro(gdat, gdatmodi, 'psfnconv')
        if gdat.typeverb > 2:
            print('Evaluating the PSF convolution kernel...')
        psfnconv = [[[] for i in gdat.indxener] for m in gdat.indxevtt]
        if gdat.typepixl == 'cart':
            
            gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binspara.angl, gmod.typemodlpsfn, strgmodl)
            fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            for mm, m in enumerate(gdat.indxevtt):
                for ii, i in enumerate(gdat.indxener):
                    if gmod.typemodlpsfn == 'singgaus':
                        sigm = psfp[i+m*gdat.numbener]
                    else:
                        sigm = fwhm[i, m] / 2.355
                    gmodstat.psfnconv[mm][ii] = AiryDisk2DKernel(sigm / gdat.sizepixl)
        
        stopchro(gdat, gdatmodi, 'psfnconv')
    
    if (gmod.typeevalpsfn == 'kern' or gmod.typeevalpsfn == 'full') and gmod.numbparaelem > 0:
        if strgmodl == 'true' or boolinit or gdat.boolmodipsfn:
            if gdat.typepixl == 'heal':
                gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binspara.angl, gmod.typemodlpsfn, strgmodl)
                gmodstat.psfnintp = sp.interpolate.interp1d(gdat.binspara.angl, gmodstat.psfn, axis=1, fill_value='extrapolate')
                fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            if gdat.typepixl == 'cart':
                if gdat.kernevaltype == 'ulip':
                    gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binspara.angl, gmod.typemodlpsfn, strgmodl)
                    gmodstat.psfnintp = sp.interpolate.interp1d(gdat.binspara.angl, gmodstat.psfn, axis=1, fill_value='extrapolate')
                    if gdat.booldiagmode:
                        if not np.isfinite(gmodstat.psfnintp(0.05)).all():
                            raise Exception('')

                if gdat.kernevaltype == 'bspx':
                    
                    gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binspara.anglcart.flatten(), gmod.typemodlpsfn, strgmodl)
                    
                    # side length of the upsampled kernel
                    gdat.numbsidekernusam = 100
                    # side length of the original kernel
                    gdat.numbsidekern = gdat.numbsidekernusam / factkernusam 
                    gdat.indxsidekern = np.arange(gdat.numbsidekern)

    	        	# pad by one row and one column
    	        	#psf = np.zeros((gdat.numbsidekernusam+1, gdat.numbsidekernusam+1))
    	        	#psf[0:gdat.numbsidekernusam, 0:gdat.numbsidekernusam] = psf0
		        	
    	        	# make design matrix for each factkernusam x factkernusam region
                    nx = factkernusam + 1
                    y, x = mgrid[0:nx, 0:nx] / float(factkernusam)
                    x = x.flatten()
                    y = y.flatten()
                    kernmatrdesi = np.array([full(nx*nx, 1), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y]).T
    	        	
                    # output np.array of coefficients
                    gmodstat.psfnintp = np.empty((gdat.numbsidekern, gdat.numbsidekern, kernmatrdesi.shape[1]))

    	        	# solve p = kernmatrdesi psfnintp for psfnintp
                    for iy in gdat.indxsidekern:
                        for ix in gdat.indxsidekern:
                            p = psf[iy*factkernusam:(iy+1)*factkernusam+1, ix*factkernusam:(ix+1)*factkernusam+1].flatten()
                            gmodstat.psfnintp[iy, ix, :] = dot(linalg.inv(dot(kernmatrdesi.T, kernmatrdesi)), dot(kernmatrdesi.T, p))
        else:
            gmodstat.psfnintp = gdat.fitt.this.psfnintp
    sbrt = dict()
    for name in gmod.listnamediff:
        sbrt[name] = []
        
    if gmod.numbparaelem > 0:
        if gmod.boolelemsbrtdfncanyy:
            sbrtdfnc = []
        if gmod.boolelemsbrtextsbgrdanyy:
            sbrtextsbgrd = []
        if gmod.boolelemdeflsubhanyy:
            deflsubh = []
        # retrieve or initialize state variable
        if gmod.boolelemsbrtdfncanyy:
            sbrtdfnc = np.zeros_like(gdat.expo)
        if gmod.boolelemdeflsubhanyy:
            deflsubh = np.zeros((gdat.numbpixl, 2))
        if gmod.boolelemsbrtextsbgrdanyy: 
            sbrtextsbgrd = np.zeros_like(gdat.expo)
        
        # element kernel evaluation
        if gmod.boolelemsbrtdfncanyy:
            initchro(gdat, gdatmodi, 'elemsbrtdfnc')
            sbrt['dfnc'] = []
            for l in gmod.indxpopl:
                if gmod.boolelemsbrtdfnc[l]:
                    for k in range(gmodstat.numbelem[l]):
                        if gmod.boolelemlght[l]:
                            varbamplextd = gmodstat.dictelem[l]['spec'][:, k]
                        if gmod.typeelem[l].startswith('clus'):
                            varbamplextd = gmodstat.dictelem[l]['nobj'][None, k]
                        if gmod.typeelem[l] == 'clusvari':
                            sbrtdfnc[0, listindxpixlelem[l][k], 0] += gmodstat.dictelem[l]['nobj'][k] / 2. / np.pi / gmodstat.dictelem[l]['gwdt'][k]**2 * \
                                np.exp(-0.5 * ((gmodstat.dictelem[l]['lgal'][k] - gdat.lgalgrid[listindxpixlelem[l][k]])**2 + \
                                    (gmodstat.dictelem[l]['bgal'][k] - gdat.bgalgrid[listindxpixlelem[l][k]])**2) / gmodstat.dictelem[l]['gwdt'][k]**2)
                            
                        if gmod.boolelempsfn[l]:
                            print('sbrtdfnc')
                            summgene(sbrtdfnc)
                            sbrtdfnc[:, listindxpixlelem[l][k], :] += retr_sbrtpnts(gdat, gmodstat.dictelem[l]['lgal'][k], \
                                                             gmodstat.dictelem[l]['bgal'][k], varbamplextd, gmodstat.psfnintp, listindxpixlelem[l][k])
                        
                        if gmod.typeelem[l].startswith('lghtline'):
                            sbrtdfnc[:, 0, 0] += gmodstat.dictelem[l]['spec'][:, k]
                        
            sbrt['dfnc'] = sbrtdfnc
            
            if gdat.booldiagmode:
                if not np.isfinite(sbrtdfnc).all():
                    raise Exception('Element delta function brightness not finite.')

            setattr(gmodstat, 'sbrtdfnc', sbrt['dfnc'])

            if gdat.booldiagmode:
                cntppntschec = retr_cntp(gdat, sbrt['dfnc'])
                numbelemtemp = 0
                for l in gmod.indxpopl:
                    if gmod.boolelemsbrtdfnc[l]:
                        numbelemtemp += np.sum(gmodstat.numbelem[l])
                if np.amin(cntppntschec) < -0.1:
                    raise Exception('Point source spectral surface brightness is not positive-definite.')
            
            stopchro(gdat, gdatmodi, 'elemsbrtdfnc')
        
        if gmod.boolelemdeflsubhanyy:
            initchro(gdat, gdatmodi, 'elemdeflsubh')
            if gdat.typeverb > 2:
                print('Perturbing subhalo deflection field')
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lens':
                    for kk, k in enumerate(indxelem[l]):
                        asca = gmodstat.dictelem[l]['asca'][k]
                        acut = gmodstat.dictelem[l]['acut'][k]
                        if gmod.typeelemspateval[l] == 'locl':
                            indxpixl = listindxpixlelem[l][kk]
                        else:
                            indxpixl = gdat.indxpixl
                        deflsubh[indxpixl, :] += retr_defl(gdat, indxpixl, \
                                                     gmodstat.dictelem[l]['lgal'][kk], gmodstat.dictelem[l]['bgal'][kk], gmodstat.dictelem[l]['defs'][kk], \
                                                     asca=asca, acut=acut)
            
                    # temp -- find out what is causing the features in the element convergence maps
                    #for kk, k in enumerate(indxelem[l]):
                    #    indxpixlpnts = retr_indxpixl(gdat, gmodstat.dictelem[l]['bgal'][kk], gmodstat.dictelem[l]['lgal'][kk])
                    #    if deflsubh[listindxpixlelem[l][kk], :]
            
            if gdat.typeverb > 2:
                print('deflsubh')
                summgene(deflsubh)
            setattr(gmodstat, 'deflsubh', deflsubh)
            
            if gdat.booldiagmode:
                if not np.isfinite(deflsubh).all():
                    raise Exception('Element deflection is not finite.')

            defl += deflsubh
            if gdat.typeverb > 2:
                print('After adding subhalo deflection to the total deflection')
                print('defl')
                summgene(defl)

            stopchro(gdat, gdatmodi, 'elemdeflsubh')

        if gmod.boolelemsbrtextsbgrdanyy:
            initchro(gdat, gdatmodi, 'elemsbrtextsbgrd')
            if strgstat == 'this':
                for l in gmod.indxpopl:
                    if gmod.typeelem[l] == 'lghtgausbgrd':
                        for k in range(gmodstat.numbelem[l]):
                            sbrtextsbgrd[:, listindxpixlelem[l][k], :] += gmodstat.dictelem[l]['spec'][:, k, None, None] / \
                                    2. / np.pi / gmodstat.dictelem[l]['gwdt'][k]**2 * \
                                    np.exp(-0.5 * ((gmodstat.dictelem[l]['lgal'][k] - gdat.lgalgrid[None, listindxpixlelem[l][k], None])**2 + \
                                                (gmodstat.dictelem[l]['bgal'][k] - gdat.bgalgrid[None, listindxpixlelem[l][k], None])**2) / gmodstat.dictelem[l]['gwdt'][k]**2)
                
                setattr(gmodstat, 'sbrtextsbgrd', sbrtextsbgrd)
            sbrt['extsbgrd'] = []
            sbrt['extsbgrd'] = sbrtextsbgrd
            
            if gdat.booldiagmode:
                cntppntschec = retr_cntp(gdat, sbrt['extsbgrd'])
                if np.amin(cntppntschec) < -0.1:
                    raise Exception('Point source spectral surface brightness is not positive-definite.')
        
            stopchro(gdat, gdatmodi, 'elemsbrtextsbgrd')
    
        if gdat.typeverb > 2:
            print('Element related state variables after perturbations...')
            if gmod.boolelemsbrtdfncanyy:
                print('sbrtdfnc')
                summgene(sbrtdfnc)
            if gmod.boolelemdeflsubhanyy:
                print('deflsubh')
                summgene(deflsubh)
            if gmod.boolelemsbrtextsbgrdanyy:
                print('sbrtextsbgrd')
                summgene(sbrtextsbgrd)
    
    if gmod.boollens:
        
        # lensed surface brightness
        initchro(gdat, gdatmodi, 'sbrtlens')
        
        if gdat.typeverb > 2:
            print('Evaluating lensed surface brightness...')
        
        if strgstat == 'this' or gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
            sbrt['bgrd'] = []
        if gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
            sbrt['bgrdgalx'] = []
        
        if gdat.numbener > 1:
            specsour = retr_spec(gdat, np.array([fluxsour]), sind=np.array([sindsour]))
            if gdat.typeverb > 2:
                print('sindsour')
                print(sindsour)
        else:
            specsour = np.array([fluxsour])
        
        if gdat.typeverb > 2:
            print('lgalsour')
            print(lgalsour)
            print('bgalsour')
            print(bgalsour)
            print('sizesour')
            print(sizesour)
            print('ellpsour')
            print(ellpsour)
            print('anglsour')
            print(anglsour)
            print('fluxsour')
            print(fluxsour)
            print('specsour')
            print(specsour)

        if gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
        
            if gdat.typeverb > 2:
                print('Interpolating the background emission...')

            sbrt['bgrdgalx'] = retr_sbrtsers(gdat, gdat.lgalgrid[indxpixlelem[0]], gdat.bgalgrid[indxpixlelem[0]], \
                                                                            lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            if gdat.typeverb > 2:
                print('sbrt[bgrdgalx]')
                summgene(sbrt['bgrdgalx'])
                print('sbrtextsbgrd')
                summgene(sbrtextsbgrd)
            sbrt['bgrd'] = sbrt['bgrdgalx'] + sbrtextsbgrd
        
            sbrt['lens'] = np.empty_like(gdat.cntpdata)
            for ii, i in enumerate(gdat.indxener):
                for mm, m in enumerate(gdat.indxevtt):
                    sbrtbgrdobjt = sp.interpolate.RectBivariateSpline(gdat.meanpara.bgalcart, gdat.meanpara.lgalcart, \
                                                            sbrt['bgrd'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    bgalprim = gdat.bgalgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 1]
                    lgalprim = gdat.lgalgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 0]
                    # temp -- T?
                    sbrt['lens'][ii, :, m] = sbrtbgrdobjt(bgalprim, lgalprim, grid=False).flatten()
        else:
            if gdat.typeverb > 2:
                print('Not interpolating the background emission...')
            
            sbrt['lens'] = retr_sbrtsers(gdat, gdat.lgalgrid - defl[gdat.indxpixl, 0], \
                                                   gdat.bgalgrid - defl[gdat.indxpixl, 1], \
                                                   lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            sbrt['bgrd'] = retr_sbrtsers(gdat, gdat.lgalgrid, \
                                                   gdat.bgalgrid, \
                                                   lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
        setattr(gmodthis, 'sbrtlens', sbrt['lens'])

        if gdat.booldiagmode:
            if not np.isfinite(sbrt['lens']).all():
                raise Exception('Lensed emission is not finite.')
            if (sbrt['lens'] == 0).all():
                raise Exception('Lensed emission is zero everynp.where.')

        stopchro(gdat, gdatmodi, 'sbrtlens')
        
    ### background surface brightness
    sbrtback = []
    # temp
    #sbrtback = np.empty((numbback, gdat.numbener, indxpixlelem[yy].size, gdat.numbevtt))
    
    # evaluate host galaxy surface brightness
    if gmod.typeemishost != 'none':
        initchro(gdat, gdatmodi, 'sbrthost')
        for e in gmod.indxsersfgrd:
            if gdat.typeverb > 2:
                print('Evaluating the host galaxy surface brightness...')
            if gdat.numbener > 1:
                spechost = retr_spec(gdat, np.array([fluxhost[e]]), sind=np.array([sindhost[e]]))
            else:
                spechost = np.array([fluxhost[e]])
            
            if gdat.typeverb > 2:
                print('lgalhost[e]')
                print(lgalhost[e] * gdat.anglfact)
                print('bgalhost[e]')
                print(bgalhost[e] * gdat.anglfact)
                print('spechost')
                print(spechost)
                print('sizehost[e]')
                print(sizehost[e])
                print('ellphost[e]')
                print(ellphost[e])
                print('anglhost[e]')
                print(anglhost[e])
                print('serihost[e]')
                print(serihost[e])
            
            sbrt['hostisf%d' % e] = retr_sbrtsers(gdat, gdat.lgalgrid, gdat.bgalgrid, lgalhost[e], \
                                                                        bgalhost[e], spechost, sizehost[e], ellphost[e], anglhost[e], serihost[e])
            
            setattr(gmodstat, 'sbrthostisf%d' % e, sbrt['hostisf%d' % e])
                
        #sbrthost = sbrt['host']
        if gdat.typeverb > 2:
            for e in gmod.indxsersfgrd:
                print('e')
                print(e)
                print('sbrt[hostisf%d]')
                summgene(sbrt['hostisf%d' % e])
        stopchro(gdat, gdatmodi, 'sbrthost')
    
    ## model emission
    initchro(gdat, gdatmodi, 'sbrtmodl')
    if gdat.typeverb > 2:
        print('Summing up the model emission...')
    
    sbrt['modlraww'] = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt))
    for name in gmod.listnamediff:
        if name.startswith('back'):
            gmod.indxbacktemp = int(name[4:8])
            
            if gdat.typepixl == 'heal' and (gmod.typeevalpsfn == 'full' or gmod.typeevalpsfn == 'conv') and not gmod.boolunifback[gmod.indxbacktemp]:
                sbrttemp = getattr(gmod, 'sbrtbackhealfull')[gmod.indxbacktemp]
            else:
                sbrttemp = gmod.sbrtbacknorm[gmod.indxbacktemp]
           
            if gmod.boolspecback[gmod.indxbacktemp]:
                sbrt[name] = sbrttemp * bacp[gmod.indxbacpback[gmod.indxbacktemp]]
            else:
                sbrt[name] = sbrttemp * bacp[gmod.indxbacpback[gmod.indxbacktemp][gdat.indxener]][:, None, None]
        
        sbrt['modlraww'] += sbrt[name]
        if gdat.booldiagmode:
            if np.amax(sbrttemp) == 0.:
                raise Exception('')

        if gdat.typeverb > 2:
            print('name')
            print(name)
            print('sbrt[name]')
            summgene(sbrt[name])
    if gdat.typeverb > 2:
        for ii, i in enumerate(gdat.indxener):
            print('ii, i')
            print(ii, i)
            for mm, m in enumerate(gdat.indxevtt):
                print('mm, m')
                print(mm, m)
                print('sbrt[modlraww][ii, :, mm]')
                summgene(sbrt['modlraww'][ii, :, mm])
    
    # convolve the model with the PSF
    if gmod.convdiffanyy and (gmod.typeevalpsfn == 'full' or gmod.typeevalpsfn == 'conv'):
        sbrt['modlconv'] = []
        # temp -- isotropic background proposals are unnecessarily entering this clause
        if gdat.typeverb > 2:
            print('Convolving the model image with the PSF...') 
        sbrt['modlconv'] = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for ii, i in enumerate(gdat.indxener):
            for mm, m in enumerate(gdat.indxevtt):
                if gdat.strgcnfg == 'pcat_ferm_igal_mock_test':
                    print('Convolving ii, i, mm, m')
                    print(ii, i, mm, m)
                if gdat.typepixl == 'cart':
                    if gdat.numbpixl == gdat.numbpixlcart:
                        sbrt['modlconv'][ii, :, mm] = convolve_fft(sbrt['modlraww'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)), \
                                                                                                                             psfnconv[mm][ii]).flatten()
                    else:
                        sbrtfull = np.zeros(gdat.numbpixlcart)
                        sbrtfull[gdat.indxpixlrofi] = sbrt['modlraww'][ii, :, mm]
                        sbrtfull = sbrtfull.reshape((gdat.numbsidecart, gdat.numbsidecart))
                        sbrt['modlconv'][ii, :, mm] = convolve_fft(sbrtfull, psfnconv[mm][ii]).flatten()[gdat.indxpixlrofi]
                    indx = np.where(sbrt['modlconv'][ii, :, mm] < 1e-50)
                    sbrt['modlconv'][ii, indx, mm] = 1e-50
                if gdat.typepixl == 'heal':
                    sbrt['modlconv'][ii, :, mm] = hp.smoothing(sbrt['modlraww'][ii, :, mm], fwhm=fwhm[i, m])[gdat.indxpixlrofi]
                    sbrt['modlconv'][ii, :, mm][np.where(sbrt['modlraww'][ii, :, mm] <= 1e-50)] = 1e-50
        
        setattr(gmodstat, 'sbrtmodlconv', sbrt['modlconv'])
        # temp -- this could be made faster -- need the copy() statement because sbrtdfnc gets added to sbrtmodl afterwards
        sbrt['modl'] = np.copy(sbrt['modlconv'])
    else:
        if gdat.typeverb > 2:
            print('Skipping PSF convolution of the model...')
        sbrt['modl'] = np.copy(sbrt['modlraww'])
    
    if gdat.typeverb > 2:
        print('sbrt[modl]')
        summgene(sbrt['modl'])

    ## add PSF-convolved delta functions to the model
    if gmod.numbparaelem > 0 and gmod.boolelemsbrtdfncanyy:
        if gdat.typeverb > 2:
            print('Adding delta functions into the model...')
            print('sbrt[dfnc]')
            summgene(sbrt['dfnc'])
        sbrt['modl'] += sbrt['dfnc']
    stopchro(gdat, gdatmodi, 'sbrtmodl')
    
    if gdat.typeverb > 2:
        print('sbrt[modl]')
        summgene(sbrt['modl'])

    ### count map
    initchro(gdat, gdatmodi, 'expo')
    cntp = dict()
    cntp['modl'] = retr_cntp(gdat, sbrt['modl'])
    
    if gdat.booldiagmode:
        setattr(gmodstat, 'cntpmodl', cntp['modl'])
    stopchro(gdat, gdatmodi, 'expo')

    # mock data specific
    if strgmodl == 'true' and strgstat == 'this':
        
        # generate count data
        cntptemp = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxevtt:
                    cntptemp[i, j, m] = np.random.poisson(cntp['modl'][i, j, m])
        setattr(gdat, 'cntpdata', cntptemp)
    
        if not gdat.boolsqzeexpo and np.amax(cntptemp) == 0:
            print('cntp[modl]')
            summgene(cntp['modl'])
            print('gdat.boolsqzeexpo')
            print(gdat.boolsqzeexpo)
            print('cntptemp')
            summgene(cntptemp)
            raise Exception('Data is zero.')
        
        proc_cntpdata(gdat)
    
    ## diagnostics
    if gdat.booldiagmode:
        frac = cntp['modl'] / np.mean(cntp['modl'])
        if np.amin(frac) < -1e-3 and np.amin(cntp['modl']) < -0.1:
            raise Exception('')
        
        indxcubebadd = np.where(cntp['modl'] < 0.)[0]
        if indxcubebadd.size > 0:
            print('Warning! Model prediction is negative. Correcting to 1e-20...')
            cntp['modl'][indxcubebadd] = 1e-20
    stopchro(gdat, gdatmodi, 'modl')

    # log-prior
    initchro(gdat, gdatmodi, 'lpri')
    if gdat.typeverb > 2:
        print('Evaluating the prior...')
        
    lpri = np.zeros(gmod.numblpri)
    if gmod.numbparaelem > 0:
        
        for l in gmod.indxpopl:
            lpri[0] -= 0.5 * gdat.priofactdoff * gmod.numbparagenrelemsing[l] * gmodstat.numbelem[l]
        
        if gdat.penalpridiff:
            sbrtdatapnts = gdat.sbrtdata - sbrt['dfnc']
            if gdat.typepixl == 'heal':
                raise Exception('')
            if gdat.typepixl == 'cart':
                psecodimdatapnts = np.empty((gdat.numbener, gdat.numbsidecarthalf, gdat.numbevtt))
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binspara.angl, gmod.typemodlpsfn, strgmodl)
                fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
                sigm = fwhm / 2.355
                psecodimdatapntsprio = np.exp(-2. * gdat.meanpara.mpolodim[None, :, None] / (0.1 / sigm[:, None, :]))
                lpridiff = 0.
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        psecdatapnts = retr_psec(gdat, sbrtdatapnts[i, :, m])
                        psecodimdatapnts[i, :, m] = retr_psecodim(gdat, psecdatapnts)
                        psecodimdatapnts[i, :, m] /= psecodimdatapnts[i, 0, m]
                        lpridiff += -0.5 * np.sum((psecodimdatapnts[i, :, m] - psecodimdatapntsprio[i, :, m])**2)
                        setattr(gmodstat, 'psecodimdatapntsen%02devt%d' % (i, m), psecodimdatapnts[i, :, m])
                        setattr(gmodstat, 'psecodimdatapntsprioen%02devt%d'% (i, m), psecodimdatapntsprio[i, :, m])
            lpri[1] = lpridiff 
            setattr(gmodstat, 'lpridiff', lpridiff)
        
        if gmod.typemodltran == 'pois':
            meanelem = gmodstat.paragenrscalfull[gmod.indxpara.meanelem]
            for l in gmod.indxpopl:
                lpri[2] += retr_lprbpois(gmodstat.numbelem[l], meanelem[l])
        
        for l in gmod.indxpopl:
            for g, (strgfeat, strgpdfn) in enumerate(zip(gmod.namepara.genrelem[l], gmod.listscalparagenrelem[l])):
                indxlpritemp = 3 + l * gmod.numbparagenrelem + g
                lpri[indxlpritemp] = retr_lprielem(gdat, strgmodl, l, g, strgfeat, strgpdfn, gmodstat.paragenrscalfull, gmodstat.dictelem, gmodstat.numbelem)
    lpritotl = np.sum(lpri)
    
    if gdat.typeverb > 1:
        print('lpritotl')
        print(lpritotl)
    
    ### log-likelihood
    initchro(gdat, gdatmodi, 'llik')
    llik = retr_llik(gdat, strgmodl, cntp['modl'])
    
    if gdat.typeverb > 2:
        print('cntp[modl]')
        summgene(cntp['modl'])
        print('np.sum(cntp[modl], (1, 2))')
        print(np.sum(cntp['modl'], (1, 2)))
        print('np.sum(gdat.cntpdata, (1, 2))')
        print(np.sum(gdat.cntpdata, (1, 2)))

    if gdat.booldiagmode:
        if not np.isfinite(llik).all():
            raise Exception('Likelihood is not finite.')
    
    gmodstat.lliktotl = np.sum(llik)
    if gdat.booldiagmode:
        if isinstance(gmodstat.lliktotl, np.ndarray):
            raise Exception('')
        if not np.isfinite(gmodstat.lliktotl).all():
            raise Exception('')

    numbdoff = gdat.numbdata - gmod.numbparagenrbase
    if gmod.numbparaelem > 0:
        for l in gmod.indxpopl:
            numbdoff -= len(gmodstat.indxparagenrfullelem[l]['full'])

    setattr(gmodstat, 'llik', llik)
    setattr(gmodstat, 'llikmean', gmodstat.lliktotl / gdat.numbdata) 
    setattr(gmodstat, 'llikcmea', gmodstat.lliktotl / (gdat.numbdata - numbdoff)) 

    if gdat.typeverb > 2:
        print('llik')
        summgene(llik)
    if gdat.typeverb > 1:
        print('gmodstat.lliktotl')
        print(gmodstat.lliktotl)
    stopchro(gdat, gdatmodi, 'llik')

    lpostotl = lpritotl + gmodstat.lliktotl
    if gdat.typeverb > 1:
        print('lpostotl')
        print(lpostotl)

    setattr(gmodstat, 'lpritotl', lpritotl) 
    setattr(gmodstat, 'gmodstat.lliktotl', gmodstat.lliktotl)
    setattr(gmodstat, 'lpostotl', lpostotl) 
    
    stopchro(gdat, gdatmodi, 'lpri')
    
    if strgstat == 'next':
        return

    initchro(gdat, gdatmodi, 'tert')
    
    setattr(gmodstat, 'lpri', lpri)
    
    if gmod.numbparaelem > 0:
        setattr(gmodstat, 'lpripena', lpri[0])
    
    dicttert = {}
    
    ## load necessary variables
        
    ## derived variables
    ## residual count map 
    cntp['resi'] = []
    cntp['resi'] = gdat.cntpdata - cntp['modl']
    
    setattr(gmodstat, 'cntpmodl', cntp['modl'])
    setattr(gmodstat, 'cntpresi', cntp['resi'])
    setattr(gmodstat, 'llik', llik)
    #if gmod.boollens:
    #    setattr(gmodstat, 'deflhost', deflhost)
    
    if gmod.boollens:
        
        setattr(gmodstat, 'defl', defl)
        for e in gmod.indxsersfgrd:
            masshostbein = massfrombein * beinhost[e]**2
            setattr(gmodstat, 'masshostisf%dbein' % e, masshostbein)
        ### sort with respect to deflection at scale radius
        if gmod.numbparaelem > 0:
            for l in gmod.indxpopl:
                if gmodstat.numbelem[l] > 0:
                    indxelemsortampl = np.argsort(gmodstat.dictelem[l][nameparaelemsort[l]])[::-1]
                    for nameparagenrelem in gmod.namepara.genrelem[l]:
                        gmodstat.dictelem[l][nameparagenrelem + 'sort'] = gmodstat.dictelem[l][nameparagenrelem][indxelemsortampl]

        deflsing = np.zeros((gdat.numbpixlcart, 2, numbdeflsingplot))
        conv = np.zeros((gdat.numbpixlcart))
        convpsec = np.zeros(((gdat.numbsidecarthalf)**2))
        convpsecodim = np.zeros((gdat.numbsidecarthalf))
        if gmod.numbparaelem > 0:
            if boolelemlens:
                gmod.indxpopllens = gmod.typeelem.index('lens')
        numbdeflsing = 2
        if gmod.numbparaelem > 0:
            if boolelemlens:
                if numbelem[indxpopllens] > 0:
                    numbdeflsing += min(numbdeflsubhplot, numbelem[indxpopllens]) 
                    numbdeflsing += 1
                for k in range(numbdeflsing):
                    indxpixltemp = gdat.indxpixlcart
                    if k == 0:
                        # temp -- should take other sersics into account
                        deflsing[indxpixltemp, :, k] = deflhost[0]
                    elif k == 1:
                        deflsing[indxpixltemp, :, k] = deflextr
                    elif k == 2:
                        deflsing[indxpixltemp, :, k] = defl - deflextr - deflhost[0]
                    else:
                        asca = gmodstat.dictelem[indxpopllens]['ascasort'][None, k-3]
                        acut = gmodstat.dictelem[indxpopllens]['acutsort'][None, k-3]
                        deflsing[listindxpixlelem[indxpopllens][k], :, k] = retr_defl(gdat, listindxpixlelem[indxpopllens][k], \
                                                      gmodstat.dictelem[indxpopllens]['lgalsort'][None, k-3], gmodstat.dictelem[indxpopllens]['bgalsort'][None, k-3], \
                                                      gmodstat.dictelem[indxpopllens]['defssort'][None, k-3], asca=asca, acut=acut)

        # convergence
        ## total
        conv[:] = retr_conv(gdat, defl) 
        convhost = np.zeros((gmod.numbsersfgrd, gdat.numbpixlcart))
        for e in gmod.indxsersfgrd:
            convhost[e, :] = retr_conv(gdat, deflhost[e]) 
        
        ### power spectrum
        #### two dimensional
        convpsec[:] = retr_psec(gdat, conv[:])
        
        #### one dimensional
        convpsecodim[:] = retr_psecodim(gdat, convpsec[:]) 
        setattr(gmodstat, 'convpsec', convpsec)
        setattr(gmodstat, 'convpsecodim', convpsecodim)
        setattr(gmodstat, 'conv', conv[...])
        for e in gmod.indxsersfgrd:
            setattr(gmodstat, 'convisf%d' % e, convhost[e, ...])
        
        ## subhalos
        if gmod.numbparaelem > 0:
            if boolelemlens:
                convelem = np.zeros((gdat.numbpixl))
                convpsecelem = np.zeros(((gdat.numbsidecarthalf)**2))
                convpsecelemodim = np.zeros((gdat.numbsidecarthalf))
                ### convergence
                convelem[:] = retr_conv(gdat, deflsubh) 
                ###  power spectrum
                ##### two dimensional
                convpsecelem[:] = retr_psec(gdat, convelem[:])
                ##### one dimensional
                convpsecelemodim[:] = retr_psecodim(gdat, convpsecelem[:]) 
                setattr(gmodstat, 'convpsecelem', convpsecelem)
                setattr(gmodstat, 'convpsecelemodim', convpsecelemodim)
                setattr(gmodstat, 'convelem', convelem[...])
                setattr(gmodstat, 'defl', defl)
        
        ### magnification
        magn = np.empty((gdat.numbpixlcart))
        histdefl = np.empty((gdat.numbdefl))
        if gmod.numbparaelem > 0 and boolelemlens:
            histdeflsubh = np.empty((gdat.numbdefl))
        deflsingmgtd = np.zeros((gdat.numbpixlcart, numbdeflsingplot))
        magn[:] = 1. / retr_invm(gdat, defl) 
        histdefl[:] = np.histogram(defl, bins=gdat.binspara.defl)[0]
        if gmod.numbparaelem > 0:
            if boolelemlens:
                histdeflsubh[:] = np.histogram(deflsubh, bins=gdat.binspara.deflsubh)[0]
        deflsingmgtd[:, :] = np.sqrt(np.sum(deflsing[...]**2, axis=1))
        if gmod.numbparaelem > 0:
            if boolelemlens:
                setattr(gmodstat, 'histdeflsubh', histdeflsubh)
        setattr(gmodstat, 'histdefl', histdefl)
        setattr(gmodstat, 'magn', magn[...])
        setattr(gmodstat, 'deflsing', deflsing[...])
        setattr(gmodstat, 'deflsingmgtd', deflsingmgtd[...])
    
    ## element related
    if gmod.numbparaelem > 0:
        if gdat.numbpixl == 1:
            for l in gmod.indxpopl:
                for k in range(gmodstat.numbelem[l]):
                    setattr(gmodstat, 'speclinepop%d%04d' % (l, k), gmodstat.dictelem[l]['spec'][:, k])
        
        if gdat.typedata == 'mock' and strgmodl == 'true' and gdat.numbpixl > 1:
            gdat.refrlgal = [[] for l in gmod.indxpopl]
            gdat.refrbgal = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                gdat.refrlgal[l] = np.tile(gmodstat.dictelem[l]['lgal'], [3] + list(np.ones(gmodstat.dictelem[l]['lgal'].ndim, dtype=int)))
                gdat.refrbgal[l] = np.tile(gmodstat.dictelem[l]['bgal'], [3] + list(np.ones(gmodstat.dictelem[l]['bgal'].ndim, dtype=int)))
    
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lghtpntspuls':
                gmodstat.dictelem[l]['per1'] = retr_per1(gmodstat.dictelem[l]['per0'], gmodstat.dictelem[l]['magf'])
        
    if gmod.numbparaelem > 0:
        if strgstat == 'this' or gdat.boolrefeforc and strgmodl == 'fitt':
            # correlate the fitting model elements with the reference elements
            if gdat.boolinforefr and not (strgmodl == 'true' and gdat.typedata == 'mock') and gdat.boolasscrefr:
                indxelemrefrasschits = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
                indxelemfittasschits = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gmod.indxpopl:
                        if gdat.refr.numbelem[q] == 0:
                            continue
                        
                        indxelemfittmatr = np.empty((gdat.refr.numbelem[q], gmodstat.numbelem[l]), dtype=int)
                        indxelemrefrmatr = np.empty((gdat.refr.numbelem[q], gmodstat.numbelem[l]), dtype=int)
                        matrdist = np.empty((gdat.refr.numbelem[q], gmodstat.numbelem[l]))
                        for k in range(gmodstat.numbelem[l]):
                            # construct a matrix of angular distances between reference and fitting elements
                            if gmod.typeelem[l].startswith('lghtline'):
                                matrdist[:, k] = abs(gdat.refrelin[q][0, :] - gmodstat.dictelem[l]['elin'][k]) / gdat.refrelin[q][0, :]
                            else:
                                matrdist[:, k] = retr_angldist(gdat, gdat.refr.dictelem[q]['lgal'][0, :], gdat.refr.dictelem[q]['bgal'][0, :], gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k])
                            indxelemrefrmatr[:, k] = np.arange(gdat.refr.numbelem[q])
                            indxelemfittmatr[:, k] = k
                        matrdist = matrdist.flatten()
                        indxelemrefrmatr = indxelemrefrmatr.flatten()
                        indxelemfittmatr = indxelemfittmatr.flatten()

                        # take only angular separations smaller than some threshold
                        indxmatrthrs = np.where(matrdist < gdat.anglassc)
                        matrdist = matrdist[indxmatrthrs]
                        indxelemrefrmatr = indxelemrefrmatr[indxmatrthrs]
                        indxelemfittmatr = indxelemfittmatr[indxmatrthrs]

                        # sort the remaining associations with respect to distance
                        indxmatrsort = np.argsort(matrdist)
                        matrdist = matrdist[indxmatrsort]
                        indxelemrefrmatr = indxelemrefrmatr[indxmatrsort]
                        indxelemfittmatr = indxelemfittmatr[indxmatrsort]
                        
                        for c in range(matrdist.size):
                            if indxelemrefrmatr[c] in indxelemrefrasschits[q][l] or indxelemfittmatr[c] in indxelemfittasschits[q][l]:
                                continue
                            indxelemrefrasschits[q][l].append(indxelemrefrmatr[c])
                            indxelemfittasschits[q][l].append(indxelemfittmatr[c])
                        
                        indxelemrefrasschits[q][l] = np.array(indxelemrefrasschits[q][l])
                        indxelemfittasschits[q][l] = np.array(indxelemfittasschits[q][l])
                setattr(gmodstat, 'indxelemrefrasschits', indxelemrefrasschits)
                setattr(gmodstat, 'indxelemfittasschits', indxelemfittasschits)
                
                indxelemrefrasscmiss = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
                indxelemfittasscfals = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gmod.indxpopl:
                        # indices of the reference elements not associated with the fitting model elements
                        if gdat.refr.numbelem[q] > 0:
                            indxelemrefrasscmiss[q][l] = np.setdiff1d(np.arange(gdat.refr.numbelem[q]), indxelemrefrasschits[q][l])
                        # indices of the fitting model elements not associated with the reference elements
                        if gmodstat.numbelem[l] > 0:
                            indxelemfittasscfals[q][l] = np.setdiff1d(np.arange(gmodstat.numbelem[l]), indxelemfittasschits[q][l])
                setattr(gmodstat, 'indxelemrefrasscmiss', indxelemrefrasscmiss)
                setattr(gmodstat, 'indxelemfittasscfals', indxelemfittasscfals)
                
                for q in gdat.indxrefr:
                    if gdat.refr.numbelem[q] == 0:
                        continue
                    for l in gmod.indxpopl:
                        # collect the associated reference element parameter for each fitting element 
                        for strgfeat in gdat.refr.namepara.elemonly[q][l]:
                            name = strgfeat + gdat.listnamerefr[q]
                            if strgfeat != 'spec' and strgfeat != 'specplot':
                                refrfeat = getattr(gdat.refr, strgfeat)
                                gmodstat.dictelem[l][name] = np.zeros(gmodstat.numbelem[l])
                                if len(refrfeat[q]) > 0 and len(indxelemrefrasschits[q][l]) > 0:
                                    gmodstat.dictelem[l][name][indxelemfittasschits[q][l]] = refrfeat[q][0, indxelemrefrasschits[q][l]]
                        
                        print('temp')
                        continue

                        # collect the error in the associated reference element amplitude
                        for strgfeat in gdat.listnameparaetotlelemcomm[q][l]:
                            refrfeat = getattr(gdat.refr, strgfeat)
                            if strgfeat == gmod.nameparagenrelemampl[l] and len(indxelemfittasschits[q][l]) > 0:
                                gmodstat.dictelem[l]['aerr' + gdat.listnamerefr[q]] = np.zeros(gmodstat.numbelem[l])
                                fittfeattemp = gmodstat.dictelem[l][strgfeat][indxelemfittasschits[q][l]]
                                refrfeattemp = refrfeat[q][0, indxelemrefrasschits[q][l]]
                                if gdat.booldiagmode:
                                    if not np.isfinite(refrfeattemp).all():
                                        raise Exception('')
                                gmodstat.dictelem[l]['aerr' + gdat.listnamerefr[q]][indxelemfittasschits[q][l]] = 100. * (fittfeattemp - refrfeattemp) / refrfeattemp
                
            if gdat.boolrefeforc and strgmodl == 'fitt':
                for l in gmod.indxpopl:
                    for strgfeat in gmod.namepara.genrelem[l]:
                        if strgfeat in gdat.refr.namepara.elem[gdat.indxrefrforc[l]]:
                            if len(indxelemrefrasschits[gdat.indxrefrforc[l]][l]) == 0:
                                continue
                            refrfeat = getattr(gdat.refr, strgfeat)[gdat.indxrefrforc[l]][0, indxelemrefrasschits[gdat.indxrefrforc[l]][l]]
                            if len(gmodstat.dictelem[l][strgfeat]) == 0:
                                continue
                            lpritotl += -2. * np.sum(1e6 * (gmodstat.dictelem[l][strgfeat][indxelemfittasschits[gdat.indxrefrforc[l]][l]] - refrfeat)**2 / refrfeat**2)

    # other tertiary variables continues
    ## number of degrees of freedom

    chi2doff = np.sum(cntp['resi']**2 / gdat.varidata) / numbdoff
    if gdat.booldiagmode:
        if not np.isfinite(cntp['resi']).all():
            raise Exception('')
        if not np.isfinite(numbdoff):
            raise Exception('')
        if not np.isfinite(chi2doff):
            raise Exception('')
    setattr(gmodstat, 'numbdoff', numbdoff)
    setattr(gmodstat, 'chi2doff', chi2doff)
    
    if gmod.boolelempsfn and gmod.numbparaelem > 0:
        gmodstat.fwhmpsfn = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            
    if gmod.numbparaelem > 0:
        
        ### derived parameters
        for l in gmod.indxpopl:

            # luminosity
            if gmod.boolelemlght[l] and 'flux' in gmod.namepara.genrelem[l]:
                for strgfeat in gmod.namepara.genrelem[l]:
                    if strgfeat.startswith('reds') and strgfeat != 'reds':
                        namerefr = strgfeat[-4:]
                        gmodstat.dictelem[l]['lumi' + namerefr] = np.zeros(gmodstat.numbelem[l]) + np.nan
                        gmodstat.dictelem[l]['dlos' + namerefr] = np.zeros(gmodstat.numbelem[l]) + np.nan
                        reds = gmodstat.dictelem[l]['reds' + namerefr]
                        indxgood = np.where(np.isfinite(gmodstat.dictelem[l]['reds' + namerefr]))[0]
                        if indxgood.size > 0:
                            # temp -- these units only work for energy units of keV
                            dlos = gdat.adisobjt(reds)
                            gmodstat.dictelem[l]['dlos' + namerefr][indxgood] = dlos
                            lumi = retr_lumi(gdat, gmodstat.dictelem[l]['flux'], dlos, reds)
                            gmodstat.dictelem[l]['lumi' + namerefr][indxgood] = lumi
        
            if gmod.typeelem[l] == 'lghtpntsagnntrue':
                gmodstat.dictelem[l]['reds'] = gdat.redsfromdlosobjt(gmodstat.dictelem[l]['dlos'])
            if gmod.typeelem[l] == 'lghtpntspuls':
                gmodstat.dictelem[l]['mass'] = full([numbelem[l]], 3.)

            if gdat.typeverb > 2:
                print('l')
                print(l)
            if gdat.boolbinsspat:
                #### radial and angular coordinates
                gmodstat.dictelem[l]['gang'] = retr_gang(gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal'])
                gmodstat.dictelem[l]['aang'] = retr_aang(gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal'])
            
            if gmod.boolelemlght[l]:
                #### number of expected counts
                if gdat.boolbinsspat:
                    gmodstat.dictelem[l]['cnts'] = retr_cntspnts(gdat, [gmodstat.dictelem[l]['lgal'], gmodstat.dictelem[l]['bgal']], gmodstat.dictelem[l]['spec'])
                else:
                    gmodstat.dictelem[l]['cnts'] = retr_cntspnts(gdat, [gmodstat.dictelem[l]['elin']], gmodstat.dictelem[l]['spec'])
            
            #### delta log-likelihood
            gmodstat.dictelem[l]['deltllik'] = np.zeros(gmodstat.numbelem[l])
            if not (strgmodl == 'true' and gdat.checprio): 
                if gdat.typeverb > 2:
                    print('Calculating log-likelihood differences when removing elements from the model.')
                for k in range(gmodstat.numbelem[l]):
                    
                    # construct gdatmodi
                    gdatmoditemp = tdpy.gdatstrt()
                    gdatmoditemp.this = tdpy.gdatstrt()
                    gdatmoditemp.next = tdpy.gdatstrt()
                    gdatmoditemp.this.indxelemfull = gmodstat.indxelemfull
                    gdatmoditemp.this.paragenrscalfull = gmodstat.paragenrscalfull
                    gdatmoditemp.this.paragenrunitfull = gmodstat.paragenrunitfull

                    prop_stat(gdat, gdatmoditemp, strgmodl, deth=True, thisindxpopl=l, thisindxelem=k)
                    proc_samp(gdat, gdatmoditemp, 'next', strgmodl)#, boolinit=boolinit)
                    
                    if gdat.booldiagmode:
                        if not np.isfinite(gmodstat.lliktotl):
                            raise Exception('')
                    
                    gdatobjttemp = retr_gdatobjt(gdat, gdatmoditemp, strgmodl)#, boolinit=boolinit)
                    nextlliktotl = gdatobjttemp.next.lliktotl
                    gmodstat.dictelem[l]['deltllik'][k] = gmodstat.lliktotl - nextlliktotl
                    
                if gdat.typeverb > 2:
                    print('deltllik calculation ended.')
    
    # more derived parameters
    if (gmod.typeevalpsfn == 'kern' or gmod.typeevalpsfn == 'full') and (strgmodl == 'true' or boolinit or gdat.boolmodipsfn):
        ### PSF FWHM
        if gdat.typepixl == 'cart':
            fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
        setattr(gmodstat, 'fwhm', fwhm)
    
    if gmod.numbparaelem > 0 and gmod.boolelemsbrtdfncanyy:
        
        if gmod.numbparaelem > 0:
            sbrt['dfnctotl'] = np.zeros_like(gdat.expo)
            sbrt['dfncsubt'] = np.zeros_like(gdat.expo)
            sbrt['dfncsupt'] = np.zeros_like(gdat.expo)
            for l in gmod.indxpopl:
                if gmod.boolcalcerrr[l]:
                    sbrt['dfncfull'] = np.zeros_like(gdat.expo)
                if gmod.boolelemsbrt[l]:
                    for k in range(gmodstat.numbelem[l]):
                        
                        # read normalization from the element dictionary
                        if gmod.boolelemlght[l]:
                            varbamplextd = gmodstat.dictelem[l]['spec'][:, k]
                        if gmod.typeelem[l].startswith('clus'):
                            varbamplextd = gmodstat.dictelem[l]['nobj'][None, k]
                        
                        # calculate imprint on the element surface brightness state variable
                        if gmod.boolelempsfn[l]:
                            sbrttemp = retr_sbrtpnts(gdat, gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                                                                    varbamplextd, gmodstat.psfnintp, listindxpixlelem[l][k])
                        indxpixltemp = listindxpixlelem[l][k]

                        if gmod.typeelem[l].startswith('lghtline'):
                            sbrttemp = gmodstat.dictelem[l]['spec'][:, k, None, None]
                        
                        # add it to the state variable depending on the significance
                        sbrt['dfnctotl'][:, indxpixltemp, :] += sbrttemp
                        if gmodstat.dictelem[l]['deltllik'][k] > 35:
                            sbrt['dfncsupt'][:, indxpixltemp, :] += sbrttemp
                        if gmodstat.dictelem[l]['deltllik'][k] < 35:
                            sbrt['dfncsubt'][:, indxpixltemp, :] += sbrttemp
                        
                        # calculate imprint without PSF truncation to calculate approximation errors
                        if gmod.boolcalcerrr[l]:
                            sbrt['dfncfull'][:, :, :] += retr_sbrtpnts(gdat, gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                                                                            varbamplextd, gmodstat.psfnintp, gdat.indxpixl)
            
                setattr(gmodstat, 'sbrtdfncsubtpop%d' % l, sbrt['dfncsubt'])
                
    if gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
        if gdat.booldiagmode:
            numbtemp = 0
            for l in gmod.indxpopl:
                if gmod.boolelemsbrtextsbgrd[l]:
                    numbtemp += np.sum(gmodstat.numbelem[l])
            if numbtemp > 0 and (sbrtextsbgrd == 0.).all():
                raise Exception('')

        sbrt['bgrdexts'] = sbrtextsbgrd

    #### count maps
    cntp = dict()
    for name in gmod.listnamegcom:
        cntp[name] = retr_cntp(gdat, sbrt[name])
        setattr(gmodstat, 'cntp' + name, cntp[name])
    
    ### spatial averages
    sbrtmean = dict()
    sbrtstdv = dict()
    for name in gmod.listnamegcom:
        sbrtmean[name], sbrtstdv[name] = retr_spatmean(gdat, sbrt[name])
        for b in gdat.indxspatmean:
            setattr(gmodstat, 'sbrt%smea%d' % (name, b), sbrtmean[name][b])
            setattr(gmodstat, 'sbrt%sstd%d' % (name, b), sbrtstdv[name][b])
    
    if gmod.numbparaelem > 0:
        if gmod.boolelemsbrtdfncanyy:
            for i in gdat.indxener:
                if 'dark' in gmod.listnamegcom:
                    fracsdenmeandarkdfncsubt = sbrtmean['dfncsubt'][0][0][i] / (sbrtmean['dfncsubt'][0][0][i] + sbrtmean['dark'][0][0][i])
                else:
                    fracsdenmeandarkdfncsubt = 1.
                setattr(gmodstat, 'fracsdenmeandarkdfncsubten%02d' % i, np.array([fracsdenmeandarkdfncsubt]))
            
            if 'dark' in gmod.listnamegcom:
                booldfncsubt = float(np.where(sbrtmean['dfncsubt'][0][0] > sbrtmean['dark'][0][0])[0].any())
            else:
                booldfncsubt = 1.
            setattr(gmodstat, 'booldfncsubt', np.array([booldfncsubt]))

    # find the 1-point function of the count maps of all emission components including the total emission
    for name in gmod.listnamegcom:
        namehistcntp = 'histcntp' + name
        for m in gdat.indxevtt:
            if gdat.numbevtt > 1:
                namehistcntp += 'evt%d' % m
            for i in gdat.indxener: 
                if gdat.numbener > 1:
                    namehistcntp += 'en%02d' % i
                
                histcntp = np.histogram(cntp[name][i, :, m], bins=gdat.binspara.cntpmodl)[0]
                setattr(gmodstat, namehistcntp, histcntp)
                
                if False and i == 0 and m == 0 and (name == 'dfnc' or name == 'dfncsubt'):
                    for strgbins in ['lowr', 'higr']:
                        strgtemp = 'histcntp' + strgbins + name + 'en%02devt%d' % (i, m)
                        if strgbins == 'lowr':
                            setattr(gmod, strgtemp, np.array([float(np.sum(histcntp[:gdat.numbtickcbar-1]))]))
                        else:
                            setattr(gmod, strgtemp, np.array([float(np.sum(histcntp[gdat.numbtickcbar-1:]))]))
            else:
                histcntp = np.histogram(cntp[name][:, 0, m], bins=gdat.binspara.cntpmodl)[0]
                setattr(gmodstat, 'histcntp' + name + 'evt%d' % m, histcntp)

    if gmod.boollens:
        if strgmodl == 'true':
            s2nr = []
            s2nr = cntp['lens'] / np.sqrt(cntp['modl'])
            setattr(gmodstat, 's2nr', s2nr)
        cntplensgrad = np.empty((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt, 2))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                cntplenstemp = np.zeros(gdat.numbpixlcart)
                cntplenstemp[gdat.indxpixlrofi] = cntp['lens'][i, :, m]
                cntplensgrad[i, :, m, :] = retr_gradmaps(gdat, cntplenstemp) * gdat.sizepixl
        
        cntplensgradmgtd = np.sqrt(np.sum(cntplensgrad**2, axis=3))
        cntplensgrad *= gdat.sizepixl
        indx = np.where(np.fabs(cntplensgrad) > 1. * gdat.sizepixl)
        cntplensgrad[indx] = np.sign(cntplensgrad[indx]) * 1. * gdat.sizepixl
        deflmgtd = np.sqrt(np.sum(defl**2, axis=1))
        setattr(gmodstat, 'deflmgtd', deflmgtd)
        setattr(gmodstat, 'cntplensgrad', cntplensgrad)
        setattr(gmodstat, 'cntplensgradmgtd', cntplensgradmgtd)

    if gmod.numbparaelem > 0:
        for l in gmod.indxpopl:
            if gmod.boolelemlght[l]:
                #### spectra
                if gdat.boolbinsspat:
                    sindcolr = [gmodstat.dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    gmodstat.dictelem[l]['specplot'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], sind=gmodstat.dictelem[l]['sind'], \
                                                                 curv=gmodstat.dictelem[l]['curv'], expc=gmodstat.dictelem[l]['expc'], \
                                                                 sindcolr=sindcolr, spectype=gmod.spectype[l], plot=True)
                
                if gdat.typedata == 'inpt':
                    if gdat.typeexpr == 'ferm':
                        # temp
                        try:
                            gmodstat.dictelem[l]['sbrt0018'] = gdat.sbrt0018objt(gmodstat.dictelem[l]['bgal'], gmodstat.dictelem[l]['lgal'])
                        except:
                            gmodstat.dictelem[l]['sbrt0018'] = gmodstat.dictelem[l]['bgal'] * 0.

            if gmod.typeelem[l] == 'lens':
                #### distance to the source
                if gmod.boollens:
                    gmodstat.dictelem[l]['diss'] = retr_angldist(gdat, gmodstat.dictelem[l]['lgal'],  gmodstat.dictelem[l]['bgal'], lgalsour, bgalsour)
                
                if gmod.boollenssubh:
                    gmodstat.dictelem[l]['deflprof'] = np.empty((gdat.numbanglfull, gmodstat.numbelem[l]))
                    gmodstat.dictelem[l]['mcut'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['rele'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['reln'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relk'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relf'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['reld'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relc'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relm'] = np.empty(gmodstat.numbelem[l])

                    # temp -- this can be placed earlier in the code
                    cntplensobjt = sp.interpolate.RectBivariateSpline(gdat.meanpara.bgalcart, gdat.meanpara.lgalcart, \
                                                            cntp['lens'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    for k in np.arange(gmodstat.numbelem[l]):
                        
                        asca = gmodstat.dictelem[l]['asca'][k]
                        acut = gmodstat.dictelem[l]['acut'][k]
                        
                        #### deflection profiles
                        gmodstat.dictelem[l]['deflprof'][:, k] = retr_deflcutf(gdat.meanpara.anglfull, gmodstat.dictelem[l]['defs'][k], asca, acut)
         
                        ### truncated mass 
                        gmodstat.dictelem[l]['mcut'][k] = retr_mcut(gdat, gmodstat.dictelem[l]['defs'][k], asca, acut, adishost, mdencrit)

                        #### dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        gmodstat.dictelem[l]['rele'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                                              gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        
                        gmodstat.dictelem[l]['relf'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                                         gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, cntpmodl=cntp['modl'][0, :, 0])
                        
                        deflelem = retr_defl(gdat, gdat.indxpixl, gmodstat.dictelem[l]['lgal'][k], \
                                                                            gmodstat.dictelem[l]['bgal'][k], gmodstat.dictelem[l]['defs'][k], asca=asca, acut=acut)
                        bgalprim = gdat.bgalgrid - deflelem[:, 1]
                        lgalprim = gdat.lgalgrid - deflelem[:, 0]
                        gmodstat.dictelem[l]['relm'][k] = np.mean(abs(cntp['lens'][0, :, 0] - cntplensobjt(bgalprim, lgalprim, grid=False).flatten()))
                        
                        
                        gmodstat.dictelem[l]['relk'][k] = gmodstat.dictelem[l]['relm'][k] / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
                        gmodstat.dictelem[l]['reln'][k] = gmodstat.dictelem[l]['rele'][k] / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
                        gmodstat.dictelem[l]['reld'][k] = retr_rele(gdat, gdat.cntpdata[0, :, 0], gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                                                                         gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        gmodstat.dictelem[l]['relc'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['lgal'][k], gmodstat.dictelem[l]['bgal'][k], \
                                                       gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, absv=False) / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
               
        ### distribution of element parameters and features
        #### calculate the model filter
        listindxelemfilt = [[[] for l in gmod.indxpopl] for namefilt in gdat.listnamefilt]
        for k, namefilt in enumerate(gdat.listnamefilt):
            for l in gmod.indxpopl:
                if namefilt == '':
                    listindxelemfilt[k][l] = np.arange(gmodstat.numbelem[l])
                if namefilt == 'imagbndr':
                    listindxelemfilt[k][l] = np.where((np.fabs(gmodstat.dictelem[l]['lgal']) < gdat.maxmgangdata) & (np.fabs(gmodstat.dictelem[l]['bgal']) < gdat.maxmgangdata))[0]
                if namefilt == 'deltllik':
                    listindxelemfilt[k][l] = np.where(gmodstat.dictelem[l]['deltllik'] > 0.5 * gmod.numbparagenrelemsing[l])[0]
                if namefilt == 'nrel':
                    listindxelemfilt[k][l] = np.where(gmodstat.dictelem[l]['reln'] > 0.3)[0]
    
        for l in gmod.indxpopl:
            # histograms of element parameters
            for namefrst in gmod.namepara.elem[l]:
               
                ## one dimensional
                if namefrst[:-4] == 'etag':
                    continue
                if namefrst == 'specplot' or namefrst == 'deflprof':
                    continue
                elif namefrst == 'spec':
                    histfrst = np.zeros((gdat.numbbinsplot, gdat.numbener))
                    for i in gdat.indxener:
                        histfrst[:, i] = np.histogram(gmodstat.dictelem[l]['spec'][i, listindxelemfilt[0][l]], gdat.binspara.spec)[0]
                elif namefrst == 'cnts':
                    histfrst = np.histogram(gmodstat.dictelem[l]['cnts'][listindxelemfilt[0][l]], gdat.binspara.cnts)[0]
                else:
                #elif not (namefrst == 'curv' and gmod.spectype[l] != 'curv' or namefrst == 'expc' \
                #                                        and gmod.spectype[l] != 'expc' or namefrst.startswith('sindarry') and \
                #                                                                                          gmod.spectype[l] != 'colr'):
                    binsfrst = getattr(gdat.binspara, namefrst)
                    #if len(gmodstat.dictelem[l][namefrst]) > 0 and len(listindxelemfilt[0][l]) > 0:
                    histfrst = np.histogram(gmodstat.dictelem[l][namefrst][listindxelemfilt[0][l]], binsfrst)[0]
                    strgvarb = 'hist' + namefrst + 'pop%d' % l
                    setattr(gmodstat, strgvarb, histfrst)
                        
                #### two dimensional
                for nameseco in gmod.namepara.elem[l]:
                    if namefrst == 'spec' or namefrst == 'specplot' or namefrst == 'deflprof' or \
                            nameseco == 'spec' or nameseco == 'specplot' or nameseco == 'deflprof':
                        continue
                    
                    if not checstrgfeat(namefrst, nameseco):
                        continue

                    binsseco = getattr(gdat.binspara, nameseco)
                    histtdim = np.histogram2d(gmodstat.dictelem[l][namefrst][listindxelemfilt[0][l]], \
                                                            gmodstat.dictelem[l][nameseco][listindxelemfilt[0][l]], [binsfrst, binsseco])[0]
            
                    setattr(gmodstat, 'hist' + namefrst + nameseco + 'pop%d' % l, histtdim)
                
            ### priors on element parameters and features
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                xdat = gmodstat.dictelem[l][nameparagenrelem]
                minm = getattr(gmod.minmpara, nameparagenrelem + 'pop%d' % l)
                maxm = getattr(gmod.maxmpara, nameparagenrelem + 'pop%d' % l)
                scal = getattr(gmod.scalpara, nameparagenrelem + 'pop%d' % l)
                booltemp = False
                if scal.startswith('expo') or scal.startswith('dexp'):
                    if scal.startswith('expo'):
                        if scal == 'expo':
                            sexp = getattr(gmod, 'gangdistsexppop%d' % l)
                        else:
                            sexp = gmodstat.paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distscal')[l]]
                        pdfn = pdfn_expo(xdat, maxm, sexp)
                    if scal.startswith('dexp'):
                        pdfn = pdfn_dnp.exp(xdat, maxm, scal)
                    booltemp = True
                if scal.startswith('self') or scal.startswith('logt'):
                    if scal.startswith('self'):
                        pdfn = 1. / (maxm - minm) + np.zeros_like(xdat)
                    else:
                        pdfn = 1. / (np.log(maxm) - np.log(minm)) + np.zeros_like(xdat)
                    booltemp = True
                # temp 
                if scal.startswith('powr'):
                    slop = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'slopprio' + nameparagenrelem + 'pop%d' % l)]
                    pdfn = pdfn_powr(xdat, minm, maxm, slop)
                    booltemp = True
                if scal.startswith('dpowslopbrek'):
                    pdfn = pdfn_dpow(xdat, minm, maxm, brek, sloplowr, slopuppr)
                    booltemp = True
                if scal == 'lnormeanstdv':
                    pdfn = pdfn_lnor(xdat, meanlnor, stdvlnor)
                    booltemp = True
                if scal.startswith('igam'):
                    cutf = getattr(gdat, 'cutf' + nameparagenrelem)
                    pdfn = pdfn_igam(xdat, slop, cutf)
                    booltemp = True
                if scal.startswith('gaus'):
                    # this does not work for mismodeling
                    meanvarb = gmodstat.paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distmean')[l]]
                    stdv = gmodstat.paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'diststdv')[l]]
                    if nameparagenrelem == 'expc' and gmod.spectype[l] == 'expc':
                        pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                    else:
                        pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                    booltemp = True
                
                # temp -- meanelem will not be defined
                #if booltemp:
                #    gmodstat.dictelem[l]['hist' + nameparagenrelem + 'prio'] = gmodstat.numbelem[l] * pdfn * np.interp(xdat, xdatplot, delt)
                
                #setattr(gmodstat, 'hist' + nameparagenrelem + 'pop%dprio' % l, gmodstat.dictelem[l]['hist' + nameparagenrelem + 'prio'])
                #if strgmodl == 'true':
                #    setattr(gmodstat, 'refrhist' + nameparagenrelem + 'pop%dprio' % l, gmodstat.dictelem[l]['hist' + nameparagenrelem + 'prio'])
    
    if gmod.numbparaelem > 0:
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lens':
                if gmodstat.numbelem[l] > 0:
                    ## total truncated mass of the subhalo as a cross check
                    # temp -- generalize
                    asca = gmodstat.dictelem[l]['asca']
                    acut = gmodstat.dictelem[l]['acut']
                    factmcutfromdefs = retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut) 
                    masssubh = np.array([np.sum(factmcutfromdefs * gmodstat.dictelem[l]['defs'])])
    
    ## derived variables as a function of other derived variables
    if gmod.numbparaelem > 0:
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lghtpntspuls'):
                massshel = np.empty(gdat.numbanglhalf)
                for k in gdat.indxanglhalf:
                    indxelemshel = np.where((gdat.binspara.anglhalf[k] < gmodstat.dictelem[l]['gang']) & (gmodstat.dictelem[l]['gang'] < gdat.binspara.anglhalf[k+1]))
                    massshel[k] = np.sum(gmodstat.dictelem[l]['mass'][indxelemshel])
                setattr(gmodstat, 'massshelpop%d' % l, massshel)
            
    if gmod.boollens or gmod.numbparaelem > 0 and gmod.boollenssubh:
        # find the host, subhalo masses and subhalo mass fraction as a function of halo-centric radius
        listnametemp = gdat.liststrgcalcmasssubh
        listnamevarbmass = []
        listnamevarbmassscal = []
        listnamevarbmassvect = []
        for e in gmod.indxsersfgrd:
            if boolllenshost:
                listnamevarbmassscal += ['masshosttotl']
                for strgtemp in listnametemp:
                    listnamevarbmassvect.append('masshostisf%d' % e + strgtemp)
                    listnamevarbmassscal.append('masshostisf%d' % e + strgtemp + 'bein')
        if gmod.numbparaelem > 0 and gmod.boollenssubh:
            listnamevarbmassscal.append('masssubhtotl')
            listnamevarbmassscal.append('fracsubhtotl')
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masssubh' + strgtemp)
                listnamevarbmassvect.append('fracsubh' + strgtemp)
                listnamevarbmassscal.append('masssubh' + strgtemp + 'bein')
                listnamevarbmassscal.append('fracsubh' + strgtemp + 'bein')

        for name in listnamevarbmassvect:
            dicttert[name] = np.zeros(gdat.numbanglhalf)
            if 'isf' in name:
                indxisfrtemp = int(name.split('isf')[1][0])
            angl = np.sqrt((gdat.meanpara.lgalcartmesh - lgalhost[indxisfrtemp])**2 + (gdat.meanpara.bgalcartmesh - bgalhost[indxisfrtemp])**2).flatten()
            for k in gdat.indxanglhalf:
                if name[4:8] == 'host':
                    convtemp = conv[:]
                if name[4:8] == 'subh':
                    convtemp = convelem[:]
                
                if name.endswith('delt'):
                    indxpixl = np.where((gdat.binspara.anglhalf[k] < angl) & (angl < gdat.binspara.anglhalf[k+1]))[0]
                    dicttert[name][k] = 1e6 * np.sum(convtemp[indxpixl]) * mdencrit * \
                                                gdat.apix * adishost**2 / 2. / np.pi * gdat.deltanglhalf[k] / gdat.meanpara.anglhalf[k]
                if name.endswith('intg'):
                    indxpixl = np.where(angl < gdat.meanpara.anglhalf[k])[0]
                    dicttert[name][k] = np.sum(convtemp[indxpixl]) * mdencrit * gdat.apix * adishost**2
                
                if name[:4] == 'frac':
                    masshosttotl = 0.
                    for e in gmod.indxsersfgrd:
                        masshosttotl += dicttert['masshostisf%d' % e + name[-4:]][k]
                    if masshosttotl != 0.:
                        dicttert['fracsubh' + name[8:]][k] = dicttert['masssubh' + name[8:]][k] / masshosttotl
            setattr(gmodstat, name, dicttert[name])
            
            # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
            dicttert[name + 'bein'] = np.interp(beinhost, gdat.meanpara.anglhalf, dicttert[name])
            setattr(gmodstat, name + 'bein', dicttert[name + 'bein'])
        
    #if gmod.numbparaelem > 0:
    #    ## copy element parameters to the global object
    #    feat = [[] for l in gmod.indxpopl]
    #    for l in gmod.indxpopl:
    #        feat[l] = dict()
    #        for strgfeat in gmod.namepara.genrelem[l]:
    #            if strgfeat[:-4] == 'etag':
    #                continue
    #            if len(gmodstat.dictelem[l][strgfeat]) > 0:
    #                if strgmodl == 'true':
    #                    shap = list(np.ones(gmodstat.dictelem[l][strgfeat].ndim, dtype=int))
    #                    feat[l][strgfeat] = np.tile(gmodstat.dictelem[l][strgfeat], [3] + shap)
    #                if strgmodl == 'fitt':
    #                    feat[l][strgfeat] = gmodstat.dictelem[l][strgfeat]
    #                
    #    #for strgfeat in gmod.namepara.elem:
    #    #    feattemp = [[] for l in gmod.indxpopl]
    #    #    for l in gmod.indxpopl:
    #    #        if strgfeat in gmod.namepara.genrelem[l]:
    #    #            if strgfeat in feat[l]:
    #    #                feattemp[l] = feat[l][strgfeat]
    #    #            else:
    #    #                feattemp[l] = np.array([])
    #    #    setattr(gmodstat, strgfeat, feattemp)
        
    # copy true state to the reference state
    #if strgmodl == 'true':
    #    for name, valu in deepcopy(gdat.__dict__).items():
    #        if name.startswith('true'):
    #            #indx = name.find('pop')
    #            #if indx != -1 and not name.endswith('pop') and name[indx+3].isdigit():
    #            #    namerefr = name.replace('pop%s' % name[indx+3], 'ref%s' % name[indx+3])
    #            #else:
    #            #    namerefr = name
    #            #namerefr = name
    #            #namerefr = namerefr.replace('true', 'refr')
    #            name = name.replace('true', 'refr')
    #            setattr(gdat, name, valu)
    
    if gmod.numbparaelem > 0 and gdat.priofactdoff != 0.:
        if strgmodl == 'true':
            for q in gdat.indxrefr:
                for strgfeat in gdat.refr.namepara.elem[q]:
                    
                    if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':
                        continue
                    reca = np.zeros(gdat.numbbinsplot) - 1.
                    
                    indxelempars = np.where(gmodstat.dictelem[q]['deltllik'] > 2.5)[0]
                    
                    refrhistpars = np.zeros(gdat.numbbinsplot) - 1.
                    
                    histparaelem = getattr(gmodstat, 'hist' + strgfeat + 'pop%d' % q)
                    indxrefrgood = np.where(histparaelem > 0)[0]
                    reca[indxrefrgood] = 0.
                    refrhistpars[indxrefrgood] = 0.
                    refrhist = getattr(gmodstat, 'hist' + strgfeat + 'pop%d' % q)

                    bins = getattr(gdat.binspara, strgfeat)
                    if len(indxelempars) > 0:
                        refrhistpars = np.histogram(gmodstat.dictelem[q][strgfeat][indxelempars], bins=bins)[0].astype(float)
                        if indxrefrgood.size > 0:
                            reca[indxrefrgood] = refrhistpars[indxrefrgood] / refrhist[indxrefrgood]
                    
                    setattr(gmodstat, 'histpars' + strgfeat + 'pop%d' % q, refrhistpars)
                    setattr(gmodstat, 'reca' + strgfeat + 'pop%d' % q, reca)
        
        print('gdat.rtagmock')
        print(gdat.rtagmock)
        if gdat.rtagmock is not None:
            if gmod.numbparaelem > 0:
                for l in gmod.indxpopl:
                    for strgfeat in gmod.namepara.genrelem[l]:
                        if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':# or strgfeat.startswith('aerr'):
                            continue
                        if strgfeat in gmod.namepara.genrelem[l]:
                            hist = getattr(gmodstat, 'hist' + strgfeat + 'pop%d' % l)
                            reca = getattr(gdat.true.this, 'reca' + strgfeat + 'pop%d' % l)
                            histcorrreca = hist / reca
                            setattr(gmodstat, 'histcorrreca' + strgfeat + 'pop%d' % l, histcorrreca)

    ### Exculusive comparison with the true state
    if strgmodl == 'fitt' and gdat.typedata == 'mock':
        if gmod.boollens:
            numbsingcomm = min(deflsing.shape[2], gmod.deflsing.shape[2])
            deflsingresi = deflsing[0, ..., :numbsingcomm] - gmod.deflsing[..., :numbsingcomm]
            deflsingresimgtd = np.sqrt(np.sum(deflsingresi**2, axis=1))
            deflsingresiperc = 100. * deflsingresimgtd / gmod.deflsingmgtd[..., :numbsingcomm]
            setattr(gmodstat, 'numbsingcomm', numbsingcomm)
            setattr(gmodstat, 'deflsingresi', deflsingresi)
            truedeflmgtd = getattr(gdat.true.this, 'deflmgtd')
            truedefl = getattr(gdat.true.this, 'defl')
            deflresi = defl - truedefl
            deflresimgtd = np.sqrt(np.sum(deflresi**2, axis=1))
            deflresiperc = 100. * deflresimgtd / truedeflmgtd
            setattr(gmodstat, 'deflresi', deflresi)
            setattr(gmodstat, 'deflresimgtd', deflresimgtd)
            if gmod.numbparaelem > 0:
                trueconvelem = getattr(gdat.true.this, 'convelem')
                convelemresi = convelem[:] - trueconvelem
                convelemresiperc = 100. * convelemresi / trueconvelem
                setattr(gmodstat, 'convelemresi', convelemresi)
                setattr(gmodstat, 'convelemresiperc', convelemresiperc)
            truemagn = getattr(gdat.true.this, 'magn')
            magnresi = magn[:] - truemagn
            magnresiperc = 100. * magnresi / truemagn
            setattr(gmodstat, 'magnresi', magnresi)
            setattr(gmodstat, 'magnresiperc', magnresiperc)
    
    if gmod.numbparaelem > 0:
    
        # correlate the catalog sample with the reference catalog
        if gdat.boolinforefr and not (strgmodl == 'true' and gdat.typedata == 'mock') and gdat.boolasscrefr:
            
            for q in gdat.indxrefr:
                for l in gmod.indxpopl:
                    if gdat.refr.numbelem[q] > 0:
                        cmpl = np.array([float(len(indxelemrefrasschits[q][l])) / gdat.refr.numbelem[q]])
                        if gdat.booldiagmode:
                            if cmpl > 1. or cmpl < 0.:
                                raise Exception('')
                    else:
                        cmpl = np.array([-1.])
                    setattr(gmodstat, 'cmplpop%dpop%d' % (l, q), cmpl)
                    if gmodstat.numbelem[l] > 0:
                        fdis = np.array([float(indxelemfittasscfals[q][l].size) / gmodstat.numbelem[l]])
                        if gdat.booldiagmode:
                            if fdis > 1. or fdis < 0.:
                                raise Exception('')
                    else:
                        fdis = np.array([-1.])
                    setattr(gmodstat, 'fdispop%dpop%d' % (q, l), fdis)
                        
            # collect the associated fitting element parameter for each reference element
            featrefrassc = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
            for q in gdat.indxrefr:
                for l in gmod.indxpopl:
                    featrefrassc[q][l] = dict()
                    for strgfeat in gdat.refr.namepara.elem[q]:
                        if not strgfeat in gmod.namepara.genrelem[l] or strgfeat in gdat.refr.namepara.elemonly[q][l]:
                            continue
                        if isinstance(gmodstat.dictelem[l][strgfeat], np.ndarray) and gmodstat.dictelem[l][strgfeat].ndim > 1:
                            continue
                        featrefrassc[q][l][strgfeat] = np.zeros(gdat.refr.numbelem[q]) + np.nan
                        if len(indxelemrefrasschits[q][l]) > 0 and len(gmodstat.dictelem[l][strgfeat]) > 0:
                            featrefrassc[q][l][strgfeat][indxelemrefrasschits[q][l]] = gmodstat.dictelem[l][strgfeat][indxelemfittasschits[q][l]]
                        name = strgfeat + 'asscpop%dpop%d' % (q, l)
                        setattr(gmodstat, name, featrefrassc[q][l][strgfeat])
            
            # completeness
            for q in gdat.indxrefr:
                if gdat.refr.numbelem[q] == 0:
                    continue
                
                l = gdat.refr.indxpoplfittassc[q]
                
                for nameparaelemfrst in gdat.refr.namepara.elem[q]:
                    
                    if nameparaelemfrst.startswith('etag'):
                        continue
                            
                    if nameparaelemfrst == 'spec' or nameparaelemfrst == 'specplot':
                        continue
                    
                    refrfeatfrst = gdat.refr.dictelem[q][nameparaelemfrst][0, :]
                    binsfeatfrst = getattr(gdat.binspara, nameparaelemfrst)
                        
                    for nameparaelemseco in gdat.refr.namepara.elem[q]:
                        if nameparaelemfrst == nameparaelemseco:
                            continue
                        
                        if nameparaelemseco.startswith('etag'):
                            continue
                        
                        if nameparaelemseco == 'spec' or nameparaelemseco == 'specplot':
                            continue
                        
                        if not checstrgfeat(nameparaelemfrst, nameparaelemseco):
                            continue
                        
                        # temp -- the size of the cmpl np.array should depend on strgmodl
                        cmpltdim = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot)) - 1.
                        
                        if len(indxelemrefrasschits[q][l]) > 0:
                            refrhistfeattdim = getattr(gdat.refr, 'hist%s%spop%d' % (nameparaelemfrst, nameparaelemseco, q))
                            refrfeatseco = gdat.refr.dictelem[q][nameparaelemseco][0, :]
                            binsfeatseco = getattr(gdat.binspara, nameparaelemseco)
                            
                            refrhistfeattdimassc = np.histogram2d(refrfeatfrst[indxelemrefrasschits[q][l]], \
                                                                  refrfeatseco[indxelemrefrasschits[q][l]], bins=(binsfeatfrst, binsfeatseco))[0]
                            indxgood = np.where(refrhistfeattdim != 0.)
                            if indxgood[0].size > 0:
                                cmpltdim[indxgood] = refrhistfeattdimassc[indxgood].astype(float) / refrhistfeattdim[indxgood]
                                if gdat.booldiagmode:
                                    if np.where((cmpltdim[indxgood] > 1.) | (cmpltdim[indxgood] < 0.))[0].size > 0:
                                        raise Exception('')
                        setattr(gmodstat, 'cmpl%s%spop%d' % (nameparaelemfrst, nameparaelemseco, q), cmpltdim)

                    cmplfrst = np.zeros(gdat.numbbinsplot) - 1.
                    if len(indxelemrefrasschits[q][l]) > 0:
                        refrhistfeatfrst = getattr(gdat.refr, 'hist' + nameparaelemfrst + 'pop%d' % q)
                        binsfeatfrst = getattr(gdat.binspara, nameparaelemfrst)
                        refrhistfeatfrstassc = np.histogram(refrfeatfrst[indxelemrefrasschits[q][l]], bins=binsfeatfrst)[0]
                        indxgood = np.where(refrhistfeatfrst != 0.)[0]
                        if indxgood.size > 0:
                            cmplfrst[indxgood] = refrhistfeatfrstassc[indxgood].astype(float) / refrhistfeatfrst[indxgood]
                            if gdat.booldiagmode:
                                if np.where((cmplfrst[indxgood] > 1.) | (cmplfrst[indxgood] < 0.))[0].size > 0:
                                    raise Exception('')
                   
                    setattr(gmodstat, 'cmpl%spop%d' % (nameparaelemfrst, q), cmplfrst)
            
            # false discovery rate
            for l in gmod.indxpopl:
                q = gmod.indxpoplrefrassc[l]
                
                for nameparaelemfrst in gmod.namepara.elem[l]:
                    
                    binsfeatfrst = getattr(gdat.binspara, nameparaelemfrst)
                    for nameparaelemseco in gmod.namepara.elem[l]:
                        
                        if not checstrgfeat(nameparaelemfrst, nameparaelemseco):
                            continue
                        
                        # temp -- the size of the fdis np.array should depend on strgmodl
                        fdistdim = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot))
                        
                        if len(indxelemrefrasschits[q][l]) > 0 and len(gmodstat.dictelem[l][nameparaelemseco]) > 0 and len(gmodstat.dictelem[l][nameparaelemfrst]) > 0: 
                            strgfeattdim = nameparaelemfrst + nameparaelemseco + 'pop%d' % l
                            fitthistfeattdim = getattr(gmodstat, 'hist' + strgfeattdim)
                            binsfeatseco = getattr(gdat.binspara, nameparaelemseco)
                            
                            fitthistfeattdimfals = np.histogram2d(gmodstat.dictelem[l][nameparaelemfrst][indxelemfittasscfals[q][l]], \
                                                  gmodstat.dictelem[l][nameparaelemseco][indxelemfittasscfals[q][l]], bins=(binsfeatfrst, binsfeatseco))[0]
                            indxgood = np.where(fitthistfeattdim != 0.)
                            if indxgood[0].size > 0:
                                fdistdim[indxgood] = fitthistfeattdimfals[indxgood].astype(float) / fitthistfeattdim[indxgood]
                                if gdat.booldiagmode:
                                    if np.where((fdistdim[indxgood] > 1.) | (fdistdim[indxgood] < 0.))[0].size > 0:
                                        raise Exception('')
                        
                        setattr(gmodstat, 'fdis%s%spop%d' % (nameparaelemfrst, nameparaelemseco, l), fdistdim)
                
                    fdisfrst = np.zeros(gdat.numbbinsplot)
                    if len(indxelemrefrasschits[q][l]) > 0 and len(gmodstat.dictelem[l][nameparaelemfrst]) > 0:
                        binsfeatfrst = getattr(gdat.binspara, nameparaelemfrst)
                        fitthistfeatfrstfals = np.histogram(gmodstat.dictelem[l][nameparaelemfrst][indxelemfittasscfals[q][l]], bins=binsfeatfrst)[0]
                        fitthistfeatfrst = getattr(gmodstat, 'hist' + nameparaelemfrst + 'pop%d' % l)
                        indxgood = np.where(fitthistfeatfrst != 0.)[0]
                        if indxgood.size > 0:
                            fdisfrst[indxgood] = fitthistfeatfrstfals[indxgood].astype(float) / fitthistfeatfrst[indxgood]
                            if gdat.booldiagmode:
                                if np.where((fdisfrst[indxgood] > 1.) | (fdisfrst[indxgood] < 0.))[0].size > 0:
                                    raise Exception('')
                    
                    setattr(gmodstat, 'fdis%spop%d' % (nameparaelemfrst, l), fdisfrst)
    
        # temp
        if strgmodl == 'true' and gdat.typeverb > 0:
            for l in gmod.indxpopl:
                for strgfeat in gmod.namepara.genrelem[l]:
                    minm = getattr(gmod.minmpara, strgfeat)
                    maxm = getattr(gmod.maxmpara, strgfeat)
                    if np.where(minm > gmodstat.dictelem[l][strgfeat])[0].size > 0 or np.where(maxm < gmodstat.dictelem[l][strgfeat])[0].size > 0:
                        print('Warning: element parameter outside the plot limits.')
                        print('l')
                        print(l)
                        print('Feature: ')
                        print(strgfeat)
                        print('Plot minmimum')
                        print(minm)
                        print('Plot maxmimum')
                        print(maxm)
                        if strgfeat == gmod.nameparagenrelemampl[l] and strgfeat in gmod.namepara.genrelem[l]:
                            gmod.indxparagenrelemtemp = gmod.namepara.genrelem[l].index(strgfeat)
                            if (gmod.listscalparagenrelem[l][gmod.indxparagenrelemtemp] != 'gaus' and not gmod.listscalparagenrelem[l][gmod.indxparagenrelemtemp].startswith('lnor')):
                                raise Exception('')
    stopchro(gdat, gdatmodi, 'tert')
    
        
def retr_lprielem(gdat, strgmodl, l, g, strgfeat, strgpdfn, paragenrscalfull, dictelem, numbelem):
    
    gmod = getattr(gdat, strgmodl)
    
    if strgpdfn == 'self':
        minmfeat = getattr(gmod.minmpara, strgfeat)
        maxmfeat = getattr(gmod.maxmpara, strgfeat)
        lpri = numbelem[l] * np.log(1. / (maxmfeat - minmfeat))
    if strgpdfn == 'logt':
        lpri = retr_lprilogtdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'gaus':
        lpri = retr_lprigausdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'dexp':
        maxmbgal = getattr(gmod, 'maxmbgal')
        gmod.indxpara.bgaldistscal = getattr(gmod.indxpara, 'bgaldistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_dnp.exp(dictelem[l]['bgal'], maxmbgal, paragenrscalfull[gmod.indxpara.bgaldistscal]))) 
    if strgpdfn == 'expo':
        maxmgang = getattr(gmod, 'maxmgang')
        gang = retr_gang(dictelem[l]['lgal'], dictelem[l]['bgal'])
        gmod.indxpara.gangdistscal = getattr(gmod.indxpara, 'gangdistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_expo(gang, maxmgang, paragenrscalfull[gmod.indxpara.gangdistscal]))) 
        lpri = -numbelem[l] * np.log(2. * pi) 
    if strgpdfn == 'tmpl':
        lpri = np.sum(lpdfspatprioobjt(dictelem[l]['bgal'], dictelem[l]['lgal'], grid=False))
    if strgpdfn == 'powr':
        lpri = retr_lpripowrdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'dpowslopbrek':
        lpri = retr_lpridpowdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'dsrcexpo':
        lpri += -np.sum(np.sqrt((dictelem[l]['lgal'] - lgalsour)**2 + (dictelem[l]['bgal'] - bgalsour)**2) / \
                                                                getattr(gmod, 'dsrcdistsexppop%d' % l))
    if strgpdfn == 'tmpl':
        if strgpdfn.endswith('cons'):
            pdfnspatpriotemp = getattr(gmod, 'pdfnspatpriotemp')
            spatdistcons = paragenrscalfull[getattr(gmod.indxpara, 'spatdistcons')]
            lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
            lpdfspatpriointp = lpdfspatprioobjt(gdat.meanpara.bgalcart, gdat.meanpara.lgalcart)
            lpdfspatpriointp = lpdfspatpriointp.T
            setattr(gmodstat, 'lpdfspatpriointp', lpdfspatpriointp)
            setattr(gmodstat, 'lpdfspatprioobjt', lpdfspatprioobjt)
        else:
            lpdfspatprioobjt = gmod.lpdfspatprioobjt
    
    return lpri


def checstrgfeat(strgfrst, strgseco):

    numbfrst = len(strgfrst)
    numbseco = len(strgseco)
    numb = min(numbfrst, numbseco)
    if strgfrst[:numb] < strgseco[:numb]:
        booltemp = True
    elif strgfrst[:numb] == strgseco[:numb]:
        if numbfrst >= numbseco:
            booltemp = False
        else:
            booltemp = True
    else:
        booltemp = False

    return booltemp


def retr_pathoutprtag(pathpcat, rtag):
    
    pathoutprtag = pathpcat + '/data/outp/' + rtag + '/'
    
    return pathoutprtag


def proc_finl(gdat=None, rtag=None, strgpdfn='post', listnamevarbproc=None, forcplot=False):
    
    gdatmock = None
    
    print('proc_finl()')

    if rtag is None:
        rtag = gdat.rtag
    
    # determine if the final-processing if nominal or tiling
    if isinstance(rtag, list):
        listrtagmodi = rtag
        rtagfinl = tdpy.retr_strgtimestmp() + rtag[0][15:] + 'tile'
        booltile = True
    else:
        listrtagmodi = [rtag]
        rtagfinl = rtag
        booltile = False
    
    # determine of the gdatfinl object is available 
    boolgdatfinl = chec_statfile(pathpcat, rtagfinl, 'gdatfinlpost')
    boolgdatfinlgood = False
    if boolgdatfinl:
        print('Final-processing has been performed previously.')
        pathoutprtag = retr_pathoutprtag(pathpcat, rtagfinl)
        path = pathoutprtag + 'gdatfinl' + strgpdfn
        try:
            gdat = readfile(path) 
            boolgdatfinlgood = True
        except:
            print('gdatfinl object is corrupted.')

    if boolgdatfinl and boolgdatfinlgood:
        # read gdatfinl
        pathoutprtag = retr_pathoutprtag(pathpcat, rtagfinl)
        path = pathoutprtag + 'gdatfinl' + strgpdfn
        gdatfinl = readfile(path) 
        
        if gdatfinl.fitt.numbparaelem > 0:
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.rtagmock is not None:
                        path = gdatfinl.pathoutprtagmock + 'gdatfinlpost'
                        gdatmock = readfile(path)
                    
    else:
        
        if booltile:
            gdatfinltile = tdpy.gdatstrt()
        
        indxrtaggood = []
        liststrgtile = []
        listrtaggood = []
        indxtiletemp = 0
        for n, rtagmodi in enumerate(listrtagmodi):
            
            # read gdatinit
            boolgdatinit = chec_statfile(pathpcat, rtagmodi, 'gdatinit')
            if not boolgdatinit:
                if booltile:
                    print('Initial global object not found. Skipping...')
                    continue
                else:
                    print('Initial global object not found. Quitting...')
                    return
            
            pathoutprtag = retr_pathoutprtag(pathpcat, rtagmodi)
            path = pathoutprtag + 'gdatinit'
            
            gdatinit = readfile(path) 
            if booltile:
                gdatfinltile = gdatinit
                gdatfinl = gdatinit
            else:
                gdatfinl = gdatinit

            pathoutprtagmodi = retr_pathoutprtag(pathpcat, rtagmodi)
            listgdatmodi = []
            for k in gdatinit.indxproc:
                path = pathoutprtagmodi + 'gdatmodi%04d' % k + strgpdfn
                listgdatmodi.append(readfile(path))
            
            # erase
            gdatdictcopy = deepcopy(gdatinit.__dict__)
            for strg, valu in gdatdictcopy.items():
                if strg.startswith('fitt.indxpara.'):
                    delattr(gdatinit, strg)

            if gdatinit.boolmockonly:
                print('Mock only run. Quitting final-processing...')
                return

            # read gdatmodi
            print('rtagmodi')
            print(rtagmodi)
            boolgdatmodi = chec_statfile(pathpcat, rtagmodi, 'gdatmodipost')
            if not boolgdatmodi:
                print('Modified global object not found. Quitting final-processing...')
                return
        
            ## list of other parameters to be flattened
            gdatinit.liststrgvarbarryflat = deepcopy(listgdatmodi[0].liststrgvarbarry)
            # temp
            #for strg in ['memoresi']:
            #    gdatinit.liststrgvarbarryflat.remove(strg)
   
            listparagenrscalfull = np.empty((gdatinit.numbsamptotl, gdatinit.fitt.maxmnumbpara))
            
            if booltile:
                gdatfinltile.pathoutprtag = retr_pathoutprtag(pathpcat, rtagfinl)
                numbsamptotlrsmp = gdatinit.numbsamptotl
                indxsamptotlrsmp = np.random.choice(gdatinit.indxsamptotl, size=gdatinit.numbsamptotl, replace=False)
            
            # aggregate samples from the chains
            if gdatinit.typeverb > 0:
                print('Reading gdatmodi objects from all processes...')
                timeinit = gdatinit.functime()
            
            if gdatinit.typeverb > 0:
                timefinl = gdatinit.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))
            
            if gdatinit.fitt.numbparaelem > 0:
                if len(getattr(listgdatmodi[0], 'list' + strgpdfn + 'gmodstat.indxelemfull')) == 0:
                    print('Found an empty element list. Skipping...')
                    continue
            
            if gdatinit.typeverb > 0:
                print('Accumulating np.arrays...')
                timeinit = gdatinit.functime()
            
            for strgvarb in gdatinit.liststrgvarbarryflat:
                for k in gdatinit.indxproc:
                    if k == 0:
                        shap = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb).shape
                        shap = [shap[0], gdatinit.numbproc] + list(shap[1:])
                        temp = np.zeros(shap) - 1
                    if len(shap) > 2:
                        temp[:, k, :] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)
                    else:
                        temp[:, k] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, temp)
            
            if gdatfinl.typeverb > 0:
                timefinl = gdatfinl.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))
            
            if gdatfinl.typeverb > 0:
                print('Accumulating lists...')
                timeinit = gdatfinl.functime()
            
            # lists of lists collected at each sample
            for strgvarb in listgdatmodi[0].liststrgvarblistsamp:
                listtemp = [[[] for k in gdatfinl.indxproc] for j in gdatfinl.indxsamp]
                for j in gdatfinl.indxsamp:      
                    for k in gdatfinl.indxproc:
                        listtemp[j][k] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)[j]
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listtemp)
            
            if gdatfinl.typeverb > 0:
                timefinl = gdatfinl.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))
            
            if not booltile:
                ## np.maximum likelihood sample 
                gdatfinl.maxmllikproc = np.empty(gdatfinl.numbproc)
                gdatfinl.indxswepmaxmllikproc = np.empty(gdatfinl.numbproc, dtype=int)
                gdatfinl.sampmaxmllikproc = np.empty((gdatfinl.numbproc, gdatfinl.fitt.maxmnumbpara))
                for k in gdatfinl.indxproc:
                    gdatfinl.maxmllikproc[k] = listgdatmodi[k].maxmllikswep
                    gdatfinl.indxswepmaxmllikproc[k] = listgdatmodi[k].indxswepmaxmllik
                    gdatfinl.sampmaxmllikproc[k, :] = listgdatmodi[k].sampmaxmllik
            
                listparagenrscalfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalfull')
                listparagenrunitfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrunitfull')

                # Gelman-Rubin test
                if gdatfinl.numbproc > 1:
                    if gdatfinl.typeverb > 0:
                        print('Computing the Gelman-Rubin TS...')
                        timeinit = gdatfinl.functime()
                    gdatfinl.gmrbparagenrscalbase = np.zeros(gdatfinl.fitt.numbparagenrbase)
                    gdatfinl.gmrbstat = np.zeros((gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt))
                    for k in gdatfinl.fitt.indxparagenrbase:
                        gdatfinl.gmrbparagenrscalbase[k] = tdpy.mcmc.gmrb_test(listparagenrscalfull[:, :, k])
                        if not np.isfinite(gdatfinl.gmrbparagenrscalbase[k]):
                            gdatfinl.gmrbparagenrscalbase[k] = 0.
                    listcntpmodl = getattr(gdatfinl, 'list' + strgpdfn + 'cntpmodl')
                    for i in gdatfinl.indxener:
                        for j in gdatfinl.indxpixl:
                            for m in gdatfinl.indxevtt:
                                gdatfinl.gmrbstat[i, j, m] = tdpy.mcmc.gmrb_test(listcntpmodl[:, :, i, j, m])
                    if gdatfinl.typeverb > 0:
                        timefinl = gdatfinl.functime()
                        print('Done in %.3g seconds.' % (timefinl - timeinit))

                # calculate the autocorrelation of the chains
                if gdatfinl.typeverb > 0:
                    print('Computing the autocorrelation of the chains...')
                    timeinit = gdatfinl.functime()
                gdatfinl.atcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt, int(gdatfinl.numbparagenrfull / 2)))
                gdatfinl.timeatcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt))
                gdatfinl.atcrpara = np.empty((gdatfinl.numbproc, gdatfinl.fitt.maxmnumbpara, int(gdatfinl.numbparagenrfull / 2)))
                gdatfinl.timeatcrpara = np.empty((gdatfinl.numbproc, gdatfinl.fitt.maxmnumbpara))
                for k in gdatfinl.indxproc:
                    gdatfinl.atcrpara[k, :, :], gdatfinl.timeatcrpara[k, :] = tdpy.mcmc.retr_timeatcr(listparagenrscalfull[:, k, :], typeverb=gdatfinl.typeverb)
                    listcntpmodl = getattr(gdatfinl, 'list' + strgpdfn + 'cntpmodl')
                    gdatfinl.atcrcntp[k, :], gdatfinl.timeatcrcntp[k, :] = tdpy.mcmc.retr_timeatcr(listcntpmodl[:, k, :, :, :], typeverb=gdatfinl.typeverb)
                timeatcrcntpmaxm = np.amax(gdatfinl.timeatcrcntp)
                gdatfinl.timeatcrcntpmaxm = np.amax(timeatcrcntpmaxm)
                
                if gdatfinl.typeverb > 0:
                    timefinl = gdatfinl.functime()
                    print('Done in %.3g seconds.' % (timefinl - timeinit))
                
                setattr(gdatfinl, 'list' + strgpdfn + 'sampproc', np.copy(getattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalfull')))

            # flatten the list chains from different walkers
            for strgvarb in listgdatmodi[0].liststrgvarblistsamp:
                listtemp = []
                listinpt = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                for j in gdatfinl.indxsamp:      
                    for k in gdatfinl.indxproc:
                        listtemp.append(listinpt[j][k])
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listtemp)
            
            # flatten the np.array chains from different walkers
            for strgvarb in gdatinit.liststrgvarbarryflat:
                inpt = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                shap = [inpt.shape[0] * inpt.shape[1]] + list(inpt.shape[2:])
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, inpt.reshape(shap))
            listparagenrscalfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalfull')
            listparagenrunitfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrunitfull')
        
            if booltile:

                liststrgtile.append(rtagmodi.split('_')[-2][-4:])
                listrtaggood.append(rtagmodi)
                indxrtaggood.append(n)
                indxtiletemp += 1
                
                if len(liststrgtile) == 1:
                    for strgfeat in gdatfinl.refrgmod.namepara.genrelemtotl:
                        refrfeattile = [[] for q in gdatfinl.indxrefr]
                        setattr(gdatfinl, 'refr' + strgfeat, refrfeattile)
                
                    for strgvarb in gdatfinl.liststrgvarbarrysamp:
                        if not strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                            listvarb = []
                            setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listvarb)
                        else:
                            hist = np.zeros_like(getattr(listgdatmodi[0], 'list' + strgpdfn + strgvarb))
                            setattr(gdatfinl, 'list' + strgpdfn + strgvarb, hist)
                
                    for name, valu in gdatfinl.__dict__.items():
                        if name.startswith('refrhist'):
                            setattr(gdatfinl, name, np.zeros_like(getattr(gdatfinl, name)))
                            
                #for strgfeat in gdatfinl.refrgmod.namepara.genrelemtotl:
                #    refrfeattile = getattr(gdatfinl, 'refr' + strgfeat)
                #    #refrfeat = getattr(gdatfinl, 'refr' + strgfeat)
                #    refrfeat = [[] for q in gdatfinl.indxrefr]
                #    for q in gdatfinl.indxrefr:
                #        if strgfeat in gdatfinl.refrgmod.namepara.genrelem[q]:
                #            refrfeat[q].append(refrfeattile[q])
                
                for strgvarb in gdatfinl.liststrgvarbarrysamp:
                    if strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                        # temp
                        if 'spec' in strgvarb:
                            continue
                        hist = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                        hist += getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                    
                for name, valu in gdatfinl.__dict__.items():
                    if name.startswith('refrhist'):
                        hist = getattr(gdatfinl, name)
                        hist += getattr(gdatfinl, name)

                print('Done with the tile number %d, run number %d...' % (indxtiletemp, n))
        
        if booltile:
            gdatfinl.pathplotrtag = gdatfinl.pathimag + rtagfinl + '/'
            make_fold(gdatfinl)
            indxrtaggood = np.array(indxrtaggood).astype(int)
            numbrtaggood = indxrtaggood.size
            numbtile = numbrtaggood
            print('Found %d tiles with run tags:' % numbrtaggood)
            for indxrtaggoodtemp in indxrtaggood:
                print(rtag[indxrtaggoodtemp])

            # np.concatenate reference elements from different tiles
            #for strgfeat in gdatfinl.refrgmod.namepara.genrelemtotl:
            #    refrfeat = getattr(gdatfinl, 'refr' + strgfeat, refrfeat)
            #    for q in gdatfinl.indxrefr:
            #        if strgfeat in gdatfinl.refrgmod.namepara.genrelem[q]:
            #            refrfeat[q] = np.concatenate(refrfeat[q], axis=1)
            
            for strgvarb in gdatfinl.liststrgvarbarrysamp:
                if not strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                    listvarb = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                    if 'assc' in strgvarb:
                        numbrefrelemtotl = 0
                        for k, varbrsmp in enumerate(listvarb):
                            numbrefrelemtotl += varbrsmp.shape[1]
                        shap = [gdatfinl.numbsamptotl, numbrefrelemtotl]
                        listvarbtemp = np.empty(shap)
                        cntr = 0
                        for k, varb in enumerate(listvarb):
                            listvarbtemp[:, cntr:cntr+varb.shape[1]] = varb
                            cntr += varb.shape[1]
                    else:
                        shap = [gdatfinl.numbsamptotl * numbtile] + list(listvarb[0].shape[1:])
                        listvarbtemp = np.empty(shap)
                        for k, varb in enumerate(listvarb):
                            listvarbtemp[k*gdatfinl.numbsamptotl:(k+1)*gdatfinl.numbsamptotl, ...] = varb
                    setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listvarbtemp)
        else:
            # np.maximum likelihood sample
            if gdatfinl.fitt.numbparaelem > 0:
                listindxelemfull = getattr(gdatfinl, 'list' + strgpdfn + 'indxelemfull')
            listllik = getattr(gdatfinl, 'list' + strgpdfn + 'llik')
            listlliktotl = getattr(gdatfinl, 'list' + strgpdfn + 'lliktotl')
            indxsamptotlmlik = np.argmax(np.sum(np.sum(np.sum(listllik, 3), 2), 1))
            
            # copy the np.maximum likelihood sample
            for strgvarb in listgdatmodi[0].liststrgvarbarrysamp:
                setattr(gdatfinl, 'mlik' + strgvarb, getattr(gdatfinl, 'list' + strgpdfn + strgvarb)[indxsamptotlmlik, ...])
            for strgvarb in listgdatmodi[0].liststrgvarblistsamp:
                setattr(gdatfinl, 'mlik' + strgvarb, getattr(gdatfinl, 'list' + strgpdfn + strgvarb)[indxsamptotlmlik])

            # temp -- dont gdatfinl.listllik and gdatfinl.listparagenrscalfull have the same dimensions?
            gdatfinl.mlikparagenrscalfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalfull')[indxsamptotlmlik, :]
            gdatfinl.mlikparagenrscalfull = getattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalfull')[indxsamptotlmlik, :]
            #if gdatfinl.fitt.numbparaelem > 0:
            #    gdatfinl.mlikindxelemfull = listindxelemfull[indxsamptotlmlik]
            gdatfinl.mlikparagenrscalbase = gdatfinl.mlikparagenrscalfull[gdatfinl.fitt.indxparagenrbase]
            for k, gmod.nameparagenrbase in enumerate(gdatfinl.fitt.nameparagenrbase):
                setattr(gdatfinl, 'mlik' + gmod.nameparagenrbase, gdatfinl.mlikparagenrscalbase[k])

            # add execution times to the chain output
            gdatfinl.timereal = np.zeros(gdatfinl.numbproc)
            gdatfinl.timeproc = np.zeros(gdatfinl.numbproc)
            for k in gdatfinl.indxproc:
                gdatfinl.timereal[k] = listgdatmodi[k].timereal
                gdatfinl.timeproc[k] = listgdatmodi[k].timeproc
        
            # find the np.maximum likelihood and posterior over the chains
            gdatfinl.indxprocmaxmllik = np.argmax(gdatfinl.maxmllikproc)
            #gdatfinl.maxmlliktotl = gdatfinl.maxmllikproc[gdatfinl.indxprocmaxmllik]
            gdatfinl.indxswepmaxmllik = gdatfinl.indxprocmaxmllik * gdatfinl.numbparagenrfull + gdatfinl.indxswepmaxmllikproc[gdatfinl.indxprocmaxmllik]
            gdatfinl.sampmaxmllik = gdatfinl.sampmaxmllikproc[gdatfinl.indxprocmaxmllik, :]
                
            if strgpdfn == 'post':
                levipost = retr_levipost(listlliktotl)
                setattr(gdatfinl, strgpdfn + 'levipost', levipost)
            
            if strgpdfn == 'prio':
                leviprio = np.log(np.mean(np.exp(listlliktotl)))
                setattr(gdatfinl, strgpdfn + 'leviprio', leviprio)
            
        # parse the sample vector
        listparagenrscalbase = listparagenrscalfull[:, gdatfinl.fitt.indxparagenrbase]
        for k, gmod.nameparagenrbase in enumerate(gdatfinl.fitt.nameparagenrbase):
            setattr(gdatfinl, 'list' + strgpdfn + gmod.nameparagenrbase, listparagenrscalbase[:, k])
        setattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalbase', listparagenrscalbase)

        if strgpdfn == 'post' and gdatfinl.checprio:
            pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
            path = pathoutprtag + 'gdatfinlprio'
            try:
                gdatprio = readfile(path)
            except:
                proc_finl(gdat=gdatfinl, strgpdfn='prio', listnamevarbproc=listnamevarbproc, forcplot=forcplot)
        else:
            gdatprio = None
        
        # post process samples
        ## bin element parameters
        if gdatfinl.typeverb > 0:
            print('Binning the probabilistic catalog spatially...')
            timeinit = gdatfinl.functime()
        
        if not booltile:
            if gdatfinl.fitt.numbparaelem > 0:
                if gdatfinl.boolbinsspat:
                    histlgalbgalelemstkd = [[] for l in gdatfinl.fittindxpopl]
                
                    listlgal = getattr(gdatfinl, 'list' + strgpdfn + 'lgal')
                    listbgal = getattr(gdatfinl, 'list' + strgpdfn + 'bgal')
                    for l in gdatfinl.fittindxpopl:
                        if gdatfinl.fitttypeelem[l] != 'lghtline':
                            histlgalbgalelemstkd[l] = np.zeros((gdatfinl.numbbgalpntsprob, gdatfinl.numblgalpntsprob, gdatfinl.numbbinsplot, numb))
                            temparry = np.concatenate([listlgal[n][l] for n in gdatfinl.indxsamptotl])
                            temp = np.empty((len(temparry), 3))
                            temp[:, 0] = temparry
                            temp[:, 1] = np.concatenate([listbgal[n][l] for n in gdatfinl.indxsamptotl])
                            temp[:, 2] = np.concatenate([getattr(gdatfinl, 'list' + strgpdfn + strgfeat)[n][l] for n in gdatfinl.indxsamptotl])
                            bins = getattr(gdatfinl, 'bins' + strgfeat)
                            histlgalbgalelemstkd[l][:, :, :, k] = np.histogramdd(temp, \
                                                                    bins=(gdatfinl.binslgalpntsprob, gdatfinl.binsbgalpntsprob, bins))[0]
                    setattr(gdatfinl, strgpdfn + 'histlgalbgalelemstkd', histlgalbgalelemstkd)

            if gdatfinl.typeverb > 0:
                timefinl = gdatfinl.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))

            ## construct a condensed catalog of elements
            if gdatfinl.boolcondcatl and gdatfinl.fitt.numbparaelem > 0:
                
                if gdatfinl.typeverb > 0:
                    print('Constructing a condensed catalog...')
                    timeinit = gdatfinl.functime()
                
                retr_condcatl(gdatfinl)
            
                if gdatfinl.typeverb > 0:
                    timefinl = gdatfinl.functime()
                    print('Done in %.3g seconds.' % (timefinl - timeinit))

            # construct lists of samples for each proposal type
            listindxproptype = getattr(gdatfinl, 'list' + strgpdfn + 'indxproptype')
            listboolpropaccp = getattr(gdatfinl, 'list' + strgpdfn + 'boolpropaccp')
            listboolpropfilt = getattr(gdatfinl, 'list' + strgpdfn + 'boolpropfilt')
            listindxsamptotlproptotl = []
            listindxsamptotlpropfilt = []
            listindxsamptotlpropaccp = []
            listindxsamptotlpropreje = []
            for n in gdatfinl.indxproptype:
                indxsampproptype = np.where(listindxproptype == gdatfinl.indxproptype[n])[0]
                listindxsamptotlproptotl.append(indxsampproptype)
                listindxsamptotlpropaccp.append(np.intersect1d(indxsampproptype, np.where(listboolpropaccp)[0]))
                listindxsamptotlpropfilt.append(np.intersect1d(indxsampproptype, np.where(listboolpropfilt)[0]))
                listindxsamptotlpropreje.append(np.intersect1d(indxsampproptype, np.where(np.logical_not(listboolpropaccp))[0]))
                if listindxsamptotlproptotl[n].size == 0:
                    accp = 0.
                else:
                    accp = float(listindxsamptotlpropaccp[n].size) / listindxsamptotlproptotl[n].size
                setattr(gdatfinl, 'accp' + gdatfinl.nameproptype[n], accp)

            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlproptotl', listindxsamptotlproptotl)
            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlpropaccp', listindxsamptotlpropaccp)
            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlpropreje', listindxsamptotlpropreje)
       
        if gdatfinl.fitt.numbparaelem > 0 and strgpdfn == 'post':
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.rtagmock is not None:
                        path = gdatfinl.pathoutprtagmock + 'gdatfinlpost'
                        gdatmock = readfile(path)
                    
        # posterior corrections
        if gdatfinl.fitt.numbparaelem > 0 and strgpdfn == 'post':

            ## perform corrections
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:

                    for gmod.namepara.genrelemvarbhist in gdatfinl.liststrgvarbhist:
                        strgvarb = gmod.namepara.genrelemvarbhist[0]

                        if gmod.namepara.genrelemvarbhist[1].startswith('aerr') or len(gmod.namepara.genrelemvarbhist[2]) > 0 and gmod.namepara.genrelemvarbhist[2].startswith('aerr'):
                            continue
                        if gmod.namepara.genrelemvarbhist[1] == 'spec' or gmod.namepara.genrelemvarbhist[1] == 'deflprof' or gmod.namepara.genrelemvarbhist[1] == 'specplot':
                            continue
                        if len(gmod.namepara.genrelemvarbhist[2]) > 0 and (gmod.namepara.genrelemvarbhist[2] == 'spec' or \
                                    gmod.namepara.genrelemvarbhist[2] == 'deflprof' or gmod.namepara.genrelemvarbhist[2] == 'specplot'):
                            continue
                        
                        ## internal correction
                        listhist = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                        
                        for qq in gdatmock.indxrefr:
                            l = int(gmod.namepara.genrelemvarbhist[3][qq].split('pop')[1][0])
                            qq = int(gmod.namepara.genrelemvarbhist[3][qq].split('pop')[2][0])
                            if gmod.namepara.genrelemvarbhist[1][-4:] in gdatfinl.listnamerefr and \
                                    (len(gmod.namepara.genrelemvarbhist[2]) == 0 or gmod.namepara.genrelemvarbhist[2][-4:] in gdatfinl.listnamerefr):
                                listhistincr = listhist
                            else:
                                if gmod.namepara.genrelemvarbhist[1][-4:] in gdatfinl.listnamerefr and len(gmod.namepara.genrelemvarbhist[2]) > 0:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostcmpl' + gmod.namepara.genrelemvarbhist[2] + 'pop%dpop%d' % (l, qq))], 2)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostfdis' + gmod.namepara.genrelemvarbhist[2] + 'pop%dpop%d' % (qq, l))], 2)
                                elif len(gmod.namepara.genrelemvarbhist[2][:-4]) > 0 and gmod.namepara.genrelemvarbhist[2][-4:] in gdatfinl.listnamerefr:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostcmpl' + gmod.namepara.genrelemvarbhist[1] + 'pop%dpop%d' % (l, qq))], 1)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostfdis' + gmod.namepara.genrelemvarbhist[1] + 'pop%dpop%d' % (qq, l))], 1)
                                else:
                                    listcmpltrue = getattr(gdatmock, 'listpostcmpl' + gmod.namepara.genrelemvarbhist[3][qq])
                                    listfdistrue = getattr(gdatmock, 'listpostfdis' + gmod.namepara.genrelemvarbhist[3][qq])
                                if len(gmod.namepara.genrelemvarbhist[2]) == 0:
                                    listcmplboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    listfdisboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    listhistboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    for k in gdatfinl.indxbinsplot:
                                        listcmplboot[:, k] = np.random.choice(listcmpltrue[:, k], size=gdatfinl.numbsampboot)
                                        listfdisboot[:, k] = np.random.choice(listfdistrue[:, k], size=gdatfinl.numbsampboot)
                                        listhistboot[:, k] = np.random.choice(listhist[:, k], size=gdatfinl.numbsampboot)
                                else:
                                    listcmplboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    listfdisboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    listhistboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    for a in gdatfinl.indxbinsplot:
                                        for b in gdatfinl.indxbinsplot:
                                            listcmplboot[:, a, b] = np.random.choice(listcmpltrue[:, a, b], size=gdatfinl.numbsampboot)
                                            listfdisboot[:, a, b] = np.random.choice(listfdistrue[:, a, b], size=gdatfinl.numbsampboot)
                                            listhistboot[:, a, b] = np.random.choice(listhist[:, a, b], size=gdatfinl.numbsampboot)
                                indxbadd = np.where(listcmplboot == -1)
                                indxbaddzero = np.where(listcmplboot == 0.)
                                listhistincr = listhistboot / listcmplboot * (1. - listfdisboot)
                                listhistincr[indxbadd] = -1.5
                                listhistincr[indxbaddzero] = 1.5
                            
                            listgdatmodi[0].liststrgchan += ['incr' + gmod.namepara.genrelemvarbhist[4][qq]]
                            setattr(gdatfinl, 'listpostincr' + gmod.namepara.genrelemvarbhist[4][qq], listhistincr)
                        
                            ## external correction
                            for q in gdatfinl.indxrefr:
                                nametemp = gmod.namepara.genrelemvarbhist[1] 
                                if len(gmod.namepara.genrelemvarbhist[2]) > 0:
                                    nametemp += gmod.namepara.genrelemvarbhist[2]
                                nametemp += 'pop%dpop%dpop%d' % (q, qq, l)
                                crexhist = getattr(gdatfinl, 'crex' + nametemp)
                                if crexhist is not None:
                                    
                                    listhistexcr = listhistincr * crexhist 
                                    
                                    if crexhist.ndim == 1 and listhistincr.ndim == 3:
                                        raise Exception('')
                                    
                                    listgdatmodi[0].liststrgchan += ['excr' + nametemp]
                                    setattr(gdatfinl, 'listpostexcr' + nametemp, listhistexcr)
                            
        # compute credible intervals
        if gdatfinl.typeverb > 0:
            print('Computing credible intervals...')
            timeinit = gdatfinl.functime()
       
        for strgchan in listgdatmodi[0].liststrgchan:
            
            if booltile:
                if strgchan in gdatfinl.liststrgvarbarryswep or strgchan in listgdatmodi[0].liststrgvarblistsamp:
                    continue
                if not (strgchan.startswith('hist') or strgchan.startswith('incr') or strgchan.startswith('excr')):
                    continue

            if gdatfinl.fitt.numbparaelem > 0 and strgchan in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                if 'spec' in strgchan:
                    continue
            if strgchan == 'spec':
                continue

            listtemp = getattr(gdatfinl, 'list' + strgpdfn + strgchan)
            
            if isinstance(listtemp, list):
            
                if booltile:
                    continue

                # ensure that transdimensional lists are not included
                # temp
                if strgchan in gdatfinl.fitt.namepara.genrelemtotl or strgchan == 'indxelemfull':
                    continue

                pctltemp = []
                pmeatemp = []
                meditemp = []
                errrtemp = []
                stdvtemp = []
                numb = len(listtemp[0])
                
                for k in range(numb):
                    if isinstance(listtemp[0][k], list):
                        continue
                    shap = [gdatfinl.numbsamptotl] + list(listtemp[0][k].shape)
                    temp = np.zeros(shap)
                    for n in gdatfinl.indxsamptotl:
                        temp[n, ...] = listtemp[n][k]
                    
                    pctltempsing = tdpy.retr_pctlvarb(temp)
                    pmeatempsing = np.mean(temp, axis=0)
                    meditempsing = pctltempsing[0, ...]
                    errrtempsing = tdpy.retr_errrvarb(pctltempsing)
                    stdvtempsing = np.std(temp)
                    
                    pctltemp.append(pctltempsing)
                    pmeatemp.append(pmeatempsing)
                    meditemp.append(meditempsing)
                    errrtemp.append(errrtempsing)
                    stdvtemp.append(stdvtempsing)
            else:
                # this is needed for finding posterior moments of features of associated reference elements
                if 'asscref' in strgchan:
                    if listtemp.ndim != 2:
                        raise Exception('')
                    pmeatemp = np.zeros(listtemp.shape[1])
                    pctltemp = np.zeros([3] + [listtemp.shape[1]])
                    # temp -- this only works for 2D listtemp
                    for k in range(listtemp.shape[1]):
                        indxassc = np.where(np.isfinite(listtemp[:, k]))[0]
                        if indxassc.size > 0:
                            pctltemp[:, k] = tdpy.retr_pctlvarb(listtemp[indxassc, k])
                            pmeatemp[k] = np.mean(listtemp[indxassc, k])
                else:
                    pctltemp = tdpy.retr_pctlvarb(listtemp)
                    pmeatemp = np.mean(listtemp, axis=0)
                
                errrtemp = tdpy.retr_errrvarb(pctltemp)
                stdvtemp = np.std(pctltemp, axis=0)
                meditemp = pctltemp[0, ...]
                
                if strgchan in gdatfinl.listnamevarbcpct:
                    cpcttemp = np.empty([gdatfinl.numbsampcpct] + [3] + list(listtemp.shape[1:]))
                    for n in gdatfinl.indxsampcpct:
                        cpcttemp[n, ...] = tdpy.retr_pctlvarb(listtemp[:n+1, ...])
            
            setattr(gdatfinl, 'pctl' + strgpdfn + strgchan, pctltemp)
            setattr(gdatfinl, 'medi' + strgpdfn + strgchan, meditemp)
            setattr(gdatfinl, 'pmea' + strgpdfn + strgchan, pmeatemp)
            setattr(gdatfinl, 'errr' + strgpdfn + strgchan, errrtemp)
            setattr(gdatfinl, 'stdv' + strgpdfn + strgchan, stdvtemp)
            if strgchan in gdatfinl.listnamevarbcpct:
                setattr(gdatfinl, 'cpct' + strgpdfn + strgchan, cpcttemp)
        
        if not booltile:
            pmealliktotl = getattr(gdatfinl, 'pmea' + strgpdfn + 'lliktotl')
            stdvlliktotl = getattr(gdatfinl, 'stdv' + strgpdfn + 'lliktotl')
            minmlliktotl = np.amin(listlliktotl)
            maxmlliktotl = np.amax(listlliktotl)
            skewlliktotl = np.mean(((listlliktotl - pmealliktotl) / stdvlliktotl)**3)
            kurtlliktotl = np.mean(((listlliktotl - pmealliktotl) / stdvlliktotl)**4)
            setattr(gdatfinl, 'minm' + strgpdfn + 'lliktotl', minmlliktotl)
            setattr(gdatfinl, 'maxm' + strgpdfn + 'lliktotl', maxmlliktotl)
            setattr(gdatfinl, 'skew' + strgpdfn + 'lliktotl', skewlliktotl)
            setattr(gdatfinl, 'kurt' + strgpdfn + 'lliktotl', kurtlliktotl)

            if strgpdfn == 'post':
                infopost = retr_infofromlevi(pmealliktotl, levipost)
                setattr(gdatfinl, strgpdfn + 'infopost', infopost)
            if strgpdfn == 'post' and gdatfinl.checprio:
                leviprio = getattr(gdatprio, 'prioleviprio')
                infoprio = retr_infofromlevi(pmealliktotl, leviprio)
                setattr(gdatfinl, strgpdfn + 'infoprio', infoprio)
            
            bcom = maxmlliktotl - pmealliktotl
            setattr(gdatfinl, strgpdfn + 'bcom', bcom)
        
        listnametemp = ['lliktotl']
        if gmod.numbparaelem > 0:
            listnametemp += ['lpripena']

        for namevarbscal in listnametemp:
            listtemp = getattr(gdatfinl, 'list' + strgpdfn + namevarbscal)
            minm = np.amin(listtemp)
            maxm = np.amax(listtemp)
            setattr(gdatfinl, 'minm' + namevarbscal, minm)
            setattr(gdatfinl, 'maxm' + namevarbscal, maxm)
            setattr(gdatfinl, 'scal' + namevarbscal, 'self')
            retr_axis(gdat, namevarbscal)
        
        if gdatfinl.checprio:
            for strgvarb in gdatfinl.listnamevarbscal:
                setp_pdfnvarb(gdatfinl, strgpdfn, strgvarb, strgvarb)
            for l0 in gdatfinl.fittindxpopl:
                for strgfeatfrst in gdatfinl.fitt.namepara.genrelem[l0]:
                    if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                        continue
                    setp_pdfnvarb(gdatfinl, strgpdfn, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l0)
                    for strgfeatseco in gdatfinl.fitt.namepara.genrelem[l0]:
                        if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                            continue
                        
                        if not checstrgfeat(strgfeatfrst, strgfeatseco):
                            continue
                                    
                        setp_pdfnvarb(gdatfinl, strgpdfn, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l0, nameseco=strgfeatseco)

            # calculate information gain
            if strgpdfn == 'post':
                for namevarbscal in gdatfinl.listnamevarbscal:
                    setp_info(gdatfinl, gdatprio, namevarbscal, namevarbscal)
                for l0 in gdatfinl.fittindxpopl:
                    for strgfeatfrst in gdatfinl.fitt.namepara.genrelem[l0]:
                        if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                            continue
                        setp_info(gdatfinl, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l0)
                        for strgfeatseco in gdatfinl.fitt.namepara.genrelem[l0]:
                            if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                                continue
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                                    
                            setp_info(gdatfinl, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l0, nameseco=strgfeatseco)

        if gdatfinl.typeverb > 0:
            timefinl = gdatfinl.functime()
            print('Done in %.3g seconds.' % (timefinl - timeinit))
        
        # flatten the np.arrays which have been collected at each sweep
        #setattr(gdat, 'list' + strgpdfn + strgpdfntemp + 'flat', getattr(gdat, 'list' + strgpdfn + strgpdfntemp + 'totl').flatten())
        if not booltile:
            # memory usage
            listmemoresi = getattr(gdatfinl, 'list' + strgpdfn + 'memoresi')
            gdatfinl.meanmemoresi = np.mean(listmemoresi, 1)
            gdatfinl.derimemoresi = (gdatfinl.meanmemoresi[-1] - gdatfinl.meanmemoresi[0]) / gdatfinl.numbswep

            gdatfinl.timerealtotl = time.time() - gdatfinl.timerealtotl
            gdatfinl.timeproctotl = time.clock() - gdatfinl.timeproctotl
            gdatfinl.timeproctotlswep = gdatfinl.timeproctotl / gdatfinl.numbswep
            
            if gdatfinl.timeatcrcntpmaxm == 0.:
                gdatfinl.timeprocnorm = 0.
            else:
                gdatfinl.timeprocnorm = gdatfinl.timeproctotlswep / gdatfinl.timeatcrcntpmaxm
   
        # write the final gdat object
        path = gdatfinl.pathoutprtag + 'gdatfinl' + strgpdfn

        if gdatfinl.typeverb > 0:
            print('Writing gdatfinl to %s...' % path)
        writfile(gdatfinl, path) 
       
        filestat = open(gdatfinl.pathoutprtag + 'stat.txt', 'a')
        filestat.write('gdatfinl%s written.\n' % strgpdfn)
        filestat.close()
   
        if not booltile:
            if gdatfinl.typeverb > 0:
                for k in gdatfinl.indxproc:
                    print('Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdatfinl.timereal[k], gdatfinl.timeproc[k]))
                print('Parent process has run in %d real seconds, %d CPU seconds.' % (gdatfinl.timerealtotl, gdatfinl.timeproctotl))
    
    print('HACKING!!')
    gdatfinl.strgpdfn = 'post'

    print('Checking whether post-processing plots already exist.')
    booltemp = chec_statfile(pathpcat, rtagfinl, 'plotfinl')
    if booltemp:
        print('Final plots already exist. Skipping...')
    else:
        if strgpdfn == 'post' and gdatfinl.checprio:
            path = pathoutprtag + 'gdatfinlprio'
            gdatprio = readfile(path)
        else:
            gdatprio = None
        
        if gdatfinl.makeplot and getattr(gdatfinl, 'makeplotfinl' + strgpdfn) or forcplot:
            plot_finl(gdatfinl, gdatprio=gdatprio, strgpdfn=strgpdfn, gdatmock=gdatmock, booltile=booltile)
            filestat = open(gdatfinl.pathoutprtag + 'stat.txt', 'a')
            filestat.write('plotfinl%s written.\n' % strgpdfn)
            filestat.close()


def retr_listgdat(listrtag, typegdat='finlpost'):
   
    listgdat = []
    for rtag in listrtag:
        pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
        path = pathoutprtag + 'gdat%s' % typegdat
        listgdat.append(readfile(path))

    return listgdat


def make_fold(gdat):

    for strgpdfn in gdat.liststrgpdfn:
        setattr(gdat, 'path' + strgpdfn, gdat.pathplotrtag + strgpdfn + '/') 
        path = getattr(gdat, 'path' + strgpdfn)

        for nameseco in ['finl', 'fram', 'anim', 'opti']:
            setattr(gdat, 'path' + strgpdfn + nameseco, path + nameseco + '/')
        
        for nameseco in ['diag', 'lpac', 'varbscal', 'cond', 'varbscalproc']:
            setattr(gdat, 'path' + strgpdfn + 'finl' + nameseco, path + 'finl/' + nameseco + '/')
        
        for n in gdat.indxproptype:
            setattr(gdat, 'path' + strgpdfn + 'finl' + gdat.nameproptype[n], path + 'finl/lpac/' + gdat.nameproptype[n] + '/')

        for namethrd in ['hist', 'trac', 'join', 'cova']:
            setattr(gdat, 'path' + strgpdfn + 'finlvarbscal' + namethrd, path + 'finl/varbscal/' + namethrd + '/')
            
        for strgphas in gdat.liststrgphas + ['init']:
            liststrgfold = getattr(gdat, 'liststrgfold' + strgphas)
            for nameseco in liststrgfold:
                if strgphas == 'init':
                    if nameseco == 'assc' or nameseco.startswith('cmpl') or nameseco.startswith('fdis'):
                        continue
                    setattr(gdat, 'path' + strgphas + nameseco[:-1], gdat.pathplotrtag + 'init/' + nameseco)
                else:
                    setattr(gdat, 'path' + strgpdfn + strgphas + nameseco[:-1], path + strgphas + '/' + nameseco)
    gdat.pathinfo = gdat.pathplotrtag + 'info/'
    
    ## make the directories 
    for attr, valu in gdat.__dict__.items():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)


def make_cmapdivg(strgcolrloww, strgcolrhigh):
    
    funccolr = mpl.colors.ColorConverter().to_rgb
    
    colrloww = funccolr(strgcolrloww)
    colrhigh = funccolr(strgcolrhigh)
    
    cmap = make_cmap([colrloww, funccolr('white'), 0.5, funccolr('white'), colrhigh])

    return cmap


def make_cmap(seq):
    
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    
    return mpl.colors.LinearSegmentedColormap('CustomMap', cdict)


def setp_pdfnvarb(gdat, strgpdfn, name, namefull, nameseco=None):
    
    if listvarb.ndim == 1:
        shaptemp = [gdat.numbbinspdfn, 1]
    else:
        shaptemp = [gdat.numbbinspdfn] + list(listvarb.shape[1:])
    pdfn = np.empty(shaptemp)
    if listvarb.ndim == 1:
        binsvarb = getattr(gdat.binspara, name)
        deltvarb = getattr(gdat, 'delt' + name)
        pdfn[:, 0] = np.histogram(listvarb, bins=binsvarb)[0].astype(float)
        pdfn[:, 0] /= np.sum(pdfn[:, 0])
        pdfn[:, 0] /= deltvarb
    else:
        binsvarb = np.linspace(0, gmod.maxmpara.numbelemtotl, 51)
        
    if listvarb.ndim == 2:
        for k in range(listvarb.shape[1]):
            pdfn[:, k] = np.histogram(listvarb[:, k], bins=binsvarb)[0].astype(float)
            pdfn[:, k] /= np.sum(pdfn[:, k])
        pdfn *= 50.
    if listvarb.ndim == 3:
        for k in range(listvarb.shape[1]):
            for m in range(listvarb.shape[2]):
                pdfn[:, k, m] = np.histogram(listvarb[:, k, m], bins=binsvarb)[0].astype(float)
                pdfn[:, k, m] /= np.sum(pdfn[:, k, m])
        pdfn *= 2500.
    pdfn[np.where(pdfn < 1e-50)[0]] = 1e-50
    
    setattr(gdat, 'pdfn' + strgpdfn + namefull, pdfn)


def setp_info(gdat, gdatprio, name, namefull, nameseco=None, namesecofull=None):
    
    listpost = getattr(gdat, 'listpost' + namefull)
    listprio = getattr(gdatprio, 'listprio' + namefull)
    pdfnpost = getattr(gdat, 'pdfnpost' + namefull)
    pdfnprio = getattr(gdatprio, 'pdfnprio' + namefull)
    if listpost.ndim == 3:
        infodens = np.empty((gdat.numbbinspdfn, listpost.shape[1], listpost.shape[2]))
        info = np.empty((listpost.shape[1], listpost.shape[2]))
        pvks = np.empty((listpost.shape[1], listpost.shape[2]))
    else:
        if listpost.ndim == 1:
            numbtemp = 1
        else:
            numbtemp = listpost.shape[1]
        infodens = np.empty((gdat.numbbinspdfn, numbtemp))
        info = np.empty(numbtemp)
        pvks = np.empty(numbtemp)
    if listpost.ndim == 1:
        listpost = listpost[:, None]
        listprio = listprio[:, None]
        deltvarb = getattr(gdat, 'delt' + name)
    else:
        if listpost.ndim == 2:
            deltvarb = 1. / 50
        else:
            deltvarb = 1. / 50**list2
    
    if listpost.ndim == 1 or listpost.ndim == 2:
        for k in range(listpost.shape[1]):
            infodens[:, k] = retr_infodens(pdfnpost[:, k], pdfnprio[:, k])
            info[k] = np.sum(infodens[:, k] * deltvarb)
            temp, pvks[k] = sp.stats.ks_2samp(listpost[:, k], listprio[:, k])
    if listpost.ndim == 3:
        for k in range(listpost.shape[1]):
            for m in range(listpost.shape[2]):
                infodens[:, k, m] = retr_infodens(pdfnpost[:, k, m], pdfnprio[:, k, m])
                info[k, m] = np.sum(infodens[:, k, m] * deltvarb)
                temp, pvks[k, m] = sp.stats.ks_2samp(listpost[:, k, m], listprio[:, k, m])
    
    setattr(gdat, 'pvks' + namefull, pvks)
    setattr(gdat, 'infodens' + namefull, infodens)
    setattr(gdat, 'info' + namefull, info)


# check the state file
def chec_statfile(pathpcat, rtag, strggdat, typeverb=1):
    
    print('Checking the state file %s for %s...' % (strggdat, rtag))
    
    pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
    
    # check the status file
    if not os.path.isfile(pathoutprtag + 'stat.txt'):
        if typeverb > 0:
            print('pathoutprtag')
            print(pathoutprtag)
            print('stat.txt not found.')
        return False

    # check the global object
    filestat = open(pathoutprtag + 'stat.txt', 'r')
    booltemp = False
    linesrch = strggdat + ' written.\n'
    for line in filestat:
        if line == linesrch:
            booltemp = True

    filestat.close()
    if not booltemp:
        if typeverb > 0:
            print('bad %s status.' % (strggdat))
        return False
    else:
        return True


def retr_los3(dlos, lgal, bgal):

    dglc = np.sqrt(8.5e3**2 + dlos**2 - 2. * dlos * 8.5e3 * np.cos(bgal) * np.cos(lgal))
    thet = np.arccos(np.sin(bgal) * dlos / dglc)
    phii = np.arcsin(np.sqrt(np.cos(bgal)**2 * dlos**2 + 8.5e3**2 - 2 * dlos * np.cos(bgal) * 8.5e3) / dglc)
    
    return dglc, thet, phii


def retr_glc3(dglc, thet, phii):

    xpos = dglc * np.sin(thet) * np.cos(phii)
    ypos = dglc * np.sin(thet) * np.sin(phii)
    zpos = dglc * np.cos(thet)
    dlos = np.sqrt(zpos**2 + xpos**2 + (8.5e3 - ypos)**2)
    lgal = np.arctan2(8.5e3 - ypos, xpos) - np.pi / 2
    bgal = np.arcsin(zpos / dlos)
   
    return dlos, lgal, bgal


def retr_lumipuls(geff, magf, per0):

    # temp -- this is bolometric luminosity np.whereas dictelem[l]['flux'] is differential!
    lumi = 9.6e33 * (geff / 0.2) * (magf / 10**8.5)**2 * (3e-3 / per0)*4

    return lumi


def retr_lumi(gdat, flux, dlos, reds=None):

    lumi = flux * 4. * np.pi * dlos**2 * gdat.prsccmtr**2 / gdat.ergsgevv
    
    # temp
    # redshift correction
    if reds is not None:
        lumi *= (1. + reds)**2

    return lumi


def retr_flux(gdat, lumi, dlos, reds=None):

    flux = lumi / 4. / np.pi / dlos**2 / gdat.prsccmtr**2 * gdat.ergsgevv
    
    # temp
    # redshift correction
    if reds is not None:
        pass

    return flux


def retr_per1(per0, magf):

    per1 = 3.3e-20 * (magf / 10**8.5)**2 * (3e-3 / per0)

    return per1


def retr_dlosgalx(lgal, bgal, dglc):

    # temp -- this is obviously wrong
    dlos = 8.5e3 - dglc

    return dlos


def retr_arryfromlist(listtemp):
    
    shap = [len(listtemp)] + list(listtemp[0].shape)
    arry = np.empty(shap)
    for k in range(len(listtemp)):
        arry[k, ...] = listtemp[k]
    
    return arry


def proc_cntpdata(gdat):

    # exclude voxels with vanishing exposure
    ## data counts
    if gdat.typedata == 'inpt':
        gdat.cntpdata = retr_cntp(gdat, gdat.sbrtdata)
    
    # data variance
    gdat.varidata = np.maximum(gdat.cntpdata, 1.)

    # correct the likelihoods for the constant data dependent factorial
    gdat.llikoffs = -sp.special.gammaln(gdat.cntpdata + 1)

    ## spatial average
    gdat.sbrtdatamean, gdat.sbrtdatastdv = retr_spatmean(gdat, gdat.cntpdata, boolcntp=True)
    
    # data count limits
    minmcntpdata = np.amin(gdat.cntpdata)
    maxmcntpdata = np.amax(gdat.cntpdata)
    minm = minmcntpdata
    maxm = maxmcntpdata
    setp_varb(gdat, 'cntpdata', minm=minm, maxm=maxm, lablroot='$C_{D}$', scal='asnh', strgmodl='plot')
    
    maxm = maxmcntpdata
    minm = 1e-1 * minmcntpdata
    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
        setp_varb(gdat, 'cntpmodl', minm=minm, maxm=maxm, strgmodl=strgmodl, scal='asnh')
    
    print('gdat.labltickmajrpara.cntpmodl')
    print(gdat.labltickmajrpara.cntpmodl)

    # residual limits
    maxm = np.ceil(maxmcntpdata * 0.1)
    minm = -np.ceil(maxmcntpdata * 0.1)
    setp_varb(gdat, 'cntpresi', minm=minm, maxm=maxm, lablroot='$C_{R}$', scal='asnh', strgmodl='plot')

    # 1-point function of the data counts
    for m in gdat.indxevtt:
        if gdat.numbpixl > 1:
            for i in gdat.indxener: 
                print('gdat.cntpdata[i, :, m]')
                summgene(gdat.cntpdata[i, :, m])
                print('gdat.binspara.cntpdata')
                summgene(gdat.binspara.cntpdata)
                histcntp = np.histogram(gdat.cntpdata[i, :, m], bins=gdat.binspara.cntpdata)[0]
                setattr(gdat, 'histcntpdataen%02devt%d' % (i, m), histcntp)
        else:
            histcntp = np.histogram(gdat.cntpdata[:, 0, m], bins=gdat.binspara.cntpdata)[0]
            setattr(gdat, 'histcntpdataevt%d' % m, histcntp)

    # obtain cartesian versions of the maps
    if gdat.typepixl == 'cart':
        ## data counts
        gdat.cntpdatacart = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt))
        gdat.cntpdatacart[:, gdat.indxpixlrofi, :] = gdat.cntpdata
        gdat.cntpdatacart = gdat.cntpdatacart.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
   

def retr_infodens(pdfnpost, pdfnprio):
    
    infodens = pdfnpost * np.log(pdfnpost / pdfnprio)

    return infodens


def retr_llik(gdat, strgmodl, cntpmodl):
   
    if gdat.liketype == 'pois':
        llik = gdat.cntpdata * np.log(cntpmodl) - cntpmodl
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.cntpdata - cntpmodl)**2 / gdat.varidata
     
    return llik


def retr_mapsgaus(gdat, lgal, bgal, spec, size, ellp, angl):
    
    rttrmatr = np.array([[np.cos(angl), -np.sin(angl)], [np.sin(angl), np.cos(angl)]])
    icovmatr = np.array([[1. / ((1. - ellp) * size)**2, 0.], [0., 1. / size**2]])

    posi = np.array([lgalgrid - lgal, bgalgrid - bgal])
    mapsgaus = flux * np.exp(-0.5 * np.sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / (1. - ellp)
        
    return mapsgaus


def retr_sbrtsers(gdat, lgalgrid, bgalgrid, lgal, bgal, spec, size, ellp, angl, seri=np.array([4.])):
   
    lgalrttr = (1. - ellp) * (np.cos(angl) * (lgalgrid - lgal) - np.sin(angl) * (bgalgrid - bgal))
    bgalrttr = np.sin(angl) * (lgalgrid - lgal) + np.cos(angl) * (bgalgrid - bgal) 
    angl = np.sqrt(lgalrttr**2 + bgalrttr**2)
    
    # interpolate pixel-convolved Sersic surface brightness
    if gdat.typesers == 'intp':

        shapinpt = angl.shape 
        inpt = np.empty(list(shapinpt) + [3])
        inpt[..., 0] = angl
        inpt[..., 1] = size
        inpt[..., 2] = seri
        
        sbrtsers = spec[:, None, None] * sp.interpolate.interpn((gdat.binspara.lgalsers, gdat.binspara.halfsers, gdat.binspara.indxsers), gdat.sersprof, inpt)[None, :, None]
    
    # evaluate directly de Vaucouleurs
    if gdat.typesers == 'vauc':
        sbrtsers = spec[:, None, None] * retr_sbrtsersnorm(angl, size)[None, :, None]
    
    return sbrtsers


def retr_sbrtsersnorm(angl, halfsers, indxsers=4.):

    ## this approximation works for 0.5  < indx < 10
    factsers = 1.9992 * indxsers - 0.3271
    
    ## surface brightness profile at the half-light radius for a 1 erg cm^-2 s^-1 A^-1 source
    if indxsers == 4.:
        sbrthalf = 1. / 7.2 / np.pi / halfsers**2
    else:
        sbrthalf = 1. / 2. / np.pi / np.exp(factsers) * factsers**(2 * indxsers) / indxsers / sp.special.gamma(2. * indxsers) / halfsers**2
                
    ## surface brightness profile
    sbrtsers = sbrthalf * np.exp(-factsers * ((angl / halfsers)**(1. / indxsers) - 1.))
    
    return sbrtsers


def copytdgu(varb):
    
    if isinstance(varb, np.ndarray):
        return np.copy(varb)
    else:
        return deepcopy(varb)


def proc_anim(rtag):
    
    pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
    
    print('Making animations of frame plots for %s...' % rtag)
    
    path = pathoutprtag + 'gdatinit'
    gdat = readfile(path)
    for strgpdfn in gdat.liststrgpdfn:
        for nameextn in gdat.liststrgfoldanim:
            
            pathframextn = gdat.pathimag + rtag + '/' + strgpdfn + '/fram/' + nameextn
            pathanimextn = gdat.pathimag + rtag + '/' + strgpdfn + '/anim/' + nameextn
        
            try:
                listfile = fnmatch.filter(os.listdir(pathframextn), '*_swep*.pdf')
            except:
                print('%s failed.' % pathframextn)
                continue
    
            listfiletemp = []
            for thisfile in listfile:
                listfiletemp.extend((thisfile.split('_')[0]).rsplit('/', 1))
            
            listname = list(set(listfiletemp))
            if len(listname) == 0:
                continue
            
            shuffle(listname)
    
            for name in listname:
                
                strgtemp = '%s*_swep*.pdf' % name
                listfile = fnmatch.filter(os.listdir(pathframextn), strgtemp)
                numbfile = len(listfile)
                liststrgextn = []
                for k in range(numbfile):
                    liststrgextn.append((listfile[k].split(name)[1]).split('_')[0])
                
                liststrgextn = list(set(liststrgextn))
                
                for k in range(len(liststrgextn)):
            
                    listfile = fnmatch.filter(os.listdir(pathframextn), name + liststrgextn[k] + '_swep*.pdf')
                    numbfile = len(listfile)
                    
                    indxfilelowr = 0
                    
                    if indxfilelowr < numbfile:
                        indxfileanim = np.arange(indxfilelowr, numbfile)
                    else:
                        continue
                        
                    indxfileanim = np.random.choice(indxfileanim, replace=False, size=indxfileanim.size)
                    
                    cmnd = 'convert -delay 20 -density 300 -quality 100 '
                    for n in range(indxfileanim.size):
                        cmnd += '%s%s ' % (pathframextn, listfile[indxfileanim[n]])
    
                    namegiff = '%s%s.gif' % (pathanimextn, name + liststrgextn[k])
                    cmnd += ' ' + namegiff
                    print('Processing %s' % namegiff)
                    if not os.path.exists(namegiff):
                        print('Run: %s, pdf: %s' % (rtag, strgpdfn))
                        print('Making %s animation...' % name)
                        os.system(cmnd)
                    else:
                        print('GIF already exists.')
                        pass
    
    pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
    filestat = open(pathoutprtag + 'stat.txt', 'a')
    filestat.write('animfinl written.\n')
    filestat.close()
    

def plot_samp(gdat, gdatmodi, strgstat, strgmodl, strgphas, strgpdfn='post', gdatmock=None, booltile=False):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    if not booltile:
    
        if strgstat != 'pdfn':
            numbelem = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                gmodstat.numbelem[l] = gmodstat.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int)
    
    if gdatmodi is not None:
        strgswep = '_%09d' % gdatmodi.cntrswep
    else:
        strgswep = ''
    
    if not booltile:
        # data count maps
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    if gdat.boolmakeframcent and (i != gdat.numbener / 2 or m != gdat.numbevtt / 2):
                        continue
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m)
            ## residual count maps
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    if gdat.boolmakeframcent and (i != gdat.numbener / 2 or m != gdat.numbevtt / 2):
                        continue
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpresi', i, m)
        
        if gdat.numbpixl > 1:
            if gmod.numbparaelem > 0:
                if gmod.boolelemlens:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelem', booltdim=True)
        
        # temp -- restrict other plots to indxmodlelemcomp
        if gdat.boolbinsener:
            for specconvunit in gdat.listspecconvunit:
                if not gmod.boolbfun:
                    plot_sbrt(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, specconvunit)
    
            if gmod.boolapplpsfn:
                plot_psfn(gdat, gdatmodi, strgstat, strgmodl)
    
    setp_indxswepsave(gdat)
    if gmod.numbparaelem > 0:
        # element parameter histograms
        if not (strgmodl == 'true' and gdat.typedata == 'inpt'):
            
            limtydat = gdat.limtydathistfeat

            for l in gmod.indxpopl:
                strgindxydat = 'pop%d' % l
                for nameparaderielemodim in gmod.namepara.derielemodim[l]:
                    if not (nameparaderielemodim == 'flux' or nameparaderielemodim == 'mcut' or \
                            nameparaderielemodim == 'deltllik' or nameparaderielemodim == 'defs' or nameparaderielemodim == 'nobj'):
                        continue
                                                                              
                    if gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt':
                        continue
                    indxydat = [l, slice(None)]
                    
                    name = nameparaderielemodim
                    namepopl = nameparaderielemodim + 'pop%d' % l
                    lablxdat = getattr(gmod.labltotlpara, namepopl)
                    scalxdat = getattr(gmod.scalpara, namepopl)
                    limtxdat = getattr(gmod.limtpara, namepopl)
                    meanxdat = getattr(gdat.meanpara, name)
                        
                    if gdat.numbpixl > 1:
                        listydattype = ['totl', 'sden']
                    else:
                        listydattype = ['totl']
                    for ydattype in listydattype:
                        
                        ## plot the surface density of elements
                        if ydattype == 'sden':
                            
                            # plot the surface density of elements only for the amplitude feature
                            if nameparaderielemodim != gmod.nameparagenrelemampl: 
                                continue
                            
                            if gdat.sdenunit == 'degr':
                                lablydat = r'$\Sigma_{%s}$ [deg$^{-2}$]' % gmod.lablelemextn[l]
                            if gdat.sdenunit == 'ster':
                                lablydat = r'$\Sigma_{%s}$ [sr$^{-2}$]' % gmod.lablelemextn[l]
                        
                        ## plot the total number of elements
                        if ydattype == 'totl':
                            lablydat = r'$N_{%s}$' % gmod.lablelemextn[l]
                    
                        if ydattype == 'totl' and not gdat.rtagmock is None:
                            listtypehist = ['hist', 'histcorrreca']
                        else:
                            listtypehist = ['hist']
                        
                        boolhistprio = not booltile
                        for typehist in listtypehist:
                            
                            if typehist == 'histcorrreca':
                                
                                if gmod.numbparaelem == 0 or gdat.priofactdoff == 0.:
                                    continue

                                if nameparaderielemodim == 'specplot' or nameparaderielemodim == 'spec' or nameparaderielemodim == 'deflprof':
                                    continue
                            
                                if not nameparaderielemodim in gmod.namepara.genrelem[l]:
                                    continue
                            
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'hist' + nameparaderielemodim + 'pop%d' % l, \
                                              'mean' + nameparaderielemodim, scalydat='logt', lablxdat=lablxdat, \
                                              lablydat=lablydat, histodim=True, ydattype=ydattype, \
                                              scalxdat=scalxdat, meanxdat=meanxdat, limtydat=limtydat, \
                                              limtxdat=limtxdat, boolhistprio=boolhistprio, \
                                              #indxydat=indxydat, strgindxydat=strgindxydat, \
                                              nameinte='histodim/', typehist=typehist)
    
    if not booltile:
        if gmod.numbparaelem > 0:
            # element parameter correlations
            for l in gmod.indxpopl:
                if strgmodl != 'true' and gdat.boolinforefr and gdat.boolasscrefr:
                    for strgfeat in gmod.namepara.derielemodim[l]:
                        if not (strgfeat == 'flux' or strgfeat == 'mass' or strgfeat == 'deltllik' or strgfeat == 'nobj') and \
                                                                                    (gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
                            continue
                        for q in gdat.indxrefr:
                            if not l in gdat.refrindxpoplassc[q]:
                                continue
                            if gdat.refr.numbelem[q] == 0:
                                continue
                            if not strgfeat in gdat.refr.namepara.elem[q] or strgfeat in gdat.refr.namepara.elemonly[q][l]:
                                continue
                            plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat)
                            plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat, plotdiff=True)
                    
        if not (gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # plots
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    if gmod.numbpopl > 1:
                        if gmod.numbparaelem > 0:
                            for l in gmod.indxpopl:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m, indxpoplplot=l)
        
            ## histograms of the number of counts per pixel
            limtxdat = [gdat.minmpara.cntpmodl, gdat.maxmpara.cntpmodl]
            for nameecom in gmod.listnameecomtotl:
                name = 'histcntp' + nameecom
                for m in gdat.indxevtt: 
                    for i in gdat.indxener:
                        if gdat.numbener > 1:
                            name += 'en%02d' % (i)
                        if gdat.numbevtt > 1:
                            name += 'evt%d' % (m)
                            
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                             name, 'meancntpdata', scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, histodim=True, \
                                             lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbener], limtxdat=limtxdat)

            
            ## highest amplitude element
            # temp
            if gmod.numbparaelem > 0:
                # completeness and false discovery rate
                if strgmodl != 'true' and gdat.boolasscrefr:
                    for strgclas in ['cmpl', 'fdis']:
                        nameinte = strgclas + 'odim/'
                        limtydat = [getattr(gdat, 'minm' + strgclas), getattr(gdat, 'maxm' + strgclas)]
                        for l in gmod.indxpopl:
                            for q in gdat.indxrefr:
                                if not l in gdat.refrindxpoplassc[q]:
                                    continue
                                if gdat.refr.numbelem[q] == 0 and strgclas == 'cmpl' or gmod.numbparaelem == 0 and strgclas == 'fdis':
                                    continue
                                if strgclas == 'cmpl':
                                    lablydat = getattr(gmod.lablpara, strgclas + 'pop%dpop%d' % (l, q))
                                    strgindxydat = 'pop%dpop%d' % (l, q)
                                else:
                                    lablydat = getattr(gmod.lablpara, strgclas + 'pop%dpop%d' % (q, l))
                                    strgindxydat = 'pop%dpop%d' % (q, l)
                                for strgfeat in gdat.refr.namepara.elem[q]:
                                    if strgfeat == 'etag':
                                        continue
                                    if strgclas == 'fdis' and not strgfeat in gmod.namepara.derielemodim[l]:
                                        continue
                                    if not strgfeat.startswith('spec') and not strgfeat.startswith('defl') \
                                                         and not strgfeat in gdat.refr.namepara.elemonly[q][l] and \
                                                         not (gdat.typedata == 'mock' and (strgfeat.endswith('pars') or strgfeat.endswith('nrel'))):
                                        
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgclas + strgfeat + strgindxydat, \
                                                  'mean' + strgfeat, lablxdat=lablxdat, \
                                                  lablydat=lablydat, \
                                                  #plottype='errr', \
                                                  scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                  omittrue=True, nameinte=nameinte)

            if gmod.numbparaelem > 0:
                alph = 0.1
                if strgmodl == 'true':
                    pathtemp = gdat.pathinit
                else:
                    if strgstat == 'this':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/fram/'
                    elif strgstat == 'mlik':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/finl/'
                    elif strgstat == 'pdfn':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/finl/'
                colr = retr_colr(gdat, strgstat, strgmodl, indxpopl=None)
        
                # transdimensional element parameters projected onto the data axes
                if not (strgstat == 'pdfn' and not gdat.boolcondcatl):
                    for l in gmod.indxpopl:
                        if gmod.typeelem[l] == 'lght':
                            # PS spectra
                            if strgstat == 'pdfn':
                                specplot = [np.empty((gdat.numbenerplot, gdat.numbstkscond))]
                                for r in gdat.indxstkscond:
                                    specplot[0][:, r] = gdat.dictglob['poststkscond'][r]['specplot'][0, :]
                            
                            listxdat = []
                            listplottype = []
                            
                            for k in range(specplot[l].shape[-1]):
                                listxdat.append(gdat.meanpara.enerplot)
                                listplottype.append('lghtline')
                            
                            for specconvunit in gdat.listspecconvunit:
                                listydat = []
                                
                                for k in range(specplot[l].shape[-1]):
                                    specplottemp = specplot[l]
                                    if strgmodl == 'true':
                                        specplottemp = np.copy(specplottemp[0, :, k])
                                    else:
                                        specplottemp = np.copy(specplottemp[:, k])
                                    if specconvunit[0] == 'en01':
                                        specplottemp *= gdat.meanpara.enerplot
                                    if specconvunit[0] == 'en02':
                                        specplottemp *= gdat.meanpara.enerplot**2
                                    if specconvunit[0] == 'en03':
                                        # temp
                                        pass
                                    listydat.append(specplottemp)
                                
                                lablydat = getattr(gmod.lablpara, 'flux' + specconvunit[0] + specconvunit[1] + 'totl')
                                strgtemp = specconvunit[0] + specconvunit[1]
                                if specconvunit[0] == 'en03':
                                    strgtemp += specconvunit[2]
                                path = pathtemp + strgstat + 'specpop%d%s%s.pdf' % (l,  strgtemp, strgswep)
                                limtydat = [np.amin(gdat.minmspec), np.amax(gdat.maxmspec)]
                                tdpy.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', \
                                                           lablxdat=gdat.lablenertotl, colr=colr, alph=alph, \
                                                           plottype=listplottype, limtxdat=[gdat.minmener, gdat.maxmener], lablydat=lablydat, \
                                                           limtydat=limtydat)
                    
                    if gmod.boollenssubh:

                        ## deflection profiles
                        if gdat.boolvariasca and gdat.boolvariacut:
                            lablxdat = gdat.labltotlpara.gang
                            if strgstat == 'pdfn':
                                deflprof = [np.empty((gdat.numbanglfull, gdat.numbstkscond))]
                                asca = [np.empty(gdat.numbstkscond)]
                                acut = [np.empty(gdat.numbstkscond)]
                                for r in gdat.indxstkscond:
                                    deflprof[0][:, r] = gdat.dictglob['poststkscond'][r]['deflprof'][0, :]
                                    asca[0][r] = gdat.dictglob['poststkscond'][r]['asca'][0]
                                    acut[0][r] = gdat.dictglob['poststkscond'][r]['acut'][0]

                            for l in range(len(deflprof)):
                                xdat = gdat.meanpara.anglfull * gdat.anglfact
                                listydat = []
                                listvlinfrst = []
                                listvlinseco = []
                                
                                if 'deflprof' in gmod.typeelem[l]:

                                    if strgmodl == 'true':
                                        deflproftemp = deflprof[l][0, :, :]
                                    else:
                                        deflproftemp = deflprof[l]
                                    
                                    for k in range(deflprof[l].shape[-1]):
                                        listydat.append(deflproftemp[:, k] * gdat.anglfact)
                                        if strgmodl == 'true':
                                            ascatemp = asca[l][0, k]
                                            acuttemp = acut[l][0, k]
                                        else:
                                            ascatemp = asca[l][k]
                                            acuttemp = acut[l][k]
                                        listvlinfrst.append(ascatemp * gdat.anglfact) 
                                        listvlinseco.append(acuttemp * gdat.anglfact)
                                    
                                    beinhost = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'paragenrscalfull', strgpdfn, indxvarb=gmod.indxpara.beinhost)
                                    listydat.append(xdat * 0. + gdat.anglfact * beinhost)
                                    
                                    path = pathtemp + strgstat + 'deflsubhpop%d%s.pdf' % (l, strgswep)
                                    limtydat = [1e-3, 1.]
                                    limtxdat = [1e-3, 1.]
                                    tdpy.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', \
                                                                        lablxdat=lablxdat, drawdiag=True, limtydat=limtydat, \
                                                                        limtxdat=limtxdat, colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', \
                                                                        listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
                        
                if gdat.typedata == 'mock':
                    # pulsar masses
                    for l in gmod.indxpopl:
                        if gmod.typeelem[l] == 'lghtpntspuls':
                            lablxdat = gdat.labltotlpara.gang
                            limtydat = [gdat.minmmassshel, gdat.maxmmassshel]
                            lablydat = gdat.lablmassshel
                            name = 'massshelpop%d' % l
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                           lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                    if gmod.boollens:
                        ## radial mass budget
                        lablxdat = gdat.lablanglfromhosttotl
                        for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                            
                            # host mass
                            for e in gmod.indxsersfgrd:
                                strgsersfgrd = 'isf%d' % e
                                limtydat = [gdat.minmmcut, getattr(gdat, 'plotmaxmmasshost' + strgsersfgrd + strgcalcmasssubh + 'bein')]
                                lablydat = getattr(gmod.lablpara, 'masshost' + strgsersfgrd + strgcalcmasssubh + 'totl')
                                name = 'masshost%s%s' % (strgsersfgrd, strgcalcmasssubh)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                          lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)
                            
                            if gmod.boolelemdeflsubhanyy:
                                # subhalo masses
                                limtydat = [gdat.minmmcut, getattr(gdat, 'plotmaxmmasssubh' + strgcalcmasssubh + 'bein')]
                                lablydat = getattr(gmod.lablpara, 'masssubh' + strgcalcmasssubh + 'totl')
                                name = 'masssubh%s' % (strgcalcmasssubh)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                          lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                                # subhalo mass fraction
                                limtydat = [1e-3, 0.1]
                                lablydat = getattr(gmod.lablpara, 'fracsubh' + strgcalcmasssubh + 'totl')
                                name = 'fracsubh%s' % (strgcalcmasssubh)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                            lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                alph = 0.1

                if gdat.boolmodipsfn and gmod.boolelempsfnanyy:
                    ## PSF radial profile
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            indxydat = [i, slice(None), m]
                            strgindxydat = 'en%02devt%d' % (i, m)
                            lablxdat = gdat.labltotlpara.gang
                            limtydat= np.array([1e-3, 1e3]) * gdat.anglfact**2
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psfn', \
                                                            'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalydat='logt', \
                                                            lablxdat=lablxdat, lablydat=r'$\mathcal{P}$', limtydat=limtydat)
        
                # internally and externally corrected element parameter histograms
                if gdat.typedata == 'inpt' and strgstat == 'pdfn' and gdat.rtagmock is not None:
                    limtydat = gdat.limtydathistfeat
                    for l in gmod.indxpopl:
                        strgindxydat = 'pop%d' % l
                        for strgfeat in gmod.namepara.derielemodim[l]:
                            if strgfeat.startswith('aerr') or strgfeat == 'specplot' or strgfeat == 'spec' or strgfeat == 'deflprof':
                                continue
                            lablydat = r'$N_{%s}$' % gmod.lablelemextn[l]
                            for namecorr in ['incr', 'excr']:
                                nameinte = namecorr + 'odim/'
                                for qq in gdatmock.indxrefr:
                                    if namecorr == 'excr':
                                        if not strgfeat in gmod.namepara.extrelem[l]:
                                            continue
                                        q = gdat.listnamerefr.index(strgfeat[-4:])
                                        if getattr(gdat, 'crex' + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)) is None:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                          lablydat=lablydat, histodim=True, ydattype='totl', \
                                                          scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                          nameinte=nameinte)
               
                                    else:
                                        if strgfeat in gmod.namepara.extrelem[l]:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%d' % (qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                      lablydat=lablydat, histodim=True, ydattype='totl', \
                                                      scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                      nameinte=nameinte)


    if not (gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
        if gmod.numbparaelem > 0:
            # element parameter correlations
            liststrgelemtdimvarb = getattr(gdat, 'liststrgelemtdimvarb' + strgphas)
            for strgelemtdimtype in gdat.liststrgelemtdimtype:
                for strgelemtdimvarb in liststrgelemtdimvarb:
                    if strgelemtdimvarb.startswith('cmpl'):
                        continue
                    for l0 in gmod.indxpopl:
                        for strgfrst in gmod.namepara.genrelem[l0]:
                            
                            if strgfrst.startswith('spec') or strgfrst == 'specplot' or strgfrst == 'deflprof':
                                continue

                            for strgseco in gmod.namepara.genrelem[l0]:
                                
                                if strgseco.startswith('spec') or strgseco == 'specplot' or strgseco == 'deflprof':
                                    continue
                                
                                if not checstrgfeat(strgfrst, strgseco):
                                    continue
                                    
                                if strgelemtdimvarb.startswith('hist'):
                                    
                                    strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%d' % l0
                                    plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                               l0, strgfrst + 'pop%d' % l0, \
                                                                                   strgseco + 'pop%d' % l0, \
                                                                                   strgtotl, strgpdfn=strgpdfn)
                                else:
                                    if booltile:
                                        continue

                                    if strgfrst.startswith('aerr') or strgseco.startswith('aerr'):
                                        continue
                                    if strgelemtdimvarb.startswith('fdis'):
                                        for q in gdat.indxrefr:
                                            strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%d' % (q, l0)
                                            plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                            l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
                                    elif strgelemtdimvarb.startswith('excr') or strgelemtdimvarb.startswith('incr'):
                                        for qq in gdatmock.indxrefr:
                                            if strgelemtdimvarb.startswith('excr'):
                                                for q in gdat.indxrefr:
                                                    if getattr(gdat, 'crex' + strgfrst + strgseco + 'pop%dpop%dpop%d' % (q, qq, l0)) is None:
                                                        continue
                                                    strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%dpop%d' % (q, qq, l0)
                                                    plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                                l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
                                            else:
                                                if strgfrst[-4:] in gdat.listnamerefr and strgseco[-4:] in gdat.listnamerefr:
                                                    continue
                                                strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%d' % (qq, l0)
                                                plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                            l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
        
            if not (gdat.typedata == 'mock' and (gmod.numbelemtotl == 0 or gmod.maxmpara.numbelemtotl == 0)):
                

                for q in gdat.indxrefr:
                    
                    if strgphas == 'init' and gdat.typedata == 'mock':
                        continue
                    print('strgpdfn')
                    print(strgpdfn)
                    raise Exception('')

                    if booltile:
                        continue
                    for l0 in gmod.indxpopl:
                        for refrstrgfrst in gdat.refr.namepara.elem[q]:
                            if refrstrgfrst == 'spec' or refrstrgfrst == 'specplot' or refrstrgfrst == 'deflprof' or refrstrgfrst == 'etag':
                                continue
                            if refrstrgfrst in gdat.refr.namepara.elemonly[q][l0]:
                                continue
                            for refrstrgseco in gdat.refr.namepara.elem[q]:
                                if refrstrgseco in gdat.refr.namepara.elemonly[q][l0]:
                                    continue
                                if refrstrgseco == 'spec' or refrstrgseco == 'specplot' or refrstrgseco == 'deflprof' or refrstrgseco == 'etag':
                                    continue
                                
                                if not checstrgfeat(refrstrgfrst, refrstrgseco):
                                    continue
                                        
                                if refrstrgfrst.startswith('aerr') or refrstrgseco.startswith('aerr') or refrstrgfrst == 'specplot' or refrstrgseco == 'specplot':
                                    continue
                                
                                strgtotl = 'cmpl' + refrstrgfrst + refrstrgseco + 'pop%dpop%d' % (l0, q)
                                
                                plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, 'bind', 'cmpl', \
                                                        q, refrstrgfrst + 'pop%d' % l0, refrstrgseco + 'pop%d' % l0, strgtotl, strgpdfn=strgpdfn)
            
    if not booltile:
        if not (gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # data and model count scatter
            for m in gdat.indxevttplot:
                if gdat.numbpixl > 1:
                    for i in gdat.indxener:
                        plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m, indxenerplot=i)
                else:
                    plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m)

            ## spatial priors
            # temp
            if gdat.numbpixl > 1:
                if gmod.numbparaelem > 0:
                    for l in gmod.indxpopl:
                        for strgfeat, strgpdfn in zip(gmod.namepara.genrelemmodu[l], gmod.liststrgpdfnmodu[l]):
                            if strgpdfn == 'tmplreln':
                                plot_genemaps(gdat, gdatmodi, 'fitt', strgpdfn, 'lpdfspatpriointp', booltdim=True)
                            if strgpdfn == 'tmplgaum':
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'lpdfspatpriointp', booltdim=True)
            
            # model count maps
            ## backgrounds
            if gdat.numbpixl > 1:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        for c in gmod.indxback:
                            if gmod.boolbfun:
                                continue
                            if not gmod.boolunifback[c]:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpback%04d' % c, i, m, strgcbar='cntpdata')
                
                ## count error
                if strgmodl != 'true':
                    if gmod.numbparaelem > 0:
                        for l in gmod.indxpopl:
                            if gmod.boolcalcerrr[l]:
                                for i in gdat.indxener:
                                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntperrr', i, -1, strgcbar='cntpresi')
                
                ## diffuse components 
                for i in gdat.indxener:
                    for k, name in enumerate(gmod.listnamediff):
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntp%s' % (name), i, strgcbar='cntpdata')
            
                ## model count maps
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpmodl', i, m, strgcbar='cntpdata')
            
                # likelihood
                if strgmodl != 'true':
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'llik', i, m, strgcbar='llikmaps')
                
                if gmod.boollens:
                    ## lensing signal to noise
                    if strgmodl == 'true':
                        for i in gdat.indxener:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 's2nr', i, -1)
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magn', booltdim=True)
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'conv', booltdim=True)
                    for i in gdat.indxener:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplens', i, strgcbar='cntpdata', booltdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplensgradmgtd', i, strgcbar='cntpdata', booltdim=True)
            
            if gdat.penalpridiff:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                                'psecodimdatapntsen%02devt%d' % (i, m), 'meanmpolodim', lablxdat='$l$', lablydat='$P_{resi}(l)$', \
                                                                                                                 limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psecodimdatapntsprioen%02devt%d' % (i, m), 'meanmpolodim', lablxdat='$l$', \
                                                                                           lablydat='$P_{prio}(l)$', limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                
            if gmod.boollens:
                indxydat = [slice(None)]
                strgindxydat = ''
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P(k)$', limtydat=[1e-1, 1e2], \
                                                                          scalxdat='logt', scalydat='logt', indxydat=indxydat, strgindxydat=strgindxydat)
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdefl', 'meandefl', \
                                                                        scal='self', lablxdat=r'$\alpha$ [arcsec]', lablydat=r'$N_{pix}$', \
                                                                                 strgindxydat=strgindxydat, indxydat=indxydat, histodim=True)
            if gmod.numbparaelem > 0 and gmod.boolelemdeflsubhanyy:
                indxydat = [slice(None)]
                strgindxydat = ''
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecelemodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P_{sub}(k)$', \
                                       strgindxydat=strgindxydat, indxydat=indxydat, limtydat=[1e-5, 1e-1], scalxdat='logt', scalydat='logt')
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdeflsubh', 'meandeflsubh', scal='self', lablxdat=r'$\alpha$ [arcsec]', \
                                       strgindxydat=strgindxydat, indxydat=indxydat, lablydat=r'$N_{pix}$', histodim=True)
            
            if gmod.boollens:
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrd', i, -1, strgcbar='cntpdata')
                    if gmod.numbparaelem > 0 and gmod.boolelemsbrtextsbgrdanyy:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdgalx', i, -1, strgcbar='cntpdata')
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdexts', i, -1, strgcbar='cntpdata')
                
                # gradient of the lens emission
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplensgrad', indxenerplot=i, indxevttplot=m)
                
        if not (gdat.boolshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            if gmod.boollens:
                # overall deflection field
                plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, multfact=0.1)
                
                # deflection field due to individual lenses
                for k in range(numbdeflsingplot):  
                    if k == 0:
                        multfact = 0.1
                    elif k == 1:
                        multfact = 1.
                    elif k >= 2:
                        multfact = 10.
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxdefl=k, multfact=multfact)
                
                # residual deflection field
                if strgmodl == 'fitt' and gdat.typedata == 'mock':
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, nameparagenrelem='resi', multfact=100.)
                    if strgstat != 'pdfn':
                        for k in range(numbsingcomm):
                            plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, nameparagenrelem='resi', indxdefl=k, multfact=100.)
                    
                    if gdat.numbpixl > 1:
                        if gmod.numbparaelem > 0:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresi', booltdim=True)
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresiperc', booltdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresi', booltdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresiperc', booltdim=True)
    

def dele_rtag(rtag):
    
    pathdata = pathpcat + '/data/outp/'
    pathimag = pathpcat + '/imag/'
    
    cmnd = 'rm -rf %s%s' % (pathdata, rtag)
    print(cmnd)
    os.system(cmnd)
    cmnd = 'rm -rf %s%s' % (pathimag, rtag)
    os.system(cmnd)
    print(cmnd)


def plot_infopvks(gdat, gdatprio, name, namefull, nameseco=None):
    
    pvks = getattr(gdat, 'pvks' + namefull)

    info = getattr(gdat, 'info' + namefull)

    path = gdat.pathinfo + 'info' + namefull

    if nameseco is not None:
       
        indxpoplfrst = int(namefull[-1])
        
        # information gain
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.pcolor(varbfrst, varbseco, info, cmap='Greys')
        plt.colorbar(imag)
        plot_sigmcont(gdat.fitt, '', axis, name, indxpoplfrst, strgseco=nameseco)
        if scalfrst == 'logt':
            axis.set_xscale('log')
        if scalseco == 'logt':
            axis.set_yscale('log')
        axis.set_xlabel(getattr(gdat.labltotlpara, name))
        axis.set_ylabel(getattr(gdat.labltotlpara, nameseco))
        axis.set_xlim(limtfrst)
        axis.set_ylim(limtseco)
        plt.tight_layout()
        plt.savefig(path)
        plt.close(figr)

        # KS test p value
        pathpvkstdim = gdat.pathinfo + 'pvks' + namefull
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.pcolor(varbfrst, varbseco, pvks, cmap='Greys')
        plt.colorbar(imag)
        plot_sigmcont(gdat.fitt, '', axis, name, indxpoplfrst, strgseco=nameseco)
        if scalfrst == 'logt':
            axis.set_xscale('log')
        if scalseco == 'logt':
            axis.set_yscale('log')
        axis.set_xlabel(getattr(gdat.labltotlpara, name))
        axis.set_ylabel(getattr(gdat.labltotlpara, nameseco))
        axis.set_xlim(limtfrst)
        axis.set_ylim(limtseco)
        plt.tight_layout()
        plt.savefig(pathpvkstdim)
        plt.close(figr)

    elif name != namefull:
        
        lablydat = '$D_{KL}$'
        lablxdat = getattr(gmod.lablpara, name + 'totl')
        xdat = getattr(gdat, 'mean' + name)
        ydat = getattr(gdat, 'info' + namefull)
        tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, scal)
        
        ydat = getattr(gdat, 'pvks' + namefull)
        pathpvks = gdat.pathinfo + 'pvks' + namefull
        tdpy.mcmc.plot_plot(pathpvks, xdat, ydat, lablxdat, '$p_{KS}$', scal)
        
    else:
        # horizontal axis
        xdat = getattr(gdat, 'mean' + name)
        lablxdat = getattr(gmod.lablpara, name + 'totl')
        
        # scaling
        scal = getattr(gdat, 'scal' + name) 
        
        # common title
        titl = '$D_{KL} = %.3g$, KS = %.3g $\sigma$' % (info, pvks)

        # DKL density
        pathdinf = gdat.pathinfo + 'dinf' + namefull
        ydat = getattr(gdat, 'infodens' + namefull)
        lablydat = r'$\rho_{D_{KL}}$'
        tdpy.mcmc.plot_plot(pathdinf, xdat, ydat, lablxdat, lablydat, scal, titl=titl)
        
        # prior and posterior PDFs
        pathpdfn = gdat.pathinfo + 'pdfn' + namefull
        lablydat = r'$P$'
        ydat = [getattr(gdat, 'pdfnpost' + namefull), getattr(gdatprio, 'pdfnprio' + namefull)]
        legd = ['$P$(%s|$D$)' % lablxdat, '$P$(%s)' % lablxdat]
        tdpy.mcmc.plot_plot(pathpdfn, xdat, ydat, lablxdat, lablydat, scal, colr=['k', 'k'], linestyl=['-', '--'], legd=legd, titl=titl)


def plot_finl(gdat=None, gdatprio=None, rtag=None, strgpdfn='post', gdatmock=None, booltile=None):
    
    if gdat.typeverb > 0:
        print('plot_finl()')
        print('Producing postprocessing plots...')

    timetotlinit = gdat.functime()
    
    gdat.strgbest = 'ML'
    
    if not booltile:
        # terms in the log-acceptance probability
        listindxsamptotlproptotl = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlproptotl')
        listindxsamptotlpropaccp = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlpropaccp')
        listindxsamptotlpropreje = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlpropreje')
        for n in gdat.indxproptype:
            pathbase = getattr(gdat, 'path' + strgpdfn + 'finl%s' % gdat.nameproptype[n])
            for k in gdat.indxtermlacp:
                varb = getattr(gdat, 'list' + strgpdfn + gdat.listnametermlacp[k])
                labl = gdat.listlabltermlacp[k]
                
                if listindxsamptotlproptotl[n].size > 0 and (varb[listindxsamptotlproptotl[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'totl'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlproptotl[n]], labl, titl=gdat.nameproptype[n] + ', Total')
                
                if listindxsamptotlpropaccp[n].size > 0 and (varb[listindxsamptotlpropaccp[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'accp'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropaccp[n]], labl, titl=gdat.nameproptype[n] + ', Accepted')
                
                if listindxsamptotlpropreje[n].size > 0 and (varb[listindxsamptotlpropreje[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'reje'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropreje[n]], labl, titl=gdat.nameproptype[n] + ', Rejected')
            
        if gdat.checprio and strgpdfn == 'post' and not booltile:
            # this works only for scalar variables -- needs to be generalized to all variables
            if gdatprio is None:
                pathoutprtag = retr_pathoutprtag(pathpcat, rtag)
                path = pathoutprtag + 'gdatfinlprio'
                gdatprio = readfile(path)

            for namevarbscal in gmod.namepara.scal:
                plot_infopvks(gdat, gdatprio, namevarbscal, namevarbscal)
            for l in gmod.indxpopl:
                for strgfeatfrst in gmod.namepara.genrelem[l]:
                    if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                        continue
                    plot_infopvks(gdat, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l)
                    for strgfeatseco in gmod.namepara.genrelem[l]:
                        if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                            continue
                        
                        if not checstrgfeat(strgfeatfrst, strgfeatseco):
                            continue
                                        
                        plot_infopvks(gdat, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l, nameseco=strgfeatseco)
        
        listparagenrscalfull = getattr(gdat, 'list' + strgpdfn + 'paragenrscalfull')
        listparagenrscalfull = getattr(gdat, 'list' + strgpdfn + 'paragenrscalfull')
        listparagenrscalbase = getattr(gdat, 'list' + strgpdfn + 'paragenrscalbase')
    
        listboolpropfilt = getattr(gdat, 'list' + strgpdfn + 'boolpropfilt')
        listmemoresi = getattr(gdat, 'list' + strgpdfn + 'memoresi')
        listindxproptype = getattr(gdat, 'list' + strgpdfn + 'indxproptype')
        listsampproc = getattr(gdat, 'list' + strgpdfn + 'sampproc')
    
        # Gelman-Rubin test
        pathdiag = getattr(gdat, 'path' + strgpdfn + 'finldiag')
        if gdat.numbproc > 1:
            if np.isfinite(gdat.gmrbstat).all():
                if gdat.typeverb > 0:
                    print('Gelman-Rubin TS...')
        
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                minm = min(np.amin(gdat.gmrbstat), np.amin(gdat.gmrbparagenrscalbase))
                maxm = max(np.amax(gdat.gmrbstat), np.amax(gdat.gmrbparagenrscalbase))
                bins = np.linspace(minm, maxm, 40)
                axis.hist(gdat.gmrbstat.flatten(), bins=bins, label='Data proj.')
                axis.hist(gdat.gmrbparagenrscalbase, bins=bins, label='Fixed dim.')
                axis.set_xlabel('PSRF')
                axis.set_ylabel('$N_{stat}$')
                plt.tight_layout()
                figr.savefig(pathdiag + 'gmrbhist.pdf')
                plt.close(figr)
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.plot(gmod.indxparagenrbase, gdat.gmrbparagenrscalbase)
                axis.set_xticklabels(gmod.labltotlpara.genrbase)
                axis.set_ylabel('PSRF')
                plt.tight_layout()
                figr.savefig(pathdiag + 'gmrbparagenrscalbase.pdf')
                plt.close(figr)
                
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        maps = gdat.gmrbstat[i, :, m]
                        path = pathdiag + 'gmrbdataen%02devt%d.pdf' % (i, m)
                        tdpy.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, typepixl=gdat.typepixl, \
                                                                   minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                   minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
            else:
                print('Inappropriate Gelman-Rubin test statistics encountered.')
    
        # plot autocorrelation
        if gdat.typeverb > 0:
            print('Autocorrelation...')
        tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrcntp[0, 0, 0, 0, :], gdat.timeatcrcntp[0, 0, 0, 0], strgextn='cntp')
        tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrpara[0, 0, :], gdat.timeatcrpara[0, 0], strgextn='para')
        print('Autocorrelation times:')
        for k, namepara in enumerate(gmod.namepara):
            print('%s %g' % (namepara, np.mean(gdat.timeatcrpara[:, k])))
        
        # plot proposal efficiency
        if gdat.typeverb > 0:
            print('Acceptance ratio...')
        numbtimemcmc = 20
        binstimemcmc = np.linspace(0., gdat.numbswep, numbtimemcmc)
        numbtick = 2
        sizefigrydat = 4. * gdat.numbproptype
        figr, axgr = plt.subplots(gdat.numbproptype, 1, figsize=(12., sizefigrydat), sharex='all')
        if gdat.numbproptype == 1:
            axgr = [axgr]
        for n, axis in enumerate(axgr):
            histtotl = axis.hist(listindxsamptotlproptotl[n], bins=binstimemcmc)[0]
            histaccp = axis.hist(listindxsamptotlpropaccp[n], bins=binstimemcmc)[0]
            axis.set_ylabel('%s' % gdat.nameproptype[n])
            if k == gdat.numbproptype - 1:
                axis.set_xlabel('$i_{samp}$')
        plt.tight_layout()
        figr.savefig(pathdiag + 'accpratiproptype.pdf')
        plt.close(figr)
   
        if gdat.typeverb > 0:
            print('Proposal execution times...')
        
        ## time performance
        #listchro = np.empty((gdat.numbswep, gdat.numbchro))
        #listchro = []
        #for k, name in enumerate(gdat.listnamechro):
        #    #listchro[:, k] = getattr(gdat, 'list' + strgpdfn + 'chro' + name).flatten() * 1e3
        #    listchro.append(getattr(gdat, 'list' + strgpdfn + 'chro' + name).flatten() * 1e3)
        #pathdiag = getattr(gdat, 'path' + strgpdfn + 'finldiag')
        #figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
        #axis.violin(listchro)
        #axis.set_yscale('log')
        #axis.set_ylabel('$t$ [ms]')
        #axis.set_xticklabels(gdat.listlablchro)
        #axis.axvline(mean(chro), ls='--', alpha=0.2, color='black')
        #figr.savefig(pathdiag + 'chro.pdf' % gdat.listnamechro[k])
        #plt.close(figr)

    # temp
    gdat.lablpmea = 'Mean'

    # posterior versions of the frame plots
    plot_samp(gdat, None, 'pdfn', 'fitt', 'finl', strgpdfn=strgpdfn, gdatmock=gdatmock, booltile=booltile)
   
    if booltile:
        return

    if gmod.numbparaelem > 0:
        if gdat.typeverb > 0:
            print('A mosaic of samples...')
    
        ## mosaic of images of posterior catalogs
        if gdat.numbpixl > 1:
            plot_mosa(gdat, strgpdfn)
    
    ## randomly selected trandimensional parameters
    if gmod.numbparaelem > 0:
        if gdat.typeverb > 0:
            print('Transdimensional parameters...')
    
        # choose the parameters based on persistence
        stdvlistsamptran = np.std(listparagenrscalfull[:, gmod.indxsamptrap], axis=0)
        indxtrapgood = np.where(stdvlistsamptran > 0.)[0]
        gmod.numbparaelemgood = indxtrapgood.size
        gmod.numbparaelemplot = min(3, gmod.numbparaelemgood)
        if gmod.numbparaelemplot > 0:
            indxtrapplot = np.sort(np.random.choice(gmod.indxsamptrap[indxtrapgood], size=gmod.numbparaelemplot, replace=False))

            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listelemfrst', listparagenrscalfull[:, gmod.indxsamptrap[:3]], [gmod.lablpara[k] for k in gmod.indxsamptrap[:3]])
            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listsamp', listparagenrscalfull[:, indxtrapplot], ['%d' % k for k in indxtrapplot])
            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listsamp', listparagenrscalfull[:, indxtrapplot], [gmod.lablpara[k] for k in indxtrapplot])
    
    if gdat.typeverb > 0:
        print('Scalar variables...')
    # scalar variables
    ## trace and marginal distribution of each parameter
    for name in gmod.namepara.scal:
        
        if gdat.typeverb > 0:
            print('Working on %s...' % name)
        scal = getattr(gdat, 'scal' + name) 
        corr = getattr(gdat, 'corr' + name)
        if corr is None:
            truepara = None
        else:
            truepara = getattr(gdat, 'corr' + name)
        
        listvarb = getattr(gdat, 'list' + strgpdfn + name)
        if listvarb.ndim != 1:
            if listvarb.shape[1] == 1:
                listvarb = listvarb[:, 0]
            else:
                raise Exception('')
        
        mlik = getattr(gdat, 'mlik' + name)
        path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscaltrac') + name
        tdpy.mcmc.plot_trac(path, listvarb, labltotl, truepara=truepara, scalpara=scal, listvarbdraw=[mlik], listlabldraw=[''], listcolrdraw=['r'])
        path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalhist') + name
        tdpy.mcmc.plot_hist(path, listvarb, labltotl, truepara=truepara, scalpara=scal, listvarbdraw=[mlik], listlabldraw=[''], listcolrdraw=['r'])
       
        for nameseco in gmod.namepara.scal:
            
            if name == nameseco:
                continue
            
            if gdat.typeverb > 0:
                print('Working on correlation of %s with %s...' % (name, nameseco))
            
            pathjoin = getattr(gdat, 'path' + strgpdfn + 'finlvarbscaljoin')
            if corrseco is None:
                trueparaseco = None
            else:
                trueparaseco = getattr(gdat, 'corr' + nameseco)
            
            if listvarbseco.ndim != 1:
                if listvarbseco.shape[1] == 1:
                    listvarbseco = listvarbseco[:, 0]
                else:
                    raise Exception('')
                
            listjoin = np.vstack((listvarb, listvarbseco)).T
    
            tdpy.mcmc.plot_grid(pathjoin, name + nameseco, listjoin, [labltotl, labltotlseco], scalpara=[scal, scalseco], truepara=[truepara, trueparaseco], \
                                                                                                join=True, listvarbdraw=[np.array([mlik, mlikseco])])

    if gdat.typeverb > 0:
        print('Fixed dimensional parameter covariance...')
    
    ### covariance
    ## overall
    path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
    truepara = gmod.corrparagenrscalbase
    mlikpara = gdat.mlikparagenrscalbase
    tdpy.mcmc.plot_grid(path, 'paragenrscalbase', listparagenrscalbase, gmod.labltotlpara.genrbasetotl, truepara=truepara, listvarbdraw=[mlikpara])
    
    # stacked posteiors binned in position and flux
    if gmod.numbparaelem > 0 and gdat.numbpixl > 1:
        liststrgbins = ['quad', 'full']
        for l in gmod.indxpopl:
            plot_histlgalbgalelemstkd(gdat, strgpdfn, l, 'cumu')
            for strgbins in liststrgbins:
                plot_histlgalbgalelemstkd(gdat, strgpdfn, l, strgbins, namepara.elemsign[l])

    if gdat.typeverb > 0:
        print('Prior and likelihood...')
    
    for strgpdfntemp in ['lpritotl', 'lliktotl']:

        if strgpdfntemp == 'lpritotl':
            labltemp = '\ln P(M)'
        if strgpdfntemp == 'lliktotl':
            labltemp = '\ln P(D|M)'
        labl = r'$%s$' % labltemp

        path = getattr(gdat, 'path' + strgpdfn + 'finl') + strgpdfntemp
        
        varb = getattr(gdat, 'list' + strgpdfn + strgpdfntemp)
        tdpy.mcmc.plot_hist(path, varb, labl)
        listvarbdraw = []
        listlabldraw = []
        listcolrdraw = []
        if gdat.typedata == 'mock':
            listvarbdraw += [getattr(gdat.true, strgpdfntemp)]
            listlabldraw += ['True model']
            listcolrdraw += [gdat.refr.colr]
        
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'list' + strgpdfn + strgpdfntemp), labl, \
                                listvarbdraw=listvarbdraw, listlabldraw=listlabldraw, listcolrdraw=listcolrdraw)
    
    # plot resident memory
    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxswep, np.mean(listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(pathdiag + 'memoresi.pdf')
    plt.close(figr)

    timetotlfinl = gdat.functime()
    if gdat.typeverb > 0:
        print('Plots and animations are produced in %.3g seconds.' % (timetotlfinl - timetotlinit))


def plot_sbrt(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, specconvunit):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        
        # plot reference spectra
        if gdat.listprefsbrtlabltotl is not None:
            for k in range(len(gdat.listprefsbrtlabltotl)):
                if gdat.listprefsbrttype[k] == 'shad':
                    factenerrefr = [[] for a in range(3)]
                    for a in range(3):
                        factenerrefr[a] = retr_factener(specconvunit[0], gdat.listprefsbrtener[k][a])
                    axis.plot(gdat.listprefsbrtener[k][0], gdat.listprefsbrtsbrt[k][0] * factenerrefr[0], color='m', label=gdat.listprefsbrtlabltotl[k])
                    enerpoly = np.empty(gdat.listprefsbrtener[k][1].size + gdat.listprefsbrtener[k][2].size)
                    enerpoly[:gdat.listprefsbrtener[k][1].size] = gdat.listprefsbrtener[k][1]
                    enerpoly[gdat.listprefsbrtener[k][1].size:] = gdat.listprefsbrtener[k][2][::-1]
                    sbrtpoly = np.empty(gdat.listprefsbrtener[k][1].size + gdat.listprefsbrtener[k][2].size)
                    sbrtpoly[:gdat.listprefsbrtener[k][1].size] = gdat.listprefsbrtsbrt[k][1] * factenerrefr[1]
                    sbrtpoly[gdat.listprefsbrtener[k][1].size:] = gdat.listprefsbrtsbrt[k][2][::-1] * factenerrefr[2][::-1]
                    axis.fill(enerpoly, sbrtpoly, color='m', alpha=0.5)
                else:
                    factenerrefr = retr_factener(specconvunit[0], gdat.listprefsbrtener[k][1])
                    axis.errorbar(gdat.listprefsbrtener[k][1], gdat.listprefsbrtsbrt[k][1] * factenerrefr, label=gdat.listprefsbrtlabltotl[k], color='m')
        
        if strgmodl == 'true':
            liststrgmodl = [strgmodl]
            listgdatobjt = [gdat]
        if strgmodl == 'fitt' and (strgstat == 'this' or strgstat == 'pdfn'):
            if gdat.typedata == 'mock':
                liststrgmodl = [strgmodl, 'true']
                listgdatobjt = [gdatobjt, gdat]
            else:
                liststrgmodl = [strgmodl]
                listgdatobjt = [gdatobjt]
        numbstrgstattemp = len(liststrgmodl)
        for a in range(numbstrgstattemp):
            
            indxploteleminit = []
            indxplotelemendd = []
                
            # number of transdimensional elements to be overplotted
            numbelemtemp = 0
            
            if gdat.numbpixl == 1 and strgstat != 'pdfn':
                if liststrgmodl[a] == 'fitt':
                    numbelem = [[] for l in gmod.indxpopl]
                    for l in gmod.indxpopl:
                        gmodstat.numbelem[l] = gmodstat.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int)
                        numbelemtemp += np.sum(gmodstat.numbelem[l])
                else:
                    for q in gdat.indxrefr:
                        numbelemtemp += np.sum(gdat.refr.numbelem[q])
                
            numbplot = numblablsbrtspec + numbelemtemp
            listydat = np.zeros((numbplot, gdat.numbener))
            listyerr = np.zeros((2, numbplot, gdat.numbener))
            
            cntr = 0
            cntrdata = cntr

            ## data
            listydat[cntr, :] = gdat.sbrtdatamean[b]
            listyerr[:, cntr, :] = gdat.sbrtdatastdv[b]
            cntr += 1
            
            for c in gmod.indxback:
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtback%04dmea%d' % (c, b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtback%04dmea%d' % (c, b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if gmod.numbparaelem > 0 and gmod.boolelemsbrtdfncanyy and not (liststrgmodl[a] == 'true' and gdat.refr.numbelemtotl == 0):
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncsubtmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncsubtmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if gmod.typeemishost != 'none':
                for e in gmod.indxsersfgrd:
                    listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrthostisf%dmea%d' % (e, b), strgpdfn)
                    if strgstat == 'pdfn':
                        listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], \
                                                                        'sbrthostisf%dmea%d' % (e, b), strgpdfn, strgmome='errr')
                    cntr += 1
            
            if gmod.boollens:
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtlensmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtlensmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if gdat.numbpixl == 1 and strgstat != 'pdfn':
                cntrline = cntr
                indxploteleminit.append(cntr)
                for l in gmod.indxpopl:
                    if liststrgmodl[a] == 'true':
                        for k in range(gmod.numbelem[l]):
                            listydat[cntr, :] = getattr(listgdatobjt[a], liststrgmodl[a] + 'spec')[l][0, :, k]
                            
                            if cntr == cntrline:
                                listlablsbrtspec = listlablsbrtspec[:cntr] + ['Lines'] + listlablsbrtspec[cntr:]
                            else:
                                listlablsbrtspec = listlablsbrtspec[:cntr] + [None] + listlablsbrtspec[cntr:]
                            
                            cntr += 1
                            if k == gmod.numbelem[l] - 1:
                                indxplotelemendd.append(k)
                    else:   
                        for k in range(gmodstat.numbelem[l]):
                            listydat[cntr, :] = getattr(listgdatobjt[a], strgstat + 'spec')[l][:, k]
                            
                            if cntr == cntrline:
                                listlablsbrtspec = listlablsbrtspec[:cntr] + ['Lines'] + listlablsbrtspec[cntr:]
                            else:
                                listlablsbrtspec = listlablsbrtspec[:cntr] + [None] + listlablsbrtspec[cntr:]
                
                            cntr += 1
                            if k == gmodstat.numbelem[l] - 1:
                                indxplotelemendd.append(k)
            ## total model
            if numblablsbrt > 1:
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtmodlmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtmodlmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
           
            if liststrgmodl[a] == 'true':
                listyerr = np.zeros((2, numbplot, gdat.numbener))
            
            # plot energy spectra of the data, background model components and total background
            if gdat.numbener > 1:
                
                listmrkr = ['o', '>', 's', 'h', '*', 'p', 'x']
                for k in range(100):
                    listmrkr.append('x')

                # determine the energy scaling factor
                if specconvunit[0] == 'en00':
                    factener = 1.
                if specconvunit[0] == 'en01':
                    factener = gdat.meanpara.ener
                if specconvunit[0] == 'en02':
                    factener = gdat.meanpara.ener**2
                if specconvunit[0] == 'en03':
                    # temp
                    pass
                    factener = 1.
                    #indxenerintv = np.where((gdat.meanpara.ener < specconvunit[4]) & (gdat.meanpara.ener > specconvunit[3]))[0]
                    #ener = np.concatenate((np.array([specconvunit[3]]), gdat.meanpara.ener[indxenerintv], np.array([specconvunit[4]])))
                    #
                    #for k in range(3):
                    #    if k == 0:
                    #        ydattemp = 
                    #    ydatminmener = np.interp(specconvunit[3], gdat.meanpara.ener, ydat)
                    #    ydatmaxmener = np.interp(specconvunit[4], gdat.meanpara.ener, ydat)
                    #    ydat = np.concatenate((np.array([ydatminmener]), ydat[indxenerintv], np.array([ydatmaxmener])))
                    #    ydat = np.trapz(ydat, gdat.meanpara.ener)
                    #
                    #yerrminmener = np.interp(specconvunit[3], gdat.meanpara.ener, yerr, axis=1)
                    #yerrmaxmener = np.interp(specconvunit[4], gdat.meanpara.ener, yerr, axis=1)
                    #ydat = np.stack((np.array([yerrminmener]), ydat[indxenerintv], np.array([yerrmaxmener])))
                    #
                    #
                    #yerr = np.trapz(yerr, gdat.meanpara.ener)


                xdat = gdat.meanpara.ener
                cntr = 0
                
                for k in range(listydat.shape[0]):
                    mrkr = listmrkr[cntr]
                    if k == cntrdata:
                        colr = 'black'
                        alph = 1.
                        linestyl = '-'
                    else:
                        colr = retr_colr(gdat, strgstat, liststrgmodl[a], indxpopl=None)
                        linestyl = '--'
                        alph = 0.5
                   
                    ydat = np.copy(listydat[k, :])
                    yerr = np.copy(listyerr[:, k, :])
                    
                    ydat *= factener
                    yerr *= factener
                    
                    if k == cntrdata and a > 0:
                        continue
                    
                    if liststrgmodl[a] == 'fitt':
                        labl = listlablsbrtspec[k]
                    else:
                        labl = None
                    
                    temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, color=colr, marker=mrkr, ls=linestyl, markersize=10, alpha=alph, label=labl)
                    for caps in listcaps:
                        caps.set_markeredgewidth(1)

                    if gdat.numbpixl == 1 and strgstat != 'pdfn':
                        if cntr != cntrline or k in indxplotelemendd:
                            cntr += 1
                    else:
                        cntr += 1

        if gdat.numbener > 1:
            axis.set_xlim([np.amin(gdat.binspara.ener), np.amax(gdat.binspara.ener)])
            
            if gdat.typeexpr == 'chan':
                factminm = 1e-1
                factmaxm = 1e2
            elif gdat.typeexpr == 'ferm':
                factminm = 1e1
                factmaxm = 1e-1
            else:
                factminm = 1e-4
                factmaxm = 1e0
            minmydat = factminm * gdat.factylimtbrt[0] * np.amax(listydat[cntrdata, :] * factener)
            maxmydat = factmaxm * gdat.factylimtbrt[1] * np.amax(listydat[cntrdata, :] * factener)
            limtydat = [minmydat, maxmydat]
            axis.set_ylim(limtydat)
            axis.set_yscale('log')
            axis.set_xlabel(gdat.lablenertotl)
            axis.set_xscale('log')
            labl = getattr(gmod.lablpara, 'sbrt' + specconvunit[0] + specconvunit[1] + 'stertotl')
            axis.set_ylabel(labl)
            make_legd(axis, numbcols=2)
            
            plt.tight_layout()
            path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'sdenmean%s%s%s' % (namespatmean, specconvunit[0], specconvunit[1]))
            figr.savefig(path)
            plt.close(figr)
        

def retr_factener(strgconvunit, ener):
    
    if strgconvunit == 'en00':
        factener = np.ones_like(ener)
    
    if strgconvunit == 'en01':
        factener = ener
    
    if strgconvunit == 'en02':
        factener = ener**2
    
    if strgconvunit == 'en03':
        # temp
        pass
        factener = np.ones_like(ener)
    
    return factener


def plot_pdfntotlflux():

    minm = 1e-9
    maxm = 10e-9
    numbvarb = 90
    numbparagenrfull = 100000
    numbbins = 40
    alph = 0.5
    
    binssing = np.linspace(minm, maxm, numbvarb + 1)
    meansing = (binssing[:-1] + binssing[1:]) / 2.
    deltsing = binssing[1:] - binssing[:-1]
    
    binsdoub = np.linspace(2. * minm, 2. * maxm, 2 * numbvarb)
    meandoub = (binsdoub[:-1] + binsdoub[1:]) / 2.
    deltdoub = binsdoub[1:] - binsdoub[:-1]
    
    bins = np.linspace(minm, 2. * maxm, 2 * numbvarb + 1)
    
    arry = np.empty((2, numbparagenrfull))
    
    minmslop = 1.5
    maxmslop = 3.
    numbslop = 4
    sloparry = np.linspace(minmslop, maxmslop, numbslop)
    for n in range(numbslop):
        slop = sloparry[n]
        for k in range(2):
            arry[k, :] = (np.random.rand(numbparagenrfull) * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
        
        totl = np.sum(arry, 0)
        
        powrprob = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop)) * meansing**(-slop)
        
        convprob = convolve(powrprob, powrprob) * deltdoub[0]
        
        indxdoub = np.where(meandoub <= maxm)[0]
        convprobpoly = polyval(polyfit(meandoub[indxdoub], convprob[indxdoub], 8), meandoub[indxdoub])
        
        figr, axis = plt.subplots()
        axis.hist(arry[k, :], bins=bins, alpha=alph, label='$f_1$ (Sampled)', color='b')
        axis.hist(totl, bins=bins, alpha=alph, label='$f_0$ (Sampled)', color='g')
        axis.plot(meansing, powrprob * numbparagenrfull * deltsing, label='$f_1$ (Analytic)', color='b')
        axis.plot(meandoub, convprob * numbparagenrfull * deltdoub[0], label='$f_0$ (Numerically convolved)', color='g')
        
        axis.plot(meandoub[indxdoub], convprobpoly * numbparagenrfull * deltdoub[indxdoub], label='$f_0$ (Fit)', color='r')
    
        axis.set_ylim([0.5, numbsamp])
        axis.set_xlabel('$f$')
        axis.set_xlim([np.amin(bins), np.amax(bins)])
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_ylabel('$N_{samp}$')
        make_legd(axis)
        plt.tight_layout()
        pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/powrpdfn/'
        figr.savefig(pathfold + 'powrpdfn%04d.pdf' % n)
        plt.close(figr)
        

def savefigr(gdat, gdatmodi, figr, path):
    
    #if gdatmodi is not None and gdat.numbproc > 1:
    #    gdatmodi.lock.acquire()
    #    print 'Process %d acquiring the lock...' % gdatmodi.indxprocwork 
    
    plt.savefig(path)
    
    #if gdatmodi is not None and gdat.numbproc > 1:
    #    gdatmodi.lock.release()
    #    print 'Process %d releasing the lock...' % gdatmodi.indxprocwork 
        

def plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, indxpoplfrst, strgfrst, \
                                                                                            strgseco, strgtotl, strgmome='pmea', strgpdfn='post'):
    
    gmod = getattr(gdat, strgmodl)
    
    sizelarg = 10
    sizesmll = 1
    
    if strgstat == 'pdfn':
        lablmome = getattr(gdat, 'labl' + strgmome)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if strgmodl == 'fitt':
        colrtemp = gmod.colrelem[indxpoplfrst]
        if strgstat == 'pdfn':
            labl = gdat.lablsampdist + ' ' + lablmome
            if strgelemtdimtype == 'bind':
                varb = getattr(gdat, strgmome + strgpdfn + strgtotl)
                varbfrst = gdat.binspara.strgfrst
                varbseco = getattr(gdat.binspara, strgseco)
                if strgtotl.startswith('hist') or strgtotl.startswith('exr') or strgtotl.startswith('incr') or np.amax(varb) <= 0.:
                    normtdim = None
                else:
                    normtdim = mpl.colors.LogNorm(0.5, vmax=np.amax(varb))
                imag = axis.pcolor(varbfrst, varbseco, varb.T, cmap='Blues', label=labl, norm=normtdim)
                make_cbar(gdat, axis, imag)
                
            else:
                if gdat.boolcondcatl:
                    varbfrst = np.zeros(gdat.numbprvlhigh)
                    varbseco = np.zeros(gdat.numbprvlhigh)
                    cntr = 0
                    for r in gdat.indxstkscond:
                        if r in gdat.indxprvlhigh:
                            varbfrst[cntr] = gdat.dictglob['poststkscond'][r][strgfrst][indxpoplfrst]
                            varbseco[cntr] = gdat.dictglob['poststkscond'][r][strgseco][indxpoplfrst]
                            cntr += 1
                    axis.scatter(varbfrst, varbseco, alpha=gdat.alphelem, color=colrtemp, label=gdat.lablparagenrscalfull)
        
        if strgstat == 'this' or strgstat == 'mlik':
            if strgelemtdimtype == 'bind':
                meanfrst = getattr(gdat.binspara, strgfrst)
                meanseco = getattr(gdat.binspara, strgseco)
                hist = getattr(gdatmodi, strgstat + strgtotl)
                if strgtotl.startswith('hist') or strgtotl.startswith('exr') or strgtotl.startswith('incr') or np.amax(hist) <= 0.:
                    normtdim = None
                else:
                    normtdim = mpl.colors.LogNorm(0.5, vmax=np.amax(hist))
                imag = axis.pcolor(meanfrst, meanseco, hist.T, cmap='Blues', label=gdat.lablparagenrscalfull, alpha=gdat.alphhist, norm=normtdim)
            else:
                varbfrst = getattr(gdatmodi.this, strgfrst)[indxpoplfrst]
                varbseco = getattr(gdatmodi.this, strgseco)[indxpoplfrst]
                if len(varbfrst) == 0 or len(varbseco) == 0:
                    varbfrst = np.array([limtfrst[0] * 0.1])
                    varbseco = np.array([limtseco[0] * 0.1])
                axis.scatter(varbfrst, varbseco, alpha=gdat.alphelem, color=colrtemp, label=gdat.lablparagenrscalfull)
    
    # reference elements
    if strgfrst[-4:] in gdat.listnamerefr:
        strgfrsttemp = strgfrst[-4:]
    else:
        strgfrsttemp = strgfrst
    if strgseco[-4:] in gdat.listnamerefr:
        strgsecotemp = strgseco[-4:]
    else:
        strgsecotemp = strgseco
    if hasattr(gdat.refr, strgfrsttemp) and hasattr(gdat.refr, strgsecotemp):
        for q in gdat.indxrefr:
            if strgfrsttemp in gdat.refr.namepara.elem[q] and strgsecotemp in gdat.refr.namepara.elem[q]:
                refrvarbfrst = getattr(gdat.refr, strgfrsttemp)[q]
                refrvarbseco = getattr(gdat.refr, strgsecotemp)[q]
                if len(refrvarbfrst) == 0 or len(refrvarbseco) == 0:
                    refrvarbfrst = np.array([limtfrst[0] * 0.1])
                    refrvarbseco = np.array([limtseco[0] * 0.1])
                axis.scatter(refrvarbfrst, refrvarbseco, alpha=gdat.alphelem, color=gdat.refr.colrelem[q], label=gdat.refr.lablelem[q], s=sizelarg)

    plot_sigmcont(gdat, strgmodl, axis, strgfrst, indxpoplfrst, strgseco=strgseco)
    
    scalfrst = getattr(gmod.scalpara, strgfrst)
    scalseco = getattr(gmod.scalpara, strgseco)

    if scalfrst == 'logt':
        axis.set_xscale('log')
    if scalseco == 'logt':
        axis.set_yscale('log')
    
    axis.set_xlabel(getattr(gmod.labltotlpara, strgfrst))
    axis.set_ylabel(getattr(gmod.labltotlpara, strgseco))
    axis.set_xlim(getattr(gmod.limtpara, strgfrst))
    axis.set_ylim(getattr(gmod.limtpara, strgseco))
    
    make_legd(axis)

    plt.tight_layout()
    if strgstat == 'pdfn':
        strgmometemp = strgmome
    else:
        strgmometemp = ''
    
    nameinte = strgelemtdimvarb + 'tdim/'
    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, '%s%s' % (strgmometemp, strgtotl), nameinte=nameinte)
    
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_sigmcont(gdat, strgmodl, axis, strgfrst, indxpoplfrst, strgseco=None):
    
    if strgfrst == 'deltllik' or strgseco == 'deltllik':
        for pval in gdat.pvalcont:
            if strgfrst == 'deltllik':
                deltlliksigm = scipy.stats.chi2.ppf(1. - pval, gmod.numbparagenrelemsing[indxpoplfrst])
                axis.axvline(deltlliksigm, ls='--', color='black', alpha=0.2) 
            if strgseco == 'deltllik':
                deltlliksigm = scipy.stats.chi2.ppf(1. - pval, gmod.numbparagenrelemsing[indxpoplfrst])
                axis.axhline(deltlliksigm, ls='--', color='black', alpha=0.2) 
    

def plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgydat, strgxdat, typehist='hist', \
                     indxrefrplot=None, indxydat=None, strgindxydat=None, indxxdat=None, strgindxxdat=None, plottype='none', \
                     meanxdat=None, \
                     scal=None, scalxdat=None, scalydat=None, limtxdat=None, limtydat=None, omittrue=False, nameinte='', \
                     lablxdat='', lablydat='', histodim=False, offslegd=None, booltdim=False, ydattype='totl', boolhistprio=True):
   
    gmod = getattr(gdat, strgmodl)
    gmodstat = getattr(gmod, strgstat)

    if strgydat[-8:-5] == 'pop':
        boolelem = True
    else:
        boolelem = False

    if scal is None:
        if scalxdat is None:
            scalxdat = 'linr'
        if scalydat is None:
            scalydat = 'linr'
    else:
        scalxdat = scal
        scalydat = scal

    if histodim:
        figrsize = (gdat.plotsize, 0.8 * gdat.plotsize)
    else:
        figrsize = (gdat.plotsize, gdat.plotsize)

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    
    if booltdim:
        xdat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgxdat, strgpdfn)
        ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn)
    else:
        xdat = getattr(gdat.meanpara, strgxdat[4:])
        if typehist == 'histcorrreca':
            ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'histcorrreca' + strgydat[4:], strgpdfn)
        else:
            ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn)
    
    if indxxdat is not None:
        xdat = xdat[indxxdat]
    if indxydat is not None:
        ydat = ydat[indxydat]
    
    xerr = np.zeros((2, xdat.size))
    
    if booltdim:
        axis.scatter(xdat, ydat, alpha=gdat.alphelem, color=colr, label=gdat.lablparagenrscalfull)
    else:
        if histodim:
            # temp
            if strgxdat[4:] in gmod.namepara.elem:
                deltxdat = getattr(gdat.deltpara, strgxdat[4:])
                binsxdat = getattr(gdat.binspara, strgxdat[4:])
            else:
                deltxdat = getattr(gdat.deltpara, strgxdat[4:])
                binsxdat = getattr(gdat.binspara, strgxdat[4:])

            xdattemp = binsxdat[:-1] + deltxdat / 2.
   
    if strgmodl == 'fitt':
        if boolelem:
            if strgydat.startswith('cmpl'):
                labl = gmod.lablelem[int(strgydat[-5])]
                colr = gmod.colrelem[int(strgydat[-5])]
            else:
                labl = gmod.lablelem[int(strgydat[-1])]
                colr = gmod.colrelem[int(strgydat[-1])]
        else:
            labl = gmod.labl
            colr = gmod.colr
        
        if strgstat == 'pdfn':
            if typehist == 'histcorrreca':
                yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'histcorrreca' + strgydat[4:], strgpdfn, strgmome='errr')
            else:
                yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn, strgmome='errr')
            if indxydat is not None:
                yerr = yerr[[slice(None)] + indxydat]
            
            # label
            if strgydat.startswith('hist'):
                ##  element distribution
                labl = gdat.lablsampdist
            else:
                ##  other
                labl = gdat.lablsampdist
            
            # draw points
            indxerrr = np.where((yerr[0, :] > 0.) | (yerr[1, :] > 0.))[0]
            if indxerrr.size > 0:
                labltemp = None
            else:
                labltemp = labl
            temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, label=labl, \
                                                                                       marker='o', ls='', markersize=5, color=colr, lw=1, capsize=5)

            # draw error-bar caps 
            if indxerrr.size > 0:
                temp, listcaps, temp = axis.errorbar(xdat[indxerrr], ydat[indxerrr], yerr=yerr[:, indxerrr], xerr=xerr[:, indxerrr], \
                                                                                      marker='o', ls='', markersize=5, color=colr, lw=1, capsize=5)
                for caps in listcaps:
                    caps.set_markeredgewidth(1)

        elif strgstat == 'this' or strgstat == 'mlik':
            
            if strgstat == 'this':
                labl = gdat.lablsamp
            else:
                labl = gdat.lablmlik

            if histodim:
                axis.bar(xdattemp, ydat, deltxdat, label=gdat.lablparagenrscalfull, alpha=0.5, linewidth=1, edgecolor=colr)
            else:
                if plottype == 'errr':
                    yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn, strgmome='errr')

                    if indxydat is not None:
                        yerr = yerr[[slice(None)] + indxydat]
                    temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, \
                                            marker='o', ls='', markersize=5, label=labl, lw=1, capsize=5, color=colr)
                    for caps in listcaps:
                        caps.set_markeredgewidth(1)
                else:
                    axis.plot(xdat, ydat, label=gdat.lablparagenrscalfull, alpha=0.5, color=colr)
    
    # reference histogram
    if not omittrue:
        for q in gdat.indxrefr:
            
            if boolelem:
                if strgydat[-12:-8] in gdat.listnamerefr:
                    name = 'refr' + strgydat[:-12] + 'pop%d' % q + strgydat[-4:]
                else:
                    name = 'refr' + strgydat[:-8] + 'pop%d' % q + strgydat[-4:]
            else:
                name = 'refr' + strgydat
            
            if not hasattr(gdat, name):
                continue
            
            ydattemp = getattr(gdat, name)
            
            ydat = ydattemp
            if indxydat is not None:
                ydat = ydat[indxydat]
            
            if strgydat[-8:-5] == 'pop':
                labl = gdat.refr.lablelem[q]
                colr = gdat.refr.colrelem[q]
            else:
                labl = gdat.refr.labl
                colr = gdat.refr.colr
    
            if histodim:
                axis.bar(xdattemp, ydat, deltxdat, color=colr, label=labl, alpha=gdat.alphhist, linewidth=1, edgecolor=colr)
            else:
                axis.plot(xdat, ydat, color=colr, label=labl, alpha=gdat.alphline)
                           
            try:
                if histodim:
                    if typehist == 'histcorrreca':
                        reca = getattr(gdat.true, 'reca' + strgydat[4:])
                    axis.plot(xdattemp, 10. * reca, color='purple', label='PTFN', alpha=gdat.alphline)
            except:
                pass

            if not boolelem:
                break
    
    # external reference histogram
    if histodim and strgydat == 'histfluxpop0':
        try:
            if gdat.listprefhistfluxlabl is not None:
                for k in range(len(gdat.listprefhistfluxlabl)):
                    if gdat.listprefhistfluxtype[k] == 'shad':
                        axis.plot(gdat.listprefhistfluxflux[k][0], gdat.listprefhistfluxhist[k][0], color='m', label=gdat.listprefhistfluxlabl[k])
                        enerpoly = np.empty(gdat.listprefhistfluxflux[k][1].size + gdat.listprefhistfluxflux[k][2].size)
                        enerpoly[:gdat.listprefhistfluxflux[k][1].size] = gdat.listprefhistfluxflux[k][1]
                        enerpoly[gdat.listprefhistfluxflux[k][1].size:] = gdat.listprefhistfluxflux[k][2][::-1]
                        sbrtpoly = np.empty(gdat.listprefhistfluxflux[k][1].size + gdat.listprefhistfluxflux[k][2].size)
                        sbrtpoly[:gdat.listprefhistfluxflux[k][1].size] = gdat.listprefhistfluxhist[k][1]
                        sbrtpoly[gdat.listprefhistfluxflux[k][1].size:] = gdat.listprefhistfluxhist[k][2][::-1]
                        axis.fill(enerpoly, sbrtpoly, color='m', alpha=0.5)
                    else:
                        axis.errorbar(gdat.listprefhistfluxflux[k], gdat.listprefhistfluxhist[k], label=gdat.listprefhistfluxlabl[k], color='m')
        except:
            pass

    if strgydat.startswith('histcntp'):
        ydattemp = getattr(gmodstat, strgydat)
        axis.bar(xdattemp, ydattemp, deltxdat, color='black', label='Data', alpha=gdat.alphhist, linewidth=1, edgecolor='black')
                
    # axis scales
    if scalxdat == 'logt':
        axis.set_xscale('log')
    if scalydat == 'logt':
        if np.where(ydat > 0.)[0].size > 0:
            axis.set_yscale('log')
    
    # axis labels
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)

    # superimpose prior on the feature
    ptch = None
    line = None

    if strgydat.startswith('hist') and strgydat != 'histdefl' and strgydat != 'histdeflelem' and boolhistprio:
        if strgydat[-8:-5] == 'pop':
            strgtemp = strgydat[4:-8]
            if strgtemp in gmod.namepara.genrelem[int(strgydat[-5])]:
                xdatprio = getattr(gmod, strgxdat + 'prio')
                if gdat.typedata == 'mock' and not omittrue:
                    for q in gdat.indxrefr:
                        if gdat.refr.numbelem[q] == 0:
                            continue
                        if strgtemp in gmod.namepara.genrelem[q]:
                            truexdatprio = getattr(gdat.true, strgxdat + 'prio')
                            trueydatsupr = getattr(gdat.true, strgydat + 'prio')
                            trueydatsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'true', strgydat + 'prio', strgpdfn)
                            axis.plot(truexdatprio, trueydatsupr, ls='-', alpha=gdat.alphline, color=gdat.refr.colrelem[q])

                if strgmodl != 'true':
                    ydatsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'fitt', strgydat + 'prio', strgpdfn)
                    if strgstat == 'pdfn':
                        yerrsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'fitt', strgydat + 'prio', strgpdfn, strgmome='errr')
                        labl = gdat.lablsampdist + ' hyper-distribution'
                        ptch, line = tdpy.plot_braz(axis, xdatprio, ydatsupr, yerr=yerrsupr, lcol='lightgrey', dcol='grey', labltotl=labltotl)
                    else:
                        axis.plot(xdatprio, ydatsupr, ls='--', alpha=gdat.alphline, color=gmod.colrelem[int(strgydat[-5])])
   
    for name, valu in gdat.refr.__dict__.items():
        if name[8:12] == 'hist' and name[12:16] == strgydat[4:] and name[16:19] == 'pop' and int(name[-1]) == indxpopltemp:
            colr = getattr(gdat, name + 'colr')
            linestyl = getattr(gdat, name + 'linestyl')
            axis.plot(valu[0, :], valu[1, :], ls=linestyl, color=colr)

    if strgydat.startswith('hist') and strgydat[4:-8] == 'deltllik':
        plot_sigmcont(gdat, strgmodl, axis, strgxdat[4:], int(strgydat[-1]))
   
    if indxydat is not None:
        strgydat += strgindxydat
    
    if indxxdat is not None:
        strgxdat += strgindxxdat
    
    if limtxdat is not None:
        axis.set_xlim(limtxdat)
    else:
        axis.set_xlim([np.amin(xdat), np.amax(xdat)])
    if limtydat is not None:
        axis.set_ylim([limtydat[0], limtydat[1]])
    else:
        axis.set_ylim([np.amin(ydat), np.amax(ydat)])
    
    if ydattype != 'totl':
        strgydat += ydattype
    
    try:
        make_legd(axis, offs=offslegd, ptch=ptch, line=line)
    except:
        print('Legend failed when')
        print('strgstat')
        print(strgstat)
        print('strgmodl')
        print(strgmodl)
        print('strgydat')
        print(strgydat)
        raise Exception('')

    plt.tight_layout()
    if typehist == 'histcorrreca':
        path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'histcorrreca' + strgydat[4:], nameinte=nameinte)
    else:
        path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, strgydat, nameinte=nameinte)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat, plotdiff=False):
    
    if plotdiff:
        figrsize = (gdat.plotsize, 0.7 * gdat.plotsize)
    else:
        figrsize = (gdat.plotsize, gdat.plotsize)
    figr, axis = plt.subplots(1, 1, figsize=figrsize)
    
    # prepare data to be plotted
    xdat = np.copy(getattr(gdat.refr, strgfeat)[q][0, :])
    xerr = tdpy.retr_errrvarb(getattr(gdat.refr, strgfeat)[q])
   
    ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgfeat + 'asscpop%dpop%d' % (q, l), strgpdfn)
    
    yerr = np.zeros((2, ydat.size))
    if strgstat == 'pdfn':
        yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgfeat + 'asscpop%dpop%d' % (q, l), strgpdfn, strgmome='errr')
    
    if plotdiff:
        ydat = 100. * (ydat - xdat) / xdat
    
    # handle the case when there is a single reference element
    if yerr.ndim == 1:
        ydat = np.array([ydat])
        yerr = yerr[:, None]
    
    # plot all associations
    if plotdiff:
        indx = np.where(ydat > -100.)[0]
    else:
        indx = np.where(ydat > 0.)[0]
    if indx.size > 0:
        axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='black')
    
    # temp -- plot associations inside the comparison area
    if plotdiff:
        axis.axhline(0., ls='--', alpha=gdat.alphline, color='black')
    else:
        axis.plot(binsplot, binsplot, ls='--', alpha=gdat.alphline, color='black')
    
    lablxdat = getattr(gmod.lablpara, strgfeat + 'refr')
    lablydat = getattr(gmod.lablpara, strgfeat + 'paragenrscalfull')
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)
    boollogtxaxi = False
    boollogtyaxi = False
    if indx.size > 0 and scal == 'logt':
        if not plotdiff:
            axis.set_yscale('log')
            boollogtyaxi = True
        axis.set_xscale('log')
        boollogtaxis = True
   
    if plotdiff:
        limtydat = np.array([-100., 100.])
    else:
        limtydat = np.array([minmplot, maxmplot])
    limtxdat = [minmplot, maxmplot]
    
    # overplot text
    if 'etag' in gdat.refr.namepara.elem[q]:
        for k in range(indx.size):
            if boollogtxaxi:
                sizexoff = 0.01 * xdat[indx[k]]
            else:
                sizexoff = 0.01 * (limtxdat[1] - limtxdat[0])
            if boollogtyaxi:
                sizeyoff = 0.01 * ydat[indx[k]]
            else:
                sizeyoff = 0.01 * (limtydat[1] - limtydat[0])
            axis.text(xdat[indx[k]] + sizexoff, ydat[indx[k]] + sizeyoff, gdat.refretag[q][indx[k]], verticalalignment='center', horizontalalignment='center', \
                                                                                                                                                        color='red', fontsize=1)

    axis.set_ylim(limtydat)
    axis.set_xlim(limtxdat)
   
    plt.tight_layout()
    if plotdiff:
        strgtype = 'diff'
    else:
        strgtype = ''
    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'scatassc' + strgfeat + '%spop%dpop%d' % (strgtype, q, l), nameinte='assc')
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxevttplot, indxenerplot=None):
    
    gmod = getattr(gdat, strgmodl)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn)
    if indxenerplot is None:
        xdat = gdat.cntpdata[:, :, indxevttplot].flatten()
        ydat = ydat[:, :, indxevttplot].flatten()
        nameplot = 'scatcntpevt%d' % (indxevttplot)
        if strgstat == 'pdfn':
            indxvarb = [slice(None), slice(None), indxevttplot]
    else:
        xdat = gdat.cntpdata[indxenerplot, :, indxevttplot]
        ydat = ydat[indxenerplot, :, indxevttplot]
        nameplot = 'scatcntpen%02devt%d' % (indxenerplot, indxevttplot)
        if strgstat == 'pdfn':
            indxvarb = [indxenerplot, slice(None), indxevttplot]
    if strgstat == 'pdfn':
        yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn, strgmome='errr', indxvarb=indxvarb)
    colr = gmod.colr

    if strgstat == 'pdfn':
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color=gmod.colr, capsize=5)
    else:
        axis.plot(xdat, ydat, marker='o', ls='', markersize=5, color=gmod.colr)
    gdat.limtcntpdata = [gdat.binspara.cntpdata[0], gdat.binspara.cntpdata[-1]]
    axis.set_xlim(gdat.limtcntpdata)
    axis.set_ylim(gdat.limtcntpdata)
    axis.set_ylabel('$k^{modl}$')
    axis.set_xlabel('$k^{data}$')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.tight_layout()

    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, nameplot)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 40
    numbfluxprox = len(gdat.indxpixlprox)
    bins = np.empty((gdat.numbprox, numbbins + 1))
    indxpixlproxsize = np.empty((numbfluxprox, gdat.numbpixlfull))
    for h in gdat.indxprox:
        for j in gdat.indxpixlfull:
            try:
                indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
            except:
                indxpixlproxsize[h, j] = gdat.numbpixlfull
        bins[h, :] = np.logspace(np.log10(np.amin(indxpixlproxsize[h, :])), np.log10(np.amax(indxpixlproxsize[h, :])), numbbins + 1)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in gdat.indxprox:
        axis.hist(indxpixlproxsize[h, :], bins=bins[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphhist)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixlfull, label='ROI', ls='--')
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    make_legd(axis)
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'init/indxprox.pdf')
    plt.close()
    
    
def plot_psfn_type():
    
    devi = np.linspace(0., 5., 100)
    y = np.zeros((x.size, 5))

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    singgaus = retr_singgaus(devi, 0.25)
    axis.plot(devi, singgaus, label='Single Gaussian')

    singking = retr_singking(devi, 0.25, 10.)
    axis.plot(devi, singking, label='Single King')

    doubgaus = retr_doubgaus(devi, 0.1, 0.25, 1.)
    axis.plot(devi, doubgaus, label='Double Gaussian')

    gausking = retr_gausking(devi, 0.1, 0.25, 1., 10.)
    axis.plot(devi, gausking, label='Gaussian + King')

    doubking = retr_doubking(devi, 0.1, 0.25, 10., 1., 5.)
    axis.plot(devi, doubking, label='Double King')

    make_legd(axis)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylim([1e-3, None])
    
    
def plot_evidtest():
    
    minmgain = -1.
    maxmgain = 5.
    minmdevi = 0.
    maxmdevi = 5.
    gain = np.linspace(minmgain, maxmgain, 100)
    devi = np.linspace(minmdevi, maxmdevi, 100)

    evid = np.log(np.sqrt(1. + np.exp(2. * gain[None, :])) * np.exp(-devi[:, None]**2 / 2. / (1. + 1. / np.exp(2. * gain[None, :]))))
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    figr.suptitle('Log-Bayesian Evidence For Lower-Dimension Model', fontsize=18)
    imag = axis.imshow(evid, extent=[minmgain, maxmgain, minmdevi, maxmdevi], cmap='winter', origin='lower')
    cset1 = plt.contourf(gain, devi, evid, cmap='winter')
    axis.set_xlabel('Information gain')
    axis.set_ylabel('Goodness of fit')
    plt.colorbar(imag, ax=axis, fraction=0.03)

    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'evidtest.pdf')
    plt.close(figr)
    
    
def plot_histlgalbgalelemstkd(gdat, strgpdfn, indxpoplplot, strgbins, strgfeat=None):
    
    if strgfeat is not None:
        numbparaplot = gdat.numbbinsplot
    else:
        numbparaplot = 1

    if strgbins == 'cumu':
        numbrows = 1
        numbcols = 1
    else:
        numbcols = 2
        if strgbins == 'full':
            numbrows = numbparaplot / 2
        else:
            numbrows = 2
    
    histlgalbgalelemstkd = getattr(gdat, strgpdfn + 'histlgalbgalelemstkd')

    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
    if numbrows == 1:
        axgr = [axgr]            
    for a, axrw in enumerate(axgr):
        if numbcols == 1:
            axrw = [axrw]
        for b, axis in enumerate(axrw):
            if strgfeat is not None:
                h = a * 2 + b
                if strgbins == 'full':
                    indxlowr = h
                    indxuppr = h + 1
                elif strgbins == 'cumu':
                    indxlowr = 0
                    indxuppr = numbparaplot
                else:
                    if h < 3:
                        indxlowr = 2 * h
                        indxuppr = 2 * (h + 1)
                    else:
                        indxlowr = 2 * h
                        indxuppr = numbparaplot
                temp = np.sum(histlgalbgalelemstkd[indxpoplplot][:, :, indxlowr:indxuppr], 2).T
            else:
                temp = np.sum(np.sum(histlgalbgalelemstkd[indxpoplplot], 2), 2).T
                
            if np.where(temp > 0.)[0].size > 0:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', \
                                                            extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
            else:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                
            if strgfeat is not None:
                bins = getattr(gdat.binspara, strgfeat)
            
            # superimpose reference elements
            for q in gdat.indxrefr:
                if gdat.refr.numbelem[q] == 0:
                    continue
                # temp -- backcomp
                reframpl = getattr(gdat.refr, gdat.refr.nameparagenrelemampl[q])
                if strgfeat in gdat.refr.namepara.elem[q]:
                    refrfeat = getattr(gdat.refr, strgfeat)[q]
                    if len(refrfeat) > 0:
                        indxelem = np.where((bins[indxlowr] < refrfeat[0, :]) & (refrfeat[0, :] < bins[indxuppr]))[0]
                    else:
                        indxelem = np.array([])
                else:
                    indxelem = np.arange(gdat.refr.numbelem[q])
                # temp -- backcomp
                try:
                    mrkrsize = retr_mrkrsize(gdat, strgmodl, reframpl[q][0, indxelem], gdat.refr.nameparagenrelemampl[q])
                except:
                    mrkrsize = retr_mrkrsize(gdat, strgmodl, reframpl[q][0, indxelem], gdat.refr.nameparagenrelemampl[q])

                if indxelem.size > 0:
                    axis.scatter(gdat.anglfact * gdat.refr.dictelem[q]['lgal'][0, indxelem], gdat.anglfact * gdat.refr.dictelem[q]['bgal'][0, indxelem], \
                                                s=mrkrsize, alpha=gdat.alphelem, marker=gdat.refrlistmrkrhits[q], lw=2, color=gdat.refr.colrelem[q])

            if a == numbrows - 1:
                axis.set_xlabel(gdat.labllgaltotl)
            else:
                axis.set_xticklabels([])
            if b == 0:
                axis.set_ylabel(gdat.lablbgaltotl)
            else:
                axis.set_yticklabels([])

            draw_frambndr(gdat, axis)
            
            if strgbins != 'cumu':
                titl = tdpy.mexp(bins[indxlowr]) + ' < $%s$ < ' % lablfeat + tdpy.mexp(bins[indxuppr])
                axis.set_title(titl)
    
    if strgfeat is not None:
        lablfeattotl = getattr(gmod.lablpara, strgfeat + 'totl')
        plt.figtext(0.5, 0.95, '%s' % lablfeattotl, ha='center', va='center')
    axiscomm = figr.add_axes([0.87, 0.2, 0.02, 0.6])
    cbar = figr.colorbar(imag, cax=axiscomm)

    plt.subplots_adjust()
    #plt.subplots_adjust(left=0.18, top=.9, right=0.82, bottom=0.15, hspace=0.08, wspace=0.08)
    if strgbins == 'cumu':
        strgtemp = ''
    else:
        strgtemp = strgfeat
    path = getattr(gdat, 'path' + strgpdfn + 'finl') + 'histlgalbgalelemstkd%s%spop%d' % (strgbins, strgtemp, indxpoplplot) + '.pdf'
    figr.savefig(path)
    plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.binspara.angl)

    figr, axgr = plt.subplots(1, 2, figsize=(2 * gdat.plotsize, gdat.plotsize))
    figr.suptitle('King Function', fontsize=20)
    for k, axis in enumerate(axgr):
        if k == 0:
            sigmlist = [0.25]
            gammlist = [1.01, 2.5, 10.]
        else:
            sigmlist = [0.1, 0.25, 1.]
            gammlist = [2.]
        for sigm in sigmlist:
            for gamm in gammlist:
                axis.plot(angl, retr_singking(angl, sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        make_legd(axis)
        axis.set_yscale('log')
        axis.set_xlabel(gdat.labltotlpara.gang)
        axis.set_xlabel(r'$\mathcal{K}$')
        
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'king.pdf')
    plt.close(figr)
    
   
def plot_intr(gdat):
    
    if gdat.typeverb > 0:
        print('Making PCAT introductory plots...')

    #plot_grap(plottype='meta', typeverb=1)
    plot_grap(plottype='lght0000', typeverb=1)
    #plot_grap(plottype='lght0001', typeverb=1)
    #plot_grap(plottype='lght0002', typeverb=1)
    #plot_grap(plottype='lght0003', typeverb=1)
    #plot_grap(plottype='lens0000', typeverb=1)
    plot_grap(plottype='lens0001', typeverb=1)
    
    with plt.xkcd():

        from matplotlib import patheffects
        mpl.rcParams['path.effects'] = [patheffects.withStroke(linewidth=0)]

        figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))

        catl = np.arange(80)
        probcatl = pss.pmf(catl, 30.) + 0.5 * pss.pmf(catl, 60.)
        axis.plot(catl, probcatl)
        axis.set_xticks([10, 30, 60])
        axis.set_xticklabels(["Crackpot's Catalog", "Best-fit catalog", "Not-so-best-fit catalog"])
        axis.set_yticks([])
        text = axis.set_title("Exploring the catalog space with Probabilistic cataloging")
        text.set_position([.5, 1.05])
        axis.set_xlabel('Catalog index')
        axis.set_ylabel("Probability")
        
        axis.tick_params(axis='x', colors='#B6E954')
        axis.tick_params(axis='y', colors='#B6E954')
        axis.spines['bottom'].set_color('#B6E954')
        axis.spines['top'].set_color('#B6E954') 
        axis.spines['right'].set_color('#B6E954')
        axis.spines['left'].set_color('#B6E954')
        axis.yaxis.label.set_color('#B6E954')
        axis.xaxis.label.set_color('#B6E954')
        axis.title.set_color('#B6E954')

        axis.set_axis_bgcolor('black')
        figr.set_facecolor('black')
        plt.tight_layout()
        figr.savefig(gdat.pathimag + 'talkintr.pdf', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_psfn(gdat, gdatmodi, strgstat, strgmodl):

    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            for k in range(gdat.numbprox + 1):
                if k == 0 or k == gdat.numbprox:
                    alph = 1.
                    colr = 'b'
                    if k == 0:
                        labl = 'Dimmest PS'
                    else:
                        labl = 'Brightest PS'
                else:
                    alph = 0.2
                    labl = None
                    colr = 'black'
                axis.plot(gdat.binspara.angl * gdat.anglfact, gdat.binspara.prox[k] * gmodstat.psfn[i, :, m], label=labl, color=colr, alpha=alph)
                axis.set_xlim([np.amin(gdat.binspara.angl) * gdat.anglfact, np.amax(gdat.binspara.angl) * gdat.anglfact])
                if k > 0:
                    axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
            axis.set_yscale('log')
            axis.set_xlabel(gdat.labltotlpara.gang)
            axis.set_ylabel(gdat.lablsbrttotl)
            
            limt = gdat.specfraceval * np.amax(gdat.binspara.prox[0] * gmodstat.psfn[i, :, m])

            if limt != 0.:
                axis.axhline(limt, color='red', ls=':', label='Flux floor')
            
            make_legd(axis)

            plt.tight_layout()
            name = 'psfn'
            if gdat.numbener > 1:
                name += 'en%02d' % i
            if gdat.numbevtt > 1:
                name += 'evt%d' % m
            figr.savefig(gdat.pathinit + name + '.pdf')
            plt.close(figr)


def plot_mosa(gdat, strgpdfn):

    # empty global object
    gdatmodi = tdpy.gdatstrt()
    
    listparagenrscalfull = getattr(gdat, 'list' + strgpdfn + 'paragenrscalfull')
    listparagenrunitfull = getattr(gdat, 'list' + strgpdfn + 'paragenrunitfull')

    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamptotl:
        indxsampmosa = np.random.choice(gdat.indxsamptotl, size=numbsampmosa, replace=False)
        for l in gmod.indxpopl:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            
                            n = indxsampmosa[numbcols*a+b]
                            gdatmodi.this.paragenrscalfull = listparagenrscalfull[n, :].flatten()
                            gdatmodi.this.paragenrunitfull = listparagenrunitfull[n, :].flatten()
                            if gmod.numbparaelem > 0:
                                gdatmodi.this.indxelemfull = getattr(gdat, 'list' + strgpdfn + 'indxelemfull')[n]
                                proc_samp(gdat, gdatmodi, 'this', 'fitt')

                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.labllgaltotl)
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.lablbgaltotl)
                            else:
                                axis.set_yticklabels([])
                            
                            imag = retr_imag(gdat, axis, gdat.cntpdata, '', 'fitt', 'cntpdata', i, m)
                            supr_fram(gdat, gdatmodi, 'this', 'fitt', axis, l)
                    
                    if gdat.boolbinsener:
                        plt.figtext(0.5, 0.93, gdat.strgener[i], ha='center', va='center')
                    axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
                    cbar = figr.colorbar(imag, cax=axiscomm)
                    cbar.set_ticks(gdat.valutickmajrpara.cntpdata)
                    cbar.set_ticklabels(gdat.labltickmajrpara.cntpdata)
                    plt.subplots_adjust()
                    #plt.subplots_adjust(left=0.1, top=.91, hspace=0.03, wspace=0.1, bottom=0.09)
                    if l == 1:
                        strg = ''
                    else:
                        strg = 'pop%d' % l
                    pathfinl = getattr(gdat, 'path' + strgpdfn + 'finl')
                    if m is None:
                        path = pathfinl + 'mosa' + strg + 'en%02dA.pdf' % (gdat.indxenerincl[i])
                    else:
                        path = pathfinl + 'mosa' + strg + 'en%02devtt%d.pdf' % (gdat.indxenerincl[i], gdat.indxevttincl[m])
                    figr.savefig(path)
                    plt.close(figr)
    else:
        if gdat.typeverb > 0:
            print('Skipping the mosaic plot...')


def plot_grap(plottype, typeverb=0):
        
    import networkx as nx

    figr, axis = plt.subplots(figsize=(6, 6))

    grap = nx.DiGraph()
    if plottype == 'meta':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'magenta']


    if plottype == 'lens0001':
        listcolr = ['olive', 'olive', 'black', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'olive', 'olive', 'olive', 'olive', 'olive', \
                                                                                                                                        r'black', 'olive', 'black']

    if plottype == 'lght0000':
        listcolr = [r'olive', r'black', r'magenta', r'magenta', 'magenta', r'magenta', r'olive', r'olive', r'black', r'olive', r'olive', r'black', r'olive']
    



    if plottype == 'lght0001':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta', 'black']

    if plottype == 'lght0002':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'olive', 'magenta', \
                                                                                                    'magenta', 'magenta', 'magenta', 'magenta', 'black']
    if plottype == 'lght0003':
        listcolr = ['black', 'black', 'black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', \
                                                                                                    'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']
    
    if plottype == 'lens0000':
        listcolr = ['olive', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', 'magenta']


    if plottype.startswith('meta'):
        grap.add_edges_from([ \
                             ('meanelem', 'numbelem'), \
                             ('modl','data'), \
                             ('psfp', 'modl'), \
                             ('feat','modl'), \
                             ('numbelem','feat'), \
                             ('amplslop', 'ampl'), \
                            ])
    
    if plottype.startswith('lght') or plottype.startswith('lens'):
        grap.add_edges_from([ \
                             ('meanelem', 'numbelem'), \
                             ('modl','data'), \
                             ('psfp', 'modl'), \
                             ('bacp', 'modl'), \
                             ('lgal','modl'), \
                             ('bgal','modl'), \
                             ('numbelem','lgal'), \
                             ('numbelem','bgal'), \
                            ])
    
    if plottype.startswith('lght'):
        grap.add_edges_from([ \
                             ('amplslop', 'ampl'), \
                             ('ampl', 'modl'), \
                             ('numbelem','ampl'), \
                             ('numbelem', 'sind'), \
                             ('sind','modl'), \
                            ])
    
    if plottype.startswith('lens'):
        grap.add_edges_from([ \
                             ('lenp', 'modl'), \
                             ('defsslop', 'defs'), \
                             ('defs', 'modl'), \
                             ('numbelem','defs'), \
                            ])
    
    if plottype == 'lens0001':
        grap.add_edges_from([ \
                             ('asca', 'modl'), \
                             ('numbelem','asca'), \
                             ('acut', 'modl'), \
                             ('numbelem','acut'), \
                            ])
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('sinddistmean', 'sind'), \
                            ])
    
    if plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('numbelem', 'expc'), \
                             ('expc', 'modl'), \
                            ])
    
    if plottype == 'lght0003':
        grap.add_edges_from([ \
                             ('spatdistcons', 'lgal'), \
                             ('spatdistcons', 'bgal'), \
                            ])
        
    labl = {}
    if plottype.startswith('lens'):
        nameelem = r'\rm{sub}'
    else:
        nameelem = r'\rm{pts}'
    if plottype.startswith('lght') and (plottype == 'lght0001' or plottype == 'lght0002'):
        labl['numbelem'] = r'$\vec{N}_{%s}$' % nameelem
        labl['meanelem'] = r'$\vec{\mu}_{%s}$' % nameelem
    else:
        labl['numbelem'] = '$N_{%s}$' % nameelem
        labl['meanelem'] = r'$\mu_{%s}$' % nameelem
    
    if plottype.startswith('lght'):
        if plottype == 'lght0000' or plottype == 'lght0003':
            labl['amplslop'] = r'$\alpha$'
        else:
            labl['amplslop'] = r'$\vec{\alpha}$'
    if plottype.startswith('lens'):
        labl['defsslop'] = r'$\beta$'
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        labl['sinddistmean'] = r'$\vec{\beta}$'
    
    if plottype == 'lght0003':
        labl['spatdistcons'] = r'$\gamma$'
    if plottype.startswith('lens'):
        labl['lenp'] = r'$\vec{\chi}$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['bacp'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{\theta_1}$'
    labl['bgal'] = r'$\vec{\theta_2}$'
    if plottype.startswith('meta'):
        labl['feat'] = r'$\vec{\xi}$'
    else:
        if plottype.startswith('lght'):
            labl['sind'] = r'$\vec{s}$'
            labl['ampl'] = r'$\vec{f}$'
        else:
            labl['defs'] = r'$\vec{\alpha_{\rm{s}}}$'
    if plottype == 'lens0001':
        labl['asca'] = r'$\vec{\theta_{\rm{s}}}$'
        labl['acut'] = r'$\vec{\theta_{\rm{c}}}$'
        
    if plottype == 'lght0002':
        labl['expc'] = r'$\vec{E_{\rm{c}}}$'
    labl['modl'] = r'$M_D$'
    labl['data'] = r'$D$'
    
    posi = nx.circular_layout(grap)
    posi['sinddistmean'] = np.array([0.4, 0.15])
    if plottype == 'lght0003':
        posi['spatdistcons'] = np.array([-0.2, 0.15])
    if plottype.startswith('lght'):
        posi['numbelem'] = np.array([0., 0.075])
        posi['meanelem'] = np.array([0., 0.15])
        posi['amplslop'] = np.array([0.2, 0.15])
    if plottype.startswith('lens'):
        posi['numbelem'] = np.array([-0.1, 0.075])
        posi['meanelem'] = np.array([-0.1, 0.15])
        posi['defsslop'] = np.array([0.1, 0.15])
    
    if plottype.startswith('lght'):
        if plottype == 'lght0002':
            posi['psfp'] = np.array([0.7, -0.0])
            posi['bacp'] = np.array([0.9, -0.0])
        else:
            posi['psfp'] = np.array([0.5, -0.0])
            posi['bacp'] = np.array([0.7, -0.0])
    if plottype == 'lens0000':
        posi['psfp'] = np.array([0.3, -0.0])
        posi['bacp'] = np.array([0.5, -0.0])
        posi['lenp'] = np.array([0.7, -0.0])
    if plottype == 'lens0001':
        posi['psfp'] = np.array([0.7, -0.0])
        posi['bacp'] = np.array([0.9, -0.0])
        posi['lenp'] = np.array([1.1, -0.0])
    posi['lgal'] = np.array([-0.3, -0.0])
    posi['bgal'] = np.array([-0.1, -0.0])
    if plottype.startswith('lght'):
        posi['ampl'] = np.array([0.1, -0.0])
        posi['sind'] = np.array([0.3, -0.0])
    if plottype == 'lght0002':
        posi['expc'] = np.array([0.5, -0.0])

    if plottype.startswith('lens'):
        posi['defs'] = np.array([0.1, -0.0])
    if plottype == 'lens0001':
        posi['asca'] = np.array([0.3, -0.0])
        posi['acut'] = np.array([0.5, -0.0])
    posi['modl'] = np.array([0., -0.075])
    posi['data'] = np.array([0., -0.15])
   
    if typeverb > 0:
        numb = max(len(grap.edges()), len(listcolr))
        for k in range(numb):
            try:
                print('%15s %15s %15s' % (grap.edges()[k][0], grap.edges()[k][1], listcolr[k]))
            except:
                print('unequal')

    size = 1000
    nx.draw(grap, posi, labels=labl, ax=axis, edgelist=[], nodelist=[])
    nx.draw_networkx_edges(grap, posi, ax=axis, labels=labl, edge_color=listcolr)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['modl', 'data'], node_color='grey', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbelem'], node_color='b', node_size=size)
    if plottype.startswith('lght'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'amplslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind'], node_color='g', node_size=size)
    if plottype == 'lght0001' or plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['sinddistmean'], node_color='r', node_size=size)
    if plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['expc'], node_color='g', node_size=size)
    if plottype == 'lght0003':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['spatdistcons'], node_color='r', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['psfp', 'bacp'], node_color='y', node_size=size)
    if plottype.startswith('lens'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'defsslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lenp'], node_color='y', node_size=size)
    if plottype == 'lens0000':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs'], node_color='g', node_size=size)
    if plottype == 'lens0001':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs', 'asca', 'acut'], node_color='g', node_size=size)
    
    pathplot = pathpcat + '/imag/'
    plt.tight_layout()
    figr.savefig(pathplot + 'grap%s.pdf' % plottype)
    plt.close(figr)


def plot_3fgl_thrs(gdat):

    path = pathpcat + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = astropy.io.fits.getdata(path, 0)

    bgalfgl3 = np.linspace(-90., 90., 481)
    lgalfgl3 = np.linspace(-180., 180., 960)

    bgalexpo = np.linspace(-90., 90., 400)
    lgalexpo = np.linspace(-180., 180., 800)

    #fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)
    fluxthrs = griddata([lgalfgl3, bgalfgl3], fluxthrs, [gdat.lgalheal])

    cntsthrs = fluxthrs * gdat.expo

    jbgal = np.where(abs(bgalexpo) < 10.)[0]
    jlgal = np.where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel(gdat.labllgaltotl)
    axis.set_ylabel(gdat.lablbgaltotl)

    imag = plt.imshow(fluxthrs[np.amin(jbgal):np.amax(jbgal)+1, np.amin(jlghprofi):np.amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'thrs.pdf')
    plt.close(figr)
    

def plot_init(gdat):
    
    print('Making initial plots...')

    gmod = gdat.fitt

    # make initial plots
    if gdat.makeplot:
        
        if gmod.numbparaelem > 0:
            for l in gmod.indxpopl:
                if (gmod.typeelemspateval[l] == 'locl' and gmod.maxmpara.numbelem[l] > 0) and gdat.numbpixl > 1:
                    plot_indxprox(gdat)
        
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if gdat.typedata == 'mock' and gmod.boollens:
                figr, axis, path = init_figr(gdat, None, 'post', 'cntpmodlraww', 'this', 'true', i, m, -1)
                imag = retr_imag(gdat, axis, gmod.cntpmodlraww, 'this', 'true', 'cntpdata', i, m, booltdim=True)
                make_cbar(gdat, axis, imag, 0, tick=gdat.valutickmajrpara.cntpdata, labltotl=gdat.lablcntpdata)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)

    if gdat.boolcorrexpo:
        gdat.lablnumbpixl = r'$N_{\rm{pix}}$'
        gdat.limtexpo = [gdat.minmpara.expo, gdat.maxmpara.expo]
        if gdat.boolbinsener:
            path = gdat.pathinit + 'expototlmean.pdf'
            tdpy.plot_gene(path, gdat.meanpara.ener, gdat.expototlmean, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, \
                                                                                            lablydat=gdat.lablexpototl, limtydat=gdat.limtexpo)
        
        for m in gdat.indxevtt:
            for i in gdat.indxener:
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.hist(gdat.expo[i, :, m], gdat.binspara.expo)
                axis.set_xlabel(gdat.labltotlpara.expo)
                axis.set_ylabel(gdat.labltotlpara.numbpixl)
                axis.set_xscale('log')
                axis.set_yscale('log')
                plt.tight_layout()
                name = 'histexpo'
                if gdat.numbener > 1:
                    name += 'en%02d' % i
                if gdat.numbevtt > 1:
                    name += 'evt%d' % m
                path = gdat.pathinit + name + '.pdf'
                figr.savefig(path)
                plt.close(figr)
            
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    figr, axis, path = init_figr(gdat, None, 'post', 'expo', '', '', i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.expo, None, None, 'expo', i, m)
                    make_cbar(gdat, axis, imag, i)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
                

def plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                        strgvarb='defl', nameparagenrelem='', indxdefl=None, indxpoplplot=-1, multfact=1., indxenerplot=None, indxevttplot=None):

    if indxdefl is not None:
        strgvarb += 'sing'
    strgvarb = strgvarb + nameparagenrelem
    
    defl = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    defl *= multfact
   
    if indxenerplot is not None:
        defl = defl[indxenerplot, :, indxevttplot, ...]

    if indxdefl is not None:
        defl = defl[..., indxdefl]
        strgvarb += '%04d' % indxdefl
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))

    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgvarb, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot)
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    draw_frambndr(gdat, axis)
  
    defllgal = defl[:, :, 0]
    deflbgal = defl[:, :, 1]
    fact = 4
    ptch = axis.quiver(gdat.anglfact * gdat.lgalgridcart[::fact, ::fact], gdat.anglfact * gdat.bgalgridcart[::fact, ::fact], \
                       gdat.anglfact * defllgal[::fact, ::fact], gdat.anglfact * deflbgal[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis)
    plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    plt.subplots_adjust()
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgvarb, indxenerplot=None, indxevttplot=-1, strgcbar=None, \
                                                                booltdim=False, indxpoplplot=-1, strgmome='pmea'):
    
    gmod = getattr(gdat, strgmodl)
    
    if strgcbar is None:
        strgcbar = strgvarb
  
    # construct the string for the map
    if strgvarb == 'cntpdata':
        strgplot = strgvarb
    else:
        if strgstat == 'post':
            strgtemp = strgmome + strgpdfn
        else:
            strgtemp = ''
        strgplot = strgtemp + strgvarb
    
    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot)
   
    maps = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    imag = retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot, indxevttplot, booltdim=booltdim)
    
    make_cbar(gdat, axis, imag, strgvarb)
    
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    if gdat.boolsuprelem:
        supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot)

    print('strgvarb')
    print(strgvarb)
    plt.tight_layout()
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def init( \
         # user interaction
         ## type of verbosity
         typeverb=1, \

         ## path in which PCAT data lives
         pathpcat=None, \
        
         # miscelleneaous
         ## type of PDF to sample from
         strgpdfn='post', \

         # data
         ## type of data
         ### 'mock': simulated data
         ### 'inpt': input data
         ### 'real': real data retrieved from databases
         typedata=None, \
         ## type of experiment
         typeexpr='user', \

         # diagnostics
         ## Boolean to enter the diagnostic mode
         booldiagmode=True, \
         ## squeeze exposure to check the low sample limit         
         boolsqzeexpo=False, \
         ### explode exposure to check the large sample limit         
         boolexplexpo=False, \
         ## squeeze proposal scale to check the acceptance ratio
         boolsqzeprop=False, \
         ## explode proposal scale to check the acceptance ratio
         boolexplprop=False, \
         ## Boolean to thin down the data
         boolthindata=False, \
         ## factor by which to thin down the data
         factthin=None, \
    
         # reference catalog
         ## Boolean to use the reference catalogs to associate
         boolasscrefr=None, \
        
         # sampling
         ## Boolean flag to make burn-in tempered
         boolburntmpr=False, \
         ## number of sweeps
         numbswep=100000, \
         ## number of samples
         numbsamp=None, \
         ## number of initial sweeps to be burned
         numbburn=None, \
        
         # output
         ## Boolean to make condensed catalog
         boolcondcatl=True, \
        
         refrlabltotl=None, \
         refrlablpopl=None, \
         fittlablpopl=None, \
    
         # numpy RNG seed
         seedtype=0, \
         ## Boolean flag to re-seed each chain separately
         boolseedchan=True, \
         ## optional deterministic seed for sampling element parameters
         seedelem=None, \
         
         indxevttincl=None, \
         indxenerincl=None, \
        
         listmask=None, \
        
         # number of samples for Bootstrap
         numbsampboot=None, \

         listnamefeatsele=None, \
         
         # type of mask for the exposure map
         typemaskexpo='ignr', \
         
         # type of exposure
         ## 'cons': constant
         ## 'file': provided in a file
         typeexpo='cons', \

         # maximum spatial distance out to which element kernel will be evaluated
         maxmangleval=None, \
         
         # initial state
         initpsfprefr=False, \
         initpsfp=None, \
        
         # evaluate the likelihood inside circles around elements
         typeelemspateval=None, \
        
         namestattrue=None, \
        
         # plotting
         ## Boolean flag to make the frame plots short
         boolshrtfram=True, \
        
         boolrefeforc=False, \
         indxrefrforc=None, \

         ## Boolean to overplot the elements
         boolsuprelem=True, \
         
         ## Boolean to plot the correlation between elements
         boolplotelemcorr=True, \
         
         ## Boolean flag to vary the PSF
         boolmodipsfn=False, \

         # name of the configuration
         strgcnfg=None, \
        
         # model
         ## number of spatial dimensions
         numbspatdims=2, \
         # hyperparameters
         fittampldisttype=None, \
         # metamodel settings
         
         ## PSF evaluation type
         ## kernel evaluation type
         kernevaltype='ulip', \

         # photometric model
        
         ## base parameters
         ### Sersic type
         typesers='vauc', \
         ## transdimensional parameters (elements)
         ### vary projected scale radius
         variasca=True, \
         ### vary projected cutoff radius
         variacut=True, \

         # prior
         penalpridiff=False, \
         priotype='logt', \
         priofactdoff=None, \

         # initialization
         ## initialization type
         inittype=None, \
        
         loadvaripara=False, \
         
         # save the state of the MCMC
         savestat=False, \
         namesavestat=None, \
         # recover the state from a previous run
         namerecostat=None, \
         forcsavestat=False, \

         # proposals
         ## Boolean flag to turn on proposals on element parameters
         boolpropcomp=True, \
         boolpropcova=True, \
         propwithsing=True, \
         # type of covariance estimation
         typeopti='none', \
         
         # modes of operation
         ## only generate and plot mock data
         boolmockonly=False, \
         ## perform an additional run sampling from the prior
         checprio=False, \

         strgexprsbrt=None, \
         anglassc=None, \
         nameexpr=None, \
         
         # likelihood dependent
         ## exposure map
         expo=None, \

         lgalprio=None, \
         bgalprio=None, \
         minmcntpdata=None, \
         strgexpo=None, \
         
         # number of processors
         numbproc=None, \
         
         # likelihood function
         liketype='pois', \
         # user-defined likelihood function
         retr_llik=None, \
         
         anlytype=None, \
         
         lgalcntr=0., \
         bgalcntr=0., \
         
         maxmangl=None, \
         
         # spatial grid
         ## type of spatial pixelization
         typepixl=None, \
         ## Boolean flag to force Cartesian spatial grid
         boolforccart=False, \
         # number of pixels on a side in the Cartesian grid
         numbsidecart=None, \
         # Nside in Healpix
         numbsideheal=256, \
         
         allwfixdtrue=True, \
         asscmetrtype='dist', \

         # plotting
         numbswepplot=None, \
         # Boolean flagt to make the frame plots only for the central energy and PSF bin
         boolmakeframcent=True, \
         makeplot=True, \
         makeplotinit=True, \
         makeplotfram=True, \
         makeplotfinlprio=True, \
         makeplotfinlpost=True, \
         
         makeplotintr=False, \
         scalmaps='asnh', \
         makeanim=True, \
         strgenerfull=None, \
         strgexprname=None, \
         strganglunit=None, \
         strganglunittext=None, \
         anglfact=None, \
            
         limtydathistfeat=None, \
            
         # model
         # emission
         ## elements

         ## PSF
         specfraceval=None, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=100, \
    
         listprefsbrtsbrt=None, \
         listprefsbrtener=None, \
         listprefsbrtlabltotl=None, \

         lablgangunit=None, \
         labllgal=None, \
         lablbgal=None, \
         lablfluxunit=None, \
         lablflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxevttfull=None, \
         binsenerfull=None, \
         asymfluxprop=False, \
         
         ## Boolean flag to make the PSF model informed
         boolpriopsfninfo=False, \
         
         ## spectral

         # lensing
         fittrelnpowr=0., \

         # temp
         margfactmodl=1., \
         maxmgangdata=None, \
        
         # proposals
         stdvprophypr=0.01, \
         stdvproppsfp=0.1, \
         stdvpropbacp=0.01, \
         stdvproplenp=1e-4, \
         stdvlgal=0.001, \
         stdvbgal=0.001, \
         stdvflux=0.001, \
         stdvspep=0.001, \
         stdvspmrsind=0.2, \
         varistdvlbhl=True, \
        
         rtagmock=None, \
        
         ## transdimensional proposal probabilities
         probtran=None, \
         probspmr=None, \
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
         # standard deviation of the Gaussian from which the angular splitting will be drawn for splits and merges
         radispmr=None, \
        
         defa=False, \
         **args \
        ):

    # preliminary setup
    # construct the global object 
    gdat = tdpy.gdatstrt()
    for attr, valu in locals().items():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)
    
    # copy all provided inputs to the global object
    for strg, valu in args.items():
        setattr(gdat, strg, valu)

    # PCAT folders
    if gdat.pathpcat is None:
        gdat.pathpcat = os.environ["PCAT_DATA_PATH"] + '/'
        
    if gdat.pathpcat[-1] != '/':
        gdat.pathpcat += '/'
    gdat.pathdata = gdat.pathpcat + 'data/'
    gdat.pathdataopti = gdat.pathdata + 'opti/'
    gdat.pathimag = gdat.pathpcat + 'imag/'
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathinpt = gdat.pathdata + 'inpt/'
        
    # list of parameter groups
    gdat.liststrggroppara = ['genrbase', 'genrelem', 'derifixd', 'derielem', 'genrelemextd', 'derielemextd', 'kind', 'full']
    
    # list of parameter features to be turned into lists
    gdat.listfeatparalist = ['minm', 'maxm', 'fact', 'scal', 'lablroot', 'lablunit', 'stdv', 'labltotl', 'name']
    # list of parameter features
    gdat.listfeatpara = gdat.listfeatparalist + ['limt', 'bins', 'delt', 'numb', 'indx', 'cmap', 'mean', 'tick', 'numbbins', 'valutickmajr', 'labltickmajr', 'valutickminr', 'labltickminr']
        
    # run tag
    gdat.strgswep = '%d' % (gdat.numbswep)
    
    ## time stamp
    gdat.strgtimestmp = tdpy.retr_strgtimestmp()
    
    ## name of the configuration function
    if gdat.strgcnfg is None:
        gdat.strgcnfg = inspect.stack()[1][3]
   
    gdat.strgvers = 'v0.3'
    if gdat.typeverb > 0:
        print('PCAT %s started at %s.' % (gdat.strgvers, gdat.strgtimestmp))
        print('Configuration %s' % gdat.strgcnfg)
    
    # string describing the number of sweeps
    gdat.strgnumbswep = '%d' % gdat.numbswep
    
    # output paths
    gdat.rtag = retr_rtag(gdat.strgcnfg, gdat.strgnumbswep)
    gdat.pathoutprtag = retr_pathoutprtag(gdat.pathpcat, gdat.rtag)

    # physical constants
    gdat.prsccmtr = 3.086e18
    gdat.ergsgevv = 624.151
    gdat.factnewtlght = 2.09e13 # Msun / pc
        
    gdat.listnamepdir = ['forw', 'reve']
    gdat.listlablpdir = ['f', 'r']
        
    # number of standard deviations around mean of Gaussian-distributed variables
    gdat.numbstdvgaus = 4.
    
    # start the timer
    gdat.timerealtotl = time.time()
    gdat.timeproctotl = time.clock()
   
    # list of parameter types
    ## 'genr': generative parameters
    ## 'deri': derived parameters
    gdat.liststrgtypepara = ['genr', 'deri']
    
    booltemp = chec_statfile(gdat.pathpcat, gdat.rtag, 'gdatmodi')
    if booltemp:
        print('gdatmodi already exists. Skipping...')
    else:
    
        # create output folder for the run
        os.system('mkdir -p %s' % gdat.pathoutprtag)

        # write the list of arguments to file
        fram = inspect.currentframe()
        listargs, temp, temp, listargsvals = inspect.getargvalues(fram)
        fileargs = open(gdat.pathoutprtag + 'cmndargs.txt', 'w')
        fileargs.write('PCAT call arguments\n')
        for args in listargs:
            fileargs.write('%s = %s\n' % (args, listargsvals[args]))
        fileargs.close()
        
        # write the list of arguments to file
        fileargs = open(gdat.pathoutprtag + 'args.txt', 'w')
        fileargs.write('PCAT call arguments\n')
        for args in listargs:
            fileargs.write('%20s %s\n' % (args, listargsvals[args]))
        fileargs.close()
        
        # defaults
        if gdat.typedata is None:
            if gdat.strgexprsbrt is None:
                gdat.typedata = 'mock'
            else:
                gdat.typedata = 'inpt'
        print('gdat.typedata')
        print(gdat.typedata)

        # list of models
        gdat.liststrgmodl = []
        if gdat.typedata == 'mock':
            gdat.liststrgmodl += ['true']
        gdat.liststrgmodl += ['fitt']
        
        gdat.refr = tdpy.gdatstrt()
        
        gdat.listgmod = []
        for strgmodl in gdat.liststrgmodl + ['refr']:
            setattr(gdat, strgmodl, tdpy.gdatstrt())
            gmod = getattr(gdat, strgmodl)
            for strgstat in ['this', 'next']:
                setattr(gmod, strgstat, tdpy.gdatstrt())
            for strgfeatpara in gdat.listfeatpara:
                setattr(gmod, strgfeatpara + 'para', tdpy.gdatstrt())
        
            gdat.listgmod += [gmod]

        for strgfeatpara in gdat.listfeatpara:
            setattr(gdat, strgfeatpara + 'para', tdpy.gdatstrt())
        
        ## number of processes
        gdat.strgproc = os.uname()[1]
        if gdat.numbproc is None:
            if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu' or gdat.strgproc == 'wise':
                gdat.numbproc = 1
            else:
                gdat.numbproc = 1
    
        if gdat.typedata == 'inpt' and gdat.rtagmock is not None:
            print('Will use %s to account for selection effects.' % gdat.rtagmock)
            gdat.pathoutprtagmock = retr_pathoutprtag(gdat.pathpcat, gdat.rtagmock)

        ## number of burned sweeps
        if gdat.numbburn is None:
            print('gdat.numbswep')
            print(gdat.numbswep)
            gdat.numbburn = int(gdat.numbswep / 10)
            print('gdat.numbburn')
            print(gdat.numbburn)
    
        # burn-in
        gdat.factburntmpr = 0.75
        gdat.numbburntmpr = gdat.factburntmpr * gdat.numbburn
        
        if (gdat.boolsqzeprop or gdat.boolexplprop) and gdat.typeopti == 'hess':
            raise Exception('')

        print('gdat.boolpriopsfninfo')
        print(gdat.boolpriopsfninfo)
        
        print('gdat.typeexpr')
        print(gdat.typeexpr)
        
        ## factor by which to thin the sweeps to get samples
        if gdat.factthin is not None and gdat.numbsamp is not None:
            raise Exception('Both factthin and numbparagenrfull cannot be provided at the same time.')
        elif gdat.factthin is None and gdat.numbsamp is None:
            gdat.factthin = int(np.ceil(1e-3 * (gdat.numbswep - gdat.numbburn)))
            gdat.numbsamp = int((gdat.numbswep - gdat.numbburn) / gdat.factthin)
        elif gdat.numbsamp is not None:
            gdat.factthin = int((gdat.numbswep - gdat.numbburn) / gdat.numbsamp)
        elif gdat.factthin is not None:
            gdat.numbsamp = int((gdat.numbswep - gdat.numbburn) / gdat.factthin)
        if not isinstance(gdat.numbsamp, int) or not isinstance(gdat.factthin, int) or \
                                        not isinstance(gdat.numbburn, int) or not isinstance(gdat.numbswep, int):
            print('gdat.numbsamp')
            print(gdat.numbsamp)
            print('gdat.factthin')
            print(gdat.factthin)
            print('gdat.numbburn')
            print(gdat.numbburn)
            print('gdat.numbswep')
            print(gdat.numbswep)
            raise Exception('Number of samples is not an integer.')

        # samples to be saved
        gdat.indxsamp = np.arange(gdat.numbsamp)
        
        # samples to be saved from all chains
        gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
        gdat.indxsamptotl = np.arange(gdat.numbsamptotl)
        gdat.numbsweptotl = gdat.numbswep * gdat.numbproc
        
        if gdat.typeverb > 0:
            print('%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % \
                                                                                                    (gdat.numbswep, gdat.numbburn, gdat.factthin))
            print('The resulting chain will contain %d samples per chain and %d samples in total.' % (gdat.numbsamp, gdat.numbsamptotl))

        if gdat.anlytype is None:
            if gdat.typeexpr == 'chan':
                gdat.anlytype = 'home'
            elif gdat.typeexpr == 'ferm':
                gdat.anlytype = 'rec8pnts'
            else:
                gdat.anlytype = 'nomi'
        
        if gdat.priofactdoff is None:
            gdat.priofactdoff = 1.
        
        # experiment defaults
        if gdat.typeexpr == 'ferm':
            gdat.lablenerunit = 'GeV'
        if gdat.typeexpr == 'chan':
            gdat.lablenerunit = 'keV'
        if gdat.typeexpr == 'gene':
            gdat.lablenerunit = ''
        if gdat.typeexpr == 'fire':
            gdat.lablenerunit = '$\mu$m^{-1}'
        
        if gdat.typeexpr == 'ferm':
            if gdat.anlytype[4:8] == 'pnts':
                bins = np.logspace(np.log10(0.3), np.log10(10.), 4)
            if gdat.anlytype[4:8] == 'back':
                bins = np.logspace(np.log10(0.3), np.log10(300.), 31)
        if gdat.typeexpr == 'chan':
            if gdat.anlytype.startswith('home'):
                bins = np.array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
            if gdat.anlytype.startswith('extr'):
                bins = np.array([0.5, 2., 8.])
            if gdat.anlytype.startswith('spec'):
                bins = np.logspace(np.log10(0.5), np.log10(10.), 21)
        if gdat.typeexpr == 'fire':
            bins = np.logspace(np.log10(1. / 2.5e-6), np.log10(1. / 0.8e-6), 31)
        if gdat.typeexpr == 'hubb':
            # temp
            #bins = np.array([500., 750, 1000.])
            bins = np.array([750, 1000.])
        if gdat.typeexpr != 'gene':
            setp_varb(gdat, 'enerfull', bins=bins)
        
        setp_varb(gdat, 'numbpixl', lablroot='$N_{pix}$')
        
        if gdat.expo is not None:
            setp_varb(gdat, 'expo', minm=np.amin(gdat.expo), maxm=np.amax(gdat.expo), lablroot='$\epsilon$', cmap='OrRd', scal='logt')
        
        # energy band string
        if gdat.strgenerfull is None:
            if gdat.typeexpr == 'tess':
                gdat.strgenerfull = ['T']
            if gdat.typeexpr == 'sdss':
                gdat.strgenerfull = ['z-band', 'i-band', 'r-band', 'g-band', 'u-band']
            if gdat.typeexpr == 'hubb':
                #gdat.strgenerfull = ['F606W', 'F814W']
                gdat.strgenerfull = ['F814W']
            if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'fire': 
                gdat.strgenerfull = []
                for i in range(len(gdat.binspara.enerfull) - 1):
                    gdat.strgenerfull.append('%.3g %s - %.3g %s' % (gdat.binspara.enerfull[i], gdat.lablenerunit, gdat.binspara.enerfull[i+1], gdat.lablenerunit))
            if gdat.typeexpr == 'gene':
                gdat.strgenerfull = ['']
        
        ## PSF class
        if gdat.indxevttfull is None:
            if gdat.typeexpr == 'ferm':
                gdat.indxevttfull = np.arange(2)
            else:
                gdat.indxevttfull = np.arange(1)
        
        if gdat.indxevttincl is None:
            if gdat.typeexpr == 'ferm':
                gdat.indxevttincl = np.array([0, 1])
            else:
                gdat.indxevttincl = np.arange(1)
        
        if gdat.indxevttincl is not None:
            gdat.evttbins = True
        else:
            gdat.evttbins = False
        if gdat.evttbins:
            gdat.numbevtt = gdat.indxevttincl.size
            gdat.numbevttfull = gdat.indxevttfull.size
        else:
            gdat.numbevtt = 1
            gdat.numbevttfull = 1
            gdat.indxevttincl = np.array([0])
        gdat.indxevtt = np.arange(gdat.numbevtt)
        
        # Boolean flag to indicate that the data are binned in energy
        if gdat.typeexpr == 'gene':
            gdat.boolbinsener = False
        else:
            gdat.boolbinsener = True
        
        if gdat.boolbinsener:
            gdat.numbenerfull = len(gdat.strgenerfull)
        else:
            gdat.numbenerfull = 1
        gdat.indxenerfull = np.arange(gdat.numbenerfull)

        if gdat.typepixl is None:
            if gdat.typeexpr == 'ferm':
                gdat.typepixl = 'heal'
            else:
                gdat.typepixl = 'cart'
        
        if gdat.boolbinsener:
            gdat.meanpara.enerfull = np.sqrt(gdat.binspara.enerfull[1:] * gdat.binspara.enerfull[:-1])
        
        setp_varb(gdat, 'boolmodipsfn', valu=False, strgmodl='fitt')
        
        # default values for model types
        print('Starting to determine the default values for model types using setp_varbvalu()...')
        if gdat.typeexpr == 'hubb':
            typeemishost = 'sers'
        else:
            typeemishost = 'none'
        setp_varb(gdat, 'typeemishost', valu=typeemishost)

        setp_varb(gdat, 'lliktotl', lablroot='$L$')
        
        ### background type
        #### template
        if gdat.typeexpr == 'ferm':
            if gdat.anlytype == 'bfun':
                gdat.ordrexpa = 10
                gdat.numbexpasing = gdat.ordrexpa**2
                gdat.numbexpa = gdat.numbexpasing * 4
                gdat.indxexpa = np.arange(gdat.numbexpa)
                typeback = ['bfun%04d' % k for k in gdat.indxexpa]
            else:
                typeback = [1., 'sbrtfdfmsmthrec8pntsnorm.fits']
        if gdat.typeexpr == 'chan':
            # particle background
            if gdat.anlytype.startswith('spec'):
                # temp -- this is fake!
                sbrtparttemp = np.array([70.04, 70.04, 12.12, 15.98, 10.79, 73.59, 73.59])
                binsenerpart = np.logspace(np.log10(0.5), np.log10(10.), 6)
                meanenerpart = np.sqrt(binsenerpart[:-1] * binsenerpart[1:])
                meanenerparttemp = np.concatenate((np.array([0.5]), meanenerpart, np.array([10.])))
                typebacktemp = interp(gdat.meanpara.enerfull, meanenerparttemp, sbrtparttemp)
            if gdat.anlytype.startswith('home') :
                typebacktemp = 1.
                #typebacktemp = np.array([70.04, 12.12, 15.98, 10.79, 73.59]) / 70.04
            if gdat.anlytype.startswith('extr'):
                #typebacktemp = 'sbrtchanback' + gdat.anlytype + '.fits'
                typebacktemp = 1.
            
            if gdat.anlytype.startswith('spec'):
                typeback = [[1e2, 2.], typebacktemp]
            else:
                typeback = [1., typebacktemp]
        
        if gdat.typeexpr == 'hubb':
            typeback = [1.]
        if gdat.typeexpr == 'tess':
            typeback = [1.]
        if gdat.typeexpr == 'gene':
            typeback = [1.]
        if gdat.typeexpr == 'fire':
            typeback = [1.]
        if gdat.typeexpr != 'user':
            setp_varb(gdat, 'typeback', valu=typeback)
        
        if gdat.typeexpr == 'hubb':
            numbsersfgrd = 1
        else:
            numbsersfgrd = 0
        setp_varb(gdat, 'numbsersfgrd', valu=numbsersfgrd)
        
        if gdat.typeexpr == 'gene':
            typeelem = ['clus']
        if gdat.typeexpr == 'ferm':
            typeelem = ['lghtpnts']
        if gdat.typeexpr == 'tess':
            typeelem = ['lghtpnts']
        if gdat.typeexpr == 'chan':
            typeelem = ['lghtpnts']
        if gdat.typeexpr == 'hubb':
            typeelem = ['lghtpnts', 'lens', 'lghtgausbgrd']
        if gdat.typeexpr == 'fire':
            typeelem = ['lghtlineabso']
        if gdat.typeexpr == 'user':
            typeelem = ['user']
        setp_varb(gdat, 'typeelem', valu=typeelem)
        print('gdat.fitt.typeelem')
        print(gdat.fitt.typeelem)

        ### PSF model
        #### angular profile
        if gdat.typeexpr == 'ferm':
            typemodlpsfn = 'doubking'
        if gdat.typeexpr == 'chan':
            typemodlpsfn = 'singking'
        if gdat.typeexpr == 'sdss':
            typemodlpsfn = 'singgaus'
        if gdat.typeexpr == 'hubb':
            typemodlpsfn = 'singgaus'
        if gdat.typeexpr == 'tess':
            typemodlpsfn = 'singgaus'
        if gdat.typeexpr == 'gene':
            typemodlpsfn = 'singgaus'
        if gdat.typeexpr == 'fire':
            typemodlpsfn = None
        if gdat.typeexpr != 'user':
            setp_varb(gdat, 'typemodlpsfn', valu=typemodlpsfn)
        
        #### background names
        listnameback = ['isot']
        if gdat.typeexpr == 'ferm':
            listnameback.append('fdfm')
        #if gdat.typeexpr == 'chan':
        #    listnameback.append('part')
        setp_varb(gdat, 'listnameback', valu=listnameback)
        
        if gdat.strgpdfn == 'prio':
            gdat.lablsampdist = 'Prior'
        if gdat.strgpdfn == 'post':
            gdat.lablsampdist = 'Posterior'

        for strgmodl in gdat.liststrgmodl:
            # set up the indices of the model
            setp_indxpara(gdat, 'init', strgmodl=strgmodl)
   
        if gdat.numbswepplot is None:
            gdat.numbswepplot = 50000
   
        gdat.numbplotfram = gdat.numbswep / gdat.numbswepplot

        #setp_varb(gdat, 'colr', valu='mediumseagreen', strgmodl='refr')
        setp_varb(gdat, 'colr', valu='b', strgmodl='fitt')
        if gdat.typedata == 'mock':
            setp_varb(gdat, 'colr', valu='g', strgmodl='true')
        
        #gdat.refr.colr = 'mediumseagreen'
        #gdat.fitt.colr = 'deepskyblue'

        gdat.minmmass = 1.
        gdat.maxmmass = 10.
        
        if gdat.checprio:
            gdat.liststrgpdfn = ['prio', 'post']
        else:
            gdat.liststrgpdfn = ['post']

        gdat.lablmass = 'M'
        gdat.minmmassshel = 1e1
        gdat.maxmmassshel = 1e5
        gdat.lablmassshel = '$M_r$' 

        gdat.lablcurv = r'\kappa'
        gdat.lablexpc = r'E_{c}'
        
        gmod.scalcurvplot = 'self'
        gmod.scalexpcplot = 'self'
        
        #gdat.minmper0 = 1e-3 
        #gdat.maxmper0 = 1e1
        #
        #gdat.minmmagf = 10**7.5
        #gdat.maxmmagf = 10**16
        
        # temp -- automatize this eventually
        #gmod.minmper0 = gdat.minmper0
        #gmod.minmper0 = gdat.minmper0
        #gmod.maxmper0 = gdat.maxmper0
        #gmod.maxmper0 = gdat.maxmper0
        #gmod.minmmagf = gdat.minmmagf
        #gmod.minmmagf = gdat.minmmagf
        #gmod.maxmmagf = gdat.maxmmagf
        #gmod.maxmmagf = gdat.maxmmagf

        gdat.fitt.listelemmrkr = ['+', '_', '3']
        gdat.true.listmrkrhits = ['x', '|', '4']
        gdat.true.listmrkrmiss = ['s', 'o', 'p']
        gdat.true.listlablmiss = ['s', 'o', 'p']
        
        # list of scalings
        gdat.listscaltype = ['self', 'logt', 'atan', 'gaus', 'pois', 'expo']
        
        # number of grids
        gdat.numbgrid = 1
        gdat.indxgrid = np.arange(gdat.numbgrid)

        if gdat.typepixl == 'heal' and gdat.boolforccart:
            raise Exception('Cartesian forcing can only used with cart typepixl')

        gdat.liststrgphas = ['fram', 'finl', 'anim']
        gdat.liststrgelemtdimtype = ['bind']
        
        # lensing
        ## list of strings indicating different methods of calculating the subhalo mass fraction
        gdat.liststrgcalcmasssubh = ['delt', 'intg']
        
        # input data
        if gdat.typedata == 'inpt':
            path = gdat.pathinpt + gdat.strgexprsbrt
            gdat.sbrtdata = astropy.io.fits.getdata(path)
                
            if gdat.typepixl == 'heal' or gdat.typepixl == 'cart' and gdat.boolforccart:
                if gdat.sbrtdata.ndim != 3:
                    raise Exception('exprsbrtdata should be a 3D numpy np.array if pixelization is HealPix.')
            else:
                if gdat.sbrtdata.ndim != 4:
                    raise Exception('exprsbrtdata should be a 4D numpy np.array if pixelization is Cartesian.')
            
            if gdat.typepixl == 'cart' and not gdat.boolforccart:
                gdat.sbrtdata = gdat.sbrtdata.reshape((gdat.sbrtdata.shape[0], -1, gdat.sbrtdata.shape[3]))
                        
            gdat.numbenerfull = gdat.sbrtdata.shape[0]
            if gdat.typepixl == 'heal':
                gdat.numbpixlfull = gdat.sbrtdata.shape[1]
            elif gdat.boolforccart:
                gdat.numbpixlfull = gdat.numbsidecart**2
            else:
                gdat.numbpixlfull = gdat.sbrtdata.shape[1] * gdat.sbrtdata.shape[2]
            gdat.numbevttfull = gdat.sbrtdata.shape[2]
            
            if gdat.typepixl == 'heal':
                # temp
                gdat.numbsidecart = 100
                gdat.numbsidecarthalf = int(gdat.numbsidecart / 2)
                gdat.numbsideheal = int(np.sqrt(gdat.numbpixlfull / 12))
    
        if gdat.typeexpr == 'hubb':
            gdat.hubbexpofact = 1.63050e-19
        
        if gdat.strgexpo is None:
            if gdat.typeexpr == 'ferm':
                gdat.strgexpo = 'expofermrec8pntsigal0256.fits'
        
        if gdat.typeexpo is None:
            if gdat.typeexpr == 'ferm':
                gdat.typeexpo = 'file'
            else:
                gdat.typeexpo = 'cons'
        
        print('strgexpo') 
        print(strgexpo)
        
        ## generative model
        # the factor to convert radians (i.e., internal angular unit of PCAT) to the angular unit that will be used in the output (i.e., plots and tables)
        if gdat.anglfact is None:
            if gdat.typeexpr == 'ferm':
                gdat.anglfact = 180. / np.pi
            if gdat.typeexpr == 'tess':
                gdat.anglfact = 60 * 180. / np.pi
            if gdat.typeexpr == 'sdss' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'hubb':
                gdat.anglfact = 3600 * 180. / np.pi
            if gdat.typeexpr == 'sche' or gdat.typeexpr == 'gene':
                gdat.anglfact = 1.
        
        if gdat.numbsidecart is not None and gdat.typepixl == 'cart' and not gdat.boolforccart and isinstance(strgexpo, str):
            raise Exception('numbsidecart argument should not be provided when strgexpo is a file name and pixelization is Cartesian.')
                    
        if gdat.typepixl == 'heal' or gdat.typepixl == 'cart' and gdat.boolforccart:
            if gdat.numbsidecart is None:
                gdat.numbsidecart = 100
        
        # exposure
        gdat.boolcorrexpo = gdat.expo is not None
        if gdat.typeexpo == 'cons':
            if gdat.typedata == 'mock':
                if gdat.numbsidecart is None:
                    gdat.numbsidecart = 100
            if gdat.typedata == 'mock':
                if gdat.typepixl == 'heal':
                    gdat.expo = np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.typepixl == 'cart':
                    gdat.expo = np.ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
            if gdat.typedata == 'inpt':
                gdat.expo = np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
        if gdat.typeexpo == 'file':
            path = gdat.pathinpt + gdat.strgexpo
            if gdat.typeverb > 0:
                print('Reading %s...' % path)
            gdat.expo = astropy.io.fits.getdata(path)
            
            if gdat.typepixl == 'cart':
                gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))
        
            if gdat.numbsidecart is None:
                # temp -- gdat.numbsidecart takes the value of the region 0
                if np.sqrt(gdat.expo.shape[1]) % 1. != 0.:
                    raise Exception('')
                gdat.numbsidecart = int(np.sqrt(gdat.expo.shape[1]))
        
        if gdat.typedata == 'mock':
            if gdat.typepixl == 'cart':
                gdat.numbpixlfull = gdat.numbsidecart**2
            if gdat.typepixl == 'heal':
                gdat.numbpixlfull = 12 * gdat.numbsideheal**2
        
        # initialization type
        if gdat.inittype is None:
            gdat.inittype = 'rand'

        if gdat.typeexpr != 'user':
            
            # Boolean flag to indicate binning in space
            gdat.boolbinsspat = gdat.numbpixlfull != 1

            print('gdat.boolbinsspat')
            print(gdat.boolbinsspat)
            
            if gdat.boolcorrexpo and np.amin(gdat.expo) == np.amax(gdat.expo) and not isinstance(gdat.strgexpo, float):
                raise Exception('Bad input exposure map.')
        
            if gdat.boolbinsspat:
                if gdat.typepixl == 'cart' and isinstance(gdat.strgexpo, float) and gdat.typedata == 'inpt':
                    if np.sqrt(gdat.sbrtdata.shape[1]) % 1. != 0.:
                        raise Exception('')
                    gdat.numbsidecart = int(np.sqrt(gdat.sbrtdata.shape[1]))
                
                gdat.numbsidecarthalf = int(gdat.numbsidecart / 2)

                if gdat.typepixl == 'cart':
                    gdat.numbpixlcart = gdat.numbsidecart**2
        
                ### spatial extent of the data
                if gdat.maxmgangdata is None:
                    if gdat.typeexpr == 'chan':
                        gdat.maxmgangdata = 0.492 / gdat.anglfact * gdat.numbsidecarthalf
                    if gdat.typeexpr == 'ferm':
                        gdat.maxmgangdata = 15. / gdat.anglfact
                    if gdat.typeexpr == 'tess':
                        gdat.maxmgangdata = 20. / gdat.anglfact
                    if gdat.typeexpr == 'hubb':
                        gdat.maxmgangdata = 2. / gdat.anglfact
                    if gdat.typeexpr == 'gene':
                        gdat.maxmgangdata = 1. / gdat.anglfact
                
                print('gdat.numbsidecart')
                print(gdat.numbsidecart)
                print('gdat.maxmgangdata')
                print(gdat.maxmgangdata)
        
                # pixelization
                if gdat.typepixl == 'cart':
                    gdat.apix = (2. * gdat.maxmgangdata / gdat.numbsidecart)**2
                if gdat.typepixl == 'heal':
                    temp, temp, temp, gdat.apix = tdpy.retr_healgrid(gdat.numbsideheal)
                gdat.sizepixl = np.sqrt(gdat.apix)
        
            # factor by which to multiply the y axis limits of the surface brightness plot
            if gdat.numbpixlfull == 1:
                gdat.factylimtbrt = [1e-4, 1e7]
            else:
                gdat.factylimtbrt = [1e-4, 1e3]

            # grid
            gdat.minmlgaldata = -gdat.maxmgangdata
            gdat.maxmlgaldata = gdat.maxmgangdata
            gdat.minmbgaldata = -gdat.maxmgangdata
            gdat.maxmbgaldata = gdat.maxmgangdata
            
            if gdat.typepixl == 'cart' and gdat.boolforccart:
                if gdat.typedata == 'inpt':
                    sbrtdatatemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                    for i in gdat.indxenerfull:
                        for m in gdat.indxevttfull:
                            sbrtdatatemp[i, :, m] = tdpy.retr_cart(gdat.sbrtdata[i, :, m], \
                                                            numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).flatten()
                    gdat.sbrtdata = sbrtdatatemp

                if gdat.boolcorrexpo:
                    expotemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                    for i in gdat.indxenerfull:
                        for m in gdat.indxevttfull:
                            expotemp[i, :, m] = tdpy.retr_cart(gdat.expo[i, :, m], \
                                                            numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).flatten()
                    gdat.expo = expotemp
        
        gdat.sdenunit = 'degr'

        gdat.factergskevv = 1.6e-9
        if gdat.typeexpr == 'ferm':
            gdat.listspecconvunit = [['en02', 'gevv']]
        if gdat.typeexpr == 'chan':
            gdat.listspecconvunit = [['en00', 'kevv'], ['en02', 'kevv'], ['en02', 'ergs'], ['en03', 'ergs', '0520', 0.5,  2.], \
                                                                                           ['en03', 'ergs', '0210',  2., 10.], \
                                                                                           ['en03', 'ergs', '0510', 0.5, 10.], \
                                                                                           ['en03', 'ergs', '0208',  2.,  8.], \
                                                                                           ['en03', 'ergs', '0508', 0.5,  8.], \
                                                                                           ['en03', 'ergs', '0207',  2.,  7.], \
                                                                                           ['en03', 'ergs', '0507', 0.5,  7.]]
        if gdat.typeexpr == 'hubb':
            gdat.listspecconvunit = [['en03', 'ergs']]
        if gdat.typeexpr == 'fire':
            gdat.listspecconvunit = [['en00', 'imum']]
        
        # temp
        #if gdat.typeexpr == 'chan' and (gdat.anlytype.startswith('home') or gdat.anlytype.startswith('extr')):
        #    gmod.lablpopl = ['AGN', 'Galaxy']

        if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'fire':
            gdat.enerdiff = True
        if gdat.typeexpr == 'hubb' or gdat.typeexpr == 'gene' or gdat.typeexpr == 'tess':
            gdat.enerdiff = False
        
        if gdat.indxenerincl is None:
            
            # default
            if gdat.boolbinsener:
                gdat.indxenerincl = np.arange(gdat.binspara.enerfull.size - 1)
            
            if gdat.typeexpr == 'ferm':
                if gdat.anlytype[4:8] == 'pnts':
                    gdat.indxenerincl = np.arange(3)
                if gdat.anlytype[4:8] == 'back':
                    gdat.indxenerincl = np.arange(30)
            if gdat.typeexpr == 'chan':
                if gdat.anlytype.startswith('home'):
                    gdat.indxenerincl = np.arange(5)
                if gdat.anlytype.startswith('extr'):
                    gdat.indxenerincl = np.arange(2)
            if gdat.typeexpr == 'hubb':
                gdat.indxenerincl = np.array([0])
                #gdat.indxenerincl = np.array([1])
                #gdat.indxenerincl = np.array([0, 1])
            if gdat.typeexpr == 'gene':
                gdat.indxenerincl = np.array([0])
        
        if gdat.indxenerincl is None:
            gdat.numbener = 1
        else:
            gdat.numbener = gdat.indxenerincl.size
        gdat.indxener = np.arange(gdat.numbener, dtype=int)
        
        if gdat.indxenerincl is None:
            gdat.indxenerincl = gdat.indxener
        
        if gdat.boolbinsener:
            gdat.indxenerinclbins = np.empty(gdat.numbener+1, dtype=int)
            gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
            gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
            gdat.indxenerpivt = 0
            gdat.numbenerplot = 100
            gdat.strgener = [gdat.strgenerfull[k] for k in gdat.indxenerincl]
            gdat.binspara.ener = gdat.binspara.enerfull[gdat.indxenerinclbins]
            gdat.meanpara.ener = np.sqrt(gdat.binspara.ener[1:] * gdat.binspara.ener[:-1])
            gdat.deltener = gdat.binspara.ener[1:] - gdat.binspara.ener[:-1]
            gdat.minmener = gdat.binspara.ener[0]
            gdat.maxmener = gdat.binspara.ener[-1]
            retr_axis(gdat, 'ener')

            gdat.limtener = [np.amin(gdat.binspara.ener), np.amax(gdat.binspara.ener)] 
        if gdat.boolbinsener: 
            if gdat.numbener > 1:
                gdat.enerpivt = gdat.meanpara.ener[gdat.indxenerpivt]
            # energy bin indices other than that of the pivot bin
            gdat.indxenerinde = np.setdiff1d(gdat.indxener, gdat.indxenerpivt)
        
            # temp
            if gdat.typeexpr == 'chan':
                gdat.edis = 0.3 * np.sqrt(gdat.binspara.ener) / 2.35
                gdat.edisintp = sp.interpolate.interp1d(gdat.binspara.ener, gdat.edis, fill_value='extrapolate')
            else:
                gdat.edis = None
                gdat.edisintp = None

        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
        
            setp_varb(gdat, 'cntpmodl', lablroot='$C_{M}$', scal='asnh', strgmodl=strgmodl)

            # number of elements
            if strgmodl == 'true':
                for l in gmod.indxpopl:
                    if gmod.typeelem[l] == 'lens':
                        numbelem = 25
                    else:
                        numbelem = 5
                    setp_varb(gdat, 'numbelem', minm=0, maxm=10, lablroot='N', scal='pois', valu=numbelem, popl=l, strgmodl=strgmodl, strgstat='this')
            if strgmodl == 'fitt':
                setp_varb(gdat, 'numbelem', minm=0, maxm=10, lablroot='N', scal='pois', popl='full', strgmodl=strgmodl)

            ## hyperparameters
            setp_varb(gdat, 'typemodltran', valu='drct', strgmodl=strgmodl)
            
            if gmod.typemodltran == 'pois':
                setp_varb(gdat, 'meanelem', minm=0.1, maxm=1000., scal='logt', popl='full', strgmodl=strgmodl)
        
            #### boolean flag background
            if gdat.typeexpr != 'user':
                if gdat.typeexpr == 'chan':
                    if gdat.numbpixlfull == 1:
                        boolspecback = [True, True]
                    else:
                        boolspecback = [False, False]
                else:
                    boolspecback = [False for k in gmod.indxback]
                setp_varb(gdat, 'boolspecback', valu=boolspecback, strgmodl=strgmodl)
        
            typeelemspateval = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                # these  element types slow down execution!
                if gmod.typeelem[l] == 'lens' or gmod.typeelem[l].startswith('lghtline') or gmod.typeelem[l] == 'clusvari' or gmod.typeelem[l] == 'lghtgausbgrd':
                    typeelemspateval[l] = 'full'
                else:
                    typeelemspateval[l] = 'locl'
            setp_varb(gdat, 'typeelemspateval', valu=typeelemspateval, strgmodl=strgmodl)
            
            gmod.minmpara.numbelem = np.empty(gmod.numbpopl, dtype=int)
            gmod.maxmpara.numbelem = np.empty(gmod.numbpopl, dtype=int)
            for l in gmod.indxpopl:
                gmod.maxmpara.numbelem[l] = int(getattr(gmod.maxmpara, 'numbelempop%d' % l))
                gmod.minmpara.numbelem[l] = int(getattr(gmod.minmpara, 'numbelempop%d' % l))
            gmod.maxmpara.numbelemtotl = np.sum(gmod.maxmpara.numbelem)
            gmod.minmpara.numbelemtotl = np.sum(gmod.minmpara.numbelem)
            
            # spatial distribution type
            typespatdist = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                typespatdist[l] = 'unif'
            setp_varb(gdat, 'typespatdist', valu=typespatdist, strgmodl=strgmodl)
            
            # flux distribution type
            typeprioflux = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                # temp -- this can assign powr to populations whose flux is not drawn from a power law!
                if gmod.typeelem[l].startswith('lght'):
                    typeprioflux[l] = 'powr'
                else:
                    typeprioflux[l] = None
            setp_varb(gdat, 'typeprioflux', valu=typeprioflux, strgmodl=strgmodl)
        
        if gdat.strgexprname is None:
            if gdat.typeexpr == 'chan':
                gdat.strgexprname = 'Chandra'
            if gdat.typeexpr == 'ferm':
                gdat.strgexprname = 'Fermi-LAT'
            if gdat.typeexpr == 'hubb':
                gdat.strgexprname = 'HST'
            if gdat.typeexpr == 'sche':
                gdat.strgexprname = 'XXXXX'
            if gdat.typeexpr == 'gene':
                gdat.strgexprname = 'TGAS-RAVE'
        
        if gdat.lablgangunit is None:
            if gdat.typeexpr == 'ferm':
                gdat.lablgangunit = '$^o$'
            if gdat.typeexpr == 'gene':
                gdat.lablgangunit = ''
            if gdat.typeexpr == 'sdss' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'hubb':
                gdat.lablgangunit = '$^{\prime\prime}$'
        
        if gdat.labllgal is None:
            if gdat.typeexpr == 'gene':
                gdat.labllgal = r'L_{z}'
            else:
                if gdat.typeexpr == 'ferm' and gdat.lgalcntr == 0 and gdat.bgalcntr == 0:
                    gdat.labllgal = r'l'
                else:
                    gdat.labllgal = r'\theta_1'
        if gdat.lablbgal is None:
            if gdat.typeexpr == 'gene':
                gdat.lablbgal = r'E_k'
            else:
                if gdat.typeexpr == 'ferm' and gdat.lgalcntr == 0 and gdat.bgalcntr == 0:
                    gdat.lablbgal = r'b'
                else:
                    gdat.lablbgal = r'\theta_2'

        if gdat.strgenerunit is None:
            if gdat.typeexpr == 'ferm':
                gdat.strgenerunit = 'GeV'
                gdat.nameenerunit = 'gevv'
            if gdat.typeexpr == 'chan':
                gdat.strgenerunit = 'keV'
                gdat.nameenerunit = 'kevv'
            if gdat.typeexpr == 'gene':
                gdat.strgenerunit = ''
                gdat.nameenerunit = ''
            if gdat.typeexpr == 'hubb':
                gdat.strgenerunit = 'erg'
                gdat.nameenerunit = 'ergs'
            if gdat.typeexpr == 'fire':
                gdat.strgenerunit = '$\mu$ m$^{-1}$'
                gdat.nameenerunit = 'imum'

        if gdat.nameexpr is None:
            if gdat.typeexpr == 'ferm':
                gdat.nameexpr = 'Fermi-LAT'
            if gdat.typeexpr == 'sdss':
                gdat.nameexpr = 'SDSS'
            if gdat.typeexpr == 'chan':
                gdat.nameexpr = 'Chandra'
            if gdat.typeexpr == 'hubb':
                gdat.nameexpr = 'HST'
            if gdat.typeexpr == 'gaia':
                gdat.nameexpr = 'Gaia'
        
        ## Lensing
        if gdat.radispmr is None:
            if gdat.typeexpr == 'ferm':
                gdat.radispmr = 0.6 / gdat.anglfact
            if gdat.typeexpr == 'hubb':
                gdat.radispmr = 0.15 / gdat.anglfact
            if gdat.typeexpr == 'tess':
                gdat.radispmr = 1. / gdat.anglfact
            if gdat.typeexpr == 'chan':
                if gdat.anlytype == 'spec':
                    gdat.radispmr = 0.1
                else:
                    gdat.radispmr = 0.2 / gdat.anglfact
            if gdat.typeexpr == 'sdss':
                gdat.radispmr = 0.5 / gdat.anglfact
            if gdat.typeexpr == 'gene':
                gdat.radispmr = 0.2
        
        print('gdat.radispmr')
        print(gdat.radispmr)

        if gdat.anglassc is None:
            gdat.anglassc = 5. * gdat.radispmr
    
        print('gdat.anglassc')
        print(gdat.anglassc)

        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
        
            if gdat.boolbinsspat:
                if gdat.typeexpr == 'chan' or gdat.typeexpr == 'sdss':
                    numbpsfpform = 0
                    gmod.numbpsfptotl = 0
                if gdat.typeexpr == 'chan':
                    retr_psfpchan(gmod)
                if gdat.typeexpr == 'ferm':
                    retr_psfpferm(gmod)
                if gdat.typeexpr == 'sdss':
                    retr_psfpsdss(gmod)
                if gdat.typeexpr == 'hubb':
                    retr_psfphubb(gmod)
                if gdat.typeexpr == 'tess':
                    retr_psfptess(gmod)
                if gdat.typeexpr == 'gene':
                    retr_psfpsdyn(gmod)
        

        # model evaluation approximation error tolerance in units of the fraction of the lowest PS flux
        if gdat.specfraceval is None:
            if gdat.typeexpr == 'ferm':
                gdat.specfraceval = 0.5
            else:
                gdat.specfraceval = 0.1

        gdat.binspara.lgalcart = np.linspace(gdat.minmlgaldata, gdat.maxmlgaldata, gdat.numbsidecart + 1)
        gdat.binspara.bgalcart = np.linspace(gdat.minmbgaldata, gdat.maxmbgaldata, gdat.numbsidecart + 1)
        gdat.meanpara.lgalcart = (gdat.binspara.lgalcart[0:-1] + gdat.binspara.lgalcart[1:]) / 2.
        gdat.meanpara.bgalcart = (gdat.binspara.bgalcart[0:-1] + gdat.binspara.bgalcart[1:]) / 2.
        
        # reference elements
        gdat.numbrefr = 0
        if gdat.typedata == 'mock':
            gdat.numbrefr = gmod.numbpopl
        if gdat.typedata == 'inpt':
            if gdat.typeexpr == 'ferm':
                gdat.numbrefr = 2
            if gdat.typeexpr == 'chan':
                gdat.numbrefr = 2
        print('gdat.numbrefr')
        print(gdat.numbrefr)

        gdat.indxrefr = np.arange(gdat.numbrefr)
        if gdat.boolasscrefr is None:
            gdat.boolasscrefr = [True for q in gdat.indxrefr]
        
        gdat.listnamerefr = [] 
        gdat.refr.nameparagenrelemampl = [[] for q in gdat.indxrefr]
        gdat.refr.namepara.elem = [[] for q in gdat.indxrefr]
        gdat.refr.namepara.elemodim = [[] for q in gdat.indxrefr]
        gdat.boolinforefr = False
        gdat.listpathwcss = []
        
        gdat.numbpixllgalshft = []
        gdat.numbpixlbgalshft = []
        gdat.refrindxpoplassc = [[] for q in gdat.indxrefr] 
        
        # temp -- this allows up to 3 reference populations
        gdat.true.colrelem = ['darkgreen', 'olivedrab', 'mediumspringgreen']
        # temp -- this allows up to 3 reference populations
        gdat.fitt.colrelem = ['royalblue', 'dodgerblue', 'navy']
        if gdat.typedata == 'mock':
            gdat.boolinforefr = True
            gdat.listnamerefr = ['moc%d' % l for l in gmod.indxpopl] 
            gdat.indxrefr = np.arange(gdat.numbrefr)
        if gdat.typedata == 'inpt':
            if gdat.typeexpr == 'ferm':
                gdat.boolinforefr = True
                retr_refrferminit(gdat)
                for q in gdat.indxrefr:
                    gdat.refrindxpoplassc[q] = gmod.indxpopl
            if gdat.typeexpr == 'chan':
                gdat.boolinforefr = True
                retr_refrchaninit(gdat)
                for q in gdat.indxrefr:
                    gdat.refrindxpoplassc[q] = gmod.indxpopl
            
            for q in gdat.indxrefr:
                if 'lgal' in gdat.refr.namepara.elem[q] and 'bgal' in gdat.refr.namepara.elem[q]:
                    gdat.refr.namepara.elem[q] += ['gang', 'aang']
                for strgfeat in gdat.refr.namepara.elem[q]:
                    setattr(gdat.refr, strgfeat, [[] for q in gdat.indxrefr])
            
            if gdat.typeexpr == 'ferm':
                retr_refrfermfinl(gdat)
            if gdat.typeexpr == 'chan':
                retr_refrchanfinl(gdat)
        
        if gdat.typeexpr == 'hubb':
            boollenshost = True
        else:
            boollenshost = False
        setp_varb(gdat, 'boollenshost', valu=boollenshost)
  
        if gdat.typeexpr == 'hubb':
            boollenssubh = True
        else:
            boollenssubh = False
        setp_varb(gdat, 'boollenssubh', valu=boollenssubh)
  
        if gdat.typeexpr == 'hubb':
            boollens = True
        else:
            boollens = False
        setp_varb(gdat, 'boollens', valu=boollens)
  
        if gdat.typeexpr == 'hubb':
            boolemishost = True
        else:
            boolemishost = False
        setp_varb(gdat, 'boolemishost', valu=boolemishost)
  
        for strgmodl in gdat.liststrgmodl:
            
            gmod = getattr(gdat, strgmodl)

            ## names of the variables for which cumulative posteriors will be plotted
            if gmod.boollenssubh:
                gmod.listnamevarbcpct = ['convelem']
            else:
                gmod.listnamevarbcpct = []
        
        # the adis in the file is kpc
        fileh5py = h5py.File(gdat.pathdata + 'inpt/adis.h5','r')
        
        gdat.redsintp = fileh5py['reds'][()]
        gdat.adisintp = fileh5py['adis'][()] * 1e6 # [pc]

        gdat.adisobjt = sp.interpolate.interp1d(gdat.redsintp, gdat.adisintp, fill_value='extrapolate')

        gdat.redsfromdlosobjt = sp.interpolate.interp1d(gdat.adisintp * gdat.redsintp, gdat.redsintp, fill_value='extrapolate')
        fileh5py.close()
        
        setp_varb(gdat, 'lgal', minm=-10., maxm=10., lablroot='$l$')
        
        for strgmodl in gdat.liststrgmodl:
            
            gmod = getattr(gdat, strgmodl)
            
            if gdat.typedata == 'mock':
                if gmod.boollenshost:
                    setp_varb(gdat, 'redshost', valu=0.2, strgmodl='true')
                    setp_varb(gdat, 'redssour', valu=1., strgmodl='true')
                
                setp_indxpara(gdat, 'finl', strgmodl='true')
                
            ### background parameters
            if gdat.typeexpr == 'chan':
                if gdat.anlytype.startswith('extr'):
                    meanbacpbac1 = 1.
                else:
                    meanbacpbac1 = 70.04
                stdvbacpbac1 = 1e-5 * meanbacpbac1
                setp_varb(gdat, 'bacp', mean=meanbacpbac1, stdv=stdvbacpbac1, back=1, scal='gaus', strgmodl='true')

                if gdat.numbpixlfull == 1:
                    bacp = [1e0, 1e2]
                    setp_varb(gdat, 'bacp', limt=bacp, back=0)
                else:
                    bacp = [1e-1, 1e3]
                    setp_varb(gdat, 'bacp', limt=bacp, ener='full', back=0)
                if gdat.numbpixlfull == 1:
                    bacp = 10.
                    setp_varb(gdat, 'bacp', valu=bacp)
                else:
                    setp_varb(gdat, 'bacp', valu=170., back=0, ener=0)
                    setp_varb(gdat, 'bacp', valu=17.4, back=0, ener=1)
                    setp_varb(gdat, 'bacp', valu=27., back=0, ener=2)
                    setp_varb(gdat, 'bacp', valu=11.8, back=0, ener=3)
                    setp_varb(gdat, 'bacp', valu=101., back=0, ener=4)
            if gdat.typeexpr == 'ferm':
                if 'ferm_bubb' in gdat.strgcnfg:
                    setp_varb(gdat, 'bacp', limt=[1e-10, 1e10], ener='full', back='full')
                else:
                    # isotropic + unresolved
                    setp_varb(gdat, 'bacp', limt=[1e-7, 1e-2], ener=0, back=0)
                    setp_varb(gdat, 'bacp', limt=[1e-9, 1e-3], ener=1, back=0)
                    setp_varb(gdat, 'bacp', limt=[1e-10, 1e-4], ener=2, back=0)
                    # diffuse
                    setp_varb(gdat, 'bacp', limt=[1e-6, 1e-2], ener=0, back=1)
                    setp_varb(gdat, 'bacp', limt=[1e-7, 1e-3], ener=1, back=1)
                    setp_varb(gdat, 'bacp', limt=[1e-8, 1e-4], ener=2, back=1)
                    # dark
                    setp_varb(gdat, 'bacp', limt=[1e-11, 1e-4], ener=0, back=2)
                    setp_varb(gdat, 'bacp', limt=[1e-11, 1e-4], ener=1, back=2)
                    setp_varb(gdat, 'bacp', limt=[1e-11, 1e-4], ener=2, back=2)

                setp_varb(gdat, 'bacp', valu=5e-6, ener=0, back=0)
                setp_varb(gdat, 'bacp', valu=5e-6, ener=0, back=0)
                setp_varb(gdat, 'bacp', valu=2e-8, ener=1, back=0)
                setp_varb(gdat, 'bacp', valu=2e-9, ener=2, back=0)
                setp_varb(gdat, 'bacp', valu=1e-5, ener=4, back=0)
                setp_varb(gdat, 'bacp', valu=7e-7, ener=0, back=1)
                setp_varb(gdat, 'bacp', valu=1e-4, ener=0, back=1)
                setp_varb(gdat, 'bacp', valu=1e-5, ener=1, back=1)
                setp_varb(gdat, 'bacp', valu=7e-7, ener=2, back=1)
                setp_varb(gdat, 'bacp', valu=3e-8, ener=4, back=1)

                # Fourier basis
                for strgmodl in gdat.liststrgmodl:
                    for c in gmod.indxback:
                        if isinstance(typeback[c], str):
                            if 'bfun' in typeback[c]:
                                setp_varb(gdat, 'bacp', limt=[1e-10, 1e10], ener='full', back=c)

            if gdat.typeexpr == 'hubb':
                bacp = [1e-10, 1e-6]
            if gdat.typeexpr == 'gene':
                setp_varb(gdat, 'bacp', minm=1e-1, maxm=1e3, valu=1e1, lablroot='$A$', scal='logt', ener=0, back=0, strgmodl=strgmodl)
            
            if gdat.typeexpr == 'fire':
                bacp = [1e-1, 1e1]
            if gdat.typeexpr == 'tess':
                bacp = [1e-1, 1e1]
                setp_varb(gdat, 'bacp', limt=bacp, ener='full', back=0)
            
            if gdat.typeexpr == 'hubb':
                bacp = 2e-7
            if gdat.typeexpr == 'chan':
                bacp = 1.
                if gdat.numbpixlfull == 1:
                    setp_varb(gdat, 'bacp', valu=bacp, back=0)
                else:
                    setp_varb(gdat, 'bacp', valu=bacp, ener='full', back=0)

                    # particle background
                    if gdat.typeexpr == 'chan':
                        bacp = 70.04
                        setp_varb(gdat, 'bacp', valu=bacp, back=1)
                    
                    # particle background
                    #if gdat.typeexpr == 'chan':
                    #    if gdat.anlytype == 'spec':
                    #        bacp = [1e-8, 1e-6]
                    #    else:
                    #        bacp = [1e-1, 1e2]
                    #    setp_varb(gdat, 'bacp', limt=bacp, back=1)
                
            ### element parameter boundaries
            #### spatial
            if gdat.boolbinsspat:
                if gdat.typeexpr == 'ferm':
                    minmgang = 1e-1 / gdat.anglfact
                else:
                    minmgang = 1e-2 / gdat.anglfact
                setp_varb(gdat, 'minmgang', valu=minmgang, popl='full', strgmodl=strgmodl)
        
            # parameter defaults
            for l in gmod.indxpopl:
                if gmod.typeelem[l].startswith('lghtline'):
                    enertemp = np.sqrt(gdat.limtener[0] * gdat.limtener[1])
                    # temp -- these should depend on population index
                    setp_varb(gdat, 'elin', limt=gdat.limtener, strgmodl=strgmodl)
                    setp_varb(gdat, 'sigm', limt=np.array([1e-1, 1e0]) * enertemp, strgmodl=strgmodl)
                    setp_varb(gdat, 'gamm', limt=np.array([1e-1, 1e0]) * enertemp, strgmodl=strgmodl)
            
            if gdat.boolbinsspat:
                minmdefs = 0.003 / gdat.anglfact
                setp_varb(gdat, 'minmdefs', valu=minmdefs, strgmodl=strgmodl)
        
            if gdat.typeexpr == 'ferm':
                setp_varb(gdat, 'curv', limt=[-1., 1.], strgmodl=strgmodl)
    
            if gdat.boolbinsspat:
                maxmdefs = 1. / gdat.anglfact
                setp_varb(gdat, 'maxmdefs', valu=maxmdefs, strgmodl=strgmodl)
        
            # true model parameters
            if gdat.typedata == 'mock':
                gmod.numbelem = np.zeros(gmod.numbpopl, dtype=int)
                if gmod.typemodltran == 'pois':
                    for l in gmod.indxpopl:
                        setattr(gdat.true.this, 'meanelempop%d' % l, getattr(gdat.true.this, 'numbelempop%d' % l))
                        gmod.numbelem[l] = getattr(gdat.true.this, 'numbelempop%d' % l)
            
                        if gmod.numbelem[l] > gmod.maxmpara.numbelem[l]:
                            raise Exception('True number of elements is larger than maximum.')

            gdat.stdvhostsour = 0.04 / gdat.anglfact
            
            ## distribution
            ### flux
            if gmod.boollenssubh:
                ### projected scale radius
                limtasca = np.array([0., 0.1]) / gdat.anglfact
                setp_varb(gdat, 'asca', minm=minmasca, maxm=maxmasca)
                ### projected cutoff radius
                limtacut = np.array([0., 2.]) / gdat.anglfact
                setp_varb(gdat, 'acut', minm=minmacut, maxm=maxmacut)

            if gdat.boolbinsspat:

                setp_varb(gdat, 'gangdisttype', valu=['self'], strgmodl=strgmodl)
        
                for l in gmod.indxpopl:
                    if gmod.typespatdist[l] == 'gangexpo':
                        setp_varb(gdat, 'maxmgang', valu=gmod.maxmlgal, strgmodl=strgmodl)
                        if gdat.typeexpr == 'ferm':
                            gangdistsexp = 5. / gdat.anglfact
                        setp_varb(gdat, 'gangdistsexp', valu=gangdistsexp, strgmodl=strgmodl, popl=l)
                    if gmod.typespatdist[l] == 'dsrcexpo':
                        if gdat.typeexpr == 'hubb':
                            dsrcdistsexp = 0.5 / gdat.anglfact
                        setp_varb(gdat, 'dsrcdistsexp', valu=dsrcdistsexp, strgmodl=strgmodl, popl=l)
        

                if strgmodl == 'true':
                    if gmod.boollenshost or boolemishost:
                        setp_varb(gdat, 'lgalhost', mean=0., stdv=gdat.stdvhostsour, strgmodl='true', isfr='full')
                        setp_varb(gdat, 'bgalhost', mean=0., stdv=gdat.stdvhostsour, strgmodl='true', isfr='full')
                    if gmod.boollens:
                        setp_varb(gdat, 'lgalsour', mean=0., stdv=gdat.stdvhostsour, strgmodl='true')
                        setp_varb(gdat, 'bgalsour', mean=0., stdv=gdat.stdvhostsour, strgmodl='true')
                if strgmodl == 'fitt':
                    if gmod.boollenshost or boolemishost:
                        setp_varb(gdat, 'lgalhost', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', isfr='full')
                        setp_varb(gdat, 'bgalhost', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', isfr='full')
                    if gmod.boollens:
                        setp_varb(gdat, 'lgalsour', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt')
                        setp_varb(gdat, 'bgalsour', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt')
            
                if gmod.boollens:
                    setp_varb(gdat, 'redshost', limt=[0., 0.4], strgmodl=strgmodl)
                    setp_varb(gdat, 'redssour', limt=[0.5, 1.5], strgmodl=strgmodl)
                    setp_varb(gdat, 'fluxsour', limt=np.array([1e-22, 1e-17]), strgmodl=strgmodl)
                    setp_varb(gdat, 'sindsour', limt=np.array([0., 4.]), strgmodl=strgmodl)
                    setp_varb(gdat, 'sizesour', limt=[0.1 / gdat.anglfact, 2. / gdat.anglfact], strgmodl=strgmodl)
                    setp_varb(gdat, 'ellpsour', limt=[0., 0.5], strgmodl=strgmodl)
                    setp_varb(gdat, 'redshost', valu=0.2, strgmodl=strgmodl)
                    setp_varb(gdat, 'redssour', valu=1., strgmodl=strgmodl)
            
                if gmod.boollenshost or boolemishost:
                    setp_varb(gdat, 'fluxhost', limt=np.array([1e-20, 1e-15]), isfr='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'sindhost', limt=np.array([0., 4.]), isfr='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'sizehost', limt=[0.1 / gdat.anglfact, 4. / gdat.anglfact], isfr='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'beinhost', limt=[0.5 / gdat.anglfact, 2. / gdat.anglfact], isfr='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'ellphost', limt=[0., 0.5], isfr='full', strgmodl=strgmodl)
                    setp_varb(gdat, 'anglhost', limt=[0., np.pi], isfr='full', strgmodl=strgmodl)
                    if strgmodl == 'fitt':
                        setp_varb(gdat, 'serihost', limt=[1., 8.], isfr='full', strgmodl=strgmodl)
                    if strgmodl == 'true':
                        setp_varb(gdat, 'serihost', valu=4., isfr='full', strgmodl=strgmodl)
                        setp_varb(gdat, 'serihost', limt=[1., 8.], isfr='full', strgmodl=strgmodl)
            
                if gmod.boollens:
                    setp_varb(gdat, 'sherextr', limt=[0., 0.1], strgmodl=strgmodl)
                    setp_varb(gdat, 'anglsour', limt=[0., np.pi], strgmodl=strgmodl)
                    setp_varb(gdat, 'sangextr', limt=[0., np.pi], strgmodl=strgmodl)
            
                # temp -- to be removed
                #gmod.factlgal = gmod.maxmlgal - gmod.minmlgal
                #gmod.factbgal = gmod.maxmbgal - gmod.minmbgal
                #gmod.minmaang = -np.pi
                #gmod.maxmaang = pi
            
            # loglikelihood difference for each element
            setp_varb(gdat, 'deltllik', lablroot='$\Delta \log L$', minm=1., maxm=100., strgmodl=strgmodl)
            setp_varb(gdat, 'deltllik', lablroot='$\Delta \log L$', minm=1., maxm=100., popl=l, strgmodl=strgmodl)
            setp_varb(gdat, 'deltllik', lablroot='$\Delta \log L$', minm=1., maxm=100., popl=l, strgmodl=strgmodl, iele='full')
            
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lens':
                    meanslop = 1.9
                    stdvslop = 0.5
                    scal = 'gaus'
                else:
                    minmslop = 0.5
                    maxmslop = 3.
                    scal = 'logt'
                if scal == 'gaus':
                    mean = meanslop
                    stdv = stdvslop
                else:
                    limt = [minmslop, maxmslop]
                
                if gmod.typeelem[l].startswith('clus'):
                    valu = 2.
                name = 'slopprio' + gmod.nameparagenrelemampl[l]

                setp_varb(gdat, name, minm=minmslop, maxm=maxmslop, scal=scal, lablroot='$\alpha$', popl=l, strgmodl=strgmodl)
                
                if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
                    setp_varb(gdat, 'gwdtslop', limt=[0.5, 4.], scal='logt', popl=l, strgmodl=strgmodl)
        
                if gdat.typeexpr != 'user':
                    if gdat.boolbinsspat:
                        setp_varb(gdat, 'spatdistcons', valu=1e-3, popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'gangslop', valu=1.1, popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'bgaldistscal', valu=2. / gdat.anglfact, popl=l, strgmodl=strgmodl)

                if gdat.typeexpr == 'ferm':
                    setp_varb(gdat, 'sloplowrprioflux', valu=1.5, popl=l)
                    setp_varb(gdat, 'slopupprprioflux', valu=2.5, popl=l)
                    setp_varb(gdat, 'brekprioflux', valu=1e-9, popl=l)
                if gmod.typeelem[l] == 'lghtpnts':
                    setp_varb(gdat, 'slopprioflux', valu=2.2, popl=l, strgmodl=strgmodl)
                if gmod.typeelem[l].startswith('lghtline'):
                    setp_varb(gdat, 'slopprioflux', valu=2., popl=l, strgmodl=strgmodl)
                if gmod.typeelem[l] == 'lens':
                    setp_varb(gdat, 'defsslop', valu=1.9, popl=l, strgmodl=strgmodl)

                if gmod.typeelem[l] == 'lens':
                    setp_varb(gdat, 'ascadistmean', valu=0.05 / gdat.anglfact, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'ascadiststdv', valu=0.04 / gdat.anglfact, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'acutdistmean', valu=1. / gdat.anglfact, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'acutdiststdv', valu=0.04 / gdat.anglfact, popl=l, strgmodl=strgmodl)
                
                if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
                    setp_varb(gdat, 'gwdtslop', valu=2., popl=l, strgmodl=strgmodl)
                
                if gdat.typeexpr == 'ferm':
                    sinddistmean = 2.15
                if gdat.typeexpr == 'chan':
                    sinddistmean = 1.
                if gdat.typeexpr == 'hubb':
                    sinddistmean = 1.
                if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'hubb':
                    setp_varb(gdat, 'sinddistmean', valu=sinddistmean, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'sinddiststdv', valu=0.5, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'curvdistmean', valu=2., popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'curvdiststdv', valu=0.2, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'expcdistmean', valu=2., popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'expcdiststdv', valu=0.2, popl=l, strgmodl=strgmodl)
            
                if gmod.typeelem[l] == 'lghtpntspuls':
                    setp_varb(gdat, 'per0distmean', valu=3e-3, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'per0diststdv', valu=0.3, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'magfdistmean', valu=10**8.5, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'magfdiststdv', valu=0.7, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'dglcslop', valu=2., popl=l, strgmodl=strgmodl)
                elif gmod.typeelem[l] == 'lghtpntsagnntrue':
                    setp_varb(gdat, 'dlosslop', valu=-2., popl=l, strgmodl=strgmodl)
                    
                    setp_varb(gdat, 'lum0sloplowr', valu=0.5, popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'lum0slopuppr', valu=1.5, popl=l, strgmodl=strgmodl)
                
            if gmod.boollenshost:
                setp_varb(gdat, 'beinhost', valu=1.5 / gdat.anglfact)
                setp_varb(gdat, 'sizesour', valu=0.3 / gdat.anglfact)
                setp_varb(gdat, 'sizehost', valu=1. / gdat.anglfact)
                setp_varb(gdat, 'ellpsour', valu=0.2)
                setp_varb(gdat, 'fluxsour', valu=1e-18)
                setp_varb(gdat, 'sindsour', valu=1.5)
                setp_varb(gdat, 'fluxhost', valu=1e-16)
                setp_varb(gdat, 'sindhost', valu=2.5)
                setp_varb(gdat, 'ellphost', valu=0.2)
                setp_varb(gdat, 'sangextr', valu=np.pi / 2.)
                setp_varb(gdat, 'serihost', valu=4.)
                
            if gdat.typeexpr != 'user':
                if gdat.boolbinsspat:
                    minm = -gdat.maxmgangdata
                    maxm = gdat.maxmgangdata
                    for l in gmod.indxpopl:
                        setp_varb(gdat, 'lgal', minm=minm, maxm=maxm, lablroot='$l$', strgmodl=strgmodl)
                        setp_varb(gdat, 'bgal', minm=minm, maxm=maxm, lablroot='$b$', strgmodl=strgmodl)
                        setp_varb(gdat, 'lgal', minm=minm, maxm=maxm, lablroot='l_{gal}', popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'bgal', minm=minm, maxm=maxm, lablroot='b_{gal}', popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'lgal', minm=minm, maxm=maxm, lablroot='l_{gal}', popl=l, iele='full', strgmodl=strgmodl)
                        setp_varb(gdat, 'bgal', minm=minm, maxm=maxm, lablroot='b_{gal}', popl=l, iele='full', strgmodl=strgmodl)
                
                minm = 0.1
                maxm = 10.
                for l in gmod.indxpopl:
                    if strgmodl == 'fitt':
                        setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', lablroot='N')
                    setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', lablroot='N', strgmodl=strgmodl)
                    setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', lablroot='N', popl=l, strgmodl=strgmodl)
                    setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', lablroot='N', popl=l, iele='full', strgmodl=strgmodl)
            
                if gdat.boolbinsspat:
                    for l in gmod.indxpopl:
                        setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, lablroot=r'$\theta$', strgmodl=strgmodl)
                        setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, lablroot=r'$\theta$', popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, lablroot=r'$\theta$', popl=l, strgmodl=strgmodl, iele='full')
                        setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, lablroot=r'$\psi$', strgmodl=strgmodl)
                        setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, lablroot=r'$\psi$', popl=l, strgmodl=strgmodl)
                        setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, lablroot=r'$\psi$', popl=l, strgmodl=strgmodl, iele='full')

        
        # copy the true model to the inference model if the inference model parameter has not been specified
        #temp = deepcopy(gdat.__dict__)
        #for strg, valu in temp.items():
        #    if strg.startswith('true') and not strg[4:].startswith('indx'):
        #        try:
        #            valumodl = getattr(gdat.fitt, strg[4:])
        #            if valumodl is None:
        #                raise
        #            if gdat.typeverb > 1:
        #                print 'Received custom input for ' + strg[4:]
        #        except:
        #            setattr(gdat.fitt, strg[4:], getattr(gdat, strg))
        
        # check inputs
        if gdat.numbburn > gdat.numbswep:
            raise Exception('Bad number of burn-in sweeps.')
        if gdat.factthin > gdat.numbswep - gdat.numbburn or gdat.factthin < 1:
            raise Exception('Bad thinning factor.')
        if gdat.typepixl == 'heal' and gdat.numbspatdims > 2:
            raise Exception('More than 2 spatial dimensions require Cartesian binning.')
        
        if gdat.defa:
            return gdat
        
        if gdat.typeverb > 0:
            if gdat.boolburntmpr:
                print('Warning: Tempered burn-in.')
    
        if gdat.typedata == 'inpt':
            gdat.minmpara.sind = -1.
            gdat.maxmpara.sind = 2.
            gdat.minmpara.curv = -1.
            gdat.maxmpara.curv = 1.
            gdat.minmpara.expc = 0.1
            gdat.maxmpara.expc = 10.

            for q in gdat.indxrefr:
                for strgfeat in gdat.refr.namepara.elem[q]:
                    if strgfeat == 'etag' or strgfeat == 'gang' or strgfeat == 'aang':
                        continue
                    refrfeat = getattr(gdat.refr, strgfeat)
                        
                    if len(refrfeat[q]) == 0 or refrfeat[q].ndim < 2:
                        raise Exception('')
        
        if gdat.typedata != 'mock':
            gdat.refr.numbelem = np.zeros(gdat.numbrefr, dtype=int)
        
        for strgmodl in gdat.liststrgmodl:
            
            # set up the indices of the fitting model
            setp_indxpara(gdat, 'finl', strgmodl=strgmodl)
            
            # construct the model
            setp_paragenrscalbase(gdat, strgmodl=strgmodl)
   
            gmod = getattr(gdat, strgmodl)

            if strgmodl == 'true':
                # transfer the true model to the reference model 
                #for strg, valu in gdat.true.__dict__.items():
                #    setattr(gdat.refr, strg, valu)
                for name in ['listmrkrmiss', 'listlablmiss', 'colr', 'colrelem', 'namepara', 'nameparagenrelemampl', 'numbelem']:
                    setattr(gdat.refr, name, getattr(gdat.true, name))
        gdat.refr.indxpoplfittassc = gdat.fitt.indxpopl
        gdat.fitt.indxpoplrefrassc = gdat.fitt.indxpopl

        # to be deleted
        # determine total label
        #for name in ['expo', 'numbpixl']:
        #    lablroot = getattr(gdat.lablrootpara, name)
        #    lablunit = getattr(gdat.lablunitpara, name)
        #    labltotl = tdpy.retr_labltotlsing(lablroot, lablunit)
        #    setattr(gdat.labltotlpara, name, labltotl)
    
        # set the reference model to true model
        # derived lens parameter minima and maxima
        print('Defining minima and maxima for derived parameters...')
        for strgmodl in gdat.liststrgmodl:
            for e in gmod.indxsersfgrd:
                strgsersfgrd = 'isf%d' % e
                setp_varb(gdat, 'masshost' + strgsersfgrd + 'bein', limt=[1e7, 1e14], strgmodl=strgmodl)
                for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                    setp_varb(gdat, 'masshost' + strgsersfgrd + strgcalcmasssubh + 'bein', limt=[1e7, 1e14], strgmodl=strgmodl)
            if gmod.numbparaelem > 0:
                if gmod.boollenssubh:
                    for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                        setp_varb(gdat, 'masssubh' + strgsersfgrd + 'bein', limt=[1e7, 1e10], strgmodl=strgmodl)
                        setp_varb(gdat, 'fracsubh' + strgsersfgrd + 'bein', limt=[0., 1.], strgmodl=strgmodl)
        
        gdat.typeelem = []
        gdat.typeelemspateval = []
        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            for typeelemtemp in gmod.typeelem:
                if not typeelemtemp in gdat.typeelem:
                    gdat.typeelem.append(typeelemtemp)
            for typeelemspatevaltemp in typeelemspateval:
                if not typeelemspatevaltemp in gdat.typeelemspateval:
                    gdat.typeelemspateval.append(typeelemspatevaltemp)
        
        for strgvarb in ['boolelempsfn']:
            varbcomm = False
            for strgmodl in gdat.liststrgmodl:
                gmod = getattr(gdat, strgmodl)
                varb = getattr(gmod, strgvarb)
                varbcomm = varbcomm or varb
            setattr(gdat, strgvarb + 'anyy', varbcomm) 

        #gdat.fitt.namepara.genrelemtagg = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
        #for q in gdat.indxrefr:
        #    for strgfeat in gdat.refr.namepara.elem[q]:
        #        for l in gmod.indxpopl:
        #            gdat.fitt.namepara.genrelemtagg[q][l].append(strgfeat + gdat.listnamerefr[q])
        
        gdat.listnamevarbstat = ['paragenrscalfull', 'paragenrunitfull', 'indxelemfull', 'lliktotl', 'llik', 'lpritotl', 'lpri']
        if gdat.typepixl == 'cart' and (gmod.typeevalpsfn == 'conv' or gmod.typeevalpsfn == 'full'):
            gdat.listnamevarbstat += ['psfnconv']
        if gmod.boolelemsbrtdfncanyy:
            gdat.listnamevarbstat += ['sbrtdfnc']
        if gmod.boolelemsbrtextsbgrdanyy:
            gdat.listnamevarbstat += ['sbrtextsbgrd']
        if gmod.boollens:
            gdat.listnamevarbstat += ['sbrtlens']
        if gmod.boollens or gmod.typeemishost != 'none':
            for e in gmod.indxsersfgrd:
                if gmod.boollens:
                    gdat.listnamevarbstat += ['deflhostisf%d' % e]
                if gmod.typeemishost != 'none':
                    gdat.listnamevarbstat += ['sbrthostisf%d' % e]
        if gmod.convdiffanyy and (gmod.typeevalpsfn == 'full' or gmod.typeevalpsfn == 'conv'):
            gdat.listnamevarbstat += ['sbrtmodlconv']
        if gmod.boolelemdeflsubhanyy:
            gdat.listnamevarbstat += ['deflsubh']
        
        # paths
        ## data
        gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
        gdat.pathprox = gdat.pathdata + 'prox/'
        ## plot
        gdat.pathplotrtag = gdat.pathimag + gdat.rtag + '/'
        gdat.pathinit = gdat.pathplotrtag + 'init/'
        gdat.pathinitintr = gdat.pathinit + 'intr/'
        
        if gdat.boolbinsspat:
            gdat.ascaglob = 0.05 / gdat.anglfact
            gdat.acutglob = 1. / gdat.anglfact
            gdat.cutfdefs = 3e-3 / gdat.anglfact

        # plotting
        gdat.lablsampdist = 'Posterior'
        gdat.lablparagenrscalfull = 'Sample'
        gdat.lablmlik = 'Maximum likelihood'
        gdat.lablmedi = 'Median'
        gdat.lablpmea = 'Mean'
        gdat.lablstdv = 'Std. dev.'
  
        # number of samples for which cumulative posterior will be calculated
        gdat.numbsampcpct = 10
        gdat.indxsampcpct = np.arange(gdat.numbsampcpct)
        
        # p value contours 
        gdat.pvalcont = [0.317, 0.0455, 2.7e-3, 6e-5, 1.3e-6]

        ## number of bins in histogram plots
        gdat.numbbinsplot = 20
        gdat.indxbinsplot = np.arange(gdat.numbbinsplot)
        
        ## number of bins in hyperprior plots
        gdat.numbbinsplotprio = 100
        # temp
        if gdat.typedata == 'inpt':
            for l in gmod.indxpopl:
                for strgpdfn in gmod.listscalparagenrelem[l]:
                    if strgpdfn.startswith('gaum') and gmod.lgalprio is None and gmod.bgalprio is None:
                        raise Exception('If typespatdist is "gaus", spatial coordinates of the prior catalog should be provided via lgalprio and bgalprio.')
        
        # temp -- have these definitions separate for all features
        # feature plotting factors and scalings
        gdat.dictglob = {}
        
        gdat.listnamechro = ['totl', 'prop', 'diag', 'save', 'plot', 'proc', 'pars', 'modl', 'llik', 'sbrtmodl']
        gdat.listlablchro = ['Total', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Process', 'Parse', 'Model', 'Likelihood', 'Total emission']
        if gmod.numbparaelem > 0:
            gdat.listnamechro += ['spec']
            gdat.listlablchro += ['Spectrum calculation']
        if gmod.boollens:
            gdat.listnamechro += ['deflzero', 'deflhost', 'deflextr', 'sbrtlens', 'sbrthost']
            gdat.listlablchro += ['Array initialization', 'Host Deflection', 'External deflection', 'Lensed emission', 'Host emission']
        if gmod.boolelemsbrtdfncanyy:
            gdat.listnamechro += ['elemsbrtdfnc']
            gdat.listlablchro += ['Dfnc S Brght']
        if gmod.boolelemdeflsubhanyy:
            gdat.listnamechro += ['elemdeflsubh']
            gdat.listlablchro += ['Subh Defl']
        if gmod.boolelemsbrtextsbgrdanyy:
            gdat.listnamechro += ['elemsbrtextsbgrd']
            gdat.listlablchro += ['Bkg Exts S Brght']
        booltemp = False
        for strgmodl in gdat.liststrgmodl:
            booltemp = booltemp or gmod.typeevalpsfn
        if booltemp or gmod.typeevalpsfn == 'full' or gmod.typeevalpsfn == 'full':
            gdat.listnamechro += ['psfnconv']
            gdat.listlablchro += ['Img for PSF Conv.']
        
        gdat.listnamechro += ['expo', 'lpri', 'tert']
        gdat.listlablchro += ['Exposure', 'Prior', 'Tertiary']
        gdat.numbchro = len(gdat.listnamechro)
        
        if gdat.typedata != 'mock':
            if gmod.boolelemlghtanyy and gdat.typeexpr == 'ferm' and gdat.maxmgangdata == 20. / gdat.anglfact:
                path = gdat.pathinpt + 'sbrt0018.png'
                gdat.sbrt0018 = sp.ndimage.imread(path, flatten=True)
                gdat.sbrt0018 -= np.amin(gdat.sbrt0018)
                gdat.sbrt0018 /= np.amax(gdat.sbrt0018)
                binslgaltemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[1])
                binsbgaltemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[0])
                gdat.sbrt0018objt = sp.interpolate.RectBivariateSpline(binsbgaltemp, binslgaltemp, gdat.sbrt0018)

        # log-prior register
        ## indices of split and merge term
        indxlprispme = -1
        ## number of elements
        numb = 0
        for l in gmod.indxpopl:
            numb += len(gmod.namepara.genrelem[l])
        
        # process index
        gdat.indxproc = np.arange(gdat.numbproc)

        if gmod.boollens or gdat.typedata == 'mock' and gmod.boollens:
            retr_axis(gdat, 'mcut')
            retr_axis(gdat, 'bein')

            # angular deviation
            gdat.numbanglhalf = 10
            gdat.indxanglhalf = np.arange(gdat.numbanglhalf)
            retr_axis(gdat, 'anglhalf')
            gdat.numbanglfull = 1000
            gdat.indxanglfull = np.arange(gdat.numbanglfull)
            gdat.minmpara.anglfull = 0.
            gdat.maxmpara.anglfull = 3. * gdat.maxmgangdata
            retr_axis(gdat, 'anglfull')
        
        # temp
        #gdat.binspara.anglcosi = np.sort(np.cos(gdat.binspara.angl))
        
        # temp
        #gdat.meshbackener = np.meshgrid(gdat.gmod.indxback, gdat.indxener, indexing='ij')
        
        # plotting
        ## the normalized offset for text annotation of point sources in the frames
        gdat.offstextimag = gdat.maxmgangdata * 0.05
        
        ## figure size
        gdat.plotsize = 6
        ## size of the images
        gdat.sizeimag = 1.3 * gdat.plotsize
        
        ## label of the models
        gdat.fitt.lablmodl = 'Model'
        if gdat.typedata == 'mock':
            gdat.refr.lablmodl = 'True'
        else:
            gdat.refr.lablmodl = 'Ref'
        
        # element parameters common between the fitting and reference models
        gdat.namepara.elemcomm = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
        for q in gdat.indxrefr:
            for l in gmod.indxpopl:
                for strgfeat in gmod.listnameparatotlelem[l]:
                    if strgfeat in gdat.refr.namepara.elem[q]:
                        gdat.namepara.elemcomm[q][l].append(strgfeat)
        
        if gdat.typedata == 'mock':
            gdat.refr.indxpopl = gdat.true.indxpopl
            gdat.refr.lablpopl = gdat.true.lablpopl

        
        for strgmodl in ['refr', 'fitt']:
        
            gmod = getattr(gdat, strgmodl)
            
            print('strgmodl')
            print(strgmodl)
            # labels of elements
            lablelem = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                lablelem[l] = gmod.lablmodl + ' ' + gmod.lablpopl[l]
            setp_varb(gdat, 'lablelem', valu=lablelem, strgmodl=strgmodl)
            
            lablelemmiss = [[] for l in gmod.indxpopl]
            lablelemhits = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                lablelemmiss[l] = gmod.lablelem[l] + ' miss'
                lablelemhits[l] = gmod.lablelem[l] + ' hit'
            setp_varb(gdat, 'lablelemmiss', valu=lablelemmiss, strgmodl=strgmodl)
            setp_varb(gdat, 'lablelemhits', valu=lablelemhits, strgmodl=strgmodl)
            
            lablhost = gmod.lablmodl + ' host'
            setp_varb(gdat, 'lablhost', valu=lablhost, strgmodl=strgmodl)
            
            lablsour = gmod.lablmodl + ' sour'
            setp_varb(gdat, 'lablsour', valu=lablsour, strgmodl=strgmodl)

        ## PSF class indices for which images will be plotted
        if gdat.numbevtt == 1:
            gdat.indxevttplot = gdat.indxevtt
        else:
            gdat.indxevttplot = np.concatenate((np.array([-1]), gdat.indxevtt))
        
        gdat.numbenerevtt = gdat.numbener * gdat.numbevtt
        
        # temp
        gdat.boolintpanglcosi = False

        if gdat.boolthindata:
            gdat.factdatathin = 10
            if gdat.typepixl != 'cart' or gdat.numbsidecart % gdat.factdatathin != 0:
                raise Exception('Cannot thin the data.')
            #gdat.indxpixlkeep = gdat.indxpixlfull[::gdat.factdatathin]
            #gdat.numbpixlkeep = gdat.indxpixlkeep.size
            gdat.indxpixlkill = np.setdiff1d(gdat.indxpixlfull, gdat.indxpixlkeep)
            gdat.numbsidecart = gdat.numbsidecart / 10
            gdat.numbsidecarthalf = int(gdat.numbsidecart / 2)
            gdat.lgalgrid = gdat.lgalgrid[gdat.indxpixlkeep]
            gdat.bgalgrid = gdat.bgalgrid[gdat.indxpixlkeep]
            gdat.indxpixlfull = gdat.indxpixlfull[gdat.indxpixlkeep]
            
        # the function to measure time
        # temp
        gdat.strgfunctime = 'clck'
        if gdat.strgfunctime == 'clck':
            gdat.functime = time.clock
        if gdat.strgfunctime == 'time':
            gdat.functime = time.time

        ## longitude
        gdat.numblgalpntsprob = gdat.numbsidepntsprob
        gdat.numbbgalpntsprob = gdat.numbsidepntsprob
        gdat.binspara.lgalpntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
        gdat.binspara.bgalpntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
        gdat.indxlgalpntsprob = np.arange(gdat.numblgalpntsprob)
        gdat.indxbgalpntsprob = np.arange(gdat.numbbgalpntsprob)

        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            if gmod.boollens or gdat.typedata == 'mock' and gmod.boollens:
                retr_axis(gdat, 'defl')
                retr_axis(gdat, 'deflsubh')

        # lensing problem setup
        ## number of deflection components to plot

        gdat.binspara.lgalcartmesh, gdat.binspara.bgalcartmesh = np.meshgrid(gdat.binspara.lgalcart, gdat.binspara.bgalcart, indexing='ij')
        gdat.meanpara.lgalcartmesh, gdat.meanpara.bgalcartmesh = np.meshgrid(gdat.meanpara.lgalcart, gdat.meanpara.bgalcart, indexing='ij')
        if gdat.typepixl == 'cart':
            gdat.sizepixl = np.sqrt(gdat.apix)
            gdat.indxsidecart = np.arange(gdat.numbsidecart)
            gdat.indxpixlrofi = np.arange(gdat.numbpixlcart)
            gdat.indxsidemesh = np.meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
            gdat.lgalgrid = gdat.meanpara.lgalcart[gdat.indxsidemesh[0].flatten()]
            gdat.bgalgrid = gdat.meanpara.bgalcart[gdat.indxsidemesh[1].flatten()]
            gdat.shapcart = (gdat.numbsidecart, gdat.numbsidecart)
            gdat.lgalgridfull = np.copy(gdat.lgalgrid)
            gdat.bgalgridfull = np.copy(gdat.bgalgrid)
            gdat.lgalgridcart = gdat.lgalgrid.reshape(gdat.shapcart)
            gdat.bgalgridcart = gdat.bgalgrid.reshape(gdat.shapcart)
            gdat.indxpent = np.meshgrid(gdat.indxener, gdat.indxsidecart, gdat.indxsidecart, gdat.indxevtt, indexing='ij')
        if gdat.typepixl == 'heal':
            lgalheal, bgalheal, gdat.numbpixlfull, gdat.apix = tdpy.retr_healgrid(gdat.numbsideheal)
            lgalheal = np.deg2rad(lgalheal)
            bgalheal = np.deg2rad(bgalheal)
   
            gdat.indxpixlrofi = np.where((np.fabs(lgalheal) < gdat.maxmgangdata) & (np.fabs(bgalheal) < gdat.maxmgangdata))[0]
            
            gdat.indxpixlrofimarg = np.where((np.fabs(lgalheal) < 1.2 * gdat.maxmgangdata) & (np.fabs(bgalheal) < 1.2 * gdat.maxmgangdata))[0]

            gdat.lgalgrid = lgalheal
            gdat.bgalgrid = bgalheal
        
        gdat.indxpixlfull = np.arange(gdat.numbpixlfull)
        if gdat.typepixl == 'cart':
            gdat.indxpixlcart = np.arange(gdat.numbpixlcart)
        
        if gdat.evttbins:
            # PSF class string
            gdat.strgevtt = []
            for m in gdat.indxevtt:
                gdat.strgevtt.append('PSF%d' % gdat.indxevttincl[m])
        
        # power spectra
        if gdat.typepixl == 'cart':
            setp_varb(gdat, 'anglodim', minm=0., maxm=1., boolinvr=True)
            setp_varb(gdat, 'mpolodim', minm=0., maxm=1.)
            #retr_axis(gdat, 'anglodim', boolinvr=True)
            #retr_axis(gdat, 'mpolodim')
                
            for strgmodl in gdat.liststrgmodl:
                gmod = getattr(gdat, strgmodl)
                
                gdat.numbwvecodim = gdat.numbsidecart
                gdat.minmanglodim = 0.
                gdat.maxmanglodim = 2. * gdat.maxmgangdata
                gdat.minmmpolodim = 0.
                gdat.maxmmpolodim = 1. / 2. / gdat.sizepixl

                if gmod.boollens or gdat.typedata == 'mock' and gmod.boollens:
                    # temp -- this should minima, maxima of adishost and the true metamodel into account
                    gdat.minmwvecodim = gdat.minmmpolodim / np.amax(gmod.adishost)
                    gdat.maxmwvecodim = gdat.maxmmpolodim / np.amin(gmod.adishost)
                    gdat.minmwlenodim = gdat.minmanglodim * np.amin(gmod.adishost)
                    gdat.maxmwlenodim = gdat.maxmanglodim * np.amax(gmod.adishost)
                    retr_axis(gdat, 'wvecodim', strgmodl=strgmodl)
                    retr_axis(gdat, 'wlenodim', strgmodl=strgmodl, boolinvr=True)
                    gdat.meanpara.wveclgal, gdat.meanpara.wvecbgal = np.meshgrid(gdat.meanpara.wvecodim, gdat.meanpara.wvecodim, indexing='ij')
                    gdat.meanpara.wvec = np.sqrt(gdat.meanpara.wveclgal**2 + gdat.meanpara.wvecbgal**2)
            gdat.meanpara.mpollgal, gdat.meanpara.mpolbgal = np.meshgrid(gdat.meanpara.mpolodim, gdat.meanpara.mpolodim, indexing='ij')
            gdat.meanpara.mpol = np.sqrt(gdat.meanpara.mpollgal**2 + gdat.meanpara.mpolbgal**2)

        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            
            # element parameter vector indices
            gmod.indxparagenrelemlgal = 0
            gmod.indxparagenrelembgal = 1
            gmod.indxparagenrelemflux = 2
            gmod.indxparagenrelemsind = 3
            gmod.indxparagenrelemcurv = 4
            gmod.indxparagenrelemexpc = 4

        # check the exposure map data structure
        if gdat.boolcorrexpo:
            booltemp = False
            if gdat.expo.ndim != 3:
                booltemp = True
            if gdat.typepixl == 'cart' and gdat.expo.shape[1] != gdat.numbpixlcart:
                booltemp = True
            if booltemp:
                raise Exception('Exposure does not have the right data structure. It should be a list of 3D np.arrays.')
            
            if gdat.boolsqzeexpo:
                gdat.expo *= 1e-10
            if gdat.boolexplexpo:
                gdat.expo *= 1e10
        
        if gdat.boolthindata:
            #gdat.expo[:, gdat.indxpixlkill, :] = 0.
            expotemp = np.copy(gdat.expo[:, gdat.indxpixlfull[::gdat.factdatathin], :])
            sbrttemp = np.copy(gdat.sbrtdata[:, gdat.indxpixlfull[::gdat.factdatathin], :])
            gdat.expo = expotemp 
            gdat.sbrtdata = sbrttemp
            
        # only include desired energy and PSF class bins
        gdat.indxcubeincl = np.meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
        
        ## exposure
        if gdat.boolcorrexpo:
            # temp -- for some reason lists of np.arrays require manual processing
            gdat.expo = gdat.expo[tuple(gdat.indxcubeincl)]
            if gdat.typedata == 'inpt':
                gdat.sbrtdata = gdat.sbrtdata[tuple(gdat.indxcubeincl)]
        
        ## backgrounds
        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            gmod.sbrtbacknormincl = [[] for c in gmod.indxback]
            for c in gmod.indxback:
                gmod.sbrtbacknormincl[c] = gmod.sbrtbacknorm[c][tuple(gdat.indxcubeincl)]
        
        # obtain cartesian versions of the maps
        #if gdat.typepixl == 'cart':
        #    gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
        #    for strgmodl in gdat.liststrgmodl:
        #        gmod.sbrtbacknormcart = []
        #        for c in getattr(gmod, 'gmod.indxback'):
        #            gmod.sbrtbacknormcart.append(gmod.sbrtbacknorm[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
        
        # mask the exposure map
        if gdat.listmask is not None:
            for mask in gdat.listmask:
                if mask[0] == 'sqre':
                    indxpixlmask = np.where((gdat.lgalgrid > mask[1]) & (gdat.lgalgrid < mask[2]) & (gdat.bgalgrid > mask[3]) & (gdat.bgalgrid < mask[4]))[0]
                if mask[0] == 'circ':
                    indxpixlmask = np.where(np.sqrt((gdat.lgalgrid - mask[1])**2 + (gdat.bgalgrid - mask[2])**2) < mask[3])[0]
                if mask[0] == 'hstr':
                    indxpixlmask = np.where((gdat.bgalgrid > mask[1]) & (gdat.bgalgrid < mask[2]))[0]
                if gdat.typemaskexpo == 'zero':
                    gdat.expo[:, indxpixlmask, :] = 0.
                if gdat.typemaskexpo == 'ignr':
                    gdat.expo[:, indxpixlmask, :] = 1e-49

        # plotting
        ## ROI
        if gdat.boolbinsspat:
            gdat.exttrofi = np.array([gdat.minmlgaldata, gdat.maxmlgaldata, gdat.minmbgaldata, gdat.maxmbgaldata])
            gdat.exttrofi *= gdat.anglfact 
            gdat.frambndrdata = gdat.maxmgangdata * gdat.anglfact

        ## marker size
        gdat.minmmrkrsize = 100
        gdat.maxmmrkrsize = 500
        ## marker line width
        gdat.mrkrlinewdth = 3
        ## marker opacity
        gdat.alphhist = 0.5
        gdat.alphline = 0.5
        gdat.alphbndr = 0.5
        gdat.alphelem = 1.
        gdat.alphmaps = 1.
        
        # number of colorbar ticks in the maps
        gdat.numbtickcbar = 11
        
        ## color bars
        gdat.minmlpdfspatpriointp = np.log(1. / 2. / gdat.maxmgangdata) - 10.
        gdat.maxmlpdfspatpriointp = np.log(1. / 2. / gdat.maxmgangdata) + 10.
        gmod.scallpdfspatpriointp = 'self'
        gdat.cmaplpdfspatpriointp = 'PuBu'
        
        gdat.minmllikmaps = -10.
        gdat.maxmllikmaps = 0.
        gmod.scalllikmaps = 'asnh'
        gdat.cmapllikmaps = 'YlGn'
        
        gdat.minmperc = 0.
        gdat.maxmperc = 1e2
        gdat.scalperc = 'asnh'
        gdat.cmapperc = 'afmhot'
        
        gdat.minmpercresi = -1e2
        gdat.maxmpercresi = 1e2
        gdat.scalpercresi = 'asnh'
        gdat.cmappercresi = 'coolwarm'
        
        gdat.scalpara.cntpdata = 'logt'
        gdat.cmappara.cntpdata = 'Greys'
        
        gdat.scalpara.cntpmodl = 'logt'
        gdat.cmappara.cntpmodl = 'Greys'
        
        gdat.scalpara.cntpresi = 'asnh'
        gdat.cmappara.cntpresi = make_cmapdivg('Red', 'Orange')

        gdat.minmconv = 1e-2
        gdat.maxmconv = 10.
        gdat.scalconv = 'logt'
        gdat.cmapconv = 'Purples'
        
        gdat.minmconvelem = 1e-4
        gdat.maxmconvelem = 1e-1
        gdat.scalconvelem = 'logt'
        gdat.cmapconvelem = 'Purples'
        
        gdat.minms2nr = 0.
        gdat.maxms2nr = 10.
        gmod.scals2nr = 'asnh'
        gdat.cmaps2nr = 'magma'
        
        gdat.minmmagn = -1e2
        gdat.maxmmagn = 1e2
        gmod.scalmagn = 'asnh'
        gdat.cmapmagn = 'BrBG'
        
        gdat.minmdeflresiperc = -100.
        gdat.maxmdeflresiperc = 100.
        gmod.scaldeflresiperc = 'self'
        gdat.cmapdeflresiperc = 'Oranges'
        
        gdat.minmconvelemresi = -0.1
        gdat.maxmconvelemresi = 0.1
        gmod.scalconvelemresi = 'self'
        gdat.cmapconvelemresi = 'PiYG'
        
        gdat.minmconvelemresiperc = -100.
        gdat.maxmconvelemresiperc = 100.
        gmod.scalconvelemresiperc = 'self'
        gdat.cmapconvelemresiperc = 'PiYG'
        
        gdat.minmmagnresi = -10.
        gdat.maxmmagnresi = 10.
        gmod.scalmagnresi = 'self'
        gdat.cmapmagnresi = 'PRGn'
        
        gdat.minmmagnresiperc = -100.
        gdat.maxmmagnresiperc = 100.
        gmod.scalmagnresiperc = 'self'
        gdat.cmapmagnresiperc = 'PRGn'
        
        gdat.lgalgrid = gdat.lgalgrid[gdat.indxpixlrofi]
        gdat.bgalgrid = gdat.bgalgrid[gdat.indxpixlrofi]
   
        if gdat.boolcorrexpo:
            if np.amax(gdat.expo) <= 0.:
                raise Exception('Bad exposure.')

        # temp
        #gdat.expo[np.where(gdat.expo < 1e-50)] = 1e-50
        
        # exclude voxels with vanishing exposure
        if gdat.boolcorrexpo:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    gdat.indxpixlrofi = np.intersect1d(gdat.indxpixlrofi, np.where(gdat.expo[i, :, m] > 0.)[0])
        
        gdat.indxcuberofi = np.meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
        gdat.numbpixl = gdat.indxpixlrofi.size
        gdat.indxpixl = np.arange(gdat.numbpixl)
        gdat.numbdata = gdat.numbener * gdat.numbevtt * gdat.numbpixl

        #gdat.lgalgridrofi = gdat.lgalgrid[gdat.indxpixlrofi]
        #gdat.bgalgridrofi = gdat.bgalgrid[gdat.indxpixlrofi]


        if gdat.typedata == 'inpt':
            gdat.sbrtdata = gdat.sbrtdata[tuple(gdat.indxcuberofi)]

        ## exposure
        if gdat.boolcorrexpo:
            gdat.expofull = np.copy(gdat.expo)
            gdat.expo = gdat.expo[tuple(gdat.indxcuberofi)]
        
            gdat.minmpara.expo = np.amin(gdat.expo[np.where(gdat.expo > 1e-100)])
            gdat.maxmpara.expo = np.amax(gdat.expo)
            gdat.minmpara.expo = np.amin(gdat.minmpara.expo)
            gdat.maxmpara.expo = np.amax(gdat.maxmpara.expo)
        
        # required to convert to an index of non-zero exposure pixels
        #if gdat.minmpara.expo > 0:
        #    gdat.indxpixlroficnvt = np.arange(gdat.numbpixlfull)
        #else:
        #    cntr = 0
        #    gdat.indxpixlroficnvt = full(gdat.numbpixlfull, -1)
        #    for j in gdat.indxpixlfull:
        #        if j in gdat.indxpixlrofi:
        #            gdat.indxpixlroficnvt[j] = cntr
        #            cntr += 1
        #
        
        ## backgrounds
        for strgmodl in gdat.liststrgmodl:
            if gdat.typepixl == 'heal':
                sbrtbackhealfull = [[] for c in gmod.indxback]
                for c in gmod.indxback:
                    sbrtbackhealfull[c] = np.copy(gmod.sbrtbacknorm[c])
            gmod.sbrtbacknormincl = [[] for c in gmod.indxback]
            for c in gmod.indxback:
                gmod.sbrtbacknormincl[c] = gmod.sbrtbacknorm[c][tuple(gdat.indxcuberofi)]
        
        if gdat.boolcorrexpo:
            gdat.expototl = []
            gdat.expototlmean = []
            gdat.expototl = np.sum(gdat.expo, axis=2)
            gdat.expototlmean = np.mean(gdat.expototl, axis=1)

        if gdat.typeelemspateval == 'locl':
            if gdat.typeexpr == 'gene':
                gdat.maxmangl = 1.
            if gdat.typeexpr == 'ferm':
                gdat.maxmangl = 20. / gdat.anglfact
            if gdat.typeexpr == 'tess':
                gdat.maxmangl = 25. / gdat.anglfact
            if gdat.typeexpr == 'chan':
                gdat.maxmangl = 15. / gdat.anglfact
            if gdat.typeexpr == 'hubb':
                gdat.maxmangl = 1. / gdat.anglfact
        else:
            gdat.maxmangl = gdat.maxmgangdata * np.sqrt(2.) * 2. * 1.1
            
        gdat.listnamespatmean = ['full']
        if gdat.typeexpr == 'ferm':
            gdat.listnamespatmean += ['innr']
        gdat.numbspatmean = len(gdat.listnamespatmean)
        gdat.indxspatmean = np.arange(gdat.numbspatmean)
        gdat.listindxcubespatmean = [[] for b in gdat.indxspatmean]
        gdat.indxcube = np.meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        for b, namespatmean in enumerate(gdat.listnamespatmean):
            if namespatmean == 'full':
                gdat.listindxcubespatmean[b] = gdat.indxcube
            if namespatmean == 'innr':
                gdat.indxpixlinnr = np.where(np.sqrt(gdat.lgalgrid**2 + gdat.bgalgrid**2) < 5. / gdat.anglfact)[0]
                gdat.listindxcubespatmean[b] = np.meshgrid(gdat.indxener, gdat.indxpixlinnr, gdat.indxevtt, indexing='ij')
        
        if gdat.numbpixl > 1:
            # store pixels as unit vectors
            gdat.xdatgrid, gdat.ydatgrid, gdat.zaxigrid = retr_unit(gdat.lgalgrid, gdat.bgalgrid)
   
            # construct a lookup table for converting HealPix pixels to ROI pixels
            if gdat.typepixl == 'heal':
                path = gdat.pathpixlcnvt + 'pixlcnvt_%09g.p' % gdat.maxmgangdata

                if os.path.isfile(path):
                    fobj = open(path, 'rb')
                    gdat.pixlcnvt = pickle.load(fobj)
                    fobj.close()
                else:
                    gdat.pixlcnvt = np.zeros(gdat.numbpixlfull, dtype=int) - 1
                    numbpixlmarg = gdat.indxpixlrofimarg.size
                    for k in range(numbpixlmarg):
                        dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                        gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
                    fobj = open(path, 'wb')
                    pickle.dump(gdat.pixlcnvt, fobj, protocol=pickle.HIGHEST_PROTOCOL)
                    fobj.close()
            
            # dummy pixel indices for full (nonlocal) element kernel evaluation 
            gdat.listindxpixl = []
            if gdat.typedata == 'mock':
                numb = max(np.sum(gmod.maxmpara.numbelem), np.sum(gmod.maxmpara.numbelem)) + 2
            else:
                numb = np.sum(gmod.maxmpara.numbelem) + 2
            for k in range(int(numb)):
                gdat.listindxpixl.append([])
                for kk in range(k):
                    gdat.listindxpixl[k].append(gdat.indxpixl)
            
            # spatial averaging setup
            # temp
            
            # temp -- check if 1000 is too much
            gdat.numbanglelem = 1000
        
        # turn off relevant proposal types
        gdat.numbprop = 5
        gdat.indxprop = np.arange(gdat.numbprop)
        
        gdat.numbstdp = gmod.numbparagenrbase - gmod.numbpopl
        cntr = 0
        for l in gmod.indxpopl:
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                setattr(gdat.fitt.inxparagenrscalelemkind, nameparagenrelem + 'pop%d' % l, gdat.numbstdp)
                cntr += 1
        gdat.numbstdp += cntr
        
        gdat.lablstdp = np.copy(np.array(gmod.labltotlpara.genrbase[gmod.numbpopl:]))
        gdat.namestdp = np.copy(np.array(gmod.nameparagenrbase[gmod.numbpopl:]))
        for l in gmod.indxpopl:
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                gdat.lablstdp = np.append(gdat.lablstdp, getattr(gdat.fitt.labltotlpara, nameparagenrelem))
                gdat.namestdp = np.append(gdat.namestdp, nameparagenrelem + 'pop%d' % l)
        gdat.namestdp = gdat.namestdp.astype(object)
        gdat.lablstdp = list(gdat.lablstdp)
        gdat.indxstdp = np.arange(gdat.numbstdp)
        gdat.indxstdpprop = gdat.indxstdp
        
        # proposal scale indices for each parameter
        indxelemfull = [list(range(gmod.maxmpara.numbelem[l])) for l in gmod.indxpopl]
        gdat.fitt.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, indxelemfull, 'fitt')
        
        gdat.indxstdppara = np.zeros(gmod.numbparagenrfull, dtype=int) - 1
        cntr = 0
        gdat.indxstdppara[gmod.numbpopl:gmod.numbparagenrbase] = gmod.indxparagenrbase[gmod.numbpopl:] - gmod.numbpopl
        if gmod.numbparaelem > 0:
            for l in gmod.indxpopl:
                for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    for indx in gdat.fitt.this.indxparagenrfullelem[l][nameparagenrelem]:
                        gdat.indxstdppara[indx] = cntr + gmod.numbparagenrbase - gmod.numbpopl
                    cntr += 1
        
        # for the fitting model, define proposal type indices
        for name, valu in gmod.indxpara.__dict__.items():
            if not name.startswith('numbelem') and name != 'dist':
                if not isinstance(valu, int):
                    continue
                indxstdp = gdat.indxstdppara[valu]
                setattr(gdat, 'indxstdp' + name, indxstdp)
        
        # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
        gmod.corr = tdpy.gdatstrt()
        for k in gmod.indxvarbscal:
            name = gmod.namepara.scal[k]
            try:
                temp = getattr(gdat.true.this, name)
            except:
                temp = None
            setattr(gmod.corr, name, temp)

        gmod.corrparagenrscalbase = np.empty(gmod.numbparagenrbase)
        for k in gmod.indxparagenrbase:
            try:
                gmod.corrparagenrscalbase[k] = getattr(gdat.true, gmod.nameparagenrbase[k])
            except:
                gmod.corrparagenrscalbase[k] = None

        for namepara in gdat.fitt.listnameparaglob:
            setattr(gdat.labltotlpara, namepara, getattr(gdat.fitt.labltotlpara, namepara))

        # set parameter features common between true and fitting models
        for strgmodl in gdat.liststrgmodl:
    
            gmod = getattr(gdat, strgmodl)
            
            for namepara in gmod.namepara.kind:
                try:
                    getattr(gdat.minmpara, namepara)
                    getattr(gdat.maxmpara, namepara)
                except:
                    try:
                        setattr(gdat.minmpara, namepara, min(getattr(gdat.fitt.minmpara, namepara), getattr(gdat.true.minmpara, namepara)))
                        setattr(gdat.maxmpara, namepara, max(getattr(gdat.fitt.maxmpara, namepara), getattr(gdat.true.maxmpara, namepara)))
                    except:
                        try:
                            setattr(gdat.minmpara, namepara, getattr(gdat.fitt.minmpara, namepara))
                            setattr(gdat.maxmpara, namepara, getattr(gdat.fitt.maxmpara, namepara))
                        except:
                            setattr(gdat.minmpara, namepara, getattr(gdat.true.minmpara, namepara))
                            setattr(gdat.maxmpara, namepara, getattr(gdat.true.minmpara, namepara))
        
        # set plot limits for each model if not already set (for Gaussian, log-normal distributions)
        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            
            for namepara in gmod.namepara.kind:
                minm = getattr(gmod.minmpara, namepara)
                maxm = getattr(gmod.maxmpara, namepara)
                limt = np.array([minm, maxm])
                setattr(gmod.limtpara, namepara, limt)
        
        # construct bins for scalar parameters
        for namevarbscal in gmod.namepara.scal:
            
            # variables with only label and scaling
            if namevarbscal == 'lliktotl' or namevarbscal == 'lpripena':
                continue

            print('temp -- place here setp_varb for all variables')
            #retr_axis(gdat, namevarbscal)
        
        gmod = gdat.fitt
        # proposal scale
        if gmod.boollens or gdat.typedata == 'mock':
            
            gdat.stdp = 1e-4 + np.zeros(gdat.numbstdp)
            
            if gmod.typemodltran == 'pois' and gmod.numbpopl > 0:
                if gmod.maxmpara.numbelem[0] > 0:
                    gdat.stdp[gdat.indxstdpmeanelempop0] = 1e-1
            
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.sigcen00evt0]] = 3e-2
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.bacpback0000en00]] = 1e-3
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.bacpback0000en00]] = 1e-1
            
            if gmod.boollens:
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.lgalsour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.bgalsour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.fluxsour]] = 1e-2
                if gdat.numbener > 1:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.sindsour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.sizesour]] = 1e-1
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.ellpsour]] = 1e-1
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.anglsour]] = 1e-1
            if gmod.typeemishost != 'none':
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.lgalhostisf0]] = 3e-4
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.bgalhostisf0]] = 3e-4
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.fluxhostisf0]] = 1e-3
                if gdat.numbener > 1:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.sindhostisf0]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.sizehostisf0]] = 3e-3
            
            if gmod.boollens:
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.beinhostisf0]] = 1e-3
            if gmod.typeemishost != 'none':
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.ellphostisf0]] = 1e-2
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.anglhostisf0]] = 1e-2
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.serihostisf0]] = 1e-2
            if gmod.boollens:
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.sherextr]] = 1e-1
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.sangextr]] = 3e-2
            
        else:
            
            if gdat.typeexpr == 'ferm':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
                
                if gmod.typemodltran == 'pois' and gmod.numbparaelem > 0:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.meanelem]] = 4e-2
                    
                    for l in gmod.indxpopl:
                        if gmod.typeprioflux[l] == 'powr':
                            gdat.stdp[gdat.indxstdppara[gmod.indxpara.sloppriofluxpop0]] = 1e-1
                        else:
                            gdat.stdp[gdat.indxstdppara[gmod.indxpara.brekpriofluxpop0]] = 1e-1
                            gdat.stdp[gdat.indxstdppara[gmod.indxpara.sloplowrpriofluxpop0]] = 1e-1
                            gdat.stdp[gdat.indxstdppara[gmod.indxpara.slopupprpriofluxpop0]] = 1e-1
                
                gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 5e-3
                gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en01')]] = 1e-2
                gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en02')]] = 3e-2
                
                if 'fdfm' in gmod.listnameback: 
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0001en00')]] = 8e-4
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0001en01')]] = 1e-3
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0001en02')]] = 2e-3
                
                if 'dark' in gmod.listnameback: 
                    gmod.indxbackdark = gmod.listnameback.index('dark')
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback%04den00' % gmod.indxbackdark)]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback%04den01' % gmod.indxbackdark)]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback%04den02' % gmod.indxbackdark)]] = 3e-2
                
                if gmod.numbparaelem > 0:
                    gdat.stdp[gdat.indxstdppop0flux] = 8e-2

                    if gmod.spectype[0] == 'colr':
                        gdat.stdp[gdat.indxstdppop0sindcolr0001] = 8e-2
                        gdat.stdp[gdat.indxstdppop0sindcolr0002] = 2e-1
            
            if gdat.typeexpr == 'chan':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
                
                if gmod.typemodltran == 'pois' and gmod.numbparaelem > 0:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.meanelem]] = 2e-1
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.sloppriofluxpop0]] = 2e-1
                
                if gmod.numbparaelem > 0 and gdat.boolbinsspat:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.psfp]] = 4e-1
                
                if gdat.indxenerincl.size == 5 and (gdat.indxenerincl == np.arange(5)).all():
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en01')]] = 3e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en02')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en03')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en04')]] = 1e-2
                elif gdat.indxenerincl.size == 2 and (gdat.indxenerincl == np.array([2])).all():
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
                
                if gmod.numbparaelem > 0:
                    if gdat.boolbinsspat:
                        gdat.stdp[gdat.fitt.inxparagenrscalelemkind.lgalpop0] = 2e-2
                        gdat.stdp[gdat.fitt.inxparagenrscalelemkind.bgalpop0] = 2e-2
                        if gdat.numbener > 1:
                            gdat.stdp[gdat.indxstdppop0sind] = 2e-1
                    gdat.stdp[gdat.indxstdppop0flux] = 2e-1
            
            if gdat.typeexpr == 'gene':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
            
                if gmod.typemodltran == 'pois' and gmod.numbparaelem > 0:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.meanelem]] = 2e-1
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.slopprionobjpop0]] = 3e-1
                    try:
                        gdat.stdp[gdat.indxstdppara[gmod.indxpara.gwdtsloppop0]] = 3e-1
                    except:
                        pass

                if gmod.typeevalpsfn != 'none' and gdat.boolmodipsfn:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.psfp]] = 4e-1
                
                gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
            
                if gmod.numbparaelem > 0:
                    gdat.stdp[gdat.fitt.inxparagenrscalelemkind.lgalpop0] = 4e-2
                    gdat.stdp[gdat.fitt.inxparagenrscalelemkind.bgalpop0] = 4e-2
                    gdat.stdp[gdat.fitt.inxparagenrscalelemkind.nobjpop0] = 3e-1
                    try:
                        gdat.stdp[gdat.indxstdppop0gwdt] = 5e-1
                    except:
                        pass

            if gdat.typeexpr == 'fire':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
        
        if gdat.boolsqzeprop:
            gdat.stdp[:]= 1e-100
        
        if gdat.boolexplprop:
            gdat.stdp[:] = 1e100

        if (gdat.stdp > 1e100).any():
            raise Exception('')
            
        if (gdat.stdp == 0).any():
            raise Exception('')
            
        if gdat.stdp.size != gdat.numbstdp or gdat.indxstdp.size != gdat.stdp.size:
            print('gdat.stdp')
            summgene(gdat.stdp)
            print('gdat.numbstdp')
            print(gdat.numbstdp)
            print('gdat.indxstdp')
            print(gdat.indxstdp)
            raise Exception('')

        if gdat.typeverb > 1:
            # temp
            for strgmodl in gdat.liststrgmodl:
                print('strgmodl')
                print(strgmodl)
                print('Fixed dimensional parameters:')
                print('%20s%25s%5s%20s%20s' % ('name', 'labltotl', 'scal', 'minm', 'maxm'))
                for k in gmod.indxparagenrbase:
                    print('%20s%25s%5s%20.6g%20.6g' % (gmod.nameparagenrbase[k], gmod.labltotlpara.genrbase[k], gmod.scalpara.genrbase[k], \
                                                                                gmod.minmpara.genrbase[k], gmod.maxmpara.genrbase[k]))
                
                print('Element parameters')
                print('%20s%20s' % ('nameparagenrelem', 'scalcomp'))
                for l in gmod.indxpopl:
                    for nameparagenrelem, scalcomp in zip(gmod.namepara.genrelem[l], gmod.listscalparagenrelem[l]):
                        print('%20s%20s' % (nameparagenrelem, scalcomp))
                
                print('%20s%20s' % ('strgmodu', 'pdfnmodu'))
                for l in gmod.indxpopl:
                    for strgmodu, pdfnmodu in zip(gmod.namepara.genrelemmodu[l], gmod.liststrgpdfnmodu[l]):
                        print('%20s%20s' % (strgmodu, pdfnmodu))
                
                print('%20s%20s' % ('strgfeat', 'pdfnprio'))
                for l in gmod.indxpopl:
                    for strgfeat, pdfnprio in zip(gmod.namepara.genrelem[l], gmod.listscalparagenrelem[l]):
                        print('%20s%20s' % (strgfeat, pdfnprio))
                
        # proposals
        # terms in the log-acceptance probability
        gdat.listnametermlacp = []
        gdat.listlabltermlacp = []
        for l in gmod.indxpopl:
            if gmod.numbpopl > 1:
                strgpopl = '%d,' % l
            else:
                strgpopl = ''
            for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                labl = getattr(gmod.lablrootpara, nameparagenrelem)
                gdat.listlabltermlacp += ['$u_{%s%s}$' % (strgpopl, labl)]
        gdat.listnametermlacp += ['ltrp']
        gdat.listlabltermlacp += [u'$\ln P(q)$']
        gdat.listnametermlacp += ['ljcb']
        gdat.listlabltermlacp += [r'$\ln \alpha_j$']
        
        gdat.numbtermlacp = len(gdat.listnametermlacp)
        gdat.indxtermlacp = np.arange(gdat.numbtermlacp)
        
        if gdat.probtran is None:
            if gmod.numbparaelem > 0:
                gdat.probtran = 0.4
            else:
                gdat.probtran = 0.
        if gdat.probspmr is None:
            if gmod.numbparaelem > 0:
                gdat.probspmr = gdat.probtran / 2.
            else:
                gdat.probspmr = 0.
        
        gdat.probbrde = 1. - gdat.probspmr

        if gdat.probbrde < 0:
            raise Exception('')
        gdat.lablproptype = ['Within']
        gdat.nameproptype = ['with']
        if gmod.numbparaelem > 0:
            gdat.lablproptype += ['Birth', 'Death', 'Split', 'Merge']
            gdat.nameproptype += ['brth', 'deth', 'splt', 'merg']
        gdat.numbproptype = len(gdat.lablproptype)
        gdat.nameproptype = np.array(gdat.nameproptype)
        cntr = tdpy.cntr()
        if gmod.numbparaelem > 0.:
            # birth
            gdat.indxproptypebrth = cntr.incr()
            # death
            gdat.indxproptypedeth = cntr.incr()
            if gdat.probspmr > 0.:
                # split
                gdat.indxproptypesplt = cntr.incr()
                # merge
                gdat.indxproptypemerg = cntr.incr()
   
        gdat.indxproptype = np.arange(gdat.numbproptype)
        gmod.indxpara.prop = np.arange(gmod.numbparagenrbase)
        gdat.numbstdpparagenrscalbase = gmod.numbparagenrbase - gmod.numbpopl
        #### filter for model elements
        gdat.listnamefilt = ['']
        if gdat.priofactdoff != 1.:
            gdat.listnamefilt += ['pars']
        #### model elements inside the image
        if gdat.boolelempsfnanyy:
            gdat.listnamefilt += ['bndr']
        #### model subhalos inside high normalized relevance region
        if 'lens' in gdat.typeelem:
            gdat.listnamefilt += ['nrel']
        
        if gdat.typedata == 'inpt':
            proc_cntpdata(gdat)
        
        # interpolated prior for models
        for strgmodl in gdat.liststrgmodl:
        
            gmod = getattr(gdat, strgmodl)
        
            lpdfprio = [None for l in gmod.indxpopl]
            lpdfprioobjt = [None for l in gmod.indxpopl]
            lpdfpriointp = [None for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                for strgfeat, strgpdfn in zip(gmod.namepara.genrelem, gmod.listscalparagenrelem):
                    if strgpdfn == 'tmplgrad':
                        pdfnpriotemp = np.empty((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                        lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                        lpdfpriointp = lpdfprioobjt(gdat.meanpara.bgalcart, gdat.meanpara.lgalcart)
            
        gdat.indxpoplcrin = 0
        if gmod.numbparaelem > 0:
            if gdat.rtagmock is not None:
                path = gdat.pathoutprtagmock + 'gdatfinlpost'
                gdatmock = readfile(path)
            gdat.liststrgvarbhist = []
            cntr = 0
            for l0 in gmod.indxpopl:
                for a, strgfeatfrst in enumerate(gmod.namepara.genrelem[l0]):
                    if strgfeatfrst == 'spec':
                        continue
                    gdat.liststrgvarbhist.append([[] for k in range(5)])
                    gdat.liststrgvarbhist[cntr][0] = 'hist' + strgfeatfrst + 'pop%d' % l
                    gdat.liststrgvarbhist[cntr][1] = strgfeatfrst
                    if gdat.rtagmock is not None:
                        # cmpl
                        gdat.liststrgvarbhist[cntr][3] = [[] for qq in gdatmock.indxrefr]
                        # fdis
                        gdat.liststrgvarbhist[cntr][4] = [[] for qq in gdatmock.indxrefr]
                        booltemp = True
                        if strgfeatfrst[-4:] in gdat.listnamerefr:
                            q = gdat.listnamerefr.index(strgfeatfrst[-4:])
                            booltemp = not strgfeatfrst in gdat.refr.namepara.elemonly[q][l]
                        if booltemp:
                            gdat.liststrgvarbhist[cntr][3][qq] = strgfeatfrst + 'pop%dpop%d' % (l, qq)
                            gdat.liststrgvarbhist[cntr][4][qq] = strgfeatfrst + 'pop%dpop%d' % (qq, l)
                    cntr += 1    
                    for b, strgfeatseco in enumerate(gmod.namepara.genrelem[l0]):
                        
                        if strgfeatseco == 'spec':
                            continue

                        if not checstrgfeat(strgfeatfrst, strgfeatseco):
                            continue
                                        
                        gdat.liststrgvarbhist.append([[] for k in range(5)])
                        gdat.liststrgvarbhist[cntr][0] = 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l0
                        gdat.liststrgvarbhist[cntr][1] = strgfeatfrst
                        gdat.liststrgvarbhist[cntr][2] = strgfeatseco
                        gdat.liststrgvarbhist[cntr][3] = [[] for qq in gdat.indxrefr]
                        gdat.liststrgvarbhist[cntr][4] = [[] for qq in gdat.indxrefr]
                        if gdat.rtagmock is not None:
                            booltempfrst = True
                            booltempseco = True
                            if strgfeatfrst[-4:] in gdat.listnamerefr:
                                q = gdat.listnamerefr.index(strgfeatfrst[-4:])
                                booltempfrst = not strgfeatfrst in gdat.refr.namepara.elemonly[q][l]
                            if strgfeatseco[-4:] in gdat.listnamerefr:
                                q = gdat.listnamerefr.index(strgfeatseco[-4:])
                                booltempseco = not strgfeatseco in gdat.refr.namepara.elemonly[q][l]
                            for qq in gdatmock.indxrefr:
                                if booltempfrst and booltempseco:
                                    gdat.liststrgvarbhist[cntr][3][qq] = strgfeatfrst + strgfeatseco + 'pop%dpop%d' % (l0, qq)
                                    gdat.liststrgvarbhist[cntr][4][qq] = strgfeatfrst + strgfeatseco + 'pop%dpop%d' % (qq, l0)
                                elif booltempfrst:
                                    gdat.liststrgvarbhist[cntr][3][qq] = strgfeatfrst + 'pop%dpop%d' % (l0, qq)
                                    gdat.liststrgvarbhist[cntr][4][qq] = strgfeatfrst + 'pop%dpop%d' % (qq, l0)
                                elif booltempseco:
                                    gdat.liststrgvarbhist[cntr][3][qq] = strgfeatseco + 'pop%dpop%d' % (l0, qq)
                                    gdat.liststrgvarbhist[cntr][4][qq] = strgfeatseco + 'pop%dpop%d' % (qq, l0)
                        cntr += 1    
        
        # selection effects
        if gdat.typedata == 'inpt' and gmod.numbparaelem > 0:
            if gdat.numbsampboot is None:
                gdat.numbsampboot = gdat.numbsamp
        
            gdat.boolcrex = False
            if gdat.rtagmock is not None:
                for qq in gdatmock.indxrefr:
                    for q in gdat.indxrefr:
                        for l in gmod.indxpopl:
                            for strgfeatfrst in gmod.namepara.genrelem[l]:
                                
                                if gdat.typeexpr == 'chan' and strgfeatfrst == 'redswo08':
                                    crex = (1. + gdat.meanpara.redswo08)**2
                                else:
                                    crex = None
                                
                                setattr(gdat, 'crex' + strgfeatfrst + 'pop%dpop%dpop%d' % (q, qq, l), crex)
                                
                                for strgfeatseco in gmod.namepara.genrelem[l]:
                                    
                                    if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                        continue
                                                
                                    if gdat.typeexpr == 'chan' and (strgfeatfrst == 'redswo08' or strgfeatseco == 'redswo08'):
                                        crex = np.empty((gdat.numbbinsplot, gdat.numbbinsplot))
                                        if strgfeatfrst == 'redswo08':
                                            crex[:, :] = (1. + gdat.meanpara.redswo08[:, None])**2
                                        else:
                                            crex[:, :] = (1. + gdat.meanpara.redswo08[None, :])**2
                                    else:
                                        crex = None
                                    
                                    setattr(gdat, 'crex' + strgfeatfrst + strgfeatseco + 'pop%dpop%dpop%d' % (q, qq, l), crex)
        
                if gdat.refr.numbelemtotl > 0:
                    for listtemp in gdat.liststrgvarbhist:
                        strgvarb = listtemp[0]
                        for qq in gdatmock.indxrefr:
                            for q in gdat.indxrefr:
                                nametemp = listtemp[1]
                                if len(listtemp[2]) > 0:
                                    nametemp += listtemp[2]
                                l = int(listtemp[4][qq].split('pop')[2][0])
                                nametemp += 'pop%dpop%dpop%d' % (q, qq, l)
                                crexhist = getattr(gdat, 'crex' + nametemp)
                                if crexhist is not None:
                                    gdat.boolcrex = True
            
            ## internal correction
            gdat.boolcrin = gdat.typedata == 'inpt' and gdat.rtagmock is not None
        
        if gmod.numbparaelem > 0:
            # variables for which two dimensional functions will be plotted
            gdat.liststrgelemtdimvarbinit = ['hist']
            gdat.liststrgelemtdimvarbfram = deepcopy(gdat.liststrgelemtdimvarbinit)
            if gdat.boolinforefr:
                gdat.liststrgelemtdimvarbfram += ['cmpl', 'fdis']
            gdat.liststrgelemtdimvarbfinl = deepcopy(gdat.liststrgelemtdimvarbfram)
            if gdat.typedata == 'inpt':
                if gdat.boolcrex:
                    gdat.liststrgelemtdimvarbfinl += ['excr']
                if gdat.boolcrin:
                    gdat.liststrgelemtdimvarbfinl += ['incr']
            gdat.liststrgelemtdimvarbanim = deepcopy(gdat.liststrgelemtdimvarbfram)
        
        gdat.liststrgfoldinit = ['']
        if gmod.numbparaelem > 0 or gdat.typedata == 'mock' and gmod.numbparaelem > 0:
            gdat.liststrgfoldinit += ['', 'histodim/', 'histtdim/', 'scattdim/', 'cmpltdim/']
        gdat.liststrgfoldfram = ['']
        if gmod.numbparaelem > 0:
            gdat.liststrgfoldfram += ['scattdim/']
        gdat.liststrgfoldfinl = ['']
        if gdat.boolinforefr and gmod.numbparaelem > 0:
            gdat.liststrgfoldfram += ['assc']
            gdat.liststrgfoldfinl += ['assc']
        gdat.liststrgfoldanim = deepcopy(gdat.liststrgfoldfram)

        if gmod.numbparaelem > 0:
            for strgdims in ['odim/', 'tdim/']:
                for strgelemtdimvarb in gdat.liststrgelemtdimvarbfram:
                    gdat.liststrgfoldfram += [strgelemtdimvarb + strgdims]
                for strgelemtdimvarb in gdat.liststrgelemtdimvarbfinl:
                    gdat.liststrgfoldfinl += [strgelemtdimvarb + strgdims]

        # make folders
        #gdat.pathprio = gdat.pathplotrtag + 'prio/'
        #gdat.pathpost = gdat.pathplotrtag + 'post/'
        make_fold(gdat)

        setp_indxswepsave(gdat)
            
        if gdat.typeopti == 'hess':
            pathopti = gdat.pathoutprtag + 'opti.h5'
            if os.path.exists(pathopti):
                thisfile = h5py.File(pathopti, 'r')
                if thisfile['stdp'][()].size == gdat.stdp.size:
                    print('Recovering the proposal scale from the previous run...')
                    gdat.stdp = thisfile['stdp'][()]
                thisfile.close()

        if gdat.rtagmock is not None:
            if gdat.typedata == 'inpt':
                path = gdat.pathoutprtagmock + 'gdatfinlpost'
                booltemp = True
                try:
                    gdatmock = readfile(path)
                except:
                    booltemp = False
                    gdat.rtagmock = None

                if booltemp:
                    numbparaelem = gdatmock.true.numbparaelem
                    if gdatmock.trueindxpopl != gmod.indxpopl:
                        raise Exception('')
                    for l in gmod.indxpopl:
                        for strgfeat in gmod.namepara.genrelem[l]:
                            if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':
                                continue

                            if strgfeat[-4:] in gdat.listnamerefr:
                                continue
                            reca = getattr(gdatmock.true, 'reca' + strgfeat + 'pop%d' % l)
                            setattr(gdat.true, 'reca' + strgfeat + 'pop%d' % l, reca)
                    gmod.namepara.genrelem = gdatmock.truegmod.namepara.genrelem
        
        if gmod.typeelemspateval[l] == 'locl' and gmod.numbparaelem > 0 or \
                            gdat.typedata == 'mock' and gmod.typeelemspateval[l] == 'locl' and gmod.numbparaelem > 0:
            gdat.numbprox = 3
            gdat.indxprox = np.arange(gdat.numbprox)
            minmparagenrscalelemampl = getattr(gdat.fitt.minmpara, gmod.nameparagenrelemampl[0])
            maxmparagenrscalelemampl = getattr(gdat.fitt.maxmpara, gmod.nameparagenrelemampl[0])
            gdat.binspara.prox = np.logspace(np.log10(minmparagenrscalelemampl), np.log10(maxmparagenrscalelemampl), gdat.numbprox + 1)
            
            # determine the maximum angle at which the contribution of the element will be computed
            if gdat.boolbinsspat:
                if gdat.maxmangleval is None:
                    if gdat.typeexpr == 'chan':
                        gdat.maxmangleval = np.array([5., 6., 9.]) / gdat.anglfact
                    elif gdat.typeexpr == 'gene':
                        gdat.maxmangleval = np.array([0.1, 0.2, 0.3]) / gdat.anglfact
                    elif gdat.typeexpr == 'ferm':
                        gdat.maxmangleval = np.array([7., 9., 15.]) / gdat.anglfact
                    else:
                        gdat.maxmangleval = np.empty(gdat.numbprox)
                        for h in gdat.indxprox:
                            if gdat.specfraceval == 0:
                                gdat.maxmangleval[h] = 3. * gdat.maxmgang
                            else:  
                                frac = min(1e-2, gdat.specfraceval * gdat.binspara.prox[0] / gdat.binspara.prox[h+1])
                                psfnwdth = retr_psfnwdth(gdat, gmodstat.psfn, frac)
                                gdat.indxmaxmangl = np.unravel_index(np.argmax(psfnwdth), psfnwdth.shape)
                                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
            
            if gdat.typeverb > 1:
                if gmod.typeelemspateval == 'locl':
                    print('maxmangleval')
                    print(gdat.anglfact * gdat.maxmangleval[l], ' [%s]' % gdat.strganglunit)

            setp_varb(gdat, 'angl', minm=0., maxm=10.)
            if gdat.boolelempsfnanyy and gdat.maxmpara.angl < np.amax(gdat.maxmangleval):
                print('gdat.maxmpara.angl')
                print(gdat.maxmpara.angl)
                print('gdat.maxmangleval')
                print(gdat.maxmangleval)
                
                raise Exception('Angular axis is too short.')

            # make a look-up table of nearby pixels for each pixel
            path = gdat.pathprox + 'indxprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.typepixl, 1e2 * np.amin(gdat.maxmangleval), \
                                                                                                                1e2 * np.amax(gdat.maxmangleval), gdat.numbprox)
            
            if gdat.typeverb > 1:
                print('gdat.typepixl')
                print(gdat.typepixl)
                print('gdat.minmlgaldata')
                print(gdat.minmlgaldata)
                print('gdat.minmbgaldata')
                print(gdat.minmbgaldata)
                print('gdat.maxmlgaldata')
                print(gdat.maxmlgaldata)
                print('gdat.maxmbgaldata')
                print(gdat.maxmbgaldata)
            if gdat.typeverb > 0:
                print('Element evaluation will be performed up to')
                if gdat.boolbinsspat:
                    print(gdat.maxmangleval * gdat.anglfact)

            if os.path.isfile(path):
                if gdat.typeverb > 0:
                    print('Previously computed nearby pixel look-up table will be used.')
                    print('Reading %s...' % path)
                fobj = open(path, 'rb')
                gdat.indxpixlprox = pickle.load(fobj)
                fobj.close()
            else:
                if gdat.typeverb > 0:
                    print('Computing the look-up table...')
                gdat.indxpixlprox = [[] for h in gdat.indxprox]
                cntrsave = -1.
                # temp
                for j in gdat.indxpixl:
                    dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
                    dist[j] = 0.
                    for h in gdat.indxprox:
                        indxpixlproxtemp = np.where(dist < gdat.maxmangleval[h])[0]
                        if indxpixlproxtemp.size > 2e4:
                            indxpixlproxtemp = -1
                            if gdat.maxmangl < np.sqrt(2.) * gdat.maxmgangdata:
                                raise Exception('Angular axis used to interpolate the PSF should be longer.')
                        
                        if indxpixlproxtemp.size < 10:
                            raise Exception('Pixel hash list should not have fewer than 10 pixels.')

                        gdat.indxpixlprox[h].append(indxpixlproxtemp)
                    cntrsave = tdpy.show_prog(j, gdat.numbpixl, cntrsave)
                fobj = open(path, 'wb')
                pickle.dump(gdat.indxpixlprox, fobj, protocol=pickle.HIGHEST_PROTOCOL)
                fobj.close()
            
            gdat.numbpixlprox = np.zeros(gdat.numbprox) 
            for h in gdat.indxprox:
                for j in gdat.indxpixl:
                    gdat.numbpixlprox[h] += len(gdat.indxpixlprox[h][j])
            gdat.numbpixlprox[h] /= len(gdat.indxpixlprox[h])
            
            if (gdat.numbpixlprox - np.mean(gdat.numbpixlprox) == 0.).all():
                raise Exception('Number of pixels in the hash lists should be different.')

        gdat.minmgang = 1e-3 * np.sqrt(2.) * gdat.maxmgangdata
        gdat.maxmgang = np.sqrt(2.) * gdat.maxmgangdata
        
        # try to pass true metamodel minima and maxima to common minima and maxima when that feature does not exist in the fitting metamodel
        if gdat.typedata == 'mock':
            for q in gdat.indxrefr:
                for strgfeat in gmod.namepara.genrelem[q]:
                    booltemp = False
                    for l in gmod.indxpopl:
                        if strgfeat in gmod.namepara.genrelem[l]:
                            booltemp = True
                    if not booltemp:
                        try:
                            setattr(gdat.minmpara, 'minm' + strgfeat + gdat.listnamerefr[q], getattr(gdat.true.minm, strgfeat))
                            setattr(gdat.maxmpara, 'maxm' + strgfeat + gdat.listnamerefr[q], getattr(gdat.true.maxm, strgfeat))
                        except:
                            pass

        ## reference spectra
        if gdat.listprefsbrtlabltotl is None:
            if gdat.typeexpr == 'chan' and gdat.boolbinsspat:
                gdat.listprefsbrtener = [[[] for k in range(3)]]
                gdat.listprefsbrtsbrt = [[[] for k in range(3)]]
                gdat.listprefsbrtlabltotl = ['Moretti+(2012)']
                gdat.listprefsbrttype = ['shad']
                
                for k, strgextn in enumerate(['', '_lower', '_higher']):
                    path = gdat.pathinpt + 'Moretti2012%s.csv' % strgextn
                    enerrefrplot = np.loadtxt(path, delimiter=',')[:, 0]
                    sbrtrefrplot = np.loadtxt(path, delimiter=',')[:, 1] / gdat.factergskevv / enerrefrplot**2 * (180. / np.pi)**2
                    gdat.listprefsbrtener[0][k] = enerrefrplot
                    gdat.listprefsbrtsbrt[0][k] = sbrtrefrplot

        # temp
        if gdat.numbener > 1:
            if gdat.enerpivt == 0.:
                raise Exception('Pivot energy cannot be zero.')
            #if gdat.typeexpr != 'fire':
            #    gdat.enerexpcfact = gdat.enerpivt - gdat.meanpara.ener
            #if gmod.numbparaelem > 0 and gdat.numbener > 1:
            #    minmsinddistmeanpop0 = getattr(gmod, 'minmsinddistmeanpop0')
            #    factspecener = (gdat.meanpara.ener / gdat.enerpivt)**(-np.sqrt(np.amin(minmsinddistmeanpop0) * np.amax(maxmsinddistmeanpop0)))
        else:
            pass
            #gdat.factspecener = np.array([1.])

        # temp -- this assumes square ROI
        if gdat.boolbinsspat:
            gdat.frambndrmodl = gdat.maxmlgaldata * gdat.anglfact
        
        if gmod.boollenshost or gdat.typedata == 'mock' and gmod.boollenshost:
            
            if gdat.typesers == 'intp':
                # construct pixel-convolved Sersic surface brightness template
                gdat.factsersusam = 10
                maxmlgal = 4. * np.sqrt(2.) * gdat.maxmlgal
                gdat.numblgalsers = int(np.ceil(maxmlgal / gdat.sizepixl))
                gdat.numblgalsersusam = (1 + gdat.numblgalsers) * gdat.factsersusam
                retr_axis(gdat, 'lgalsers')
                retr_axis(gdat, 'lgalsersusam')
                retr_axis(gdat, 'bgalsersusam')
                
                gdat.numbhalfsers = 20
                gdat.numbindxsers = 20
                    
                retr_axis(gdat, 'halfsers')
                retr_axis(gdat, 'indxsers')
                
                gdat.binspara.lgalsersusammesh, gdat.binspara.bgalsersusammesh = np.meshgrid(gdat.binspara.lgalsersusam, gdat.binspara.bgalsersusam, indexing='ij')
                gdat.binspara.radisersusam = np.sqrt(gdat.binspara.lgalsersusammesh**2 + gdat.binspara.bgalsersusammesh**2)
                 
                gdat.sersprofcntr = np.empty((gdat.numblgalsers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
                gdat.sersprof = np.empty((gdat.numblgalsers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
                
                for n in range(gdat.numbindxsers + 1):
                    for k in range(gdat.numbhalfsers + 1):
                        
                        profusam = retr_sbrtsersnorm(gdat.binspara.radisersusam, gdat.binspara.halfsers[k], indxsers=gdat.binspara.indxsers[n])
        
                        ## take the pixel average
                        indxbgallowr = gdat.factsersusam * (gdat.numblgalsers + 1) / 2
                        indxbgaluppr = gdat.factsersusam * (gdat.numblgalsers + 3) / 2
                        for a in range(gdat.numblgalsers):
                            indxlgallowr = gdat.factsersusam * a
                            indxlgaluppr = gdat.factsersusam * (a + 1) + 1
                            gdat.sersprofcntr[a, k, n] = profusam[(indxlgallowr+indxlgaluppr)/2, 0]
                            gdat.sersprof[a, k, n] = np.mean(profusam[indxlgallowr:indxlgaluppr, :])
                
                temp, indx = unique(gdat.binspara.lgalsers, return_index=True)
                gdat.binspara.lgalsers = gdat.binspara.lgalsers[indx]
                gdat.sersprof = gdat.sersprof[indx, :, :]
                gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]
        
                indx = np.argsort(gdat.binspara.lgalsers)
                gdat.binspara.lgalsers = gdat.binspara.lgalsers[indx]
                gdat.sersprof = gdat.sersprof[indx, :, :]
                gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]

        #for strg, valu in gmod.cmappara.__dict__.items():
        #    retr_ticklabl(gdat, strg)
                
        # generate true data
        if gdat.typedata == 'mock':
            
            if gdat.typeverb > 0:
                print('Generating mock data...')

            if gdat.seedtype == 'rand':
                np.random.seed()
            else:
                if gdat.typeverb > 0:
                    print('Setting the seed for the RNG to %d...' % gdat.seedtype)
                np.random.seed(gdat.seedtype)
        
            ## unit sample vector
            gdat.true.this.paragenrunitfull = np.random.rand(gdat.true.numbparagenrfull)
            gdat.true.this.paragenrscalfull = np.zeros(gdat.true.numbparagenrfull)
            
            if gdat.true.numbparaelem > 0:
                gdat.true.this.numbelempopl = np.empty(gdat.true.maxmpara.numbelem[l], dtype=int)
                for l in gdat.true.indxpopl:
                    gdat.true.this.paragenrunitfull[gdat.true.indxpara.numbelem[l]] = getattr(gdat.true.this, 'numbelempop%d' % l)
                    gdat.true.this.numbelempopl[l] = getattr(gdat.true.this, 'numbelempop%d' % l)

                gdat.true.this.indxelemfull = [[] for l in gdat.true.indxpopl]
                for l in gdat.true.indxpopl:
                    gdat.true.this.indxelemfull[l] = list(range(gdat.true.numbelem[l]))
                gdat.true.this.indxparagenrfullelem = retr_indxparagenrfullelem(gdat, gdat.true.this.indxelemfull, 'true')
            else:
                gdat.true.this.indxelemfull = []
                gdat.true.this.indxparagenrfullelem = None

            if gdat.true.numbparaelem > 0:
                if gdat.seedelem is None:
                    np.random.seed()
                else:
                    np.random.seed(gdat.seedelem)
                gdat.true.this.paragenrunitfull[gdat.true.numbparagenrbase:] = np.random.rand(gdat.true.numbparagenrelemtotl)
            
            gdat.true.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'true', gdat.true.this.paragenrunitfull, gdat.true.this.indxparagenrfullelem)
            
            # impose true values (valu)
            for k in gdat.true.indxparagenr:
                
                if gdat.true.numbparaelem > 0 and (k in gdat.true.indxpara.numbelem or \
                                            gdat.true.typemodltran == 'pois' and k in gdat.true.indxpara.meanelem):
                        continue
        
                # assume the true PSF
                if gdat.true.typeevalpsfn != 'none' and gdat.numbpixl > 1 and k in gdat.true.indxpara.psfp:
                    gdat.true.this.paragenrscalfull[k] = gdat.true.psfpexpr[k-gdat.true.indxpara.psfp[0]]
                else:
                    ## read input mock model parameters
                    try:
                        # impose user-defined true parameter
                        gdat.true.this.paragenrscalfull[k] = getattr(gdat.true, gdat.true.namepara.genrscalfull[k])
                    except:
                        pass
        
            if gdat.typeverb > 0:
                show_paragenrscalfull(gdat, None, strgmodl='true')

            if gmod.boollenshost:
                proc_samp(gdat, None, 'this', 'true', boolinit=True)
            
            #for strgmodl in gdat.liststrgmodl:
            #    gmod = getattr(gdat, strgmodl)
            #    print('gmod.minmpara.numbelempop0')
            #    print(gmod.minmpara.numbelempop0)
            #    print('gmod.minmpara.numbelem')
            #    print(gmod.minmpara.numbelem)
            #raise Exception('')
        
            # construct bins for element parameters of the true model
            for strgmodl in ['true']:
                
                gmod = getattr(gdat, strgmodl)

                # list of names for element parameters, concatenated across all populations
                for l in gmod.indxpopl:
                    if gmod.maxmpara.numbelem[l] > 0:
                        # temp -- does not cover the case when different populations have parameters with the same name
                        for strgfeat in gmod.listnameparaglob:
                        #for strgfeat in gmod.namepara.genrelem[l]:
                            if strgfeat[:-4] == 'etag':
                                continue
                            #retr_axis(gdat, strgfeat)
                            #if strgfeat in gmod.namepara.elem:
                            #    retr_axis(gdat, strgfeat + 'prio')
        
            proc_samp(gdat, None, 'this', 'true', boolinit=True)
        
            # transfer current state of the true model to the reference model 
            for strg, valu in gdat.true.this.__dict__.items():
                if strg == 'dictelem':
                    # modify the current state of the element parameters of the true model to include uncertainty
                    valutemp = [[] for l in gdat.true.indxpopl]
                    for l in gdat.true.indxpopl:
                        valutemp[l] = dict()
                        for nameparaelem in gdat.true.this.dictelem[l]:
                            valutemp[l][nameparaelem] = np.zeros((3, gdat.true.this.dictelem[l][nameparaelem].size))
                            valutemp[l][nameparaelem][0, :] = gdat.true.this.dictelem[l][nameparaelem]
                else:
                    valutemp = valu
                setattr(gdat.refr, strg, valutemp)
            if gdat.makeplot and gdat.makeplotinit:
                plot_samp(gdat, None, 'this', 'true', 'init')
            
        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            print('gmod.minmpara.numbelempop0')
            print(gmod.minmpara.numbelempop0)
            print('gmod.minmpara.numbelem')
            print(gmod.minmpara.numbelem)
        
        ## initialization
        gdat.fitt.this = tdpy.gdatstrt()
        gdat.fitt.next = tdpy.gdatstrt()
        init_stat(gdat)
        
        # process the parameter vector
        proc_samp(gdat, None, 'this', 'fitt', boolinit=True)
    
        #liststrgcbar = ['llikmaps', 'perc', 'percresi', 'expo', 'lpdfspatpriointp', 'conv', 'magn', 'deflcomp', 'resiconvelem', 'resimagn']
        #for strgcbar in liststrgcbar:
        #    retr_ticklabl(gdat, strgcbar)
        
        # temp
        #for strgmodl in gdat.liststrgmodl:
        #    for namesele in gdat.listnamesele:
        #        for namefeat in gdat.listnamefeatsele:
        #            for strglimt in gdat.liststrglimt:
        #                try:
        #                    getattr(gdat, strglimt + namefeat + namesele)
        #                except:
        #                    setattr(gdat, strglimt + namefeat + namesele, getattr(gdat, strglimt + namefeat))

        # construct bins for element parameters of the fitting model
        #for strgmodl in ['fitt']:
        #    
        #    gmod = getattr(gdat, strgmodl)

        #    # list of names for element parameters, concatenated across all populations
        #    for l in gmod.indxpopl:
        #        if gmod.maxmpara.numbelem[l] > 0:
        #            # temp -- does not cover the case when different populations have parameters with the same name
        #            for strgfeat in gmod.listnameparaglob:
        #            #for strgfeat in gmod.namepara.genrelem[l]:
        #                if strgfeat[:-4] == 'etag':
        #                    continue
        #                #retr_axis(gdat, strgfeat)
        #                #if strgfeat in gmod.namepara.elem:
        #                #    retr_axis(gdat, strgfeat + 'prio')
        
        gdat.numbbinspdfn = 50
                
        # scalar variable setup continued
        for strgbins in ['lowr', 'higr']:
            for strgecom in ['dfnc', 'dfncsubt']:
                setattr(gdat, 'scalhistcntp' + strgbins + strgecom + 'en00evt0', 'self')
                setattr(gdat, 'minmhistcntp' + strgbins + strgecom + 'en00evt0', 0.)
                setattr(gdat, 'maxmhistcntp' + strgbins + strgecom + 'en00evt0', gdat.numbpixl)
                setattr(gdat, 'facthistcntp' + strgbins + strgecom + 'en00evt0', 1.)
        for i in gdat.indxener:
            setattr(gdat, 'scalfracsdenmeandarkdfncsubten%02d' % i, 'self')
            setattr(gdat, 'minmfracsdenmeandarkdfncsubten%02d' % i, 0.)
            setattr(gdat, 'maxmfracsdenmeandarkdfncsubten%02d' % i, 1.)
            setattr(gdat, 'factfracsdenmeandarkdfncsubten%02d' % i, 1.)
        
        gmod.scalbooldfncsubt = 'self'
        gdat.minmbooldfncsubt = -0.5
        gdat.maxmbooldfncsubt = 1.5
        gdat.factbooldfncsubt = 1.

        #sys.stdout = logg(gdat)
        #gdat.log.close()

        # initial plots
        if gdat.makeplot and gdat.makeplotinit:
            plot_init(gdat)

        if gdat.typeverb > 0:
            sizetotl = 0.
            for root, dirs, listfile in os.walk(gdat.pathoutp):
                for thisfile in listfile:
                    sizetotl += os.path.getsize(root + '/' + thisfile) / 2**30
            if sizetotl > 10.:
                print('Warning: PCAT data path size is %d GB' % sizetotl)

        if gdat.typedata == 'inpt':
        
            ## rotate element coordinates to the ROI center
            if gdat.typepixl == 'heal' and (gdat.lgalcntr != 0. or gdat.bgalcntr != 0.):
                for q in gdat.indxrefr:
                    for l in gmod.indxpopl:
                        rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='ZYX')
                        gdat.refr.dictelem[q]['bgal'][0, :], gdat.refrlgal[0, :] = rttr(pi / 2. - gdat.refrbgal[0, :], gdat.refrlgal[0, :])
                        gdat.refr.dictelem[q]['bgal'][0, :] = pi / 2. - gdat.refrbgal[0, :]

            ## assign zero to nonspecified uncertainties for the reference element features
            for q in gdat.indxrefr:
                for strgfeat in gdat.refr.namepara.elem[q]:
                    if strgfeat == 'gang' or strgfeat == 'aang':
                        continue
                    if strgfeat == 'etag':
                        continue
                    refrfeat = getattr(gdat.refr, strgfeat)
                    if refrfeat[q].ndim == 1:
                        refrfeat[q] = np.tile(refrfeat[q], (3, 1)) 
            
        # temp
        #if gdat.refr.numbelem > 0:
        #    gdat.refrfluxbrgt, gdat.refrfluxbrgtassc = retr_fluxbrgt(gdat, gdat.refrlgal, gdat.refrbgal, gdat.refrflux[0, :])
        
        print('gdat.liketype')
        print(gdat.liketype)

        print('Data settings')
        print('gdat.numbener')
        print(gdat.numbener)
        print('gdat.numbevtt')
        print(gdat.numbevtt)

        print('Model settings')
        print('gdat.fitt.numbpopl')
        print(gdat.fitt.numbpopl)
        print('gdat.fitt.numbparagenrbase')
        print(gdat.fitt.numbparagenrbase)
        
        for strgmodl in gdat.liststrgmodl:
            for l in gmod.indxpopl:
                for strgfeat, strgpdfn in zip(gmod.namepara.genrelemmodu[l], gmod.liststrgpdfnmodu[l]):
                    if strgpdfn == 'tmpl':
                        if gdat.lgalprio is None or gdat.bgalprio is None:
                            gdat.lgalprio = np.concatenate((gmod.lgal))
                            gdat.bgalprio = np.concatenate((gmod.bgal))
                        gdat.numbspatprio = gdat.lgalprio.size
        
                        # spatial template for the catalog prior
                        # temp -- this should move outside the if
                        gdat.pdfnspatpriotemp = np.zeros((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                        for k in range(gdat.numbspatprio):
                            gdat.pdfnspatpriotemp[:] += 1. / np.sqrt(2. * np.pi) / gdat.stdvspatprio * \
                                                                exp(-0.5 * (gdat.binspara.lgalcartmesh - gdat.lgalprio[k])**2 / gdat.stdvspatprio**2) * \
                                                                exp(-0.5 * (gdat.binspara.bgalcartmesh - gdat.bgalprio[k])**2 / gdat.stdvspatprio**2)
                        gdat.pdfnspatpriotemp /= np.amax(gdat.pdfnspatpriotemp)
        
        if gdat.typedata == 'inpt':

            # rotate reference elements to the spatial coordinate system of PCAT
            # temp -- this does not rotate the uncertainties!

            if gdat.typeverb > 0:
                print('Rotating the reference elements...')
            for q in gdat.indxrefr:
                # temp -- this should depend on q
                if len(gdat.listpathwcss) > 0:
                    listhdun = ap.io.fits.open(gdat.listpathwcss)
                    wcso = ap.wcs.WCS(listhdun[0].header)
                    skycobjt = ap.coordinates.SkyCoord("galactic", l=gdat.refr.dictelem[q]['lgal'][0, :] * 180. / pi, b=gdat.refr.dictelem[q]['bgal'][0, :] * 180. / pi, unit='deg')
                    rasc = skycobjt.fk5.ra.degree
                    decl = skycobjt.fk5.dec.degree
                    lgal, bgal = wcso.wcs_world2pix(rasc, decl, 0)
                    lgal -= gdat.numbpixllgalshft + gdat.numbsidecarthalf
                    bgal -= gdat.numbpixlbgalshft + gdat.numbsidecarthalf
                    lgal *= gdat.sizepixl
                    bgal *= gdat.sizepixl
                    gdat.refr.dictelem[q]['lgal'][0, :] = bgal
                    gdat.refr.dictelem[q]['bgal'][0, :] = lgal

            ## preprocess reference element features
            for q in gdat.indxrefr:
                # temp -- this should depend on q
                # temp -- this does not properly calculate uncertainties
                gdat.refrgang[q] = np.zeros((3, gdat.refr.dictelem[q]['lgal'].shape[1]))
                gdat.refraang[q] = np.zeros((3, gdat.refr.dictelem[q]['lgal'].shape[1]))
                gdat.refrgang[q][:, :] = retr_gang(gdat.refr.dictelem[q]['lgal'][0, :], gdat.refr.dictelem[q]['bgal'][0, :])[None, :]
                gdat.refraang[q][:, :] = retr_aang(gdat.refr.dictelem[q]['lgal'][0, :], gdat.refr.dictelem[q]['bgal'][0, :])[None, :]

            # save all reference element features
            for strgfeat in gdat.refr.namepara.elemtotl:
                refrfeattotl = [[] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for strgfeat in gdat.refr.namepara.elem[q]:
                        refrfeat = getattr(gdat.refr, strgfeat)
                        for l in gmod.indxpopl:
                            if len(refrfeat[q]) > 0:
                                refrfeattotl[q] = refrfeat[q]
                setattr(gdat.refr, strgfeat + 'totl', refrfeattotl)
            
            # find the reference elements inside the ROI
            gdat.indxrefrpntsrofi = [[] for q in gdat.indxrefr]
            for q in gdat.indxrefr:
                gdat.indxrefrpntsrofi[q] = np.where((np.fabs(gdat.refr.dictelem[q]['lgal'][0, :]) < gdat.maxmgangdata) & \
                                                                        (np.fabs(gdat.refr.dictelem[q]['bgal'][0, :]) < gdat.maxmgangdata))[0]
            for strgfeat in gdat.refr.namepara.elemtotl:
                refrfeat = getattr(gdat.refr, strgfeat)
                refrfeatrofi = [[] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    if len(refrfeat[q]) > 0:
                        refrfeatrofi[q] = refrfeat[q][..., gdat.indxrefrpntsrofi[q]]
                setattr(gdat.refr, strgfeat, refrfeatrofi)
            
            # temp -- gdat.refr.numbelem is defined twice, one before and one after the filter. The initial definition is needed for strgfeat definitions.
            gdat.refr.numbelem = [[] for q in gdat.indxrefr]
            gdat.refr.numbelemtotl = 0
            for q in gdat.indxrefr:
                gdat.refr.numbelem[q] = 0
                gdat.refr.numbelem[q] = gdat.refr.dictelem[q]['lgal'].shape[1]
                gdat.refr.numbelem[q] = np.sum(gdat.refr.numbelem[q])
                gdat.refr.numbelemtotl += np.sum(gdat.refr.numbelem[q]) 
            
            ## check that all reference element features are finite
            for q in gdat.indxrefr:
                for strgfeat in gdat.refr.namepara.elem[q]:
                    if strgfeat == 'etag':
                        continue
                    refrfeat = getattr(gdat.refr, strgfeat)
                    if len(refrfeat[q]) > 0:
                        indxbadd = np.where(np.logical_not(np.isfinite(refrfeat[q])))
                        if indxbadd[0].size > 0:
                            refrfeat[q][indxbadd] = 0.
                            if gdat.typeverb > 0:
                                print('Warning: Provided reference element feature is not finite. Defaulting to 0...')
                        
                        if refrfeat[q].size == 0:
                            print('Warning! A reference element feature has length zero!')
                            print('strgfeat')
                            print(strgfeat)
                        else:
                            if np.amin(refrfeat[q]) == 0. and np.amax(refrfeat[q]) == 0.:
                                print('Warning! A reference element feature is all np.zeros!')
                                raise Exception('')
            
            ## element feature indices ordered with respect to the amplitude variable
            refrfeatsort = [[] for q in gdat.indxrefr]
            if not (gdat.typedata == 'mock' and gmod.numbparaelem == 0):
                for q in gdat.indxrefr:
                    refrparagenrscalelemampl = getattr(gdat.refr, gdat.refr.nameparagenrelemampl[q])
                    if len(refrparagenrscalelemampl[q]) > 0:
                        indxelem = np.argsort(refrparagenrscalelemampl[q][0, :])[::-1]
                        for strgfeat in gdat.refr.namepara.elem[q]:
                            refrfeat = getattr(gdat.refr, strgfeat)
                            if len(refrfeat[q]) > 0:
                                refrfeatsort[q] = refrfeat[q][..., indxelem]
                setattr(gdat.refr, strgfeat, refrfeatsort)
            
            # bin reference element features
            for q in gdat.indxrefr:
                for strgfeatfrst in gdat.refr.namepara.elem[q]:
                    if strgfeatfrst.startswith('etag'):
                        continue
                    refrfeatfrst = getattr(gdat.refr, strgfeatfrst)
                    if len(refrfeatfrst[q]) > 0:
                        binsfrst = getattr(gdat.binspara, strgfeatfrst)
                        hist = np.histogram(refrfeatfrst[q][0, :], binsfrst)[0]
                        setattr(gdat.refr, 'hist' + strgfeatfrst + 'pop%d' % q, hist)
                        for strgfeatseco in gdat.refr.namepara.elem[q]:
                            if strgfeatseco.startswith('etag'):
                                continue
                            refrfeatseco = getattr(gdat.refr, strgfeatseco)
                            
                            strgfeattdim = strgfeatfrst + strgfeatseco + 'pop%d' % q
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                            
                            if len(refrfeatseco[q]) > 0:
                                binsseco = getattr(gdat.binspara, strgfeatseco)
                                hist = np.histogram2d(refrfeatfrst[q][0, :], refrfeatseco[q][0, :], bins=(binsfrst, binsseco))[0]
                                setattr(gdat.refr, 'hist' + strgfeattdim, hist)
            
        if gmod.numbparaelem > 0:
            # plot settings
            ## upper limit of histograms
            if gdat.limtydathistfeat is None:
                gdat.limtydathistfeat = [0.5, max(100., 10**np.ceil(np.log10(gdat.refr.numbelemtotl)))]
                #gdat.limtydathistfeat = [0.5, max(100., 10**np.ceil(np.log10(gmod.maxmpara.numbelemtotl)))]

        # initial plots
        if gdat.makeplot and gdat.makeplotinit:
            # problem-specific plots
            if gdat.makeplotintr:
                plot_intr(gdat)
                #plot_pert()
                #plot_king(gdat)
                plot_lens(gdat)
                #plot_3fgl_thrs(gdat)
                #if gdat.typeexpr == 'ferm':
                #    plot_fgl3(gdat)
        
        # find the pixels at which data count maps have local maxima
        if gdat.typepixl == 'cart':
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    # temp
                    gdat.indxxdatmaxm, gdat.indxydatmaxm = tdpy.retr_indximagmaxm(gdat.cntpdatacart[i, :, m])
    
        if not gdat.boolsqzeexpo and np.amax(gdat.cntpdata) < 1.:
            raise Exception('Data counts per pixel is less than 1.')
        
        # check the data
        if (np.fabs(gdat.cntpdata - np.round(gdat.cntpdata)) > 1e-3).any():
            raise Exception('')
        if np.amin(gdat.cntpdata) < 0.:
            raise Exception('')
    
        # list of variables for which the posterior is collected at each proposal
        gdat.liststrgvarbarryswep = ['memoresi', 'accpprob', 'boolpropfilt', 'boolpropaccp', 'indxproptype', 'amplpert']
        for namechro in gdat.listnamechro:
            gdat.liststrgvarbarryswep += ['chro' + namechro]
        gdat.liststrgvarbarryswep += ['ltrp']
        if gdat.probtran > 0.:
            for l in gmod.indxpopl:
                gdat.liststrgvarbarryswep += ['auxiparapop%d' % l]
        gdat.liststrgvarbarryswep += ['ljcb']
    
        # write the numpy RNG state to file
        with open(gdat.pathoutprtag + 'stat.p', 'wb') as thisfile:
        	pickle.dump(np.random.get_state(), thisfile)
        
        # process lock for simultaneous plotting
        lock = mp.Manager().Lock()
        
        if gdat.typeverb > 0:
            print('Writing the global state to the disc before spawning workers...')
        path = gdat.pathoutprtag + 'gdatinit'
        writfile(gdat, path) 
        gdat.filestat = open(gdat.pathoutprtag + 'stat.txt', 'w')
        gdat.filestat.write('gdatinit written.\n')
        gdat.filestat.close()
        
        # exit before running the sampler
        if gdat.boolmockonly:
            if gdat.typeverb > 0:
                print('Mock dataset is generated. Quitting...')
            return gdat.rtag
        
        # perform an initial run, sampling from the prior
        if gdat.checprio:
            
            if gdat.typeverb > 0:
                print('Sampling from the prior...')
            
            ## perform sampling
            worksamp(gdat, lock, strgpdfn='prio')
            
            ## post process the samples
            proc_finl(gdat=gdat, strgpdfn='prio')
            
        if gdat.typeverb > 0:
            print('Sampling from the posterior...')
        
        # run the sampler
        worksamp(gdat, lock)
    
    # post process the samples
    proc_finl(gdat=gdat)
    
    # make animations
    if gdat.makeanim and gdat.numbplotfram > 1:
        proc_anim(gdat.rtag)

    if gdat.typeverb > 0:
        print('The output is at ' + gdat.pathoutprtag)
        if gdat.makeplot:
            print('The plots are at ' + gdat.pathplotrtag)
        print('PCAT has run successfully. Returning to the OS...')

    return gdat.rtag


def initarry( \
             dictvarbvari, \
             dictvarb, \
             listnamecnfgextn, \
             forcneww=False, \
             forcprev=False, \
             strgpara=False, \
             
             # Boolean flag to execute the runs in parallel
             boolexecpara=True, \
             
             strgcnfgextnexec=None, \
             listnamevarbcomp=[], \
             listscalvarbcomp=[], \
             listlablvarbcomp=[], \
             listtypevarbcomp=[], \
             listpdfnvarbcomp=[], \
             listgdatvarbcomp=[], \
             
             # parameter name, axis label, tick values and scaling of the input variable changed across PCAT runs
             namexaxivari=None, \
             lablxaxivari=None, \
             tickxaxivari=None, \
             scalxaxivari=None, \
            ):
    
    print('Running PCAT in array mode...')
    
    numbiter = len(dictvarbvari)
    indxiter = np.arange(numbiter) 
    
    cntrcomp = 0
    
    if boolexecpara:
        cntrproc = 0

    listrtag = []
    listpridchld = []
    for k, strgcnfgextn in enumerate(listnamecnfgextn):
        
        if strgcnfgextnexec is not None:
            if strgcnfgextn != strgcnfgextnexec:
                continue
        
        strgcnfg = inspect.stack()[1][3] + '_' + strgcnfgextn
    
        dictvarbtemp = deepcopy(dictvarb)
        for strgvarb, valu in dictvarbvari[strgcnfgextn].items():
            dictvarbtemp[strgvarb] = valu
        dictvarbtemp['strgcnfg'] = strgcnfg
    
        listrtagprev = retr_listrtagprev(strgcnfg, gdat.pathpcat)
        cntrcomp += 1

        if (not forcneww and strgcnfgextnexec is None or forcprev and strgcnfgextnexec is not None) and len(listrtagprev) > 0:
            print('Found at least one previous run with the configuration %s' % strgcnfg)
            print('Skipping...')
            listrtag.append(listrtagprev[-1])
        else:
            if len(listrtagprev) > 0:
                print('Found at least one previous run. But, repeating the run anways...')
            else:
                print('Did not find any previous run.')
            if boolexecpara and strgcnfgextnexec is None:
                cntrproc += 1
                prid = os.fork()
                if prid > 0:
                    listpridchld.append(prid)
                else:
                    print('Forking a child process to run the configuration extension...')
                    rtag = init(**dictvarbtemp)
                    os._exit(0)
            else:
                print('Calling the main PCAT function without forking a child...')
                listrtag.append(init(**dictvarbtemp))
    
    if boolexecpara and strgcnfgextnexec is None:
        for prid in listpridchld:
            os.waitpid(prid, 0)
        if cntrproc > 0:
            print('Exiting before comparion plots because of parallel execution...')
            return
    
    if cntrcomp == 0:
        print('Found no runs...')

    print('Final-processing run outputs...')
    for rtag in listrtag:
        print(rtag)
        proc_finl(rtag=rtag, strgpdfn='post')
        proc_anim(rtag)
    
    strgtimestmp = tdpy.retr_strgtimestmp()
    
    if strgcnfgextnexec is not None or namexaxivari is None: 
        return
    
    print('Making plots to compare the output of different PCAT runs...')
     
    if 'boolmockonly' in dictvarb and dictvarb['boolmockonly']:
        listgdat = retr_listgdat(listrtag, typegdat='init')
    else:
        listgdat = retr_listgdat(listrtag)
    
    numbgdat = len(listgdat)

    for namevarbscal in listgdat[0].listnamevarbscal:
        booltemp = True
        for k in range(1, numbgdat - 1):
            if not namevarbscal in listgdat[k].listnamevarbscal:
                booltemp = False
        if booltemp:
            if namevarbscal in listnamevarbcomp:
                raise Exception('')
            listnamevarbcomp += [namevarbscal]
            listscalvarbcomp += [getattr(listgdat[0], 'scal' + namevarbscal)]
            listlablvarbcomp += [getattr(listgdat[0], 'labl' + namevarbscal + 'totl')]
            listtypevarbcomp += ['pctl']
            listpdfnvarbcomp += ['post']
            listgdatvarbcomp += ['post']
    
    # add others to the variable list
    listnamevarbcomp += ['lliktotl', 'lliktotl', 'infopost', 'bcom', 'lliktotl', 'lliktotl', 'lliktotl', 'levipost']
    listscalvarbcomp += ['self', 'self', 'self', 'self', 'self', 'self', 'self', 'self']
    listlablvarbcomp += ['$\ln P(D|M_{min})$', '$\ln P(D|M_{max})$', '$D_{KL}$', '$\eta_B$', '$\sigma_{P(D|M)}$', r'$\gamma_{P(D|M)}$', \
                                                                                                                    r'$\kappa_{P(D|M)}$', '$\ln P_H(D)$']
    listtypevarbcomp += ['minm', 'maxm', '', '', 'stdv', 'skew', 'kurt', '']
    listpdfnvarbcomp += ['post', 'post', 'post', 'post', 'post', 'post', 'post', 'post']
    listgdatvarbcomp += ['post', 'post', 'post', 'post', 'post', 'post', 'post', 'post']
    
    arrytemp = np.array([len(listnamevarbcomp), len(listscalvarbcomp), len(listlablvarbcomp), len(listtypevarbcomp), len(listpdfnvarbcomp), len(listgdatvarbcomp)])
    if (arrytemp - np.mean(arrytemp) != 0.).all():
        raise Exception('')

    # add log-evidence to the variable list, if prior is also sampled
    booltemp = True
    for k in range(numbgdat):
        if not listgdat[k].checprio:
            booltemp = False
    
    if booltemp:
        listgdatprio = retr_listgdat(listrtag, typegdat='finlprio')
        
        listnamevarbcomp += ['leviprio']
        listscalvarbcomp += ['self']
        listlablvarbcomp += ['$\ln P_{pr}(D)$']
        listtypevarbcomp += ['']
        listpdfnvarbcomp += ['prio']
        listgdatvarbcomp += ['prio']
    
    # time stamp
    strgtimestmp = tdpy.retr_strgtimestmp()
    
    dictoutp = dict()
    liststrgvarbtotl = []
    for (typevarbcomp, pdfnvarbcomp, namevarbcomp) in zip(listtypevarbcomp, listpdfnvarbcomp, listnamevarbcomp):
        strgtemp = typevarbcomp + pdfnvarbcomp + namevarbcomp
        liststrgvarbtotl.append(strgtemp)
        dictoutp[strgtemp] = [[] for k in range(numbiter)]
    
    for k in indxiter:
        for a, strgvarbtotl in enumerate(liststrgvarbtotl):
            if listgdatvarbcomp[a] == 'prio':
                gdattemp = listgdatprio[k]
            else:
                gdattemp = listgdat[k]
            dictoutp[strgvarbtotl][k] = getattr(gdattemp, strgvarbtotl)

    pathbase = '%s/imag/%s_%s/' % (gdat.pathpcat, strgtimestmp, inspect.stack()[1][3])
    cmnd = 'mkdir -p %s' % pathbase 
    os.system(cmnd)
    cmnd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%smrgd.pdf' % pathbase
    for strgvarbtotl, varboutp in dictoutp.items():
        
        figr, axis = plt.subplots(figsize=(6, 6))
        ydat = np.empty(numbiter)
        yerr = np.zeros((2, numbiter))
            
        indxlist = liststrgvarbtotl.index(strgvarbtotl)
        
        if listscalvarbcomp is None:
            scalyaxi = getattr(listgdat[0], 'scal' + listnamevarbcomp[indxlist])
        else:
            scalyaxi = listscalvarbcomp[indxlist]
        
        lablyaxi = listlablvarbcomp[indxlist]
        
        try:
            if listtypevarbcomp[indxlist] == 'pctl':
                trueyaxi = getattr(listgdat[0], 'true' + listnamevarbcomp[indxlist])
            else:
                trueyaxi = getattr(listgdat[0], 'true' + listtypevarbcomp[indxlist] + listnamevarbcomp[indxlist])
        except:
            trueyaxi = None
        
        for k in indxiter:
            
            if isinstance(varboutp[k], list) or isinstance(varboutp[k], np.ndarray) and varboutp[k].ndim > 2:
                raise Exception('')
            elif isinstance(varboutp[k], float):
                ydat[k] = varboutp[k]
            else:
                if listtypevarbcomp[indxlist] != 'pctl':
                    yerr[:, k] = 0.
                if varboutp[k].ndim == 2:
                    if varboutp[k].shape[1] != 1:
                        raise Exception('varboutp format is wrong.')
                    varboutp[k] = varboutp[k][:, 0]
                    if listtypevarbcomp[indxlist] == 'pctl':
                        yerr[:, k] = getattr(listgdat[k], 'errr' + listpdfnvarbcomp[indxlist] + listnamevarbcomp[indxlist])[:, 0]
                else:
                    if listtypevarbcomp[indxlist] == 'pctl':
                        yerr[:, k] = getattr(listgdat[k], 'errr' + listpdfnvarbcomp[indxlist] + listnamevarbcomp[indxlist])
                ydat[k] = varboutp[k][0]
        
        axis.errorbar(indxiter+1., ydat, yerr=yerr, color='b', ls='', markersize=15, marker='o', lw=3)
        indxrtagyerr = np.where((yerr[0, :] > 0.) | (yerr[1, :] > 0.))[0]
        if indxrtagyerr.size > 0:
            temp, listcaps, temp = axis.errorbar(indxiter[indxrtagyerr]+1., ydat[indxrtagyerr], yerr=yerr[:, indxrtagyerr], \
                                                                                color='b', ls='', capsize=15, markersize=15, marker='o', lw=3)
            for caps in listcaps:
                caps.set_markeredgewidth(3)
        
        if trueyaxi is not None:
            axis.axhline(trueyaxi, ls='--', color='g')
        
        if lablxaxivari is None:
            lablxaxivari = getattr(listgdat[0], 'labl' + namexaxivari + 'totl')
        
        if scalxaxivari is None:
            scalxaxivari = getattr(listgdat[0], 'scal' + namexaxivari)
        
        axis.set_xlabel(lablxaxivari)
        axis.set_xticks(indxiter+1.)
        axis.set_xticklabels(tickxaxivari)
        
        axis.set_ylabel(lablyaxi)
        if scalyaxi == 'logt':
            axis.set_yscale('log')
        plt.tight_layout()
        
        pathfull = '%s%s_%s_%s.pdf' % (pathbase, strgtimestmp, inspect.stack()[1][3], liststrgvarbtotl[indxlist])
        print('Writing to %s...' % pathfull)
        plt.savefig(pathfull)
        plt.close(figr)
    
        cmnd += ' %s' % pathfull

    print(cmnd)
    os.system(cmnd)

    print('Making animations...')
    for rtag in listrtag:
        print('Working on %s...' % rtag)
        proc_anim(rtag=rtag)
    
    print('Compiling run plots...')
    cmnd = 'python comp_rtag.py'
    for rtag in listrtag: 
        cmnd += ' %s' % rtag
    os.system(cmnd)

    return listrtag


def retr_rtag(strgcnfg, strgnumbswep):
    
    rtag = strgcnfg + '_' + strgnumbswep
    
    return rtag


class logg(object):
    
    def __init__(self, gdat):
        self.terminal = sys.stdout
        gdat.pathstdo = gdat.pathoutprtag + 'stdo.txt'
        self.log = open(gdat.pathstdo, 'a')
        pathlink = gdat.pathplotrtag + 'stdo.txt'
        os.system('ln -s %s %s' % (gdat.pathstdo, pathlink))
    
    def write(self, strg):
        self.terminal.write(strg)
        self.log.write(strg)  

    def flush(self):
        pass


def worktrac(pathoutprtag, lock, strgpdfn, indxprocwork):
	
    try:
        return work(pathoutprtag, lock, strgpdfn, indxprocwork)
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def opti_hess(gdat, gdatmodi):
    
    gmod = gdat.fitt

    if gmod.numbparaelem > 0:
        cntr = 0
        for l in gmod.indxpopl:
            for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                gdatmodi.indxparastdp[gmod.numbparagenrbase-gmod.numbpopl+cntr] = np.concatenate(gdatmodi.this.indxparagenrfullelem[nameparagenrelem])
                cntr += 1
    
    if gmod.numbparaelem > 0:
        gdatmodi.next.indxelemfull = gdatmodi.this.indxelemfull
        gdatmodi.next.indxparagenrfullelem = gdatmodi.this.indxparagenrfullelem
    else:
        gdatmodi.next.indxparagenrfullelem = None

    gdatmodi.stdpmatr = np.zeros((gdat.numbstdp, gdat.numbstdp)) 
    gdatmodi.hess = np.zeros((gdat.numbstdp, gdat.numbstdp)) 
    deltlpos = np.zeros((3, 3))
    diffpara = np.empty(gdat.numbstdp)
    for k, indxparatemp in enumerate(gdatmodi.indxparastdp):
        if len(indxparatemp) == 0:
            diffpara[k] = 0.
        else:
            diffpara[k] = min(min(np.amin(gdatmodi.this.paragenrunitfull[indxparatemp]) * 0.9, np.amin(1. - gdatmodi.this.paragenrunitfull[indxparatemp]) * 0.9), 1e-5)

    #gdatmodi.this.sampunitsave = np.copy(gdatmodi.this.paragenrunitfull)
    
    #if gmod.numbparaelem > 0:
    #    gdatmodi.dictmodi = [[] for l in gmod.indxpopl]
    #    for l in gmod.indxpopl:
    #        gdatmodi.dictmodi[l] = dict()
    #        gdatmodi.dictmodi[l][gmod.nameparagenrelemampl[l] + 'indv'] = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[gmod.nameparagenrelemampl[l]][l]]
    #        for nameparagenrelem in gmod.namepara.genrelem[l]:
    #            gdatmodi.dictmodi[l]['stdv' + nameparagenrelem + 'indv'] = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l][nameparagenrelem]]
    #if gmod.numbparaelem > 0:
    #    gdatmodi.this.indxparagenrfullelemconc = np.concatenate([gdatmodi.this.indxparagenrfullelem[l]['full'] for l in gmod.indxpopl])
    #if gdat.boolpropcomp:
    #    indxsamptranprop = gdatmodi.this.indxparagenrfullelemconc
    #else:
    #    indxsamptranprop = []
    
    deltlpos[1, 1] = gdatmodi.this.lliktotl
    for indxstdpfrst in gdat.indxstdpprop:
        for indxstdpseco in gdat.indxstdpprop:
            
            if indxstdpfrst > indxstdpseco:
                continue
            
            if indxstdpfrst == indxstdpseco:
                
                #if gmod.numbparaelem > 0:
                #    if k in gdatmodi.this.indxparagenrfullelemconc:
                #        indxtrapmoditemp = k - gmod.indxparagenrfulleleminit
                #        indxpoplmoditemp = np.array([np.amin(np.where(indxtrapmoditemp // gmod.numbparagenrelemcumr == 0))])
                #        numbparapoplinittemp = indxtrapmoditemp - gmod.numbparagenrelemcuml[indxpoplmoditemp[0]]
                #        indxelemmoditemp = [numbparapoplinittemp // gmod.numbparagenrelemsing[indxpoplmoditemp[0]]]
                #        gmod.indxparagenrelemmoditemp = numbparapoplinittemp % gmod.numbparagenrelemsing[indxpoplmoditemp[0]]
                #        nameparagenrelem = gmod.namepara.genrelem[indxpoplmoditemp[0]][gmod.indxparagenrelemmoditemp] 
                #        indxsampampltemp = k - gmod.indxparagenrelemmoditemp + gmod.indxparagenrelemampl[indxpoplmoditemp[0]]
                #        #amplfact = gdatmodi.this.paragenrscalfull[indxsampampltemp] / getattr(gdat, 'minm' + gmod.nameparagenrelemampl[indxpoplmoditemp[0]])
                #        stdv = 1. / np.sqrt(gdatmodi.hess[indxstdpfrst, indxstdpseco])
                #        gdatmodi.stdpmatr[indxstdpfrst, indxstdpseco] += stdv
                #        gdatmodi.dictmodi[indxpoplmoditemp[0]]['stdv' + nameparagenrelem + 'indv'][indxelemmoditemp[0]] = stdv
                #        gdatmodi.dictmodi[indxpoplmoditemp[0]][gmod.nameparagenrelemampl[indxpoplmoditemp[0]] + 'indv'][indxelemmoditemp[0]] = \
                #                                                                                             gdatmodi.this.paragenrscalfull[indxsampampltemp]
                
                if len(gdatmodi.indxparastdp[indxstdpseco]) == 0:
                    continue
                
                for a in range(2):
                    gdatmodi.next.paragenrunitfull = np.copy(gdatmodi.this.paragenrunitfull)
                    if a == 0:
                        gdatmodi.next.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] -= diffpara[indxstdpseco]
                    if a == 1:
                        gdatmodi.next.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] += diffpara[indxstdpseco]
                    
                    gdatmodi.next.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gdatmodi.next.paragenrunitfull, gdatmodi.next.indxparagenrfullelem)
                    
                    proc_samp(gdat, gdatmodi, 'next', 'fitt')
                    if a == 0:
                        deltlpos[0, 1] = gdatmodi.next.lliktotl
                    if a == 1:
                        deltlpos[2, 1] = gdatmodi.next.lliktotl
                
                gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / diffpara[indxstdpseco]**2 * np.fabs(deltlpos[0, 1] + \
                                                                                                        deltlpos[2, 1] - 2. * deltlpos[1, 1])
            else:
                # temp
                continue

                for a in range(4):
                    gdatmodi.this.paragenrunitfull = np.copy(gdatmodi.this.sampunitsave)
                    if a == 0:
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpfrst]] -= diffpara
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] -= diffpara
                    if a == 1:
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpfrst]] += diffpara
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] += diffpara
                    if a == 2:
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpfrst]] -= diffpara
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] += diffpara
                    if a == 3:
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpfrst]] += diffpara
                        gdatmodi.this.paragenrunitfull[gdatmodi.indxparastdp[indxstdpseco]] -= diffpara
                    proc_samp(gdat, gdatmodi, 'this', 'fitt')
                    if a == 0:
                        deltlpos[0, 0] = gdatmodi.this.lpostotl
                    if a == 1:
                        deltlpos[2, 2] = gdatmodi.this.lpostotl
                    if a == 2:
                        deltlpos[1, 2] = gdatmodi.this.lpostotl
                    if a == 3:
                        deltlpos[2, 1] = gdatmodi.this.lpostotl
                gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / diffpara**2 * \
                                                                                (deltlpos[2, 2] + deltlpos[0, 0] - deltlpos[1, 2] - deltlpos[2, 1])
            
            if not np.isfinite(gdatmodi.hess[indxstdpfrst, indxstdpseco]):
                raise Exception('')
            if gdat.booldiagmode and not np.isfinite(gdatmodi.next.paragenrscalfull).all():
                raise Exception('')
            if gdatmodi.hess[indxstdpfrst, indxstdpseco] == 0.:
                raise Exception('')

    gdatmodi.hess[np.where(gdatmodi.hess == 0)] = 10.

    # temp
    #gdatmodi.stdpmatr = np.sqrt(linalg.inv(gdatmodi.hess))
    numbdoffefff = gmod.numbparagenrbase
    if gmod.numbparaelem > 0:
        numbdoffefff += gmod.numbparagenrelem * 10
    gdatmodi.stdpmatr = np.sqrt(1. / gdatmodi.hess) / np.sqrt(numbdoffefff)
    
    if (gdatmodi.stdpmatr == 0).any():
        raise Exception('')
    
    gdatmodi.stdp = gdatmodi.stdpmatr[gdat.indxstdp, gdat.indxstdp]
    

def worksamp(gdat, lock, strgpdfn='post'): 
    
    pathorig = gdat.pathoutprtag + 'stat.txt'
    pathlink = gdat.pathplotrtag + 'stat.txt'
    os.system('ln -s %s %s' % (pathorig, pathlink))
    
    if gdat.numbproc == 1:
        worktrac(gdat.pathoutprtag, lock, strgpdfn, 0)
    else:
        if gdat.typeverb > 0:
            print('Forking the sampler...')

        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(worktrac, gdat.pathoutprtag, lock, strgpdfn)
        pool.map(workpart, gdat.indxproc)

        pool.close()
        pool.join()
    
    gdat.filestat = open(gdat.pathoutprtag + 'stat.txt', 'a')
    gdat.filestat.write('gdatmodi%s written.\n' % strgpdfn)
    gdat.filestat.close()


def work(pathoutprtag, lock, strgpdfn, indxprocwork):
    
    print('Worker #%d' % indxprocwork)
    
    # read the initial global object, gdatinit
    path = pathoutprtag + 'gdatinit'
    gdat = readfile(path) 
    
    gmod = gdat.fitt
    
    # define time functions
    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for this chain
    if gdat.boolseedchan:
        np.random.seed(indxprocwork + 1000)

    # construct a global object for the walker
    gdatmodi = tdpy.gdatstrt()
    gdatmodi.this = tdpy.gdatstrt()
    gdatmodi.next = tdpy.gdatstrt()
    gdatmodi.indxprocwork = indxprocwork
    
    gdatmodi.this = gdat.fitt.this

    # path of gdatmodi
    gdatmodi.pathgdatmodi = gdat.pathoutprtag + 'gdatmodi%04d' % gdatmodi.indxprocwork + gdat.strgpdfn
    
    print('Determining the parameter indices of the fitting model with only the floating parameters...')

    gdatmodi.booldone = False
    gdatmodi.lock = lock
    gdatmodi.indxprocwork = indxprocwork
    
    # plotting factors for scalar variables
    for name in gmod.namepara.scal:
        if name in gmod.nameparagenrbase:
            gmod.indxpara.temp = np.where(gmod.nameparagenrbase == name)[0]
    
    # find the list of variables for which the posterior will be calculated
    if not gdat.boolmockonly:
        
        if gdat.typeverb > 1:
            print('gdatmodi.this.paragenrunitfull')
            print(gdatmodi.this.paragenrunitfull)
            show_paragenrscalfull(gdat, gdatmodi)
        proc_samp(gdat, gdatmodi, 'this', 'fitt')
        
        gdat.liststrgvarbarrysamp = []
        gdat.liststrgvarblistsamp = []
        for strg, valu in gdatmodi.this.__dict__.items():
            if not strg in gdat.liststrgvarbarryswep:
                if isinstance(valu, np.ndarray) or isinstance(valu, float):
                    gdat.liststrgvarbarrysamp.append(strg)
                elif isinstance(valu, list) and strg != 'indxparagenrfullelem' and strg != 'psfnconv' and \
                                                                     strg != 'trueindxelemasscmiss' and strg != 'trueindxelemasschits':
                    gdat.liststrgvarblistsamp.append(strg)
        if gdat.typeverb == 2:
            print('gdat.liststrgvarbarrysamp')
            print(gdat.liststrgvarbarrysamp)
            print('gdat.liststrgvarblistsamp')
            print(gdat.liststrgvarblistsamp)
        
        gdat.liststrgvarblistswep = []
        if gdat.typeverb == 2:
            print('gdat.liststrgvarblistswep')
            print(gdat.liststrgvarblistswep)

        gdat.liststrgvarblist = gdat.liststrgvarblistsamp + gdat.liststrgvarblistswep

        gdatmodi.liststrgchan = gdat.liststrgvarbarryswep + ['paragenrscalbase'] + gmod.namepara.scal
        
        if gdat.typeverb == 2:
            print('gdatmodi.liststrgchan')
            print(gdatmodi.liststrgchan)
    
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
    
    ## sample index
    gdatmodi.cntrswep = 0
   
    if gdat.booldiagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    
    # definitions required for the initial sample
    gdatmodi.this.boolpropfilt = True
    gdatmodi.this.boolpropaccp = True
    
    # dummy definitions required for logs
    gdatmodi.this.indxproptype = np.zeros(1, dtype=int)
    for l in gmod.indxpopl:
        setattr(gdatmodi.this, 'auxiparapop%d' % l, np.zeros(gmod.numbparagenrelemsing[l]))
    gdatmodi.this.lpri = np.zeros(gmod.numblpri)
    gdatmodi.this.lpau = np.zeros(1)
    gdatmodi.this.ltrp = np.zeros(1)
    gdatmodi.this.ljcb = np.zeros(1)
    gdatmodi.this.accpprob = np.zeros(1)
    gdatmodi.this.memoresi = np.zeros(1)
    gdatmodi.this.amplpert = np.zeros(1)
    
    # make sure the first sample derived variables are generated on gdatmodi
    proc_samp(gdat, gdatmodi, 'this', 'fitt')
    
    # log the initial state
    if False and gdat.typeverb > 1:
        tdpy.show_memo(gdatmodi, 'gdatmodi')
    
    for k, name in enumerate(gdat.listnamechro):
        setattr(gdatmodi.this, 'chro' + name, 0.)
    
    gdatmodi.stdp = np.copy(gdat.stdp)
    
    # indices of parameters corresping to each proposal scale
    gdatmodi.indxparastdp = [[] for k in gdat.indxstdp]
    for k in gmod.indxparagenrbase:
        if k < gmod.numbpopl:
            continue
        gdatmodi.indxparastdp[k-gmod.numbpopl] = [k]
    
    workdict = {}
    # list of variable names with type numpy array
    for strgvarb in gdat.liststrgvarbarry:
        valu = getattr(gdatmodi.this, strgvarb)
        if strgvarb in gdat.liststrgvarbarryswep:
            if isinstance(valu, dict):
                shap = [gdat.numbswep, len(valu.keys())]
            elif isinstance(valu, float) or isinstance(valu, bool):
                shap = [gdat.numbswep, 1]
            else:
                shap = [gdat.numbswep] + list(valu.shape)
        else:
            shap = [gdat.numbsamp] + list(valu.shape)
        workdict['list' + gdat.strgpdfn + strgvarb] = np.zeros(shap)
   
    # list of variable names with type list
    for strgvarb in gdat.liststrgvarblist:
        workdict['list' + gdat.strgpdfn + strgvarb] = []
    
    ## saved state of the sample index used for logging progress status
    gdatmodi.percswepsave = -1
   
    # store the initial sample as the best fit sample
    gdatmodi.maxmllikswep = np.sum(gdatmodi.this.llik)
    gdatmodi.indxswepmaxmllik = -1 
    gdatmodi.sampmaxmllik = np.copy(gdatmodi.this.paragenrscalfull)
    
    if gdat.typeverb > 0:
        print('Sampling...')
        print('gdat.stdp')
        for k in gdat.indxstdp:
            print('%04d %s %g' % (k, gdat.namestdp[k], gdat.stdp[k]))

    gdatmodi.this.stdp = np.copy(gdat.stdp)

    gdatmodi.optidone = False 
    
    while gdatmodi.cntrswep < gdat.numbswep:
        
        initchro(gdat, gdatmodi, 'totl')
        
        # Boolean flag to indicate burn-in phase
        gdatmodi.boolburn = gdatmodi.cntrswep < gdat.numbburn
        
        # temp
        if gdat.typeopti == 'hess' and gdatmodi.cntrswep % gdat.numbstdp * 4 == 0 and gdatmodi.cntrswep < gdat.numbburn:
            if gdat.typeverb > 0:
                print('Optimizing proposal scale...')
            opti_hess(gdat, gdatmodi)
            
            if (gdatmodi.stdpmatr[gdat.indxstdp, gdat.indxstdp] < 0.5).all():
                gdatmodi.optidone = True
           
        if gdat.typeopti == 'hess' and gdatmodi.cntrswep == gdat.numbburn:
            path = gdat.pathoutprtag + 'opti.h5'
            if gdat.typeverb > 0:
                print('Writing the estimated covariance matrix to %s...' % path)
            thisfile = h5py.File(path, 'w')
            thisfile.create_dataset('stdp', data=gdatmodi.stdp)
            thisfile.close()
            
            if gdat.makeplot:
                
                xdat = gdat.indxstdp
                ydat = gdatmodi.stdp
                
                pathopti = getattr(gdat, 'path' + gdat.strgpdfn + 'opti')
                path = pathopti + 'stdv%d.pdf' % gdatmodi.indxprocwork
                tdpy.plot_gene(path, xdat, ydat, scalydat='logt', \
                                lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist', limtydat=[np.amin(ydat) / 2., 2. * np.amax(ydat)])
                
                # plot uncertainties of element parameters as a function of amplitude parameter
                if gmod.numbparaelem > 0:
                    for l in gmod.indxpopl:
                        for nameparagenrelem in gmod.namepara.genrelem[l]:
                            path = pathopti + 'stdv' + nameparagenrelem + 'pop%d.pdf' % l
                            xdat = [gdatmodi.dictmodi[l][gmod.nameparagenrelemampl[l] + 'indv'], meanplot]
                            
                            if nameparagenrelem == gmod.nameparagenrelemampl[l]:
                                ydat = [gdatmodi.dictmodi[l]['stdv' + nameparagenrelem + 'indv'], \
                                                                gdatmodi.stdp[getattr(gdat, 'indxstdp' + nameparagenrelem)] / (meanplot / minm)**2.]
                            else:
                                ydat = [gdatmodi.dictmodi[l]['stdv' + nameparagenrelem + 'indv'], \
                                                                gdatmodi.stdp[getattr(gdat, 'indxstdp' + nameparagenrelem)] / (meanplot / minm)**0.5]
                            lablxdat = getattr(gmod.lablpara, gmod.nameparagenrelemampl[l] + 'totl')
                            scalxdat = getattr(gdat, 'scal' + gmod.nameparagenrelemampl[l] + 'plot')
                            limtxdat = np.array(getattr(gdat, 'limt' + gmod.nameparagenrelemampl[l]))
                            tdpy.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                             lablydat=r'$\sigma_{%s}$' % getattr(gmod.lablpara, nameparagenrelem), plottype=['scat', 'lghtline'])
                            #tdpy.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                            #                                 lablydat=r'$\sigma_{%s}$%s' % (getattr(gmod.lablpara, nameparagenrelem), getattr(gmod.lablpara, nameparagenrelem + 'unit')), plottype=['scat', 'lghtline'])
                            
                            tdpy.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                             lablydat=r'$\sigma_{%s}$' % getattr(gmod.lablpara, nameparagenrelem), plottype=['scat', 'lghtline'])


        if gdat.typeverb > 1:
            print('-' * 10)
            print('Sweep %d' % gdatmodi.cntrswep)

        # decide whether to make a frame
        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and \
                                                gdatmodi.indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                and gdat.makeplotfram and gdat.makeplot
        
        # decide whether to make a log
        boollogg = False
        if gdat.typeverb > 0:
            gdatmodi.this.percswep = 5 * int(20. * gdatmodi.cntrswep / gdat.numbswep) 
            if gdatmodi.this.percswep > gdatmodi.percswepsave or thismakefram:
                gdatmodi.percswepsave = gdatmodi.this.percswep
                minmswepintv = max(0, gdatmodi.cntrswep - 1000)
                maxmswepintv = gdatmodi.cntrswep + 1
                if maxmswepintv > minmswepintv:
                    boollogg = True
        
        # propose the next sample
        if gdat.typeverb > 1:        
            print('-----')
            print('thislliktotl')
            print(gdatmodi.this.lliktotl)
            print('thislpostotl')
            print(gdatmodi.this.lpostotl)
            print('Proposing...')
        
        if gdat.boolburntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
            gdatmodi.this.facttmpr = ((gdatmodi.cntrswep + 1.) / gdat.numbburntmpr)**4
            gdatmodi.this.tmprfactstdv = 1. / gdatmodi.this.facttmpr
            #gdatmodi.this.tmprlposelem = -1000. * (1. - gdatmodi.this.facttmpr) * np.concatenate(gdatmodi.this.indxparagenrfullelem['full']).size
            gdatmodi.this.tmprlposelem = 0.
        else:
            gdatmodi.this.tmprfactstdv = 1.
            gdatmodi.this.tmprlposelem = 0. 
        
        # temp -- this can be faster
        for l in gmod.indxpopl:
            setattr(gdatmodi.this, 'auxiparapop%d' % l, np.empty(gmod.numbparagenrelemsing[l]))

        if gdat.typeverb > 1:
            show_paragenrscalfull(gdat, gdatmodi)
        
        # make a proposal
        initchro(gdat, gdatmodi, 'prop')
        prop_stat(gdat, gdatmodi, 'fitt')
        stopchro(gdat, gdatmodi, 'prop')

        if gdat.booldiagmode:
        
            for k in gmod.indxparagenrbase:
                if gmod.scalpara.genrbase[k] == 'logt' and gdatmodi.this.paragenrscalfull[k] < 0.:
                    raise Exception('')

            if not np.isfinite(gdatmodi.next.paragenrscalfull).all():
                raise Exception('')
        
        if gdat.typeverb > 1:
            show_paragenrscalfull(gdat, gdatmodi, strgstat='next')
    
        if (thismakefram or gdat.boolsave[gdatmodi.cntrswep] or boollogg):
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this', 'fitt')
        
        # diagnostics
        if gdat.booldiagmode:
            
            initchro(gdat, gdatmodi, 'diag')
            
            indxsampbadd = np.where((gdatmodi.this.paragenrunitfull[gmod.numbpopl:] > 1.) | (gdatmodi.this.paragenrunitfull[gmod.numbpopl:] < 0.))[0] + 1
            if indxsampbadd.size > 0:
                raise Exception('Unit sample vector went outside [0,1].')
            
            if not np.isfinite(gdatmodi.this.lliktotl):
                raise Exception('Log-likelihood is infinite!')
    
            #indxsampclos = np.where((gdatmodi.this.paragenrscalfull < 0.01) & (gdatmodi.this.paragenrscalfull % 1. != 0.))[0]
            #indxsampclos = list(indxsampclos)
            #for indxparagenrfulltemp in indxsampclos:
            #    for l in gmod.indxpopl:
            #        if not indxparagenrfulltemp in gdatmodi.this.indxparagenrfullelem[l]['full']:
            #            indxsampclos.remove(indxparagenrfulltemp)
            #indxsampclos = np.array(indxsampclos)
            #if indxsampclos.size > 0:
            #    print 'Warning! State is too close to 0!'
            #    print gmod.namepara[indxsampclos]

            #indxsampclos = np.where((gdatmodi.this.paragenrscalfull > 0.99) & (gdatmodi.this.paragenrscalfull % 1. != 0.))[0]
            #indxsampclos = list(indxsampclos)
            #for indxparagenrfulltemp in indxsampclos:
            #    for l in gmod.indxpopl:
            #        if not indxparagenrfulltemp in gdatmodi.this.indxparagenrfullelem[l]['full']:
            #            indxsampclos.remove(indxparagenrfulltemp)
            #indxsampclos = np.array(indxsampclos)
            #if indxsampclos.size > 0:
            #    print 'Warning! State is too close to 1!'
            #    print gmod.namepara[indxsampclos]

            if gdatmodi.cntrswep == 0:
                gdatmodi.this.lliktotlprev = gdatmodi.this.lliktotl
            
            lliktotldiff = gdatmodi.this.lliktotl - gdatmodi.this.lliktotlprev

            if gdatmodi.this.lliktotl - gdatmodi.this.lliktotlprev < -10.:
                raise Exception('loglikelihood drop is very unlikely!')
            gdatmodi.this.lliktotlprev = gdatmodi.this.lliktotl
       
            for strgstat in ['this', 'next']:
                for strgvarb in ['paragenrscalfull', 'paragenrunitfull']:
                    varb = getattr(getattr(gdatmodi, strgstat), strgvarb)
                    if not np.isfinite(varb).all():
                        raise Exception('Sample vector is not finite.')
            
            if gmod.numbparaelem > 0:
                if gmod.boolelemsbrtdfncanyy:
                    thissbrtdfnc = getattr(gdatmodi.this, 'sbrtdfnc')
                    frac = np.amin(thissbrtdfnc) / np.mean(thissbrtdfnc)
                    cntppntschec = retr_cntp(gdat, thissbrtdfnc)
                    if np.amin(cntppntschec) < -0.1 and frac < -1e-3:
                        raise Exception('thissbrtdfnc went negative by %.3g percent.' % (100. * frac))
                    
            # check the population index
            if (gdatmodi.this.cntpmodl <= 0.).any() or not (np.isfinite(gdatmodi.this.cntpmodl)).all():
                raise Exception('Current flux model is not positive')

            if gmod.numbparaelem > 0:
                for l in gmod.indxpopl:
                    if gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]] != len(gdatmodi.this.indxelemfull[l]):
                        raise Exception('Number of elements is inconsistent with the element index list.')

                    if gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]] != len(gdatmodi.this.indxelemfull[l]):
                        raise Exception('Number of elements is inconsistent across data structures.')
                    
                    for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                        if gmod.listscalparagenrelem[l][k] == 'gaus' or gmod.listscalparagenrelem[l][k] == 'igam' \
                                                                                            or gmod.listscalparagenrelem[l][k] == 'expo':
                            continue
                        comp = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l][nameparagenrelem]]
                        minm = getattr(gdat.fitt.minmpara, nameparagenrelem)
                        maxm = getattr(gdat.fitt.maxmpara, nameparagenrelem)
                        indxtemp = np.where((comp < minm) | (comp > maxm))[0]
                        if indxtemp.size > 0:
                            raise Exception('A component of an element went outside the prior range.')
        
            stopchro(gdat, gdatmodi, 'diag')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            initchro(gdat, gdatmodi, 'save')
        
            if gdat.savestat:
                
                if gdat.namesavestat is not None:
                    strgcnfg = gdat.namesavestat
                else:
                    strgcnfg = gdat.strgcnfg
                path = gdat.pathoutp + 'stat_' + strgcnfg + '.h5'
                
                booltemp = False
                if os.path.isfile(path) and gdatmodi.indxprocwork == 0:
                    thisfilechec = h5py.File(path, 'r')
                    if thisfilechec['lliktotl'][...] > gdatmodi.this.lliktotl:
                        if gdat.typeverb > 0:
                            print('Not saving the state to %s because loglikelihood is lower...' % path)
                            print('Likelihood in the file:')
                            print(thisfilechec['lliktotl'][...])
                    else:
                        booltemp = True
                    thisfilechec.close()
                else:
                    booltemp = True
                if gdat.forcsavestat:
                    booltemp = True
                if booltemp:
                    if gdatmodi.indxprocwork > 0:
                        continue
                    if gdat.typeverb > 0:
                        print('Saving the state to %s...' % path)
        
                    thisfile = h5py.File(path, 'w')
                    thisfile.create_dataset('lliktotl', data=gdatmodi.this.lliktotl)
                    for gmod.nameparagenrbase in gmod.nameparagenrbase:
                        valu = gdatmodi.this.paragenrscalfull[gmod.indxparagenrbase]
                        thisfile.create_dataset(gmod.nameparagenrbase, data=valu)
                    if gmod.numbparaelem > 0:
                        for l in gmod.indxpopl:
                            for nameparagenrelem in gmod.namepara.genrelem[l]:
                                comp = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrfullelem[l][nameparagenrelem]]
                                for k in np.arange(comp.size):
                                    name = nameparagenrelem + 'pop%d%04d' % (l, k)
                                    thisfile.create_dataset(name, data=comp[k])
                    thisfile.close()
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]
            
            # fill the sample lists
            for strgvarb in gdat.liststrgvarbarrysamp:
                valu = getattr(gdatmodi.this, strgvarb)
                workdict['list' + gdat.strgpdfn + strgvarb][indxsampsave, ...] = valu
            for strgvarb in gdat.liststrgvarblistsamp:
                workdict['list' + gdat.strgpdfn + strgvarb].append(deepcopy(getattr(gdatmodi.this, strgvarb)))
            stopchro(gdat, gdatmodi, 'save')

        # plot the current sample
        if thismakefram:
            
            initchro(gdat, gdatmodi, 'plot')
            
            writfile(gdatmodi, gdatmodi.pathgdatmodi) 

            if gdat.typeverb > 0:
                print('Process %d is in queue for making a frame.' % gdatmodi.indxprocwork)
            
            if gdat.numbproc > 1:
                gdatmodi.lock.acquire()
            
            if gdat.typeverb > 0:
                print('Process %d started making a frame.' % gdatmodi.indxprocwork)
            
            plot_samp(gdat, gdatmodi, 'this', 'fitt', 'fram')
            
            if gdat.typeverb > 0:
                print('Process %d finished making a frame.' % gdatmodi.indxprocwork)
        
            if gdat.numbproc > 1:
                gdatmodi.lock.release()
        
            stopchro(gdat, gdatmodi,  'plot')
    
        # determine the acceptance probability
        if gdatmodi.this.boolpropfilt:
            
            initchro(gdat, gdatmodi, 'proc')
            proc_samp(gdat, gdatmodi, 'next', 'fitt')
            stopchro(gdat, gdatmodi, 'proc')
        
            calc_probprop(gdat, gdatmodi)
            
            if gdat.booldiagmode:
                if not gdatmodi.this.indxproptype > 2 and gdatmodi.this.ljcb != 0.:
                    raise Exception('log Jacobian can only be be nonzero when a split or merge is proposed.')
                if not gdatmodi.this.indxproptype > 2 and gdatmodi.this.ltrp != 0.:
                    raise Exception('log ratio proposal probability can only be be nonzero when a split or merge is proposed.')
           
            # evaluate the acceptance probability
            gdatmodi.this.deltlpostotl = gdatmodi.next.lpostotl - gdatmodi.this.lpostotl
            gdatmodi.this.accplprb = gdatmodi.this.deltlpostotl + gdatmodi.this.tmprlposelem - gdatmodi.this.lpau + gdatmodi.this.ltrp + gdatmodi.this.ljcb
            gdatmodi.this.accpprob[0] = np.exp(gdatmodi.this.accplprb)
            if gdat.typeverb > 1:
                print('gdatmodi.this.lpritotl')
                print(gdatmodi.this.lpritotl)
                print('gdatmodi.next.lpritotl')
                print(gdatmodi.next.lpritotl)
                print('gdatmodi.this.lliktotl')
                print(gdatmodi.this.lliktotl)
                print('gdatmodi.next.lliktotl')
                print(gdatmodi.next.lliktotl)
                print('gdatmodi.this.lpostotl')
                print(gdatmodi.this.lpostotl)
                print('gdatmodi.next.lpostotl')
                print(gdatmodi.next.lpostotl)
                
                print('gdatmodi.this.deltlpostotl')
                print(gdatmodi.this.deltlpostotl)
                print('gdatmodi.this.tmprlposelem')
                print(gdatmodi.this.tmprlposelem)
                print('gdatmodi.this.lpau')
                print(gdatmodi.this.lpau)
                print('gdatmodi.this.ltrp')
                print(gdatmodi.this.ltrp)
                print('gdatmodi.this.ljcb')
                print(gdatmodi.this.ljcb)
            
                print('gdatmodi.this.accplprb')
                print(gdatmodi.this.accplprb)
        else:
            gdatmodi.this.accpprob[0] = 0.
    
        # accept or reject the proposal
        booltemp = gdatmodi.this.accpprob[0] >= np.random.rand()
        
        if gdat.booldiagmode:
            if gdatmodi.this.indxproptype == 0:
                if gdat.boolsqzeprop and not booltemp:
                    raise Exception('')

        if booltemp:
            if gdat.typeverb > 1:
                print('Accepted.')
            
            # update the current state
            updt_stat(gdat, gdatmodi)

            # check if the accepted sample has maximal likelihood
            if gdatmodi.this.lliktotl > gdatmodi.maxmllikswep:
                gdatmodi.maxmllikswep = gdatmodi.this.lliktotl
                gdatmodi.indxswepmaxmllik = gdatmodi.cntrswep
                gdatmodi.sampmaxmllik = np.copy(gdatmodi.this.paragenrscalfull)
            
            # register the sample as accepted
            gdatmodi.this.boolpropaccp = True

        # reject the sample
        else:

            if gdat.typeverb > 1:
                print('Rejected.')

            gdatmodi.this.boolpropaccp = False
            
        ## variables to be saved for each sweep
        for strg in gdat.liststrgvarbarryswep:
            workdict['list' + gdat.strgpdfn + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi.this, strg)
        
        workdict['list' + gdat.strgpdfn + 'accpprob'][gdatmodi.cntrswep, 0] = gdatmodi.this.accpprob[0]
        
        # log the progress
        if boollogg:
            
            print('--------------')
            print('Sweep number %d' % gdatmodi.cntrswep)
            print('%3d%% completed.' % gdatmodi.this.percswep)
            print('%30s %50s %10s' % ('Prop', 'Accp rate', 'Scale'))
            
            indxswepintv = np.arange(minmswepintv, maxmswepintv)
            for k in gdat.indxproptype:
                indxswepprop = indxswepintv[np.where(workdict['list' + gdat.strgpdfn + 'indxproptype'][indxswepintv, 0] == k)]
                boolproptype = workdict['list' + gdat.strgpdfn + 'indxproptype'][indxswepintv, 0] == k
                boolaccp = workdict['list' + gdat.strgpdfn + 'boolpropaccp'][indxswepintv, 0] == 1
                numbaccp = np.where(boolaccp & boolproptype)[0].size
                numbtotl = np.where(boolproptype)[0].size
                if numbtotl > 0:
                    percaccp = 100. * numbaccp / float(numbtotl)
                else:
                    percaccp = 0.
                if k in gdat.indxstdp:
                    strgstdp = '%.3g' % gdat.stdp[k]
                else:
                    strgstdp = ''
                print('%30s %50s' % (gdat.lablproptype[k], 'acceptance rate: %3d%% (%5d out of %5d)' % (percaccp, numbaccp, numbtotl)))
                
            if gdat.boolburntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
                print('Tempered burn-in')
                print('gdatmodi.this.facttmpr')
                print(gdatmodi.this.facttmpr)
            print 
            numbpara = gmod.numbparagenrbase
            if gmod.numbparaelem > 0:
                for l in gmod.indxpopl:
                    numbpara += gdatmodi.this.indxparagenrfullelem[l]['full'].size
            if gmod.numbparaelem > 0:
                print('Number of elements:')
                for l in gmod.indxpopl:
                    print(gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int))
            print('Current number of parameters:')
            print(numbpara)
            print('gdatmodi.this.numbdoff')
            print(gdatmodi.this.numbdoff)
            for attr, valu in gdatmodi.__dict__.items():
                if isinstance(valu, np.ndarray):
                    if 8 * valu.size * gdat.numbsamptotl > 1e9:
                        print('Warning! %s has total length %d and size %s' % (attr, valu.size * gdat.numbsamptotl, \
                                                                                        tdpy.retr_strgmemo(8 * valu.size * gdat.numbsamptotl)))
            if gmod.numbparaelem > 0:
                if gmod.typemodltran == 'pois':
                    print('Mean number of elements:')
                    print(gdatmodi.this.paragenrscalfull[gmod.indxpara.meanelem])
                for l in gmod.indxpopl:
                    if gmod.nameparagenrelemampl[l] == 'flux' and gmod.typeprioflux[l] == 'powr' or gmod.nameparagenrelemampl[l] != 'flux':
                        print('Log-slope of the amplitude parameter distribution, population %d:' % l)
                        indxparagenrbase = getattr(gmod.indxpara, 'slopprio' + gmod.nameparagenrelemampl[l] + 'pop%d' % l)
                        print(gdatmodi.this.paragenrscalfull[indxparagenrbase])
                    else:
                        print('Flux distribution break:')
                        print(gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'brek' + gmod.nameparagenrelemampl[l] + 'pop%d' % l)])
                        print('Flux distribution lower slope:')
                        print(gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'sloplowr' + gmod.nameparagenrelemampl[l] + 'pop%d' % l)])
                        print('Flux distribution upper slope:')
                        print(gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'slopuppr' + gmod.nameparagenrelemampl[l] + 'pop%d' % l)])
            print('Backgrounds')
            print(gdatmodi.this.paragenrscalfull[gmod.indxpara.bacp])
            if gmod.numbparaelem > 0:
                print('Log-prior penalization term: ')
                print(gdatmodi.this.lpripena)
                print('Completeness')
                for q in gdat.indxrefr:
                    if gdat.refr.numbelem[q] == 0:
                        continue
                    l = gdat.refr.indxpoplfittassc[q]
                    print('Reference Population %d, Fitting Population %d' % (q, l))
                    #print('Total:')
                    #print(getattr(gdatmodi.this, 'cmpl' + namevarb))
                    print('Binned in significance feature:')
                    print(getattr(gdatmodi.this, 'cmpl' + gdat.refr.namepara.elemsign[q] + 'pop%d' % q))
                print('False discovery rate')
                for l in gmod.indxpopl:
                    if gdat.fitt.this.numbelem[l] == 0:
                        continue
                    q = gmod.indxpoplrefrassc[l]
                    print('Fitting population %d, Reference Population %d' % (l, q))
                    #print('Total:')
                    #print(getattr(gdatmodi.this, 'fdis' + namevarb))
                    print('Binned in significance feature:')
                    print(getattr(gdatmodi.this, 'fdis' + gdat.fitt.namepara.elemsign[l] + 'pop%d' % l))
    
            print('gdatmodi.this.lliktotl')
            print(gdatmodi.this.lliktotl)
            print('Chi2 per degree of freedom')
            print(gdatmodi.this.chi2doff)
        
        # save the execution time for the sweep
        stopchro(gdat, gdatmodi, 'totl')
        
        if boollogg:
            print('Chronometers: ')
            for k, name in enumerate(gdat.listnamechro):
                #for name, valu in gdat.indxchro.items():
                    #if valu == k:
                thischro = getattr(gdatmodi.this, 'chro' + name)
                print('%s: %.3g msec' % (name, thischro * 1e3))
                booltemp = False
                for l in gmod.indxpopl:
                    if gmod.typeelemspateval[l] == 'locl' and gmod.maxmpara.numbelem[l] > 0:
                        booltemp = True
                if name == 'llik' and gdat.numbpixl > 1 and gmod.numbparaelem > 0 and booltemp:
                    print('%.3g per pixel' % (thischro * 1e3 / np.amin(gdat.numbpixlprox)))
            print 

        if gdat.typeverb > 1:
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
        
        # update the sweep counter
        gdatmodi.cntrswep += 1
        
    for strgvarb in gdat.liststrgvarbarry + gdat.liststrgvarblistsamp:
        valu = workdict['list' + gdat.strgpdfn + strgvarb]
        setattr(gdatmodi, 'list' + gdat.strgpdfn + strgvarb, valu)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    delattr(gdatmodi, 'lock')
    
    gdatmodi.booldone = True

    writfile(gdatmodi, gdatmodi.pathgdatmodi) 

