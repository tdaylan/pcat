# plotting
import matplotlib as mpl
mpl.use('agg')
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

import aspendos
import chalcedon

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
            meanener = gdat.bctrpara.enerplot
        else:
            meanener = gdat.bctrpara.ener

        if gmod.spectype == 'gaus':
            spec = 1. / edis[None, :] / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.bctrpara.ener[:, None] - elin[None, :]) / edis[None, :])**2)
        if gmod.spectype == 'voig':
            args = (gdat.bctrpara.ener[:, None] + 1j * gamm[None, :]) / np.sqrt(2.) / sigm[None, :]
            spec = 1. / sigm[None, :] / np.sqrt(2. * pi) * flux[None, :] * real(scipy.special.wofz(args))
        if gmod.spectype == 'edis':
            edis = edisintp(elin)[None, :]
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.bctrpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'pvoi':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.bctrpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'lore':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.bctrpara.ener[:, None] - elin[None, :]) / edis)**2)
        if gmod.spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :])
        if gmod.spectype == 'colr':
            if plot:
                spec = np.zeros((gdat.numbenerplot, flux.size))
            else:
                spec = np.empty((gdat.numbener, flux.size))
                for i in gdat.indxener:
                    if i < gdat.indxenerpivt:
                        spec[i, :] = flux * (gdat.bctrpara.ener[i] / gdat.enerpivt)**(-sindcolr[i])
                    elif i == gdat.indxenerpivt:
                        spec[i, :] =  flux
                    else:
                        spec[i, :] = flux * (gdat.bctrpara.ener[i] / gdat.enerpivt)**(-sindcolr[i-1])
        if gmod.spectype == 'curv':
            spec = flux[None, :] * meanener[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if gmod.spectype == 'expc':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :]) * np.exp(-(meanener - gdat.enerpivt)[:, None] / expc[None, :])
    
    return spec


### find the surface brightness due to one point source
def retr_sbrtpnts(gdat, xpos, ypos, spec, psfnintp, indxpixlelem):
    
    # calculate the distance to all pixels from each point source
    dist = retr_angldistunit(gdat, xpos, ypos, indxpixlelem)
    
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

    wdth = np.zeros((gdat.numbener, gdat.numbdqlt))
    for i in gdat.indxener:
        for m in gdat.indxdqlt:
            psfntemp = psfn[i, :, m]
            indxanglgood = np.argsort(psfntemp)
            intpwdth = max(frac * np.amax(psfntemp), np.amin(psfntemp))
            if intpwdth >= np.amin(psfntemp[indxanglgood]) and intpwdth <= np.amax(psfntemp[indxanglgood]):
                wdthtemp = sp.interpolate.interp1d(psfntemp[indxanglgood], gdat.blimpara.angl[indxanglgood], fill_value='extrapolate')(intpwdth)
            else:
                wdthtemp = 0.
            wdth[i, m] = wdthtemp
                        
    return wdth


def samp_xposyposfromtmpl(gdat, probtmpl):
    
    indxpixldraw = np.random.choice(gdat.indxpixl, p=probtmpl)
    xpos = gdat.xposgrid[indxpixldraw] + randn(gdat.sizepixl)
    ypos = gdat.yposgrid[indxpixldraw] + randn(gdat.sizepixl)
    
    return xpos, ypos


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
        maxmparagenrscalbase = gmod.maxmparagenrscalbase[thisindxparagenrbase]

        if scalparagenrbase == 'self':
            paragenrscalbaseunit = cdfn_self(paragenrscalbase, listminmparagenrscalbase, maxmparagenrscalbase)
        elif scalparagenrbase == 'logt':
            paragenrscalbaseunit = cdfn_logt(paragenrscalbase, listminmparagenrscalbase, maxmparagenrscalbase)

        elif scalparagenrbase == 'atan':
            gmod.listmaxmparagenrscalbase = gmod.listmaxmparagenrscalbase[thisindxparagenrbase]
            paragenrscalbaseunit = cdfn_atan(paragenrscalbase, listminmparagenrscalbase, gmod.listmaxmparagenrscalbase)
    
    elif scalparagenrbase == 'gaus' or scalparagenrbase == 'eerr':
        gmod.meanpara.genrbasescal = gmod.meanpara.genrbasescal[thisindxparagenrbase]
        gmod.stdvpara.genrbasescal = gmod.stdvpara.genrbasescal[thisindxparagenrbase]
        if scalparagenrbase == 'eerr':
            gmod.cdfnlistminmparagenrscalbaseunit = gmod.cdfnlistminmparagenrscalbaseunit[thisindxparagenrbase]
            gmod.listparagenrscalbaseunitdiff = gmod.listparagenrscalbaseunitdiff[thisindxparagenrbase]
            paragenrscalbaseunit = cdfn_eerr(paragenrscalbase, gmod.meanpara.genrbasescal, gmod.stdvpara.genrbasescal, \
                                                                            gmod.cdfnlistminmparagenrscalbaseunit, gmod.listparagenrscalbaseunitdiff)
        else:
            paragenrscalbaseunit = cdfn_gaus(paragenrscalbase, gmod.meanpara.genrbasescal, gmod.stdvpara.genrbasescal)

    elif scalparagenrbase == 'pois':
        paragenrscalbaseunit = paragenrscalbase
    
    if gdat.booldiag:
        if paragenrscalbaseunit == 0:
            print('Warning. CDF is zero.')

    return paragenrscalbaseunit


def icdf_paragenrscalfull(gdat, strgmodl, paragenrunitfull, indxparagenrfullelem):
    
    gmod = getattr(gdat, strgmodl)

    # tobechanged
    # temp -- change zeros to empty
    paragenrscalfull = np.empty_like(paragenrunitfull)
    for scaltype in gdat.listscaltype:
        listindxparagenrbasescal = gmod.indxpara.genrbasescal[scaltype]
        if len(listindxparagenrbasescal) == 0:
            continue
        paragenrscalfull[listindxparagenrbasescal] = icdf_paragenrscalbase(gdat, strgmodl, paragenrunitfull[listindxparagenrbasescal], \
                                                                                                            scaltype, listindxparagenrbasescal)
    
        if gdat.booldiag and not np.isfinite(paragenrscalfull).all():
            print('')
            print('')
            print('')
            raise Exception('')

    if indxparagenrfullelem is not None:
        for l in gmod.indxpopl:
            for g in gmod.indxparagenrelemsing[l]:
                indxparagenrfulltemp = indxparagenrfullelem[l][gmod.namepara.genrelem[l][g]]
                if indxparagenrfulltemp.size == 0:
                    continue
                
                paragenrscalfull[indxparagenrfulltemp] = icdf_trap(gdat, strgmodl, paragenrunitfull[indxparagenrfulltemp], paragenrscalfull, \
                                                                                        gmod.scalpara.genrelem[l][g], gmod.namepara.genrelem[l][g], l)
    
                if gdat.booldiag:
                    if not np.isfinite(paragenrscalfull[indxparagenrfulltemp]).all():
                        print('')
                        print('')
                        print('')
                        raise Exception('')

    if gdat.booldiag and not np.isfinite(paragenrscalfull).all():
        print('')
        print('')
        print('')
        print('paragenrscalfull')
        print(paragenrscalfull)
        raise Exception('')
    
    return paragenrscalfull

    
def icdf_paragenrscalbase(gdat, strgmodl, paragenrunitbase, scaltype, indxparagenrbasescal):
    
    gmod = getattr(gdat, strgmodl)
    
    if scaltype == 'self' or scaltype == 'logt' or scaltype == 'atan':
        minmparagenrscalbase = gmod.minmpara.genrbase[indxparagenrbasescal]
        maxmparagenrscalbase = gmod.maxmpara.genrbase[indxparagenrbasescal]

    if scaltype == 'self':
        paragenrscalbase = tdpy.icdf_self(paragenrunitbase, minmparagenrscalbase, maxmparagenrscalbase)
    elif scaltype == 'logt':
        paragenrscalbase = tdpy.icdf_logt(paragenrunitbase, minmparagenrscalbase, maxmparagenrscalbase)
    elif scaltype == 'atan':
        listmaxmparagenrscalbase = gmod.listmaxmparagenrscalbase[indxparagenrbasescal]
        paragenrscalbase = tdpy.icdf_atan(paragenrunitbase, minmparagenrscalbase, listmaxmparagenrscalbase)
    elif scaltype == 'gaus' or scaltype == 'eerr':
        listmeanparagenrscalbase = gmod.meanpara.genrbase[indxparagenrbasescal]
        liststdvparagenrscalbase = gmod.stdvpara.genrbase[indxparagenrbasescal]
        if scaltype == 'eerr':
            cdfnminmparagenrscalbaseunit = gmod.cdfnminmparagenrscalbaseunit[indxparagenrbasescal]
            listparagenrscalbaseunitdiff = gmod.listparagenrscalbaseunitdiff[indxparagenrbasescal]
            paragenrscalbase = tdpy.icdf_eerr(paragenrunitbase, listmeanparagenrscalbase, liststdvparagenrscalbase, \
                                                                cdfnminmparagenrscalbaseunit, listparagenrscalbaseunitdiff)
        else:
            paragenrscalbase = tdpy.icdf_gaus(paragenrunitbase, listmeanparagenrscalbase, liststdvparagenrscalbase)
    elif scaltype == 'drct':
        paragenrscalbase = paragenrunitbase
    
    if gdat.booldiag:
        if not np.isfinite(paragenrscalbase).all():
            print('')
            print('')
            print('')
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
        if gdat.booldiag:
            if not np.isfinite(slop):
                print('')
                print('')
                print('')
                raise Exception('')
            if maxm < minm:
                print('')
                print('')
                print('')
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
        minm = getattr(gmod.minmpara, nameparagenrelem)
        maxm = getattr(gmod.maxmpara, nameparagenrelem)
        icdf = tdpy.icdf_self(cdfn, minm, maxm)
    
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
        slop = paragenrscalfull[getattr(gmod.indxpara, 'slopprio%s%d' % (nameparagenrelem, l))[l]]
        cutf = getattr(gdat, 'cutf' + nameparagenrelem)
        icdf = tdpy.icdf_igam(cdfn, slop, cutf)
    
    if scalcomp == 'gaus':
        distmean = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'distmean')[l]]
        diststdv = paragenrscalfull[getattr(gmod.indxpara, nameparagenrelem + 'diststdv')[l]]
        icdf = tdpy.icdf_gaus(cdfn, distmean, diststdv)
    
    if gdat.booldiag:
        if not np.isfinite(icdf).all():
            print('')
            print('')
            print('')
            print('icdf')
            print(icdf)
            raise Exception('')

    return icdf


def cdfn_trap(gdat, gdatmodi, strgmodl, icdf, indxpoplthis):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    
    scaltemp = gmod.scalpara.genrelem[indxpoplthis]
    
    cdfn = np.empty_like(icdf)
    for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[indxpoplthis]):
        
        if scaltemp[k] == 'self' or scaltemp[k] == 'dexp' or scaltemp[k] == 'expo' \
                                                                        or scaltemp[k] == 'powr' or scaltemp[k] == 'dpowslopbrek':
            minm = getattr(gdat.fitt.minmpara, nameparagenrelem)
            if scaltemp[k] == 'powr':
                maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                slop = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'slop')[indxpoplthis]]
                cdfn[k] = cdfn_powr(icdf[k], minm, maxm, slop)
            elif scaltemp[k] == 'dpowslopbrek':
                maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                brek = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'distbrek')[indxpoplthis]]
                sloplowr = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'sloplowr')[indxpoplthis]]
                slopuppr = gdatobjt.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'slopuppr')[indxpoplthis]]
                cdfn[k] = cdfn_dpow(icdf[k], minm, maxm, brek, sloplowr, slopuppr)
            else:
                maxm = getattr(gdat.fitt.maxmpara, nameparagenrelem)
                cdfn[k] = cdfn_self(icdf[k], minm, maxm)
        if scaltemp[k] == 'lnormeanstdv':
            distmean = gdatmodi.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'diststdv')[indxpoplthis]]
            cdfn[k] = cdfn_lnor(icdf[k], distmean, slop)
        if scaltemp[k] == 'igam':
            slop = gdatmodi.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'slop')[indxpoplthis]]
            cutf = getattr(gdat, 'cutf' + nameparagenrelem)
            cdfn[k] = cdfn_igam(icdf[k], slop, cutf)
        if scaltemp[k] == 'gaus':
            distmean = gdatmodi.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'diststdv')[indxpoplthis]]
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
        gdatmodi.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gdatmodi.this.indxelemfull, 'fitt')


def initcompfromstat(gdat, gdatmodi, namerefr):
    
    for l in gmod.indxpopl:
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
            minm = getattr(gdat.fitt.minmpara, nameparagenrelem)
            maxm = getattr(gdat.fitt.maxmpara, nameparagenrelem)
            try:
                comp = getattr(gdat, namerefr + nameparagenrelem)[l][0, :]
                if gmod.scalpara.genrelem[l][g] == 'self' or gmod.scalpara.genrelem[l][g] == 'logt':
                    if gmod.scalpara.genrelem[l][g] == 'self':
                        compunit = cdfn_self(comp, minm, maxm)
                    if gmod.scalpara.genrelem[l][g] == 'logt':
                        compunit = cdfn_logt(comp, minm, maxm)
                if gmod.scalpara.genrelem[l][g] == 'expo':
                    scal = getattr(gdat.fitt, 'gangdistsexp')
                    maxm = getattr(gdat.fitt.maxm, nameparagenrelem)
                    compunit = cdfn_expo(icdf, maxm, scal)
                if gmod.scalpara.genrelem[l][g] == 'powr' or gmod.scalpara.genrelem[l][g] == 'igam':
                    slop = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'slop')[l]]
                    if gmod.scalpara.genrelem[l][g] == 'powr':
                        compunit = cdfn_powr(comp, minm, maxm, slop)
                    if gmod.scalpara.genrelem[l][g] == 'igam':
                        cutf = getattr(gdat, 'cutf' + nameparagenrelem)
                        compunit = cdfn_igam(comp, slop, cutf)
                if gmod.scalpara.genrelem[l][g] == 'dpowslopbrek':
                    brek = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'distbrek')[l]]
                    sloplowr = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'sloplowr')[l]]
                    slopuppr = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'slopuppr')[l]]
                    compunit = cdfn_powr(comp, minm, maxm, brek, sloplowr, slopuppr)
                if gmod.scalpara.genrelem[l][g] == 'gaus':
                    distmean = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'distmean')[l]]
                    diststdv = gdatmodi.this.paragenrscalfull[getattr(gdat.fitt.indxpara, 'genrbase' + nameparagenrelem + 'diststdv')[l]]
                    compunit = cdfn_gaus(comp, distmean, diststdv)
            except:
                if gdat.typeverb > 0:
                    print('Initialization from the reference catalog failed for %s. Sampling randomly...' % nameparagenrelem)
                compunit = np.random.rand(gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]].astype(int))
            gdatmodi.this.paragenrunitfull[gdatmodi.this.indxparagenrelemfull[l][nameparagenrelem]] = compunit


### find the set of pixels in proximity to a position on the map
def retr_indxpixlelemconc(gdat, strgmodl, dictelem, l):
    
    gmod = getattr(gdat, strgmodl)
    
    xpos = dictelem[l]['xpos']
    ypos = dictelem[l]['ypos']
    varbampl = dictelem[l][gmod.nameparagenrelemampl[l]]
    
    if gmod.typeelemspateval[l] == 'locl':
        listindxpixlelem = [[] for k in range(xpos.size)]
        for k in range(xpos.size):
            indxpixlpnts = retr_indxpixl(gdat, ypos[k], xpos[k])
            
            indxfluxproxtemp = np.digitize(varbampl[k], gdat.blimpara.prox)
            if indxfluxproxtemp > 0:
                indxfluxproxtemp -= 1
            if indxfluxproxtemp == gdat.blimpara.prox.size - 1:
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
def retr_angldistunit(gdat, xpos, ypos, indxpixlelem, retranglcosi=False):
   
    if gdat.typepixl == 'heal':
        xdat, ydat, zaxi = retr_unit(xpos, ypos)
        anglcosi = gdat.xdatgrid[indxpixlelem] * xdat + gdat.ydatgrid[indxpixlelem] * ydat + gdat.zaxigrid[indxpixlelem] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = np.arccos(anglcosi)
            return angldist
    
    else:
        angldist = np.sqrt((xpos - gdat.xposgrid[indxpixlelem])**2 + (ypos - gdat.yposgrid[indxpixlelem])**2)
        
        return angldist
    

### find the pixel index of a point on the map
def retr_indxpixl(gdat, ypos, xpos):

    if gdat.typepixl == 'heal':
        indxpixl = gdat.pixlcnvt[hp.ang2pix(gdat.numbsideheal, np.pi / 2. - ypos, xpos)]
        if gdat.booldiag:
            if (indxpixl == -1).any():  
                print('')
                print('')
                print('')
                raise Exception('pixlcnvt went negative!')

    if gdat.typepixl == 'cart':
        indxlgcr = np.floor(gdat.numbsidecart * (xpos - gdat.minmxposdata) / 2. / gdat.maxmgangdata).astype(int)
        indxbgcr = np.floor(gdat.numbsidecart * (ypos - gdat.minmyposdata) / 2. / gdat.maxmgangdata).astype(int)

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
        path = gdat.pathinit + nameinte + strgplot + '.%s' % gdat.typefileplot
    elif strgstat == 'pdfn' or strgstat == 'mlik':
        path = gdat.pathplotcnfg + strgpdfn + '/finl/' + nameinte + strgstat + strgplot + '.%s' % gdat.typefileplot
    elif strgstat == 'this':
        path = gdat.pathplotcnfg + strgpdfn + '/fram/' + nameinte + strgstat + strgplot + '_swep%09d.%s' % (gdatmodi.cntrswep, gdat.typefileplot)
    
    return path


### determine the marker size
def retr_mrkrsize(gdat, strgmodl, compampl, nameparagenrelemampl):
    
    gmod = getattr(gdat, strgmodl)
    minm = getattr(gdat.minmpara, nameparagenrelemampl) 
    maxm = getattr(gdat.maxmpara, nameparagenrelemampl)
    mrkrsize = (np.sqrt(compampl) - np.sqrt(minm)) / (np.sqrt(maxm) - np.sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


## experiment specific
def retr_psfphubb(gdat, gmod):

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
    
    fermscal = np.zeros((gdat.numbdqlt, numbpsfpscal))
    fermform = np.zeros((gdat.numbener, gdat.numbdqlt, numbpsfpform))
    
    strgpara = ['score', 'gcore', 'stail', 'gtail', 'ntail']
    for m in gdat.indxdqlt:
        if gdat.anlytype.startswith('rec8'):
            irfn = astropy.io.fits.getdata(path, 1 + 3 * gdat.indxdqltincl[m])
            fermscal[m, :] = astropy.io.fits.getdata(path, 2 + 3 * gdat.indxdqltincl[m])['PSFSCALE']
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
            fermform[:, m, k] = sp.interpolate.interp1d(enerirfn, np.mean(irfn[strgpara[k]].squeeze(), axis=0), fill_value='extrapolate')(gdat.bctrpara.ener)
    # convert N_tail to f_core
    for m in gdat.indxdqlt:
        for i in gdat.indxener:
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = np.sqrt((fermscal[None, :, 0] * (10. * gdat.bctrpara.ener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gmod.psfpexpr = np.zeros(gdat.numbener * gdat.numbdqlt * numbpsfpform)
    for m in gdat.indxdqlt:
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
       
    gdat.refr.namepara.elem[0] += ['xpos', 'ypos', 'flux', 'sind', 'otyp', 'lumi']
    gdat.refr.namepara.elem[1] += ['xpos', 'ypos', 'magt', 'reds', 'otyp']


def retr_refrchanfinl(gdat):
    
    booltemp = False
    if gdat.anlytype.startswith('extr'):
        if gdat.numbsidecart == 300:
            gdat.numbpixlxposshft[0] = 1490
            gdat.numbpixlyposshft[0] = 1430
        else:
            booltemp = True
    elif gdat.anlytype.startswith('home'):
        gdat.numbpixlxposshft[0] = 0
        gdat.numbpixlyposshft[0] = 0
    
        if gdat.numbsidecart == 600:
            pass
        elif gdat.numbsidecart == 100:
            indxtile = int(gdat.anlytype[-4:])
            numbsidecntr = int(gdat.anlytype[8:12])
            numbtileside = numbsidecntr / gdat.numbsidecart
            indxtilexaxi = indxtile // numbtileside
            indxtileyaxi = indxtile % numbtileside
            gdat.numbpixlxposshft[0] += indxtilexaxi * gdat.numbsidecart
            gdat.numbpixlyposshft[0] += indxtileyaxi * gdat.numbsidecart
        elif gdat.numbsidecart == 300:
            gdat.numbpixlxposshft[0] += 150
            gdat.numbpixlyposshft[0] += 150
        else:
            booltemp = True
    else:
        booltemp = True

    if booltemp:
        print('')
        print('')
        print('')
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
    xposchan = hdun[1].data['_Glon'] / 180. * pi
    yposchan = hdun[1].data['_Glat'] / 180. * pi
    fluxchansoft = hdun[1].data['SFlux']
    fluxchanhard = hdun[1].data['HFlux']
    objttypechan = hdun[1].data['Otype']
    gdat.refrlumi[0][0] = hdun[1].data['Lx']
    
    # position
    gdat.refr.dictelem[0]['xpos'] = xposchan
    gdat.refr.dictelem[0]['ypos'] = yposchan

    # spectra
    gdat.refrspec = [[np.zeros((3, gdat.numbener, xposchan.size))]]
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
    objttypechantemp = np.zeros(xposchan.size) - 1.
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
    gdat.refrxpos[1] = np.deg2rad(data['_Glon'])
    gdat.refrxpos[1] = ((gdat.refrxpos[1] - pi) % (2. * pi)) - pi
    gdat.refrypos[1] = np.deg2rad(data['_Glat'])
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
    for name in ['xpos', 'ypos', 'sind', 'otyp', 'lumi', 'magt', 'reds']:
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
   
    gdat.refr.namepara.elem[0] += ['xpos', 'ypos', 'flux', 'sind', 'curv', 'expc', 'tvar', 'etag', 'styp', 'sindcolr0001', 'sindcolr0002']
    gdat.refr.namepara.elem[1] += ['xpos', 'ypos', 'flux0400', 'per0', 'per1']


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
    
    gdat.refr.dictelem[0]['xpos'] = np.deg2rad(fgl3['glon'])
    gdat.refr.dictelem[0]['xpos'] = np.pi - ((gdat.refr.dictelem[0]['xpos'] - np.pi) % (2. * np.pi))
    gdat.refr.dictelem[0]['ypos'] = np.deg2rad(fgl3['glat'])
    
    gdat.refr.numbelemfull = gdat.refr.dictelem[0]['xpos'].size

    gdat.refrspec = [np.empty((3, gdat.numbener, gdat.refr.dictelem[0]['xpos'].size))]
    gdat.refrspec[0][0, :, :] = np.stack((fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    
    fgl3specstdvtemp = np.stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.refrspec[0][np.where(np.isfinite(gdat.refrspec[0]) == False)] = 0.
    
    gdat.refrflux[0] = gdat.refrspec[0][:, gdat.indxenerpivt, :]
    gdat.refrsindcolr0001[0] = -np.log(gdat.refrspec[0][:, 1, :] / gdat.refrflux[0]) / np.log(gdat.bctrpara.ener[1] / gdat.enerpivt)
    gdat.refrsindcolr0002[0] = -np.log(gdat.refrspec[0][:, 2, :] / gdat.refrflux[0]) / np.log(gdat.bctrpara.ener[2] / gdat.enerpivt)
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = np.deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3xposstdv = fgl3axisstdv * abs(np.cos(fgl3anglstdv))
    fgl3yposstdv = fgl3axisstdv * abs(np.sin(fgl3anglstdv))

    gdat.refretag[0] = np.zeros(gdat.refr.dictelem[0]['xpos'].size, dtype=object)
    for k in range(gdat.refr.dictelem[0]['xpos'].size):
        gdat.refretag[0][k] = '%s, %s, %s' % (fgl3['Source_Name'][k], fgl3['CLASS1'][k], fgl3['ASSOC1'][k])
    gdat.refrtvar[0] = fgl3['Variability_Index']
    
    gdat.refrstyp[0] = np.zeros_like(gdat.refr.dictelem[0]['xpos']) - 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PowerLaw        ')] = 0
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'LogParabola     ')] = 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLExpCutoff     ')] = 2
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLSuperExpCutoff')] = 3
    indx = np.where(gdat.refrstyp[0] == -1)[0]
    if indx.size > 0:
        print('')
        print('')
        print('')
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
   
    gdat.refrxpos[1] = np.deg2rad(data['glon'])
    gdat.refrxpos[1] = ((gdat.refrxpos[1] - np.pi) % (2. * np.pi)) - np.pi
    gdat.refrypos[1] = np.deg2rad(data['glat'])
    
    gdat.refrper0[1] = data['P0']
    gdat.refrper1[1] = data['P1']
    gdat.refrflux0400[1] = data['S400']
    #gdat.refrdism[1] = data['DM']
    #gdat.refrdlos[1] = data['Dist']

    # error budget
    for name in ['xpos', 'ypos', 'per0', 'per1', 'flux0400', 'tvar', 'styp']:
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


def retr_xposypos(gang, aang):
    
    xpos = gang * np.cos(aang)
    ypos = gang * np.sin(aang)

    return xpos, ypos


def retr_gang(xpos, ypos):
    
    gang = np.arccos(np.cos(xpos) * np.cos(ypos))

    return gang


def retr_aang(xpos, ypos):

    aang = np.arctan2(ypos, xpos)

    return aang


def show_paragenrscalfull(gdat, gdatmodi, strgstat='this', strgmodl='fitt', indxsampshow=None):
    
    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    print('strgmodl: ' + strgmodl)
    print('strgstat: ' + strgstat)
    print('%5s %20s %30s %30s %15s' % ('index', 'namepara', 'paragenrunitfull', 'paragenrscalfull', 'scalpara'))
    for k in gmod.indxpara.genrfull:
        
        if indxsampshow is not None and not k in indxsampshow:
            continue
        
        if gdat.booldiag:
            if gmod.namepara.genr[k] == '0.0':
                print('')
                print('')
                print('')
                raise Exception('')
            if gmod.scalpara.genr[k] == '0.0':
                print('')
                print('')
                print('')
                raise Exception('')
        
        if gmod.numbpopl > 0:
            
            booltemp = False
            for l in gmod.indxpopl:
                if k == gmod.indxparagenrelemsing[l][0]:
                    booltemp = True
            if booltemp:
                print('')
        print('%5d %20s %30g %30g %15s' % (k, gmod.namepara.genr[k], gmodstat.paragenrunitfull[k], gmodstat.paragenrscalfull[k], gmod.scalpara.genr[k]))
    

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
    
    if gmod.numbpopl > 0:
        if gdat.booldiag:
            for l in gmod.indxpopl:
                if len(gmodthis.indxelemfull[l]) > len(set(gmodthis.indxelemfull[l])):
                    print('')
                    print('')
                    print('')
                    raise Exception('Repeating entry in the element index list!')

        thisindxparagenrfullelem = retr_indxparagenrelemfull(gdat, gmodthis.indxelemfull, strgmodl)
        setattr(gmodthis, 'indxparagenrfullelem', thisindxparagenrfullelem)
    else:
        thisindxparagenrfullelem = None
    
    gdatmodi.this.boolpropfilt = True 

    # index of the population in which a transdimensional proposal will be attempted
    if gmod.numbpopl > 0:
        if thisindxpopl is None:
            gdatmodi.indxpopltran = np.random.choice(gmod.indxpopl)
        else:
            gdatmodi.indxpopltran = thisindxpopl
        numbelemtemp = gmodthis.paragenrscalfull[gmod.indxpara.numbelem[gdatmodi.indxpopltran]]
    
    # forced death or birth does not check for the prior on the dimensionality on purpose!
    if gmod.numbpopl > 0 and (deth or brth or np.random.rand() < gdat.probtran) and \
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
        
        if gdat.booldiag and (gdatmodi.stdp > 1e2).any():
            print('')
            print('')
            print('')
            raise Exception('')

        thisindxparagenrfullelemconc = []
        for l in gmod.indxpopl:
            thisindxparagenrfullelemconc.append(thisindxparagenrfullelem[l]['full'])

        # get the indices of the current parameter vector
        if gmod.numbpopl > 0:
            thisindxsampfull = np.concatenate([gmod.indxpara.genrbasepert] + thisindxparagenrfullelemconc)
        else:
            thisindxsampfull = gmod.indxpara.genrbasepert
        
        thisstdp = gdatmodi.stdp[gdat.indxstdppara[thisindxsampfull]]
        if not np.isfinite(thisstdp).all():
            print('')
            print('')
            print('')
            raise Exception('')
        gdatmodi.this.indxproptype = 0
    
    if gdat.booldiag:
        if gdat.probspmr == 0 and gdatmodi.this.indxproptype > 2:
            print('')
            print('')
            print('')
            raise Exception('')

    if gdat.typeverb > 1:
        print('gdatmodi.this.indxproptype')
        print(gdatmodi.this.indxproptype)

    if gdatmodi.this.indxproptype == 0:
        gmodnext.paragenrunitfull = np.copy(gmodthis.paragenrunitfull)
        if gmod.numbpopl > 0:
            gmodnext.indxelemfull = gmodthis.indxelemfull
    if gdatmodi.this.indxproptype > 0:
        gmodnext.paragenrunitfull = np.copy(gmodthis.paragenrunitfull)
        gmodnext.paragenrscalfull = np.copy(gmodthis.paragenrscalfull)
        if gmod.numbpopl > 0:
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
        
        if gdat.booldiag:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(gmodnext.paragenrunitfull).all():
                raise Exception('')

        indxsamplowr = np.where(gmodnext.paragenrunitfull[gmod.numbpopl:] < 0.)[0]
        if indxsamplowr.size > 0:
            gmodnext.paragenrunitfull[gmod.numbpopl+indxsamplowr] = abs(gmodnext.paragenrunitfull[gmod.numbpopl+indxsamplowr]) % 1.
        
        if gdat.booldiag:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

        indxsampuppr = np.where(gmodnext.paragenrunitfull[gmod.numbpopl:] > 1.)[0]
        if indxsampuppr.size > 0:
            gmodnext.paragenrunitfull[gmod.numbpopl+indxsampuppr] = (gmodnext.paragenrunitfull[gmod.numbpopl+indxsampuppr] - 1.) % 1.
        
        if gdat.booldiag:
            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 1).any():
                raise Exception('')

            if (gmodnext.paragenrunitfull[gmod.numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(gmodnext.paragenrunitfull).all():
                raise Exception('')

        gmodnext.paragenrscalfull = icdf_paragenrscalfull(gdat, strgmodl, gmodnext.paragenrunitfull, thisindxparagenrfullelem)

        if gdat.booldiag:
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
        gdatmodi.indxpara.genrfullelemaddd = retr_indxparagenrelem(gmod, gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.append(gdatmodi.indxpara.genrfullelemaddd)
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
        indxparagenrfullelemdeth = retr_indxparagenrelem(gmod, gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
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
        gdatmodi.indxpara.genrfullelemfrst = retr_indxparagenrelem(gmod, l, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.insert(0, gdatmodi.indxpara.genrfullelemfrst)
        
        # sample indices for the second element
        gdatmodi.indxsampseco = gdatmodi.indxpara.genrfullelemaddd
        
        # take the parent element parameters
        for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            gdatmodi.comppare[k] = np.copy(gmodthis.paragenrscalfull[thisindxparagenrfullelem[gdatmodi.indxpopltran][nameparagenrelem][gdatmodi.indxelemfullmodi[0]]])
        
        # draw the auxiliary parameters
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if gmod.boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.this.auxipara[g] = np.random.randn() * gdat.radispmr
            elif g == gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.this.auxipara[g] = np.random.rand()
            else:
                gdatmodi.this.auxipara[g] = icdf_trap(gdat, strgmodl, np.random.rand(), gmodthis.paragenrscalfull, gmod.scalpara.genrelem[gdatmodi.indxpopltran][g], \
                                                                                                     gmod.namepara.genrelem[gdatmodi.indxpopltran][g], l)

        # determine the new parameters
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.this.auxipara[1]) * gdatmodi.this.auxipara[0]
        else:
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.this.auxipara[2]) * gdatmodi.this.auxipara[0]
            gdatmodi.compfrst[1] = gdatmodi.comppare[1] + (1. - gdatmodi.this.auxipara[2]) * gdatmodi.this.auxipara[1]
        gdatmodi.compfrst[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] = gdatmodi.this.auxipara[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] * \
                                                                                                        gdatmodi.comppare[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]
        if gmod.typeelem[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.this.auxipara[1] * gdatmodi.this.auxipara[0]
        else:
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.this.auxipara[2] * gdatmodi.this.auxipara[0]
            gdatmodi.compseco[1] = gdatmodi.comppare[1] - gdatmodi.this.auxipara[2] * gdatmodi.this.auxipara[1]
        gdatmodi.compseco[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] = (1. - gdatmodi.this.auxipara[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]) * \
                                                                                                        gdatmodi.comppare[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]
        for g in range(gmod.numbparagenrelemsing[gdatmodi.indxpopltran]):
            if not gmod.boolcompposi[gdatmodi.indxpopltran][g] and g != gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]:
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
            if np.fabs(gdatmodi.compfrst[0]) > maxmxpos or np.fabs(gdatmodi.compseco[0]) > maxmxpos or \
                                                                    np.fabs(gdatmodi.compfrst[1]) > maxmypos or np.fabs(gdatmodi.compseco[1]) > maxmypos:
                gdatmodi.this.boolpropfilt = False
        if gdatmodi.compfrst[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] < getattr(gmod.minmpara, gmod.nameparagenrelemampl[gdatmodi.indxpopltran]) or \
           gdatmodi.compseco[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] < getattr(gmod.minmpara, gmod.nameparagenrelemampl[gdatmodi.indxpopltran]):
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
        if gdat.booldiag:
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
        gdatmodi.indxpara.genrfullelemfrst = retr_indxparagenrelem(gmod, l, gdatmodi.mergindxelemfrst)
        gdatmodi.indxsamptran.append(gdatmodi.indxpara.genrfullelemfrst)

        ## second element index to be merged
        gdatmodi.mergindxelemseco = gmodthis.indxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullmergseco]
       
        ## second element
        gdatmodi.indxpara.genrfullelemseco = retr_indxparagenrelem(gmod, l, gdatmodi.mergindxelemseco)
        gdatmodi.indxsamptran.append(gdatmodi.indxpara.genrfullelemseco)
        
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
        gdatmodi.this.auxipara[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] / \
                                        (gdatmodi.compfrst[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] + gdatmodi.compseco[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]) 
        for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[gdatmodi.indxpopltran]):
            if not gmod.boolcompposi[gdatmodi.indxpopltran][g] and g != gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]:
                gdatmodi.this.auxipara[g] = gdatmodi.compseco[g]

        # merged element
        gdatmodi.comppare[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] + \
                                                                                                gdatmodi.compseco[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]
        if gdatmodi.comppare[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]] > getattr(gdat, 'maxm' + gmod.nameparagenrelemampl[gdatmodi.indxpopltran]):
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
                gdatmodi.comppare[g] = gdatmodi.compfrst[g] + (1. - gdatmodi.this.auxipara[gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]]) * \
                                                                                            (gdatmodi.compseco[g] - gdatmodi.compfrst[g])
            elif g == gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]:
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
            print('indxparagenrfullelemfrst: ', gdatmodi.indxpara.genrfullelemfrst)
            print('indxparagenrfullelemseco: ', gdatmodi.indxpara.genrfullelemseco)

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
            print('xposfrst: ', gdat.anglfact * gdatmodi.compfrst[0])
            print('yposfrst: ', gdat.anglfact * gdatmodi.compfrst[1])
            print('amplfrst: ', gdatmodi.compfrst[2])
            print('xposseco: ', gdat.anglfact * gdatmodi.compseco[0])
            print('yposseco: ', gdat.anglfact * gdatmodi.compseco[1])
            print('amplseco: ', gdatmodi.compseco[2])
            print('xpospare: ', gdat.anglfact * gdatmodi.comppare[0])
            print('ypospare: ', gdat.anglfact * gdatmodi.comppare[1])
            print('fluxpare: ', gdatmodi.comppare[2])
            print('auxipara[0][0]: ', gdat.anglfact * gdatmodi.this.auxipara[0])
            print('auxipara[0][1]: ', gdat.anglfact * gdatmodi.this.auxipara[1])
            print('auxipara[0][2]: ', gdatmodi.this.auxipara[2])
                
    if gmod.numbpopl > 0 and gdatmodi.this.indxproptype > 0 and gdatmodi.this.boolpropfilt:
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
    
    if gmod.numbpopl > 0:
        if gdatmodi.this.indxproptype == 0:
            indxparagenrfullelem = thisindxparagenrfullelem
        else:
            indxparagenrfullelem = retr_indxparagenrelemfull(gdat, gmodnext.indxelemfull, strgmodl)
    if gdat.typeverb > 1:
        print('gdatmodi.indxsampmodi')
        print(gdatmodi.indxsampmodi)
        if gmod.numbpopl > 0:
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
                                                                            gmod.scalpara.genrelem[gdatmodi.indxpopltran][g], \
                                                                            gmod.namepara.genrelem[gdatmodi.indxpopltran][g], gdatmodi.indxpopltran)

    if gdat.booldiag:
        if gmod.numbpopl > 0:
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
        lpautemp = 0.5 * gdat.factpriodoff * gmod.numbparagenrelemsing[gdatmodi.indxpopltran]
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
            elif g != gmod.indxpara.genrelemampl[gdatmodi.indxpopltran]:
                dictelemtemp[0][nameparagenrelem] = gdatmodi.this.auxipara[g]
                gdatmodi.this.lpau += retr_lprielem(gdat, 'fitt', gdatmodi.indxpopltran, g, \
                                            gmod.namepara.genrelem[gdatmodi.indxpopltran][g], gmod.scalpara.genrelem[gdatmodi.indxpopltran][g], \
                                            gdatmodi.this.paragenrscalfull, dictelemtemp, [1])
        if gdatmodi.this.indxproptype == 4:
            gdatmodi.this.lpau *= -1.

    if gdatmodi.this.indxproptype > 2 and gdatmodi.this.boolpropfilt:
        ## the ratio of the probability of the reverse and forward proposals, and
        if gdatmodi.this.indxproptype == 3:
            gdatmodi.this.probmergtotl = retr_probmerg(gdat, gdatmodi, gdatmodi.next.paragenrscalfull, gdatmodi.next.indxparagenrelemfull, gdatmodi.indxpopltran, 'pair', \
                                                                               typeelem=gmod.typeelem)
            gdatmodi.this.ltrp = np.log(gdatmodi.this.numbelem[gdatmodi.indxpopltran] + 1) + np.log(gdatmodi.this.probmergtotl)

        else:
            gdatmodi.this.probmergtotl = retr_probmerg(gdat, gdatmodi, gdatmodi.this.paragenrscalfull, gdatmodi.this.indxparagenrelemfull, gdatmodi.indxpopltran, 'pair', \
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


def retr_indxparagenrelemfull(gdat, indxelemfull, strgmodl):
    '''
    Compute the generative parameter indices of the occupied (full) elements
    '''
    
    gmod = getattr(gdat, strgmodl)
    
    ## element parameters
    if gmod.numbpopl > 0:
        indxparagenrfullelem = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            indxparagenrfulltemp = gmod.numbparagenrbase + gmod.numbparagenrelemcuml[l] + np.array(indxelemfull[l], dtype=int) * gmod.numbparagenrelemsing[l]
            
            cntr = tdpy.cntr()
            
            indxparagenrfullelem[l] = dict()
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                indxparagenrfullelem[l][nameparagenrelem] = indxparagenrfulltemp + cntr.incr()
            indxparagenrfullelem[l]['full'] = np.repeat(indxparagenrfulltemp, gmod.numbparagenrelemsing[l]) + np.tile(gmod.indxparagenrelemsing[l], len(indxelemfull[l]))
        
        if gdat.booldiag:
            
            for l in gmod.indxpopl:
                if len(indxparagenrfullelem[l]['full']) > 0:
                    
                    if np.amax(indxparagenrfullelem[l]['full']) > gmod.numbparagenrelempopl[l] + gmod.numbparagenrbase or \
                       gmod.this.numbelempop0 < gmod.minmpara.numbelempop0 or \
                       gmod.this.numbelempop0 > gmod.maxmpara.numbelempop0:
                        print('')
                        print('')
                        print('')
                        print('l')
                        print(l)
                        print('strgmodl')
                        print(strgmodl)
                        print('indxparagenrfullelem[l][full]')
                        summgene(indxparagenrfullelem[l]['full'])
                        print('gmod.numbparagenrelempopl[l]')
                        print(gmod.numbparagenrelempopl[l])
                        print('gmod.numbparagenrbase')
                        print(gmod.numbparagenrbase)
                        print('gmod.numbparagenrelemcuml[l]')
                        print(gmod.numbparagenrelemcuml[l])
                        print('np.array(indxelemfull[l], dtype=int)')
                        print(np.array(indxelemfull[l], dtype=int))
                        print('gmod.numbparagenrelemsing[l]')
                        print(gmod.numbparagenrelemsing[l])
                        print('indxelemfull')
                        print(indxelemfull)
                        print('gmod.this.numbelempop0')
                        print(gmod.this.numbelempop0)
                        print('gmod.minmpara.numbelempop0')
                        print(gmod.minmpara.numbelempop0)
                        print('gmod.maxmpara.numbelempop0')
                        print(gmod.maxmpara.numbelempop0)
                        if np.amax(indxparagenrfullelem[l]['full']) > gmod.numbparagenrelempopl[l] + gmod.numbparagenrbase:
                            raise Exception('np.amax(indxparagenrfullelem[l][full]) > gmod.numbparagenrelempopl[l] + gmod.numbparagenrbase')
                        elif gmod.this.numbelempop0 < gmod.minmpara.numbelempop0:
                            raise Exception('gmod.this.numbelempop0 < gmod.minmpara.numbelempop0')
                        else:
                            raise Exception('gmod.this.numbelempop0 > gmod.maxmpara.numbelempop0')
        
    else:
        indxparagenrfullelem = None
    
    return indxparagenrfullelem
    

def retr_weigmergodim(gdat, elin, elinothr):
    
    weigmerg = np.exp(-0.5 * ((elin - elinothr) / gdat.radispmr)**2)
    
    return weigmerg


def retr_weigmergtdim(gdat, xpos, xposothr, ypos, yposothr):
    
    weigmerg = np.exp(-0.5 * (((xpos - xposothr) / gdat.radispmr)**2 + ((ypos - yposothr) / gdat.radispmr)**2))
    
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
            xpostotl = paragenrscalfull[indxparagenrfullelem['xpos'][indxpopltran]]
            ypostotl = paragenrscalfull[indxparagenrfullelem['ypos'][indxpopltran]]
            xpos = xpostotl[gdatmodi.indxelemfullmodi[0]]
            ypos = ypostotl[gdatmodi.indxelemfullmodi[0]]
            xposothr = np.concatenate((xpostotl[:gdatmodi.indxelemfullmodi[0]], xpostotl[gdatmodi.indxelemfullmodi[0]+1:]))
            yposothr = np.concatenate((ypostotl[:gdatmodi.indxelemfullmodi[0]], ypostotl[gdatmodi.indxelemfullmodi[0]+1:]))
            weigmerg = retr_weigmergtdim(gdat, xpos, xposothr, ypos, yposothr)
        listweigmerg.append(weigmerg) 

    # determine the probability of merging the second element given the first element
    if strgtype == 'seco':
        probmerg = listweigmerg[0] / np.sum(listweigmerg[0])
    
    # determine the probability of merging the pair
    if strgtype == 'pair':
        if gmod.typeelem[indxpopltran].startswith('lghtline'):
            weigpair = retr_weigmergtdim(gdat, elin, elintotl[gdatmodi.indxelemfullmodi[1]])
        else:
            weigpair = retr_weigmergtdim(gdat, xpos, xpostotl[gdatmodi.indxelemfullmodi[1]], ypos, ypostotl[gdatmodi.indxelemfullmodi[1]])
        probmerg = weigpair / np.sum(listweigmerg[0]) + weigpair / np.sum(listweigmerg[1])
        
    if gdat.booldiag:
        if not np.isfinite(probmerg).all():
            raise Exception('Merge probability is infinite.')

    return probmerg

    
def retr_indxparagenrelem(gmod, l, u):

    indxparagenrelem = gmod.numbparagenrbase + gmod.numbparagenrelemcuml[l] + u * gmod.numbparagenrelemsing[l] + gmod.indxparagenrelemsing[l]

    return indxparagenrelem


def gang_detr():

    gang, aang, xpos, ypos = sympy.symbols('gang aang xpos ypos')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])


def retr_psfn(gdat, psfp, indxenertemp, thisangl, typemodlpsfn, strgmodl):

    gmod = getattr(gdat, strgmodl)
    
    indxpsfpinit = gmod.numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxdqlt[None, :])
    
    if gdat.typeexpr == 'ferm':
        scalangl = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(gdat.blimpara.angl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
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
        fact = 2. * np.pi * np.trapz(psfnnorm * np.sin(gdat.blimpara.angl[None, :, None]), gdat.blimpara.angl, axis=1)[:, None, :]
        psfn /= fact

    return psfn


def retr_unit(xpos, ypos):

    xdat = np.cos(ypos) * np.cos(xpos)
    ydat = -np.cos(ypos) * np.sin(xpos)
    zaxi = np.sin(ypos)

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
        indxmpol = np.where((gdat.bctrpara.mpol > gdat.blimpara.mpolodim[k]) & (gdat.bctrpara.mpol < gdat.blimpara.mpolodim[k+1]))
        psecodim[k] = np.mean(psec[indxmpol])
    # temp
    #psecodim *= gdat.bctrpara.mpolodim**2
    
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
    arrystks = np.zeros((numbstks, gmod.numbparagenrelem))
    for n in gdat.indxsamptotl:
        indxparagenrfullelem = retr_indxparagenrelemfull(gdat, gdat.listpostindxelemfull[n], 'fitt') 
        for l in gmod.indxpopl:
            for k in np.arange(len(gdat.listpostindxelemfull[n][l])):
                for m, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    arrystks[indxstks[n][l][k], m] = gdat.listpostparagenrscalfull[n, gmodstat.indxparagenrelemfull[l][nameparagenrelem][k]]

    if gdat.typeverb > 0:
        print('Constructing the distance matrix for %d stacked samples...' % arrystks.shape[0])
        timeinit = gdat.functime()
    
    gdat.distthrs = np.empty(gmod.numbparagenrelem)
    for k, nameparagenrelem in enumerate(gmod.namepara.elem):
       # temp
       l = 0
       gdat.distthrs[k] = gdat.stdp[getattr(gdat, 'indxstdppop%d' % l + nameparagenrelem)]
    
    # construct lists of samples for each proposal type
    listdisttemp = [[] for k in range(gmod.numbparagenrelem)]
    indxstksrows = [[] for k in range(gmod.numbparagenrelem)]
    indxstkscols = [[] for k in range(gmod.numbparagenrelem)]
    thisperc = 0
    cntr = 0
    for k in gmod.indxpara.genrelemtotl:
        for n in range(numbstks):
            dist = np.fabs(arrystks[n, k] - arrystks[:, k])
            indxstks = np.where(dist < gdat.distthrs[k])[0]
            if indxstks.size > 0:
                for j in indxstks:
                    cntr += 1
                    listdisttemp[k].append(dist[j])
                    indxstksrows[k].append(n)
                    indxstkscols[k].append(j)
            
            nextperc = np.floor(100. * float(k * numbstks + n) / numbstks / gmod.numbparagenrelem)
            if nextperc > thisperc:
                thisperc = nextperc
            if cntr > 1e6:
                break
        
        listdisttemp[k] = np.array(listdisttemp[k])
        indxstksrows[k] = np.array(indxstksrows[k])
        indxstkscols[k] = np.array(indxstkscols[k])

        if cntr > 1e6:
            break
    
    listdist = [[] for k in range(gmod.numbparagenrelem)]
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
            indxindx = np.where((listdist[0][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmxpos < gdat.anglassc) & \
                             (listdist[1][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmypos < gdat.anglassc))[0]
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
                    for k in gmod.indxpara.genrelemtotl:
                        temp = listdist[k][indxstkscntr, indxstkstemp].tonp.array()[0]
                        totl = totl + temp**2

                    indxleft = np.argsort(totl)[0]
                    
                    indxstksthis = indxstkstemp[indxleft]
                
                    thisbool = True
                    for k in gmod.indxpara.genrelemtotl:
                        if listdist[k][indxstkscntr, indxstksthis] > gdat.distthrs[k]:
                            thisbool = False

                    if thisbool:
                        indxstksassc[cntr].append(indxstksthis)
                        indxstksleft.remove(indxstksthis)
            
                # temp
                #if gdat.boolmakeplot:
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

                #    plot_genemaps(gdat, gdatmodi, 'this', 'cntpdata', strgpdfn, indxenerplot=0, indxdqltplot=0, cond=True)
                
            cntr += 1
        
    gdat.dictglob['poststkscond'] = []
    gdat.dictglob['liststkscond'] = []
    # for each condensed element
    for r in range(len(indxstksassc)): 
        gdat.dictglob['liststkscond'].append([])
        gdat.dictglob['liststkscond'][r] = {}
        gdat.dictglob['poststkscond'].append([])
        gdat.dictglob['poststkscond'][r] = {}
        for strgfeat in gmod.namepara.genr.elem:
            gdat.dictglob['liststkscond'][r][strgfeat] = []

        # for each associated sample associated with the central stacked sample 
        for k in range(len(indxstksassc[r])):
            indxsamptotlcntr = indxtupl[indxstksassc[r][k]][0]
            indxpoplcntr = indxtupl[indxstksassc[r][k]][1]
            indxelemcntr = indxtupl[indxstksassc[r][k]][2]
            
            for strgfeat in gmod.namepara.genr.elem:
                temp = getattr(gdat, 'list' + strgfeat)
                if temp[indxsamptotlcntr][indxpoplcntr].size > 0:
                    temp = temp[indxsamptotlcntr][indxpoplcntr][..., indxelemcntr]
                    gdat.dictglob['liststkscond'][r][strgfeat].append(temp)

    for r in range(len(gdat.dictglob['liststkscond'])):
        for strgfeat in gmod.namepara.genr.elem:
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
    setp_varb(gdat, 'prvl')
    gdat.histprvl = np.histogram(gdat.prvl, bins=gdat.blimpara.prvl)[0]
    if gdat.boolmakeplot:
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
    
    if gdat.boolbindspat:
        xpos = listposi[0]
        ypos = listposi[1]
        indxpixlpnts = retr_indxpixl(gdat, ypos, xpos)
    else:
        elin = listposi[0]
        indxpixlpnts = np.zeros_like(elin, dtype=int)
    for k in range(spec.shape[1]):
        cnts[:, k] += spec[:, k] * gdat.expototl[:, indxpixlpnts[k]]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None]
    cnts = np.sum(cnts, axis=0)

    return cnts


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


def retr_rele(gdat, maps, xpos, ypos, defs, asca, acut, indxpixlelem, absv=True, cntpmodl=None):
    
    grad = retr_gradmaps(gdat, maps)
        
    defl = chalcedon.retr_defl(gdat.xposgrid, gdat.yposgrid, indxpixlelem, xpos, ypos, defs, asca=asca, acut=acut)

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


def setp_indxpara_init(gdat, strgmodl):
    
    print('setp_indxpara_init(): Initial setup for model %s...' % (strgmodl))

    gmod = getattr(gdat, strgmodl)
    
    if strgmodl == 'fitt':
        gmod.lablmodl = 'Model'
    if strgmodl == 'true':
        gmod.lablmodl = 'True'
    
    # transdimensional element populations
    gmod.numbpopl = len(gmod.typeelem)
    gmod.indxpopl = np.arange(gmod.numbpopl)
    
    gmod.namepara.genrelem = [[] for l in gmod.indxpopl]
    gmod.scalpara.genrelem = [[] for l in gmod.indxpopl]
    gmod.namepara.derielemodim = [[] for l in gmod.indxpopl]
    
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
    
    # name of the generative element parameter used for the amplitude
    gmod.nameparagenrelemampl = [[] for l in gmod.indxpopl]
    gmod.indxpara.genrelemampl = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l] == 'lghtpntspuls':
            gmod.nameparagenrelemampl[l] = 'per0'
            gmod.indxpara.genrelemampl[l] = 2
        elif gmod.typeelem[l] == 'lghtpntsagnntrue':
            gmod.nameparagenrelemampl[l] = 'lum0'
            gmod.indxpara.genrelemampl[l] = 2
        elif gmod.typeelem[l].startswith('lghtline'):
            gmod.nameparagenrelemampl[l] = 'flux'
            gmod.indxpara.genrelemampl[l] = 1
        elif gmod.typeelem[l].startswith('lghtpnts'):
            gmod.nameparagenrelemampl[l] = 'flux'
            gmod.indxpara.genrelemampl[l] = 2
        elif gmod.typeelem[l].startswith('lghtgausbgrd'):
            gmod.nameparagenrelemampl[l] = 'flux'
            gmod.indxpara.genrelemampl[l] = 2
        if gmod.typeelem[l] == 'lens':
            gmod.nameparagenrelemampl[l] = 'defs'
            gmod.indxpara.genrelemampl[l] = 2
        if gmod.typeelem[l].startswith('clus'):
            gmod.nameparagenrelemampl[l] = 'nobj'
            gmod.indxpara.genrelemampl[l] = 2
        if gmod.typeelem[l] == 'lens':
            gmod.nameparagenrelemampl[l] = 'defs'
        if gmod.typeelem[l] == 'clus':
            gmod.nameparagenrelemampl[l] = 'nobj'
        if len(gmod.nameparagenrelemampl[l]) == 0:
            raise Exception('Amplitude feature undefined.')
    
    # galaxy components
    gmod.indxsersfgrd = np.arange(gmod.numbsersfgrd)


def setp_indxpara_finl(gdat, strgmodl):
    
    print('setp_indxpara_finl(): Building parameter indices for model %s...' % (strgmodl))
        
    gmod = getattr(gdat, strgmodl)
    
    # element setup
    ## Boolean flag to calculate the approximation error due to the finite PSF kernel
    boolcalcerrr = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        # if the PSF evaluation is local and the number of pixels is small, do it
        if gmod.typeelemspateval[l] == 'locl' and gdat.numbpixlfull < 1e5:
            boolcalcerrr[l] = True
        else:
            boolcalcerrr[l] = False
    setp_varb(gdat, 'boolcalcerrr', valu=boolcalcerrr, strgmodl=strgmodl)
    
    ## minimum number of elements summed over all populations
    gmod.minmpara.numbelemtotl = np.amin(gmod.minmpara.numbelem)
            
    ## maximum number of elements summed over all populations
    gmod.maxmpara.numbelemtotl = np.sum(gmod.maxmpara.numbelem) 

    ## name of the element feature used to sort the elements, for each population
    gmod.nameparaelemsort = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lght'):
            gmod.nameparaelemsort[l] = 'flux'
        if gmod.typeelem[l] == 'lens':
            gmod.nameparaelemsort[l] = 'defs'
        if gmod.typeelem[l].startswith('clus'):
            gmod.nameparaelemsort[l] = 'nobj'
    
    ## subscript used to denote the elements of each population
    gmod.lablelemsubs = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gdat.numbgrid > 1:
            if gmod.typeelem[l] == 'lghtpnts':
                gmod.lablelemsubs[l] = r'\rm{fps}'
            if gmod.typeelem[l] == 'lghtgausbgrd':
                gmod.lablelemsubs[l] = r'\rm{bgs}'
        else:
            if gmod.typeelem[l].startswith('lghtpntspuls'):
                gmod.lablelemsubs[l] = r'\rm{pul}'
            if gmod.typeelem[l].startswith('lghtpntsagnn'):
                gmod.lablelemsubs[l] = r'\rm{agn}'
            elif gmod.typeelem[l] == 'lghtpnts':
                gmod.lablelemsubs[l] = r'\rm{pts}'
        if gmod.typeelem[l] == 'lens':
            gmod.lablelemsubs[l] = r'\rm{sub}'
        if gmod.typeelem[l].startswith('clus'):
            gmod.lablelemsubs[l] = r'\rm{cls}'
        if gmod.typeelem[l].startswith('lghtline'):
            gmod.lablelemsubs[l] = r'\rm{lin}'
    
    ## indices of element populations for each coordinate grid
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
    
    ## indices of the coordinate grids for each population
    indxgridpopl = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        for y in gdat.indxgrid:
            if l in gmod.indxpoplgrid[y]:
                indxgridpopl[l] = y
    
    ## Boolean flag to calculate the surface brightness due to elements
    gmod.boolcalcelemsbrt = False
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lghtpnts'):
            gmod.boolcalcelemsbrt = True
    
    ## Boolean flag to calculate the surface brightness due to elements in the background coordinate grid
    if 'lghtgausbgrd' in gmod.typeelem:
        gmod.boolcalcelemsbrtbgrd = True
    else:
        gmod.boolcalcelemsbrtbgrd = False
    
    ## Boolean flag indicating, for each element population, whether elements are light sources
    gmod.boolelemlght = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lght'):
            gmod.boolelemlght[l] = True
        else:
            gmod.boolelemlght[l] = False
    
    ## Boolean flag indicating whether any population of elements is a light source
    gmod.boolelemlghtanyy = True in gmod.boolelemlght
    
    ## Boolean flag indicating whether the surface brightness on the background grid is lensed to due elements
    gmod.boolelemlens = False
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lens'):
            gmod.boolelemlens = True
    
    ## Boolean flag indicating, for each population of elements, whether the elements are delta-function light sources
    gmod.boolelemsbrtdfnc = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.maxmpara.numbelem[l] > 0 and (gmod.typeelem[l].startswith('lght') and not gmod.typeelem[l].endswith('bgrd') or gmod.typeelem[l].startswith('clus')):
            gmod.boolelemsbrtdfnc[l] = True
        else:
            gmod.boolelemsbrtdfnc[l] = False
    ## Boolean flag indicating whether any population of elements has delta-function light sources
    gmod.boolelemsbrtdfncanyy = True in gmod.boolelemsbrtdfnc

    ## Boolean flag indicating, for each population of elements, whether the elements are subhalos deflecting the light from the background grid
    gmod.boolelemdeflsubh = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l] == 'lens':
            gmod.boolelemdeflsubh[l] = True
        else:
            gmod.boolelemdeflsubh[l] = False
    ## Boolean flag indicating whether any population of elements has elements that are subhalos deflecting the light from the background grid
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
    
    if gdat.boolbindspat:
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
        if gdat.typeexpr == 'gmix':
            minmflux = 0.1
            maxmflux = 100.
        if gdat.typeexpr.startswith('HST_WFC3'):
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
            if gdat.typeexpr.startswith('HST_WFC3'):
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
            setp_varb(gdat, 'yposdistscal', limt=[0.5 / gdat.anglfact, 5. / gdat.anglfact], popl='full', strgmodl=strgmodl)
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
        gmod.sbrtbacknorm[c] = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
        if gmod.typeback[c] == 'data':
            gmod.sbrtbacknorm[c] = np.copy(gdat.sbrtdata)
            gmod.sbrtbacknorm[c][np.where(gmod.sbrtbacknorm[c] == 0.)] = 1e-100
        elif isinstance(gmod.typeback[c], float):
            gmod.sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull)) + gmod.typeback[c]
        elif isinstance(gmod.typeback[c], list) and isinstance(gmod.typeback[c], float):
            gmod.sbrtbacknorm[c] = retr_spec(gdat, np.array([gmod.typeback[c]]), sind=np.array([gmod.typeback[c]]))[:, 0, None, None]
        elif isinstance(gmod.typeback[c], np.ndarray) and gmod.typeback[c].ndim == 1:
            gmod.sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull)) + gmod.typeback[c][:, None, None]
        elif gmod.typeback[c].startswith('bfunfour') or gmod.typeback[c].startswith('bfunwfou'):
            indxexpatemp = int(gmod.typeback[c][8:]) 
            indxterm = indxexpatemp // ordrexpa**2
            indxexpaxdat = (indxexpatemp % ordrexpa**2) // ordrexpa + 1
            indxexpaydat = (indxexpatemp % ordrexpa**2) % ordrexpa + 1
            if namebfun == 'bfunfour':
                ampl = 1.
                func = gdat.bctrpara.yposcart 
            if namebfun == 'bfunwfou':
                functemp = np.exp(-0.5 * (gdat.bctrpara.yposcart / (1. / gdat.anglfact))**2)
                ampl = np.sqrt(functemp)
                func = functemp
            argsxpos = 2. * np.pi * indxexpaxdat * gdat.bctrpara.xposcart / gdat.maxmgangdata
            argsypos = 2. * np.pi * indxexpaydat * func / gdat.maxmgangdata
            if indxterm == 0:
                termfrst = np.sin(argsxpos)
                termseco = ampl * np.sin(argsypos)
            if indxterm == 1:
                termfrst = np.sin(argsxpos)
                termseco = ampl * np.cos(argsypos)
            if indxterm == 2:
                termfrst = np.cos(argsxpos)
                termseco = ampl * np.sin(argsypos)
            if indxterm == 3:
                termfrst = np.cos(argsxpos)
                termseco = ampl * np.cos(argsypos)
            gmod.sbrtbacknorm[c] = (termfrst[None, :] * termseco[:, None]).flatten()[None, :, None] * \
                                                        np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
            
        else:
            path = gdat.pathinpt + gmod.typeback[c]
            gmod.sbrtbacknorm[c] = astropy.io.fits.getdata(path)
            
            if gdat.typepixl == 'cart':
                if not gdat.boolforccart:
                    if gmod.sbrtbacknorm[c].shape[2] != gdat.numbsidecart:
                        raise Exception('Provided background template must have the chosen image dimensions.')
                
                gmod.sbrtbacknorm[c] = gmod.sbrtbacknorm[c].reshape((gmod.sbrtbacknorm[c].shape[0], -1, gmod.sbrtbacknorm[c].shape[-1]))
    
            if gdat.typepixl == 'cart' and gdat.boolforccart:
                sbrtbacknormtemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
                for i in gdat.indxenerfull:
                    for m in gdat.indxdqltfull:
                        sbrtbacknormtemp[i, :, m] = tdpy.retr_cart(gmod.sbrtbacknorm[c][i, :, m], \
                                                numbsidexpos=gdat.numbsidecart, numbsideypos=gdat.numbsidecart, \
                                                minmxpos=gdat.anglfact*gdat.minmxposdata, maxmxpos=gdat.anglfact*gdat.maxmxposdata, \
                                                minmypos=gdat.anglfact*gdat.minmyposdata, maxmypos=gdat.anglfact*gdat.maxmyposdata).flatten()
                gmod.sbrtbacknorm[c] = sbrtbacknormtemp

        # determine spatially uniform background templates
        for i in gdat.indxenerfull:
            for m in gdat.indxdqltfull:
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
        raise Exception('At least one background template must be positive everywhere.')
    
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
                for m in gdat.indxdqlt:
                    meansigc = gmod.psfpexpr[i * gmod.numbpsfptotl + m * gmod.numbpsfptotl * gdat.numbener]
                    stdvsigc = meansigc * 0.1
                    setp_varb(gdat, 'sigcen%02devt%d' % (i, m), mean=meansigc, stdv=stdvsigc, labl=['$\sigma$', ''], scal='gaus', strgmodl=strgmodl)
                    
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
            if gdat.typeexpr == 'gmix':
                minmsigm = 0.01 / gdat.anglfact
                maxmsigm = 0.1 / gdat.anglfact
            if gdat.typeexpr == 'ferm':
                minmsigm = 0.1
                maxmsigm = 10.
            if gdat.typeexpr.startswith('HST_WFC3'):
                minmsigm = 0.01 / gdat.anglfact
                maxmsigm = 0.1 / gdat.anglfact
            if gdat.typeexpr == 'chan':
                minmsigm = 0.1 / gdat.anglfact
                maxmsigm = 2. / gdat.anglfact
            minmgamm = 1.5
            maxmgamm = 20.
            setp_varb(gdat, 'sigc', valu= 0.05/gdat.anglfact, minm=minmsigm, maxm=maxmsigm, labl=['$\sigma_c$', ''], ener='full', dqlt='full', strgmodl=strgmodl, strgstat='this')
            setp_varb(gdat, 'sigt', minm=minmsigm, maxm=maxmsigm, ener='full', dqlt='full', strgmodl=strgmodl)
            setp_varb(gdat, 'gamc', minm=minmgamm, maxm=maxmgamm, ener='full', dqlt='full', strgmodl=strgmodl)
            setp_varb(gdat, 'gamt', minm=minmgamm, maxm=maxmgamm, ener='full', dqlt='full', strgmodl=strgmodl)
            
        setp_varb(gdat, 'psff', minm=0., maxm=1., ener='full', dqlt='full', strgmodl=strgmodl)
 
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
        if gmod.boolelemsbrt[l] and gmod.maxmpara.numbelem[l] > 0:
            if not 'dfnc' in listnameecom:
                listnameecom += ['dfnc']
            if not 'dfncsubt' in listnameecom:
                listnameecom += ['dfncsubt']
    gmod.listnameecomtotl = listnameecom + ['modl']
    
    for c in gmod.indxback:
        setp_varb(gdat, 'cntpback%04d' % c, labl=['$C_{%d}$' % c, ''], minm=1., maxm=100., scal='logt', strgmodl=strgmodl)
    
    gmod.listnamegcom = deepcopy(gmod.listnameecomtotl)
    if gmod.boollens:
        gmod.listnamegcom += ['bgrd']
        if gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
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
    
    if gdat.typeverb > 0:
        if strgmodl == 'true':
            strgtemp = 'true'
        if strgmodl == 'fitt':
            strgtemp = 'fitting'
        print('Building elements for the %s model...' % strgtemp)
    
    # element parameters that correlate with the statistical significance of the element
    gmod.namepara.elemsign = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lght'):
            gmod.namepara.elemsign[l] = 'flux'
        if gmod.typeelem[l] == 'lens':
            gmod.namepara.elemsign[l] = 'defs'
        if gmod.typeelem[l].startswith('clus'):
            gmod.namepara.elemsign[l] = 'nobj'
    
    # define the names and scalings of element parameters
    for l in gmod.indxpopl:
        
        # 'temp' these scale definitions are already done in setp_varb (which allows user-defined values) and redundant here. Should be removed.
        if gmod.typeelem[l].startswith('lghtline'):
            gmod.namepara.genrelem[l] = ['elin']
            gmod.scalpara.genrelem[l] = ['logt']
        elif gmod.typespatdist[l] == 'diskscal':
            gmod.namepara.genrelem[l] = ['xpos', 'ypos']
            gmod.scalpara.genrelem[l] = ['self', 'dexp']
        elif gmod.typespatdist[l] == 'gangexpo':
            gmod.namepara.genrelem[l] = ['gang', 'aang']
            gmod.scalpara.genrelem[l] = ['expo', 'self']
        elif gmod.typespatdist[l] == 'glc3':
            gmod.namepara.genrelem[l] = ['dglc', 'thet', 'phii']
            gmod.scalpara.genrelem[l] = ['powr', 'self', 'self']
        else:
            gmod.namepara.genrelem[l] = ['xpos', 'ypos']
            gmod.scalpara.genrelem[l] = ['self', 'self']
        
        # amplitude
        if gmod.typeelem[l] == 'lghtpntsagnntrue':
            gmod.namepara.genrelem[l] += ['lum0']
            gmod.scalpara.genrelem[l] += ['dpowslopbrek']
        elif gmod.typeelem[l] == 'lghtpntspuls':
            gmod.namepara.genrelem[l] += ['per0']
            gmod.scalpara.genrelem[l] += ['lnormeanstdv']
        elif gmod.typeelem[l].startswith('lght'):
            gmod.namepara.genrelem[l] += ['flux']
            gmod.scalpara.genrelem[l] += [gmod.typeprioflux[l]]
        elif gmod.typeelem[l] == 'lens':
            gmod.namepara.genrelem[l] += ['defs']
            gmod.scalpara.genrelem[l] += ['powr']
        elif gmod.typeelem[l].startswith('clus'):
            gmod.namepara.genrelem[l] += ['nobj']
            gmod.scalpara.genrelem[l] += ['powr']
       
        # shape
        if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
            gmod.namepara.genrelem[l] += ['gwdt']
            gmod.scalpara.genrelem[l] += ['powr']
        if gmod.typeelem[l] == 'lghtlinevoig':
            gmod.namepara.genrelem[l] += ['sigm']
            gmod.scalpara.genrelem[l] += ['logt']
            gmod.namepara.genrelem[l] += ['gamm']
            gmod.scalpara.genrelem[l] += ['logt']
        
        # others
        if gmod.typeelem[l] == 'lghtpntspuls':
            gmod.namepara.genrelem[l] += ['magf']
            gmod.scalpara.genrelem[l] += ['lnormeanstdv']
            gmod.namepara.genrelem[l] += ['geff']
            gmod.scalpara.genrelem[l] += ['self']
        elif gmod.typeelem[l] == 'lghtpntsagnntrue':
            gmod.namepara.genrelem[l] += ['dlos']
            gmod.scalpara.genrelem[l] += ['powr']

        if gdat.numbener > 1 and gmod.typeelem[l].startswith('lghtpnts'):
            if gmod.spectype[l] == 'colr':
                for i in gdat.indxener:
                    if i == 0:
                        continue
                    gmod.namepara.genrelem[l] += ['sindcolr%04d' % i]
                    gmod.scalpara.genrelem[l] += ['self']
            else:
                gmod.namepara.genrelem[l] += ['sind']
                gmod.scalpara.genrelem[l] += ['self']
                if gmod.spectype[l] == 'curv':
                    gmod.namepara.genrelem[l] += ['curv']
                    gmod.scalpara.genrelem[l] += ['self']
                if gmod.spectype[l] == 'expc':
                    gmod.namepara.genrelem[l] += ['expc']
                    gmod.scalpara.genrelem[l] += ['self']
        if gmod.typeelem[l] == 'lens':
            if gdat.variasca:
                gmod.namepara.genrelem[l] += ['asca']
                gmod.scalpara.genrelem[l] += ['self']
            if gdat.variacut:
                gmod.namepara.genrelem[l] += ['acut']
                gmod.scalpara.genrelem[l] += ['self']
    
    # names of element parameters for each scaling
    gmod.namepara.genrelemscal = [{} for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        for scaltype in gdat.listscaltype:
            gmod.namepara.genrelemscal[l][scaltype] = []
            for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                if scaltype == gmod.scalpara.genrelem[l][k]:
                    gmod.namepara.genrelemscal[l][scaltype].append(nameparagenrelem)

    # variables for which whose marginal distribution and pair-correlations will be plotted
    for l in gmod.indxpopl:
        gmod.namepara.derielemodim[l] = ['deltllik']
        
        if gmod.typeelem[l] == 'lens':
            gmod.namepara.derielemodim[l] += ['deflprof']
        
        if gdat.boolbindspat:
            # to be deleted
            #if not 'xpos' in gmod.namepara.derielemodim[l]:
            #    gmod.namepara.derielemodim[l] += ['xpos']
            #if not 'ypos' in gmod.namepara.derielemodim[l]:
            #    gmod.namepara.derielemodim[l] += ['ypos']
            
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
            gmod.namepara.derielemodim[l] += ['mcut', 'distsour', 'rele']#, 'reln', 'relk', 'relf', 'relm', 'reld', 'relc']
    
        #for k in range(len(gmod.namepara.derielemodim[l])):
        #    gmod.namepara.derielemodim[l][k] += 'pop%d' % l
        
        # check later
        # temp
        #if strgmodl == 'fitt':
        #    for q in gdat.indxrefr: 
        #        if gmod.nameparagenrelemampl[l] in gdat.refr.namepara.elem[q]:
        #            gmod.namepara.derielemodim[l].append('aerr' + gdat.listnamerefr[q])
    
    # derived parameters
    #gmod.listnameparaderitotl = [temptemp for temp in gmod.namepara.deri.elem for temptemp in temp]
    ##gmod.listnameparaderitotl += gmod.namepara.scal
    #
    #for namediff in gmod.listnamediff:
    #    gmod.listnameparaderitotl += ['cntp' + namediff]
    #
    #if gdat.typeverb > 1:
    #    print('gmod.listnameparaderitotl')
    #    print(gmod.listnameparaderitotl)

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

    gmod.namepara.elemodim = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        gmod.namepara.elemodim[l] += gmod.namepara.genrelem[l] + gmod.namepara.derielemodim[l]
    
    if gdat.typeverb > 1:
        print('gmod.namepara.derielemodim')
        print(gmod.namepara.derielemodim)
        print('gmod.namepara.elemodim')
        print(gmod.namepara.elemodim)
        
    if gdat.booldiag:
        if np.unique(np.array(gmod.namepara.derielemodim[l])).size != len(gmod.namepara.derielemodim[l]):
            print('gmod.namepara.derielemodim')
            print(gmod.namepara.derielemodim)
            raise Exception('')
        
        for name in gmod.namepara.derielemodim[l]:
            if name in gmod.namepara.genrelem[l]:
                print('gmod.namepara.derielemodim')
                print(gmod.namepara.derielemodim)
                print('gmod.namepara.genrelem')
                print(gmod.namepara.genrelem)
                raise Exception('')
    
    # defaults
    gmod.liststrgpdfnmodu = [[] for l in gmod.indxpopl]
    gmod.namepara.genrelemmodu = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        if gmod.typeelem[l].startswith('lght'): 
            if gdat.typeexpr == 'ferm' and gdat.xposcntr == 0.:
                if l == 1:
                    gmod.liststrgpdfnmodu[l] += ['tmplnfwp']
                    gmod.namepara.genrelemmodu[l] += ['xposypos']
                if l == 2:
                    gmod.liststrgpdfnmodu[l] += ['tmplnfwp']
                    gmod.namepara.genrelemmodu[l] += ['xposypos']
    
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
    
    #gmod.namepara.genr.elemeval = [[] for l in gmod.indxpopl]
    #for l in gmod.indxpopl:
    #    if gmod.typeelem[l].startswith('clus'):
    #        gmod.namepara.genr.elemeval[l] = ['xpos', 'ypos', 'nobj']
    #    if gmod.typeelem[l] == 'clusvari':
    #        gmod.namepara.genr.elemeval[l] += ['gwdt']
    #    if gmod.typeelem[l] == 'lens':
    #        gmod.namepara.genr.elemeval[l] = ['xpos', 'ypos', 'defs', 'asca', 'acut']
    #    if gmod.typeelem[l].startswith('lghtline'):
    #        gmod.namepara.genr.elemeval[l] = ['elin', 'spec']
    #    elif gmod.typeelem[l] == 'lghtgausbgrd':
    #        gmod.namepara.genr.elemeval[l] = ['xpos', 'ypos', 'gwdt', 'spec']
    #    elif gmod.typeelem[l].startswith('lght'):
    #        gmod.namepara.genr.elemeval[l] = ['xpos', 'ypos', 'spec']
    
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
    
    for l in gmod.indxpopl:
        setp_varb(gdat, 'slopprionobjpop%d' % l, labl=['$\alpha_{%d}$' % l, ''], strgmodl=strgmodl, popl=l)

    if strgmodl == 'true':
        gmod.indxpoplassc = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            if gmod.numbpopl == 3 and gmod.typeelem[1] == 'lens':
                gmod.indxpoplassc[l] = [l]
            else:
                gmod.indxpoplassc[l] = gmod.indxpopl

    # number of element parameters
    if gmod.numbpopl > 0:
        gmod.numbparagenrelemsing = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparagenrelempopl = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparagenrelemcuml = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparagenrelemcumr = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparaderielempopl = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparaderielemsing = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparaelempopl = np.zeros(gmod.numbpopl, dtype=int)
        gmod.numbparaelemsing = np.zeros(gmod.numbpopl, dtype=int)
        for l in gmod.indxpopl:
            # number of generative element parameters for a single element of a specific population
            gmod.numbparagenrelemsing[l] = len(gmod.namepara.genrelem[l])
            # number of derived element parameters for a single element of a specific population
            gmod.numbparaderielemsing[l] = len(gmod.namepara.derielemodim[l])
            # number of element parameters for a single element of a specific population
            gmod.numbparaelemsing[l] = len(gmod.namepara.elem[l])
            # number of generative element parameters for all elements of a specific population
            gmod.numbparagenrelempopl[l] = gmod.numbparagenrelemsing[l] * gmod.maxmpara.numbelem[l]
            # number of generative element parameters up to the beginning of a population
            gmod.numbparagenrelemcuml[l] = np.sum(gmod.numbparagenrelempopl[:l])
            # number of generative element parameters up to the end of a population
            gmod.numbparagenrelemcumr[l] = np.sum(gmod.numbparagenrelempopl[:l+1])
            # number of derived element parameters for all elements of a specific population
            gmod.numbparaderielempopl[l] = gmod.numbparaderielemsing[l] * gmod.maxmpara.numbelem[l]
            # number of element parameters for all elements of a specific population
            gmod.numbparaelempopl[l] = gmod.numbparaelemsing[l] * gmod.maxmpara.numbelem[l]
        # number of generative element parameters summed over all populations
        gmod.numbparagenrelem = np.sum(gmod.numbparagenrelempopl)
        # number of derived element parameters summed over all populations
        gmod.numbparaderielem = np.sum(gmod.numbparaderielempopl)
        # number of element parameters summed over all populations
        gmod.numbparaelemtotl = np.sum(gmod.numbparaelempopl)
    
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
        gmod.numblpri = 1
        #+ gmod.numbparagenrelempopl * gmod.numbpopl
    else:
        gmod.numblpri = 0
    if gdat.boolpenalpridiff:
        gmod.numblpri += 1
    indxlpri = np.arange(gmod.numblpri)

    # append the population tags to element parameter names
    #for l in gmod.indxpopl:
    #    gmod.namepara.genrelem[l] = [gmod.namepara.genrelem[l][g] + 'pop%d' % l for g in gmod.indxparagenrelemsing[l]]
    
    # Boolean flag indicating if a generative element parameter is "positional"
    gmod.boolcompposi = [[] for l in gmod.indxpopl]
    for l in gmod.indxpopl:
        gmod.boolcompposi[l] = np.zeros(gmod.numbparagenrelemsing[l], dtype=bool)
        if gmod.typeelem[l].startswith('lghtline'):
            gmod.boolcompposi[l][0] = True
        else:
            gmod.boolcompposi[l][0] = True
            gmod.boolcompposi[l][1] = True
    
    # flattened list of element parameters
    for strgtypepara in ['genrelem', 'derielemodim', 'elemodim']:
        nameparaelem = getattr(gmod.namepara, strgtypepara)
        setattr(gmod.namepara, strgtypepara + 'flat', [])
        nameparaelemflat = getattr(gmod.namepara, strgtypepara + 'flat')
        for l in gmod.indxpopl:
            for nameparaelemodim in nameparaelem[l]:
                nameparaelemflat.append(nameparaelemodim + 'pop%d' % l)
    
    #gmod.numbparaelem = np.empty(gmod.numbpopl, dtype=int)
    #for l in gmod.indxpopl:
    #    gmod.numbparaelem[l] = len(gmod.namepara.elem[l])
    
    gmod.numbdeflsingplot = gdat.numbdeflsubhplot
    if gmod.numbpopl > 0:
        gmod.numbdeflsingplot += 3

    gmod.convdiffanyy = True in convdiff

    cntr = tdpy.cntr()
    
    # collect lensing here
    if gmod.boollens:
        if gmod.boollenshost:
            setp_varb(gdat, 'redshost', valu=0.2, minm=0., maxm=0.4, strgmodl=strgmodl)
            setp_varb(gdat, 'redssour', valu=1., minm=0.5, maxm=1.5, strgmodl=strgmodl)

        gmod.adislens = gdat.adisobjt(gmod.redshost)
        gmod.adissour = gdat.adisobjt(gmod.redssour)
        gmod.adislenssour = gmod.adissour - (1. + gmod.redshost) / (1. + gmod.redssour) * gmod.adislens
        gmod.ratimassbeinsqrd = chalcedon.retr_ratimassbeinsqrd(gmod.adissour, gmod.adislens, gmod.adislenssour)
        gmod.mdencrit = chalcedon.retr_mdencrit(gmod.adissour, gmod.adislens, gmod.adislenssour)
    
        gmod.bctrpara.adislens = gdat.adisobjt(gmod.bctrpara.redshost)
        
    
    # base parameter indices (indxpara) are being defined here
    # define parameter indices
    print('Defining the indices of individual parameters...')
    if gmod.numbpopl > 0:

        # number of elements
        for l in gmod.indxpopl:
            indx = cntr.incr()
            setattr(gmod.indxpara, 'numbelempop%d' % l, indx)
        
        # hyperparameters
        ## mean number of elements
        if gmod.typemodltran == 'pois':
            #gmod.indxpara.meanelem = np.empty(gmod.numbpopl, dtype=int)
            for l in gmod.indxpopl:
                if gmod.maxmpara.numbelem[l] > 0:
                    indx = cntr.incr()
                    print('meanelempop%d % l')
                    print('meanelempop%d' % l)
                    print('indx')
                    print(indx)
                    setattr(gmod.indxpara, 'meanelempop%d' % l, indx)
                    #gmod.indxpara.meanelem[l] = indx

        ## parameters parametrizing priors on element parameters
        liststrgvarb = []
        for l in gmod.indxpopl:
            if gmod.maxmpara.numbelem[l] > 0:
                for strgpdfnelemgenr, strgfeat in zip(gmod.scalpara.genrelem[l], gmod.namepara.genrelem[l]):
                    liststrgvarb = []
                    if strgpdfnelemgenr == 'expo' or strgpdfnelemgenr == 'dexp':
                        liststrgvarb += [strgfeat + 'distscal']
                    if strgpdfnelemgenr == 'powr':
                        liststrgvarb += ['slopprio' + strgfeat]
                    if strgpdfnelemgenr == 'dpow':
                        liststrgvarb += [strgfeat + 'distbrek', strgfeat + 'sloplowr', strgfeat + 'slopuppr']
                    if strgpdfnelemgenr == 'gausmean' or strgpdfnelemgenr == 'lnormean':
                        liststrgvarb += [strgfeat + 'distmean']
                    if strgpdfnelemgenr == 'gausstdv' or strgpdfnelemgenr == 'lnorstdv':
                        liststrgvarb += [strgfeat + 'diststdv']
                    if strgpdfnelemgenr == 'gausmeanstdv' or strgpdfnelemgenr == 'lnormeanstdv':
                        liststrgvarb += [strgfeat + 'distmean', strgfeat + 'diststdv']
                    for strgvarb in liststrgvarb:
                        indx = cntr.incr()
                        setattr(gmod.indxpara, strgvarb + 'pop%d' % l, indx)
            #setattr(gmod.indxpara, strgvarb, np.zeros(gmod.numbpopl, dtype=int) - 1)

        # can potentially be used to define arrays of parameter indices for a given hyperparameter across all populations
        ## this is where 'slopprionobj'-like definitions happen
        ## turned off
        #for l in gmod.indxpopl:
        #    if gmod.maxmpara.numbelem[l] > 0:
        #        for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
        #            
        #            if gmod.scalpara.genrelem[l][k] == 'self':
        #                continue
        #            indx = cntr.incr()

        #            if gmod.scalpara.genrelem[l][k] == 'dpow':
        #                for nametemp in ['brek', 'sloplowr', 'slopuppr']:
        #                    strg = '%s' % nametemp + nameparagenrelem
        #                    setattr(gmod.indxpara, strg, indx)
        #                    setattr(gmod.indxpara, strg, indx)
        #            else:
        #                if gmod.scalpara.genrelem[l][k] == 'expo' or gmod.scalpara.genrelem[l][k] == 'dexp':
        #                    strghypr = 'scal'
        #                if gmod.scalpara.genrelem[l][k] == 'powr':
        #                    strghypr = 'slop'
        #                if gmod.scalpara.genrelem[l][k] == 'gausmean' or gmod.scalpara.genrelem[l][k] == 'gausmeanstdv' or \
        #                                gmod.scalpara.genrelem[l][k] == 'lnormean' or gmod.scalpara.genrelem[l][k] == 'lnormeanstdv':
        #                    strghypr = 'mean'
        #                if gmod.scalpara.genrelem[l][k] == 'gausstdv' or gmod.scalpara.genrelem[l][k] == 'gausmeanstdv' or \
        #                                gmod.scalpara.genrelem[l][k] == 'lnorstdv' or gmod.scalpara.genrelem[l][k] == 'lnormeanstdv':
        #                    strghypr = 'stdv'
        #                strg = strghypr + 'prio' + nameparagenrelem
        #                print('strg')
        #                print(strg)
        #                setattr(gmod.indxpara, strg, indx)
        
        #raise Exception('')
    
    # parameter groups 
    ## number element elements for all populations
    if gmod.numbpopl > 0:
        gmod.indxpara.numbelem = np.empty(gmod.numbpopl, dtype=int)
        for l in gmod.indxpopl:
            gmod.indxpara.numbelem[l] = indx
    
    ## PSF parameters
    if gmod.typeevalpsfn != 'none':
        for m in gdat.indxdqlt:
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

        gmod.numbpsfptotldqlt = gdat.numbdqlt * gmod.numbpsfptotl
        gmod.numbpsfptotlener = gdat.numbener * gmod.numbpsfptotl
        numbpsfp = gmod.numbpsfptotl * gdat.numbener * gdat.numbdqlt
        indxpsfpform = np.arange(numbpsfpform)
        indxpsfptotl = np.arange(gmod.numbpsfptotl)
   
        gmod.indxpara.psfp = np.sort(gmod.indxpara.psfp)

    ## background parameters
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
        gmod.indxpara.xpossour = cntr.incr()
        gmod.indxpara.ypossour = cntr.incr()
        gmod.indxpara.fluxsour = cntr.incr()
        if gdat.numbener > 1:
            gmod.indxpara.sindsour = cntr.incr()
        gmod.indxpara.sizesour = cntr.incr()
        gmod.indxpara.ellpsour = cntr.incr()
        gmod.indxpara.anglsour = cntr.incr()
    if gmod.typeemishost != 'none' or gmod.boollens:
        for e in gmod.indxsersfgrd: 
            if gmod.typeemishost != 'none':
                setattr(gmod.indxpara, 'xposhostisf%d' % e, cntr.incr())
                setattr(gmod.indxpara, 'yposhostisf%d' % e, cntr.incr())
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
    if gdat.typeexpr.startswith('HST_WFC3'):
        gmod.listnamecomplens = ['hostlght', 'hostlens', 'sour', 'extr']
        for namecomplens in gmod.listnamecomplens:
            setattr(gmod, 'liststrg' + namecomplens, [])
            setattr(gmod.indxpara, namecomplens, [])
    if gmod.boollens or gmod.typeemishost != 'none':
        gmod.liststrghostlght += ['xposhost', 'yposhost', 'ellphost', 'anglhost']
        gmod.liststrghostlens += ['xposhost', 'yposhost', 'ellphost', 'anglhost']
    if gmod.typeemishost != 'none':
        gmod.liststrghostlght += ['fluxhost', 'sizehost', 'serihost']
        if gdat.numbener > 1:
            gmod.liststrghostlght += ['sindhost']
    if gmod.boollens:
        gmod.liststrghostlens += ['beinhost']
        gmod.liststrgextr += ['sherextr', 'sangextr']
        gmod.liststrgsour += ['xpossour', 'ypossour', 'fluxsour', 'sizesour', 'ellpsour', 'anglsour']
        if gdat.numbener > 1:
            gmod.liststrgsour += ['sindsour']
    
    for strg, valu in gmod.__dict__.items():
        
        if isinstance(valu, list) or isinstance(valu, np.ndarray):
            continue
        
        if gdat.typeexpr.startswith('HST_WFC3'):
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
        gmod.indxpara.host = gmod.indxpara.hostlght + gmod.indxpara.hostlens
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
        


def retr_strgcalc(gdat):

    ## list of strings indicating different methods of calculating the subhalo mass fraction
    gdat.liststrgcalcmasssubh = ['delt', 'intg']
        


def plot_lens(gdat):
    
    if gmod.boolelemdeflsubh:
        xdat = gdat.blimpara.angl[1:] * gdat.anglfact
        lablxdat = gdat.labltotlpara.gang
        
        listdeflscal = np.array([4e-2, 4e-2, 4e-2]) / gdat.anglfact
        listanglscal = np.array([0.05, 0.1, 0.05]) / gdat.anglfact
        listanglcutf = np.array([1.,    1.,  10.]) / gdat.anglfact
        listasym = [False, False, False]
        listydat = []
        for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
            listydat.append(chalcedon.retr_deflcutf(gdat.blimpara.angl[1:], deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
        
        for scalxdat in ['self', 'logt']:
            path = gdat.pathinitintr + 'deflcutf' + scalxdat + '.%s' % gdat.typefileplot
            tdpy.plot_gene(path, xdat, listydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, \
                                                             lablydat=r'$\alpha_n$ [$^{\prime\prime}$]', limtydat=[1e-3, 1.5e-2], limtxdat=[None, 2.])
       
        # pixel-convoltuion of the Sersic profile
        # temp -- y axis labels are wrong, should be per solid angle
        xdat = gdat.blimpara.xpossers * gdat.anglfact
        for n in range(gdat.numbindxsers + 1):
            for k in range(gdat.numbhalfsers + 1):
                if k != 5:
                    continue
                path = gdat.pathinitintr + 'sersprofconv%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, gdat.sersprof[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                #path = gdat.pathinitintr + 'sersprofcntr%04d%04d.%s' % (n, k, gdat.typefileplot)
                #tdpy.plot_gene(path, xdat, gdat.sersprofcntr[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], scalxdat='logt', \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
       
        xdat = gdat.blimpara.angl * gdat.anglfact
        listspec = np.array([1e-19, 1e-18, 1e-18, 1e-18]) / gdat.anglfact
        listsize = np.array([0.3, 1., 1., 1.]) / gdat.anglfact
        listindx = np.array([4., 2., 4., 10.])
        listydat = []
        listlabl = []
        for spec, size, indx in zip(listspec, listsize, listindx):
            listydat.append(spec * retr_sbrtsersnorm(gdat.blimpara.angl, size, indxsers=indx))
            listlabl.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
        path = gdat.pathinitintr + 'sersprof.%s' % gdat.typefileplot
        tdpy.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, \
                                                                                                   listlegd=listlegd, listhlin=1e-7, limtydat=[1e-8, 1e0])
    
        #valulevl = np.linspace(7.5, 9., 5)
        valulevl = [7.0, 7.3, 7.7, 8., 8.6]
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        cont = axis.contour(gdat.blimpara.redshost, gdat.blimpara.redssour, minmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=20, fmt='%.3g')
        axis.set_xlabel(r'$z_{\rm{hst}}$')
        axis.set_ylabel(r'$z_{\rm{src}}$')
        axis.set_title(r'$M_{c,min}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsminm.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        valulevl = np.linspace(9., 11., 20)
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=9., vmax=11.)
        cont = axis.contour(gdat.blimpara.redshost, gdat.blimpara.redssour, maxmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=15, fmt='%.3g')
        axis.set_xlabel('$z_{hst}$')
        axis.set_ylabel('$z_{src}$')
        axis.set_title(r'$M_{c,max}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsmaxm.%s' % gdat.typefileplot
        plt.colorbar(imag) 
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adislens * gdat.sizepixl * 1e-3)
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adislens * 2. * gdat.maxmgangdata * 1e-3)
        axis.set_xlabel('$z_h$')
        axis.set_yscale('log')
        axis.set_ylabel(r'$\lambda$ [kpc]')
        path = gdat.pathinitintr + 'wlenreds.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        fracacutasca = np.logspace(-1., 2., 20)
        mcut = aspendos.retr_mcutfrommscl(fracacutasca)
        axis.lognp.log(fracacutasca, mcut)
        axis.set_xlabel(r'$\tau_n$')
        axis.set_ylabel(r'$M_{c,n} / M_{0,n}$')
        axis.axhline(1., ls='--')
        path = gdat.pathinitintr + 'mcut.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
       

def setup_pcat(gdat):
    '''
    Routines common among sample() and wrappers()
    These should be made feedable to sample()
    '''
    
    # folders
    if gdat.pathbase is None:
        gdat.pathbase = os.environ["PCAT_DATA_PATH"] + '/'
        
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathvisu = gdat.pathbase + 'visuals/'
    
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathinpt = gdat.pathdata + 'inpt/'
        
    # run tag
    # list of parameter features to be turned into lists
    gdat.liststrgfeatparalist = ['minm', 'maxm', 'scal', 'lablroot', 'lablunit', 'labl', 'labltotl', 'name', 'mean', 'stdv']
    
    # list of parameter features
    gdat.liststrgfeatpara = gdat.liststrgfeatparalist + ['limt', 'numbbins', 'blim', 'delt', 'indxbins', 'cmap', 'bctr', 'tick', 'numbbins', 'valutickmajr', \
                                                                                                                # index of the parameter in the parameter vector
                                                                                                                'indx', \

                                                                                                                # to be deleted
                                                                                                                #'numb', \
                                                                                                            'labltickmajr', 'valutickminr', 'labltickminr']
    # list of scalings
    gdat.listscaltype = ['self', 'logt', 'atan', 'gaus', 'pois', 'expo', 'drct']
        
    # number of standard deviations around Gaussian mean used to plot variables
    gdat.numbstdvgaus = 4.
    

def setup_pcat_model(gdat):
    '''
    Routines that act on the models common among sample() and wrappers()
    These should be made feedable to sample()
    '''
    
    #setp_varb(gdat, 'colr', valu='mediumseagreen', strgmodl='refr')
    setp_varb(gdat, 'colr', valu='b', strgmodl='fitt')
    if gdat.typedata == 'simu':
        setp_varb(gdat, 'colr', valu='g', strgmodl='true')
        

def init_image( \
              
         ## type of experiment
         ### 'gmix': two-dimensional Gaussian Mixture Model (GMM)
         typeexpr='gmix', \
        
         # Boolean flag to bin the data in energy
         boolbinsener=None, \
        
         dicttrue=None, \
         dictfitt=None, \

         indxdqltincl=None, \
         indxenerincl=None, \
        
         listmask=None, \
        
         anlytype=None, \
         
         strgexprsbrt=None, \
         
         # type of mask for the exposure map
         typemaskexpo='ignr', \
         
         # type of exposure
         ## 'cons': constant
         ## 'file': provided in a file
         typeexpo='cons', \

         # modes of operation
         ## only generate and plot simulated data
         boolsimuonly=False, \
         
         ### Sersic type
         typesers='vauc', \
         ## transdimensional parameters (elements)
         ### vary projected scale radius
         variasca=True, \
         ### vary projected cutoff radius
         variacut=True, \

         # maximum spatial distance out to which element kernel will be evaluated
         maxmangleval=None, \
         
         ## exposure map
         expo=None, \
        
         strgcnfgsimu=None, \
        
         # initial state
         initpsfprefr=False, \
         initpsfp=None, \
        
         # diagnostics
         ## Boolean to turn on diagnostic mode
         booldiag=True, \
         
         # Boolean flag to make the frame plots only for the central energy and PSF bin
         boolmakeframcent=True, \
         
         # strings for energy bins
         strgenerfull=None, \
         
         strgexprname=None, \
         
         strganglunit=None, \
         
         strganglunittext=None, \
         
         anglfact=None, \
            
         limtydathistfeat=None, \
            
         ## optional deterministic seed for sampling element parameters
         typeseedelem=None, \
         
         # model
         # emission
         ## elements

         ## Boolean to thin down the data
         boolthindata=False, \
         
         ## PSF
         specfraceval=None, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=100, \
    
         listprefsbrtsbrt=None, \
         listprefsbrtener=None, \
         listprefsbrtlabltotl=None, \

         lablgangunit=None, \
         lablxpos=None, \
         lablypos=None, \
         lablfluxunit=None, \
         lablflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxdqltfull=None, \
         binsenerfull=None, \
         asymfluxprop=False, \
         
         ## Boolean flag to make the PSF model informed
         boolpriopsfninfo=False, \
         
         ## spectral

         # temp
         margfactmodl=1., \
         maxmgangdata=None, \
        
         # proposals
         stdvprophypr=0.01, \
         stdvproppsfp=0.1, \
         stdvpropbacp=0.01, \
         stdvproplenp=1e-4, \
         stdvxpos=0.001, \
         stdvypos=0.001, \
         stdvflux=0.001, \
         stdvspep=0.001, \
         stdvspmrsind=0.2, \
         varistdvlbhl=True, \
        
         # hyperparameters
         fittampldisttype=None, \
         # metamodel settings
         
         ## PSF evaluation type
         ## kernel evaluation type
         kernevaltype='ulip', \

         # photometric model
        
         # random state
         ## seed for numpy random number generator
         typeseed=0, \
         
         boolpenalpridiff=False, \
         loadvaripara=False, \
         
         anglassc=None, \
         nameexpr=None, \
         
         xposprio=None, \
         yposprio=None, \
         minmcntpdata=None, \
         strgexpo=None, \
         
         # likelihood function
         liketype='pois', \
         
         xposcntr=0., \
         yposcntr=0., \
         
         maxmangl=None, \
         
         ## Boolean flag to make the frame plots short
         boolmakeshrtfram=False, \
        
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

         boolmakeplot=True, \
         boolmakeplotinit=True, \
         boolmakeplotfram=True, \
         boolmakeplotfinlprio=True, \
         boolmakeplotfinlpost=True, \
         
         ## Boolean to overplot the elements
         boolplotelem=True, \
         
         boolmakeplotintr=False, \
         scalmaps='asnh', \
         
         ## file type of the plot
         typefileplot='png', \
         
         
         # arguments common among sample() and wrappers
         
         # modes of operation
         ## perform an additional run sampling from the prior
         checprio=False, \

         # name of the configuration
         strgcnfg=None, \
        
         # elements
         ## reference element catalog
         ## Boolean flag to associate elements to those in reference element catalogs
         boolasscrefr=None, \
         
         # data
         ## type of data
         ### 'simu': simulated data
         ### 'inpt': input data
         ### 'real': real data retrieved from databases
         typedata=None, \
         
         # user interaction
         ## type of verbosity
         typeverb=1, \

         ## base path where visuals and data will be written
         pathbase=None, \
        
         ## transdimensional proposal probabilities
         probtran=None, \
         probspmr=None, \
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
         # standard deviation of the Gaussian from which the angular splitting will be drawn for splits and merges
         radispmr=None, \
         
         # factor by which to multiply the prior for each additional degree of freedom in the model
         factpriodoff=None, \

         # Boolean flag to force the posterior towards the reference
         boolrefeforc=False, \
         # index of the reference catalog towards which the posterior will be forced
         indxrefrforc=None, \

         **dictpcat

              ):
   
    # preliminary setup
    # construct the global object 
    gdat = tdpy.gdatstrt()
    
    # copy locals (inputs) to the global object
    dictinpt = dict(locals())
    for attr, valu in dictinpt.items():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)

    print('gdat.typeexpr')
    print(gdat.typeexpr)
        
    dictpcatinpt = retr_dictpcatinpt()
    
    listlablcnfg = ['', '', '', '', '']
    listnamecnfgextn = ['truevlow', 'trueloww', 'nomi', 'truehigh', 'truevhig']
    dictpcatinptvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictpcatinptvari[namecnfgextn] = retr_dictpcatinpt()
    
    # defaults
    if gdat.typedata is None:
        if gdat.strgexprsbrt is None:
            gdat.typedata = 'simu'
        else:
            gdat.typedata = 'inpt'
    print('gdat.typedata')
    print(gdat.typedata)

    # list of models
    gdat.liststrgmodl = []
    if gdat.typedata == 'simu':
        gdat.liststrgmodl += ['true']
    gdat.liststrgmodl += ['fitt']
    
    for strgmodl in gdat.liststrgmodl:
        for strgstat in ['this', 'next']:
            setattr(gdat, strgstat, tdpy.gdatstrt())
        
    setup_pcat(gdat)

    for strgmodl in gdat.liststrgmodl + ['refr']:
        setattr(gdat, strgmodl, tdpy.gdatstrt())
        gmod = getattr(gdat, strgmodl)
        for strgfeatpara in gdat.liststrgfeatpara:
            setattr(gmod, strgfeatpara + 'para', tdpy.gdatstrt())
    
    # copy user-provided inputs into gmod
    for strgfeatpara in gdat.liststrgfeatpara:
        setattr(gdat, strgfeatpara + 'para', tdpy.gdatstrt())
        
    if gdat.dicttrue is not None:
        for name in gdat.dicttrue:
            setattr(gdat.true, name, gdat.dicttrue[name])
    if gdat.dictfitt is not None:
        for name in gdat.dictfitt:
            setattr(gdat.fitt, name, gdat.dictfitt[name])
    
    ## name of the configuration function
    if gdat.strgcnfg is None:
        gdat.strgcnfg = inspect.stack()[1][3]
   
    if gdat.typedata == 'inpt' and gdat.strgcnfgsimu is not None:
        print('Will use %s to account for selection effects.' % gdat.strgcnfgsimu)
        gdat.pathoutpcnfgsimu = retr_pathoutpcnfg(gdat.pathbase, gdat.strgcnfgsimu)

    if gdat.anlytype is None:
        if gdat.typeexpr == 'chan':
            gdat.anlytype = 'home'
        elif gdat.typeexpr == 'ferm':
            gdat.anlytype = 'rec8pnts'
        else:
            gdat.anlytype = 'nomi'
    
    # experiment defaults
    if gdat.typeexpr == 'ferm':
        gdat.lablenerunit = 'GeV'
    if gdat.typeexpr == 'chan':
        gdat.lablenerunit = 'keV'
    if gdat.typeexpr == 'gmix':
        gdat.lablenerunit = ''
    if gdat.typeexpr == 'fire':
        gdat.lablenerunit = '$\mu$m^{-1}'
    
    if gdat.typeexpr == 'ferm':
        if gdat.anlytype[4:8] == 'pnts':
            blim = np.logspace(np.log10(0.3), np.log10(10.), 4)
        if gdat.anlytype[4:8] == 'back':
            blim = np.logspace(np.log10(0.3), np.log10(300.), 31)
    elif gdat.typeexpr == 'chan':
        if gdat.anlytype.startswith('home'):
            blim = np.array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
        if gdat.anlytype.startswith('extr'):
            blim = np.array([0.5, 2., 8.])
        if gdat.anlytype.startswith('spec'):
            blim = np.logspace(np.log10(0.5), np.log10(10.), 21)
    elif gdat.typeexpr == 'fire':
        blim = np.logspace(np.log10(1. / 2.5e-6), np.log10(1. / 0.8e-6), 31)
    elif gdat.typeexpr.startswith('HST_WFC3'):
        # temp
        #blim = np.array([500., 750, 1000.])
        blim = np.array([750, 1000.])
    else:
        print('')
        print('')
        print('')
        print('gdat.typeexpr')
        print(gdat.typeexpr)
        raise Exception('gdat.typeexpr is not defined.')

    if gdat.typeexpr != 'gmix':
        setp_varb(gdat, 'enerfull', blim=blim)
    
    setp_varb(gdat, 'numbpixl', labl=['$N_{pix}$', ''])
    
    if gdat.expo is not None:
        setp_varb(gdat, 'expo', minm=np.amin(gdat.expo), maxm=np.amax(gdat.expo), labl=['$\epsilon$', ''], cmap='OrRd', scal='logt')
    
    # string indicating the energy band
    if gdat.strgenerfull is None:
        if gdat.typeexpr == 'tess':
            gdat.strgenerfull = ['T']
        if gdat.typeexpr == 'sdss':
            gdat.strgenerfull = ['z-band', 'i-band', 'r-band', 'g-band', 'u-band']
        if gdat.typeexpr.startswith('HST_WFC3'):
            #gdat.strgenerfull = ['F606W', 'F814W']
            gdat.strgenerfull = ['F814W']
        if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'fire': 
            gdat.strgenerfull = []
            for i in range(len(gdat.blimpara.enerfull) - 1):
                gdat.strgenerfull.append('%.3g %s - %.3g %s' % (gdat.blimpara.enerfull[i], gdat.lablenerunit, gdat.blimpara.enerfull[i+1], gdat.lablenerunit))
        if gdat.typeexpr == 'gmix':
            gdat.strgenerfull = ['']
    
    ## PSF class
    if gdat.indxdqltfull is None:
        if gdat.typeexpr == 'ferm':
            gdat.indxdqltfull = np.arange(2)
        else:
            gdat.indxdqltfull = np.arange(1)
    
    if gdat.indxdqltincl is None:
        if gdat.typeexpr == 'ferm':
            gdat.indxdqltincl = np.array([0, 1])
        else:
            gdat.indxdqltincl = np.arange(1)
    
    if gdat.indxdqltincl is not None:
        gdat.boolbinddqlt = True
    else:
        gdat.boolbinddqlt = False
    
    # number of data quality bins
    if gdat.boolbinddqlt:
        gdat.numbdqlt = gdat.indxdqltincl.size
        gdat.numbdqltfull = gdat.indxdqltfull.size
    else:
        gdat.numbdqlt = 1
        gdat.numbdqltfull = 1
        gdat.indxdqltincl = np.array([0])
    gdat.indxdqlt = np.arange(gdat.numbdqlt)
    
    # Boolean flag to indicate that the data are binned in energy
    if gdat.boolbinsener is None:
        if gdat.typeexpr == 'gmix':
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
        gdat.bctrpara.enerfull = np.sqrt(gdat.blimpara.enerfull[1:] * gdat.blimpara.enerfull[:-1])
    
    ## Boolean flag to vary the PSF
    setp_varb(gdat, 'boolmodipsfn', valu=False, strgmodl='fitt')
    
    # default values for model types
    print('Starting to determine the default values for model types using setp_varbvalu()...')
    if gdat.typeexpr.startswith('HST_WFC3'):
        typeemishost = 'sers'
    else:
        typeemishost = 'none'
    setp_varb(gdat, 'typeemishost', valu=typeemishost)

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
            typebacktemp = interp(gdat.bctrpara.enerfull, meanenerparttemp, sbrtparttemp)
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
    
    if gdat.typeexpr.startswith('HST_WFC3'):
        typeback = [1.]
    if gdat.typeexpr == 'tess':
        typeback = [1.]
    if gdat.typeexpr == 'gmix':
        typeback = [1.]
    if gdat.typeexpr == 'fire':
        typeback = [1.]
    
    setp_varb(gdat, 'typeback', valu=typeback)
    
    if gdat.typeexpr.startswith('HST_WFC3'):
        numbsersfgrd = 1
    else:
        numbsersfgrd = 0
    setp_varb(gdat, 'numbsersfgrd', valu=numbsersfgrd)
    
    if gdat.typeexpr == 'gmix':
        typeelem = ['clus']
    if gdat.typeexpr == 'ferm':
        typeelem = ['lghtpnts']
    if gdat.typeexpr == 'tess':
        typeelem = ['lghtpnts']
    if gdat.typeexpr == 'chan':
        typeelem = ['lghtpnts']
    if gdat.typeexpr.startswith('HST_WFC3'):
        typeelem = ['lens']
        #typeelem = ['lghtpnts', 'lens', 'lghtgausbgrd']
    if gdat.typeexpr == 'fire':
        typeelem = ['lghtlineabso']
    setp_varb(gdat, 'typeelem', valu=typeelem)
    
    if gdat.typeexpr == 'gmix':
        legdelem = ['Cluster']
    if gdat.typeexpr.startswith('HST_WFC3'):
        legdelem = ['Subhalo']
    setp_varb(gdat, 'legdelem', valu=legdelem)

    ### PSF model
    #### angular profile
    if gdat.typeexpr == 'ferm':
        typemodlpsfn = 'doubking'
    if gdat.typeexpr == 'chan':
        typemodlpsfn = 'singking'
    if gdat.typeexpr == 'sdss':
        typemodlpsfn = 'singgaus'
    if gdat.typeexpr.startswith('HST_WFC3'):
        typemodlpsfn = 'singgaus'
    if gdat.typeexpr == 'tess':
        typemodlpsfn = 'singgaus'
    if gdat.typeexpr == 'gmix':
        typemodlpsfn = 'singgaus'
    if gdat.typeexpr == 'fire':
        typemodlpsfn = None
    setp_varb(gdat, 'typemodlpsfn', valu=typemodlpsfn)
    
    #### background names
    listnameback = ['isot']
    if gdat.typeexpr == 'ferm':
        listnameback.append('fdfm')
    #if gdat.typeexpr == 'chan':
    #    listnameback.append('part')
    setp_varb(gdat, 'listnameback', valu=listnameback)
        
    defn_tdim(gdat)
    
    for strgmodl in gdat.liststrgmodl:
        # set up the indices of the model
        setp_indxpara_init(gdat, strgmodl=strgmodl)
   
    gdat.numbdeflsubhplot = 2
    #gdat.refr.colr = 'mediumseagreen'
    #gdat.fitt.colr = 'deepskyblue'

    if gdat.factpriodoff is None:
        gdat.factpriodoff = 1.
        
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
        gdat.numbdqltfull = gdat.sbrtdata.shape[2]
        
        if gdat.typepixl == 'heal':
            # temp
            gdat.numbsidecart = 100
            gdat.numbsidecarthalf = int(gdat.numbsidecart / 2)
            gdat.numbsideheal = int(np.sqrt(gdat.numbpixlfull / 12))
    
    if gdat.typeexpr.startswith('HST_WFC3'):
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
        if gdat.typeexpr == 'sdss' or gdat.typeexpr == 'chan' or gdat.typeexpr.startswith('HST'):
            gdat.anglfact = 3600 * 180. / np.pi
        if gdat.typeexpr == 'sche' or gdat.typeexpr == 'gmix':
            gdat.anglfact = 1.
    
    if gdat.numbsidecart is not None and gdat.typepixl == 'cart' and not gdat.boolforccart and isinstance(strgexpo, str):
        raise Exception('numbsidecart argument should not be provided when strgexpo is a file name and pixelization is Cartesian.')
                
    if gdat.typepixl == 'heal' or gdat.typepixl == 'cart' and gdat.boolforccart:
        if gdat.numbsidecart is None:
            gdat.numbsidecart = 100
    
    # exposure
    gdat.boolcorrexpo = gdat.expo is not None
    if gdat.typeexpo == 'cons':
        if gdat.typedata == 'simu':
            if gdat.numbsidecart is None:
                gdat.numbsidecart = 100
        if gdat.typedata == 'simu':
            if gdat.typepixl == 'heal':
                gdat.expo = np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
            if gdat.typepixl == 'cart':
                gdat.expo = np.ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbdqltfull))
        if gdat.typedata == 'inpt':
            gdat.expo = np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
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
    
    if gdat.typedata == 'simu':
        if gdat.typepixl == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2
        if gdat.typepixl == 'heal':
            gdat.numbpixlfull = 12 * gdat.numbsideheal**2
        
   
    # Boolean flag to indicate binning in space
    gdat.boolbindspat = gdat.numbpixlfull != 1

    print('gdat.boolbindspat')
    print(gdat.boolbindspat)
    
    if gdat.boolcorrexpo and np.amin(gdat.expo) == np.amax(gdat.expo) and not isinstance(gdat.strgexpo, float):
        raise Exception('Bad input exposure map.')
    
    if gdat.boolbindspat:
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
            if gdat.typeexpr.startswith('HST_WFC3'):
                gdat.maxmgangdata = 2. / gdat.anglfact
            if gdat.typeexpr == 'gmix':
                gdat.maxmgangdata = 1. / gdat.anglfact
        
        print('gdat.numbsidecart')
        print(gdat.numbsidecart)
        print('gdat.maxmgangdata')
        print(gdat.maxmgangdata)
    
        # pixelization
        if gdat.typepixl == 'cart':
            if gdat.typeexpr == 'HST_WFC3_UVIS':
                gdat.sizepixl = 0.04 # [arcsec]
            if gdat.typeexpr == 'HST_WFC3_IR':
                gdat.sizepixl = 0.13 # [arcsec]
            gdat.apix = gdat.sizepixl**2
        if gdat.typepixl == 'heal':
            temp, temp, temp, gdat.apix = tdpy.retr_healgrid(gdat.numbsideheal)
            gdat.sizepixl = np.sqrt(gdat.apix)
    
    # factor by which to multiply the y axis limits of the surface brightness plot
    if gdat.numbpixlfull == 1:
        gdat.factylimtbrt = [1e-4, 1e7]
    else:
        gdat.factylimtbrt = [1e-4, 1e3]

    # grid
    gdat.minmxposdata = -gdat.maxmgangdata
    gdat.maxmxposdata = gdat.maxmgangdata
    gdat.minmyposdata = -gdat.maxmgangdata
    gdat.maxmyposdata = gdat.maxmgangdata
    
    if gdat.typepixl == 'cart' and gdat.boolforccart:
        if gdat.typedata == 'inpt':
            sbrtdatatemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
            for i in gdat.indxenerfull:
                for m in gdat.indxdqltfull:
                    sbrtdatatemp[i, :, m] = tdpy.retr_cart(gdat.sbrtdata[i, :, m], \
                                                    numbsidexpos=gdat.numbsidecart, numbsideypos=gdat.numbsidecart, \
                                                    minmxpos=gdat.anglfact*gdat.minmxposdata, maxmxpos=gdat.anglfact*gdat.maxmxposdata, \
                                                    minmypos=gdat.anglfact*gdat.minmyposdata, maxmypos=gdat.anglfact*gdat.maxmyposdata).flatten()
            gdat.sbrtdata = sbrtdatatemp

        if gdat.boolcorrexpo:
            expotemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbdqltfull))
            for i in gdat.indxenerfull:
                for m in gdat.indxdqltfull:
                    expotemp[i, :, m] = tdpy.retr_cart(gdat.expo[i, :, m], \
                                                    numbsidexpos=gdat.numbsidecart, numbsideypos=gdat.numbsidecart, \
                                                    minmxpos=gdat.anglfact*gdat.minmxposdata, maxmxpos=gdat.anglfact*gdat.maxmxposdata, \
                                                    minmypos=gdat.anglfact*gdat.minmyposdata, maxmypos=gdat.anglfact*gdat.maxmyposdata).flatten()
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
    if gdat.typeexpr.startswith('HST_WFC3'):
        gdat.listspecconvunit = [['en03', 'ergs']]
    if gdat.typeexpr == 'fire':
        gdat.listspecconvunit = [['en00', 'imum']]
    
    # temp
    #if gdat.typeexpr == 'chan' and (gdat.anlytype.startswith('home') or gdat.anlytype.startswith('extr')):
    #    gmod.lablpopl = ['AGN', 'Galaxy']

    if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr == 'fire':
        gdat.enerdiff = True
    if gdat.typeexpr.startswith('HST_WFC3') or gdat.typeexpr == 'gmix' or gdat.typeexpr == 'tess':
        gdat.enerdiff = False
    
    if gdat.indxenerincl is None:
        
        # default
        if gdat.boolbinsener:
            gdat.indxenerincl = np.arange(gdat.blimpara.enerfull.size - 1)
        
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
        if gdat.typeexpr.startswith('HST_WFC3'):
            gdat.indxenerincl = np.array([0])
            #gdat.indxenerincl = np.array([1])
            #gdat.indxenerincl = np.array([0, 1])
        if gdat.typeexpr == 'gmix':
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
        gdat.blimpara.ener = gdat.blimpara.enerfull[gdat.indxenerinclbins]
        gdat.bctrpara.ener = np.sqrt(gdat.blimpara.ener[1:] * gdat.blimpara.ener[:-1])
        gdat.deltener = gdat.blimpara.ener[1:] - gdat.blimpara.ener[:-1]
        gdat.minmener = gdat.blimpara.ener[0]
        gdat.maxmener = gdat.blimpara.ener[-1]
        setp_varb(gdat, 'ener')

        gdat.limtener = [np.amin(gdat.blimpara.ener), np.amax(gdat.blimpara.ener)] 
    if gdat.boolbinsener: 
        if gdat.numbener > 1:
            gdat.enerpivt = gdat.bctrpara.ener[gdat.indxenerpivt]
        # energy bin indices other than that of the pivot bin
        gdat.indxenerinde = np.setdiff1d(gdat.indxener, gdat.indxenerpivt)
    
        # temp
        if gdat.typeexpr == 'chan':
            gdat.edis = 0.3 * np.sqrt(gdat.blimpara.ener) / 2.35
            gdat.edisintp = sp.interpolate.interp1d(gdat.blimpara.ener, gdat.edis, fill_value='extrapolate')
        else:
            gdat.edis = None
            gdat.edisintp = None

    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
    
        setp_varb(gdat, 'cntpmodl', labl=['$C_{M}$', ''], scal='asnh', strgmodl=strgmodl)
        
        if isinstance(gdat.fitt.maxmpara, dict):
            raise Exception('')

        if isinstance(gdat.true.maxmpara, dict):
            raise Exception('')

        for strgstat in ['this', 'next']:
            setattr(gmod, strgstat, tdpy.gdatstrt())
        
        # number of elements
        if strgmodl == 'true':
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lens':
                    numbelem = 2
                else:
                    numbelem = 2
                setp_varb(gdat, 'numbelem', minm=0, maxm=3, labl=['N', ''], scal='drct', valu=numbelem, popl=l, strgmodl=strgmodl, strgstat='this')
        if strgmodl == 'fitt':
            setp_varb(gdat, 'numbelem', minm=0, maxm=3, labl=['N', ''], scal='drct', popl='full', strgmodl=strgmodl)
        
        # total number of elements summed over populations
        setp_varb(gdat, 'numbelemtotl', maxm=10, labl=['$N_{tot}$', ''], scal='drct', strgmodl=strgmodl)

        ## hyperparameters
        setp_varb(gdat, 'typemodltran', valu='drct', strgmodl=strgmodl)
        
        if gmod.typemodltran == 'pois':
            setp_varb(gdat, 'meanelem', minm=0.1, maxm=1000., scal='logt', popl='full', strgmodl=strgmodl)
    
        #### boolean flag background
        if gdat.typeexpr == 'chan':
            if gdat.numbpixlfull == 1:
                boolspecback = [True, True]
            else:
                boolspecback = [False, False]
        else:
            boolspecback = [False for k in gmod.indxback]
        setp_varb(gdat, 'boolspecback', valu=boolspecback, strgmodl=strgmodl)
    
        # type of the spatial extent of the model evaluation due to elements
        typeelemspateval = [[] for l in gmod.indxpopl]
        for l in gmod.indxpopl:
            # these  element types slow down execution!
            if gmod.typeelem[l] == 'lens' or gmod.typeelem[l].startswith('lghtline') or gmod.typeelem[l] == 'clusvari' or gmod.typeelem[l] == 'lghtgausbgrd':
                typeelemspateval[l] = 'full'
            else:
                typeelemspateval[l] = 'locl'
        setp_varb(gdat, 'typeelemspateval', valu=typeelemspateval, strgmodl=strgmodl)
        
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
        if gdat.typeexpr.startswith('HST_WFC3'):
            gdat.strgexprname = 'HST'
        if gdat.typeexpr == 'sche':
            gdat.strgexprname = 'XXXXX'
        if gdat.typeexpr == 'gmix':
            gdat.strgexprname = 'TGAS-RAVE'
    
    if gdat.lablgangunit is None:
        if gdat.typeexpr == 'ferm':
            gdat.lablgangunit = '$^o$'
        if gdat.typeexpr == 'gmix':
            gdat.lablgangunit = ''
        if gdat.typeexpr == 'sdss' or gdat.typeexpr == 'chan' or gdat.typeexpr.startswith('HST'):
            gdat.lablgangunit = '$^{\prime\prime}$'
    
    if gdat.lablxpos is None:
        if gdat.typeexpr == 'gmix':
            gdat.lablxpos = r'L_{z}'
        else:
            if gdat.typeexpr == 'ferm' and gdat.xposcntr == 0 and gdat.yposcntr == 0:
                gdat.lablxpos = r'l'
            else:
                gdat.lablxpos = r'\theta_1'
    if gdat.lablypos is None:
        if gdat.typeexpr == 'gmix':
            gdat.lablypos = r'E_k'
        else:
            if gdat.typeexpr == 'ferm' and gdat.xposcntr == 0 and gdat.yposcntr == 0:
                gdat.lablypos = r'b'
            else:
                gdat.lablypos = r'\theta_2'

    if gdat.strgenerunit is None:
        if gdat.typeexpr == 'ferm':
            gdat.strgenerunit = 'GeV'
            gdat.nameenerunit = 'gevv'
        if gdat.typeexpr == 'chan':
            gdat.strgenerunit = 'keV'
            gdat.nameenerunit = 'kevv'
        if gdat.typeexpr == 'gmix':
            gdat.strgenerunit = ''
            gdat.nameenerunit = ''
        if gdat.typeexpr.startswith('HST_WFC3'):
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
        if gdat.typeexpr.startswith('HST_WFC3'):
            gdat.nameexpr = 'HST'
        if gdat.typeexpr == 'gaia':
            gdat.nameexpr = 'Gaia'
    
    ## Lensing
    if gdat.radispmr is None:
        if gdat.typeexpr == 'ferm':
            gdat.radispmr = 0.6 / gdat.anglfact
        if gdat.typeexpr.startswith('HST_WFC3'):
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
        if gdat.typeexpr == 'gmix':
            gdat.radispmr = 0.2
    
    print('gdat.radispmr')
    print(gdat.radispmr)

    if gdat.anglassc is None:
        gdat.anglassc = 5. * gdat.radispmr
    
    print('gdat.anglassc')
    print(gdat.anglassc)

    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
    
        if gdat.boolbindspat:
            if gdat.typeexpr == 'chan' or gdat.typeexpr == 'sdss':
                numbpsfpform = 0
                gmod.numbpsfptotl = 0
            if gdat.typeexpr == 'chan':
                retr_psfpchan(gmod)
            if gdat.typeexpr == 'ferm':
                retr_psfpferm(gmod)
            if gdat.typeexpr == 'sdss':
                retr_psfpsdss(gmod)
            if gdat.typeexpr.startswith('HST_WFC3'):
                retr_psfphubb(gdat, gmod)
            if gdat.typeexpr == 'tess':
                retr_psfptess(gmod)
            if gdat.typeexpr == 'gmix':
                retr_psfpsdyn(gmod)
    

    # model evaluation approximation error tolerance in units of the fraction of the lowest PS flux
    if gdat.specfraceval is None:
        if gdat.typeexpr == 'ferm':
            gdat.specfraceval = 0.5
        else:
            gdat.specfraceval = 0.1

    if gdat.boolbindspat:
        gdat.blimpara.xposcart = np.linspace(gdat.minmxposdata, gdat.maxmxposdata, gdat.numbsidecart + 1)
        gdat.blimpara.yposcart = np.linspace(gdat.minmyposdata, gdat.maxmyposdata, gdat.numbsidecart + 1)
        gdat.bctrpara.xposcart = (gdat.blimpara.xposcart[0:-1] + gdat.blimpara.xposcart[1:]) / 2.
        gdat.bctrpara.yposcart = (gdat.blimpara.yposcart[0:-1] + gdat.blimpara.yposcart[1:]) / 2.
    
    # reference elements
    gdat.numbrefr = 0
    if gdat.typedata == 'simu':
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
    
    gdat.numbpixlxposshft = []
    gdat.numbpixlyposshft = []
    gdat.refrindxpoplassc = [[] for q in gdat.indxrefr] 
    
    # temp -- this allows up to 3 reference populations
    gdat.true.colrelem = ['darkgreen', 'olivedrab', 'mediumspringgreen']
    # temp -- this allows up to 3 reference populations
    gdat.fitt.colrelem = ['royalblue', 'dodgerblue', 'navy']
    if gdat.typedata == 'simu':
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
            if 'xpos' in gdat.refr.namepara.elem[q] and 'ypos' in gdat.refr.namepara.elem[q]:
                gdat.refr.namepara.elem[q] += ['gang', 'aang']
            for strgfeat in gdat.refr.namepara.elem[q]:
                setattr(gdat.refr, strgfeat, [[] for q in gdat.indxrefr])
        
        if gdat.typeexpr == 'ferm':
            retr_refrfermfinl(gdat)
        if gdat.typeexpr == 'chan':
            retr_refrchanfinl(gdat)
    
    if gdat.typeexpr.startswith('HST_WFC3'):
        boollenshost = True
    else:
        boollenshost = False
    setp_varb(gdat, 'boollenshost', valu=boollenshost)
  
    ## Boolean flag to turn on deflection due to elements
    if gdat.typeexpr.startswith('HST_WFC3'):
        boollenssubh = True
    else:
        boollenssubh = False
    setp_varb(gdat, 'boollenssubh', valu=boollenssubh)
  
    if gdat.typeexpr.startswith('HST_WFC3'):
        boollens = True
    else:
        boollens = False
    setp_varb(gdat, 'boollens', valu=boollens)
  
    if gdat.typeexpr.startswith('HST_WFC3'):
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
    
    #setp_varb(gdat, 'xpos', minm=-10., maxm=10., labl=['$l$', ''])
    
    for strgmodl in gdat.liststrgmodl:
        
        gmod = getattr(gdat, strgmodl)
        
        ## group the maximum number of elements for each population into an array 'temp' this repeats the generic process and need to be removed
        gmod.minmpara.numbelem = np.empty(gmod.numbpopl, dtype=int)
        gmod.maxmpara.numbelem = np.zeros(gmod.numbpopl, dtype=int)
        for l in gmod.indxpopl:
            gmod.minmpara.numbelem[l] = getattr(gmod.minmpara, 'numbelempop%d' % l)
            gmod.maxmpara.numbelem[l] = getattr(gmod.maxmpara, 'numbelempop%d' % l)
    
    if gdat.typedata == 'simu':
        setp_indxpara_finl(gdat, strgmodl='true')
            
    for strgmodl in gdat.liststrgmodl:
        
        gmod = getattr(gdat, strgmodl)
        
        if gdat.typeexpr.startswith('HST_WFC3'):
            minm = 0.1
            maxm = 10.
            
            setp_varb(gdat, 'defs', minm=minm, maxm=maxm, scal='powr', labl=[r'$\alpha$', ''], strgmodl=strgmodl)
            setp_varb(gdat, 'mcut', minm=minm, maxm=maxm, scal='powr', labl=['$m_c$', ''], strgmodl=strgmodl)
            setp_varb(gdat, 'asca', minm=minm, maxm=maxm, scal='self', labl=['$\theta_s$', ''], strgmodl=strgmodl)
            setp_varb(gdat, 'acut', minm=minm, maxm=maxm, scal='self', labl=['$\theta_c$', ''], strgmodl=strgmodl)
            setp_varb(gdat, 'rele', minm=minm, maxm=maxm, scal='self', labl=['$R$', ''], strgmodl=strgmodl)
            ## distance to the source
            setp_varb(gdat, 'distsour', minm=minm, maxm=maxm, scal='powr', labl=['$\delta \theta_S$', ''], strgmodl=strgmodl)
            ## relevance
            #setp_varb(gdat, 'rele', minm=minm, maxm=maxm, scal='powr', labl=['$R_{%d}$', ''], strgmodl=strgmodl)
            
            # fields defined on a specific grid
            ## total deflection
            setp_varb(gdat, 'defl', minm=gdat.maxmgangdata/1e4, maxm=gdat.maxmgangdata, numbbins=10, scal='powr', labl=[r'$\alpha$', ''], strgmodl=strgmodl)
            ## subhalo deflection
            setp_varb(gdat, 'deflsubh', minm=gdat.maxmgangdata/1e4, maxm=gdat.maxmgangdata, numbbins=10, scal='powr', labl=['$\alpha_s$', ''], strgmodl=strgmodl)
            ## deflection profile of an individual subhalo
            setp_varb(gdat, 'deflprof', minm=minm, maxm=maxm, scal='powr', labl=['$\alpha(r)$', ''], strgmodl=strgmodl)
            
            for l in gmod.indxpopl:
                setp_varb(gdat, 'defs', minm=minm, maxm=maxm, scal='powr', labl=[r'$\alpha$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'defs', minm=minm, maxm=maxm, scal='powr', labl=[r'$\alpha$', ''], popl=l, iele='full', strgmodl=strgmodl)
        
                setp_varb(gdat, 'asca', minm=minm, maxm=maxm, scal='self', labl=['$\theta_s$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'asca', minm=minm, maxm=maxm, scal='self', labl=['$\theta_s$', ''], popl=l, iele='full', strgmodl=strgmodl)
        
                setp_varb(gdat, 'acut', minm=minm, maxm=maxm, scal='self', labl=['$\theta_c$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'acut', minm=minm, maxm=maxm, scal='self', labl=['$\theta_c$', ''], popl=l, iele='full', strgmodl=strgmodl)
        
                setp_varb(gdat, 'mcut', minm=minm, maxm=maxm, scal='powr', labl=['$m_c$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'mcut', minm=minm, maxm=maxm, scal='powr', labl=['$m_c$', ''], popl=l, strgmodl=strgmodl, iele='full')
        
                setp_varb(gdat, 'rele', minm=minm, maxm=maxm, scal='powr', labl=['$R_{%d}$' % l, ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'rele', minm=minm, maxm=maxm, scal='powr', labl=['$R_{%d}$' % l, ''], popl=l, strgmodl=strgmodl, iele='full')
        
                setp_varb(gdat, 'deflprof', minm=minm, maxm=maxm, scal='powr', labl=['$\alpha_{%d}(r)$' % l, ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'deflprof', minm=minm, maxm=maxm, scal='powr', labl=['$\alpha_{%d}(r)$' % l, ''], popl=l, strgmodl=strgmodl, iele='full')
                
                setp_varb(gdat, 'distsour', minm=minm, maxm=maxm, scal='powr', labl=['$\delta \theta_S$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'distsour', minm=minm, maxm=maxm, scal='powr', labl=['$\delta \theta_S$', ''], popl=l, strgmodl=strgmodl, iele='full')
        
        ### background parameters
        if gdat.typeexpr == 'chan':
            if gdat.anlytype.startswith('extr'):
                meanbacpbac1 = 1.
            else:
                meanbacpbac1 = 70.04
            stdvbacpbac1 = 1e-5 * meanbacpbac1
            setp_varb(gdat, 'bacp', mean=meanbacpbac1, stdv=stdvbacpbac1, back=1, scal='gaus', strgmodl='true')

            if gdat.numbpixlfull == 1:
                setp_varb(gdat, 'bacp', minm=1., maxm=100., back=0)
            else:
                setp_varb(gdat, 'bacp', minm=0.1, maxm=1000., ener='full', back=0)
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

        if gdat.typeexpr.startswith('HST_WFC3'):
            setp_varb(gdat, 'bacp', minm=1e-1, maxm=1e3, valu=1e1, labl=['$A$', ''], scal='logt', ener=0, back=0, strgmodl=strgmodl, strgstat='this')
        if gdat.typeexpr == 'gmix':
            setp_varb(gdat, 'bacp', minm=1e-1, maxm=1e3, valu=1e1, labl=['$A$', ''], scal='logt', ener=0, back=0, strgmodl=strgmodl)
        
        if gdat.typeexpr == 'fire':
            bacp = [1e-1, 1e1]
        if gdat.typeexpr == 'tess':
            bacp = [1e-1, 1e1]
            setp_varb(gdat, 'bacp', limt=bacp, ener='full', back=0)
        
        # temp
        if gmod.boollens:
            gmod.indxisfr = np.arange(1)

        if gdat.typeexpr.startswith('HST_WFC3'):
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
        if gdat.boolbindspat:
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
        
        if gdat.boolbindspat:
            minmdefs = 0.003 / gdat.anglfact
            setp_varb(gdat, 'minmdefs', valu=minmdefs, strgmodl=strgmodl)
    
        if gdat.typeexpr == 'ferm':
            setp_varb(gdat, 'curv', limt=[-1., 1.], strgmodl=strgmodl)
    
        if gdat.boolbindspat:
            maxmdefs = 1. / gdat.anglfact
            setp_varb(gdat, 'maxmdefs', valu=maxmdefs, strgmodl=strgmodl)
    
        # true model parameters
        if gdat.typedata == 'simu':
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
            minm = 0.
            maxm = 0.1 / gdat.anglfact
            setp_varb(gdat, 'asca', minm=minm, maxm=maxm, strgmodl=strgmodl, popl=l)
            
            ### projected cutoff radius
            limtacut = np.array([0., 2.]) / gdat.anglfact
            minm = 0.
            maxm = 2. / gdat.anglfact
            setp_varb(gdat, 'acut', minm=minm, maxm=maxm, strgmodl=strgmodl, popl=l)

        if gdat.boolbindspat:

            setp_varb(gdat, 'gangdisttype', valu=['self'], strgmodl=strgmodl)
    
            for l in gmod.indxpopl:
                if gmod.typespatdist[l] == 'gangexpo':
                    setp_varb(gdat, 'maxmgang', valu=gmod.maxmxpos, strgmodl=strgmodl)
                    if gdat.typeexpr == 'ferm':
                        gangdistsexp = 5. / gdat.anglfact
                    setp_varb(gdat, 'gangdistsexp', valu=gangdistsexp, strgmodl=strgmodl, popl=l)
                if gmod.typespatdist[l] == 'dsrcexpo':
                    if gdat.typeexpr.startswith('HST_WFC3'):
                        dsrcdistsexp = 0.5 / gdat.anglfact
                    setp_varb(gdat, 'dsrcdistsexp', valu=dsrcdistsexp, strgmodl=strgmodl, popl=l)
    
            if gmod.boollens:
                setp_varb(gdat, 'xpossour', valu=0., mean=0., stdv=gdat.stdvhostsour, labl=['$x_{S}$', 'arcsec'], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'ypossour', valu=0., mean=0., stdv=gdat.stdvhostsour, labl=['$y_{S}$', 'arcsec'], strgmodl=strgmodl, strgstat='this')
                
                setp_varb(gdat, 'sherextr', limt=[0., 0.1], labl=['$\rho_{ext}$', ''], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'anglsour', limt=[0., np.pi], labl=['$\phi_{S}$', 'degree'], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'sangextr', limt=[0., np.pi], labl=['$\phi_{ext}$', ''], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'sangextr', labl=['$\phi_{ext}$', ''], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'sizesour', limt=[0.1 / gdat.anglfact, 2. / gdat.anglfact], \
                                                                        labl=['$R_{S}$', 'arcsec'], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'ellpsour', limt=[0., 0.5], labl=['$\epsilon_{S}$', ''], strgmodl=strgmodl, strgstat='this')
            
                setp_varb(gdat, 'fluxsour', valu=1e-20, limt=np.array([1e-22, 1e-17]), labl=['$f_{S}$', 'erg/s'], strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'sindsour', limt=np.array([0., 4.]), strgmodl=strgmodl, strgstat='this')
        
            if gmod.boollenshost or gmod.boolemishost:
                for e in gmod.indxsersfgrd:
                    setp_varb(gdat, 'xposhost', labl=['$x_{H%d}$' % e, 'arcsec'], isfr=e)
                    setp_varb(gdat, 'yposhost', labl=['$y_{H%d}$' % e, 'arcsec'], isfr=e)
                    setp_varb(gdat, 'fluxhost', labl=['$f_{H%d}$' % e, 'erg/s'], isfr=e)
                    setp_varb(gdat, 'sizehost', labl=['$R_{H%d}$' % e, 'arcsec'], isfr=e)
                    setp_varb(gdat, 'beinhost', labl=['$\theta_{E,H%d}$' % e, 'arcsec'], isfr=e)
                    setp_varb(gdat, 'serihost', labl=['$n_{Ser,H%d}$' % e, ''], isfr=e)
            
                    setp_varb(gdat, 'sindhost', limt=np.array([0., 4.]), isfr=e, strgmodl=strgmodl)
                    setp_varb(gdat, 'beinhost', limt=[0.5 / gdat.anglfact, 2. / gdat.anglfact], isfr=e, strgmodl=strgmodl)
                    setp_varb(gdat, 'ellphost', limt=[0., 0.5], labl=['$\epsilon_{H%d}$' % e, ''], isfr=e, strgmodl=strgmodl)
                    setp_varb(gdat, 'anglhost', limt=[0., np.pi], labl=['$\phi_{H%d}$' % e, 'degree'], isfr=e, strgmodl=strgmodl, strgstat='this')
                    setp_varb(gdat, 'xposhost', valu=0., mean=0., stdv=gdat.stdvhostsour, strgmodl='true', isfr=e, strgstat='this')
                    setp_varb(gdat, 'yposhost', valu=0., mean=0., stdv=gdat.stdvhostsour, strgmodl='true', isfr=e, strgstat='this')
                
            if gmod.boolemishost:
                for e in gmod.indxsersfgrd:
                    setp_varb(gdat, 'fluxhost', valu=1e-16, limt=np.array([1e-20, 1e-15]), strgmodl=strgmodl, strgstat='this', isfr=e)
                    setp_varb(gdat, 'sizehost', valu=1. / gdat.anglfact, limt=[0.1 / gdat.anglfact, 4. / gdat.anglfact], strgmodl=strgmodl, strgstat='this', isfr=e)
                    setp_varb(gdat, 'serihost', valu=4., limt=[1., 8.], strgmodl=strgmodl, strgstat='this', isfr=e)
            
            if gmod.boollenshost:
                setp_varb(gdat, 'beinhost', valu=1.5 / gdat.anglfact, isfr='full', strgmodl=strgmodl, strgstat='this')
                setp_varb(gdat, 'fluxsour', valu=1e-18, strgmodl='true', strgstat='this')
                setp_varb(gdat, 'sindsour', valu=1.5, strgmodl='true', strgstat='this')
                setp_varb(gdat, 'sindhost', valu=2.5, strgmodl='true', isfr='full', strgstat='this')
            
            if strgmodl == 'fitt':
                if gmod.boollenshost or gmod.boolemishost:
                    setp_varb(gdat, 'xposhost', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', isfr='full')
                    setp_varb(gdat, 'yposhost', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', isfr='full')
                if gmod.boollens:
                    setp_varb(gdat, 'xpossour', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt')
                    setp_varb(gdat, 'ypossour', limt=[-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt')
        
            # temp -- to be removed
            #gmod.factxpos = gmod.maxmxpos - gmod.minmxpos
            #gmod.factypos = gmod.maxmypos - gmod.minmypos
            #gmod.minmaang = -np.pi
            #gmod.maxmaang = pi
        
        # hyperparameters
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lens':
                scal = 'gaus'
                minm = None
                maxm = None
                mean = 1.9
                stdv = 0.5
            elif gmod.typeelem[l].startswith('clus'):
                scal = 'logt'
                minm = 0.5
                maxm = 2.
                valu = 1.5
                mean = None
                stdv = None
            else:
                scal = 'logt'
                minm = 0.5
                maxm = 3.
                mean = None
                stdv = None
            name = 'slopprio' + gmod.nameparagenrelemampl[l]
            
            setp_varb(gdat, name, minm=minm, maxm=maxm, scal=scal, mean=mean, stdv=stdv, labl=[r'$\alpha$', ''], popl=l, strgmodl=strgmodl)
            
            # below this line is unnecessary to be deleted
            if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
                setp_varb(gdat, 'gwdtslop', limt=[0.5, 4.], scal='logt', popl=l, strgmodl=strgmodl)
    
            if gdat.boolbindspat:
                setp_varb(gdat, 'spatdistcons', valu=1e-3, popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'gangslop', valu=1.1, popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'yposdistscal', valu=2. / gdat.anglfact, popl=l, strgmodl=strgmodl)

            if gdat.typeexpr == 'ferm':
                setp_varb(gdat, 'sloplowrprioflux', valu=1.5, popl=l)
                setp_varb(gdat, 'slopupprprioflux', valu=2.5, popl=l)
                setp_varb(gdat, 'brekprioflux', valu=1e-9, popl=l)
            if gmod.typeelem[l] == 'lghtpnts':
                setp_varb(gdat, 'slopprioflux', valu=2.2, popl=l, strgmodl=strgmodl)
            if gmod.typeelem[l].startswith('lghtline'):
                setp_varb(gdat, 'slopprioflux', valu=2., popl=l, strgmodl=strgmodl)
            if gmod.typeelem[l] == 'lens':
                setp_varb(gdat, 'sloppriodefs', valu=1.9, popl=l, strgmodl=strgmodl, strgstat='this')
            
            if gmod.typeelem[l] == 'lens':
                setp_varb(gdat, 'meanprioasca', valu=0.05 / gdat.anglfact, popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'stdvprioasca', valu=0.04 / gdat.anglfact, popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'meanprioacut', valu=1. / gdat.anglfact, popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'stdvprioacut', valu=0.04 / gdat.anglfact, popl=l, strgmodl=strgmodl)
            
            if gmod.typeelem[l] == 'lghtgausbgrd' or gmod.typeelem[l] == 'clusvari':
                setp_varb(gdat, 'gwdtslop', valu=2., popl=l, strgmodl=strgmodl)
            
            if gdat.typeexpr == 'ferm':
                sinddistmean = 2.15
            if gdat.typeexpr == 'chan':
                sinddistmean = 1.
            if gdat.typeexpr.startswith('HST_WFC3'):
                sinddistmean = 1.
            if gdat.typeexpr == 'ferm' or gdat.typeexpr == 'chan' or gdat.typeexpr.startswith('HST'):
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
            
    # copy the true model to the fitting model if the fitting model parameter has not been specified
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
    if gdat.typepixl == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')
    
    # list of names of variables to be calculated
    retr_strgcalc(gdat)

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
    
    if gdat.typedata != 'simu':
        gdat.refr.numbelem = np.zeros(gdat.numbrefr, dtype=int)
    
    setup_pcat_model(gdat)

    # set the reference model to true model
    gdat.refr.labl = 'True'
    print('Setting the remaining few parameters in the reference model to those in the true model...')
    for strgmodl in gdat.liststrgmodl:
        if strgmodl == 'true':
            
            #for strg, valu in gdat.true.__dict__.items():
            #    setattr(gdat.refr, strg, valu)
            
            for name in ['listmrkrmiss', 'listlablmiss', 'colr', 'colrelem', 'namepara', 'nameparagenrelemampl', 'numbelem']:
                setattr(gdat.refr, name, getattr(gdat.true, name))
    
    # set up the indices of the fitting model
    setp_indxpara_finl(gdat, strgmodl='fitt')
        
    for strgmodl in gdat.liststrgmodl:
        
        gmod = getattr(gdat, strgmodl)
        
        minmredshost = 0.01
        maxmredshost = 0.4
        minmredssour = 0.01
        maxmredssour = 2.
        numbreds = 200
        
        asca = 0.1 / gdat.anglfact
        acut = 1. / gdat.anglfact
    
        # to be deleted
        #minmmass = np.zeros((numbreds + 1, numbreds + 1))
        #maxmmass = np.zeros((numbreds + 1, numbreds + 1))
        #for k, redshost in enumerate(gmod.blimpara.redshost):
        #    for n, redssour in enumerate(gmod.blimpara.redssour):
        #        if redssour > redshost:
        #            gmod.adislens = gdat.adisobjt(redshost)
        #            gmod.adissour = gdat.adisobjt(redssour)
        #            gmod.adislenssour = gmod.adissour - (1. + redshost) / (1. + redssour) * gmod.adislens
        #            factmcutfromdefs = chalcedon.retr_factmcutfromdefs(gmod.adissour, gmod.adislens, gmod.adislenssour, asca, acut)
        #            minmmass[n, k] = np.log10(factmcutfromdefs * gmod.minmpara.defs)
        #            maxmmass[n, k] = np.log10(factmcutfromdefs * gmod.maxmpara.defs)
        
        if gdat.boolbindspat:
            minm = -gdat.maxmgangdata
            maxm = gdat.maxmgangdata
            for l in gmod.indxpopl:
                if gdat.typeexpr == 'ferm':
                    lablxpos = '$l$'
                    lablypos = '$b$'
                    lablunitxpos = 'degree'
                    lablunitypos = 'degree'
                if gdat.typeexpr.startswith('HST_WFC3'):
                    lablxpos = '$x$'
                    lablypos = '$y$'
                    lablunitxpos = 'arcsec'
                    lablunitypos = 'arcsec'
                setp_varb(gdat, 'xpos', minm=minm, maxm=maxm, labl=[lablxpos, lablunitxpos], strgmodl=strgmodl)
                setp_varb(gdat, 'ypos', minm=minm, maxm=maxm, labl=[lablypos, lablunitypos], strgmodl=strgmodl)
                setp_varb(gdat, 'xpos', minm=minm, maxm=maxm, labl=[lablxpos, lablunitxpos], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'ypos', minm=minm, maxm=maxm, labl=[lablypos, lablunitypos], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'xpos', minm=minm, maxm=maxm, labl=[lablxpos, lablunitxpos], popl=l, iele='full', strgmodl=strgmodl)
                setp_varb(gdat, 'ypos', minm=minm, maxm=maxm, labl=[lablypos, lablunitypos], popl=l, iele='full', strgmodl=strgmodl)
        
        if gdat.typeexpr == 'gmix':
            minm = 0.1
            maxm = 10.
            for l in gmod.indxpopl:
                setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', labl=['N', ''], strgmodl=strgmodl)
                setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', labl=['N', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'nobj', minm=minm, maxm=maxm, scal='powr', labl=['N', ''], popl=l, iele='full', strgmodl=strgmodl)
        
        






        if gdat.boolbindspat:
            for l in gmod.indxpopl:
                setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, labl=[r'$\theta$', ''], strgmodl=strgmodl)
                setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, labl=[r'$\theta$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'aang', minm=-np.pi, maxm=np.pi, labl=[r'$\theta$', ''], popl=l, strgmodl=strgmodl, iele='full')
                setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, labl=[r'$\psi$', ''], strgmodl=strgmodl)
                setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, labl=[r'$\psi$', ''], popl=l, strgmodl=strgmodl)
                setp_varb(gdat, 'gang', minm=0, maxm=gdat.maxmgangdata, labl=[r'$\psi$', ''], popl=l, strgmodl=strgmodl, iele='full')

        # loglikelihood difference for each element
        setp_varb(gdat, 'deltllik', labl=['$\Delta \log L$', ''], minm=1., maxm=100., strgmodl=strgmodl)
        setp_varb(gdat, 'deltllik', labl=['$\Delta \log L$', ''], minm=1., maxm=100., popl=l, strgmodl=strgmodl)
        setp_varb(gdat, 'deltllik', labl=['$\Delta \log L$', ''], minm=1., maxm=100., popl=l, strgmodl=strgmodl, iele='full')
        
    # construct the fitting model
    setp_paragenrscalbase(gdat, strgmodl='fitt')
    
    if gdat.typedata == 'simu':
        # construct the true model
        setp_paragenrscalbase(gdat, strgmodl='true')
   
    gdat.refr.indxpoplfittassc = gdat.fitt.indxpopl
    gdat.fitt.indxpoplrefrassc = gdat.fitt.indxpopl

    print('Defining minima and maxima for derived parameters...')
    
    # derived lens parameter minima and maxima
    for strgmodl in gdat.liststrgmodl:
        for e in gmod.indxsersfgrd:
            strgsersfgrd = 'isf%d' % e
            # scalar
            setp_varb(gdat, 'masshost%stotl' % strgsersfgrd, limt=[1e7, 1e14], strgmodl=strgmodl, labl=[r'$M_{\rm{hst,%d}}$' % e, ''])
            setp_varb(gdat, 'masshost%sbein' % strgsersfgrd, limt=[1e7, 1e14], strgmodl=strgmodl, labl=[r'$M_{\rm{hst,E,%d}}$' % e, ''])
            
            # defined on a grid
            setp_varb(gdat, 'masshost%sdelt' % strgsersfgrd, limt=[1e7, 1e14], strgmodl=strgmodl, labl=[r'$dM_{\rm{hst,%d}}$' % e, ''])
            setp_varb(gdat, 'masshost%sintg' % strgsersfgrd, limt=[1e7, 1e14], strgmodl=strgmodl, labl=[r'$M_{\rm{hst,%d<}}$' % e, ''])
            
            if gmod.numbpopl > 0:
                if gmod.boollenssubh:
                    # to be deleted
                    #setattr(gmod.lablrootpara, 'masssubh%sintg', r'$M_{\rm{sub}}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sdelt', r'$\rho_{\rm{sub}}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sintgbein', r'$M_{\rm{sub,E}}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sdeltbein', r'$\rho_{\rm{sub,E}}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sintgunit', '$10^9 M_{\odot}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sdeltunit', '$M_{\odot}$/kpc')
                    #setattr(gmod.lablrootpara, 'masssubh%sintgbeinunit', '$10^9 M_{\odot}$')
                    #setattr(gmod.lablrootpara, 'masssubh%sdeltbeinunit', '$M_{\odot}$/kpc')
                    #setattr(gmod.lablrootpara, 'fracsubh%sintg', r'f_{\rm{sub}}')
                    #setattr(gmod.lablrootpara, 'fracsubh%sdelt', r'f_{\rho,\rm{sub}}')
                    #setattr(gmod.lablrootpara, 'fracsubh%sintgbein', r'$f_{\rm{sub,E}}$')
                    #setattr(gmod.lablrootpara, 'fracsubh%sdeltbein', r'$f_{\rho,\rm{sub,E}}$')
                    
                    # scalar
                    setp_varb(gdat, 'masssubh%sbein' % strgsersfgrd, limt=[1e7, 1e10], strgmodl=strgmodl, labl=[r'$M_{\rm{sub,E,%d}}$' % e, ''])
                    setp_varb(gdat, 'fracsubh%sbein' % strgsersfgrd, limt=[0., 1.], strgmodl=strgmodl, labl=[r'$f_{\rm{sub,E,%d}}$' % e, ''])
                   
                    # defined on a grid
                    setp_varb(gdat, 'masssubh%sdelt' % strgsersfgrd, limt=[1e7, 1e10], strgmodl=strgmodl, labl=[r'$dM_{\rm{sub,%d}}$' % e, ''])
                    setp_varb(gdat, 'fracsubh%sdelt' % strgsersfgrd, limt=[0., 1.], strgmodl=strgmodl, labl=[r'$df_{\rm{sub,%d}}$' % e, ''])
    
                    setp_varb(gdat, 'masssubh%sintg' % strgsersfgrd, limt=[1e7, 1e10], strgmodl=strgmodl, labl=[r'$M_{\rm{sub,%d<}}$' % e, ''])
                    setp_varb(gdat, 'fracsubh%sintg' % strgsersfgrd, limt=[0., 1.], strgmodl=strgmodl, labl=[r'$f_{\rm{sub,%d<}}$' % e, ''])
    
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
    gdat.pathplotcnfg = gdat.pathvisu + gdat.strgcnfg + '/'
    gdat.pathinit = gdat.pathplotcnfg + 'init/'
    gdat.pathinitintr = gdat.pathinit + 'intr/'
    
    if gdat.boolbindspat:
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
            for strgpdfn in gmod.scalpara.genrelem[l]:
                if strgpdfn.startswith('gaum') and gmod.xposprio is None and gmod.yposprio is None:
                    raise Exception('If typespatdist is "gaus", spatial coordinates of the prior catalog should be provided via xposprio and yposprio.')
    
    # temp -- have these definitions separate for all features
    # feature plotting factors and scalings
    gdat.dictglob = {}
    
    gdat.listnamechro = ['totl', 'prop', 'diag', 'save', 'plot', 'proc', 'pars', 'modl', 'llik', 'sbrtmodl']
    gdat.listlablchro = ['Total', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Process', 'Parse', 'Model', 'Likelihood', 'Total emission']
    if gmod.numbpopl > 0:
        gdat.listnamechro += ['spec']
        gdat.listlablchro += ['Spectrum calculation']
    if gmod.boollens:
        gdat.listnamechro += ['deflzero', 'deflhost', 'deflextr', 'sbrtlens', 'sbrthost']
        gdat.listlablchro += ['Array initialization', 'Lens Host Deflection', 'External deflection', 'Lensed emission', 'Lens Host emission']
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
    
    if gdat.typedata != 'simu':
        if gmod.boolelemlghtanyy and gdat.typeexpr == 'ferm' and gdat.maxmgangdata == 20. / gdat.anglfact:
            path = gdat.pathinpt + 'sbrt0018.png'
            gdat.sbrt0018 = sp.ndimage.imread(path, flatten=True)
            gdat.sbrt0018 -= np.amin(gdat.sbrt0018)
            gdat.sbrt0018 /= np.amax(gdat.sbrt0018)
            binsxpostemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[1])
            binsypostemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[0])
            gdat.sbrt0018objt = sp.interpolate.RectBivariateSpline(binsypostemp, binsxpostemp, gdat.sbrt0018)

    # log-prior register
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    numb = 0
    for l in gmod.indxpopl:
        numb += len(gmod.namepara.genrelem[l])
    
    if gmod.boollens or gdat.typedata == 'simu' and gmod.boollens:
        setp_varb(gdat, 'mcut')
        
        # 'temp' turn these back on as necessary
        ##setp_varb(gdat, 'bein')

        ## angular deviation
        setp_varb(gdat, 'anglhalf', minm=0., maxm=3*gdat.maxmgangdata, labl=['$\theta$', ''], numbbins=1000)
        setp_varb(gdat, 'anglfull', minm=0., maxm=3*gdat.maxmgangdata, numbbins=1000)
        
    setp_varb(gdat, 'anglfromhost', minm=0., maxm=3*gdat.maxmgangdata, numbbins=1000, labl=['$\theta_{\rm{0,hst}}$', ''])
    
    # temp
    #gdat.blimpara.anglcosi = np.sort(np.cos(gdat.blimpara.angl))
    
    # temp
    #gdat.meshbackener = np.meshgrid(gdat.gmod.indxback, gdat.indxener, indexing='ij')
    
    # plotting
    ## the normalized offset for text annotation of point sources in the frames
    gdat.offstextimag = gdat.maxmgangdata * 0.05
    
    ## figure size
    gdat.plotsize = 4.5
    ## size of the images
    gdat.sizeimag = 1.3 * gdat.plotsize
    
    ## label of the models
    gdat.fitt.lablmodl = 'Model'
    if gdat.typedata == 'simu':
        gdat.refr.lablmodl = 'True'
    else:
        gdat.refr.lablmodl = 'Ref'
    
    # element parameters common between the fitting and reference models
    # 'temp' check later what this is useful for 
    #gdat.namepara.elemcomm = [[[] for l in gmod.indxpopl] for q in gdat.indxrefr]
    #for q in gdat.indxrefr:
    #    for l in gmod.indxpopl:
    #        for strgfeat in gmod.listnameparatotlelem[l]:
    #            if strgfeat in gdat.refr.namepara.elem[q]:
    #                gdat.namepara.elemcomm[q][l].append(strgfeat)
    
    if gdat.typedata == 'simu':
        gdat.refr.indxpopl = gdat.true.indxpopl
        gdat.refr.lablpopl = gdat.true.lablpopl

    
    for strgmodl in ['refr', 'fitt']:
    
        gmod = getattr(gdat, strgmodl)
        
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
    if gdat.numbdqlt == 1:
        gdat.indxdqltplot = gdat.indxdqlt
    else:
        gdat.indxdqltplot = np.concatenate((np.array([-1]), gdat.indxdqlt))
    
    gdat.numbenerdqlt = gdat.numbener * gdat.numbdqlt
    
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
        gdat.xposgrid = gdat.xposgrid[gdat.indxpixlkeep]
        gdat.yposgrid = gdat.yposgrid[gdat.indxpixlkeep]
        gdat.indxpixlfull = gdat.indxpixlfull[gdat.indxpixlkeep]
        
    # the function to measure time
    # temp
    gdat.strgfunctime = 'clck'
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

    ## longitude
    gdat.numbxpospntsprob = gdat.numbsidepntsprob
    gdat.numbypospntsprob = gdat.numbsidepntsprob
    gdat.blimpara.xpospntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.blimpara.ypospntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.indxxpospntsprob = np.arange(gdat.numbxpospntsprob)
    gdat.indxypospntsprob = np.arange(gdat.numbypospntsprob)

    # lensing problem setup
    ## number of deflection components to plot

    gdat.blimpara.xposcartmesh, gdat.blimpara.yposcartmesh = np.meshgrid(gdat.blimpara.xposcart, gdat.blimpara.yposcart, indexing='ij')
    gdat.bctrpara.xposcartmesh, gdat.bctrpara.yposcartmesh = np.meshgrid(gdat.bctrpara.xposcart, gdat.bctrpara.yposcart, indexing='ij')
    if gdat.typepixl == 'cart':
        gdat.sizepixl = np.sqrt(gdat.apix)
        gdat.indxsidecart = np.arange(gdat.numbsidecart)
        gdat.indxpixlrofi = np.arange(gdat.numbpixlcart)
        gdat.indxsidemesh = np.meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
        gdat.xposgrid = gdat.bctrpara.xposcart[gdat.indxsidemesh[0].flatten()]
        gdat.yposgrid = gdat.bctrpara.yposcart[gdat.indxsidemesh[1].flatten()]
        gdat.shapcart = (gdat.numbsidecart, gdat.numbsidecart)
        gdat.xposgridfull = np.copy(gdat.xposgrid)
        gdat.yposgridfull = np.copy(gdat.yposgrid)
        gdat.xposgridcart = gdat.xposgrid.reshape(gdat.shapcart)
        gdat.yposgridcart = gdat.yposgrid.reshape(gdat.shapcart)
        gdat.indxpent = np.meshgrid(gdat.indxener, gdat.indxsidecart, gdat.indxsidecart, gdat.indxdqlt, indexing='ij')
    if gdat.typepixl == 'heal':
        xposheal, yposheal, gdat.numbpixlfull, gdat.apix = tdpy.retr_healgrid(gdat.numbsideheal)
        xposheal = np.deg2rad(xposheal)
        yposheal = np.deg2rad(yposheal)
   
        gdat.indxpixlrofi = np.where((np.fabs(xposheal) < gdat.maxmgangdata) & (np.fabs(yposheal) < gdat.maxmgangdata))[0]
        
        gdat.indxpixlrofimarg = np.where((np.fabs(xposheal) < 1.2 * gdat.maxmgangdata) & (np.fabs(yposheal) < 1.2 * gdat.maxmgangdata))[0]

        gdat.xposgrid = xposheal
        gdat.yposgrid = yposheal
    
    gdat.indxpixlfull = np.arange(gdat.numbpixlfull)
    if gdat.typepixl == 'cart':
        gdat.indxpixlcart = np.arange(gdat.numbpixlcart)
    
    if gdat.boolbinddqlt:
        # PSF class string
        gdat.strgdqlt = []
        for m in gdat.indxdqlt:
            gdat.strgdqlt.append('PSF%d' % gdat.indxdqltincl[m])
    
    # power spectra
    if gdat.typepixl == 'cart':
        setp_varb(gdat, 'wvecodim', minm=0., maxm=1., boolinvr=True)
        setp_varb(gdat, 'wlenodim', minm=0., maxm=1., boolinvr=True)
        setp_varb(gdat, 'anglodim', minm=0., maxm=1., boolinvr=True)
        setp_varb(gdat, 'mpolodim', minm=0., maxm=1.)
        gdat.numbmpolodim = gdat.bctrpara.mpolodim.size
        gdat.indxmpolodim = np.arange(gdat.numbmpolodim)
        #setp_varb(gdat, 'anglodim', boolinvr=True)
        #setp_varb(gdat, 'mpolodim')
            
        for strgmodl in gdat.liststrgmodl:
            gmod = getattr(gdat, strgmodl)
            
            gdat.numbwvecodim = gdat.numbsidecart
            gdat.minmpara.anglodim = 0.
            gdat.maxmpara.anglodim = 2. * gdat.maxmgangdata
            gdat.minmpara.mpolodim = 0.
            gdat.maxmpara.mpolodim = 1. / 2. / gdat.sizepixl

            if gmod.boollens or gdat.typedata == 'simu' and gmod.boollens:
                # temp -- this should minima, maxima of adislens and the true metamodel into account
                gdat.minmpara.wvecodim = gdat.minmpara.mpolodim / np.amax(gmod.adislens)
                gdat.maxmpara.wvecodim = gdat.maxmpara.mpolodim / np.amin(gmod.adislens)
                gdat.minmpara.wlenodim = gdat.minmpara.anglodim * np.amin(gmod.adislens)
                gdat.maxmpara.wlenodim = gdat.maxmpara.anglodim * np.amax(gmod.adislens)
                setp_varb(gdat, 'wvecodim', strgmodl=strgmodl)
                setp_varb(gdat, 'wlenodim', strgmodl=strgmodl, boolinvr=True)
                gdat.bctrpara.wvecxpos, gdat.bctrpara.wvecypos = np.meshgrid(gdat.bctrpara.wvecodim, gdat.bctrpara.wvecodim, indexing='ij')
                gdat.bctrpara.wvec = np.sqrt(gdat.bctrpara.wvecxpos**2 + gdat.bctrpara.wvecypos**2)
        gdat.bctrpara.mpolxpos, gdat.bctrpara.mpolypos = np.meshgrid(gdat.bctrpara.mpolodim, gdat.bctrpara.mpolodim, indexing='ij')
        gdat.bctrpara.mpol = np.sqrt(gdat.bctrpara.mpolxpos**2 + gdat.bctrpara.mpolypos**2)

    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
        
        # element parameter vector indices
        gmod.indxpara.genrelemxpos = 0
        gmod.indxpara.genrelemypos = 1
        gmod.indxpara.genrelemflux = 2
        gmod.indxpara.genrelemsind = 3
        gmod.indxpara.genrelemcurv = 4
        gmod.indxpara.genrelemexpc = 4

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
    gdat.indxcubeincl = np.meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxdqltincl, indexing='ij')
    
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
    #    gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbdqlt))
    #    for strgmodl in gdat.liststrgmodl:
    #        gmod.sbrtbacknormcart = []
    #        for c in getattr(gmod, 'gmod.indxback'):
    #            gmod.sbrtbacknormcart.append(gmod.sbrtbacknorm[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbdqlt)))
    
    # mask the exposure map
    if gdat.listmask is not None:
        for mask in gdat.listmask:
            if mask[0] == 'sqre':
                indxpixlmask = np.where((gdat.xposgrid > mask[1]) & (gdat.xposgrid < mask[2]) & (gdat.yposgrid > mask[3]) & (gdat.yposgrid < mask[4]))[0]
            if mask[0] == 'circ':
                indxpixlmask = np.where(np.sqrt((gdat.xposgrid - mask[1])**2 + (gdat.yposgrid - mask[2])**2) < mask[3])[0]
            if mask[0] == 'hstr':
                indxpixlmask = np.where((gdat.yposgrid > mask[1]) & (gdat.yposgrid < mask[2]))[0]
            if gdat.typemaskexpo == 'zero':
                gdat.expo[:, indxpixlmask, :] = 0.
            if gdat.typemaskexpo == 'ignr':
                gdat.expo[:, indxpixlmask, :] = 1e-49

    # plotting
    ## ROI
    if gdat.boolbindspat:
        gdat.exttrofi = np.array([gdat.minmxposdata, gdat.maxmxposdata, gdat.minmyposdata, gdat.maxmyposdata])
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
    
    gdat.scalpara.cntpresi = 'asnh'
    gdat.cmappara.cntpresi = make_cmapdivg('Red', 'Orange')

    setp_varb(gdat, 'conv', minm=1e-2, maxm=10., labl=['$\kappa$', ''], cmap='Purples', scal='logt')
    setp_varb(gdat, 's2nr', minm=0., maxm=10., labl=['SNR', ''], cmap='magma', scal='asnh')
    setp_varb(gdat, 'magn', minm=-1e2, maxm=1e2, labl=['$\mu$', ''], cmap='BrBG', scal='asnh')
    
    gdat.minmdeflresiperc = -100.
    gdat.maxmdeflresiperc = 100.
    gmod.scaldeflresiperc = 'self'
    gdat.cmapdeflresiperc = 'Oranges'
    
    setp_varb(gdat, 'convelem', minm=1e-4, maxm=1e-1, labl=['$C_{el}$', ''], cmap='Purples', scal='logt')
    setp_varb(gdat, 'convelemresi', minm=-0.1, maxm=0.1, labl=['$C_{el} - $C_{el}$', ''], cmap='PiYG', scal='self')
    setp_varb(gdat, 'convelemresiperc', minm=-100., maxm=100., labl=['$C_{el} - $C_{el}$', ''], cmap='PiYG', scal='self')
    
    gdat.minmmagnresi = -10.
    gdat.maxmmagnresi = 10.
    gmod.scalmagnresi = 'self'
    gdat.cmapmagnresi = 'PRGn'
    
    gdat.minmmagnresiperc = -100.
    gdat.maxmmagnresiperc = 100.
    gmod.scalmagnresiperc = 'self'
    gdat.cmapmagnresiperc = 'PRGn'
    
    gdat.xposgrid = gdat.xposgrid[gdat.indxpixlrofi]
    gdat.yposgrid = gdat.yposgrid[gdat.indxpixlrofi]
   
    if gdat.boolcorrexpo:
        if np.amax(gdat.expo) <= 0.:
            raise Exception('Bad exposure.')

    # temp
    #gdat.expo[np.where(gdat.expo < 1e-50)] = 1e-50
    
    # exclude voxels with vanishing exposure
    if gdat.boolcorrexpo:
        for i in gdat.indxener:
            for m in gdat.indxdqlt:
                gdat.indxpixlrofi = np.intersect1d(gdat.indxpixlrofi, np.where(gdat.expo[i, :, m] > 0.)[0])
    
    gdat.indxcuberofi = np.meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxdqlt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = np.arange(gdat.numbpixl)
    gdat.numbdata = gdat.numbener * gdat.numbdqlt * gdat.numbpixl

    #gdat.xposgridrofi = gdat.xposgrid[gdat.indxpixlrofi]
    #gdat.yposgridrofi = gdat.yposgrid[gdat.indxpixlrofi]


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
        if gdat.typeexpr == 'gmix':
            gdat.maxmangl = 1.
        if gdat.typeexpr == 'ferm':
            gdat.maxmangl = 20. / gdat.anglfact
        if gdat.typeexpr == 'tess':
            gdat.maxmangl = 25. / gdat.anglfact
        if gdat.typeexpr == 'chan':
            gdat.maxmangl = 15. / gdat.anglfact
        if gdat.typeexpr.startswith('HST_WFC3'):
            gdat.maxmangl = 1. / gdat.anglfact
    else:
        gdat.maxmangl = gdat.maxmgangdata * np.sqrt(2.) * 2. * 1.1
        
    gdat.listnamespatmean = ['full']
    if gdat.typeexpr == 'ferm':
        gdat.listnamespatmean += ['innr']
    gdat.numbspatmean = len(gdat.listnamespatmean)
    gdat.indxspatmean = np.arange(gdat.numbspatmean)
    gdat.listindxcubespatmean = [[] for b in gdat.indxspatmean]
    gdat.indxcube = np.meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxdqlt, indexing='ij')
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        if namespatmean == 'full':
            gdat.listindxcubespatmean[b] = gdat.indxcube
        if namespatmean == 'innr':
            gdat.indxpixlinnr = np.where(np.sqrt(gdat.xposgrid**2 + gdat.yposgrid**2) < 5. / gdat.anglfact)[0]
            gdat.listindxcubespatmean[b] = np.meshgrid(gdat.indxener, gdat.indxpixlinnr, gdat.indxdqlt, indexing='ij')
    
    if gdat.numbpixl > 1:
        # store pixels as unit vectors
        gdat.xdatgrid, gdat.ydatgrid, gdat.zaxigrid = retr_unit(gdat.xposgrid, gdat.yposgrid)
   
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
                    dist = retr_angldistunit(gdat, xposheal[gdat.indxpixlrofimarg[k]], yposheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                    gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
                fobj = open(path, 'wb')
                pickle.dump(gdat.pixlcnvt, fobj, protocol=pickle.HIGHEST_PROTOCOL)
                fobj.close()
        
        # dummy pixel indices for full (nonlocal) element kernel evaluation 
        gdat.listindxpixl = []
        if gdat.typedata == 'simu':
            numb = max(gdat.true.maxmpara.numbelemtotl, gdat.fitt.maxmpara.numbelemtotl)
        else:
            numb = gdat.fitt.maxmpara.numbelemtotl
        numb += 2
        for k in range(int(numb)):
            gdat.listindxpixl.append([])
            for kk in range(k):
                gdat.listindxpixl[k].append(gdat.indxpixl)
        
        # spatial averaging setup
        # temp
        
        # temp -- check if 1000 is too much
        gdat.numbanglelem = 1000
    
    for namepara in gdat.fitt.namepara.glob:
        setattr(gdat.labltotlpara, namepara, getattr(gdat.fitt.labltotlpara, namepara))

    # set parameter features common between true and fitting models
    for strgmodl in gdat.liststrgmodl:
    
        gmod = getattr(gdat, strgmodl)
        
        # 'temp' fix this part later
        #for namepara in gmod.namepara.kind:
        #    
        #    try:
        #        getattr(gdat.minmpara, namepara)
        #        getattr(gdat.maxmpara, namepara)
        #    except:
        #        try:
        #            setattr(gdat.minmpara, namepara, min(getattr(gdat.fitt.minmpara, namepara), getattr(gdat.true.minmpara, namepara)))
        #            setattr(gdat.maxmpara, namepara, max(getattr(gdat.fitt.maxmpara, namepara), getattr(gdat.true.maxmpara, namepara)))
        #        except:
        #            try:
        #                setattr(gdat.minmpara, namepara, getattr(gdat.fitt.minmpara, namepara))
        #                setattr(gdat.maxmpara, namepara, getattr(gdat.fitt.maxmpara, namepara))
        #            except:
        #                setattr(gdat.minmpara, namepara, getattr(gdat.true.minmpara, namepara))
        #                setattr(gdat.maxmpara, namepara, getattr(gdat.true.minmpara, namepara))
    
    # set plot limits for each model if not already set (for Gaussian, log-normal distributions)
    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
        
        for namepara in gmod.namepara.kind:
            if hasattr(gmod.minmpara, namepara):
                minm = getattr(gmod.minmpara, namepara)
                maxm = getattr(gmod.maxmpara, namepara)
                limt = np.array([minm, maxm])
                setattr(gmod.limtpara, namepara, limt)
    
    if gdat.checprio:
        gdat.liststrgpdfn = ['prio', 'post']
    else:
        gdat.liststrgpdfn = ['post']

    if gdat.typeverb > 1:
        # temp
        for strgmodl in gdat.liststrgmodl:
            print('strgmodl')
            print(strgmodl)
            print('Fixed dimensional parameters:')
            print('%20s%25s%5s%20s%20s' % ('name', 'labltotl', 'scal', 'minm', 'maxm'))
            for k in gmod.indxpara.genrbase:
                print('%20s%25s%5s%20.6g%20.6g' % (gmod.namepara.genrbase[k], gmod.labltotlpara.genrbase[k], gmod.scalpara.genrbase[k], \
                                                                              gmod.minmpara.genrbase[k], gmod.maxmpara.genrbase[k]))
            
            print('Element parameters')
            print('%20s%20s' % ('nameparagenrelem', 'scalcomp'))
            for l in gmod.indxpopl:
                for nameparagenrelem, scalcomp in zip(gmod.namepara.genrelem[l], gmod.scalpara.genrelem[l]):
                    print('%20s%20s' % (nameparagenrelem, scalcomp))
            
            print('%20s%20s' % ('strgmodu', 'pdfnmodu'))
            for l in gmod.indxpopl:
                for strgmodu, pdfnmodu in zip(gmod.namepara.genrelemmodu[l], gmod.liststrgpdfnmodu[l]):
                    print('%20s%20s' % (strgmodu, pdfnmodu))
            
            print('%20s%20s' % ('strgfeat', 'pdfnprio'))
            for l in gmod.indxpopl:
                for strgfeat, pdfnprio in zip(gmod.namepara.genrelem[l], gmod.scalpara.genrelem[l]):
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
        if gmod.numbpopl > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
    if gdat.probspmr is None:
        if gmod.numbpopl > 0:
            gdat.probspmr = gdat.probtran / 2.
        else:
            gdat.probspmr = 0.
    
    gdat.probbrde = 1. - gdat.probspmr

    if gdat.probbrde < 0:
        raise Exception('')
    gdat.lablproptype = ['Within']
    gdat.nameproptype = ['with']
    if gmod.numbpopl > 0:
        gdat.lablproptype += ['Birth', 'Death', 'Split', 'Merge']
        gdat.nameproptype += ['brth', 'deth', 'splt', 'merg']
    gdat.numbproptype = len(gdat.lablproptype)
    gdat.nameproptype = np.array(gdat.nameproptype)
    cntr = tdpy.cntr()
    if gmod.numbpopl > 0.:
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
    if gdat.factpriodoff != 1.:
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
            
            for strgfeat, strgpdfn in zip(gmod.namepara.genr, gmod.scalpara.genr):
                if strgpdfn == 'tmplgrad':
                    pdfnpriotemp = np.empty((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                    lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                    lpdfpriointp = lpdfprioobjt(gdat.bctrpara.yposcart, gdat.bctrpara.xposcart)
        
    gdat.indxpoplcrin = 0
    if gmod.numbpopl > 0:
        if gdat.strgcnfgsimu is not None:
            path = gdat.pathoutpcnfgsimu + 'gdatfinlpost'
            gdatsimu = readfile(path)
        gdat.liststrgvarbhist = []
        cntr = 0
        for l0 in gmod.indxpopl:
            for a, strgfeatfrst in enumerate(gmod.namepara.genrelem[l0]):
                if strgfeatfrst == 'spec':
                    continue
                gdat.liststrgvarbhist.append([[] for k in range(5)])
                gdat.liststrgvarbhist[cntr][0] = 'hist' + strgfeatfrst + 'pop%d' % l
                gdat.liststrgvarbhist[cntr][1] = strgfeatfrst
                if gdat.strgcnfgsimu is not None:
                    # cmpl
                    gdat.liststrgvarbhist[cntr][3] = [[] for qq in gdatsimu.indxrefr]
                    # fdis
                    gdat.liststrgvarbhist[cntr][4] = [[] for qq in gdatsimu.indxrefr]
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
                    if gdat.strgcnfgsimu is not None:
                        booltempfrst = True
                        booltempseco = True
                        if strgfeatfrst[-4:] in gdat.listnamerefr:
                            q = gdat.listnamerefr.index(strgfeatfrst[-4:])
                            booltempfrst = not strgfeatfrst in gdat.refr.namepara.elemonly[q][l]
                        if strgfeatseco[-4:] in gdat.listnamerefr:
                            q = gdat.listnamerefr.index(strgfeatseco[-4:])
                            booltempseco = not strgfeatseco in gdat.refr.namepara.elemonly[q][l]
                        for qq in gdatsimu.indxrefr:
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
    if gdat.typedata == 'inpt' and gmod.numbpopl > 0:
        if gdat.numbsampboot is None:
            gdat.numbsampboot = gdat.numbsamp
    
        gdat.boolcrex = False
        if gdat.strgcnfgsimu is not None:
            for qq in gdatsimu.indxrefr:
                for q in gdat.indxrefr:
                    for l in gmod.indxpopl:
                        for strgfeatfrst in gmod.namepara.genrelem[l]:
                            
                            if gdat.typeexpr == 'chan' and strgfeatfrst == 'redswo08':
                                crex = (1. + gdat.bctrpara.redswo08)**2
                            else:
                                crex = None
                            
                            setattr(gdat, 'crex' + strgfeatfrst + 'pop%dpop%dpop%d' % (q, qq, l), crex)
                            
                            for strgfeatseco in gmod.namepara.genrelem[l]:
                                
                                if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                    continue
                                            
                                if gdat.typeexpr == 'chan' and (strgfeatfrst == 'redswo08' or strgfeatseco == 'redswo08'):
                                    crex = np.empty((gdat.numbbinsplot, gdat.numbbinsplot))
                                    if strgfeatfrst == 'redswo08':
                                        crex[:, :] = (1. + gdat.bctrpara.redswo08[:, None])**2
                                    else:
                                        crex[:, :] = (1. + gdat.bctrpara.redswo08[None, :])**2
                                else:
                                    crex = None
                                
                                setattr(gdat, 'crex' + strgfeatfrst + strgfeatseco + 'pop%dpop%dpop%d' % (q, qq, l), crex)
    
            if gdat.refr.numbelemtotl > 0:
                for listtemp in gdat.liststrgvarbhist:
                    strgvarb = listtemp[0]
                    for qq in gdatsimu.indxrefr:
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
        gdat.boolcrin = gdat.typedata == 'inpt' and gdat.strgcnfgsimu is not None
    
    if gmod.numbpopl > 0:
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
    if gmod.numbpopl > 0 or gdat.typedata == 'simu' and gmod.numbpopl > 0:
        gdat.liststrgfoldinit += ['', 'histodim/', 'histtdim/', 'scattdim/', 'cmpltdim/']
    gdat.liststrgfoldfram = ['']
    if gmod.numbpopl > 0:
        gdat.liststrgfoldfram += ['scattdim/']
    gdat.liststrgfoldfinl = ['']
    if gdat.boolinforefr and gmod.numbpopl > 0:
        gdat.liststrgfoldfram += ['assc']
        gdat.liststrgfoldfinl += ['assc']
    gdat.liststrgfoldanim = deepcopy(gdat.liststrgfoldfram)

    if gmod.numbpopl > 0:
        for strgdims in ['odim/', 'tdim/']:
            for strgelemtdimvarb in gdat.liststrgelemtdimvarbfram:
                gdat.liststrgfoldfram += [strgelemtdimvarb + strgdims]
            for strgelemtdimvarb in gdat.liststrgelemtdimvarbfinl:
                gdat.liststrgfoldfinl += [strgelemtdimvarb + strgdims]

    # make folders
    #gdat.pathprio = gdat.pathplotcnfg + 'prio/'
    #gdat.pathpost = gdat.pathplotcnfg + 'post/'
    make_fold(gdat)

    if gdat.strgcnfgsimu is not None:
        if gdat.typedata == 'inpt':
            path = gdat.pathoutpcnfgsimu + 'gdatfinlpost'
            booltemp = True
            gdatsimu = readfile(path)

            if booltemp:
                numbparaelem = gdatsimu.true.numbparaelem
                if gdatsimu.trueindxpopl != gmod.indxpopl:
                    raise Exception('')
                for l in gmod.indxpopl:
                    for strgfeat in gmod.namepara.genrelem[l]:
                        if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':
                            continue

                        if strgfeat[-4:] in gdat.listnamerefr:
                            continue
                        reca = getattr(gdatsimu.true, 'reca' + strgfeat + 'pop%d' % l)
                        setattr(gdat.true, 'reca' + strgfeat + 'pop%d' % l, reca)
                gmod.namepara.genr.elem = gdatsimu.truegmod.namepara.genr.elem
    
    setp_varb(gdat, 'angl', minm=0., maxm=10., numbbins=10)
    
    if gmod.typeelemspateval[l] == 'locl' and gmod.numbpopl > 0:
        gdat.numbprox = 3
        gdat.indxprox = np.arange(gdat.numbprox)
        minmparagenrscalelemampl = getattr(gdat.fitt.minmpara, gmod.nameparagenrelemampl[0])
        maxmparagenrscalelemampl = getattr(gdat.fitt.maxmpara, gmod.nameparagenrelemampl[0])
        gdat.blimpara.prox = np.logspace(np.log10(minmparagenrscalelemampl), np.log10(maxmparagenrscalelemampl), gdat.numbprox + 1)
        
        # determine the maximum angle at which the contribution of the element will be computed
        if gdat.boolbindspat:
            if gdat.maxmangleval is None:
                if gdat.typeexpr == 'chan':
                    gdat.maxmangleval = np.array([5., 6., 9.]) / gdat.anglfact
                elif gdat.typeexpr == 'gmix':
                    gdat.maxmangleval = np.array([0.1, 0.2, 0.3]) / gdat.anglfact
                elif gdat.typeexpr == 'ferm':
                    gdat.maxmangleval = np.array([7., 9., 15.]) / gdat.anglfact
                else:
                    gdat.maxmangleval = np.empty(gdat.numbprox)
                    for h in gdat.indxprox:
                        if gdat.specfraceval == 0:
                            gdat.maxmangleval[h] = 3. * gdat.maxmgang
                        else:  
                            frac = min(1e-2, gdat.specfraceval * gdat.blimpara.prox[0] / gdat.blimpara.prox[h+1])
                            psfnwdth = retr_psfnwdth(gdat, gmodstat.psfn, frac)
                            gdat.indxmaxmangl = np.unravel_index(np.argmax(psfnwdth), psfnwdth.shape)
                            gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.typeverb > 1:
            if gmod.typeelemspateval == 'locl':
                print('maxmangleval')
                print(gdat.anglfact * gdat.maxmangleval[l], ' [%s]' % gdat.strganglunit)

        if gdat.boolelempsfnanyy and gdat.maxmpara.angl < np.amax(gdat.maxmangleval):
            print('')
            print('')
            print('')
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
            print('gdat.minmxposdata')
            print(gdat.minmxposdata)
            print('gdat.minmyposdata')
            print(gdat.minmyposdata)
            print('gdat.maxmxposdata')
            print(gdat.maxmxposdata)
            print('gdat.maxmyposdata')
            print(gdat.maxmyposdata)
        if gdat.typeverb > 0:
            print('Element evaluation will be performed up to')
            if gdat.boolbindspat:
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
                dist = retr_angldistunit(gdat, gdat.xposgrid[j], gdat.yposgrid[j], gdat.indxpixl)
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
    if gdat.typedata == 'simu':
        for q in gdat.indxrefr:
            for strgfeat in gmod.namepara.genrelem[q]:
                booltemp = False
                for l in gmod.indxpopl:
                    if strgfeat in gmod.namepara.genrelem[l]:
                        booltemp = True
                if not booltemp:
                    setattr(gdat.minmpara, 'minm' + strgfeat + gdat.listnamerefr[q], getattr(gdat.true.minm, strgfeat))
                    setattr(gdat.maxmpara, 'maxm' + strgfeat + gdat.listnamerefr[q], getattr(gdat.true.maxm, strgfeat))

    ## reference spectra
    if gdat.listprefsbrtlabltotl is None:
        if gdat.typeexpr == 'chan' and gdat.boolbindspat:
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
        #    gdat.enerexpcfact = gdat.enerpivt - gdat.bctrpara.ener
        #if gmod.numbpopl > 0 and gdat.numbener > 1:
        #    minmsinddistmeanpop0 = getattr(gmod, 'minmsinddistmeanpop0')
        #    factspecener = (gdat.bctrpara.ener / gdat.enerpivt)**(-np.sqrt(np.amin(minmsinddistmeanpop0) * np.amax(maxmsinddistmeanpop0)))
    else:
        pass
        #gdat.factspecener = np.array([1.])

    # temp -- this assumes square ROI
    if gdat.boolbindspat:
        gdat.frambndrmodl = gdat.maxmxposdata * gdat.anglfact
    
    if gmod.boollenshost or gdat.typedata == 'simu' and gmod.boollenshost:
        
        if gdat.typesers == 'intp':
            # construct pixel-convolved Sersic surface brightness template
            gdat.factsersusam = 10
            maxmxpos = 4. * np.sqrt(2.) * gdat.maxmxpos
            gdat.numbxpossers = int(np.ceil(maxmxpos / gdat.sizepixl))
            gdat.numbxpossersusam = (1 + gdat.numbxpossers) * gdat.factsersusam
            setp_varb(gdat, 'xpossers')
            setp_varb(gdat, 'xpossersusam')
            setp_varb(gdat, 'ypossersusam')
            
            gdat.numbhalfsers = 20
            gdat.numbindxsers = 20
                
            setp_varb(gdat, 'halfsers')
            setp_varb(gdat, 'indxsers')
            
            gdat.blimpara.xpossersusammesh, gdat.blimpara.ypossersusammesh = np.meshgrid(gdat.blimpara.xpossersusam, gdat.blimpara.ypossersusam, indexing='ij')
            gdat.blimpara.radisersusam = np.sqrt(gdat.blimpara.xpossersusammesh**2 + gdat.blimpara.ypossersusammesh**2)
             
            gdat.sersprofcntr = np.empty((gdat.numbxpossers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
            gdat.sersprof = np.empty((gdat.numbxpossers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
            
            for n in range(gdat.numbindxsers + 1):
                for k in range(gdat.numbhalfsers + 1):
                    
                    profusam = retr_sbrtsersnorm(gdat.blimpara.radisersusam, gdat.blimpara.halfsers[k], indxsers=gdat.blimpara.indxsers[n])
    
                    ## take the pixel average
                    indxyposlowr = gdat.factsersusam * (gdat.numbxpossers + 1) / 2
                    indxyposuppr = gdat.factsersusam * (gdat.numbxpossers + 3) / 2
                    for a in range(gdat.numbxpossers):
                        indxxposlowr = gdat.factsersusam * a
                        indxxposuppr = gdat.factsersusam * (a + 1) + 1
                        gdat.sersprofcntr[a, k, n] = profusam[(indxxposlowr+indxxposuppr)/2, 0]
                        gdat.sersprof[a, k, n] = np.mean(profusam[indxxposlowr:indxxposuppr, :])
            
            temp, indx = unique(gdat.blimpara.xpossers, return_index=True)
            gdat.blimpara.xpossers = gdat.blimpara.xpossers[indx]
            gdat.sersprof = gdat.sersprof[indx, :, :]
            gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]
    
            indx = np.argsort(gdat.blimpara.xpossers)
            gdat.blimpara.xpossers = gdat.blimpara.xpossers[indx]
            gdat.sersprof = gdat.sersprof[indx, :, :]
            gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]

    #for strg, valu in gmod.cmappara.__dict__.items():
    #    retr_ticklabl(gdat, strg)
            
    # types of two dimensional distributions
    ## 'bind' and 'scat'
    gdat.liststrgelemtdimtype = ['bind']

    # generate true data
    if gdat.typedata == 'simu':
        
        if gdat.typeverb > 0:
            print('Generating simulated data...')

        if gdat.typeseed == 'rand':
            np.random.seed()
        else:
            if gdat.typeverb > 0:
                print('Setting the seed for the RNG to %d...' % gdat.typeseed)
            np.random.seed(gdat.typeseed)
    
        gdat.true.this.indxpara = tdpy.gdatstrt()
        
        ## unit sample vector
        gdat.true.this.paragenrunitfull = np.random.rand(gdat.true.numbparagenr)
        gdat.true.this.paragenrscalfull = np.zeros(gdat.true.numbparagenr)
        
        if gdat.typeverb > 0:
            show_paragenrscalfull(gdat, None, strgmodl='true')

        if gdat.true.numbpopl > 0:
            gdat.true.this.numbelempopl = np.empty(gdat.true.maxmpara.numbelem[l], dtype=int)
            for l in gdat.true.indxpopl:
                gdat.true.this.paragenrunitfull[gdat.true.indxpara.numbelem[l]] = getattr(gdat.true.this, 'numbelempop%d' % l)
                gdat.true.this.numbelempopl[l] = getattr(gdat.true.this, 'numbelempop%d' % l)

            gdat.true.this.indxelemfull = [[] for l in gdat.true.indxpopl]
            for l in gdat.true.indxpopl:
                gdat.true.this.indxelemfull[l] = list(range(gdat.true.numbelem[l]))

            gdat.true.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gdat.true.this.indxelemfull, 'true')
        else:
            gdat.true.this.indxelemfull = []
            gdat.true.this.indxparagenrelemfull = None

        if gdat.true.numbpopl > 0:
            if gdat.typeseedelem is None:
                np.random.seed()
            else:
                np.random.seed(gdat.typeseedelem)
            gdat.true.this.paragenrunitfull[gdat.true.numbparagenrbase:] = np.random.rand(gdat.true.numbparagenrelem)
        
        gdat.true.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'true', gdat.true.this.paragenrunitfull, gdat.true.this.indxparagenrelemfull)
        
        # impose true values (valu)
        for k in gdat.true.indxpara.genr:
            
            if gdat.true.numbpopl > 0 and (k in gdat.true.indxpara.numbelem or \
                                        gdat.true.typemodltran == 'pois' and k in gdat.true.indxpara.meanelem):
                    continue
    
            # assume the true PSF
            if gdat.true.typeevalpsfn != 'none' and gdat.numbpixl > 1 and k in gdat.true.indxpara.psfp:
                gdat.true.this.paragenrscalfull[k] = gdat.true.psfpexpr[k-gdat.true.indxpara.psfp[0]]
            elif hasattr(gdat.true.this, gdat.true.namepara.genrscalfull[k]):
                ## read input simulated model parameters
                # impose user-defined true parameter
                gdat.true.this.paragenrscalfull[k] = getattr(gdat.true.this, gdat.true.namepara.genrscalfull[k])
    
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
                    for strgfeat in gmod.namepara.glob:
                    #for strgfeat in gmod.namepara.genrelem[l]:
                        if strgfeat[:-4] == 'etag':
                            continue
                        #setp_varb(gdat, strgfeat)
                        #if strgfeat in gmod.namepara.elem:
                        #    setp_varb(gdat, strgfeat + 'prio')
    
        proc_samp(gdat, None, 'this', 'true', boolinit=True)
    
        # set the reference model to true model
        print('Setting the reference model to the true model...')
        for strg, valu in gdat.true.this.__dict__.items():
            if strg == 'dictelem':
                # modify the current state of the element parameters of the true model to include uncertainty
                valutemp = [[] for l in gdat.true.indxpopl]
                for l in gdat.true.indxpopl:
                    valutemp[l] = dict()
                    for nameparaelem in gdat.true.this.dictelem[l]:
                        valutemp[l][nameparaelem] = np.zeros([3] + list(gdat.true.this.dictelem[l][nameparaelem].shape))
                        valutemp[l][nameparaelem][0, ...] = gdat.true.this.dictelem[l][nameparaelem]
            else:
                valutemp = valu
            setattr(gdat.refr, strg, valutemp)
        
        if gdat.boolmakeplot and gdat.boolmakeplotinit:
            plot_samp(gdat, None, 'this', 'true', 'init')
        
    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
        print('gmod.minmpara.numbelempop0')
        print(gmod.minmpara.numbelempop0)
        print('gmod.minmpara.numbelem')
        print(gmod.minmpara.numbelem)
        
    # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
    gmod.corr = tdpy.gdatstrt()
    for k in gmod.indxvarbscal:
        name = gmod.namepara.scal[k]
        temp = getattr(gdat.true.this, name)
        setattr(gmod.corr, name, temp)

    gmod.corrparagenrscalbase = np.empty(gmod.numbparagenrbase)
    for k in gmod.indxpara.genrbase:
        gmod.corrparagenrscalbase[k] = getattr(gdat.true, gmod.namepara.genrbase[k])

    dictglob = sample( \
                      **dictpcat, \
                     )










# new parameters to define
## relnpowr = 0.










def setp_paragenrscalbase(gdat, strgmodl='fitt'):
    '''
    Setup labels and scales for base parameters
    '''
    
    print('setp_paragenrscalbase(): Building the %s model base paremeter names and scales...' % strgmodl)
    gmod = getattr(gdat, strgmodl)
    
    # list of labels for background components
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
    gdat.numblablsbrt = 0
    for l in gmod.indxpopl:
        if gmod.boolelemsbrt[l]:
            listlablsbrt.append(gmod.lablpopl[l])
            listlablsbrt.append(gmod.lablpopl[l] + ' subt')
            gdat.numblablsbrt += 2
    if gmod.boollens:
        listlablsbrt.append('Source')
        gdat.numblablsbrt += 1
    if gmod.typeemishost != 'none':
        for e in gmod.indxsersfgrd:
            listlablsbrt.append('Lens Host %d' % e)
            gdat.numblablsbrt += 1
    if gmod.numbpopl > 0:
        if 'clus' in gmod.typeelem or 'clusvari' in gmod.typeelem:
            listlablsbrt.append('Uniform')
            gdat.numblablsbrt += 1
    
    listlablsbrtspec = ['Data']
    listlablsbrtspec += deepcopy(listlablsbrt)
    if len(listlablsbrt) > 1:
        listlablsbrtspec.append('Total Model')
    
    gdat.numblablsbrtspec = len(listlablsbrtspec)
    
    # number of generative parameters per element, depends on population
    #numbparaelem = gmod.numbparagenrelempopl + numbparaelemderi

    # maximum total number of parameters
    #numbparagenrfull = gmod.numbparagenrbase + gmod.numbpara.totl.elem
    
    #numbparaelemkind = gmod.numbparagenrbase
    #for l in gmod.indxpopl:
    #    numbparaelemkind += gmod.numbparagenrelemsing[l]
    
    #nameparagenrbase
    #gmod.namepara.genr.elem
    
    #listnameparaderifixd
    #listnameparaderielem
    
    #gmod.namepara.genrelemextd = gmod.namepara.genr.elem * maxm.numbelem
    #listnameparaderielemextd = gmod.namepara.genr.elem * maxm.numbelem
    
    #
    ## stack
    ## gmod.listnameparastck
    #gmod.listnameparastck = np.zeros(gmod.maxmnumbpara, dtype=object)
    #gmod.listscalparastck = np.zeros(gmod.maxmnumbpara, dtype=object)
    #
    #gmod.listnameparastck[gmod.indxpara.genrbase] = gmod.namepara.genr.base
    #gmod.listscalparastck[gmod.indxpara.genrbase] = gmod.listscalparagenrbase
    #for k in range(gmod.numbpara.totl.elem):
    #    for l in gmod.indxpopl:  
    #        if k >= gmod.numbparagenrelemcuml[l]:
    #            indxpopltemp = l
    #            indxelemtemp = (k - gmod.numbparagenrelemcuml[indxpopltemp]) // gmod.numbparagenrelemsing[indxpopltemp]
    #            gmod.indxpara.genrelemtemp = (k - gmod.numbparagenrelemcuml[indxpopltemp]) % gmod.numbparagenrelemsing[indxpopltemp]
    #            break
    #    gmod.listnameparastck[gmod.numbparagenrbase+k] = '%spop%d%04d' % (gmod.namepara.genrelem[indxpopltemp][gmod.indxpara.genrelemtemp], indxpopltemp, indxelemtemp)
    #    gmod.listscalparastck[gmod.numbparagenrbase+k] = gmod.scalpara.genrelem[indxpopltemp][gmod.indxpara.genrelemtemp]
    #
    #
    #if np.where(gmod.listscalpara == 0)[0].size > 0:
    #    print('gmod.listscalpara[gmod.indxpara.genrbase]')
    #    print(gmod.listscalpara[gmod.indxpara.genrbase])
    #    raise Exception('')
    #
    ## labels and scales for variables
    if gmod.boollens:
        for e in gmod.indxsersfgrd:
            # scalar
            setp_varb(gdat, 'masshostisf%dbein' % e, labl=[r'$M_{\rm{hst,%d,C}}$' % e, ''])
            
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
    
    gdat.lablfactpriodoff = r'$\alpha_{p}$'
    gmod.scalfactpriodoff = 'self'

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
    gdat.lablxposunit = gdat.lablgangunit
    gdat.lablyposunit = gdat.lablgangunit
   
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
            setp_varb(gdat, 'fdispop%dpop%d' % (l, q), minm=0., maxm=1., labl=['$F_{%d%d}$' % (l, q), ''], strgmodl=strgmodl)
            setp_varb(gdat, 'cmplpop%dpop%d' % (l, q), minm=0., maxm=1., labl=['$C_{%d%d}$' % (l, q), ''], strgmodl=strgmodl)
                    
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
    if gdat.typeexpr.startswith('HST_WFC3'):
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

    if gdat.boolbindspat:
        gdat.minmbein = 0.
        gdat.maxmbein = 1. / gdat.anglfact
    
    # scalar variables
    if gdat.boolbindspat:
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
    gmod.namepara.genrbase = []
    for name, k in gmod.indxpara.__dict__.items():
        
        #print('name')
        #print(name)
        #print('k')
        #print(k)
        
        if np.isscalar(k):
            gmod.namepara.genrbase.append(name)
    
    gmod.numbparagenrbase = len(gmod.namepara.genrbase)
    gmod.indxpara.genrbase = np.arange(gmod.numbparagenrbase)
    
    if gdat.booldiag:
        for name in gmod.namepara.genrbase:
            if 'pop0pop0' in name:
                raise Exception('')

        if gmod.numbparagenrbase == 0:
            print('gmod.namepara.genrbase')
            print(gmod.namepara.genrbase)
            raise Exception('')

        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')

    # base parameters to be perturbed
    gmod.indxpara.genrbasepert = gmod.indxpara.genrbase[gmod.numbpopl:]
    ## list of scalar variable names
    gmod.namepara.scal = list(gmod.namepara.genrbase) 
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
            gmod.namepara.derielemflat.append(gmod.namepara.derielemodim[l][k] + 'pop%d' % l)
            for d in range(gmod.maxmpara.numbelem[l]):
                gmod.namepara.derielemextd[l][k].append(gmod.namepara.derielemodim[l][k] + 'pop%d' % l + '%04d' % d)
                gmod.namepara.derielemextdflat.append(gmod.namepara.derielemextd[l][k][d])

    if gdat.booldiag:
        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')

    # list of element parameter names (derived and generative), counting label-degenerate element parameters only once 
    # to be deleted
    #gmod.namepara.elem = [[] for l in gmod.indxpopl]
    #for l in gmod.indxpopl:
    #    gmod.namepara.genr.elem[l].extend(gmod.namepara.genrelem[l])
    #    gmod.namepara.deri.elem[l].extend(gmod.namepara.derielemodim[l])
    
    gmod.namepara.elemflat = []
    for l in gmod.indxpopl:
        gmod.namepara.elemflat.extend(gmod.namepara.elem[l])

    gmod.namepara.genrelemdefa = deepcopy(gmod.namepara.elemflat)
    if gdat.boolbinsener and gmod.boolelemlghtanyy:
        for strgfeat in ['sind', 'curv', 'expc'] + ['sindcolr%04d' % i for i in gdat.indxenerinde]:
            if not strgfeat in gmod.namepara.genrelemdefa:
                gmod.namepara.genrelemdefa.append(strgfeat)

    if gdat.booldiag:
        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')

    # list of flattened generative element parameter names, counting label-degenerate element parameters only once
    gmod.namepara.genrelemkind = gmod.namepara.genrelemflat + gmod.namepara.derielemflat
    gmod.numbparagenrelemkind = len(gmod.namepara.genrelemkind)
    
    gmod.numbparagenrelemextdflat = len(gmod.namepara.genrelemextdflat)
    gmod.indxpara.genrelemextdflat = np.arange(gmod.numbparagenrelemextdflat)
    
    gmod.namepara.deribase = ['numbelemtotl']

    gmod.namepara.base = gmod.namepara.genrbase + gmod.namepara.deribase

    # list of parameter names (derived and generative), counting label-degenerate element parameters only once, element lists flattened
    gmod.namepara.kind = gmod.namepara.base + gmod.namepara.genrelemflat + gmod.namepara.derielemflat
    
    if gdat.booldiag:
        if np.unique(np.array(gmod.namepara.kind)).size != len(gmod.namepara.kind):
            print('')
            print('')
            print('')
            print('gmod.namepara.kind')
            print(gmod.namepara.kind)
            print('gmod.namepara.base')
            print(gmod.namepara.base)
            print('gmod.namepara.genrelemflat')
            print(gmod.namepara.genrelemflat)
            print('gmod.namepara.derielemflat')
            print(gmod.namepara.derielemflat)
            raise Exception('np.unique(np.array(gmod.namepara.kind)).size != len(gmod.namepara.kind)')

    gmod.numbparakind = len(gmod.namepara.kind)
    gmod.indxparakind = np.arange(gmod.numbparakind)

    # list of generative parameter names, separately including all label-degenerate element parameters, element lists flattened
    gmod.namepara.genrscalfull = gmod.namepara.genrbase + gmod.namepara.genrelemextdflat
    gmod.namepara.genrscalfull = np.array(gmod.namepara.genrscalfull)
    gmod.numbparagenr = len(gmod.namepara.genrscalfull)
    gmod.indxpara.genrfull = np.arange(gmod.numbparagenr)

    # list of generative parameter names, counting label-degenerate element parameters only once, element lists flattened
    gmod.listnameparagenrscal = gmod.namepara.genrbase + gmod.namepara.genrelemflat
    
    # list of parameter names (derived and generative), element lists flattened
    gmod.namepara.para = gmod.namepara.base + gmod.namepara.genrelemextdflat + gmod.namepara.derielemextdflat
    
    # to be deleted
    #for e in gmod.indxsersfgrd:
    #    strgsersfgrd = 'isf%d' % e
    #    gmod.namepara.scal += ['masshost' + strgsersfgrd + 'bein']
    #    for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
    #        gmod.namepara.scal += ['masshost' + strgsersfgrd + strgcalcmasssubh + 'bein']
    #if gmod.numbpopl > 0:
    #    if gmod.boollenssubh:
    #        for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
    #            gmod.namepara.scal += ['masssubh' + strgcalcmasssubh + 'bein', 'fracsubh' + strgcalcmasssubh + 'bein'] 
    
    if gmod.numbpopl > 0:
        gmod.namepara.scal += ['lpripena']
    if False and gmod.boolelemsbrtdfncanyy:
        for strgbins in ['lowr', 'higr']:
            gmod.namepara.scal += ['histcntp%sdfncen00evt0' % strgbins]
            gmod.namepara.scal += ['histcntp%sdfncsubten00evt0' % strgbins]
        for i in gdat.indxener:
            gmod.namepara.scal += ['fracsdenmeandarkdfncsubten%02d' % i]
        gmod.namepara.scal += ['booldfncsubt']
    if gmod.numbpopl > 0:
        for q in gdat.indxrefr:
            if gdat.boolasscrefr[q]:
                for l in gmod.indxpopl:
                    gmod.namepara.scal += ['cmplpop%dpop%d' % (l, q)]
                    gmod.namepara.scal += ['fdispop%dpop%d' % (q, l)]
    
    gmod.numbvarbscal = len(gmod.namepara.scal)
    gmod.indxvarbscal = np.arange(gmod.numbvarbscal)
    
    if gdat.booldiag:
        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')

    # determine total label
    gmod.namepara.glob = gmod.namepara.kind + gmod.namepara.genrelemextdflat + gmod.namepara.derielemextdflat
    gmod.namepara.glob += ['cntpmodl']
    for l in gmod.indxpopl:
        for g in gmod.indxparagenrelemsing[l]:
            if not gmod.namepara.genrelem[l][g] in gmod.namepara.glob:
                gmod.namepara.glob.append(gmod.namepara.genrelem[l][g])
                gmod.namepara.glob.append(gmod.namepara.derielemodim[l][g])
    
    if gdat.typeverb > 1:
        print('gmod.namepara.base')
        print(gmod.namepara.base)
        print('gmod.namepara.kind')
        print(gmod.namepara.kind)
        print('gmod.namepara.genrelemextdflat')
        print(gmod.namepara.genrelemextdflat)
        print('gmod.namepara.derielemextdflat')
        print(gmod.namepara.derielemextdflat)
        print('gmod.namepara.glob')
        print(gmod.namepara.glob)
        print('gmod.namepara.genrelem')
        print(gmod.namepara.genrelem)
        print('gmod.namepara.derielemodim')
        print(gmod.namepara.derielemodim)
        print('gmod.indxparagenrelemsing')
        print(gmod.indxparagenrelemsing)
        print('gmod.namepara.glob')
        print(gmod.namepara.glob)
    
    for name in gmod.namepara.glob:
        
        # set default scaling to 'self'
        if not hasattr(gmod.scalpara, name):
            setattr(gmod.scalpara, name, 'self')
        
        # set default number of bins to 10
        if not hasattr(gmod.numbbinspara, name):
            setattr(gmod.numbbinspara, name, 10)
        
    # define fact
    # to be deleted
    #for l in gmod.indxpopl:
    #    for k in gmod.indxparakind:
    #        name = gmod.namepara.kind[k]
    #        scal = getattr(gmod.scalpara, name)
    #        if scal == 'self' or scal == 'logt':
    #            minm = getattr(gmod.minmpara, name)
    #            maxm = getattr(gmod.maxmpara, name)
    #            if scal == 'self':
    #                fact = maxm - minm
    #            if scal == 'logt':
    #                fact = np.log(maxm / minm)
    #            
    #            if fact == 0:
    #                print('name')
    #                print(name)
    #                raise Exception('')
    #            setattr(gmod.factpara, name, fact)

    if gmod.numbpopl > 0:
        gmod.indxpara.genrelem = gmod.numbparagenrbase + np.arange(gmod.numbparagenrelem)
    
    print('gmod.namepara.genrelem')
    print(gmod.namepara.genrelem)
    gmod.namepara.genr = np.concatenate((gmod.namepara.genrbase, gmod.namepara.genrelemextdflat))
    
    # array of indices of all (base and transdimensional) generative parameters
    gmod.indxpara.genr = np.concatenate((gmod.indxpara.genrbase, gmod.indxpara.genrelem))

    if gdat.booldiag:
        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')
    
    ## arrays of parameter features (e.g., gmod.minmpara, maxm, labl, scal, etc.)
    print('Constructing arrays of parameter features such as minima, maxima, labels, scalings, etc...')
    for strgfeatpara in gdat.liststrgfeatparalist:
        
        if strgfeatpara == 'name':
            continue

        gmodtypefeat = getattr(gmod, strgfeatpara + 'para')
                
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
        #            listfeat[l][k] = getattr(gmodtypefeat, listname[l][k] + strgextn)
        #            listfeatflat.append(listfeat[l][k])
        #    setattr(gmodtypefeat, strgtypepara + 'elem', listfeat)
        #    setattr(gmodtypefeat, strgtypepara + 'elemflat', listfeatflat)
        
        if gdat.booldiag:
            if gmod.numbparagenrbase == 0:
                raise Exception('')
    
        # parameter subgroups (sgrp)
        # gdat.liststrgsgrppara is whatever follows gdat.liststrgfeatparalist after the dot (e.g., genrbase, genr)
        gdat.liststrgsgrppara = list(gmod.namepara.__dict__.keys())

        for strgsgrppara in gdat.liststrgsgrppara:
            
            if not strgsgrppara in ['genrbase', 'genr', 'numbelem', 'deribase', 'kind']:
                continue

            listtemp = getattr(gmod.namepara, strgsgrppara)
                
            if isinstance(listtemp[0], list):
                numbiter = len(listtemp)
                listname = listtemp
            else:
                numbiter = 1
                listname = [listtemp]

            for k in range(numbiter):
                
                feat = [0. for name in listname[k]]
                for n, name in enumerate(listname[k]):
                    
                    if hasattr(gmodtypefeat, name):
                        feat[n] = getattr(gmodtypefeat, name)
                    
                feat = np.array(feat)

                setattr(gmodtypefeat, strgsgrppara, feat)
    
    # dictionaries
    gmod.indxpara.kindscal = {}
    gmod.indxpara.genrbasescal = {}
    for scaltype in gdat.listscaltype:
        gmod.indxpara.kindscal[scaltype] = np.where(scaltype == gmod.scalpara.kind)[0]
        gmod.indxpara.genrbasescal[scaltype] = np.where(scaltype == gmod.scalpara.genrbase)[0]
    
    for strgfeatpara in gdat.liststrgfeatparalist:
        gmodtypefeat = getattr(gmod, strgfeatpara + 'para')
        if not isinstance(gmodtypefeat.genrbase, np.ndarray):
            print('TURNING gmodtypefeat.genrbase INTO NUMPY ARRAY. THIS SHOULD HAVE BEEN DONE BEFORE.')
            gmodtypefeat.genrbase = np.array(gmodtypefeat.genrbase)

        dicttemp = dict()
        for scaltype in gdat.listscaltype:
            dicttemp[scaltype] = gmodtypefeat.genrbase[gmod.indxpara.genrbasescal[scaltype]]
        setattr(gmodtypefeat, 'genrbasescal', dicttemp)

    if gdat.booldiag:
        if len(gmod.namepara.genrbase) == 0:
            raise Exception('')

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
    gmod.indxpara.genrbasescal = dict()
    for scaltype in gdat.listscaltype:
        gmod.indxpara.genrbasescal[scaltype] = np.where(np.array(gmod.scalpara.genrbase) == scaltype)[0]
    
    if gdat.booldiag:
        if np.where(gmod.scalpara.genrbase == 0)[0].size > 0:
            print('')
            print('')
            print('')
            print('gmod.scalpara.genrbase')
            print(gmod.scalpara.genrbase)
            raise Exception('np.where(gmod.scalpara.genrbase == 0)[0].size > 0')
    

def retr_liststrgcnfgprev(strgcnfg, pathbase):
    
    # list of PCAT run plot outputs
    pathvisu = pathbase + 'visuals/'
    liststrgcnfg = fnmatch.filter(os.listdir(pathvisu), '2*')
    
    liststrgcnfgprev = []
    for strgcnfg in liststrgcnfg:
        strgstat = pathbase + '/data/outp/' + strgcnfg
        
        if chec_statfile(pathbase, strgcnfg, 'gdatmodipost', typeverb=0) and strgcnfg + '_' + strgcnfg[16:].split('_')[-1] == strgcnfg[16:]:
            liststrgcnfgprev.append(strgcnfg) 
    
    liststrgcnfgprev.sort()

    return liststrgcnfgprev


def make_legd(axis, offs=None, loca=1, numbcols=1, ptch=None, line=None):
   
    hand, labl = axis.get_legend_handles_labels()
    legd = axis.legend(hand, labl, fancybox=True, frameon=True, bbox_to_anchor=offs, bbox_transform=axis.transAxes, ncol=numbcols, loc=loca, labelspacing=1, handlelength=2)
    legd.get_frame().set_fill(True)
    legd.get_frame().set_facecolor('white')


def setp_namevarbsing(gdat, gmod, strgmodl, strgvarb, popl, ener, dqlt, back, isfr, iele):
    
    if popl == 'full':
        indxpopltemp = gmod.indxpopl
    elif popl != 'none':
        indxpopltemp = [popl]
    
    if ener == 'full':
        indxenertemp = gdat.indxener
    elif ener != 'none':
        indxenertemp = [ener]
    
    if dqlt == 'full':
        indxdqlttemp = gdat.indxdqlt
    elif dqlt != 'none':
        indxdqlttemp = [dqlt]
    
    if isfr == 'full':
        indxisfrtemp = gmod.indxisfr
    elif isfr != 'none':
        indxisfrtemp = [isfr]
    
    if back == 'full':
        gmod.indxbacktemp = gmod.indxback
    elif isinstance(back, int):
        gmod.indxbacktemp = np.array([back])
    
    liststrgvarb = []
    if iele != 'none':
        for l in gmod.indxpopl:
            if iele == 'full':
                listiele = np.arange(gmod.maxmpara.numbelemtotl)
            else:
                listiele = [iele]
            for k in listiele:
                liststrgvarb.append(strgvarb + 'pop%d%04d' % (l, k))
    
    if popl != 'none' and ener == 'none' and dqlt == 'none' and back == 'none' and iele == 'none':
        for l in indxpopltemp:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    
    if popl == 'none' and ener == 'none' and dqlt == 'none' and back == 'none' and isfr != 'none':
        for e in indxisfrtemp:
            liststrgvarb.append(strgvarb + 'isf%d' % e)
    
    if popl == 'none' and ener != 'none' and dqlt != 'none' and back == 'none':
        for i in indxenertemp:
            for m in indxdqlttemp:
                liststrgvarb.append(strgvarb + 'en%02devt%d' % (i, m))
    
    if popl == 'none' and ener != 'none' and dqlt == 'none' and back != 'none':
        for c in gmod.indxbacktemp:
            for i in indxenertemp:
                liststrgvarb.append(strgvarb + 'back%04den%02d' % (c, i))
    
    if popl == 'none' and ener == 'none' and dqlt == 'none' and back != 'none':
        for c in gmod.indxbacktemp:
            liststrgvarb.append(strgvarb + 'back%04d' % c)
    
    if popl == 'none' and ener != 'none' and dqlt == 'none' and back == 'none':
        for i in indxenertemp:
            liststrgvarb.append(strgvarb + 'en%02d' % i)
    
    if popl == 'none' and ener == 'none' and dqlt == 'none' and back == 'none' and isfr == 'none':
        liststrgvarb.append(strgvarb)
    
    if gdat.booldiag:
        for strgvarb in liststrgvarb:
            if liststrgvarb.count(strgvarb) != 1:
                print('')
                print('')
                print('')
                print('liststrgvarb')
                print(liststrgvarb)
                print('popl')
                print(popl)
                print('ener')
                print(ener)
                print('dqlt')
                print(dqlt)
                print('back')
                print(back)
                print('isfr')
                print(isfr)
                print('iele')
                print(iele)
                raise Exception('liststrgvarb.count(strgvarb) != 1')
    
    return liststrgvarb


def setp_varb(gdat, \
              strgvarbbase, \
              
              valu=None, \
              
              minm=None, \
              maxm=None, \
              limt=None, \
              
              mean=None, \
              stdv=None, \
              
              scal=None, \
              
              bctr=None, \
              blim=None, \
              numbbins=None, \
              
              labl=None, \
              cmap=None, \
              
              popl='none', \
              ener='none', \
              dqlt='none', \
              back='none', \
              isfr='none', \
              iele='none', \
              
              boolinvr=False, \
              strgmodl=None, \
              strgstat=None, \
             ):
    '''
    Set up variable values across all models (true and fitting) as well as all populations, energy bins, 
    event bins, background components, and Sersic components 
    '''
    
    if limt is not None and (minm is not None or maxm is not None): 
        print('')
        print('')
        print('')
        raise Exception('limt is not None and (minm is not None or maxm is not None)')
    
    if limt is not None:
        minm = limt[0]
        maxm = limt[1]
    
    if valu is not None:
        if minm is not None:
            if minm > valu:
                print('')
                print('')
                print('')
                raise Exception('minm > valu')
        if maxm is not None:
            if maxm < valu:
                print('')
                print('')
                print('')
                raise Exception('maxm < valu')

    # determine the list of models
    if strgmodl is None:
        if gdat.typedata == 'simu':
            liststrgmodl = ['true', 'fitt', 'plot']
        else:
            liststrgmodl = ['fitt', 'plot']
    else:
        if strgmodl == 'true' or strgmodl == 'plot' or strgmodl == 'refr':
            liststrgmodl = [strgmodl]
        else:
            liststrgmodl = ['fitt', 'plot']
    
    for strgmodl in liststrgmodl:
        
        if strgmodl == 'plot':
            gmod = gdat.fitt
            gmodoutp = gdat
        else:
            gmod = getattr(gdat, strgmodl)
            gmodoutp = gmod

        # get the list of names of the variable
        liststrgvarbnone = setp_namevarbsing(gdat, gmod, strgmodl, strgvarbbase, popl, ener, dqlt, back, isfr, 'none')
        
        if iele != 'none':
            liststrgvarb = setp_namevarbsing(gdat, gmod, strgmodl, strgvarbbase, popl, ener, dqlt, back, isfr, iele)
        else:
            liststrgvarb = liststrgvarbnone
        
        # set the values of each variable in the list
        for strgvarb in liststrgvarb:
            if minm is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, minm, 'minmpara')
            
            if maxm is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, maxm, 'maxmpara')
            
            if mean is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, mean, 'meanpara')
            
            if stdv is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, stdv, 'stdvpara')
            
            if bctr is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, bctr, 'bctrpara')
            
            if blim is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, blim, 'blimpara')
            
            if cmap is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, cmap, 'cmappara')
            
            if valu is not None:
                if strgstat is None:
                    setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, valu, '')
                elif strgstat == 'this':
                    setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, valu, 'this')
                else:
                    raise Exception('')

            if scal is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, scal, 'scalpara')
                
            if labl is not None:
                
                if len(labl) != 2 or not isinstance(labl[0], str) and isinstance(labl[1], str):
                    print('')
                    print('')
                    print('')
                    raise Exception('labl input to setp_varb() has an issue.')

                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, labl[0], 'lablrootpara')
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, labl[1], 'lablunitpara')
                labltotl = tdpy.retr_labltotlsing(labl[0], labl[1])
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, labltotl, 'labltotlpara')
            
            if numbbins is not None:
                setp_varbcore(gdat, strgmodl, gmodoutp, strgvarb, numbbins, 'numbbinspara')
                
            # create limtpara, binspara, indxpara, meanpara, deltpara, etc.
            if minm is not None and maxm is not None or mean is not None and stdv is not None:
                
                if numbbins is None:
                    if isinstance(minm, int) and isinstance(maxm, int):
                        numbbins = maxm - minm + 1
                    else:
                        numbbins = 10

                if gdat.booldiag:
                    if numbbins is None:
                        print('')
                        print('')
                        print('')
                        raise Exception('numbbins is None')

                # determine minima and maxima for Gaussian or log-Gaussian distributed parameters
                if mean is not None:
                    minm = mean - gdat.numbstdvgaus * stdv
                    maxm = mean + gdat.numbstdvgaus * stdv
                
                # uniformly-distributed
                if scal == 'self' or scal == 'drct' or scal == 'pois' or scal == 'gaus' or scal is None:
                    binsunif = np.linspace(minm, maxm, numbbins + 1)
                # 'temp' fix below, every scal should be different
                elif scal == 'logt' or scal == 'powr':
                    binsunif = np.linspace(np.log10(minm), np.log10(maxm), numbbins + 1)
                    if gdat.booldiag:
                        if minm <= 0.:
                            print('')
                            print('')
                            print('')
                            print('minm')
                            print(minm)
                            print('scal')
                            print(scal)
                            raise Exception('minm <= 0 but scal == logt or scal == powr.')
                elif scal == 'asnh':
                    binsunif = np.linspace(np.arcsinh(minm), np.arcsinh(maxm), numbbins + 1)
                else:
                    print('scal')
                    print(scal)
                    raise Exception('')
                
                if boolinvr:
                    binsunif = binsunif[::-1]
                
                bctrunif = (binsunif[1:] + binsunif[:-1]) / 2.
                
                if scal == 'self' or scal == 'drct' or scal == 'pois' or scal == 'gaus' or scal is None:
                    bctr = bctrunif
                    blim = binsunif
                    minmunif = minm
                    maxmunif = maxm
                if scal == 'logt' or scal == 'powr':
                    bctr = 10**bctrunif
                    blim = 10**binsunif
                    minmunif = np.log10(minm)
                    maxmunif = np.log10(maxm)
                if scal == 'asnh':
                    bctr = np.sinh(bctrunif)
                    blim = np.sinh(binsunif)
                    minmunif = np.arcsinh(minm)
                    maxmunif = np.arcsinh(maxm)
                
                indxbins = np.arange(numbbins)

                delt = np.diff(blim) 
                limt = np.array([minm, maxm]) 
                
                # 'self' is not yet defined
                if (scal == 'asnh' or scal == 'logt' or scal == 'powr'):
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

                setattr(gmodoutp.indxbinspara, strgvarb, indxbins)
                setattr(gmodoutp.numbbinspara, strgvarb, numbbins)
                setattr(gmodoutp.limtpara, strgvarb, limt)
                setattr(gmodoutp.blimpara, strgvarb, blim)
                setattr(gmodoutp.bctrpara, strgvarb, bctr)
                setattr(gmodoutp.deltpara, strgvarb, delt)
       

def retr_ticklabltemp_TOBEDELETED(gdat, strgcbar):
    
    minm = getattr(gdat.minmpara, strgcbar)
    maxm = getattr(gdat.maxmpara, strgcbar)
    scal = getattr(gdat.scalpara, strgcbar)
    numb = gdat.numbtickcbar - 1
    setp_varb(gdat, strgcbar, numb=numb)

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



def setp_varb_TOBEDELETED(gdat, strgvarb, strgmodl=None, boolinvr=False):
    
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
        if not hasattr(gdattemp.minmpara, strgvarb) or \
           not hasattr(gdattemp.maxmpara, strgvarb) or \
           not hasattr(gdattemp.numbbinspara, strgvarb) or \
           not hasattr(gdattemp.scalpara, strgvarb):
            continue
            
        minm = getattr(gdattemp.minmpara, strgvarb)
        maxm = getattr(gdattemp.maxmpara, strgvarb)
        numb = getattr(gdattemp.numbbinspara, strgvarb)
        scal = getattr(gdattemp.scalpara, strgvarb)

        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            binsscal = np.linspace(minm, maxm, numb + 1)
        # 'temp' fix 'powr', it's wrong
        elif scal == 'logt' or scal == 'powr':
            binsscal = np.linspace(np.log10(minm), np.log10(maxm), numb + 1)
            if gdat.booldiag:
                if minm <= 0.:
                    raise Exception('')

        elif scal == 'asnh':
            binsscal = np.linspace(np.arcsinh(minm), np.arcsinh(maxm), numb + 1)
        else:
            print('')
            print('scal')
            print(scal)
            raise Exception('')

        if boolinvr:
            binsscal = binsscal[::-1]
        
        meanvarbscal = (binsscal[1:] + binsscal[:-1]) / 2.
        
        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            meanvarb = meanvarbscal
            blim = binsscal
        # 'temp' fix 'powr', it's wrong
        if scal == 'logt' or scal == 'powr':
            meanvarb = 10**meanvarbscal
            blim = 10**binsscal
        if scal == 'asnh':
            meanvarb = np.sinh(meanvarbscal)
            blim = np.sinh(binsscal)

        delt = np.diff(bins) 
        limt = np.array([np.amin(bins), np.amax(bins)]) 
        
        setattr(gdattemp.limtpara, strgvarb, limt)
        setattr(gdattemp.blimpara, strgvarb, bins)
        setattr(gdattemp.bctrpara, strgvarb, meanvarb)
        setattr(gdattemp.deltpara, strgvarb, delt)


def setp_varbcore(gdat, strgmodl, gdattemptemp, strgvarbtemp, valu, strgtake):
    
    if strgtake == '':
        gdattemp = gdattemptemp
    else:
        gdattemp = getattr(gdattemptemp, strgtake)

    if hasattr(gdattemp, strgvarbtemp):
        valutemp = getattr(gdattemp, strgvarbtemp)
        if gdat.typeverb > 0:
            print('Received custom value for variable %s, model %s, and feature %s: %s' % (strgvarbtemp, strgmodl, strgtake, valutemp))
    else:
        if strgvarbtemp.startswith('beinh') or strgvarbtemp.startswith('cntpmod'):
            print('Setting default value for variable %s, model %s, and feature %s: %s' % (strgvarbtemp, strgmodl, strgtake, valu))
        setattr(gdattemp, strgvarbtemp, valu)


def intp_sinc(gdat, xpos, ypos):

    intpsinc = 4. * gdat.numbsidepsfn**2 * np.sum(gdat.temppsfn * sinc(gdat.numbsidepsfn * (gdat.gridpsfnxpos + xpos) - gdat.gridpsfnxpos) * \
                                                                                      sinc(gdat.numbsidepsfn * (gdat.gridpsfnypos + ypos) - gdat.gridpsfnypos))

    return intpsinc


def retr_fluxbrgt(gdat, xpos, ypos, flux):

    if xpos.size == 0:
        fluxbrgt = np.array([0.])
        fluxbrgtassc = np.array([0.])
    else:
        indxbrgt = np.argmax(flux)
        fluxbrgt = flux[indxbrgt]

    return fluxbrgt, fluxbrgtassc


def init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxdqltplot, indxpoplplot):
    
    figrsize = (gdat.plotsize, gdat.plotsize)
    figr, axis = plt.subplots(figsize=figrsize)
    
    nameplot = strgplot

    if gdat.numbener > 1:
        nameplot += 'en%02d' % gdat.indxenerincl[indxenerplot]
    
    if gdat.numbener > 1:
        if indxdqltplot == -1:
            nameplot += 'evtA'
        else:
            nameplot += 'evt%d' % gdat.indxdqltincl[indxdqltplot]
    
    if gdat.fitt.numbpopl > 1:
        if indxpoplplot == -1:
            nameplot += 'popA'
        else:
            nameplot += 'pop%d' % indxpoplplot

    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, nameplot)
    
    axis.set_xlabel(gdat.fitt.labltotlpara.xpospop0)
    axis.set_ylabel(gdat.fitt.labltotlpara.ypospop0)
    titl = ''
    if indxenerplot is not None and gdat.numbener > 1 and strgplot.endswith('cnts'):
        titl = gdat.strgener[indxenerplot]
    if indxdqltplot is not None and gdat.numbdqlt > 1 and strgplot.endswith('cnts'):
        titl += ' ' + gdat.strgdqlt[indxdqltplot]
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


def retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot=None, indxdqltplot=-1, booltdim=False, imag=None):
    
    draw_frambndr(gdat, axis)
    
    # take the relevant energy and PSF bins
    if indxenerplot is not None:
        if indxdqltplot == -1:
            maps = np.sum(maps[indxenerplot, ...], axis=1)
        else:
            maps = maps[indxenerplot, :, indxdqltplot]
    
    # project the map to 2D
    if gdat.typepixl == 'heal':
        maps = tdpy.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                    minmxpos=gdat.anglfact*gdat.minmxposdata, maxmxpos=gdat.anglfact*gdat.maxmxposdata, \
                                                                    minmypos=gdat.anglfact*gdat.minmyposdata, maxmypos=gdat.anglfact*gdat.maxmyposdata)
    
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
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05, aspect=15)
    
    if hasattr(gdat.valutickmajrpara, strgvarb):
        valutickmajr = getattr(gdat.valutickmajrpara, strgvarb)
        labltickmajr = getattr(gdat.labltickmajrpara, strgvarb)
        cbar.set_ticks(valutickmajr)
        cbar.set_ticklabels(labltickmajr)
    
    #print('valutickmajr')
    #print(valutickmajr)
    #print('labltickmajr')
    #print(labltickmajr)
    #print('strgvarb')
    #print(strgvarb)
    #print('')
    
    return cbar


def make_legdmaps(gdat, strgstat, strgmodl, axis, mosa=False, assc=False):
    
    gmod = getattr(gdat, strgmodl)
    
    # transdimensional elements
    if strgmodl == 'fitt' and (strgstat == 'pdfn' and gdat.boolcondcatl or strgstat == 'this') and gmod.numbpopl > 0:
        for l in gmod.indxpopl:
            colr = retr_colr(gdat, strgstat, strgmodl, l)
            if strgstat == 'pdfn':
                labl = 'Condensed %s' % (gmod.legdpopl[l])
            else:
                labl = 'Sample %s' % (gmod.legdpopl[l])
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
                                                           label='%s Lens Host' % gmod.lablmodl, marker='s', lw=gdat.mrkrlinewdth, color=gmod.colr)
    
    if gdat.typedata == 'simu':
        if gmod.boollens:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Source' % gdat.refr.labl, marker='>', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
        
        if gmod.typeemishost != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Lens Host' % gdat.refr.labl, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
    
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
                xpos = np.copy(gdat.refr.dictelem[q]['xpos'][0, :])
                ypos = np.copy(gdat.refr.dictelem[q]['ypos'][0, :])
                numbelem = int(gdat.refr.numbelem[q])
                
                if gdatmodi is not None and gmod.numbpopl > 0 and assc:   
                    ### hit
                    indx = gdatmodi.this.indxelemrefrasschits[q][l]
                    if indx.size > 0:
                        axis.scatter(gdat.anglfact * xpos[indx], gdat.anglfact * ypos[indx], s=mrkrsize[indx], alpha=gdat.alphelem, label=gdat.refr.lablhits, \
                                                                      marker=gdat.refrlistmrkrhits[q], lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
                    ### missed
                    indx = gdatmodi.this.indxelemrefrasscmiss[q][l]
                else:
                    indx = np.arange(xpos.size)
                
                if indx.size > 0: 
                    axis.scatter(gdat.anglfact * xpos[indx], gdat.anglfact * ypos[indx], s=mrkrsize[indx], alpha=gdat.alphelem, facecolor='none', \
                                                             label=gdat.refr.listlablmiss, marker=gdat.refr.listmrkrmiss[q], \
                                                             lw=gdat.mrkrlinewdth, color=gdat.refr.colrelem[q])
        
            sizexoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            sizeyoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            if 'etag' in gdat.refr.namepara.elem[q]:
                for k in range(indx.size):
                    axis.text(gdat.anglfact * xpos[indx[k]] + sizexoff, gdat.anglfact * ypos[indx[k]] + sizeyoff, gdat.refretag[q][indx[k]], \
                                                                                            verticalalignment='center', horizontalalignment='center', \
                                                                                                             color='red', fontsize=1)

    # temp -- generalize this to input refrxposhost vs.
    if gdat.typedata == 'simu':
        ## host galaxy position
        if gmod.typeemishost != 'none':
            for e in gmod.indxsersfgrd:
                xposhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'xposhostisf%d' % (e))]
                yposhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'yposhostisf%d' % (e))]
                axis.scatter(gdat.anglfact * xposhost, gdat.anglfact * yposhost, facecolor='none', alpha=0.7, \
                                             label='%s Lens Host %d' % (gdat.refr.labl, e), s=300, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refr.colr)
        if gmod.boollens:
            ## host galaxy Einstein radius
            for e in gmod.indxsersfgrd:
                truexposhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'xposhostisf%d' % (e))]
                trueyposhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'yposhostisf%d' % (e))]
                truebeinhost = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % (e))]
                axis.add_patch(plt.Circle((gdat.anglfact * truexposhost, \
                                           gdat.anglfact * trueyposhost), \
                                           gdat.anglfact * truebeinhost, \
                                           edgecolor=gdat.refr.colr, facecolor='none', lw=gdat.mrkrlinewdth))
            
        if gmod.boollens:
            ## source galaxy position
            axis.scatter(gdat.anglfact * gmodstat.paragenrscalfull[gmod.indxpara.xpossour], \
                                                    gdat.anglfact * gmodstat.paragenrscalfull[gmod.indxpara.ypossour], \
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
            if gmod.numbpopl > 0:
                colr = retr_colr(gdat, strgstat, strgmodl, l)
                mrkrsize = retr_mrkrsize(gdat, strgmodl, gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[gmod.nameparagenrelemampl[l]][l]], gmod.nameparagenrelemampl[l])
                if 'xpos' in gdatmodi.this.indxparagenrelemfull:
                    xpos = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l]['xpos']]
                    ypos = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l]['ypos']]
                else:
                    gang = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l]['gang']]
                    aang = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l]['aang']]
                    xpos, ypos = retr_xposypos(gang, aang)
                axis.scatter(gdat.anglfact * xpos, gdat.anglfact * ypos, s=mrkrsize, alpha=gdat.alphelem, label='Sample', marker=gmod.listelemmrkr[l], \
                                                                                                                           lw=gdat.mrkrlinewdth, color=colr)

            ## source
            if gmod.boollens:
                xpossour = gdatmodi.this.paragenrscalfull[gmod.indxpara.xpossour]
                ypossour = gdatmodi.this.paragenrscalfull[gmod.indxpara.ypossour]
                axis.scatter(gdat.anglfact * xpossour, gdat.anglfact * ypossour, facecolor='none', \
                                                  alpha=gdat.alphelem, \
                                                  label='%s Source' % gmod.lablpara, s=300, marker='<', lw=gdat.mrkrlinewdth, color=gmod.colr)
    
            if gmod.typeemishost != 'none':
                ## host
                xposhost = [[] for e in gmod.indxsersfgrd]
                yposhost = [[] for e in gmod.indxsersfgrd]
                for e in gmod.indxsersfgrd:
                    xposhost[e] = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'xposhostisf%d' % (e))]
                    yposhost[e] = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'yposhostisf%d' % (e))]
                    axis.scatter(gdat.anglfact * xposhost[e], gdat.anglfact * yposhost[e], facecolor='none', \
                                                     alpha=gdat.alphelem, \
                                                     label='%s Lens Host' % gmod.lablpara, s=300, marker='s', lw=gdat.mrkrlinewdth, color=gmod.colr)
                    if gmod.boollens:
                        beinhost = gdatmodi.this.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % (e))]
                        axis.add_patch(plt.Circle((gdat.anglfact * xposhost[e], gdat.anglfact * yposhost[e]), \
                                                       gdat.anglfact * beinhost, edgecolor=gmod.colr, facecolor='none', \
                                                       lw=gdat.mrkrlinewdth, ls='--'))
                
    # temp
    if strgstat == 'pdfn' and gdat.boolcondcatl and gmod.numbpopl > 0:
        xpos = np.zeros(gdat.numbprvlhigh)
        ypos = np.zeros(gdat.numbprvlhigh)
        ampl = np.zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                xpos[cntr] = gdat.dictglob['poststkscond'][r]['xpos'][0]
                ypos[cntr] = gdat.dictglob['poststkscond'][r]['ypos'][0]
                # temp -- this does not allow sources with different spectra to be assigned to the same stacked sample
                ampl[cntr] = gdat.dictglob['poststkscond'][r][gmod.nameparagenrelemampl[l]][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, strgmodl, ampl, gmod.nameparagenrelemampl[l])
        
        colr = retr_colr(gdat, strgstat, strgmodl, l)
        axis.scatter(gdat.anglfact * xpos, gdat.anglfact * ypos, s=mrkrsize, \
                                    label='Condensed', marker=gmod.listelemmrkr[l], color='black', lw=gdat.mrkrlinewdth)
        for r in gdat.indxstkscond:
            xpos = np.array([gdat.dictglob['liststkscond'][r]['xpos']])
            ypos = np.array([gdat.dictglob['liststkscond'][r]['ypos']])
            axis.scatter(gdat.anglfact * xpos, gdat.anglfact * ypos, s=mrkrsize, \
                                                marker=gmod.listelemmrkr[l], color='black', alpha=0.1, lw=gdat.mrkrlinewdth)


def retr_colr(gdat, strgstat, strgmodl, indxpopl=None):
    
    if strgmodl == 'true':
        if indxpopl is None:
            colr = gdat.refr.colr
        else:
            colr = gdat.refr.colrelem[indxpopl]
    if strgmodl == 'fitt':
        if strgstat == 'this' or strgstat == 'pdfn':
            gmod = getattr(gdat, strgmodl)
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
    
    fluxpare, xpospare, ypospare, fluxauxi, xposauxi, yposauxi = sympy.symbols('fluxpare xpospare ypospare fluxauxi xposauxi yposauxi')
    
    matr = sympy.Matrix([[ fluxpare,      fluxauxi, 0,            0, 0,            0], \
                         [-fluxpare, 1  - fluxauxi, 0,            0, 0,            0], \
                         [-xposauxi,             0, 1, 1 - fluxauxi, 0,            0], \
                         [-xposauxi,             0, 1,    -fluxauxi, 0,            0], \
                         [-yposauxi,             0, 0,            0, 1, 1 - fluxauxi], \
                         [-yposauxi,             0, 0,            0, 1,    -fluxauxi]])

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

def retr_angldist(gdat, xposfrst, yposfrst, xposseco, yposseco):
    
    # temp -- heal does not work when the dimension of xposfrst is 1
    if gdat.typepixl == 'heal':
        dir1 = np.array([xposfrst, yposfrst])
        dir2 = np.array([xposseco, yposseco])
        angldist = hp.rotator.angdist(dir1, dir2)
    else:
        angldist = np.sqrt((xposfrst - xposseco)**2 + (yposfrst - yposseco)**2)

    return angldist


def retr_deflextr(gdat, indxpixlelem, sher, sang):
    
    factcosi = sher * np.cos(2. * sang)
    factsine = sher * np.cos(2. * sang)
    deflxpos = factcosi * gdat.xposgrid[indxpixlelem] + factsine * gdat.yposgrid[indxpixlelem]
    deflypos = factsine * gdat.xposgrid[indxpixlelem] - factcosi * gdat.yposgrid[indxpixlelem]
    
    return np.vstack((deflxpos, deflypos)).T 


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
    gmod.this.paragenrunitfull = np.random.rand(gmod.numbparagenr)
    gmod.this.paragenrscalfull = np.empty(gmod.numbparagenr)

    ## impose user-specified initial state
    ### number of elements
    ## create dummy indxparagenrfullelem 
    gmod.this.indxparagenrelemfull = None
    if gmod.numbpopl > 0:
        if gdat.inittype == 'refr':
            for l in gmod.indxpopl:
                gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]] = gmod.paragenrunitfull[gmod.indxpara.numbelem[l]]
        else:
            for l in gmod.indxpopl:
                if gmod.typemodltran == 'pois':
                    meanelemtemp = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, \
                                                gmod.this.indxparagenrelemfull)[gmod.indxpara.meanelem[l]]
                
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
    
    if gdat.booldiag:
        if gdat.typedata == 'simu' and gdat.inittype == 'refr':
            for l in gmod.indxpopl:
                if gmod.paragenrunitfull[gmod.indxpara.numbelem[l]] > gmod.maxmpara.numbelem[l]:
                    raise Exception('')

    if gmod.numbpopl > 0:
        gmod.this.indxelemfull = []
        for l in gmod.indxpopl:
            gmod.this.indxelemfull.append(list(range(gmod.this.paragenrunitfull[gmod.indxpara.numbelem[l]].astype(int))))
        gmod.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gmod.this.indxelemfull, 'fitt')

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
                        if attr.startswith('xpospop'):
                            gmod.indxpopl = int(attr[7])
                            if gmod.indxpopl > maxmindxpopl:
                                maxmindxpopl = gmod.indxpopl
                numbpoplinpt = maxmindxpopl + 1
                
                if numbpoplinpt != gmod.numbpopl:
                    print('State file and fitting metamodel have different number of populations.')
                
                # find the number of elements provided
                cntr = np.zeros(gmod.numbpoplinpt, dtype=int)
                for attr in thisfile:
                    if attr.startswith('xpospop'):
                        gmod.indxpopl = int(attr[7])
                        cntr[indxpopl] += 1
                if gdat.typeverb > 0:
                    print('Number of elements found:')
                    print(cntr)

            for attr in thisfile:
                for k, namepara in enumerate(gmod.namepara.genrbase):
                    if namepara == attr:
                        if gmod.namepara.genr.base.startswith('numbelem'):
                            try:
                                indxpopltemp = int(gmod.namepara.genrbase[-1])
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
                        if (gmod.numbparaelem == 0 or gmod.numbpopl > 0 and not k in gmod.indxpara.numbelem) and \
                                        (not np.isfinite(gmod.this.paragenrunitfull[k]) or gmod.this.paragenrunitfull[k] < 0. or \
                                        gmod.this.paragenrunitfull[k] > 1.):
                            raise Exception('CDF of the retreived state parameter is bad.')
            if gmod.numbpopl > 0:
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
            
            gmod.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gmod.this.indxelemfull, 'fitt')
            gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrelemfull)
            
            if (gmod.this.paragenrunitfull == 0).all():
                raise Exception('Bad initialization.')
    
            if gmod.numbpopl > 0 and gmod.this.indxparagenrelemfull is not None:
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
        if gdat.typedata == 'simu':
            for k, namepara in enumerate(gmod.namepara.genrbase):
                if not (gdat.inittype == 'pert' and gmod.namepara.genr.base.startswith('numbelem')):
                    gmod.indxpara.true = k
                    gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gmodstat.paragenrscalfull[gmod.indxpara.true], k)
        if gmod.numbpopl > 0:
            gmod.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gmod.this.indxelemfull, 'fitt')
        if gdat.typeverb > 1:
            show_paragenrscalfull(gdat, gdatmodi)
        if gmod.this.indxparagenrelemfull is not None:
            print('Initializing elements from the reference element parameters...')
            show_paragenrscalfull(gdat, gdatmodi)
            gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrelemfull)
            show_paragenrscalfull(gdat, gdatmodi)
            initcompfromstat(gdat, gdatmodi, 'refr')
        gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrelemfull)
    
    ## impose user-specified individual initial values
    for k, namepara in enumerate(gmod.namepara.genrbase):
        if namepara.startswith('numbelem'):
            continue
        if gdat.inittype == 'reco' or  gdat.inittype == 'refr' or gdat.inittype == 'pert':
            getattr(gdat, 'init' + gmod.namepara.genrbase)
            print('Conflicting initial state arguments detected, init keyword takes precedence.')
            
        initvalu = getattr(gdat, 'init' + gmod.namepara.genrbase)
        gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', initvalu, k)
        if gdat.typeverb > 0:
            print('Received initial condition for %s: %.3g' % (gmod.namepara.genr.base, initvalu))
    
    ## PSF
    if gdat.initpsfp is not None:
        print('Initializing the metamodel PSF from the provided initial state...')
        if gdat.initpsfp.size != gmod.indxpara.psfp.size:
            raise Exception('')
        for k, namepara in enumerate(gmod.namepara.genrbase):
            if k in gmod.indxpara.psfp:
                gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gdat.initpsfp[k-gmod.indxpara.psfp[0]], k)
    if gdat.initpsfprefr:
        print('Initializing the metamodel PSF from the reference state...')
        for k, namepara in enumerate(gmod.namepara.genrbase):
            if k in gmod.indxpara.psfp:
                gmod.this.paragenrunitfull[k] = cdfn_paragenrscalbase(gdat.fitt, '', gmod.psfpexpr[k-gmod.indxpara.psfp[0]], k)

    if gdat.inittype == 'rand' or gdat.inittype == 'reco' and not boolinitreco:
        if gdat.typeverb > 0:
            print('Initializing from a random state...')
        gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrelemfull)
    
    if gmod.numbpopl > 0:
        gmod.this.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gmod.this.indxelemfull, 'fitt')

    # check the initial unit sample vector for bad entries
    if gmod.numbpopl > 0:
        indxsampdiff = np.setdiff1d(gmod.indxpara.genrfull, gmod.indxpara.numbelem)
        
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
    
    gmod.this.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gmod.this.paragenrunitfull, gmod.this.indxparagenrelemfull)
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
   

def initchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi.this, 'chro' + name, gdat.functime())
    

def stopchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi.this, 'chro' + name, gdat.functime() - getattr(gdatmodi.this, 'chro' + name))


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
    lpdfspatprioobjt = sp.interpolate.RectBivariateSpline(gdat.blimpara.yposcart, gdat.blimpara.xposcart, lpdfspatprio)
    
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

    if gmod.typeevalpsfn != 'none' and (strgmodl == 'true' or boolinit or gdat.boolmodipsfn):
        psfp = gmodstat.paragenrscalfull[gmod.indxpara.psfp]
        if gdat.booldiag:
            if np.where(psfp == 0)[0].size == psfp.size:
                raise Exception('')
        setattr(gmodstat, 'psfp', psfp)
    bacp = gmodstat.paragenrscalfull[gmod.indxpara.bacp]
   
    if gmod.numbpopl > 0:
        
        # temp -- this may slow down execution
        gmodstat.indxparagenrelemfull = retr_indxparagenrelemfull(gdat, gmodstat.indxelemfull, strgmodl)
        
        # check if all active generative parameters are finite
        if gdat.booldiag:
            indxparatemp = []
            for l in gmod.indxpopl:
                indxparatemp.append(gmodstat.indxparagenrelemfull[l]['full'])
            indxparatemp.append(gmod.indxpara.genrbase)
            gmodstat.indxpara.genrfull = np.concatenate(indxparatemp)
            if not np.isfinite(gmodstat.paragenrscalfull[gmodstat.indxpara.genrfull]).all():
                raise Exception('')

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
                gmodstat.dictelem[l][nameparagenrelem] = gmodstat.paragenrscalfull[gmodstat.indxparagenrelemfull[l][nameparagenrelem]]
                if gdat.booldiag:
                    if ((abs(gmodstat.paragenrscalfull[gmodstat.indxparagenrelemfull[l][nameparagenrelem]]) < 1e-100 ) & \
                        (abs(gmodstat.paragenrscalfull[gmodstat.indxparagenrelemfull[l][nameparagenrelem]]) > 0.)).any():
                        raise Exception('')

                    if gmodstat.numbelem[l] != len(gmodstat.dictelem[l][nameparagenrelem]):
                        print('l')
                        print(l)
                        print('gmodstat.numbelem[l]')
                        print(gmodstat.numbelem[l])
                        print('gmodstat.dictelem[l]')
                        print(gmodstat.dictelem[l])
                        print('gmodstat.dictelem[l][nameparagenrelem]')
                        print(gmodstat.dictelem[l][nameparagenrelem])
                        print('nameparagenrelem')
                        print(nameparagenrelem)
                        raise Exception('')
    
        if gdat.boolbinsener:
            if gdat.typeverb > 2:
                print('Calculating element spectra...')
            initchro(gdat, gdatmodi, 'spec')
            for l in gmod.indxpopl:
                if gmod.typeelem[l].startswith('lght'):
                    for strgfeat in gmod.namepara.genrelem[l]:
                        sindcolr = [gmodstat.dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                        gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], \
                                                    sind=gmodstat.dictelem[l]['sind'], curv=gmodstat.dictelem[l]['curv'], \
                                                    expc=gmodstat.dictelem[l]['expc'], sindcolr=sindcolr, spectype=gmod.spectype[l])
                        if gmod.typeelem[l].startswith('lghtline'):
                            if gmod.typeelem[l] == 'lghtlinevoig':
                                gmodstat.dictelem[l]['spec'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], \
                                                                                elin=gmodstat.dictelem[l]['elin'], sigm=gmodstat.dictelem[l]['sigm'], \
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
    
        if gdat.booldiag:
            for l in gmod.indxpopl:
                for g, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                    if (gmod.scalpara.genrelem[l][g] != 'gaus' and not gmod.scalpara.genrelem[l][g].startswith('lnor')) and  \
                       (gmod.scalpara.genrelem[l][g] != 'expo' and (gmodstat.dictelem[l][nameparagenrelem] < getattr(gmod.minmpara, nameparagenrelem)).any()) or \
                                        (gmodstat.dictelem[l][nameparagenrelem] > getattr(gmod.maxmpara, nameparagenrelem)).any():
                        
                        print('')
                        print('')
                        print('')
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
                        print('gmod.scalpara.genrelem[l][g]')
                        print(gmod.scalpara.genrelem[l][g])
                        raise Exception('')
           
        # calculate element spectra
        # temp
        if gdat.booldiag:
            for l in gmod.indxpopl:
                if gmod.typeelem[l] == 'lens':
                    if gdat.variasca:
                        indx = np.where(gmodstat.paragenrscalfull[gmodstat.indxparagenrelemfull[l]['acut']] < 0.)[0]
                        if indx.size > 0:
                            raise Exception('')
                    if gdat.variacut:
                        indx = np.where(gmodstat.paragenrscalfull[gmodstat.indxparagenrelemfull[l]['asca']] < 0.)[0]
                        if indx.size > 0:
                            raise Exception('')
    
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lght'):
                    
                # evaluate horizontal and vertical position for elements whose position is a power law in image-centric radius
                if gmod.typespatdist[l] == 'glc3':
                    gmodstat.dictelem[l]['dlos'], gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos'] = retr_glc3(gmodstat.dictelem[l]['dglc'], \
                                                                                                        gmodstat.dictelem[l]['thet'], gmodstat.dictelem[l]['phii'])
                
                if gmod.typespatdist[l] == 'gangexpo':
                    gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos'], = retr_xposypos(gmodstat.dictelem[l]['gang'], \
                                                                                                        gmodstat.dictelem[l]['aang'])
                    
                    if gdat.booldiag:
                        if gmodstat.numbelem[l] > 0:
                            if np.amin(gmodstat.dictelem[l]['xpos']) < gmod.minmxpos or \
                               np.amax(gmodstat.dictelem[l]['xpos']) > gmod.maxmxpos or \
                               np.amin(gmodstat.dictelem[l]['ypos']) < gmod.minmypos or \
                               np.amax(gmodstat.dictelem[l]['ypos']) > gmod.maxmypos:
                                raise Exception('Bad coordinates!')

                if gmod.typespatdist[l] == 'los3':
                    gmodstat.dictelem[l]['dglc'], gmodstat.dictelem[l]['thet'], gmodstat.dictelem[l]['phii'] = retr_los3(gmodstat.dictelem[l]['dlos'], \
                                                                                                        gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos'])

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
        xpossour = gmodstat.paragenrscalfull[gmod.indxpara.xpossour]
        ypossour = gmodstat.paragenrscalfull[gmod.indxpara.ypossour]
    
    if gdat.typeverb > 2:
        print('Evaluating the likelihood...')
    
    # process a sample vector and the occupancy list to calculate secondary variables
    if gmod.boollens:
        gmodstat.fluxsour = gmodstat.paragenrscalfull[gmod.indxpara.fluxsour]
        if gdat.numbener > 1:
            gmodstat.sindsour = gmodstat.paragenrscalfull[gmod.indxpara.sindsour]
        gmodstat.sizesour = gmodstat.paragenrscalfull[gmod.indxpara.sizesour]
        gmodstat.ellpsour = gmodstat.paragenrscalfull[gmod.indxpara.ellpsour]
        gmodstat.anglsour = gmodstat.paragenrscalfull[gmod.indxpara.anglsour]
    if gmod.typeemishost != 'none':
        gmodstat.xposhost = [[] for e in gmod.indxsersfgrd]
        gmodstat.yposhost = [[] for e in gmod.indxsersfgrd]
        gmodstat.fluxhost = [[] for e in gmod.indxsersfgrd]
        if gdat.numbener > 1:
            gmodstat.sindhost = [[] for e in gmod.indxsersfgrd]
        gmodstat.sizehost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            gmodstat.xposhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'xposhostisf%d' % e)]
            gmodstat.yposhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'yposhostisf%d' % e)]
            gmodstat.fluxhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'fluxhostisf%d' % e)]
            if gdat.numbener > 1:
                gmodstat.sindhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sindhostisf%d' % e)]
            gmodstat.sizehost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'sizehostisf%d' % e)]
    if gmod.boollens:
        gmodstat.beinhost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            gmodstat.beinhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'beinhostisf%d' % e)]
    if gmod.typeemishost != 'none':
        gmodstat.ellphost = [[] for e in gmod.indxsersfgrd]
        gmodstat.anglhost = [[] for e in gmod.indxsersfgrd]
        gmodstat.serihost = [[] for e in gmod.indxsersfgrd]
        for e in gmod.indxsersfgrd:
            gmodstat.ellphost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'ellphostisf%d' % e)]
            gmodstat.anglhost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'anglhostisf%d' % e)]
            gmodstat.serihost[e] = gmodstat.paragenrscalfull[getattr(gmod.indxpara, 'serihostisf%d' % e)]
    if gmod.boollens:
        numbpixltemp = gdat.numbpixlcart
        defl = np.zeros((numbpixltemp, 2))
        
    # determine the indices of the pixels over which element kernels will be evaluated
    if gdat.boolbindspat:
        if gmod.numbpopl > 0:
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
                print('xposhost[e]')
                print(gmodstat.xposhost[e])
                print('yposhost[e]')
                print(gmodstat.yposhost[e])
                print('beinhost[e]')
                print(gmodstat.beinhost[e])
                print('gmodstat.ellphost[e]')
                print(gmodstat.ellphost[e])
                print('gmodstat.anglhost[e]')
                print(gmodstat.anglhost[e])

            deflhost[e] = chalcedon.retr_defl(gdat.xposgrid, gdat.yposgrid, indxpixlmiss, gmodstat.xposhost[e], gmodstat.yposhost[e], \
                                                                        gmodstat.beinhost[e], ellp=gmodstat.ellphost[e], angl=gmodstat.anglhost[e])
             
            if gdat.booldiag:
                if not np.isfinite(deflhost[e]).all():
                    print('')
                    print('')
                    print('')
                    print('gdat.xposgrid')
                    summgene(gdat.xposgrid)
                    print('gdat.yposgrid')
                    summgene(gdat.yposgrid)
                    print('indxpixlmiss')
                    summgene(indxpixlmiss)
                    print('xposhost[e]')
                    print(gmodstat.xposhost[e])
                    print('yposhost[e]')
                    print(gmodstat.yposhost[e])
                    print('beinhost[e]')
                    print(gmodstat.beinhost[e])
                    print('gmodstat.ellphost[e]')
                    print(gmodstat.ellphost[e])
                    print('gmodstat.anglhost[e]')
                    print(gmodstat.anglhost[e])
                    print('deflhost[e]')
                    print(deflhost[e])
                    summgene(deflhost[e])
                    raise Exception('not np.isfinite(deflhost[e]).all()')
        
            if gdat.booldiag:
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
        psfnconv = [[[] for i in gdat.indxener] for m in gdat.indxdqlt]
        if gdat.typepixl == 'cart':
            
            gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.blimpara.angl, gmod.typemodlpsfn, strgmodl)
            fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            for mm, m in enumerate(gdat.indxdqlt):
                for ii, i in enumerate(gdat.indxener):
                    if gmod.typemodlpsfn == 'singgaus':
                        sigm = psfp[i+m*gdat.numbener]
                    else:
                        sigm = fwhm[i, m] / 2.355
                    psfnconv[mm][ii] = AiryDisk2DKernel(sigm / gdat.sizepixl)
        
        stopchro(gdat, gdatmodi, 'psfnconv')
    
    if (gmod.typeevalpsfn == 'kern' or gmod.typeevalpsfn == 'full') and gmod.numbpopl > 0:
        if strgmodl == 'true' or boolinit or gdat.boolmodipsfn:
            if gdat.typepixl == 'heal':
                gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.blimpara.angl, gmod.typemodlpsfn, strgmodl)
                gmodstat.psfnintp = sp.interpolate.interp1d(gdat.blimpara.angl, gmodstat.psfn, axis=1, fill_value='extrapolate')
                fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            if gdat.typepixl == 'cart':
                if gdat.kernevaltype == 'ulip':
                    gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.blimpara.angl, gmod.typemodlpsfn, strgmodl)
                    gmodstat.psfnintp = sp.interpolate.interp1d(gdat.blimpara.angl, gmodstat.psfn, axis=1, fill_value='extrapolate')
                    if gdat.booldiag:
                        if not np.isfinite(gmodstat.psfnintp(0.05)).all():
                            raise Exception('')

                if gdat.kernevaltype == 'bspx':
                    
                    gmodstat.psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.blimpara.anglcart.flatten(), gmod.typemodlpsfn, strgmodl)
                    
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
        
    if gmod.numbpopl > 0:
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
                                np.exp(-0.5 * ((gmodstat.dictelem[l]['xpos'][k] - gdat.xposgrid[listindxpixlelem[l][k]])**2 + \
                                    (gmodstat.dictelem[l]['ypos'][k] - gdat.yposgrid[listindxpixlelem[l][k]])**2) / gmodstat.dictelem[l]['gwdt'][k]**2)
                            
                        if gmod.boolelempsfn[l]:
                            sbrtdfnc[:, listindxpixlelem[l][k], :] += retr_sbrtpnts(gdat, gmodstat.dictelem[l]['xpos'][k], \
                                                             gmodstat.dictelem[l]['ypos'][k], varbamplextd, gmodstat.psfnintp, listindxpixlelem[l][k])
                        
                        if gmod.typeelem[l].startswith('lghtline'):
                            sbrtdfnc[:, 0, 0] += gmodstat.dictelem[l]['spec'][:, k]
                        
            sbrt['dfnc'] = sbrtdfnc
            
            if gdat.booldiag:
                if not np.isfinite(sbrtdfnc).all():
                    raise Exception('Element delta function brightness not finite.')

            setattr(gmodstat, 'sbrtdfnc', sbrt['dfnc'])

            if gdat.booldiag:
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
                        deflsubh[indxpixl, :] += chalcedon.retr_defl(gdat.xposgrid, gdat.yposgrid, indxpixl, \
                                                     gmodstat.dictelem[l]['xpos'][kk], gmodstat.dictelem[l]['ypos'][kk], gmodstat.dictelem[l]['defs'][kk], \
                                                     asca=asca, acut=acut)
            
                    # temp -- find out what is causing the features in the element convergence maps
                    #for kk, k in enumerate(indxelem[l]):
                    #    indxpixlpnts = retr_indxpixl(gdat, gmodstat.dictelem[l]['ypos'][kk], gmodstat.dictelem[l]['xpos'][kk])
                    #    if deflsubh[listindxpixlelem[l][kk], :]
            
            if gdat.typeverb > 2:
                print('deflsubh')
                summgene(deflsubh)
            setattr(gmodstat, 'deflsubh', deflsubh)
            
            if gdat.booldiag:
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
                                    np.exp(-0.5 * ((gmodstat.dictelem[l]['xpos'][k] - gdat.xposgrid[None, listindxpixlelem[l][k], None])**2 + \
                                    (gmodstat.dictelem[l]['ypos'][k] - gdat.yposgrid[None, listindxpixlelem[l][k], None])**2) / gmodstat.dictelem[l]['gwdt'][k]**2)
                
                setattr(gmodstat, 'sbrtextsbgrd', sbrtextsbgrd)
            sbrt['extsbgrd'] = []
            sbrt['extsbgrd'] = sbrtextsbgrd
            
            if gdat.booldiag:
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
        
        if strgstat == 'this' or gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
            sbrt['bgrd'] = []
        if gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
            sbrt['bgrdgalx'] = []
        
        if gdat.numbener > 1:
            specsour = retr_spec(gdat, np.array([gmodstat.fluxsour]), sind=np.array([gmodstat.sindsour]))
            if gdat.typeverb > 2:
                print('gmodstat.sindsour')
                print(gmodstat.sindsour)
        else:
            specsour = np.array([gmodstat.fluxsour])
        
        if gdat.typeverb > 2:
            print('xpossour')
            print(xpossour)
            print('ypossour')
            print(ypossour)
            print('gmodstat.sizesour')
            print(gmodstat.sizesour)
            print('gmodstat.ellpsour')
            print(gmodstat.ellpsour)
            print('gmodstat.anglsour')
            print(gmodstat.anglsour)
            print('gmodstat.fluxsour')
            print(gmodstat.fluxsour)
            print('specsour')
            print(specsour)

        if gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
        
            if gdat.typeverb > 2:
                print('Interpolating the background emission...')

            sbrt['bgrdgalx'] = retr_sbrtsers(gdat, gdat.xposgrid[indxpixlelem[0]], gdat.yposgrid[indxpixlelem[0]], \
                                                                            xpossour, ypossour, specsour, gmodstat.sizesour, gmodstat.ellpsour, gmodstat.anglsour)
            
            if gdat.typeverb > 2:
                print('sbrt[bgrdgalx]')
                summgene(sbrt['bgrdgalx'])
                print('sbrtextsbgrd')
                summgene(sbrtextsbgrd)
            sbrt['bgrd'] = sbrt['bgrdgalx'] + sbrtextsbgrd
        
            sbrt['lens'] = np.empty_like(gdat.cntpdata)
            for ii, i in enumerate(gdat.indxener):
                for mm, m in enumerate(gdat.indxdqlt):
                    sbrtbgrdobjt = sp.interpolate.RectBivariateSpline(gdat.bctrpara.yposcart, gdat.bctrpara.xposcart, \
                                                            sbrt['bgrd'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    yposprim = gdat.yposgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 1]
                    xposprim = gdat.xposgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 0]
                    # temp -- T?
                    sbrt['lens'][ii, :, m] = sbrtbgrdobjt(yposprim, xposprim, grid=False).flatten()
        else:
            if gdat.typeverb > 2:
                print('Not interpolating the background emission...')
            
            sbrt['lens'] = retr_sbrtsers(gdat, gdat.xposgrid - defl[gdat.indxpixl, 0], \
                                                   gdat.yposgrid - defl[gdat.indxpixl, 1], \
                                                   xpossour, ypossour, specsour, gmodstat.sizesour, gmodstat.ellpsour, gmodstat.anglsour)
            
            sbrt['bgrd'] = retr_sbrtsers(gdat, gdat.xposgrid, \
                                                   gdat.yposgrid, \
                                                   xpossour, ypossour, specsour, gmodstat.sizesour, gmodstat.ellpsour, gmodstat.anglsour)
            
        setattr(gmodstat, 'sbrtlens', sbrt['lens'])

        if gdat.booldiag:
            if not np.isfinite(sbrt['lens']).all():
                raise Exception('Lensed emission is not finite.')
            if (sbrt['lens'] == 0).all():
                raise Exception('Lensed emission is zero everywhere.')

        stopchro(gdat, gdatmodi, 'sbrtlens')
        
    ### background surface brightness
    sbrtback = []
    # temp
    #sbrtback = np.empty((numbback, gdat.numbener, indxpixlelem[yy].size, gdat.numbdqlt))
    
    # evaluate host galaxy surface brightness
    if gmod.typeemishost != 'none':
        initchro(gdat, gdatmodi, 'sbrthost')
        for e in gmod.indxsersfgrd:
            if gdat.typeverb > 2:
                print('Evaluating the host galaxy surface brightness...')
            if gdat.numbener > 1:
                spechost = retr_spec(gdat, np.array([gmodstat.fluxhost[e]]), sind=np.array([gmodstat.sindhost[e]]))
            else:
                spechost = np.array([gmodstat.fluxhost[e]])
            
            if gdat.typeverb > 2:
                print('xposhost[e]')
                print(gmodstat.xposhost[e] * gdat.anglfact)
                print('yposhost[e]')
                print(gmodstat.yposhost[e] * gdat.anglfact)
                print('spechost')
                print(spechost)
                print('gmodstat.sizehost[e]')
                print(gmodstat.sizehost[e])
                print('gmodstat.ellphost[e]')
                print(gmodstat.ellphost[e])
                print('gmodstat.anglhost[e]')
                print(gmodstat.anglhost[e])
                print('gmodstat.serihost[e]')
                print(gmodstat.serihost[e])
            
            sbrt['hostisf%d' % e] = retr_sbrtsers(gdat, gdat.xposgrid, gdat.yposgrid, gmodstat.xposhost[e], \
                                                                        gmodstat.yposhost[e], spechost, gmodstat.sizehost[e], gmodstat.ellphost[e], gmodstat.anglhost[e], gmodstat.serihost[e])
            
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
    
    sbrt['modlraww'] = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbdqlt))
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
        if gdat.booldiag:
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
            for mm, m in enumerate(gdat.indxdqlt):
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
        sbrt['modlconv'] = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbdqlt))
        for ii, i in enumerate(gdat.indxener):
            for mm, m in enumerate(gdat.indxdqlt):
                if gdat.strgcnfg == 'pcat_ferm_igal_simu_test':
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
    if gmod.numbpopl > 0 and gmod.boolelemsbrtdfncanyy:
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
    
    if gdat.booldiag:
        setattr(gmodstat, 'cntpmodl', cntp['modl'])
    stopchro(gdat, gdatmodi, 'expo')

    # simulated-data specific
    if strgmodl == 'true' and strgstat == 'this':
        
        # generate count data
        cntptemp = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbdqlt))
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxdqlt:
                    cntptemp[i, j, m] = np.random.poisson(cntp['modl'][i, j, m])
        setattr(gdat, 'cntpdata', cntptemp)
        
        if np.amax(cntptemp) == 0:
            print('')
            print('')
            print('')
            print('cntptemp')
            summgene(cntptemp)
            print('sbrt[modl]')
            summgene(sbrt['modl'])
            print('gdat.apix')
            print(gdat.apix)
            print('gdat.expo')
            summgene(gdat.expo)
            print('cntp[modl][i, j, m]')
            summgene(cntp['modl'][i, j, m])
            raise Exception('Generated data has zero counts everywhere.')

        print('Will process the true model...')
        proc_cntpdata(gdat)
    
    ## diagnostics
    if gdat.booldiag:
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
    if gmod.numbpopl > 0:
        
        # parsimony
        for l in gmod.indxpopl:
            lpri[0] -= 0.5 * gdat.factpriodoff * gmod.numbparagenrelemsing[l] * gmodstat.numbelem[l]
        
    # power spectrum of the model
    if gdat.boolpenalpridiff:
        sbrtdatapnts = gdat.sbrtdata - sbrt['dfnc']
        if gdat.typepixl == 'heal':
            raise Exception('')
        if gdat.typepixl == 'cart':
            psecodimdatapnts = np.empty((gdat.numbener, gdat.numbsidecarthalf, gdat.numbdqlt))
            psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.blimpara.angl, gmod.typemodlpsfn, strgmodl)
            fwhm = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            sigm = fwhm / 2.355
            psecodimdatapntsprio = np.exp(-2. * gdat.bctrpara.mpolodim[None, :, None] / (0.1 / sigm[:, None, :]))
            lpridiff = 0.
            for i in gdat.indxener:
                for m in gdat.indxdqlt:
                    psecdatapnts = retr_psec(gdat, sbrtdatapnts[i, :, m])
                    psecodimdatapnts[i, :, m] = retr_psecodim(gdat, psecdatapnts)
                    psecodimdatapnts[i, :, m] /= psecodimdatapnts[i, 0, m]
                    lpridiff += -0.5 * np.sum((psecodimdatapnts[i, :, m] - psecodimdatapntsprio[i, :, m])**2)
                    setattr(gmodstat, 'psecodimdatapntsen%02devt%d' % (i, m), psecodimdatapnts[i, :, m])
                    setattr(gmodstat, 'psecodimdatapntsprioen%02devt%d'% (i, m), psecodimdatapntsprio[i, :, m])
        lpri[1] = lpridiff 
        setattr(gmodstat, 'lpridiff', lpridiff)
        
    if gmod.numbpopl > 0:
        # Poission prior on the number of elements
        if gmod.typemodltran == 'pois':
            meanelem = gmodstat.paragenrscalfull[gmod.indxpara.meanelem]
            for l in gmod.indxpopl:
                lpri[2] += retr_lprbpois(gmodstat.numbelem[l], meanelem[l])
        
        for l in gmod.indxpopl:
            for g, (strgfeat, strgpdfn) in enumerate(zip(gmod.namepara.genrelem[l], gmod.scalpara.genrelem[l])):
                indxlpritemp = 0# + l * gmod.numbparagenrelempopl + g
                lpri[indxlpritemp] = retr_lprielem(gdat, strgmodl, l, g, strgfeat, strgpdfn, gmodstat.paragenrscalfull, gmodstat.dictelem, gmodstat.numbelem)
    lpritotl = np.sum(lpri)
    
    if gdat.typeverb > 1:
        print('lpritotl')
        print(lpritotl)
    
    ### log-likelihood
    initchro(gdat, gdatmodi, 'llik')
    llik = retr_llik_bind(gdat, strgmodl, cntp['modl'])
    
    if gdat.typeverb > 2:
        print('cntp[modl]')
        summgene(cntp['modl'])
        print('np.sum(cntp[modl], (1, 2))')
        print(np.sum(cntp['modl'], (1, 2)))
        print('np.sum(gdat.cntpdata, (1, 2))')
        print(np.sum(gdat.cntpdata, (1, 2)))

    if gdat.booldiag:
        if not np.isfinite(llik).all():
            raise Exception('Likelihood is not finite.')
    
    gmodstat.lliktotl = np.sum(llik)
    if gdat.booldiag:
        if isinstance(gmodstat.lliktotl, np.ndarray):
            raise Exception('')
        if not np.isfinite(gmodstat.lliktotl).all():
            raise Exception('')

    numbdoff = gdat.numbdata - gmod.numbparagenrbase
    if gmod.numbpopl > 0:
        for l in gmod.indxpopl:
            numbdoff -= len(gmodstat.indxparagenrelemfull[l]['full'])

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
    
    if gmod.numbpopl > 0:
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
            masshostbein = gmod.ratimassbeinsqrd * gmodstat.beinhost[e]**2
            setattr(gmodstat, 'masshostisf%dbein' % e, masshostbein)
        ### sort with respect to deflection at scale radius
        if gmod.numbpopl > 0:
            for l in gmod.indxpopl:
                if gmodstat.numbelem[l] > 0:
                    indxelemsortampl = np.argsort(gmodstat.dictelem[l][gmod.nameparaelemsort[l]])[::-1]
                    for nameparagenrelem in gmod.namepara.genrelem[l]:
                        gmodstat.dictelem[l][nameparagenrelem + 'sort'] = gmodstat.dictelem[l][nameparagenrelem][indxelemsortampl]

        deflsing = np.zeros((gdat.numbpixlcart, 2, gmod.numbdeflsingplot))
        conv = np.zeros((gdat.numbpixlcart))
        convpsec = np.zeros(((gdat.numbsidecarthalf)**2))
        convpsecodim = np.zeros((gdat.numbsidecarthalf))
        if gmod.numbpopl > 0:
            if gmod.boolelemlens:
                gmod.indxpopllens = gmod.typeelem.index('lens')
        numbdeflsing = 2
        if gmod.numbpopl > 0:
            if gmod.boolelemlens:
                if gmodstat.numbelem[gmod.indxpopllens] > 0:
                    numbdeflsing += min(gdat.numbdeflsubhplot, gmodstat.numbelem[gmod.indxpopllens]) 
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
                        asca = gmodstat.dictelem[gmod.indxpopllens]['ascasort'][None, k-3]
                        acut = gmodstat.dictelem[gmod.indxpopllens]['acutsort'][None, k-3]
                        deflsing[listindxpixlelem[gmod.indxpopllens][k], :, k] = chalcedon.retr_defl(gdat.xposgrid, gdat.yposgrid, listindxpixlelem[gmod.indxpopllens][k], \
                                                      gmodstat.dictelem[gmod.indxpopllens]['xpossort'][None, k-3], gmodstat.dictelem[gmod.indxpopllens]['ypossort'][None, k-3], \
                                                      gmodstat.dictelem[gmod.indxpopllens]['defssort'][None, k-3], asca=asca, acut=acut)

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
        if gmod.numbpopl > 0:
            if gmod.boolelemlens:
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
        deflsingmgtd = np.zeros((gdat.numbpixlcart, gmod.numbdeflsingplot))
        magn = 1. / retr_invm(gdat, defl) 
        histdefl = np.histogram(defl, bins=gmod.blimpara.defl)[0]
        deflsingmgtd = np.sqrt(np.sum(deflsing[...]**2, axis=1))
        if gmod.numbpopl > 0:
            if gmod.boolelemlens:
                histdeflsubh = np.histogram(deflsubh, bins=gdat.blimpara.deflsubh)[0]
                setattr(gmodstat, 'histdeflsubh', histdeflsubh)
        setattr(gmodstat, 'histdefl', histdefl)
        setattr(gmodstat, 'magn', magn)
        setattr(gmodstat, 'deflsing', deflsing)
        setattr(gmodstat, 'deflsingmgtd', deflsingmgtd)
    
    ## element related
    if gmod.numbpopl > 0:
        if gdat.numbpixl == 1:
            for l in gmod.indxpopl:
                for k in range(gmodstat.numbelem[l]):
                    setattr(gmodstat, 'speclinepop%d%04d' % (l, k), gmodstat.dictelem[l]['spec'][:, k])
        
        if gdat.typedata == 'simu' and strgmodl == 'true' and gdat.numbpixl > 1:
            gdat.refrxpos = [[] for l in gmod.indxpopl]
            gdat.refrypos = [[] for l in gmod.indxpopl]
            for l in gmod.indxpopl:
                gdat.refrxpos[l] = np.tile(gmodstat.dictelem[l]['xpos'], [3] + list(np.ones(gmodstat.dictelem[l]['xpos'].ndim, dtype=int)))
                gdat.refrypos[l] = np.tile(gmodstat.dictelem[l]['ypos'], [3] + list(np.ones(gmodstat.dictelem[l]['ypos'].ndim, dtype=int)))
    
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lghtpntspuls':
                gmodstat.dictelem[l]['per1'] = retr_per1(gmodstat.dictelem[l]['per0'], gmodstat.dictelem[l]['magf'])
        
    if gmod.numbpopl > 0:
        if strgstat == 'this' or gdat.boolrefeforc and strgmodl == 'fitt':
            # correlate the fitting model elements with the reference elements
            if gdat.boolinforefr and not (strgmodl == 'true' and gdat.typedata == 'simu') and gdat.boolasscrefr:
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
                                matrdist[:, k] = retr_angldist(gdat, gdat.refr.dictelem[q]['xpos'][0, :], gdat.refr.dictelem[q]['ypos'][0, :], gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k])
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
                        
                        # collect the error in the associated reference element amplitude
                        if len(indxelemfittasschits[q][l]) > 0:
                            gmodstat.dictelem[l]['aerr' + gdat.listnamerefr[q]] = np.zeros(gmodstat.numbelem[l])
                            fittfeattemp = gmodstat.dictelem[l][gmod.nameparagenrelemampl[l]][indxelemfittasschits[q][l]]
                            refrfeattemp = gdat.refr.dictelem[q][gmod.nameparagenrelemampl[l]][0, indxelemrefrasschits[q][l]]
                            if gdat.booldiag:
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
    if gdat.booldiag:
        if not np.isfinite(cntp['resi']).all():
            raise Exception('')
        if not np.isfinite(numbdoff):
            raise Exception('')
        if not np.isfinite(chi2doff):
            raise Exception('')
    setattr(gmodstat, 'numbdoff', numbdoff)
    setattr(gmodstat, 'chi2doff', chi2doff)
    
    if gmod.boolelempsfn and gmod.numbpopl > 0:
        gmodstat.fwhmpsfn = 2. * retr_psfnwdth(gdat, gmodstat.psfn, 0.5)
            
    if gmod.numbpopl > 0:
        
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
            if gdat.boolbindspat:
                #### radial and angular coordinates
                gmodstat.dictelem[l]['gang'] = retr_gang(gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos'])
                gmodstat.dictelem[l]['aang'] = retr_aang(gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos'])
            
            if gmod.boolelemlght[l]:
                #### number of expected counts
                if gdat.boolbindspat:
                    gmodstat.dictelem[l]['cnts'] = retr_cntspnts(gdat, [gmodstat.dictelem[l]['xpos'], gmodstat.dictelem[l]['ypos']], gmodstat.dictelem[l]['spec'])
                else:
                    gmodstat.dictelem[l]['cnts'] = retr_cntspnts(gdat, [gmodstat.dictelem[l]['elin']], gmodstat.dictelem[l]['spec'])
            
            #### delta log-likelihood
            gmodstat.dictelem[l]['deltllik'] = np.zeros(gmodstat.numbelem[l])
            if not (strgmodl == 'true' and gdat.checprio): 
                if gdat.typeverb > 2:
                    print('Calculating log-likelihood differences when removing elements from the model.')
                for k in range(gmodstat.numbelem[l]):
                    
                    # construct gdatmodi
                    gdatmoditemp = None
                    gdat.true.next.indxpara = tdpy.gdatstrt()

                    gdatmoditemp = tdpy.gdatstrt()
                    gdatmoditemp.this = tdpy.gdatstrt()
                    gdatmoditemp.next = tdpy.gdatstrt()
                    gdatmoditemp.this.indxpara = tdpy.gdatstrt()
                    gdatmoditemp.next.indxpara = tdpy.gdatstrt()
                    gdatmoditemp.this.indxelemfull = gmodstat.indxelemfull
                    gdatmoditemp.this.paragenrscalfull = gmodstat.paragenrscalfull
                    gdatmoditemp.this.paragenrunitfull = gmodstat.paragenrunitfull

                    prop_stat(gdat, gdatmoditemp, strgmodl, deth=True, thisindxpopl=l, thisindxelem=k)
                    proc_samp(gdat, gdatmoditemp, 'next', strgmodl)#, boolinit=boolinit)
                    
                    if gdat.booldiag:
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
    
    if gmod.numbpopl > 0 and gmod.boolelemsbrtdfncanyy:
        
        if gmod.numbpopl > 0:
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
                            sbrttemp = retr_sbrtpnts(gdat, gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
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
                            sbrt['dfncfull'][:, :, :] += retr_sbrtpnts(gdat, gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
                                                                                                            varbamplextd, gmodstat.psfnintp, gdat.indxpixl)
            
                setattr(gmodstat, 'sbrtdfncsubtpop%d' % l, sbrt['dfncsubt'])
                
    if gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
        if gdat.booldiag:
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
    
    if gmod.numbpopl > 0:
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
        for m in gdat.indxdqlt:
            if gdat.numbdqlt > 1:
                namehistcntp += 'evt%d' % m
            for i in gdat.indxener: 
                if gdat.numbener > 1:
                    namehistcntp += 'en%02d' % i
                
                histcntp = np.histogram(cntp[name][i, :, m], bins=gdat.blimpara.cntpmodl)[0]
                setattr(gmodstat, namehistcntp, histcntp)
                
                if False and i == 0 and m == 0 and (name == 'dfnc' or name == 'dfncsubt'):
                    for strgbins in ['lowr', 'higr']:
                        strgtemp = 'histcntp' + strgbins + name + 'en%02devt%d' % (i, m)
                        if strgbins == 'lowr':
                            setattr(gmod, strgtemp, np.array([float(np.sum(histcntp[:gdat.numbtickcbar-1]))]))
                        else:
                            setattr(gmod, strgtemp, np.array([float(np.sum(histcntp[gdat.numbtickcbar-1:]))]))
            else:
                histcntp = np.histogram(cntp[name][:, 0, m], bins=gdat.blimpara.cntpmodl)[0]
                setattr(gmodstat, 'histcntp' + name + 'evt%d' % m, histcntp)

    if gmod.boollens:
        if strgmodl == 'true':
            s2nr = []
            s2nr = cntp['lens'] / np.sqrt(cntp['modl'])
            setattr(gmodstat, 's2nr', s2nr)
        cntplensgrad = np.empty((gdat.numbener, gdat.numbpixlcart, gdat.numbdqlt, 2))
        for i in gdat.indxener:
            for m in gdat.indxdqlt:
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

    if gmod.numbpopl > 0:
        for l in gmod.indxpopl:
            if gmod.boolelemlght[l]:
                #### spectra
                if gdat.boolbindspat:
                    sindcolr = [gmodstat.dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    gmodstat.dictelem[l]['specplot'] = retr_spec(gdat, gmodstat.dictelem[l]['flux'], sind=gmodstat.dictelem[l]['sind'], \
                                                                 curv=gmodstat.dictelem[l]['curv'], expc=gmodstat.dictelem[l]['expc'], \
                                                                 sindcolr=sindcolr, spectype=gmod.spectype[l], plot=True)
                
                if gdat.typedata == 'inpt':
                    if gdat.typeexpr == 'ferm':
                        # temp
                        try:
                            gmodstat.dictelem[l]['sbrt0018'] = gdat.sbrt0018objt(gmodstat.dictelem[l]['ypos'], gmodstat.dictelem[l]['xpos'])
                        except:
                            gmodstat.dictelem[l]['sbrt0018'] = gmodstat.dictelem[l]['ypos'] * 0.

            if gmod.typeelem[l] == 'lens':
                #### distance to the source
                if gmod.boollens:
                    gmodstat.dictelem[l]['distsour'] = retr_angldist(gdat, gmodstat.dictelem[l]['xpos'],  gmodstat.dictelem[l]['ypos'], xpossour, ypossour)
                
                if gmod.boollenssubh:
                    gmodstat.dictelem[l]['deflprof'] = np.empty((gdat.numbbinspara.anglfull, gmodstat.numbelem[l]))
                    gmodstat.dictelem[l]['mcut'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['rele'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['reln'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relk'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relf'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['reld'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relc'] = np.empty(gmodstat.numbelem[l])
                    gmodstat.dictelem[l]['relm'] = np.empty(gmodstat.numbelem[l])

                    # temp -- this can be placed earlier in the code
                    cntplensobjt = sp.interpolate.RectBivariateSpline(gdat.bctrpara.yposcart, gdat.bctrpara.xposcart, \
                                                            cntp['lens'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    for k in np.arange(gmodstat.numbelem[l]):
                        
                        asca = gmodstat.dictelem[l]['asca'][k]
                        acut = gmodstat.dictelem[l]['acut'][k]
                        
                        #### deflection profiles
                        gmodstat.dictelem[l]['deflprof'][:, k] = chalcedon.retr_deflcutf(gdat.bctrpara.anglfull, gmodstat.dictelem[l]['defs'][k], asca, acut)
         
                        ### truncated mass 
                        gmodstat.dictelem[l]['mcut'][k] = aspendos.retr_mcut(gdat, gmodstat.dictelem[l]['defs'][k], asca, acut, gmod.adislens, gmod.mdencrit)

                        #### relevance, the dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        gmodstat.dictelem[l]['rele'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
                                                                              gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        
                        #gmodstat.dictelem[l]['relf'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
                        #                                                 gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, cntpmodl=cntp['modl'][0, :, 0])
                        #
                        #deflelem = chalcedon.retr_defl(gdat.xposgrid, gdat.yposgrid, gdat.indxpixl, gmodstat.dictelem[l]['xpos'][k], \
                        #                                                    gmodstat.dictelem[l]['ypos'][k], gmodstat.dictelem[l]['defs'][k], asca=asca, acut=acut)
                        #yposprim = gdat.yposgrid - deflelem[:, 1]
                        #xposprim = gdat.xposgrid - deflelem[:, 0]
                        #gmodstat.dictelem[l]['relm'][k] = np.mean(abs(cntp['lens'][0, :, 0] - cntplensobjt(yposprim, xposprim, grid=False).flatten()))
                        #
                        #
                        #gmodstat.dictelem[l]['relk'][k] = gmodstat.dictelem[l]['relm'][k] / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
                        #gmodstat.dictelem[l]['reln'][k] = gmodstat.dictelem[l]['rele'][k] / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
                        #gmodstat.dictelem[l]['reld'][k] = retr_rele(gdat, gdat.cntpdata[0, :, 0], gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
                        #                                                                                 gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        #gmodstat.dictelem[l]['relc'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], gmodstat.dictelem[l]['xpos'][k], gmodstat.dictelem[l]['ypos'][k], \
                        #                               gmodstat.dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, absv=False) / gmodstat.dictelem[l]['defs'][k] * gdat.sizepixl
               
        ### distribution of element parameters and features
        #### calculate the model filter
        listindxelemfilt = [[[] for l in gmod.indxpopl] for namefilt in gdat.listnamefilt]
        for k, namefilt in enumerate(gdat.listnamefilt):
            for l in gmod.indxpopl:
                if namefilt == '':
                    listindxelemfilt[k][l] = np.arange(gmodstat.numbelem[l])
                if namefilt == 'imagbndr':
                    listindxelemfilt[k][l] = np.where((np.fabs(gmodstat.dictelem[l]['xpos']) < gdat.maxmgangdata) & (np.fabs(gmodstat.dictelem[l]['ypos']) < gdat.maxmgangdata))[0]
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
                        histfrst[:, i] = np.histogram(gmodstat.dictelem[l]['spec'][i, listindxelemfilt[0][l]], gdat.blimpara.spec)[0]
                elif namefrst == 'cnts':
                    histfrst = np.histogram(gmodstat.dictelem[l]['cnts'][listindxelemfilt[0][l]], gdat.blimpara.cnts)[0]
                else:
                #elif not (namefrst == 'curv' and gmod.spectype[l] != 'curv' or namefrst == 'expc' \
                #                                        and gmod.spectype[l] != 'expc' or namefrst.startswith('sindarry') and \
                #                                                                                          gmod.spectype[l] != 'colr'):
                    blimfeatfrst = getattr(gdat.blimpara, namefrst)
                    #if len(gmodstat.dictelem[l][namefrst]) > 0 and len(listindxelemfilt[0][l]) > 0:
                    histfrst = np.histogram(gmodstat.dictelem[l][namefrst][listindxelemfilt[0][l]], blimfeatfrst)[0]
                    strgvarb = 'hist' + namefrst + 'pop%d' % l
                    setattr(gmodstat, strgvarb, histfrst)
                        
                #### two dimensional
                for nameseco in gmod.namepara.elem[l]:
                    if namefrst == 'spec' or namefrst == 'specplot' or namefrst == 'deflprof' or \
                            nameseco == 'spec' or nameseco == 'specplot' or nameseco == 'deflprof':
                        continue
                    
                    if not checstrgfeat(namefrst, nameseco):
                        continue

                    blimfeatseco = getattr(gdat.blimpara, nameseco)
                    histtdim = np.histogram2d(gmodstat.dictelem[l][namefrst][listindxelemfilt[0][l]], \
                                                            gmodstat.dictelem[l][nameseco][listindxelemfilt[0][l]], [blimfeatfrst, blimfeatseco])[0]
            
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
                            sexp = gmodstat.paragenrscalfull[getattr(gmod.indxpara.genrbase, nameparagenrelem + 'distscal')[l]]
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
    
    if gmod.numbpopl > 0:
        for l in gmod.indxpopl:
            if gmod.typeelem[l] == 'lens':
                if gmodstat.numbelem[l] > 0:
                    ## total truncated mass of the subhalo as a cross check
                    # temp -- generalize
                    asca = gmodstat.dictelem[l]['asca']
                    acut = gmodstat.dictelem[l]['acut']
                    factmcutfromdefs = chalcedon.retr_factmcutfromdefs(gmod.adissour, gmod.adislens, gmod.adislenssour, asca, acut) 
                    masssubh = np.array([np.sum(factmcutfromdefs * gmodstat.dictelem[l]['defs'])])
    
    ## derived variables as a function of other derived variables
    if gmod.numbpopl > 0:
        for l in gmod.indxpopl:
            if gmod.typeelem[l].startswith('lghtpntspuls'):
                massshel = np.empty(gdat.numbanglhalf)
                for k in gdat.indxanglhalf:
                    indxelemshel = np.where((gdat.blimpara.anglhalf[k] < gmodstat.dictelem[l]['gang']) & (gmodstat.dictelem[l]['gang'] < gdat.blimpara.anglhalf[k+1]))
                    massshel[k] = np.sum(gmodstat.dictelem[l]['mass'][indxelemshel])
                setattr(gmodstat, 'massshelpop%d' % l, massshel)
            
    if gmod.boollens or gmod.numbpopl > 0 and gmod.boollenssubh:
        # find the host, subhalo masses and subhalo mass fraction as a function of halo-centric radius
        listnametemp = gdat.liststrgcalcmasssubh
        listnamevarbmass = []
        listnamevarbmassscal = []
        listnamevarbmassvect = []
        for e in gmod.indxsersfgrd:
            listnamevarbmassscal += ['masshosttotl']
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masshostisf%d' % e + strgtemp)
                listnamevarbmassscal.append('masshostisf%d' % e + strgtemp + 'bein')
        if gmod.numbpopl > 0 and gmod.boollenssubh:
            listnamevarbmassscal.append('masssubhtotl')
            listnamevarbmassscal.append('fracsubhtotl')
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masssubh' + strgtemp)
                listnamevarbmassvect.append('fracsubh' + strgtemp)
                listnamevarbmassscal.append('masssubh' + strgtemp + 'bein')
                listnamevarbmassscal.append('fracsubh' + strgtemp + 'bein')
        
        for name in listnamevarbmassvect:
            
            #dicttert[name] = np.zeros(gdat.numbanglhalf)
            dicttert[name] = np.zeros(gdat.numbbinspara.anglhalf)
            if 'isf' in name:
                indxisfrtemp = int(name.split('isf')[1][0])
            angl = np.sqrt((gdat.bctrpara.xposcartmesh - gmodstat.xposhost[indxisfrtemp])**2 + (gdat.bctrpara.yposcartmesh - gmodstat.yposhost[indxisfrtemp])**2).flatten()
            for k in gdat.indxbinspara.anglhalf:
                if name[4:8] == 'host':
                    convtemp = conv[:]
                if name[4:8] == 'subh':
                    convtemp = convelem[:]
                
                if name.endswith('delt'):
                    indxpixl = np.where((gdat.blimpara.anglhalf[k] < angl) & (angl < gdat.blimpara.anglhalf[k+1]))[0]
                    dicttert[name][k] = 1e6 * np.sum(convtemp[indxpixl]) * gmod.mdencrit * \
                                                gdat.apix * gmod.adislens**2 / 2. / np.pi * gdat.deltpara.anglhalf[k] / gdat.bctrpara.anglhalf[k]
                if name.endswith('intg'):
                    indxpixl = np.where(angl < gdat.bctrpara.anglhalf[k])[0]
                    dicttert[name][k] = np.sum(convtemp[indxpixl]) * gmod.mdencrit * gdat.apix * gmod.adislens**2
                
                if name[:4] == 'frac':
                    masshosttotl = 0.
                    for e in gmod.indxsersfgrd:
                        masshosttotl += dicttert['masshostisf%d' % e + name[-4:]][k]
                    if masshosttotl != 0.:
                        dicttert['fracsubh' + name[8:]][k] = dicttert['masssubh' + name[8:]][k] / masshosttotl
            setattr(gmodstat, name, dicttert[name])
            
            # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
            dicttert[name + 'bein'] = np.interp(gmodstat.beinhost, gdat.bctrpara.anglhalf, dicttert[name])
            setattr(gmodstat, name + 'bein', dicttert[name + 'bein'])
        
    #if gmod.numbpopl > 0:
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
    
    if gmod.numbpopl > 0 and gdat.factpriodoff != 0.:
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

                    blim = getattr(gdat.blimpara, strgfeat)
                    if len(indxelempars) > 0:
                        refrhistpars = np.histogram(gmodstat.dictelem[q][strgfeat][indxelempars], blim=blim)[0].astype(float)
                        if indxrefrgood.size > 0:
                            reca[indxrefrgood] = refrhistpars[indxrefrgood] / refrhist[indxrefrgood]
                    
                    setattr(gmodstat, 'histpars' + strgfeat + 'pop%d' % q, refrhistpars)
                    setattr(gmodstat, 'reca' + strgfeat + 'pop%d' % q, reca)
        
        print('gdat.strgcnfgsimu')
        print(gdat.strgcnfgsimu)
        if gdat.strgcnfgsimu is not None:
            if gmod.numbpopl > 0:
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
    if strgmodl == 'fitt' and gdat.typedata == 'simu':
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
            if gmod.numbpopl > 0:
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
    
    if gmod.numbpopl > 0:
    
        # correlate the catalog sample with the reference catalog
        if gdat.boolinforefr and not (strgmodl == 'true' and gdat.typedata == 'simu') and gdat.boolasscrefr:
            
            for q in gdat.indxrefr:
                for l in gmod.indxpopl:
                    if gdat.refr.numbelem[q] > 0:
                        cmpl = np.array([float(len(indxelemrefrasschits[q][l])) / gdat.refr.numbelem[q]])
                        if gdat.booldiag:
                            if cmpl > 1. or cmpl < 0.:
                                raise Exception('')
                    else:
                        cmpl = np.array([-1.])
                    setattr(gmodstat, 'cmplpop%dpop%d' % (l, q), cmpl)
                    if gmodstat.numbelem[l] > 0:
                        fdis = np.array([float(indxelemfittasscfals[q][l].size) / gmodstat.numbelem[l]])
                        if gdat.booldiag:
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
                    blimfeatfrst = getattr(gdat.blimpara, nameparaelemfrst)
                        
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
                            blimfeatseco = getattr(gdat.blimpara, nameparaelemseco)
                            
                            refrhistfeattdimassc = np.histogram2d(refrfeatfrst[indxelemrefrasschits[q][l]], \
                                                                  refrfeatseco[indxelemrefrasschits[q][l]], bins=(blimfeatfrst, blimfeatseco))[0]
                            indxgood = np.where(refrhistfeattdim != 0.)
                            if indxgood[0].size > 0:
                                cmpltdim[indxgood] = refrhistfeattdimassc[indxgood].astype(float) / refrhistfeattdim[indxgood]
                                if gdat.booldiag:
                                    if np.where((cmpltdim[indxgood] > 1.) | (cmpltdim[indxgood] < 0.))[0].size > 0:
                                        raise Exception('')
                        setattr(gmodstat, 'cmpl%s%spop%d' % (nameparaelemfrst, nameparaelemseco, q), cmpltdim)

                    cmplfrst = np.zeros(gdat.numbbinsplot) - 1.
                    if len(indxelemrefrasschits[q][l]) > 0:
                        refrhistfeatfrst = getattr(gdat.refr, 'hist' + nameparaelemfrst + 'pop%d' % q)
                        blimfeatfrst = getattr(gdat.blimpara, nameparaelemfrst)
                        refrhistfeatfrstassc = np.histogram(refrfeatfrst[indxelemrefrasschits[q][l]], blim=blimfeatfrst)[0]
                        indxgood = np.where(refrhistfeatfrst != 0.)[0]
                        if indxgood.size > 0:
                            cmplfrst[indxgood] = refrhistfeatfrstassc[indxgood].astype(float) / refrhistfeatfrst[indxgood]
                            if gdat.booldiag:
                                if np.where((cmplfrst[indxgood] > 1.) | (cmplfrst[indxgood] < 0.))[0].size > 0:
                                    raise Exception('')
                   
                    setattr(gmodstat, 'cmpl%spop%d' % (nameparaelemfrst, q), cmplfrst)
            
            # false discovery rate
            for l in gmod.indxpopl:
                q = gmod.indxpoplrefrassc[l]
                
                for nameparaelemfrst in gmod.namepara.elem[l]:
                    
                    blimfeatfrst = getattr(gdat.blimpara, nameparaelemfrst)
                    for nameparaelemseco in gmod.namepara.elem[l]:
                        
                        if not checstrgfeat(nameparaelemfrst, nameparaelemseco):
                            continue
                        
                        # temp -- the size of the fdis np.array should depend on strgmodl
                        fdistdim = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot))
                        
                        if len(indxelemrefrasschits[q][l]) > 0 and len(gmodstat.dictelem[l][nameparaelemseco]) > 0 and len(gmodstat.dictelem[l][nameparaelemfrst]) > 0: 
                            strgfeattdim = nameparaelemfrst + nameparaelemseco + 'pop%d' % l
                            fitthistfeattdim = getattr(gmodstat, 'hist' + strgfeattdim)
                            blimfeatseco = getattr(gdat.blimpara, nameparaelemseco)
                            
                            fitthistfeattdimfals = np.histogram2d(gmodstat.dictelem[l][nameparaelemfrst][indxelemfittasscfals[q][l]], \
                                                  gmodstat.dictelem[l][nameparaelemseco][indxelemfittasscfals[q][l]], bins=(blimfeatfrst, blimfeatseco))[0]
                            indxgood = np.where(fitthistfeattdim != 0.)
                            if indxgood[0].size > 0:
                                fdistdim[indxgood] = fitthistfeattdimfals[indxgood].astype(float) / fitthistfeattdim[indxgood]
                                if gdat.booldiag:
                                    if np.where((fdistdim[indxgood] > 1.) | (fdistdim[indxgood] < 0.))[0].size > 0:
                                        raise Exception('')
                        
                        setattr(gmodstat, 'fdis%s%spop%d' % (nameparaelemfrst, nameparaelemseco, l), fdistdim)
                
                    fdisfrst = np.zeros(gdat.numbbinsplot)
                    if len(indxelemrefrasschits[q][l]) > 0 and len(gmodstat.dictelem[l][nameparaelemfrst]) > 0:
                        blimfeatfrst = getattr(gdat.blimpara, nameparaelemfrst)
                        fitthistfeatfrstfals = np.histogram(gmodstat.dictelem[l][nameparaelemfrst][indxelemfittasscfals[q][l]], blim=blimfeatfrst)[0]
                        fitthistfeatfrst = getattr(gmodstat, 'hist' + nameparaelemfrst + 'pop%d' % l)
                        indxgood = np.where(fitthistfeatfrst != 0.)[0]
                        if indxgood.size > 0:
                            fdisfrst[indxgood] = fitthistfeatfrstfals[indxgood].astype(float) / fitthistfeatfrst[indxgood]
                            if gdat.booldiag:
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
                            gmod.indxpara.genrelemtemp = gmod.namepara.genrelem[l].index(strgfeat)
                            if (gmod.scalpara.genrelem[l][gmod.indxpara.genrelemtemp] != 'gaus' and not gmod.scalpara.genrelem[l][gmod.indxpara.genrelemtemp].startswith('lnor')):
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
        maxmypos = getattr(gmod, 'maxmypos')
        gmod.indxpara.yposdistscal = getattr(gmod.indxpara, 'yposdistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_dnp.exp(dictelem[l]['ypos'], maxmypos, paragenrscalfull[gmod.indxpara.yposdistscal]))) 
    if strgpdfn == 'expo':
        maxmgang = getattr(gmod, 'maxmgang')
        gang = retr_gang(dictelem[l]['xpos'], dictelem[l]['ypos'])
        gmod.indxpara.gangdistscal = getattr(gmod.indxpara, 'gangdistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_expo(gang, maxmgang, paragenrscalfull[gmod.indxpara.gangdistscal]))) 
        lpri = -numbelem[l] * np.log(2. * pi) 
    if strgpdfn == 'tmpl':
        lpri = np.sum(lpdfspatprioobjt(dictelem[l]['ypos'], dictelem[l]['xpos'], grid=False))
    if strgpdfn == 'powr':
        lpri = retr_lpripowrdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'dpowslopbrek':
        lpri = retr_lpridpowdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, paragenrscalfull, l)
    if strgpdfn == 'dsrcexpo':
        lpri += -np.sum(np.sqrt((dictelem[l]['xpos'] - xpossour)**2 + (dictelem[l]['ypos'] - ypossour)**2) / \
                                                                getattr(gmod, 'dsrcdistsexppop%d' % l))
    if strgpdfn == 'tmpl':
        if strgpdfn.endswith('cons'):
            pdfnspatpriotemp = getattr(gmod, 'pdfnspatpriotemp')
            spatdistcons = paragenrscalfull[getattr(gmod.indxpara, 'spatdistcons')]
            lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
            lpdfspatpriointp = lpdfspatprioobjt(gdat.bctrpara.yposcart, gdat.bctrpara.xposcart)
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


def retr_pathoutpcnfg(pathbase, strgcnfg):
    
    pathoutpcnfg = pathbase + '/data/outp/' + strgcnfg + '/'
    
    return pathoutpcnfg


def proc_finl(gdat=None, strgcnfg=None, strgpdfn='post', listnamevarbproc=None, forcplot=False):
    
    gdatsimu = None
    
    print('proc_finl()')

    if strgcnfg is None:
        strgcnfg = gdat.strgcnfg
    
    # determine if the final-processing if nominal or tiling
    if isinstance(strgcnfg, list):
        liststrgcnfgmodi = strgcnfg
        strgcnfgfinl = tdpy.retr_strgtimestmp() + strgcnfg[0][15:] + 'tile'
        booltile = True
    else:
        liststrgcnfgmodi = [strgcnfg]
        strgcnfgfinl = strgcnfg
        booltile = False
    
    # determine of the gdatfinl object is available 
    boolgdatfinl = chec_statfile(pathbase, strgcnfgfinl, 'gdatfinlpost')
    boolgdatfinlgood = False
    if boolgdatfinl:
        print('Final-processing has been performed previously.')
        pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfgfinl)
        path = pathoutpcnfg + 'gdatfinl' + strgpdfn
        gdat = readfile(path) 
        boolgdatfinlgood = True

    if boolgdatfinl and boolgdatfinlgood:
        # read gdatfinl
        pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfgfinl)
        path = pathoutpcnfg + 'gdatfinl' + strgpdfn
        gdatfinl = readfile(path) 
        
        if gdatfinl.fitt.numbpopl > 0:
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.strgcnfgsimu is not None:
                        path = gdatfinl.pathoutpcnfgsimu + 'gdatfinlpost'
                        gdatsimu = readfile(path)
                    
    else:
        
        if booltile:
            gdatfinltile = tdpy.gdatstrt()
        
        indxstrgcnfggood = []
        liststrgtile = []
        liststrgcnfggood = []
        indxtiletemp = 0
        for n, strgcnfgmodi in enumerate(liststrgcnfgmodi):
            
            # read gdatinit
            boolgdatinit = chec_statfile(pathbase, strgcnfgmodi, 'gdatinit')
            if not boolgdatinit:
                if booltile:
                    print('Initial global object not found. Skipping...')
                    continue
                else:
                    print('Initial global object not found. Quitting...')
                    return
            
            pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfgmodi)
            path = pathoutpcnfg + 'gdatinit'
            
            gdatinit = readfile(path) 
            if booltile:
                gdatfinltile = gdatinit
                gdatfinl = gdatinit
            else:
                gdatfinl = gdatinit

            pathoutpcnfgmodi = retr_pathoutpcnfg(pathbase, strgcnfgmodi)
            listgdatmodi = []
            for k in gdatinit.indxproc:
                path = pathoutpcnfgmodi + 'gdatmodi%04d' % k + strgpdfn
                listgdatmodi.append(readfile(path))
            
            # erase
            gdatdictcopy = deepcopy(gdatinit.__dict__)
            for strg, valu in gdatdictcopy.items():
                if strg.startswith('fitt.indxpara.'):
                    delattr(gdatinit, strg)

            if gdatinit.boolsimuonly:
                print('Mock only run. Quitting final-processing...')
                return

            # read gdatmodi
            print('strgcnfgmodi')
            print(strgcnfgmodi)
            boolgdatmodi = chec_statfile(pathbase, strgcnfgmodi, 'gdatmodipost')
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
                gdatfinltile.pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfgfinl)
                numbsamptotlrsmp = gdatinit.numbsamptotl
                indxsamptotlrsmp = np.random.choice(gdatinit.indxsamptotl, size=gdatinit.numbsamptotl, replace=False)
            
            # aggregate samples from the chains
            if gdatinit.typeverb > 0:
                print('Reading gdatmodi objects from all processes...')
                timeinit = gdatinit.functime()
            
            if gdatinit.typeverb > 0:
                timefinl = gdatinit.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))
            
            if gdatinit.fitt.numbpopl > 0:
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
                    gdatfinl.gmrbstat = np.zeros((gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbdqlt))
                    for k in gdatfinl.fitt.indxpara.genr.base:
                        gdatfinl.gmrbparagenrscalbase[k] = tdpy.mcmc.gmrb_test(listparagenrscalfull[:, :, k])
                        if not np.isfinite(gdatfinl.gmrbparagenrscalbase[k]):
                            gdatfinl.gmrbparagenrscalbase[k] = 0.
                    listcntpmodl = getattr(gdatfinl, 'list' + strgpdfn + 'cntpmodl')
                    for i in gdatfinl.indxener:
                        for j in gdatfinl.indxpixl:
                            for m in gdatfinl.indxdqlt:
                                gdatfinl.gmrbstat[i, j, m] = tdpy.mcmc.gmrb_test(listcntpmodl[:, :, i, j, m])
                    if gdatfinl.typeverb > 0:
                        timefinl = gdatfinl.functime()
                        print('Done in %.3g seconds.' % (timefinl - timeinit))

                # calculate the autocorrelation of the chains
                if gdatfinl.typeverb > 0:
                    print('Computing the autocorrelation of the chains...')
                    timeinit = gdatfinl.functime()
                gdatfinl.atcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbdqlt, int(gdatfinl.numbparagenrfull / 2)))
                gdatfinl.timeatcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbdqlt))
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

                liststrgtile.append(strgcnfgmodi.split('_')[-2][-4:])
                liststrgcnfggood.append(strgcnfgmodi)
                indxstrgcnfggood.append(n)
                indxtiletemp += 1
                
                if len(liststrgtile) == 1:
                    for strgfeat in gdatfinl.refrgmod.namepara.genrelem:
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
                            
                #for strgfeat in gdatfinl.refrgmod.namepara.genrelem:
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
            gdatfinl.pathplotcnfg = gdatfinl.pathvisu + strgcnfgfinl + '/'
            make_fold(gdatfinl)
            indxstrgcnfggood = np.array(indxstrgcnfggood).astype(int)
            numbstrgcnfggood = indxstrgcnfggood.size
            numbtile = numbstrgcnfggood
            print('Found %d tiles with run tags:' % numbstrgcnfggood)
            for indxstrgcnfggoodtemp in indxstrgcnfggood:
                print(strgcnfg[indxstrgcnfggoodtemp])

            # np.concatenate reference elements from different tiles
            #for strgfeat in gdatfinl.refrgmod.namepara.genrelem:
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
            if gdatfinl.fitt.numbpopl > 0:
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
            #if gdatfinl.fitt.numbpopl > 0:
            #    gdatfinl.mlikindxelemfull = listindxelemfull[indxsamptotlmlik]
            gdatfinl.mlikparagenrscalbase = gdatfinl.mlikparagenrscalfull[gdatfinl.fitt.indxpara.genr.base]
            for k, namepara in enumerate(gdatfinl.fitt.namepara.genrbase):
                setattr(gdatfinl, 'mlik' + gmod.namepara.genr.base, gdatfinl.mlikparagenrscalbase[k])

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
        listparagenrscalbase = listparagenrscalfull[:, gdatfinl.fitt.indxpara.genr.base]
        for k, namepara in enumerate(gdatfinl.fitt.namepara.genrbase):
            setattr(gdatfinl, 'list' + strgpdfn + gmod.namepara.genr.base, listparagenrscalbase[:, k])
        setattr(gdatfinl, 'list' + strgpdfn + 'paragenrscalbase', listparagenrscalbase)

        if strgpdfn == 'post' and gdatfinl.checprio:
            pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
            path = pathoutpcnfg + 'gdatfinlprio'
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
            if gdatfinl.fitt.numbpopl > 0:
                if gdatfinl.boolbindspat:
                    histxposyposelemstkd = [[] for l in gdatfinl.fittindxpopl]
                
                    listxpos = getattr(gdatfinl, 'list' + strgpdfn + 'xpos')
                    listypos = getattr(gdatfinl, 'list' + strgpdfn + 'ypos')
                    for l in gdatfinl.fittindxpopl:
                        if gdatfinl.fitttypeelem[l] != 'lghtline':
                            histxposyposelemstkd[l] = np.zeros((gdatfinl.numbypospntsprob, gdatfinl.numbxpospntsprob, gdatfinl.numbbinsplot, numb))
                            temparry = np.concatenate([listxpos[n][l] for n in gdatfinl.indxsamptotl])
                            temp = np.empty((len(temparry), 3))
                            temp[:, 0] = temparry
                            temp[:, 1] = np.concatenate([listypos[n][l] for n in gdatfinl.indxsamptotl])
                            temp[:, 2] = np.concatenate([getattr(gdatfinl, 'list' + strgpdfn + strgfeat)[n][l] for n in gdatfinl.indxsamptotl])
                            blim = getattr(gdatfinl, 'blim' + strgfeat)
                            histxposyposelemstkd[l][:, :, :, k] = np.histogramdd(temp, \
                                                                    bins=(gdatfinl.binsxpospntsprob, gdatfinl.binsypospntsprob, bins))[0]
                    setattr(gdatfinl, strgpdfn + 'histxposyposelemstkd', histxposyposelemstkd)

            if gdatfinl.typeverb > 0:
                timefinl = gdatfinl.functime()
                print('Done in %.3g seconds.' % (timefinl - timeinit))

            ## construct a condensed catalog of elements
            if gdatfinl.boolcondcatl and gdatfinl.fitt.numbpopl > 0:
                
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
       
        if gdatfinl.fitt.numbpopl > 0 and strgpdfn == 'post':
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.strgcnfgsimu is not None:
                        path = gdatfinl.pathoutpcnfgsimu + 'gdatfinlpost'
                        gdatsimu = readfile(path)
                    
        # posterior corrections
        if gdatfinl.fitt.numbpopl > 0 and strgpdfn == 'post':

            ## perform corrections
            if gdatfinl.typedata == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:

                    for gmod.namepara.genr.elemvarbhist in gdatfinl.liststrgvarbhist:
                        strgvarb = gmod.namepara.genr.elemvarbhist[0]

                        if gmod.namepara.genr.elemvarbhist[1].startswith('aerr') or len(gmod.namepara.genr.elemvarbhist[2]) > 0 and gmod.namepara.genr.elemvarbhist[2].startswith('aerr'):
                            continue
                        if gmod.namepara.genr.elemvarbhist[1] == 'spec' or gmod.namepara.genr.elemvarbhist[1] == 'deflprof' or gmod.namepara.genr.elemvarbhist[1] == 'specplot':
                            continue
                        if len(gmod.namepara.genr.elemvarbhist[2]) > 0 and (gmod.namepara.genr.elemvarbhist[2] == 'spec' or \
                                    gmod.namepara.genr.elemvarbhist[2] == 'deflprof' or gmod.namepara.genr.elemvarbhist[2] == 'specplot'):
                            continue
                        
                        ## internal correction
                        listhist = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                        
                        for qq in gdatsimu.indxrefr:
                            l = int(gmod.namepara.genr.elemvarbhist[3][qq].split('pop')[1][0])
                            qq = int(gmod.namepara.genr.elemvarbhist[3][qq].split('pop')[2][0])
                            if gmod.namepara.genr.elemvarbhist[1][-4:] in gdatfinl.listnamerefr and \
                                    (len(gmod.namepara.genr.elemvarbhist[2]) == 0 or gmod.namepara.genr.elemvarbhist[2][-4:] in gdatfinl.listnamerefr):
                                listhistincr = listhist
                            else:
                                if gmod.namepara.genr.elemvarbhist[1][-4:] in gdatfinl.listnamerefr and len(gmod.namepara.genr.elemvarbhist[2]) > 0:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatsimu, 'listpostcmpl' + gmod.namepara.genr.elemvarbhist[2] + 'pop%dpop%d' % (l, qq))], 2)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatsimu, 'listpostfdis' + gmod.namepara.genr.elemvarbhist[2] + 'pop%dpop%d' % (qq, l))], 2)
                                elif len(gmod.namepara.genr.elemvarbhist[2][:-4]) > 0 and gmod.namepara.genr.elemvarbhist[2][-4:] in gdatfinl.listnamerefr:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatsimu, 'listpostcmpl' + gmod.namepara.genr.elemvarbhist[1] + 'pop%dpop%d' % (l, qq))], 1)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatsimu, 'listpostfdis' + gmod.namepara.genr.elemvarbhist[1] + 'pop%dpop%d' % (qq, l))], 1)
                                else:
                                    listcmpltrue = getattr(gdatsimu, 'listpostcmpl' + gmod.namepara.genr.elemvarbhist[3][qq])
                                    listfdistrue = getattr(gdatsimu, 'listpostfdis' + gmod.namepara.genr.elemvarbhist[3][qq])
                                if len(gmod.namepara.genr.elemvarbhist[2]) == 0:
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
                            
                            listgdatmodi[0].liststrgchan += ['incr' + gmod.namepara.genr.elemvarbhist[4][qq]]
                            setattr(gdatfinl, 'listpostincr' + gmod.namepara.genr.elemvarbhist[4][qq], listhistincr)
                        
                            ## external correction
                            for q in gdatfinl.indxrefr:
                                nametemp = gmod.namepara.genr.elemvarbhist[1] 
                                if len(gmod.namepara.genr.elemvarbhist[2]) > 0:
                                    nametemp += gmod.namepara.genr.elemvarbhist[2]
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

            if gdatfinl.fitt.numbpopl > 0 and strgchan in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
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
        if gmod.numbpopl > 0:
            listnametemp += ['lpripena']

        for namevarbscal in listnametemp:
            listtemp = getattr(gdatfinl, 'list' + strgpdfn + namevarbscal)
            minm = np.amin(listtemp)
            maxm = np.amax(listtemp)
            setattr(gdatfinl, 'minm' + namevarbscal, minm)
            setattr(gdatfinl, 'maxm' + namevarbscal, maxm)
            setattr(gdatfinl, 'scal' + namevarbscal, 'self')
            setp_varb(gdat, namevarbscal)
        
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
        path = gdatfinl.pathoutpcnfg + 'gdatfinl' + strgpdfn

        if gdatfinl.typeverb > 0:
            print('Writing gdatfinl to %s...' % path)
        writfile(gdatfinl, path) 
       
        filestat = open(gdatfinl.pathoutpcnfg + 'stat.txt', 'a')
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
    booltemp = chec_statfile(pathbase, strgcnfgfinl, 'plotfinl')
    if booltemp:
        print('Final plots already exist. Skipping...')
    else:
        if strgpdfn == 'post' and gdatfinl.checprio:
            path = pathoutpcnfg + 'gdatfinlprio'
            gdatprio = readfile(path)
        else:
            gdatprio = None
        
        if gdatfinl.boolmakeplot and getattr(gdatfinl, 'boolmakeplotfinl' + strgpdfn) or forcplot:
            plot_finl(gdatfinl, gdatprio=gdatprio, strgpdfn=strgpdfn, gdatsimu=gdatsimu, booltile=booltile)
            filestat = open(gdatfinl.pathoutpcnfg + 'stat.txt', 'a')
            filestat.write('plotfinl%s written.\n' % strgpdfn)
            filestat.close()


def retr_listgdat(liststrgcnfg, typegdat='finlpost'):
   
    listgdat = []
    for strgcnfg in liststrgcnfg:
        pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
        path = pathoutpcnfg + 'gdat%s' % typegdat
        listgdat.append(readfile(path))

    return listgdat


def make_fold(gdat):

    gdat.liststrgphas = ['fram', 'finl', 'anim']
        
    for strgpdfn in gdat.liststrgpdfn:
        setattr(gdat, 'path' + strgpdfn, gdat.pathplotcnfg + strgpdfn + '/') 
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
                    setattr(gdat, 'path' + strgphas + nameseco[:-1], gdat.pathplotcnfg + 'init/' + nameseco)
                else:
                    setattr(gdat, 'path' + strgpdfn + strgphas + nameseco[:-1], path + strgphas + '/' + nameseco)
    gdat.pathinfo = gdat.pathplotcnfg + 'info/'
    
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
        binsvarb = getattr(gdat.blimpara, name)
        deltvarb = getattr(gdat, 'delt' + name)
        pdfn[:, 0] = np.histogram(listvarb, blim=blimvarb)[0].astype(float)
        pdfn[:, 0] /= np.sum(pdfn[:, 0])
        pdfn[:, 0] /= deltvarb
    else:
        binsvarb = np.linspace(0, gmod.maxmpara.numbelemtotl, 51)
        
    if listvarb.ndim == 2:
        for k in range(listvarb.shape[1]):
            pdfn[:, k] = np.histogram(listvarb[:, k], blim=blimvarb)[0].astype(float)
            pdfn[:, k] /= np.sum(pdfn[:, k])
        pdfn *= 50.
    if listvarb.ndim == 3:
        for k in range(listvarb.shape[1]):
            for m in range(listvarb.shape[2]):
                pdfn[:, k, m] = np.histogram(listvarb[:, k, m], blim=blimvarb)[0].astype(float)
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
def chec_statfile(pathbase, strgcnfg, strggdat, typeverb=1):
    
    print('Checking the state file %s for %s...' % (strggdat, strgcnfg))
    
    pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
    
    # check the status file
    if not os.path.isfile(pathoutpcnfg + 'stat.txt'):
        if typeverb > 0:
            print('pathoutpcnfg')
            print(pathoutpcnfg)
            print('stat.txt not found.')
        return False

    # check the global object
    filestat = open(pathoutpcnfg + 'stat.txt', 'r')
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


def retr_los3(dlos, xpos, ypos):

    dglc = np.sqrt(8.5e3**2 + dlos**2 - 2. * dlos * 8.5e3 * np.cos(ypos) * np.cos(xpos))
    thet = np.arccos(np.sin(ypos) * dlos / dglc)
    phii = np.arcsin(np.sqrt(np.cos(ypos)**2 * dlos**2 + 8.5e3**2 - 2 * dlos * np.cos(ypos) * 8.5e3) / dglc)
    
    return dglc, thet, phii


def retr_glc3(dglc, thet, phii):

    xpos = dglc * np.sin(thet) * np.cos(phii)
    ypos = dglc * np.sin(thet) * np.sin(phii)
    zpos = dglc * np.cos(thet)
    dlos = np.sqrt(zpos**2 + xpos**2 + (8.5e3 - ypos)**2)
    xpos = np.arctan2(8.5e3 - ypos, xpos) - np.pi / 2
    ypos = np.arcsin(zpos / dlos)
   
    return dlos, xpos, ypos


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


def retr_dlosgalx(xpos, ypos, dglc):

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
    setp_varb(gdat, 'cntpdata', minm=minm, maxm=maxm, labl=['$C_{D}$', ''], scal='asnh', strgmodl='plot', cmap='Greys')
    
    maxm = maxmcntpdata
    minm = 1e-1 * minmcntpdata
    for strgmodl in gdat.liststrgmodl:
        gmod = getattr(gdat, strgmodl)
        setp_varb(gdat, 'cntpmodl', minm=minm, maxm=maxm, strgmodl=strgmodl, cmap='Greys')
    
    # residual limits
    maxm = np.ceil(maxmcntpdata * 0.1)
    minm = -np.ceil(maxmcntpdata * 0.1)
    setp_varb(gdat, 'cntpresi', minm=minm, maxm=maxm, labl=['$C_{R}$', ''], scal='asnh', strgmodl='plot')

    # 1-point function of the data counts
    for m in gdat.indxdqlt:
        if gdat.numbpixl > 1:
            for i in gdat.indxener: 
                histcntp = np.histogram(gdat.cntpdata[i, :, m], bins=gdat.blimpara.cntpdata)[0]
                setattr(gdat, 'histcntpdataen%02devt%d' % (i, m), histcntp)
        else:
            histcntp = np.histogram(gdat.cntpdata[:, 0, m], bins=gdat.blimpara.cntpdata)[0]
            setattr(gdat, 'histcntpdataevt%d' % m, histcntp)

    # obtain cartesian versions of the maps
    if gdat.typepixl == 'cart':
        ## data counts
        gdat.cntpdatacart = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbdqlt))
        gdat.cntpdatacart[:, gdat.indxpixlrofi, :] = gdat.cntpdata
        gdat.cntpdatacart = gdat.cntpdatacart.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbdqlt))
   

def retr_infodens(pdfnpost, pdfnprio):
    
    infodens = pdfnpost * np.log(pdfnpost / pdfnprio)

    return infodens


def retr_llik_unbd(gdat, strgmodl, cntpmodl):
    
    if gdat.liketype == 'pois':
        llik = np.exp(-(xpos - xpossamp)**2) + np.exp(-(yposgrid - ypossamp)**2)
    if gdat.liketype == 'gaus':
        llik = np.exp(-(xpos - xpossamp)**2) + np.exp(-(yposgrid - ypossamp)**2)
    
    return llik


def retr_llik_bind(gdat, strgmodl, cntpmodl):
   
    if gdat.liketype == 'pois':
        llik = gdat.cntpdata * np.log(cntpmodl) - cntpmodl
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.cntpdata - cntpmodl)**2 / gdat.varidata
     
    return llik


def retr_mapsgaus(gdat, xpos, ypos, spec, size, ellp, angl):
    
    rttrmatr = np.array([[np.cos(angl), -np.sin(angl)], [np.sin(angl), np.cos(angl)]])
    icovmatr = np.array([[1. / ((1. - ellp) * size)**2, 0.], [0., 1. / size**2]])

    posi = np.array([xposgrid - xpos, yposgrid - ypos])
    mapsgaus = flux * np.exp(-0.5 * np.sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / (1. - ellp)
        
    return mapsgaus


def retr_sbrtsers(gdat, xposgrid, yposgrid, xpos, ypos, spec, size, ellp, angl, seri=np.array([4.])):
   
    xposrttr = (1. - ellp) * (np.cos(angl) * (xposgrid - xpos) - np.sin(angl) * (yposgrid - ypos))
    yposrttr = np.sin(angl) * (xposgrid - xpos) + np.cos(angl) * (yposgrid - ypos) 
    angl = np.sqrt(xposrttr**2 + yposrttr**2)
    
    # interpolate pixel-convolved Sersic surface brightness
    if gdat.typesers == 'intp':

        shapinpt = angl.shape 
        inpt = np.empty(list(shapinpt) + [3])
        inpt[..., 0] = angl
        inpt[..., 1] = size
        inpt[..., 2] = seri
        
        sbrtsers = spec[:, None, None] * sp.interpolate.interpn((gdat.blimpara.xpossers, gdat.blimpara.halfsers, gdat.blimpara.indxsers), gdat.sersprof, inpt)[None, :, None]
    
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


def proc_anim(strgcnfg):
    
    pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
    
    print('Making animations of frame plots for %s...' % strgcnfg)
    
    path = pathoutpcnfg + 'gdatinit'
    gdat = readfile(path)
    for strgpdfn in gdat.liststrgpdfn:
        for nameextn in gdat.liststrgfoldanim:
            
            pathframextn = gdat.pathvisu + strgcnfg + '/' + strgpdfn + '/fram/' + nameextn
            pathanimextn = gdat.pathvisu + strgcnfg + '/' + strgpdfn + '/anim/' + nameextn
        
            try:
                listfile = fnmatch.filter(os.listdir(pathframextn), '*_swep*.%s' % gdat.typefileplot)
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
                
                strgtemp = '%s*_swep*.%s' % (name, gdat.typefileplot)
                listfile = fnmatch.filter(os.listdir(pathframextn), strgtemp)
                numbfile = len(listfile)
                liststrgextn = []
                for k in range(numbfile):
                    liststrgextn.append((listfile[k].split(name)[1]).split('_')[0])
                
                liststrgextn = list(set(liststrgextn))
                
                for k in range(len(liststrgextn)):
            
                    listfile = fnmatch.filter(os.listdir(pathframextn), name + liststrgextn[k] + '_swep*.%s' % gdat.typefileplot)
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
                        print('Run: %s, pdf: %s' % (strgcnfg, strgpdfn))
                        print('Making %s animation...' % name)
                        os.system(cmnd)
                    else:
                        print('GIF already exists.')
                        pass
    
    pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
    filestat = open(pathoutpcnfg + 'stat.txt', 'a')
    filestat.write('animfinl written.\n')
    filestat.close()
    

def plot_samp(gdat, gdatmodi, strgstat, strgmodl, strgphas, strgpdfn='post', gdatsimu=None, booltile=False):
    
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
    
    print('Plotting routine plot_samp() called...')

    if not booltile:
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxdqlt:
                    if not gdat.boolmakeframcent or gdat.boolmakeframcent and (i == int(gdat.numbener / 2) and m == int(gdat.numbdqlt / 2)):
                        ## data count maps
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m)
                        ## residual count maps
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpresi', i, m)
        
        if gdat.numbpixl > 1:
            if gmod.numbpopl > 0:
                if gmod.boolelemlens:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelem', booltdim=True)
        
        # temp -- restrict other plots to indxmodlelemcomp
        if gdat.boolbinsener:
            for specconvunit in gdat.listspecconvunit:
                if not gmod.boolbfun:
                    plot_sbrt(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, specconvunit)
    
            if gmod.boolapplpsfn:
                for l in gmod.indxpopl:
                    if gmod.typeelemspateval[l] == 'locl':
                        plot_psfn(gdat, gdatmodi, strgstat, strgmodl)
    
    if gmod.numbpopl > 0:
        # element parameter histograms
        if not (strgmodl == 'true' and gdat.typedata == 'inpt'):
            
            limtydat = gdat.limtydathistfeat

            for l in gmod.indxpopl:
                strgindxydat = 'pop%d' % l
                for nameparaderielemodim in gmod.namepara.derielemodim[l]:
                    if not (nameparaderielemodim == 'flux' or nameparaderielemodim == 'mcut' or \
                            nameparaderielemodim == 'deltllik' or nameparaderielemodim == 'defs' or nameparaderielemodim == 'nobj'):
                        continue
                                                                              
                    if gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt':
                        continue
                    indxydat = [l, slice(None)]
                    
                    name = nameparaderielemodim
                    namepopl = nameparaderielemodim + 'pop%d' % l
                    lablxdat = getattr(gmod.labltotlpara, namepopl)
                    scalxdat = getattr(gmod.scalpara, namepopl)
                    limtxdat = getattr(gmod.limtpara, namepopl)
                    meanxdat = getattr(gdat.bctrpara, name)
                        
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
                                lablydat = r'$\Sigma_{%s}$ [deg$^{-2}$]' % gmod.lablelemsubs[l]
                            if gdat.sdenunit == 'ster':
                                lablydat = r'$\Sigma_{%s}$ [sr$^{-2}$]' % gmod.lablelemsubs[l]
                        
                        ## plot the total number of elements
                        if ydattype == 'totl':
                            lablydat = r'$N_{%s}$' % gmod.lablelemsubs[l]
                    
                        if ydattype == 'totl' and not gdat.strgcnfgsimu is None:
                            listtypehist = ['hist', 'histcorrreca']
                        else:
                            listtypehist = ['hist']
                        
                        boolhistprio = not booltile
                        for typehist in listtypehist:
                            
                            if typehist == 'histcorrreca':
                                
                                if gmod.numbparaelem == 0 or gdat.factpriodoff == 0.:
                                    continue

                                if nameparaderielemodim == 'specplot' or nameparaderielemodim == 'spec' or nameparaderielemodim == 'deflprof':
                                    continue
                            
                                if not nameparaderielemodim in gmod.namepara.genrelem[l]:
                                    continue
                            
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'hist' + nameparaderielemodim + 'pop%d' % l, \
                                              'bctr' + nameparaderielemodim, scalydat='logt', lablxdat=lablxdat, \
                                              lablydat=lablydat, histodim=True, ydattype=ydattype, \
                                              scalxdat=scalxdat, meanxdat=meanxdat, limtydat=limtydat, \
                                              limtxdat=limtxdat, boolhistprio=boolhistprio, \
                                              #indxydat=indxydat, strgindxydat=strgindxydat, \
                                              nameinte='histodim/', typehist=typehist)
    
    if not booltile:
        if gmod.numbpopl > 0:
            # element parameter correlations
            for l in gmod.indxpopl:
                if strgmodl != 'true' and gdat.boolinforefr and gdat.boolasscrefr:
                    for strgfeat in gmod.namepara.derielemodim[l]:
                        if not (strgfeat == 'flux' or strgfeat == 'mass' or strgfeat == 'deltllik' or strgfeat == 'nobj') and \
                                                                                    (gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
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
                    
        if not (gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # plots
            for i in gdat.indxener:
                for m in gdat.indxdqlt:
                    if gmod.numbpopl > 1:
                        if gmod.numbpopl > 0:
                            for l in gmod.indxpopl:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m, indxpoplplot=l)
        
            ## histograms of the number of counts per pixel
            limtxdat = [gdat.minmpara.cntpmodl, gdat.maxmpara.cntpmodl]
            for nameecom in gmod.listnameecomtotl:
                name = 'histcntp' + nameecom
                for m in gdat.indxdqlt: 
                    for i in gdat.indxener:
                        if gdat.numbener > 1:
                            name += 'en%02d' % (i)
                        if gdat.numbdqlt > 1:
                            name += 'evt%d' % (m)
                            
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                             name, 'bctrcntpdata', scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, histodim=True, \
                                             lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbener], limtxdat=limtxdat)

            
            ## highest amplitude element
            # temp
            if gmod.numbpopl > 0:
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
                                                         not (gdat.typedata == 'simu' and (strgfeat.endswith('pars') or strgfeat.endswith('nrel'))):
                                        
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgclas + strgfeat + strgindxydat, \
                                                  'bctr' + strgfeat, lablxdat=lablxdat, \
                                                  lablydat=lablydat, \
                                                  #plottype='errr', \
                                                  scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                  omittrue=True, nameinte=nameinte)

            if gmod.numbpopl > 0:
                alph = 0.1
                if strgmodl == 'true':
                    pathtemp = gdat.pathinit
                else:
                    if strgstat == 'this':
                        pathtemp = gdat.pathplotcnfg + strgpdfn + '/fram/'
                    elif strgstat == 'mlik':
                        pathtemp = gdat.pathplotcnfg + strgpdfn + '/finl/'
                    elif strgstat == 'pdfn':
                        pathtemp = gdat.pathplotcnfg + strgpdfn + '/finl/'
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
                                listxdat.append(gdat.bctrpara.enerplot)
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
                                        specplottemp *= gdat.bctrpara.enerplot
                                    if specconvunit[0] == 'en02':
                                        specplottemp *= gdat.bctrpara.enerplot**2
                                    if specconvunit[0] == 'en03':
                                        # temp
                                        pass
                                    listydat.append(specplottemp)
                                
                                lablydat = getattr(gmod.lablpara, 'flux' + specconvunit[0] + specconvunit[1] + 'totl')
                                strgtemp = specconvunit[0] + specconvunit[1]
                                if specconvunit[0] == 'en03':
                                    strgtemp += specconvunit[2]
                                path = pathtemp + strgstat + 'specpop%d%s%s.%s' % (l,  strgtemp, strgswep, gdat.typefileplot)
                                limtydat = [np.amin(gdat.minmspec), np.amax(gdat.maxmspec)]
                                tdpy.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', \
                                                           lablxdat=gdat.lablenertotl, colr=colr, alph=alph, \
                                                           plottype=listplottype, limtxdat=[gdat.minmener, gdat.maxmener], lablydat=lablydat, \
                                                           limtydat=limtydat)
                    
                    if gmod.boollenssubh:

                        ## deflection profiles
                        lablxdat = gdat.labltotlpara.gang
                        if strgstat == 'pdfn':
                            deflprof = [np.empty((gdat.numbanglfull, gdat.numbstkscond))]
                            asca = [np.empty(gdat.numbstkscond)]
                            acut = [np.empty(gdat.numbstkscond)]
                            for r in gdat.indxstkscond:
                                deflprof[0][:, r] = gdat.dictglob['poststkscond'][r]['deflprof'][0, :]
                                asca[0][r] = gdat.dictglob['poststkscond'][r]['asca'][0]
                                acut[0][r] = gdat.dictglob['poststkscond'][r]['acut'][0]
                        
                        for l in gmod.indxpopl:
                            if strgmodl == 'true':
                                deflprof = gdat.true.this.dictelem[l]['deflprof']
                            else:
                                deflprof = gmodstat.dictelem[l]['deflprof']

                            xdat = gdat.bctrpara.anglfull * gdat.anglfact
                            listydat = []
                            listvlinfrst = []
                            listvlinseco = []
                            
                            if 'deflprof' in gmod.typeelem[l]:

                                if strgmodl == 'true':
                                    deflproftemp = deflprof[l][0, :, k]
                                else:
                                    deflproftemp = deflprof[l][:, k]
                                
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
                                
                                path = pathtemp + strgstat + 'deflsubhpop%d%s.%s' % (l, strgswep, gdat.typefileplot)
                                limtydat = [1e-3, 1.]
                                limtxdat = [1e-3, 1.]
                                tdpy.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', \
                                                                    lablxdat=lablxdat, drawdiag=True, limtydat=limtydat, \
                                                                    limtxdat=limtxdat, colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', \
                                                                    listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
                        
                if gdat.typedata == 'simu':
                    # pulsar masses
                    for l in gmod.indxpopl:
                        if gmod.typeelem[l] == 'lghtpntspuls':
                            lablxdat = gdat.labltotlpara.gang
                            limtydat = [gdat.minmmassshel, gdat.maxmmassshel]
                            lablydat = gdat.lablmassshel
                            name = 'massshelpop%d' % l
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctranglhalf', scalydat='logt', \
                                                                           lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                    if gmod.boollens:
                        ## radial mass budget
                        lablxdat = gdat.labltotlpara.anglfromhost
                        for strgcalcmasssubh in gdat.liststrgcalcmasssubh:
                            
                            # host mass
                            for e in gmod.indxsersfgrd:
                                strgsersfgrd = 'isf%d' % e
                                
                                limtydat = [gdat.minmmcut, getattr(gdat.maxmpara, 'masshost%s%s' % (strgsersfgrd, strgcalcmasssubh))]
                                lablydat = getattr(gmod.labltotlpara, 'masshost' + strgsersfgrd + strgcalcmasssubh)
                                name = 'masshost%s%s' % (strgsersfgrd, strgcalcmasssubh)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctranglhalf', scalydat='logt', \
                                                                          lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)
                            
                                if gmod.boolelemdeflsubhanyy:
                                    # subhalo masses
                                    limtydat = [gdat.minmmcut, getattr(gdat.maxmpara, 'masssubh%s' % strgsersfgrd + strgcalcmasssubh)]
                                    lablydat = getattr(gmod.labltotlpara, 'masssubh%s' % strgsersfgrd + strgcalcmasssubh)
                                    name = 'masssubh%s' % (strgcalcmasssubh)
                                    plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctranglhalf', scalydat='logt', \
                                                                              lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                                    # subhalo mass fraction
                                    limtydat = [1e-3, 0.1]
                                    lablydat = getattr(gmod.labltotlpara, 'fracsubh%s' % strgsersfgrd + strgcalcmasssubh)
                                    name = 'fracsubh%s' % (strgcalcmasssubh)
                                    plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctranglhalf', scalydat='logt', \
                                                                lablxdat=lablxdat, lablydat=lablydat, limtydat=limtydat)

                alph = 0.1

                if gdat.boolmodipsfn and gmod.boolelempsfnanyy:
                    ## PSF radial profile
                    for i in gdat.indxener:
                        for m in gdat.indxdqlt:
                            indxydat = [i, slice(None), m]
                            strgindxydat = 'en%02devt%d' % (i, m)
                            lablxdat = gdat.labltotlpara.gang
                            limtydat= np.array([1e-3, 1e3]) * gdat.anglfact**2
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psfn', \
                                                            'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalydat='logt', \
                                                            lablxdat=lablxdat, lablydat=r'$\mathcal{P}$', limtydat=limtydat)
        
                # internally and externally corrected element parameter histograms
                if gdat.typedata == 'inpt' and strgstat == 'pdfn' and gdat.strgcnfgsimu is not None:
                    limtydat = gdat.limtydathistfeat
                    for l in gmod.indxpopl:
                        strgindxydat = 'pop%d' % l
                        for strgfeat in gmod.namepara.derielemodim[l]:
                            if strgfeat.startswith('aerr') or strgfeat == 'specplot' or strgfeat == 'spec' or strgfeat == 'deflprof':
                                continue
                            lablydat = r'$N_{%s}$' % gmod.lablelemsubs[l]
                            for namecorr in ['incr', 'excr']:
                                nameinte = namecorr + 'odim/'
                                for qq in gdatsimu.indxrefr:
                                    if namecorr == 'excr':
                                        if not strgfeat in gmod.namepara.extrelem[l]:
                                            continue
                                        q = gdat.listnamerefr.index(strgfeat[-4:])
                                        if getattr(gdat, 'crex' + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)) is None:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctr' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                          lablydat=lablydat, histodim=True, ydattype='totl', \
                                                          scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                          nameinte=nameinte)
               
                                    else:
                                        if strgfeat in gmod.namepara.extrelem[l]:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%d' % (qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'bctr' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                      lablydat=lablydat, histodim=True, ydattype='totl', \
                                                      scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                      nameinte=nameinte)


    if not (gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
        if gmod.numbpopl > 0:
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
                                        for qq in gdatsimu.indxrefr:
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
        
            if not (gdat.typedata == 'simu' and (gmod.maxmpara.numbelemtotl == 0 or gmod.maxmpara.numbelemtotl == 0)):
                

                for q in gdat.indxrefr:
                    
                    if strgphas == 'init' and gdat.typedata == 'simu':
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
        if not (gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # data and model count scatter
            for m in gdat.indxdqltplot:
                if gdat.numbpixl > 1:
                    for i in gdat.indxener:
                        plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m, indxenerplot=i)
                else:
                    plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m)

            ## spatial priors
            # temp
            if gdat.numbpixl > 1:
                if gmod.numbpopl > 0:
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
                    for m in gdat.indxdqlt:
                        for c in gmod.indxback:
                            if gmod.boolbfun:
                                continue
                            if not gmod.boolunifback[c]:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpback%04d' % c, i, m, strgcbar='cntpdata')
                
                ## count error
                if strgmodl != 'true':
                    if gmod.numbpopl > 0:
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
                    for m in gdat.indxdqlt:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpmodl', i, m, strgcbar='cntpdata')
            
                # likelihood
                if strgmodl != 'true':
                    for i in gdat.indxener:
                        for m in gdat.indxdqlt:
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
            
            if gdat.boolpenalpridiff:
                for i in gdat.indxener:
                    for m in gdat.indxdqlt:
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                                'psecodimdatapntsen%02devt%d' % (i, m), 'bctrmpolodim', lablxdat='$l$', lablydat='$P_{resi}(l)$', \
                                                                                                                 limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psecodimdatapntsprioen%02devt%d' % (i, m), 'bctrmpolodim', lablxdat='$l$', \
                                                                                           lablydat='$P_{prio}(l)$', limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                
            if gmod.boollens:
                indxydat = [slice(None)]
                strgindxydat = ''
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecodim', 'bctrwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P(k)$', limtydat=[1e-1, 1e2], \
                                                                          scalxdat='logt', scalydat='logt', indxydat=indxydat, strgindxydat=strgindxydat)
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdefl', 'bctrdefl', \
                                                                        scal='self', lablxdat=r'$\alpha$ [arcsec]', lablydat=r'$N_{pix}$', \
                                                                                 strgindxydat=strgindxydat, indxydat=indxydat, histodim=True)
            if gmod.numbpopl > 0 and gmod.boolelemdeflsubhanyy:
                indxydat = [slice(None)]
                strgindxydat = ''
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecelemodim', 'bctrwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P_{sub}(k)$', \
                                       strgindxydat=strgindxydat, indxydat=indxydat, limtydat=[1e-5, 1e-1], scalxdat='logt', scalydat='logt')
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdeflsubh', 'bctrdeflsubh', scal='self', lablxdat=r'$\alpha$ [arcsec]', \
                                       strgindxydat=strgindxydat, indxydat=indxydat, lablydat=r'$N_{pix}$', histodim=True)
            
            if gmod.boollens:
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrd', i, -1, strgcbar='cntpdata')
                    if gmod.numbpopl > 0 and gmod.boolelemsbrtextsbgrdanyy:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdgalx', i, -1, strgcbar='cntpdata')
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdexts', i, -1, strgcbar='cntpdata')
                
                # gradient of the lens emission
                for i in gdat.indxener:
                    for m in gdat.indxdqlt:
                        plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplensgrad', indxenerplot=i, indxdqltplot=m)
                
        if not (gdat.boolmakeshrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            if gmod.boollens:
                # overall deflection field
                plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, multfact=0.1)
                
                # deflection field due to individual lenses
                for k in range(gmod.numbdeflsingplot):  
                    if k == 0:
                        multfact = 0.1
                    elif k == 1:
                        multfact = 1.
                    elif k >= 2:
                        multfact = 10.
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxdefl=k, multfact=multfact)
                
                # residual deflection field
                if strgmodl == 'fitt' and gdat.typedata == 'simu':
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, nameparagenrelem='resi', multfact=100.)
                    if strgstat != 'pdfn':
                        for k in range(numbsingcomm):
                            plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, nameparagenrelem='resi', indxdefl=k, multfact=100.)
                    
                    if gdat.numbpixl > 1:
                        if gmod.numbpopl > 0:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresi', booltdim=True)
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresiperc', booltdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresi', booltdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresiperc', booltdim=True)
    

def delete_strgcnfg(strgcnfg):
    
    pathdata = pathbase + 'data/outp/'
    pathvisu = pathbase + 'visuals/'
    
    cmnd = 'rm -rf %s%s' % (pathdata, strgcnfg)
    print(cmnd)
    os.system(cmnd)
    cmnd = 'rm -rf %s%s' % (pathvisu, strgcnfg)
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
        print('Writing to %s...' % path)
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
        print('Writing to %s...' % pathpvkstdim)
        plt.savefig(pathpvkstdim)
        plt.close(figr)

    elif name != namefull:
        
        lablydat = '$D_{KL}$'
        lablxdat = getattr(gmod.lablpara, name + 'totl')
        xdat = getattr(gdat, 'bctr' + name)
        ydat = getattr(gdat, 'info' + namefull)
        tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, scal)
        
        ydat = getattr(gdat, 'pvks' + namefull)
        pathpvks = gdat.pathinfo + 'pvks' + namefull
        tdpy.mcmc.plot_plot(pathpvks, xdat, ydat, lablxdat, '$p_{KS}$', scal)
        
    else:
        # horizontal axis
        xdat = getattr(gdat, 'bctr' + name)
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


def plot_finl(gdat=None, gdatprio=None, strgcnfg=None, strgpdfn='post', gdatsimu=None, booltile=None):
    
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
            pathfinlprop = getattr(gdat, 'path' + strgpdfn + 'finl%s' % gdat.nameproptype[n])
            for k in gdat.indxtermlacp:
                varb = getattr(gdat, 'list' + strgpdfn + gdat.listnametermlacp[k])
                labl = gdat.listlabltermlacp[k]
                
                if listindxsamptotlproptotl[n].size > 0 and (varb[listindxsamptotlproptotl[n]] != 0.).any():
                    path = pathfinlprop + gdat.listnametermlacp[k] + 'totl'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlproptotl[n]], labl, titl=gdat.nameproptype[n] + ', Total')
                
                if listindxsamptotlpropaccp[n].size > 0 and (varb[listindxsamptotlpropaccp[n]] != 0.).any():
                    path = pathfinlprop + gdat.listnametermlacp[k] + 'accp'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropaccp[n]], labl, titl=gdat.nameproptype[n] + ', Accepted')
                
                if listindxsamptotlpropreje[n].size > 0 and (varb[listindxsamptotlpropreje[n]] != 0.).any():
                    path = pathfinlprop + gdat.listnametermlacp[k] + 'reje'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropreje[n]], labl, titl=gdat.nameproptype[n] + ', Rejected')
            
        if gdat.checprio and strgpdfn == 'post' and not booltile:
            # this works only for scalar variables -- needs to be generalized to all variables
            if gdatprio is None:
                pathoutpcnfg = retr_pathoutpcnfg(pathbase, strgcnfg)
                path = pathoutpcnfg + 'gdatfinlprio'
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
                blim = np.linspace(minm, maxm, 40)
                axis.hist(gdat.gmrbstat.flatten(), blim=blim, label='Data proj.')
                axis.hist(gdat.gmrbparagenrscalbase, blim=blim, label='Fixed dim.')
                axis.set_xlabel('PSRF')
                axis.set_ylabel('$N_{stat}$')
                plt.tight_layout()
                path = pathdiag + 'gmrbhist.%s' % gdat.typefileplot
                print('Writing to %s...' % path)
                figr.savefig(path)
                plt.close(figr)
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.plot(gmod.indxpara.genrbase, gdat.gmrbparagenrscalbase)
                axis.set_xticklabels(gmod.labltotlpara.genr.base)
                axis.set_ylabel('PSRF')
                plt.tight_layout()
                path = pathdiag + 'gmrbparagenrscalbase.%s' % gdat.typefileplot
                print('Writing to %s...' % path)
                figr.savefig(path)
                plt.close(figr)
                
                for i in gdat.indxener:
                    for m in gdat.indxdqlt:
                        maps = gdat.gmrbstat[i, :, m]
                        path = pathdiag + 'gmrbdataen%02devt%d.%s' % (i, m, gdat.typefileplot)
                        tdpy.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, typepixl=gdat.typepixl, \
                                                                   minmxpos=gdat.anglfact*gdat.minmxpos, maxmxpos=gdat.anglfact*gdat.maxmxpos, \
                                                                   minmypos=gdat.anglfact*gdat.minmypos, maxmypos=gdat.anglfact*gdat.maxmypos)
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
            histtotl = axis.hist(listindxsamptotlproptotl[n], blim=blimtimemcmc)[0]
            histaccp = axis.hist(listindxsamptotlpropaccp[n], blim=blimtimemcmc)[0]
            axis.set_ylabel('%s' % gdat.nameproptype[n])
            if k == gdat.numbproptype - 1:
                axis.set_xlabel('$i_{samp}$')
        plt.tight_layout()
        path = pathdiag + 'accpratiproptype.%s' % gdat.typefileplot
        print('Writing to %s...' % path)
        figr.savefig(path)
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
        #path = pathdiag + 'chro.%s' % (gdat.listnamechro[k], gdat.typefileplot)
        #print('Writing to %s...' % path)
        #figr.savefig(path)
        #plt.close(figr)

    # temp
    gdat.lablpmea = 'Mean'

    # posterior versions of the frame plots
    plot_samp(gdat, None, 'pdfn', 'fitt', 'finl', strgpdfn=strgpdfn, gdatsimu=gdatsimu, booltile=booltile)
   
    if booltile:
        return

    if gmod.numbpopl > 0:
        if gdat.typeverb > 0:
            print('A mosaic of samples...')
    
        ## mosaic of images of posterior catalogs
        if gdat.numbpixl > 1:
            plot_mosa(gdat, strgpdfn)
    
    ## randomly selected trandimensional parameters
    if gmod.numbpopl > 0:
        if gdat.typeverb > 0:
            print('Transdimensional parameters...')
    
        # choose the parameters based on persistence
        stdvlistsamptran = np.std(listparagenrscalfull[:, gmod.indxsamptrap], axis=0)
        indxtrapgood = np.where(stdvlistsamptran > 0.)[0]
        gmod.numbpara.totl.elemgood = indxtrapgood.size
        gmod.numbpara.totl.elemplot = min(3, gmod.numbpara.totl.elemgood)
        if gmod.numbpara.totl.elemplot > 0:
            indxtrapplot = np.sort(np.random.choice(gmod.indxsamptrap[indxtrapgood], size=gmod.numbpara.totl.elemplot, replace=False))

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
    tdpy.mcmc.plot_grid(path, 'paragenrscalbase', listparagenrscalbase, gmod.labltotlpara.genr.basetotl, truepara=truepara, listvarbdraw=[mlikpara])
    
    # stacked posteiors binned in position and flux
    if gmod.numbpopl > 0 and gdat.numbpixl > 1:
        liststrgbins = ['quad', 'full']
        for l in gmod.indxpopl:
            plot_histxposyposelemstkd(gdat, strgpdfn, l, 'cumu')
            for strgbins in liststrgbins:
                plot_histxposyposelemstkd(gdat, strgpdfn, l, strgbins, namepara.elemsign[l])

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
        if gdat.typedata == 'simu':
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
    path = pathdiag + 'memoresi.%s' % gdat.typefileplot
    print('Writing to %s...' % path)
    figr.savefig(path)
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
            if gdat.typedata == 'simu':
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
                
            numbplot = gdat.numblablsbrtspec + numbelemtemp
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
            
            if gmod.numbpopl > 0 and gmod.boolelemsbrtdfncanyy and not (liststrgmodl[a] == 'true' and gdat.refr.numbelemtotl == 0):
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
            if gdat.numblablsbrt > 1:
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
                    factener = gdat.bctrpara.ener
                if specconvunit[0] == 'en02':
                    factener = gdat.bctrpara.ener**2
                if specconvunit[0] == 'en03':
                    # temp
                    pass
                    factener = 1.
                    #indxenerintv = np.where((gdat.bctrpara.ener < specconvunit[4]) & (gdat.bctrpara.ener > specconvunit[3]))[0]
                    #ener = np.concatenate((np.array([specconvunit[3]]), gdat.bctrpara.ener[indxenerintv], np.array([specconvunit[4]])))
                    #
                    #for k in range(3):
                    #    if k == 0:
                    #        ydattemp = 
                    #    ydatminmener = np.interp(specconvunit[3], gdat.bctrpara.ener, ydat)
                    #    ydatmaxmener = np.interp(specconvunit[4], gdat.bctrpara.ener, ydat)
                    #    ydat = np.concatenate((np.array([ydatminmener]), ydat[indxenerintv], np.array([ydatmaxmener])))
                    #    ydat = np.trapz(ydat, gdat.bctrpara.ener)
                    #
                    #yerrminmener = np.interp(specconvunit[3], gdat.bctrpara.ener, yerr, axis=1)
                    #yerrmaxmener = np.interp(specconvunit[4], gdat.bctrpara.ener, yerr, axis=1)
                    #ydat = np.stack((np.array([yerrminmener]), ydat[indxenerintv], np.array([yerrmaxmener])))
                    #
                    #
                    #yerr = np.trapz(yerr, gdat.bctrpara.ener)


                xdat = gdat.bctrpara.ener
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
            axis.set_xlim([np.amin(gdat.blimpara.ener), np.amax(gdat.blimpara.ener)])
            
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
            print('Writing to %s...' % path)
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
    
    blim = np.linspace(minm, 2. * maxm, 2 * numbvarb + 1)
    
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
        axis.hist(arry[k, :], blim=blim, alpha=alph, label='$f_1$ (Sampled)', color='b')
        axis.hist(totl, blim=blim, alpha=alph, label='$f_0$ (Sampled)', color='g')
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
        pathfold = os.environ["TDGU_DATA_PATH"] + 'visuals/powrpdfn/'
        path = pathfold + 'powrpdfn%04d.%s' % (n, gdat.typefileplot)
        print('Writing to %s...' % path)
        figr.savefig(path)
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
                varbfrst = gdat.blimpara.strgfrst
                varbseco = getattr(gdat.blimpara, strgseco)
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
                meanfrst = getattr(gdat.blimpara, strgfrst)
                meanseco = getattr(gdat.blimpara, strgseco)
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
        xdat = getattr(gdat.bctrpara, strgxdat[4:])
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
                binsxdat = getattr(gdat.blimpara, strgxdat[4:])
            else:
                deltxdat = getattr(gdat.deltpara, strgxdat[4:])
                binsxdat = getattr(gdat.blimpara, strgxdat[4:])

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
                           
            if histodim:
                if typehist == 'histcorrreca':
                    reca = getattr(gdat.true, 'reca' + strgydat[4:])
                axis.plot(xdattemp, 10. * reca, color='purple', label='PTFN', alpha=gdat.alphline)

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
                if gdat.typedata == 'simu' and not omittrue:
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
    
    #print('Legend failed when')
    #print('strgstat')
    #print(strgstat)
    #print('strgmodl')
    #print(strgmodl)
    #print('strgydat')
    #print(strgydat)
    #raise Exception('')
    
    make_legd(axis, offs=offslegd, ptch=ptch, line=line)

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


def plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxdqltplot, indxenerplot=None):
    
    gmod = getattr(gdat, strgmodl)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn)
    if indxenerplot is None:
        xdat = gdat.cntpdata[:, :, indxdqltplot].flatten()
        ydat = ydat[:, :, indxdqltplot].flatten()
        nameplot = 'scatcntpevt%d' % (indxdqltplot)
        if strgstat == 'pdfn':
            indxvarb = [slice(None), slice(None), indxdqltplot]
    else:
        xdat = gdat.cntpdata[indxenerplot, :, indxdqltplot]
        ydat = ydat[indxenerplot, :, indxdqltplot]
        nameplot = 'scatcntpen%02devt%d' % (indxenerplot, indxdqltplot)
        if strgstat == 'pdfn':
            indxvarb = [indxenerplot, slice(None), indxdqltplot]
    if strgstat == 'pdfn':
        yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn, strgmome='errr', indxvarb=indxvarb)
    colr = gmod.colr

    if strgstat == 'pdfn':
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color=gmod.colr, capsize=5)
    else:
        axis.plot(xdat, ydat, marker='o', ls='', markersize=5, color=gmod.colr)
    gdat.limtcntpdata = [gdat.blimpara.cntpdata[0], gdat.blimpara.cntpdata[-1]]
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
    blim = np.empty((gdat.numbprox, numbbins + 1))
    indxpixlproxsize = np.empty((numbfluxprox, gdat.numbpixlfull))
    for h in gdat.indxprox:
        for j in gdat.indxpixlfull:
            indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
        bins[h, :] = np.logspace(np.log10(np.amin(indxpixlproxsize[h, :])), np.log10(np.amax(indxpixlproxsize[h, :])), numbbins + 1)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in gdat.indxprox:
        axis.hist(indxpixlproxsize[h, :], blim=blim[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphhist)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixlfull, label='ROI', ls='--')
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    make_legd(axis)
    plt.tight_layout()
    figr.savefig(gdat.pathplotcnfg + 'init/indxprox.%s' % gdat.typefileplot)
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
    figr.savefig(gdat.pathplotcnfg + 'evidtest.%s' % gdat.typefileplot)
    plt.close(figr)
    
    
def plot_histxposyposelemstkd(gdat, strgpdfn, indxpoplplot, strgbins, strgfeat=None):
    
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
    
    histxposyposelemstkd = getattr(gdat, strgpdfn + 'histxposyposelemstkd')

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
                temp = np.sum(histxposyposelemstkd[indxpoplplot][:, :, indxlowr:indxuppr], 2).T
            else:
                temp = np.sum(np.sum(histxposyposelemstkd[indxpoplplot], 2), 2).T
                
            if np.where(temp > 0.)[0].size > 0:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', \
                                                            extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
            else:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                
            if strgfeat is not None:
                blim = getattr(gdat.blimpara, strgfeat)
            
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
                mrkrsize = retr_mrkrsize(gdat, strgmodl, reframpl[q][0, indxelem], gdat.refr.nameparagenrelemampl[q])

                if indxelem.size > 0:
                    axis.scatter(gdat.anglfact * gdat.refr.dictelem[q]['xpos'][0, indxelem], gdat.anglfact * gdat.refr.dictelem[q]['ypos'][0, indxelem], \
                                                s=mrkrsize, alpha=gdat.alphelem, marker=gdat.refrlistmrkrhits[q], lw=2, color=gdat.refr.colrelem[q])

            if a == numbrows - 1:
                axis.set_xlabel(gdat.lablxpostotl)
            else:
                axis.set_xticklabels([])
            if b == 0:
                axis.set_ylabel(gdat.lablypostotl)
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
    path = getattr(gdat, 'path' + strgpdfn + 'finl') + 'histxposyposelemstkd%s%spop%d' % (strgbins, strgtemp, indxpoplplot) + '.%s' % gdat.typefileplot
    figr.savefig(path)
    plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.blimpara.angl)

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
    figr.savefig(gdat.pathplotcnfg + 'king.%s' % gdat.typefileplot)
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
        figr.savefig(gdat.pathvisu + 'talkintr.%s' % gdat.typefileplot, facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_psfn(gdat, gdatmodi, strgstat, strgmodl):

    gmod = getattr(gdat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgmodl)
    gmodstat = getattr(gdatobjt, strgstat)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for i in gdat.indxener:
        for m in gdat.indxdqlt:
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
                axis.plot(gdat.blimpara.angl * gdat.anglfact, gdat.blimpara.prox[k] * gmodstat.psfn[i, :, m], label=labl, color=colr, alpha=alph)
                axis.set_xlim([np.amin(gdat.blimpara.angl) * gdat.anglfact, np.amax(gdat.blimpara.angl) * gdat.anglfact])
                if k > 0:
                    axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
            axis.set_yscale('log')
            axis.set_xlabel(gdat.labltotlpara.gang)
            axis.set_ylabel(gdat.lablsbrttotl)
            
            limt = gdat.specfraceval * np.amax(gdat.blimpara.prox[0] * gmodstat.psfn[i, :, m])

            if limt != 0.:
                axis.axhline(limt, color='red', ls=':', label='Flux floor')
            
            make_legd(axis)

            plt.tight_layout()
            name = 'psfn'
            if gdat.numbener > 1:
                name += 'en%02d' % i
            if gdat.numbdqlt > 1:
                name += 'evt%d' % m
            figr.savefig(gdat.pathinit + name + '.%s' % gdat.typefileplot)
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
                for m in gdat.indxdqltplot:
                    
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            
                            n = indxsampmosa[numbcols*a+b]
                            gdatmodi.this.paragenrscalfull = listparagenrscalfull[n, :].flatten()
                            gdatmodi.this.paragenrunitfull = listparagenrunitfull[n, :].flatten()
                            if gmod.numbpopl > 0:
                                gdatmodi.this.indxelemfull = getattr(gdat, 'list' + strgpdfn + 'indxelemfull')[n]
                                proc_samp(gdat, gdatmodi, 'this', 'fitt')

                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.lablxpostotl)
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.lablypostotl)
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
                        path = pathfinl + 'mosa' + strg + 'en%02dA.%s' % (gdat.indxenerincl[i], gdat.typefileplot)
                    else:
                        path = pathfinl + 'mosa' + strg + 'en%02ddqlt%d.%s' % (gdat.indxenerincl[i], gdat.indxdqltincl[m], gdat.typefileplot)
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
                             ('xpos','modl'), \
                             ('ypos','modl'), \
                             ('numbelem','xpos'), \
                             ('numbelem','ypos'), \
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
                             ('spatdistcons', 'xpos'), \
                             ('spatdistcons', 'ypos'), \
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
    labl['xpos'] = r'$\vec{\theta_1}$'
    labl['ypos'] = r'$\vec{\theta_2}$'
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
    posi['xpos'] = np.array([-0.3, -0.0])
    posi['ypos'] = np.array([-0.1, -0.0])
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
            print('%15s %15s %15s' % (grap.edges()[k][0], grap.edges()[k][1], listcolr[k]))

    size = 1000
    nx.draw(grap, posi, labels=labl, ax=axis, edgelist=[], nodelist=[])
    nx.draw_networkx_edges(grap, posi, ax=axis, labels=labl, edge_color=listcolr)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['modl', 'data'], node_color='grey', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbelem'], node_color='b', node_size=size)
    if plottype.startswith('lght'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'amplslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['xpos', 'ypos', 'ampl', 'sind'], node_color='g', node_size=size)
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
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['xpos', 'ypos', 'defs'], node_color='g', node_size=size)
    if plottype == 'lens0001':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['xpos', 'ypos', 'defs', 'asca', 'acut'], node_color='g', node_size=size)
    
    pathplot = pathbase + 'visuals/'
    plt.tight_layout()
    figr.savefig(pathplot + 'grap%s.%s' % (plottype, gdat.typefileplot))
    plt.close(figr)


def plot_3fgl_thrs(gdat):

    path = pathbase + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = astropy.io.fits.getdata(path, 0)

    yposfgl3 = np.linspace(-90., 90., 481)
    xposfgl3 = np.linspace(-180., 180., 960)

    yposexpo = np.linspace(-90., 90., 400)
    xposexpo = np.linspace(-180., 180., 800)

    #fluxthrs = interp2d(xposfgl3, yposfgl3, fluxthrs)(xposexpo, yposexpo)
    fluxthrs = griddata([xposfgl3, yposfgl3], fluxthrs, [gdat.xposheal])

    cntsthrs = fluxthrs * gdat.expo

    jypos = np.where(abs(yposexpo) < 10.)[0]
    jxpos = np.where(abs(xposexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel(gdat.lablxpostotl)
    axis.set_ylabel(gdat.lablypostotl)

    imag = plt.imshow(fluxthrs[np.amin(jypos):np.amax(jypos)+1, np.amin(jlghprofi):np.amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    figr.savefig(gdat.pathplotcnfg + 'thrs.%s' % gdat.typefileplot)
    plt.close(figr)
    

def plot_init(gdat):
    
    print('Making initial plots...')

    gmod = gdat.fitt

    # make initial plots
    if gdat.boolmakeplot:
        
        if gmod.numbpopl > 0:
            for l in gmod.indxpopl:
                if (gmod.typeelemspateval[l] == 'locl' and gmod.maxmpara.numbelem[l] > 0) and gdat.numbpixl > 1:
                    plot_indxprox(gdat)
        
    for i in gdat.indxener:
        for m in gdat.indxdqlt:
            if gdat.typedata == 'simu' and gmod.boollens:
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
            path = gdat.pathinit + 'expototlmean.%s' % gdat.typefileplot
            tdpy.plot_gene(path, gdat.bctrpara.ener, gdat.expototlmean, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, \
                                                                                            lablydat=gdat.lablexpototl, limtydat=gdat.limtexpo)
        
        for m in gdat.indxdqlt:
            for i in gdat.indxener:
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.hist(gdat.expo[i, :, m], gdat.blimpara.expo)
                axis.set_xlabel(gdat.labltotlpara.expo)
                axis.set_ylabel(gdat.labltotlpara.numbpixl)
                axis.set_xscale('log')
                axis.set_yscale('log')
                plt.tight_layout()
                name = 'histexpo'
                if gdat.numbener > 1:
                    name += 'en%02d' % i
                if gdat.numbdqlt > 1:
                    name += 'evt%d' % m
                path = gdat.pathinit + name + '.%s' % gdat.typefileplot
                figr.savefig(path)
                plt.close(figr)
            
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxdqlt:
                    figr, axis, path = init_figr(gdat, None, 'post', 'expo', '', '', i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.expo, None, None, 'expo', i, m)
                    make_cbar(gdat, axis, imag, i)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
                

def plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                        strgvarb='defl', nameparagenrelem='', indxdefl=None, indxpoplplot=-1, multfact=1., indxenerplot=None, indxdqltplot=None):

    if indxdefl is not None:
        strgvarb += 'sing'
    strgvarb = strgvarb + nameparagenrelem
    
    defl = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    defl *= multfact
   
    if indxenerplot is not None:
        defl = defl[indxenerplot, :, indxdqltplot, ...]

    if indxdefl is not None:
        defl = defl[..., indxdefl]
        strgvarb += '%04d' % indxdefl
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))

    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgvarb, strgstat, strgmodl, indxenerplot, indxdqltplot, indxpoplplot)
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    draw_frambndr(gdat, axis)
  
    deflxpos = defl[:, :, 0]
    deflypos = defl[:, :, 1]
    fact = 4
    ptch = axis.quiver(gdat.anglfact * gdat.xposgridcart[::fact, ::fact], gdat.anglfact * gdat.yposgridcart[::fact, ::fact], \
                       gdat.anglfact * deflxpos[::fact, ::fact], gdat.anglfact * deflypos[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis)
    plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    plt.subplots_adjust()
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgvarb, indxenerplot=None, indxdqltplot=-1, strgcbar=None, \
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
    
    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxdqltplot, indxpoplplot)
   
    maps = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    imag = retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot, indxdqltplot, booltdim=booltdim)
    
    make_cbar(gdat, axis, imag, strgvarb)
    
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    if gdat.boolplotelem:
        supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot)

    plt.tight_layout()
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def defn_tdim(gdat):
    
    # two dimensional
    gdat.fitt.listelemmrkr = ['+', '_', '3']
    if gdat.typedata == 'simu':
        gdat.true.listmrkrhits = ['x', '|', '4']
        gdat.true.listmrkrmiss = ['s', 'o', 'p']
        gdat.true.listlablmiss = ['s', 'o', 'p']
        
    # number of grids
    gdat.numbgrid = 1
    gdat.indxgrid = np.arange(gdat.numbgrid)

    if gdat.typepixl == 'heal' and gdat.boolforccart:
        raise Exception('Cartesian forcing can only used with cart typepixl')


def sample( \
         
         # dictionary defining the model to sample from
         #dictmodl=None, \
    
         #dicttrue=None, \

         # miscelleneaous
         ## type of PDF to sample from
         strgpdfn='post', \

         # diagnostics
         ## Boolean to turn on diagnostic mode
         booldiag=True, \
         ## squeeze exposure to check the low sample limit         
         boolsqzeexpo=False, \
         ### explode exposure to check the large sample limit         
         boolexplexpo=False, \
         ## squeeze proposal scale to check the acceptance ratio
         boolsqzeprop=False, \
         ## explode proposal scale to check the acceptance ratio
         boolexplprop=False, \
         ## factor by which to thin down the data
         factthin=None, \
    
         ## file type of the plot
         typefileplot='png', \
         
         # name of the element feature used for selection
         listnamefeatsele=None, \
         
         # number of processors
         numbproc=None, \
         
         refrlabltotl=None, \
         refrlablpopl=None, \
         fittlablpopl=None, \
    
         # sampling
         ## Boolean flag to make burn-in tempered
         boolburntmpr=False, \
         ## number of sweeps
         numbswep=100000, \
         ## number of samples
         numbsamp=None, \
         ## number of initial sweeps to be burned
         numbburn=None, \
            
         # number of samples for Bootstrap
         numbsampboot=None, \

         # output
         ## Boolean flag to make condensed catalog of elements
         boolcondcatl=True, \

         # plotting
         numbswepplot=None, \
         
         # random state
         ## seed for numpy random number generator
         typeseed=0, \
         ## Boolean flag to re-seed each chain separately
         boolseedchan=True, \
         ## optional deterministic seed for sampling element parameters
         typeseedelem=None, \
         
         # plotting
         # Boolean flag to force the posterior towards the reference
         boolrefeforc=False, \
         # index of the reference catalog towards which the posterior will be forced
         indxrefrforc=None, \

         # model
         ## number of spatial dimensions
         numbspatdims=2, \
         
         # initialization
         ## initialization type
         inittype=None, \
        
         # prior
         priotype='logt', \
         
         # saving and recovering the state of the sampler
         # Boolean flag to save the state of the MCMC
         boolsavestat=False, \
         namesavestat=None, \
         # recover the state from a previous run
         namerecostat=None, \
         forcsavestat=False, \
         namestattrue=None, \
        
         # user defined function to plot
         plot_func=None, \

         # proposals
         ## Boolean flag to turn on proposals on element parameters
         boolpropcomp=True, \
         boolpropcova=True, \
         propwithsing=True, \
         # type of covariance estimation
         typeopti='none', \
         
         # user-defined likelihood function
         retr_llik=None, \


         # arguments common among sample() and wrappers
         
         # modes of operation
         ## perform an additional run sampling from the prior
         checprio=False, \

         # name of the configuration
         strgcnfg=None, \
        
         # elements
         ## reference element catalog
         ## Boolean flag to associate elements to those in reference element catalogs
         boolasscrefr=None, \
         
         # data
         ## type of data
         ### 'simu': simulated data
         ### 'inpt': input data
         ### 'real': real data retrieved from databases
         typedata=None, \
         
         # user interaction
         ## type of verbosity
         typeverb=1, \

         ## base path where visuals and data will be written
         pathbase=None, \
        
         ## transdimensional proposal probabilities
         probtran=None, \
         probspmr=None, \
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
         # standard deviation of the Gaussian from which the angular splitting will be drawn for splits and merges
         radispmr=None, \
         
         # factor by which to multiply the prior for each additional degree of freedom in the model
         factpriodoff=None, \

        ):

    # preliminary setup
    # construct the global object 
    gdat = tdpy.gdatstrt()
    for attr, valu in locals().items():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)
    
    # check inputs
    if gdat.numbburn > gdat.numbswep:
        raise Exception('Bad number of burn-in sweeps.')
    if gdat.factthin > gdat.numbswep - gdat.numbburn or gdat.factthin < 1:
        raise Exception('Bad thinning factor.')
    
    setup_pcat(gdat)

    gdat.fitt = tdpy.gdatstrt()
    
    print('gdat.liststrgfeatpara')
    print(gdat.liststrgfeatpara)
    for strgfeatpara in gdat.liststrgfeatpara:
        setattr(gdat.fitt, strgfeatpara + 'para', tdpy.gdatstrt())
        setattr(gdat, strgfeatpara + 'para', tdpy.gdatstrt())
    
    gdat.strgswep = '%d' % (gdat.numbswep)
    
    if gdat.typeverb > 0:
        if gdat.boolburntmpr:
            print('Warning: Tempered burn-in.')
    
    ## time stamp
    gdat.strgtimestmp = tdpy.retr_strgtimestmp()
    
    gdat.strgvers = 'v0.3'
    if gdat.typeverb > 0:
        print('PCAT %s started at %s.' % (gdat.strgvers, gdat.strgtimestmp))
        print('Configuration %s' % gdat.strgcnfg)
    
    # string describing the number of sweeps
    gdat.strgnumbswep = '%d' % gdat.numbswep
    
    # output paths
    gdat.strgcnfg = '%d' % gdat.strgnumbswep
    gdat.pathoutpcnfg = retr_pathoutpcnfg(gdat.pathbase, gdat.strgcnfg)

    # physical constants
    gdat.prsccmtr = 3.086e18
    gdat.ergsgevv = 624.151
        
    gdat.listnamepdir = ['forw', 'reve']
    gdat.listlablpdir = ['f', 'r']
        
    # start the timer
    gdat.timerealtotl = time.time()
    gdat.timeproctotl = time.clock()
   
    # list of parameter types
    ## 'genr': generative parameters
    ## 'deri': derived parameters
    gdat.liststrgtypepara = ['genr', 'deri']
    
    # list of parameter groups
    ## 'base': base parameters
    ## 'elem': element parameters
    gdat.liststrggroppara = ['base', 'elem', 'full']
    
    booltemp = chec_statfile(gdat.pathbase, gdat.strgcnfg, 'gdatmodi')
    if booltemp:
        print('gdatmodi already exists. Skipping...')
    else:
    
        # create output folder for the run
        os.system('mkdir -p %s' % gdat.pathoutpcnfg)

        # write the list of arguments to file
        fram = inspect.currentframe()
        listargs, temp, temp, listargsvals = inspect.getargvalues(fram)
        fileargs = open(gdat.pathoutpcnfg + 'cmndargs.txt', 'w')
        fileargs.write('PCAT call arguments\n')
        for args in listargs:
            fileargs.write('%s = %s\n' % (args, listargsvals[args]))
        fileargs.close()
        
        # write the list of arguments to file
        fileargs = open(gdat.pathoutpcnfg + 'args.txt', 'w')
        fileargs.write('PCAT call arguments\n')
        for args in listargs:
            fileargs.write('%20s %s\n' % (args, listargsvals[args]))
        fileargs.close()
        
        gdat.refr = tdpy.gdatstrt()
        
        for strgstat in ['this', 'next']:
            setattr(gdat, strgstat, tdpy.gdatstrt())
        
        ## number of processes
        gdat.strgproc = os.uname()[1]
        if gdat.numbproc is None:
            gdat.numbproc = 1
    
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
        
        if gdat.boolsqzeprop:
            gdat.stdp[:]= 1e-100
        
        if gdat.boolexplprop:
            gdat.stdp[:] = 1e100

        if (gdat.boolsqzeprop or gdat.boolexplprop) and gdat.typeopti == 'hess':
            raise Exception('')

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

        if gdat.factpriodoff is None:
            gdat.factpriodoff = 1.
        
        print('gdat.fitt.lablrootpara')
        print(gdat.fitt.lablrootpara)
        setp_varb(gdat, 'lliktotl', labl=['$L$', ''], strgmodl='fitt')
        
        if gdat.strgpdfn == 'prio':
            gdat.lablsampdist = 'Prior'
        if gdat.strgpdfn == 'post':
            gdat.lablsampdist = 'Posterior'

        if gdat.numbswepplot is None:
            gdat.numbswepplot = 50000
   
        gdat.numbplotfram = gdat.numbswep / gdat.numbswepplot

        if gdat.checprio:
            gdat.liststrgpdfn = ['prio', 'post']
        else:
            gdat.liststrgpdfn = ['post']

        # initialization type
        if gdat.inittype is None:
            gdat.inittype = 'rand'

        # process index
        gdat.indxproc = np.arange(gdat.numbproc)

        setp_indxswepsave(gdat)
        
        if gdat.typeopti == 'hess':
            pathopti = gdat.pathoutpcnfg + 'opti.h5'
            if os.path.exists(pathopti):
                thisfile = h5py.File(pathopti, 'r')
                if thisfile['stdp'][()].size == gdat.stdp.size:
                    print('Recovering the proposal scale from the previous run...')
                    gdat.stdp = thisfile['stdp'][()]
                thisfile.close()

        
        

        # turn off relevant proposal types
        gdat.numbprop = 5
        gdat.indxprop = np.arange(gdat.numbprop)
        
        gdat.numbstdp = gmod.numbparagenrbase - gmod.numbpopl
        cntr = 0
        for l in gmod.indxpopl:
            for nameparagenrelem in gmod.namepara.genrelem[l]:
                #setattr(gdat.fitt.indxpara.genrelemkind, nameparagenrelem + 'pop%d' % l, gdat.numbstdp)
                cntr += 1
        gdat.numbstdp += cntr
        
        gdat.lablstdp = np.copy(np.array(gmod.labltotlpara.genrbase[gmod.numbpopl:]))
        gdat.namestdp = np.copy(np.array(gmod.namepara.genrbase[gmod.numbpopl:]))
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
        gdat.fitt.this.indxparagenrfullelem = retr_indxparagenrelemfull(gdat, indxelemfull, 'fitt')
        
        gdat.indxstdppara = np.zeros(gmod.numbparagenr, dtype=int) - 1
        cntr = 0
        gdat.indxstdppara[gmod.numbpopl:gmod.numbparagenrbase] = gmod.indxpara.genrbase[gmod.numbpopl:] - gmod.numbpopl
        if gmod.numbpopl > 0:
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
    

        gmod = gdat.fitt
        # proposal scale
        if gmod.boollens or gdat.typedata == 'simu':
            
            gdat.stdp = 1e-4 + np.zeros(gdat.numbstdp)
            
            if gmod.typemodltran == 'pois' and gmod.numbpopl > 0:
                if gmod.maxmpara.numbelem[0] > 0:
                    gdat.stdp[gdat.indxstdpmeanelempop0] = 1e-1
            
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.sigcen00evt0]] = 3e-2
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.bacpback0000en00]] = 1e-3
            gdat.stdp[gdat.indxstdppara[gmod.indxpara.bacpback0000en00]] = 1e-1
            
            if gmod.boollens:
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.xpossour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.ypossour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.fluxsour]] = 1e-2
                if gdat.numbener > 1:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.sindsour]] = 1e-3
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.sizesour]] = 1e-1
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.ellpsour]] = 1e-1
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.anglsour]] = 1e-1
            if gmod.typeemishost != 'none':
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.xposhostisf0]] = 3e-4
                gdat.stdp[gdat.indxstdppara[gmod.indxpara.yposhostisf0]] = 3e-4
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
                
                if gmod.typemodltran == 'pois' and gmod.numbpopl > 0:
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
                
                if gmod.numbpopl > 0:
                    gdat.stdp[gdat.indxstdppop0flux] = 8e-2

                    if gmod.spectype[0] == 'colr':
                        gdat.stdp[gdat.indxstdppop0sindcolr0001] = 8e-2
                        gdat.stdp[gdat.indxstdppop0sindcolr0002] = 2e-1
            
            if gdat.typeexpr == 'chan':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
                
                if gmod.typemodltran == 'pois' and gmod.numbpopl > 0:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.meanelem]] = 2e-1
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.sloppriofluxpop0]] = 2e-1
                
                if gmod.numbpopl > 0 and gdat.boolbindspat:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.psfp]] = 4e-1
                
                if gdat.indxenerincl.size == 5 and (gdat.indxenerincl == np.arange(5)).all():
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en01')]] = 3e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en02')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en03')]] = 2e-2
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en04')]] = 1e-2
                elif gdat.indxenerincl.size == 2 and (gdat.indxenerincl == np.array([2])).all():
                    gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
                
                if gmod.numbpopl > 0:
                    if gdat.boolbindspat:
                        gdat.stdp[gdat.fitt.indxpara.genrelemkind.xpospop0] = 2e-2
                        gdat.stdp[gdat.fitt.indxpara.genrelemkind.ypospop0] = 2e-2
                        if gdat.numbener > 1:
                            gdat.stdp[gdat.indxstdppop0sind] = 2e-1
                    gdat.stdp[gdat.indxstdppop0flux] = 2e-1
            
            if gdat.typeexpr == 'gmix':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
            
                if gmod.typemodltran == 'pois' and gmod.numbpopl > 0:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.meanelem]] = 2e-1
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.slopprionobjpop0]] = 3e-1
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.gwdtsloppop0]] = 3e-1

                if gmod.typeevalpsfn != 'none' and gdat.boolmodipsfn:
                    gdat.stdp[gdat.indxstdppara[gmod.indxpara.psfp]] = 4e-1
                
                gdat.stdp[gdat.indxstdppara[getattr(gmod.indxpara, 'bacpback0000en00')]] = 2e-2
            
                if gmod.numbpopl > 0:
                    gdat.stdp[gdat.fitt.indxpara.genrelemkind.xpospop0] = 4e-2
                    gdat.stdp[gdat.fitt.indxpara.genrelemkind.ypospop0] = 4e-2
                    gdat.stdp[gdat.fitt.indxpara.genrelemkind.nobjpop0] = 3e-1
                    gdat.stdp[gdat.indxstdppop0gwdt] = 5e-1

            if gdat.typeexpr == 'fire':
                gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
        
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








        ## initialization
        gdat.fitt.this = tdpy.gdatstrt()
        gdat.fitt.this.indxpara = tdpy.gdatstrt()
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
        #            for strgfeat in gmod.namepara.glob:
        #            #for strgfeat in gmod.namepara.genrelem[l]:
        #                if strgfeat[:-4] == 'etag':
        #                    continue
        #                #setp_varb(gdat, strgfeat)
        #                #if strgfeat in gmod.namepara.elem:
        #                #    setp_varb(gdat, strgfeat + 'prio')
        
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
        if gdat.boolmakeplot and gdat.boolmakeplotinit:
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
            if gdat.typepixl == 'heal' and (gdat.xposcntr != 0. or gdat.yposcntr != 0.):
                for q in gdat.indxrefr:
                    for l in gmod.indxpopl:
                        rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.xposcntr), rad2deg(gdat.yposcntr), 0.], deg=True, eulertype='ZYX')
                        gdat.refr.dictelem[q]['ypos'][0, :], gdat.refrxpos[0, :] = rttr(pi / 2. - gdat.refrypos[0, :], gdat.refrxpos[0, :])
                        gdat.refr.dictelem[q]['ypos'][0, :] = pi / 2. - gdat.refrypos[0, :]

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
        #    gdat.refrfluxbrgt, gdat.refrfluxbrgtassc = retr_fluxbrgt(gdat, gdat.refrxpos, gdat.refrypos, gdat.refrflux[0, :])
        
        print('gdat.liketype')
        print(gdat.liketype)

        print('Data settings')
        print('gdat.numbener')
        print(gdat.numbener)
        print('gdat.numbdqlt')
        print(gdat.numbdqlt)

        print('Model settings')
        print('gdat.fitt.numbpopl')
        print(gdat.fitt.numbpopl)
        print('gdat.fitt.numbparagenrbase')
        print(gdat.fitt.numbparagenrbase)
        
        for strgmodl in gdat.liststrgmodl:
            for l in gmod.indxpopl:
                for strgfeat, strgpdfn in zip(gmod.namepara.genrelemmodu[l], gmod.liststrgpdfnmodu[l]):
                    if strgpdfn == 'tmpl':
                        if gdat.xposprio is None or gdat.yposprio is None:
                            gdat.xposprio = np.concatenate((gmod.xpos))
                            gdat.yposprio = np.concatenate((gmod.ypos))
                        gdat.numbspatprio = gdat.xposprio.size
        
                        # spatial template for the catalog prior
                        # temp -- this should move outside the if
                        gdat.pdfnspatpriotemp = np.zeros((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                        for k in range(gdat.numbspatprio):
                            gdat.pdfnspatpriotemp[:] += 1. / np.sqrt(2. * np.pi) / gdat.stdvspatprio * \
                                                                exp(-0.5 * (gdat.blimpara.xposcartmesh - gdat.xposprio[k])**2 / gdat.stdvspatprio**2) * \
                                                                exp(-0.5 * (gdat.blimpara.yposcartmesh - gdat.yposprio[k])**2 / gdat.stdvspatprio**2)
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
                    skycobjt = ap.coordinates.SkyCoord("galactic", l=gdat.refr.dictelem[q]['xpos'][0, :] * 180. / pi, \
                                                                    b=gdat.refr.dictelem[q]['ypos'][0, :] * 180. / pi, unit='deg')
                    rasc = skycobjt.fk5.ra.degree
                    decl = skycobjt.fk5.dec.degree
                    xpos, ypos = wcso.wcs_world2pix(rasc, decl, 0)
                    xpos -= gdat.numbpixlxposshft + gdat.numbsidecarthalf
                    ypos -= gdat.numbpixlyposshft + gdat.numbsidecarthalf
                    xpos *= gdat.sizepixl
                    ypos *= gdat.sizepixl
                    gdat.refr.dictelem[q]['xpos'][0, :] = ypos
                    gdat.refr.dictelem[q]['ypos'][0, :] = xpos

            ## preprocess reference element features
            for q in gdat.indxrefr:
                # temp -- this should depend on q
                # temp -- this does not properly calculate uncertainties
                gdat.refrgang[q] = np.zeros((3, gdat.refr.dictelem[q]['xpos'].shape[1]))
                gdat.refraang[q] = np.zeros((3, gdat.refr.dictelem[q]['xpos'].shape[1]))
                gdat.refrgang[q][:, :] = retr_gang(gdat.refr.dictelem[q]['xpos'][0, :], gdat.refr.dictelem[q]['ypos'][0, :])[None, :]
                gdat.refraang[q][:, :] = retr_aang(gdat.refr.dictelem[q]['xpos'][0, :], gdat.refr.dictelem[q]['ypos'][0, :])[None, :]

            # save all reference element features
            for strgfeat in gdat.refr.namepara.elem.totl:
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
                gdat.indxrefrpntsrofi[q] = np.where((np.fabs(gdat.refr.dictelem[q]['xpos'][0, :]) < gdat.maxmgangdata) & \
                                                                        (np.fabs(gdat.refr.dictelem[q]['ypos'][0, :]) < gdat.maxmgangdata))[0]
            for strgfeat in gdat.refr.namepara.elem.totl:
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
                gdat.refr.numbelem[q] = gdat.refr.dictelem[q]['xpos'].shape[1]
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
            if not (gdat.typedata == 'simu' and gmod.numbparaelem == 0):
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
                        blimfeatfrst = getattr(gdat.blimpara, strgfeatfrst)
                        hist = np.histogram(refrfeatfrst[q][0, :], blimfeatfrst)[0]
                        setattr(gdat.refr, 'hist' + strgfeatfrst + 'pop%d' % q, hist)
                        for strgfeatseco in gdat.refr.namepara.elem[q]:
                            if strgfeatseco.startswith('etag'):
                                continue
                            refrfeatseco = getattr(gdat.refr, strgfeatseco)
                            
                            strgfeattdim = strgfeatfrst + strgfeatseco + 'pop%d' % q
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                            
                            if len(refrfeatseco[q]) > 0:
                                blimfeatseco = getattr(gdat.blimpara, strgfeatseco)
                                hist = np.histogram2d(refrfeatfrst[q][0, :], refrfeatseco[q][0, :], bins=(blimfeatfrst, blimfeatseco))[0]
                                setattr(gdat.refr, 'hist' + strgfeattdim, hist)
            
        if gmod.numbpopl > 0:
            # plot settings
            ## upper limit of histograms
            if gdat.limtydathistfeat is None:
                gdat.limtydathistfeat = [0.5, max(100., 10**np.ceil(np.log10(gdat.refr.numbelemtotl)))]
                #gdat.limtydathistfeat = [0.5, max(100., 10**np.ceil(np.log10(gmod.maxmpara.numbelemtotl)))]

        # initial plots
        if gdat.boolmakeplot and gdat.boolmakeplotinit:
            # problem-specific plots
            if gdat.boolmakeplotintr:
                plot_intr(gdat)
                #plot_pert()
                #plot_king(gdat)
                plot_func(gdat)
                #plot_3fgl_thrs(gdat)
                #if gdat.typeexpr == 'ferm':
                #    plot_fgl3(gdat)
        
        # find the pixels at which data count maps have local maxima
        if gdat.typepixl == 'cart':
            for i in gdat.indxener:
                for m in gdat.indxdqlt:
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
        with open(gdat.pathoutpcnfg + 'stat.p', 'wb') as thisfile:
        	pickle.dump(np.random.get_state(), thisfile)
        
        # process lock for simultaneous plotting
        lock = mp.Manager().Lock()
        
        if gdat.typeverb > 0:
            print('Writing the global state to the disc before spawning workers...')
        path = gdat.pathoutpcnfg + 'gdatinit'
        writfile(gdat, path) 
        gdat.filestat = open(gdat.pathoutpcnfg + 'stat.txt', 'w')
        gdat.filestat.write('gdatinit written.\n')
        gdat.filestat.close()
        
        # exit before running the sampler
        if gdat.boolsimuonly:
            if gdat.typeverb > 0:
                print('Mock dataset is generated. Quitting...')
            return gdat.strgcnfg
        
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
        proc_anim(gdat.strgcnfg)

    if gdat.typeverb > 0:
        print('The output is at ' + gdat.pathoutpcnfg)
        if gdat.boolmakeplot:
            print('The plots are at ' + gdat.pathplotcnfg)
        print('PCAT has run successfully. Returning to the OS...')

    return gdat.strgcnfg


def retr_dictpcatinpt():
    
    dictpcatinpt = dict()
    dictpcatinpt['dicttrue'] = dict()
    dictpcatinpt['dictfitt'] = dict()
    dictpcatinpt['dictfitt']['maxmpara'] = tdpy.gdatstrt()
    dictpcatinpt['dicttrue']['maxmpara'] = tdpy.gdatstrt()
    
    return dictpcatinpt


def sample_parallel( \
             dictpcatinptvari, \
             listlablcnfg, \
             forcneww=False, \
             forcprev=False, \
             strgpara=False, \
             
             # input dictionary to PCAT
             dictpcatinpt=None, \
            
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
    
    numbiter = len(dictpcatinptvari)
    indxiter = np.arange(numbiter) 
    
    pathbase = os.environ["PCAT_DATA_PATH"] + '/'
    
    cntrcomp = 0
    
    if boolexecpara:
        cntrproc = 0

    liststrgcnfg = []
    listpridchld = []
    listnamecnfgextn = list(dictpcatinptvari.keys())
    print('listnamecnfgextn')
    print(listnamecnfgextn)
    for k, strgcnfgextn in enumerate(listnamecnfgextn):
        
        if strgcnfgextnexec is not None:
            if strgcnfgextn != strgcnfgextnexec:
                continue
        
        strgcnfg = inspect.stack()[1][3] + '_' + strgcnfgextn
        
        if dictpcatinpt is None:
            dictvarbtemp = dict()
        else:
            dictvarbtemp = deepcopy(dictpcatinpt)
        
        for strgvarb, valu in dictpcatinptvari[strgcnfgextn].items():
            dictvarbtemp[strgvarb] = valu
        dictvarbtemp['strgcnfg'] = strgcnfg
        
        liststrgcnfgprev = retr_liststrgcnfgprev(strgcnfg, pathbase)
        cntrcomp += 1

        if (not forcneww and strgcnfgextnexec is None or forcprev and strgcnfgextnexec is not None) and len(liststrgcnfgprev) > 0:
            print('Found at least one previous run with the configuration %s' % strgcnfg)
            print('Skipping...')
            liststrgcnfg.append(liststrgcnfgprev[-1])
        else:
            if len(liststrgcnfgprev) > 0:
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
                    strgcnfg = sample(**dictvarbtemp)
                    os._exit(0)
            else:
                print('Calling the main PCAT function without forking a child...')
                liststrgcnfg.append(sample(**dictvarbtemp))
    
    if boolexecpara and strgcnfgextnexec is None:
        for prid in listpridchld:
            os.waitpid(prid, 0)
        if cntrproc > 0:
            print('Exiting before comparion plots because of parallel execution...')
            return
    
    if cntrcomp == 0:
        print('Found no runs...')

    print('Final-processing run outputs...')
    for strgcnfg in liststrgcnfg:
        print('strgcnfg')
        print(strgcnfg)
        proc_finl(strgcnfg=strgcnfg, strgpdfn='post')
        proc_anim(strgcnfg)
    
    strgtimestmp = tdpy.retr_strgtimestmp()
    
    if strgcnfgextnexec is not None or namexaxivari is None: 
        return
    
    print('Making plots to compare the output of different PCAT runs...')
     
    if 'boolsimuonly' in dictvarb and dictvarb['boolsimuonly']:
        listgdat = retr_listgdat(liststrgcnfg, typegdat='init')
    else:
        listgdat = retr_listgdat(liststrgcnfg)
    
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
        listgdatprio = retr_listgdat(liststrgcnfg, typegdat='finlprio')
        
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

    pathvisuinsp = '%s/visu/%s_%s/' % (gdat.pathbase, strgtimestmp, inspect.stack()[1][3])
    cmnd = 'mkdir -p %s' % pathvisuinsp 
    os.system(cmnd)
    cmnd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%smrgd.%s' % (pathvisuinsp, gdat.typefileplot)
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
        
        if listtypevarbcomp[indxlist] == 'pctl':
            trueyaxi = getattr(listgdat[0], 'true' + listnamevarbcomp[indxlist])
        else:
            trueyaxi = getattr(listgdat[0], 'true' + listtypevarbcomp[indxlist] + listnamevarbcomp[indxlist])
        
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
        indxstrgcnfgyerr = np.where((yerr[0, :] > 0.) | (yerr[1, :] > 0.))[0]
        if indxstrgcnfgyerr.size > 0:
            temp, listcaps, temp = axis.errorbar(indxiter[indxstrgcnfgyerr]+1., ydat[indxstrgcnfgyerr], yerr=yerr[:, indxstrgcnfgyerr], \
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
        
        pathfull = '%s%s_%s_%s.%s' % (pathvisuinsp, strgtimestmp, inspect.stack()[1][3], liststrgvarbtotl[indxlist], gdat.typefileplot)
        print('Writing to %s...' % pathfull)
        plt.savefig(pathfull)
        plt.close(figr)
    
        cmnd += ' %s' % pathfull

    print(cmnd)
    os.system(cmnd)

    print('Making animations...')
    for strgcnfg in liststrgcnfg:
        print('Working on %s...' % strgcnfg)
        proc_anim(strgcnfg=strgcnfg)
    
    print('Compiling run plots...')
    cmnd = 'python comp_strgcnfg.py'
    for strgcnfg in liststrgcnfg: 
        cmnd += ' %s' % strgcnfg
    os.system(cmnd)

    return liststrgcnfg


class logg(object):
    
    def __init__(self, gdat):
        self.terminal = sys.stdout
        gdat.pathstdo = gdat.pathoutpcnfg + 'stdo.txt'
        self.log = open(gdat.pathstdo, 'a')
        pathlink = gdat.pathplotcnfg + 'stdo.txt'
        os.system('ln -s %s %s' % (gdat.pathstdo, pathlink))
    
    def write(self, strg):
        self.terminal.write(strg)
        self.log.write(strg)  

    def flush(self):
        pass


def worktrac(pathoutpcnfg, lock, strgpdfn, indxprocwork):
	
    try:
        return work(pathoutpcnfg, lock, strgpdfn, indxprocwork)
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def opti_hess(gdat, gdatmodi):
    
    gmod = gdat.fitt

    if gmod.numbpopl > 0:
        cntr = 0
        for l in gmod.indxpopl:
            for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                gdatmodi.indxparastdp[gmod.numbparagenrbase-gmod.numbpopl+cntr] = np.concatenate(gdatmodi.this.indxparagenrelemfull[nameparagenrelem])
                cntr += 1
    
    if gmod.numbpopl > 0:
        gdatmodi.next.indxelemfull = gdatmodi.this.indxelemfull
        gdatmodi.next.indxparagenrelemfull = gdatmodi.this.indxparagenrelemfull
    else:
        gdatmodi.next.indxparagenrelemfull = None

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
    
    #if gmod.numbpopl > 0:
    #    gdatmodi.dictmodi = [[] for l in gmod.indxpopl]
    #    for l in gmod.indxpopl:
    #        gdatmodi.dictmodi[l] = dict()
    #        gdatmodi.dictmodi[l][gmod.nameparagenrelemampl[l] + 'indv'] = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[gmod.nameparagenrelemampl[l]][l]]
    #        for nameparagenrelem in gmod.namepara.genrelem[l]:
    #            gdatmodi.dictmodi[l]['stdv' + nameparagenrelem + 'indv'] = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l][nameparagenrelem]]
    #if gmod.numbpopl > 0:
    #    gdatmodi.this.indxparagenrelemfullconc = np.concatenate([gdatmodi.this.indxparagenrelemfull[l]['full'] for l in gmod.indxpopl])
    #if gdat.boolpropcomp:
    #    indxsamptranprop = gdatmodi.this.indxparagenrelemfullconc
    #else:
    #    indxsamptranprop = []
    
    deltlpos[1, 1] = gdatmodi.this.lliktotl
    for indxstdpfrst in gdat.indxstdpprop:
        for indxstdpseco in gdat.indxstdpprop:
            
            if indxstdpfrst > indxstdpseco:
                continue
            
            if indxstdpfrst == indxstdpseco:
                
                #if gmod.numbpopl > 0:
                #    if k in gdatmodi.this.indxparagenrelemfullconc:
                #        indxtrapmoditemp = k - gmod.numbpara.genr.base
                #        indxpoplmoditemp = np.array([np.amin(np.where(indxtrapmoditemp // gmod.numbparagenrelemcumr == 0))])
                #        numbparapoplinittemp = indxtrapmoditemp - gmod.numbparagenrelemcuml[indxpoplmoditemp[0]]
                #        indxelemmoditemp = [numbparapoplinittemp // gmod.numbparagenrelemsing[indxpoplmoditemp[0]]]
                #        gmod.indxpara.genrelemmoditemp = numbparapoplinittemp % gmod.numbparagenrelemsing[indxpoplmoditemp[0]]
                #        nameparagenrelem = gmod.namepara.genrelem[indxpoplmoditemp[0]][gmod.indxpara.genrelemmoditemp] 
                #        indxsampampltemp = k - gmod.indxpara.genrelemmoditemp + gmod.indxpara.genrelemampl[indxpoplmoditemp[0]]
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
                    
                    gdatmodi.next.paragenrscalfull = icdf_paragenrscalfull(gdat, 'fitt', gdatmodi.next.paragenrunitfull, gdatmodi.next.indxparagenrelemfull)
                    
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
            if gdat.booldiag and not np.isfinite(gdatmodi.next.paragenrscalfull).all():
                raise Exception('')
            if gdatmodi.hess[indxstdpfrst, indxstdpseco] == 0.:
                raise Exception('')

    gdatmodi.hess[np.where(gdatmodi.hess == 0)] = 10.

    # temp
    #gdatmodi.stdpmatr = np.sqrt(linalg.inv(gdatmodi.hess))
    numbdoffefff = gmod.numbparagenrbase
    if gmod.numbpopl > 0:
        numbdoffefff += gmod.numbparagenrelempopl * 10
    gdatmodi.stdpmatr = np.sqrt(1. / gdatmodi.hess) / np.sqrt(numbdoffefff)
    
    if (gdatmodi.stdpmatr == 0).any():
        raise Exception('')
    
    gdatmodi.stdp = gdatmodi.stdpmatr[gdat.indxstdp, gdat.indxstdp]
    

def worksamp(gdat, lock, strgpdfn='post'): 
    
    pathorig = gdat.pathoutpcnfg + 'stat.txt'
    pathlink = gdat.pathplotcnfg + 'stat.txt'
    os.system('ln -s %s %s' % (pathorig, pathlink))
    
    if gdat.numbproc == 1:
        worktrac(gdat.pathoutpcnfg, lock, strgpdfn, 0)
    else:
        if gdat.typeverb > 0:
            print('Forking the sampler...')

        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(worktrac, gdat.pathoutpcnfg, lock, strgpdfn)
        pool.map(workpart, gdat.indxproc)

        pool.close()
        pool.join()
    
    gdat.filestat = open(gdat.pathoutpcnfg + 'stat.txt', 'a')
    gdat.filestat.write('gdatmodi%s written.\n' % strgpdfn)
    gdat.filestat.close()


def work(pathoutpcnfg, lock, strgpdfn, indxprocwork):
    
    print('Worker #%d' % indxprocwork)
    
    # read the initial global object, gdatinit
    path = pathoutpcnfg + 'gdatinit'
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
    gdatmodi.this.indxpara = tdpy.gdatstrt()
    gdatmodi.next.indxpara = tdpy.gdatstrt()
    gdatmodi.indxprocwork = indxprocwork
    
    gdatmodi.this = gdat.fitt.this

    # path of gdatmodi
    gdatmodi.pathgdatmodi = gdat.pathoutpcnfg + 'gdatmodi%04d' % gdatmodi.indxprocwork + gdat.strgpdfn
    
    print('Determining the parameter indices of the fitting model with only the floating parameters...')

    gdatmodi.booldone = False
    gdatmodi.lock = lock
    gdatmodi.indxprocwork = indxprocwork
    
    # find the list of variables for which the posterior will be calculated
    if not gdat.boolsimuonly:
        
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
   
    if gdat.booldiag:
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
    for k in gmod.indxpara.genrbase:
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
    
    # 'temp' something is wrong with this 
    #gdatmodi.sampmaxmllik = np.copy(gdatmodi.this.paragenrscalfull)
    
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
            path = gdat.pathoutpcnfg + 'opti.h5'
            if gdat.typeverb > 0:
                print('Writing the estimated covariance matrix to %s...' % path)
            thisfile = h5py.File(path, 'w')
            thisfile.create_dataset('stdp', data=gdatmodi.stdp)
            thisfile.close()
            
            if gdat.boolmakeplot:
                
                xdat = gdat.indxstdp
                ydat = gdatmodi.stdp
                
                pathopti = getattr(gdat, 'path' + gdat.strgpdfn + 'opti')
                path = pathopti + 'stdv%d.%s' % (gdatmodi.indxprocwork, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, ydat, scalydat='logt', \
                                lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist', limtydat=[np.amin(ydat) / 2., 2. * np.amax(ydat)])
                
                # plot uncertainties of element parameters as a function of amplitude parameter
                if gmod.numbpopl > 0:
                    for l in gmod.indxpopl:
                        for nameparagenrelem in gmod.namepara.genrelem[l]:
                            path = pathopti + 'stdv' + nameparagenrelem + 'pop%d.%s' % (l, gdat.typefileplot)
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
                            #                                 lablydat=r'$\sigma_{%s}$%s' % (getattr(gmod.lablpara, nameparagenrelem), \
                            #                                 getattr(gmod.lablpara, nameparagenrelem + 'unit')), plottype=['scat', 'lghtline'])
                            
                            tdpy.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                             lablydat=r'$\sigma_{%s}$' % getattr(gmod.lablpara, nameparagenrelem), plottype=['scat', 'lghtline'])


        if gdat.typeverb > 1:
            print('-' * 10)
            print('Sweep %d' % gdatmodi.cntrswep)

        # decide whether to make a frame
        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and \
                                                gdatmodi.indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                and gdat.boolmakeplotfram and gdat.boolmakeplot
        
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
            #gdatmodi.this.tmprlposelem = -1000. * (1. - gdatmodi.this.facttmpr) * np.concatenate(gdatmodi.this.indxparagenrelemfull['full']).size
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

        if gdat.booldiag:
        
            for k in gmod.indxpara.genrbase:
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
        if gdat.booldiag:
            
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
            #        if not indxparagenrfulltemp in gdatmodi.this.indxparagenrelemfull[l]['full']:
            #            indxsampclos.remove(indxparagenrfulltemp)
            #indxsampclos = np.array(indxsampclos)
            #if indxsampclos.size > 0:
            #    print 'Warning! State is too close to 0!'
            #    print gmod.namepara[indxsampclos]

            #indxsampclos = np.where((gdatmodi.this.paragenrscalfull > 0.99) & (gdatmodi.this.paragenrscalfull % 1. != 0.))[0]
            #indxsampclos = list(indxsampclos)
            #for indxparagenrfulltemp in indxsampclos:
            #    for l in gmod.indxpopl:
            #        if not indxparagenrfulltemp in gdatmodi.this.indxparagenrelemfull[l]['full']:
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
            
            if gmod.numbpopl > 0:
                if gmod.boolelemsbrtdfncanyy:
                    thissbrtdfnc = getattr(gdatmodi.this, 'sbrtdfnc')
                    frac = np.amin(thissbrtdfnc) / np.mean(thissbrtdfnc)
                    cntppntschec = retr_cntp(gdat, thissbrtdfnc)
                    if np.amin(cntppntschec) < -0.1 and frac < -1e-3:
                        raise Exception('thissbrtdfnc went negative by %.3g percent.' % (100. * frac))
                    
            # check the population index
            if (gdatmodi.this.cntpmodl <= 0.).any() or not (np.isfinite(gdatmodi.this.cntpmodl)).all():
                raise Exception('Current flux model is not positive')

            if gmod.numbpopl > 0:
                for l in gmod.indxpopl:
                    if gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]] != len(gdatmodi.this.indxelemfull[l]):
                        raise Exception('Number of elements is inconsistent with the element index list.')

                    if gdatmodi.this.paragenrscalfull[gmod.indxpara.numbelem[l]] != len(gdatmodi.this.indxelemfull[l]):
                        raise Exception('Number of elements is inconsistent across data structures.')
                    
                    for k, nameparagenrelem in enumerate(gmod.namepara.genrelem[l]):
                        if gmod.scalpara.genrelem[l][k] == 'gaus' or gmod.scalpara.genrelem[l][k] == 'igam' \
                                                                                            or gmod.scalpara.genrelem[l][k] == 'expo':
                            continue
                        comp = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l][nameparagenrelem]]
                        minm = getattr(gdat.fitt.minmpara, nameparagenrelem)
                        maxm = getattr(gdat.fitt.maxmpara, nameparagenrelem)
                        indxtemp = np.where((comp < minm) | (comp > maxm))[0]
                        if indxtemp.size > 0:
                            raise Exception('A component of an element went outside the prior range.')
        
            stopchro(gdat, gdatmodi, 'diag')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            initchro(gdat, gdatmodi, 'save')
        
            if gdat.boolsavestat:
                
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
                    for gmod.namepara.genr.base in gmod.namepara.genr.base:
                        valu = gdatmodi.this.paragenrscalfull[gmod.indxpara.genrbase]
                        thisfile.create_dataset(gmod.namepara.genr.base, data=valu)
                    if gmod.numbpopl > 0:
                        for l in gmod.indxpopl:
                            for nameparagenrelem in gmod.namepara.genrelem[l]:
                                comp = gdatmodi.this.paragenrscalfull[gdatmodi.this.indxparagenrelemfull[l][nameparagenrelem]]
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
            
            if gdat.booldiag:
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
        
        if gdat.booldiag:
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
            
            numbpara = gmod.numbparagenrbase
            if gmod.numbpopl > 0:
                for l in gmod.indxpopl:
                    numbpara += gdatmodi.this.indxparagenrelemfull[l]['full'].size
            if gmod.numbpopl > 0:
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
            if gmod.numbpopl > 0:
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
            if gmod.numbpopl > 0:
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
                if name == 'llik' and gdat.numbpixl > 1 and gmod.numbpopl > 0 and booltemp:
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

