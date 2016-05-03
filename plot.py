
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

def plot_look(globdata):

    indxpixlproxsize = zeros(globdata.numbpixl)
    figr, axis = plt.subplots(figsize=(10, 6))
    for h in range(numbspecprox):
        for j in ipixl:
            indxpixlproxsize[j] = indxpixlprox[h][j].size
        binspixlsize = logspace(log10(amin(indxpixlproxsize)), log10(amax(indxpixlproxsize)), 100)
        hist, bins = histogram(indxpixlproxsize, binspixlsize)
        mean = sqrt(bins[:-1] * bins[1:])
        axis.loglog(mean, hist, label='Flux bin %d' % h)
    axis.set_title("Number of pixels in the pixel lookup tables")
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    axis.legend()
    
    plt.savefig(globdata.plotpath + 'look.png')
    plt.close()
    
    
def plot_psfn_type():
    
    devi = linspace(0., 5., 100)
    y = zeros((x.size, 5))

    figr, axis = plt.subplots(figsize=(10, 6))
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

    axis.legend()
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylim([1e-3, None])
    plt.show()
    
    
def plot_minmspecinfo(minmspecarry, listinfo, listlevi):
    
    figr, axis = plt.subplots(figsize=(14, 10))
    ax_ = axis.twinx()
    axis.plot(minmspecarry, listinfo, label='Relative entropy')
    axis.legend(bbox_to_anchor=[0.3, 1.08], loc=2)
    
    ax_.plot(minmspecarry, listlevi, label='Log-evidence', color='g')
    ax_.legend(bbox_to_anchor=[0.7, 1.08])

    axis.set_ylabel('$D_{KL}$ [nats]')
    ax_.set_ylabel(r'$\log P(D)$ [nats]')
    axis.set_xlabel('$f_{min}$ [1/cm$^2$/s/GeV]')
    axis.set_xscale('log')
    figr.savefig(os.environ["PNTS_TRAN_DATA_PATH"] + '/png/minmspecinfo.png')
    plt.close(figr)
    
    
def plot_evidtest():
    
    minmgain = -1.
    maxmgain = 5.
    minmdevi = 0.
    maxmdevi = 5.
    gain = linspace(minmgain, maxmgain, 100)
    devi = linspace(minmdevi, maxmdevi, 100)

    evid = log(sqrt(1. + exp(2. * gain[None, :])) *                exp(-devi[:, None]**2 / 2. / (1. + 1. / exp(2. * gain[None, :]))))
    
    figr, axis = plt.subplots(figsize=(7, 7))
    figr.suptitle('Log-Bayesian Evidence For Lower-Dimension Model', fontsize=18)
    imag = axis.imshow(evid, extent=[minmgain, maxmgain, minmdevi, maxmdevi], cmap='winter', origin='lower')
    cset1 = plt.contourf(gain, devi, evid, cmap='winter')
    axis.set_xlabel('Information gain')
    axis.set_ylabel('Goodness of fit')
    plt.colorbar(imag, ax=axis, fraction=0.03)
    #figr.subplots_adjust(top=0.8)

    plt.savefig(globdata.plotpath + 'evidtest_' + rtag + '.png')
    plt.close(figr)
    
    
def plot_pntsprob(globdata, pntsprobcart, ptag, full=False, cumu=False):
    
    if cumu:
        numbcols = 1
    else:
        numbcols = 2
        
    if cumu:
        numbrows = 1
    elif full:
        numbrows = numbspec / 2
    else:
        numbrows = 2
        
    if exprtype == 'ferm':
        strgvarb = '$f$'
        strgunit = ' [1/cm$^2$/s/GeV]'
    if exprtype == 'sdss':
        strgvarb = '$C$'
        strgunit = ' [counts]'
    titl = strgvarb + strgunit

    for l in indxpopl:
        for i in indxenerfdfn:
            figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * 7, numbrows * 7), sharex='all', sharey='all')
            if numbrows == 1:
                axgr = [axgr]            
            figr.suptitle(titl, fontsize=18)
            for a, axrw in enumerate(axgr):
                if numbcols == 1:
                    axrw = [axrw]
                for b, axis in enumerate(axrw):
                    h = a * 2 + b

                    if h < 3 or full:
                        imag = axis.imshow(pntsprobcart[:, :, l, i, h], origin='lower', cmap='Reds',                                            norm=mpl.colors.LogNorm(vmin=0.01, vmax=1), extent=globdata.exttrofi)
                    else:
                        imag = axis.imshow(sum(pntsprobcart[:, :, l, i, 3:], 2), origin='lower', cmap='Reds',                                            norm=mpl.colors.LogNorm(vmin=0.01, vmax=1), extent=globdata.exttrofi)

                    # vmin=0.01, vmax=1
                
                    plt.colorbar(imag, fraction=0.05, ax=axis)

                    # superimpose true PS
                    if trueinfo:
                        if h < 3 or full:
                            indxpnts = where((binsspec[i, h] < truespec[l][0, i, :]) &                                           (truespec[l][0, i, :] < binsspec[i, h+1]))[0]
                        else:
                            indxpnts = where(binsspec[i, 3] < truespec[l][0, i, :])[0]

                        mar1 = axis.scatter(truelgal[l][indxpnts],                                           truebgal[l][indxpnts],                                           s=100, alpha=0.5, marker='x', lw=2, color='g')
                        
                    axis.set_xlabel(longlabl)
                    axis.set_ylabel(latilabl)
                    axis.set_xlim([frambndrmarg, -frambndrmarg])
                    axis.set_ylim([-frambndrmarg, frambndrmarg])
                    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
                    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
                    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
                    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
                    if not cumu:
                        if h < 3 or full:
                            axis.set_title(tdpy_util.util.mexp(binsspec[i, h]) + ' $<$ ' + strgvarb + ' $<$ ' + tdpy_util.util.mexp(binsspec[i, h+1]))
                        else:
                            axis.set_title(tdpy_util.util.mexp(binsspec[i, h]) + ' $<$ ' + strgvarb)
            figr.savefig(globdata.plotpath + 'pntsbind' + ptag + '%d%d' % (l, indxenerincl[i]) + '_' + rtag + '.png')
            plt.close(figr)
       
    
def plot_king(globdata):

    figr, axgr = plt.subplots(1, 2, figsize=(12, 6))
    figr.suptitle('King Function', fontsize=20)
    for k, axis in enumerate(axgr):
        if k == 0:
            sigmlist = [0.25]
            gammlist = [1.01, 2.5, 10.]
            lloc = 3
        else:
            sigmlist = [0.1, 0.25, 1.]
            gammlist = [2.]
            lloc = 1
        for sigm in sigmlist:
            for gamm in gammlist:
                axis.plot(rad2deg(globdata.angldisp), retr_king(sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        axis.legend(loc=lloc)
        axis.set_yscale('log')
        axis.set_xlabel(r'$\theta$ ' + strganglunit)
        axis.set_xlabel(r'$\mathcal{K}(\theta)$')
        plt.figtext(0.7, 0.7, '$\mathcal{K}(\theta) = \frac{1}{2\pi\sigma^2}(1-\frac{1}{\gamma}(\frac{x^2}{2\sigma^2\gamma})^{-\gamma})$')
        
    figr.subplots_adjust()
    plt.savefig(globdata.plotpath + 'king.png')
    plt.close(figr)
    
    
def plot_psfn(globdata, thispsfn):
    
    if exprtype == 'sdss':
        globdata.angldisptemp = rad2deg(globdata.angldisp) * 3600.
    if exprtype == 'ferm':
        globdata.angldisptemp = rad2deg(globdata.angldisp)

    with sns.color_palette("Blues", globdata.numbevtt):

        figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener,                                   figsize=(7 * globdata.numbener, 7 * globdata.numbevtt))
        figr.suptitle(r'Point Spread Function, d$P$/d$\Omega$ [1/sr]', fontsize=20)
        if globdata.numbevtt == 1:
            axgr = [axgr]
        for m, axrw in enumerate(axgr):
            if globdata.numbener == 1:
                axrw = [axrw]
            for i, axis in enumerate(axrw):
                axis.plot(globdata.angldisptemp, thispsfn[i, :, m], label='Sample')
                if trueinfo and datatype == 'mock':
                    axis.plot(globdata.angldisptemp, truepsfn[i, :, m], label='Mock', color='g', ls='--')
                if exprtype == 'ferm':
                    axis.plot(globdata.angldisptemp, fermpsfn[i, :, m], label='Fermi-LAT', color='r', ls='-.', alpha=0.4)
                axis.set_yscale('log')
                if m == globdata.numbevtt - 1:
                    axis.set_xlabel(r'$\theta$ ' + strganglunit)
                if i == 0 and exprtype == 'ferm':
                    axis.set_ylabel(evttstrg[m])
                if m == 0:
                    axis.set_title(enerstrg[i])  
                if i == globdata.numbener - 1 and m == globdata.numbevtt - 1:
                    axis.legend(loc=2)
                indxsamp = indxsamppsfipara[i*nformpara+m*globdata.numbener*nformpara]
                if psfntype == 'singgaus':
                    strg = r'$\sigma = %.3g$ ' % rad2deg(thissampvarb[indxsamp]) + strganglunit
                elif psfntype == 'singking':
                    strg = r'$\sigma = %.3g$ ' % rad2deg(thissampvarb[indxsamp]) + strganglunit + '\n'
                    strg += r'$\gamma = %.3g$' % thissampvarb[indxsamp+1]
                elif psfntype == 'doubgaus':
                    strg = r'$f = %.3g$' % thissampvarb[indxsamp] + '\n'
                    if exprtype == 'sdss':
                        paratemp = rad2deg(thissampvarb[indxsamp+1]) * 3600.
                    if exprtype == 'ferm':
                        paratemp = rad2deg(thissampvarb[indxsamp+1])
                    strg += r'$\sigma = %.3g$ ' % paratemp + strganglunit + '\n'
                    if exprtype == 'sdss':
                        paratemp = rad2deg(thissampvarb[indxsamp+2]) * 3600.
                    if exprtype == 'ferm':
                        paratemp = rad2deg(thissampvarb[indxsamp+2])
                    strg += r'$\sigma = %.3g$ ' % paratemp + strganglunit
                elif psfntype == 'gausking':
                    strg = r'$f_G = %.3g$' % thissampvarb[indxsamp] + '\n'
                    strg += r'$\sigma_G = %.3g$ ' % rad2deg(thissampvarb[indxsamp+1]) + strganglunit + '\n'
                    strg += r'$\sigma_K = %.3g$ ' % rad2deg(thissampvarb[indxsamp+2]) + strganglunit + '\n'
                    strg += r'$\gamma = %.3g$' % thissampvarb[indxsamp+3]
                elif psfntype == 'doubking':
                    strg = r'$f_c = %.3g$' % thissampvarb[indxsamp] + '\n'
                    strg += r'$\sigma_c = %.3g$ ' % rad2deg(thissampvarb[indxsamp+1]) + strganglunit + '\n'
                    strg += r'$\gamma_c = %.3g$' % thissampvarb[indxsamp+2] + '\n'
                    strg += r'$\sigma_t = %.3g$ ' % rad2deg(thissampvarb[indxsamp+3]) + strganglunit + '\n'
                    strg += r'$\gamma_t = %.3g$' % thissampvarb[indxsamp+4]
                axis.text(0.75, 0.75, strg, va='center', ha='center', transform=axis.transAxes, fontsize=18)
                
                if exprtype == 'ferm':
                    axis.set_ylim([1e0, 1e6])
                if exprtype == 'sdss':
                    axis.set_ylim([1e7, 1e11])

        plt.savefig(globdata.plotpath + 'psfnprof_' + rtag + '_%09d.png' % j)
        plt.close(figr)
    
    
def plot_fwhm(globdata, thisfwhm):
    
    figr, axis = plt.subplots()

    tranfwhm = transpose(thisfwhm)
    imag = axis.imshow(rad2deg(tranfwhm), origin='lower', extent=[binsener[0], binsener[-1], 0, 4],                      cmap='BuPu', interpolation='none', alpha=0.4)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.set_xscale('log')
    axis.set_ylabel('PSF Class')
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_yticks([0.5, 1.5, 2.5, 3.5])
    axis.set_yticklabels(['0', '1', '2', '3'])
    axis.set_xticks(binsener)
    axis.set_xticklabels(['%.2g' % binsener[i] for i in indxener])
    axis.set_title('PSF FWHM')
    for i in indxener:
        for m in indxevtt:
            axis.text(meanener[i], indxevttincl[m]+0.5, r'$%.3g^\circ$' % rad2deg(tranfwhm[m, i]), ha='center', va='center', fontsize=14)

    figr.subplots_adjust(bottom=0.2)
    plt.savefig(globdata.plotpath + 'fwhmcnts_' + rtag + '_%09d.png' % j)
    plt.close(figr)
    
    
def plot_backcntsmean(globdata, backcntsmean):
    
    figr, axis = plt.subplots()

    tranumbbackcntsrofimean = transpose(backcntsmean)
    
    imag = axis.imshow(tranumbbackcntsrofimean, origin='lower', extent=[binsener[0], binsener[-1], 0, 4],                      cmap='BuPu', interpolation='none', alpha=0.4)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.set_xscale('log')
    axis.set_ylabel('PSF Class')
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_yticks([0.5, 1.5, 2.5, 3.5])
    axis.set_yticklabels(['0', '1', '2', '3'])
    axis.set_title('Mean FDM counts inside a PSF FWHM')
    for i in indxener:
        for m in indxevtt:
            axis.text(meanener[i], indxevttincl[m]+0.5, '%.3g' % tranumbbackcntsrofimean[m, i], ha='center', va='center')
            
    figr.subplots_adjust(bottom=0.2)
    plt.savefig(globdata.plotpath + 'backcnts_' + rtag + '_%09d.png' % j)
    plt.close(figr)
    
    
def plot_datacntshist(globdata):

    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(7 * globdata.numbener, 7 * globdata.numbevtt))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            datacntstemp = datacnts[i, :, m]
            
            maxmdatacntstemp = amax(datacntstemp)
            minmdatacntstemp = amin(datacntstemp[where(datacntstemp > 0.)])
            binscntstemp = linspace(minmdatacntstemp, maxmdatacntstemp, 20)
            meancntstemp = (binscntstemp[1:] + binscntstemp[0:-1]) / 2.
            diffcntstemp = binscntstemp[1:] - binscntstemp[0:-1]
            
            datacntshist = axis.hist(datacntstemp, binscntstemp, color='b')[0]

            init = [meancntstemp[argmax(datacntshist)], 1.]
            
            global datacntshistglob
            datacntshistglob = diffcntstemp * sum(datacntshist)
    
            axis.set_xlim([amin(binscntstemp), amax(binscntstemp)])
            if m == globdata.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            #axis.set_xscale('log')
            if m == 0:
                axis.set_title(binsenerstrg[i])
            if i == 0 and exprtype == 'ferm':
                axis.set_ylabel(evttstrg[m])
        
    figr.subplots_adjust(wspace=0.3, hspace=0.2)
    plt.savefig(globdata.plotpath + 'datacntshist' + rtag + '.png')
    plt.close(figr)
    
    
def plot_intr():
    
    with plt.xkcd():

        from matplotlib import patheffects
        mpl.rcParams['path.effects'] = [patheffects.withStroke(linewidth=0)]

        figr, axis = plt.subplots(figsize=(14, 6))

        catl = arange(80)
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
        figr.subplots_adjust(bottom=0.15, top=0.9)
        plt.savefig(os.environ["PNTS_TRAN_DATA_PATH"] + '/png/talkintr.png', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_heal(globdata, heal, rofi=True, titl=''):
    
    if rofi:
        healtemp = copy(heal)
        heal = zeros(globdata.numbpixlheal)
        heal[jpixlrofi] = healtemp

    cart = tdpy_util.util.retr_cart(heal, minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    imag = plt.imshow(cart, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.title(titl)
    plt.show()
    
    
def plot_3fgl_thrs():

    expoheal = sum(sum(expoheal, 2)[1:3, :], axis=0)

    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/catl/3fgl_thrs.fits'
    fluxthrs = pf.getdata(path, 0)

    bgalfgl3 = linspace(-90., 90., 481)
    lgalfgl3 = linspace(-180., 180., 960)

    bgalexpo = linspace(-90., 90., 400)
    lgalexpo = linspace(-180., 180., 800)

    fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)

    cntsthrs = fluxthrs * expo

    jbgal = where(abs(bgalexpo) < 10.)[0]
    jlgal = where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)

    axis.set_title('3FGL Detection Flux Threshold [1/cm$^2$/s], 1.0 GeV - 10. GeV')
    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.savefig(globdata.plotpath + 'thrs_' + rtag + '.png')
    plt.close(figr)
    
    
def plot_fgl3(globdata):
    
    figr, axis = plt.subplots()
    bins = logspace(log10(amin(globdata.fgl3timevari[where(globdata.fgl3timevari > 0.)[0]])),                     log10(amax(globdata.fgl3timevari)), 100)
    axis.hist(globdata.fgl3timevari, bins=bins, label='All', log=True)
    axis.hist(globdata.fgl3timevari[globdata.indxfgl3rofi], bins=bins, label='ROI', log=True)
    axis.axvline(72.44, ls='--', alpha=0.5, color='black')
    axis.set_xlabel('3FGL time variability index')
    axis.set_xscale('log')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(globdata.plotpath + 'fgl3timevari.png')
    plt.close(figr)
    

    figr, axis = plt.subplots()
    indxfgl3scut = where(isfinite(globdata.fgl3scut))[0]
    bins = linspace(amin(globdata.fgl3scut[indxfgl3scut]),                     amax(globdata.fgl3scut[indxfgl3scut]), 100)
    axis.hist(globdata.fgl3scut[indxfgl3scut], bins=bins, label='All', log=True)
    axis.hist(globdata.fgl3scut[intersect1d(indxfgl3scut, globdata.indxfgl3rofi)],               bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral cutoff')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(globdata.plotpath + 'fgl3scut.png')
    plt.close(figr)
    
    figr, axis = plt.subplots()
    indxfgl3scur = where(isfinite(globdata.fgl3scur))[0]
    bins = linspace(amin(globdata.fgl3scur[indxfgl3scur]),                     amax(globdata.fgl3scur[indxfgl3scur]), 100)
    axis.hist(globdata.fgl3scur[indxfgl3scur], bins=bins, label='All', log=True)
    axis.hist(globdata.fgl3scur[intersect1d(indxfgl3scur, globdata.indxfgl3rofi)],               bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral curvature')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(globdata.plotpath + 'fgl3scur.png')
    plt.close(figr)
    
    figr, axis = plt.subplots()
    indxfgl3sind = where(isfinite(globdata.fgl3sind))[0]
    bins = linspace(amin(globdata.fgl3sind[indxfgl3sind]),                     amax(globdata.fgl3sind[indxfgl3sind]), 100)
    axis.hist(globdata.fgl3sind[indxfgl3sind], bins=bins, label='All', log=True)
    axis.hist(globdata.fgl3sind[intersect1d(indxfgl3sind, globdata.indxfgl3rofi)],               bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral index')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(globdata.plotpath + 'fgl3sind.png')
    plt.close(figr)
    
    strgfgl3spectype = ['LogParabola', 'PLExpCutoff', 'PLSuperExpCutoff', 'PowerLaw']

def make_anim():

    listname = ['errrcnts0A', 'datacnts0A', 'resicnts0A', 'modlcnts0A', 'histspec',         'scatspec', 'psfnprof', 'compfrac0', 'compfracspec', 'scatpixl']
    
    #print os.listdir(globdata.plotpath)
    for name in listname:
    
        strg = '%s*0.png' % name
        listfile = fnmatch.filter(os.listdir(globdata.plotpath), strg)[int(numbburn/plotperd):]
        
        print fnmatch.filter(os.listdir(globdata.plotpath), strg)
        print listfile

        nfile = len(listfile)
        jfile = choice(arange(nfile), replace=False, size=nfile)

        cmnd = 'convert -delay 20 '
        for k in range(nfile):
            cmnd += '%s ' % listfile[jfile[k]]
        cmnd += ' %s.gif' % name
        os.system(cmnd)

        
def plot_samp(globdata):

    global thisresicnts, errrmodlcnts
    
    thisresicnts = datacnts - thismodlcnts

    thispsfn = retr_psfn(thissampvarb[indxsamppsfipara], indxener, globdata.angldisp, psfntype=modlpsfntype)
    if pixltype == 'cart':
        thispsfn = thispsfn.reshape((globdata.numbener, -1, globdata.numbevtt))
            

    thisfwhm = retr_fwhm(thispsfn)
    
    plot_psfn(thispsfn)
    
    global thisbackcntsmean
    thisbackcntsmean = empty((globdata.numbener, globdata.numbevtt))
    for c in indxback:
        thisbackcntsmean += mean(thissampvarb[indxsampnormback[c, :, None, None]] * backflux[c] * expo *                                     diffener[:, None, None] * pi * thisfwhm[:, None, :]**2 / 4., 1)
    
    thiscnts = []
    for l in indxpopl:
        ppixl = retr_pixl(thissampvarb[thisindxsampbgal[l]], thissampvarb[thisindxsamplgal[l]])
        cntstemp = thissampvarb[thisindxsampspec[l]][:, :, None] * expo[:, ppixl, :] * diffener[:, None, None]
        thiscnts.append(cntstemp)

    global jtruepntsbias, jtruepntsmiss, specmtch
    indxmodl, jtruepntsbias, jtruepntsmiss = pair_catl(truelgal[l],                          truebgal[l],                          truespec[l][0, :, :],                          thissampvarb[thisindxsamplgal[l]],                          thissampvarb[thisindxsampbgal[l]],                          thissampvarb[thisindxsampspec[l]])

    thisspecmtch = thissampvarb[thisindxsampspec[l]][:, indxmodl]
    thisspecmtch[:, jtruepntsmiss] = 0.


    if thissampvarb[indxsampnumbpnts[l]] > 1:
        for l in indxpopl:
            if globdata.colrprio:
                plot_histsind(globdata, l)
            plot_scatpixl(globdata, l)
            if trueinfo:
                plot_scatspec(globdata, l, thisspecmtch=thisspecmtch)
            plot_histspec(globdata, l)
            plot_histcnts(globdata, l, thiscnts)
            plot_compfrac(globdata)

    for i in indxener:
        
        plot_datacnts(globdata, i, None)
        #plot_catl(globdata, i, None, thiscnts)
        plot_modlcnts(globdata, i, None)
        plot_resicnts(globdata, i, None, thisresicnts)

        #for m in indxevtt:
        #    plot_datacnts(globdata, i, m)
        #    plot_catl(globdata, i, m, thiscnts)
        #    plot_modlcnts(globdata, i, m)
        #    plot_resicnts(globdata, i, m, thisresicnts)
        
    #if globdata.numbener == 3:
    #    plot_datacnts(globdata, None, None)
        
    #plot_fwhm(globdata, thisfwhm)
    #plot_backcntsmean(globdata, thisbackcntsmean)
    
    tempsampvarb, tempppixl, tempcnts, temppntsflux,         tempmodlflux, tempmodlcnts = pars_samp(thisindxpntsfull, drmcsamp[:, 0])
    errrmodlcnts = thismodlcnts - tempmodlcnts
    
    for i in indxener:
        plot_errrcnts(globdata, i, None, errrmodlcnts)
    
    
def plot_histcnts(globdata, l, thiscnts):

    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(7 * globdata.numbener, 7 * globdata.numbevtt))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            if trueinfo:
                truehist = axis.hist(truecnts[l][i, :, m], binscnts[i, :], alpha=0.5, color='g', log=True, label=truelabl)
                if datatype == 'mock':
                    axis.hist(fgl3cnts[i, :, m], binscnts[i, :], alpha=0.5, color='red', log=True, label='3FGL')
            axis.hist(thiscnts[l][i, :, m], binscnts[i, :], alpha=0.5, color='b', log=True, label='Sample')
            if m == globdata.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            axis.set_xscale('log')
            if trueinfo:
                axis.set_ylim([0.1, 1e3])
            if m == 0:
                axis.set_title(binsenerstrg[i])
            if i == 0 and exprtype == 'ferm':
                axis.set_ylabel(evttstrg[m])
            if m == globdata.numbevtt / 2 and i == globdata.numbener / 2:
                axis.legend()
        
    figr.subplots_adjust(wspace=0.3)
    plt.savefig(globdata.plotpath + 'histcnts%d_' % l + rtag + '_%09d.png' % j)
    plt.close(figr)
    
def plot_datacnts(globdata, pener, pevtt, nextstat=False):

    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        if pener == None:
            axis.set_title('')
        else:
            axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)

    # plot the model count map
    if pevtt == None:
        if pener == None:
            imag = axis.imshow(sum(datacntscart, axis=3), origin='lower', extent=globdata.exttrofi, interpolation='none')
        else:
            imag = axis.imshow(sum(datacntscart[:, :, pener, :], axis=2), origin='lower',                              interpolation='none', cmap='Reds', extent=globdata.exttrofi)
    else:
        imag = axis.imshow(datacntscart[:, :, pener, pevtt], origin='lower', interpolation='none',                          cmap='Reds', extent=globdata.exttrofi)
    
    if pevtt != None or pener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
    
    # superimpose catalogs
    for l in indxpopl:

        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal[jtruepntsmiss], bgal[jtruepntsmiss], s=mrkrsize[jtruepntsmiss],                        alpha=mrkralph, label=truelabl + ', missed', marker='x', linewidth=2, color='g')
            axis.scatter(lgal[jtruepntsbias], bgal[jtruepntsbias], s=mrkrsize[jtruepntsbias],                        alpha=mrkralph, label=truelabl + ', biased', marker='o', linewidth=2, color='g')
            indxpnts = setdiff1d(arange(truenumbpnts, dtype=int), concatenate((jtruepntsbias, jtruepntsmiss)))
            axis.scatter(lgal[indxpnts], bgal[indxpnts], s=mrkrsize[indxpnts],                        alpha=mrkralph, label=truelabl + ', hit', marker='D', linewidth=2, color='g')
            for l in indxpopl:
                if jtruepntstimevari[l].size > 0:
                    axis.scatter(lgal[jtruepntstimevari[l]], bgal[jtruepntstimevari[l]], s=100,                                label=truelabl + ', variable', marker='*', linewidth=2, color='y')

        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
            
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

    if nextstat:
        
        if False:
            print 'indxprop'
            print indxprop
            print 'reje'
            print reje
        
        if indxprop == indxpropsplt:
            print 'thislgal'
            print thislgal
            print 'thisbgal'
            print thisbgal
            print 'thisspec'
            print thisspec
            print 'nextlgal0'
            print nextlgal0
            print 'nextlgal1'
            print nextlgal1
            print 'nextbgal0'
            print nextbgal0
            print 'nextbgal1'
            print nextbgal1
            print 'nextspec0'
            print nextspec0
            print 'nextspec1'
            print nextspec1
            
        if indxprop == indxpropmerg:
            print 'thislgal0'
            print thislgal0
            print 'thisbgal0'
            print thisbgal0
            print 'thisspec0'
            print thisspec0
            print 'thislgal1'
            print thislgal1
            print 'thisbgal1'
            print thisbgal1
            print 'thisspec1'
            print thisspec1
            print 'nextlgal'
            print nextlgal
            print 'nextbgal'
            print nextbgal
            print 'nextspec'
            print nextspec  
        
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(abs(modispec[pener, k]), pener)
            
            if exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=mrkralph,                        marker='o', linewidth=2, color=colr)
            
            if False:
                print 'modilgal[k]'
                print modilgal[k]
                print 'modibgal[k]'
                print modibgal[k]
                print
                
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == indxpropsplt or indxprop == indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        if pener == None:
            path = globdata.plotpath + 'datacntsAA_' + rtag + '_%09d.png' % j
        else:
            path = globdata.plotpath + 'datacnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = globdata.plotpath + 'datacnts%d%d_' % (pener, indxevttincl[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_modlcnts(globdata, pener, pevtt):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        modlcntstemp = sum(thismodlcnts[pener, :, :], axis=1)
    else:
        modlcntstemp = thismodlcnts[pener, :, pevtt]
    if pixltype == 'heal':
        modlcntstemp = tdpy_util.util.retr_cart(modlcntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                            minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        modlcntstemp = modlcntstemp.reshape((nsidecart, nsidecart)).T
    modlcntstemp[where(modlcntstemp > datacntssatu[pener])] = datacntssatu[pener]
    
    imag = plt.imshow(modlcntstemp, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    
    # superimpose catalogs
    for l in indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = truelgal[l]
            bgal = truebgal[l]
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = globdata.plotpath + 'modlcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = globdata.plotpath + 'modlcnts%d%d_' % (pener, indxevttincl[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_resicnts(globdata, pener, pevtt, resicnts, nextstat=False):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        resicntstemp = sum(resicnts[pener, :, :], axis=1)
    else:
        resicntstemp = resicnts[pener, :, pevtt]
    if pixltype == 'heal':
        resicntstemp = tdpy_util.util.retr_cart(resicntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                            minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        resicntstemp = resicntstemp.reshape((nsidecart, nsidecart))
    resicntstemp[where(resicntstemp > resicntssatu[pener])] = resicntssatu[pener]
    resicntstemp[where(resicntstemp < -resicntssatu[pener])] = -resicntssatu[pener]
    
    imag = axis.imshow(resicntstemp, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    # superimpose catalogs
    for l in indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

        
    if nextstat:
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(abs(modispec[pener, k]), pener)
            
            
            if exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=mrkralph,                        marker='o', linewidth=2, color=colr)

            
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == indxpropsplt or indxprop == indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

            
    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = globdata.plotpath + 'resicnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = globdata.plotpath + 'resicnts%d%d_' % (pener, indxevttincl[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    
    plt.close(figr)
    
    
def plot_errrcnts(globdata, pener, pevtt, errrcntsrofi):

    if pevtt == None:
        errrcntstemp = sum(errrcntsrofi[pener, :, :], axis=1)
    else:
        errrcntstemp = errrcntsrofi[pener, :, pevtt]
    
    if pixltype == 'heal':
        errrcntstemp = tdpy_util.util.retr_cart(errrcntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal,                                            minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        errrcntstemp = errrcntstemp.reshape((nsidecart, nsidecart))
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)

    # plot the error count map
    imag = axis.imshow(errrcntstemp, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    if pevtt == None:
        path = globdata.plotpath + 'errrcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = globdata.plotpath + 'errrcnts%d%d_' % (pener, indxevttincl[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_catl(globdata, pener, pevtt, thiscnts):
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        

    # superimpose catalogs
    for l in indxpopl:
        
        # model catalog
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=300, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=300, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

    
    
    if trueinfo:
        for l in indxpopl:
            numbpnts = int(truenumbpnts[l])
            for a in range(numbpnts):
                if pevtt == None:
                    cnts = sum(truecnts[l][pener, a, :])
                    sigm = sqrt(sum(truesigm[l][pener, a, :]**2))
                else:
                    cnts = truecnts[l][pener, a, pevtt]
                    sigm = truesigm[l][pener, a, pevtt]
                axis.text(truelgal[l][a] + 0.7, truelgal[l][a] - 0.7,                         '%d/%.2f' % (cnts, sigm), color='g', fontsize=13)

    for l in indxpopl:
        numbpnts = int(thissampvarb[indxsampnumbpnts[l]])
        for a in range(numbpnts):
            if pevtt == None:
                cnts = sum(thiscnts[l][pener, a, :])
            else:
                cnts = thiscnts[l][pener, a, pevtt]
            axis.text(thissampvarb[thisindxsamplgal[l][a]] - 0.5, thissampvarb[thisindxsampbgal[l][a]] + 0.3,                     '%d' % cnts, color='b', fontsize=13)
    
    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % mean(thisbackcntsmean[pener, :]), fontsize=18)
    else:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % thisbackcntsmean[pener, pevtt], fontsize=18)
        
    if pevtt == None:
        path = globdata.plotpath + 'catlcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = globdata.plotpath + 'catlcnts%d%d_' % (pener, indxevttincl[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    
    plt.close(figr)

def plot_topo():
    
    figr, axis = plt.subplots()
    datacnts = zeros(globdata.numbpixlheal)
    datacnts[jpixlrofi] = datacnts[0,:,3]
    testcart = tdpy_util.util.retr_cart(datacnts, minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    imag = axis.imshow(testcart, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.scatter(mocklgal, mockbgal, c='g')
    axis.scatter(lgal, bgal, c='b', s=50)
    axis.set_xlim(extent[0:2])
    axis.set_ylim(extent[2:4])
    axis.set_title('PSF3 Mock count map')
    plt.show()
    
    data = linspace(0., 1., numbbins)

    # parameters of interest
    indxpnts = choice(arange(numbpnts), size=3, replace=False)
    jpara = sort(concatenate((indxsampnormback[1, :], indxsamplgal[indxpnts], indxsampbgal[indxpnts], indxsampspec[0, indxpnts])))
    
    
    labl = ['A', '$\psi_1$', '$\phi_1$', '$f_1$', '$\psi_2$', '$\phi_2$', '$f_2$', '$\psi_3$', '$\phi_3$', '$f_3$']
    labllist = list(itertools.combinations(labl, 2))
    indxlist = list(itertools.combinations(jpara, 2))
    ncomb = len(labllist)
    

    llik = zeros((ncomb, numbbins, numbbins))        
    indxsamp = arange(5, npara, dtype=int)
    xlog = ones(ncomb, dtype=bool)
    ylog = ones(ncomb, dtype=bool)
    globdata.exttrofi = zeros((ncomb, 4))
    thiscntr = -1
    for k in range(ncomb):
        
        for i in range(2):
            if indxlist[k][i] == indxsampnormback[1, :]:
                globdata.exttrofi[k, 2 * i] = minmnfdm
                globdata.exttrofi[k, 2 * i + 1] = maxmnfdm
            elif where(indxlist[k][i] == indxsamplgal[indxpnts])[0].size > 0:
                globdata.exttrofi[k, 2 * i] = minmgang
                globdata.exttrofi[k, 2 * i + 1] = maxmgang
            elif where(indxlist[k][i] == indxsampbgal[indxpnts])[0].size > 0:
                globdata.exttrofi[k, 2 * i] = 0.
                globdata.exttrofi[k, 2 * i + 1] = 2. * pi
                if i == 0:
                    xlog[k] = False
                else:
                    ylog[k] = False
            else:
                globdata.exttrofi[k, 2 * i] = hyprmflx[0]
                globdata.exttrofi[k, 2 * i + 1] = maxmflux

            

    for k in range(ncomb):

        # set all the other parameters to random values
        indxsamprand = setdiff1d(indxsamp, array(indxlist[k]))
        samp[indxsamprand] = rand()

        # scan over the two parameters of interest
        for n in range(numbbins):
            for m in range(numbbins):
                samp[indxlist[k][0]] = data[n]
                samp[indxlist[k][1]] = data[m]
                llik[k, n, m] = retr_llik(samp)

        thiscntr = tdpy_util.util.show_prog(k, ncomb, thiscntr)
            
    
    llik *= 1e-7
    numbrows = 15
    numbcols = 3
    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(24, 90))
    figr.suptitle('Likelihood Topology', fontsize=40)
    for x, axrw in enumerate(axgr):
        for y, axis in enumerate(axrw):  
            k = x * numbcols + y
            im = axis.imshow(llik[k, :, :], origin='lower', cmap='Reds', extent=globdata.exttrofi[k,:], aspect='auto')
            figr.colorbar(im, ax=axis, fraction=0.04)
            axis.set_xlabel(labllist[k][0])
            axis.set_ylabel(labllist[k][1])
            if xlog[k]:
                axis.set_xscale('log')
            if ylog[k]:
                axis.set_yscale('log')
                  
    figr.subplots_adjust(top=0.97, hspace=0.4, wspace=0.4)
    plt.show()
    
    
def plot_pntsdiff():

    tempindxpntsfull = []
    tempnumbpnts = 10
    tempindxpntsfull.append(range(tempnumbpnts))
    tempsamp = rand(maxmsampsize)
    tempsamp[indxsampnumbpnts] = array([tempnumbpnts])
    tempsamp[indxsampfdfnnorm] = cdfn_logt(array([tempnumbpnts]), minmfdfnnorm, factfdfnnorm)
    tempsamp[indxsampfdfnslop] = cdfn_atan(array([1.5]), minmfdfnslop, factfdfnslop)
    tempsamp[indxsamppsfipara] = 0.5
    for c in indxback:
        tempsamp[indxsampnormback[c, 0]] = cdfn_logt(array([1.]), minmnormback[c], factnormback[c])
    tempsampvarb, tempppixl, tempcnts,         temppntsflux, tempflux, tempcntsrofi = pars_samp(tempindxpntsfull, tempsamp)


    pntscnts = tempcnts[0, :, 0] * expo[0, :, 0] * apix * diffener[0]
    isotcnts = isotflux[0, :, 0] * expo[0, :, 0] * apix * diffener[0]
    
    totlcnts0 = isotcnts
    totlcnts1 = pntscnts / 1e6 + isotcnts
    
    nrept = 100
    totlcntsheal0 = zeros((nrept, globdata.numbpixlheal))
    totlcntsheal1 = zeros((nrept, globdata.numbpixlheal))
    for k in range(nrept):
        for n in range(jpixlrofi.size):
            totlcntsheal0[k, jpixlrofi[n]] = poisson(totlcnts0[n])
            totlcntsheal1[k, jpixlrofi[n]] = poisson(totlcnts1[n])
          
    maxmcnts = max(amax(totlcntsheal0), amax(totlcntsheal1))
    binscnts = linspace(0., maxmcnts, maxmcnts + 2)
    diffcnts = binscnts[1:] - binscnts[:-1]
    meancnts = (binscnts[1:] + binscnts[:-1]) / 2.
    
    hist0 = empty((nrept, maxmcnts + 1))
    hist1 = empty((nrept, maxmcnts + 1))
    for k in range(nrept):
        hist0[k, :] = histogram(totlcntsheal0[k, jpixlrofi], binscnts)[0].astype(float)
        hist0[k, :] *= 1. / sum(hist0[k, :]) / diffcnts
        hist1[k, :] = histogram(totlcntsheal1[k, jpixlrofi], binscnts)[0].astype(float)
        hist1[k, :] *= 1. / sum(hist1[k, :]) / diffcnts
        
    
    totlcntscart0 = tdpy_util.util.retr_cart(totlcntsheal0[0, :], minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    totlcntscart1 = tdpy_util.util.retr_cart(totlcntsheal1[0, :], minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    
    fig = plt.figure(figsize=(12, 12))
    axis = figr.add_subplot(221)
    imag = axis.imshow(totlcntscart0, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    if exprtype == 'ferm':
        axis.set_xlim([globdata.maxmlgal, globdata.minmlgal])
        axis.set_ylim([globdata.minmbgal, globdata.maxmbgal])
    else:
        axis.set_xlim(array([globdata.maxmlgal, globdata.minmlgal]) * 3600.)
        axis.set_ylim(array([globdata.minmbgal, globdata.maxmbgal]) * 3600.)
    
    axis.set_title('Isotropic')

    
    axis = figr.add_subplot(222)
    imag = axis.imshow(totlcntscart1, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.set_xlabel(r'$x$ [$^\circ$]')
    axis.set_xlim([globdata.maxmlgal, globdata.minmlgal])
    axis.set_ylim([globdata.minmbgal, globdata.maxmbgal])
    axis.set_title('Isotropic + Unresolved PS')
    
    axis.scatter(tempsampvarb[trueindxsamplgal], tempsampvarb[trueindxsampbgal],                s=50, alpha=0.8, marker='x', color='g', linewidth=2)
    
    axis = figr.add_subplot(212)


    tdpy_util.util.plot_braz(ax, meancnts, hist0,  lcol='lightgreen',                         alpha=0.3, dcol='green', mcol='darkgreen', labl='Isotropic')
    tdpy_util.util.plot_braz(ax, meancnts, hist1, lcol='lightblue',                         alpha=0.3, dcol='blue', mcol='darkblue', labl='Isotropic + Unresolved PS')

    axis.set_title('Count PDF')
    axis.set_xlabel('$k$')
    axis.set_ylabel('d$P$/d$k$')
    axis.set_yscale('log')
    axis.legend()

    plt.savefig(globdata.plotpath + 'pntsdiff.png')
    plt.close(figr)
    


