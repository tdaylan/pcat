
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

from util import *

# In[ ]:

def plot_post(pathprobcatl):
         
    hdun = pf.open(pathprobcatl)

    globdata = globdatastrt()
    
    globdata.numbener = hdun[0].header['numbener']
    globdata.numbevtt = hdun[0].header['numbevtt']
    globdata.numbpopl = hdun[0].header['numbpopl']
    globdata.numbpsfipara = hdun[0].header['numbpsfipara']
    globdata.numbformpara = hdun[0].header['numbformpara']
    
    globdata.indxener = arange(globdata.numbener)
    globdata.indxevtt = arange(globdata.numbevtt)
    globdata.indxpopl = arange(globdata.numbpopl)
    
    globdata.numbsamp = hdun[0].header['numbsamp']
    globdata.numbburn = hdun[0].header['numbburn']
    globdata.numbswep = hdun[0].header['numbswep']
    globdata.factthin = hdun[0].header['factthin']
    globdata.numbpopl = hdun[0].header['numbpopl']
    globdata.numbproc = hdun[0].header['numbproc']

    globdata.maxmgang = hdun[0].header['maxmgang']
    globdata.lgalcntr = hdun[0].header['lgalcntr']
    globdata.bgalcntr = hdun[0].header['bgalcntr']
    globdata.numbspec = hdun[0].header['numbspec']

    globdata.minmlgal = hdun[0].header['minmlgal']
    globdata.maxmlgal = hdun[0].header['maxmlgal']
    globdata.minmbgal = hdun[0].header['minmbgal']
    globdata.maxmbgal = hdun[0].header['maxmbgal']

    globdata.datatype = hdun[0].header['datatype']
    globdata.regitype = hdun[0].header['regitype']
    globdata.psfntype = hdun[0].header['psfntype']
    if globdata.datatype == 'mock':
        globdata.mockpsfntype = hdun[0].header['mockpsfntype']
    globdata.exprtype = hdun[0].header['exprtype']
    globdata.pixltype = hdun[0].header['pixltype']

    if globdata.pixltype == 'heal':
        globdata.numbsideheal = hdun[0].header['numbsideheal']
    else:
        globdata.numbsidecart = hdun[0].header['numbsidecart']

    globdata.trueinfo = hdun[0].header['trueinfo']
    globdata.colrprio = hdun[0].header['colrprio']
    
    globdata.margsize = hdun[0].header['margsize']
    globdata.strgtime = hdun[0].header['strgtime']

    globdata.liststrgbackflux = []
    k = 0
    while True:
        try:
            globdata.liststrgbackflux.append(hdun[0].header['strgbackflux%04d' % k])
            k += 1
        except:
            break

    globdata.levi = hdun[0].header['levi']
    globdata.info = hdun[0].header['info']
    
    globdata.indxenerincl = hdun['indxenerincl'].data
    globdata.indxevttincl = hdun['indxevttincl'].data

    if globdata.datatype == 'mock':
        globdata.mocknumbpnts = hdun['mocknumbpnts'].data
        globdata.mockfdfnslop = hdun['mockfdfnslop'].data
        globdata.mocknormback = hdun['mocknormback'].data
        globdata.mockpsfntype = hdun['mockpsfntype'].data

    globdata.maxmnumbpnts = hdun['maxmnumbpnts'].data
        
    listspechist = hdun['spechist'].data
    pntsprob = hdun['pntsprob'].data
    
    listllik = hdun['llik'].data
    listlpri = hdun['lpri'].data
    
    gmrbstat = hdun['gmrbstat'].data
    listmodlcnts = hdun['modlcnts'].data
    
    # prior boundaries
    if False:
        globdata.minmfdfnnorm = hdun['minmfdfnnorm'].data
        globdata.maxmfdfnnorm = hdun['maxmfdfnnorm'].data
        globdata.minmfdfnslop = hdun['minmfdfnslop'].data
        globdata.maxmfdfnslop = hdun['maxmfdfnslop'].data
    else:
        globdata.minmfdfnnorm = 1.
        globdata.maxmfdfnnorm = 100.
        globdata.minmfdfnslop = 1.5
        globdata.maxmfdfnslop = 2.5

    if globdata.colrprio:
        globdata.minmsind = hdun['minmsind'].data
        globdata.maxmsind = hdun['maxmsind'].data
        globdata.binssind = hdun['binssind'].data
        globdata.meansind = hdun['meansind'].data
        
    globdata.minmspec = hdun['minmspec'].data
    globdata.maxmspec = hdun['maxmspec'].data
    globdata.binsspec = hdun['binsspec'].data
    globdata.meanspec = hdun['meanspec'].data
    
    # bins
    globdata.binsener = hdun['binsener'].data
    globdata.diffener = hdun['diffener'].data
    globdata.meanener = hdun['meanener'].data
    globdata.indxenerfdfn = hdun['indxenerfdfn'].data
    
    # utilities
    globdata.probprop = hdun['probprop'].data
    globdata.listindxprop = hdun['indxprop'].data
    globdata.listchro = hdun['chro'].data
    globdata.listaccp = hdun['accp'].data
    globdata.listindxsampmodi = hdun['sampmodi'].data
    globdata.listauxipara = hdun['auxipara'].data
    globdata.listlaccfrac = hdun['laccfrac'].data
    globdata.listnumbpair = hdun['numbpair'].data
    globdata.listjcbnfact = hdun['jcbnfact'].data
    globdata.listcombfact = hdun['combfact'].data
    globdata.listpntsfluxmean = hdun['listpntsfluxmean'].data

    # posterior distributions
    globdata.listnumbpnts = hdun['numbpnts'].data
    globdata.listfdfnnorm = hdun['fdfnnorm'].data
    globdata.listfdfnslop = hdun['fdfnslop'].data
    globdata.listpsfipara = hdun['psfipara'].data
    globdata.listnormback = hdun['normback'].data


    # setup the sampler
    init(globdata) 

    # truth globdata.information
    if globdata.trueinfo:
        globdata.truenumbpnts = hdun['truenumbpnts'].data
        globdata.truelgal = []
        globdata.truebgal = []
        globdata.truespec = []
        globdata.truesind = []
        globdata.truespec = []
        for l in globdata.indxpopl:
            globdata.truelgal.append(hdun['truelgalpopl%d' % l].data)
            globdata.truebgal.append(hdun['truebgalpopl%d' % l].data)
            globdata.truespec.append(hdun['truespecpopl%d' % l].data)
            if globdata.colrprio:
                globdata.truespec.append(hdun['truesindpopl%d' % l].data)
        globdata.indxtruepntstimevari = hdun['indxtruepntstimevari'].data
        globdata.truefdfnslop = hdun['truefdfnslop'].data
        globdata.truenormback = hdun['truenormback'].data
        globdata.truepsfipara = hdun['truepsfipara'].data
        
    listlgal = []
    listbgal = []
    listspec = []
    if globdata.colrprio:
        listsind = []
    for l in globdata.indxpopl:
        listlgal.append(hdun['lgalpopl%d' % l].data)
        listbgal.append(hdun['bgalpopl%d' % l].data)
        listspec.append(hdun['specpopl%d' % l].data)
        if globdata.colrprio:
            listsind.append(hdun['sindpopl%d' % l].data)

    #globdata.rtag = retr_rtag(globdata, None)

    #if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
    #    plotfold = '/n/pan/www/tansu/png/pcat/'
    #else:
    #    plotfold = os.environ["PCAT_DATA_PATH"] + '/png/'
    #globdata.plotpath = plotfold + globdata.strgtime + '_' + globdata.rtag + '/'
    #cmnd = 'mkdir -p ' + globdata.plotpath
    #os.system(cmnd)

    # get data cube indices
    #retr_indxcube(globdata)

    # load exposure map
    #retr_expo(globdata)

    # get strings            
    #retr_strgangl(globdata)

    # PSF model
    #retr_psfimodl(globdata)

    # proposals
    #globdata.probprop = None
    #retr_propmodl(globdata)
    #numbprop = len(globdata.probprop)
    #globdata.maxmgangmarg = globdata.maxmgang + globdata.margsize
    #globdata.exttrofi = [globdata.minmlgal, globdata.maxmlgal, globdata.minmbgal, globdata.maxmbgal]
    #if globdata.exprtype == 'sdss':
    #    globdata.exttrofi *= 3600.
    #    globdata.frambndr = globdata.maxmgang * 3600.
    #    globdata.frambndrmarg = globdata.maxmgangmarg * 3600.
    #else:
    #    globdata.frambndr = globdata.maxmgang
    #    globdata.frambndrmarg = globdata.maxmgangmarg
        
    #globdata.strgfluxunit = retr_strgfluxunit(globdata)

    # energy bin string
    #globdata.enerstrg, globdata.binsenerstrg = retr_enerstrg(globdata)
    
    #if globdata.pixltype == 'heal':
    #    retr_pixlcnvt(globdata)

        
    # Gelman-Rubin test
    if globdata.numbproc > 1:
        print 'Making the Gelman-Rubin TS plot...'
        tim0 = time.time()
        
        figr, axis = plt.subplots()
        axis.hist(gmrbstat, bins=linspace(1., amax(gmrbstat), 40))
        axis.set_title('Gelman-Rubin Convergence Test')
        axis.set_xlabel('PSRF')
        axis.set_ylabel('$N_{pix}$')
        figr.savefig(globdata.plotpath + 'gmrbdist_' + globdata.rtag + '.png')
        plt.close(figr)
        
        tim1 = time.time()
        print 'Done in %.3g seconds.' % (tim1 - tim0)

        
    print 'Calculating proposal and acceptance rates...'
    tim0 = time.time()

    mcmctimebins = linspace(0., globdata.numbswep, 100)
    figr, axgr = plt.subplots(numbprop, 1, figsize=(10, 15), sharex='all')
    for g, axis in enumerate(axgr):
        axis.hist(where(globdata.listindxprop == g)[0], bins=mcmctimebins)
        axis.hist(where((globdata.listindxprop == g) & (globdata.listaccp == True))[0], bins=mcmctimebins)
        axis.set_ylabel('%d' % g)
        if g == numbprop - 1:
            axis.set_xlabel('$i_{samp}$')
    figr.savefig(globdata.plotpath + 'propeffi_' + globdata.rtag + '.png')
    figr.subplots_adjust(hspace=0.)
    plt.close(figr)
     

    indxsampsplt = where(globdata.listindxprop == globdata.indxpropsplt)[0]
    indxsampmerg = where(globdata.listindxprop == globdata.indxpropmerg)[0]
            
    listname = ['laccfrac', 'numbpair', 'combfact', 'jcbnfact']
    listvarb = [globdata.listlaccfrac, globdata.listnumbpair, globdata.listcombfact, globdata.listjcbnfact]
    for k in range(4):
        figr, axis = plt.subplots()
        axis.hist(listvarb[k][indxsampsplt])
        axis.hist(listvarb[k][indxsampmerg])
        axis.set_ylabel('$N_{samp}$')
        axis.set_xlabel(listname[k])
        figr.subplots_adjust(bottom=0.2)
        figr.savefig(globdata.plotpath + listname[k] + globdata.rtag + '.png')
        plt.close(figr)
    
    
    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Calculating proposal execution times...'
    tim0 = time.time()

    retr_strgprop(globdata)
    
    binstime = logspace(log10(amin(globdata.listchro[where(globdata.listchro > 0)] * 1e3)), log10(amax(globdata.listchro * 1e3)), 50)
    with sns.color_palette("nipy_spectral", numbprop):
        figr, axis = plt.subplots(figsize=(14, 12))
        
        axis.hist(globdata.listchro[where(globdata.listchro[:, 0] > 0)[0], 0] * 1e3, binstime, facecolor='none', log=True, lw=1, ls='--', edgecolor='black')
        for g in range(numbprop):
            indxlistchro = where((globdata.listindxprop == g) & (globdata.listchro[:, 0] > 0))[0]
            axis.hist(globdata.listchro[indxlistchro, 0] * 1e3, binstime, edgecolor='none', log=True, alpha=0.3, label=globdata.strgprop[g])
        axis.set_xlabel('$t$ [ms]')
        axis.set_xscale('log')
        axis.set_xlim([amin(binstime), amax(binstime)])
        axis.set_ylim([0.5, None])
        axis.legend(loc=2)
        figr.savefig(globdata.plotpath + 'chroprop_' + globdata.rtag + '.png')
        plt.close(figr)

    labl = ['Total', 'Proposal', 'Prior', 'Likelihood']
    figr, axcl = plt.subplots(2, 1, figsize=(14, 10))
    for k in range(1, 4):
        axcl[0].hist(globdata.listchro[where(globdata.listchro[:, k] > 0)[0], k] * 1e3, binstime, log=True, alpha=0.5, label=labl[k])
    axcl[1].hist(globdata.listchro[where(globdata.listchro[:, 0] > 0)[0], 0] * 1e3, binstime, log=True, label=labl[0], color='black')
    axcl[1].set_title(r'$\langle t \rangle$ = %.3g ms' % mean(globdata.listchro[where(globdata.listchro[:, 0] > 0)[0], 0] * 1e3))
    axcl[0].set_xlim([amin(binstime), amax(binstime)])
    axcl[1].set_xlabel('$t$ [ms]')
    axcl[0].set_xscale('log')
    axcl[1].set_xscale('log')
    axcl[0].set_ylim([0.5, None])
    axcl[1].set_ylim([0.5, None])
    axcl[0].legend(loc=1)
    axcl[1].legend(loc=2)
    figr.savefig(globdata.plotpath + 'chrototl_' + globdata.rtag + '.png')
    plt.close(figr)

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Parsing the sample bundle and making frames...'
    tim0 = time.time()

    
    if False:
        if maxmnumbpnts[0] == 4 and globdata.numbener == 1 and globdata.numbpopl == 1:
            for k in range(maxmnumbpnts[0]):
                listpostdist = zeros((globdata.numbsamp, 3 * k)) 
                listpostdist[:, 0*k:1*k] = listlgalpnts[0][k]
                listpostdist[:, 1*k:2*k] = listbgalpnts[0][k]
                listpostdist[:, 2*k:3*k] = listspecpnts[0][k][0, :]
                path = globdata.plotpath + 'postdist_%d_' % k + globdata.rtag
                if globdata.trueinfo:
                    truepara = truesampvarb[globdata.indxsamppsfipara[ipsfipara]]
                else:
                    truepara = None
                tdpy.util.plot_mcmc(listpostdist, parastrgpsfipara[ipsfipara],                                     truepara=truepara, path=path, numbbins=globdata.numbbins, quan=True)
                
                
            
    # flux match with the true catalog
    if globdata.trueinfo:
        for l in globdata.indxpopl:
            
            discspecmtch = zeros(globdata.truenumbpnts) + globdata.numbsamp
            listindxmodl = []
            for k in range(globdata.numbsamp):
                indxmodl, indxtruepntsbias, indxtruepntsmiss = pair_catl(globdata, l,                                      listlgal[l][k, :],                                      listbgal[l][k, :],                                      listspec[l][k, :, :])
                listindxmodl.append(indxmodl)
                discspecmtch[indxtruepntsmiss] -= 1.
            discspecmtch /= globdata.numbsamp

            
            postspecmtch = zeros((3, globdata.numbener, globdata.truenumbpnts[l]))
            for i in globdata.indxener:
                listspecmtch = zeros((globdata.numbsamp, globdata.truenumbpnts[l]))
                for k in range(globdata.numbsamp):
                    indxpntstrue = where(listindxmodl[k] >= 0)[0]
                    listspecmtch[k, indxpntstrue] = listspec[l][k][i, listindxmodl[k][indxpntstrue]]
                postspecmtch[0, i, :] = percentile(listspecmtch, 16., axis=0)
                postspecmtch[1, i, :] = percentile(listspecmtch, 50., axis=0)
                postspecmtch[2, i, :] = percentile(listspecmtch, 84., axis=0)
            
            plot_scatspec(globdata, l, postspecmtch=postspecmtch)

            # store the comparison
            path = os.environ["PCAT_DATA_PATH"] + '/pcatcomp_popl%d_' % l + globdata.rtag + '.fits'
            compbund = stack((globdata.truespec[l], postspecmtch))

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Making posterior distribution plots...'

    if globdata.pixltype == 'heal':
        
        reso = 60. / globdata.numbsideheal
        numbbinslgcr = int((globdata.maxmlgal - globdata.minmlgal) / reso)
        numbbinsbgcr = int((globdata.maxmbgal - globdata.minmbgal) / reso)

        pntsprobcart = zeros((numbbinslgcr, numbbinsbgcr, globdata.numbpopl, globdata.numbener, globdata.numbspec))
        for l in globdata.indxpopl:
            for i in globdata.indxener:
                for h in range(globdata.numbspec):
                    pntsprobcart[:, :, l, i, h] = tdpy.util.retr_cart(pntsprob[l, i, :, h], 
                                                                      indxpixlrofi=globdata.indxpixlrofi, \
                                                                      numbsideinpt=globdata.numbsideheal, \
                                                                      minmlgal=globdata.minmlgal, \
                                                                      maxmlgal=globdata.maxmlgal, \
                                                                      minmbgal=globdata.minmbgal, \
                                                                      maxmbgal=globdata.maxmbgal, \
                                                                      reso=reso)
    else:
        pntsprobcart = pntsprob.reshape((globdata.numbpopl, globdata.numbener, globdata.numbsidecart, globdata.numbsidecart, globdata.numbspec))
        pntsprobcart = swapaxes(swapaxes(pntsprobcart, 0, 2), 1, 3)
        
    # stacked posteiors binned in position and flux
    plot_pntsprob(globdata, pntsprobcart, ptag='quad')
    plot_pntsprob(globdata, pntsprobcart, ptag='full', full=True)
    plot_pntsprob(globdata, pntsprobcart, ptag='cumu', cumu=True)
  
                                
    if globdata.exprtype == 'ferm':
        retr_fgl3(globdata)

    # flux distribution
    for l in globdata.indxpopl:
        plot_histspec(globdata, l, listspechist=listspechist[:, l, :, :])

    # fraction of emission components
    postpntsfluxmean = retr_postvarb(globdata.listpntsfluxmean)
    postnormback = retr_postvarb(globdata.listnormback)
    plot_compfrac(globdata, postpntsfluxmean=postpntsfluxmean, postnormback=postnormback)


    # PSF parameters
    path = globdata.plotpath + 'psfipara_' + globdata.rtag
    if globdata.psfntype == 'singgaus' or globdata.psfntype == 'singking':
        globdata.listpsfipara[:, jpsfiparainit] = rad2deg(globdata.listpsfipara[:, jpsfiparainit])
        if globdata.trueinfo:
            globdata.truepsfipara[jpsfiparainit] = rad2deg(globdata.truepsfipara[jpsfiparainit])
    elif globdata.psfntype == 'doubgaus' or globdata.psfntype == 'gausking':
        globdata.listpsfipara[:, jpsfiparainit+1] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+1])
        globdata.listpsfipara[:, jpsfiparainit+2] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+2])
        if globdata.trueinfo:
            globdata.truepsfipara[jpsfiparainit+1] = rad2deg(globdata.truepsfipara[jpsfiparainit+1])
            globdata.truepsfipara[jpsfiparainit+2] = rad2deg(globdata.truepsfipara[jpsfiparainit+2])
    elif globdata.psfntype == 'doubking':
        globdata.listpsfipara[:, jpsfiparainit+1] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+1])
        globdata.listpsfipara[:, jpsfiparainit+3] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+3])
        if globdata.trueinfo:
            globdata.truepsfipara[jpsfiparainit+1] = rad2deg(globdata.truepsfipara[jpsfiparainit+1])
            globdata.truepsfipara[jpsfiparainit+3] = rad2deg(globdata.truepsfipara[jpsfiparainit+3])
    if globdata.trueinfo and psfntype == 'doubking':
        truepara = globdata.truepsfipara
    else:
        truepara = array([None] * globdata.numbpsfipara)
    tdpy.util.plot_mcmc(globdata.listpsfipara, strgpsfipara, truepara=truepara,                         nplot=globdata.numbformpara, path=path, numbbins=globdata.numbbins, quan=True, ntickbins=3)
    
    for k in range(globdata.numbpsfipara):
        path = globdata.plotpath + 'psfipara%d_' % k + globdata.rtag + '.png'
        tdpy.util.plot_trac(globdata.listpsfipara[:, k], strgpsfipara[k], path=path, quan=True)
    
    
    # log-likelihood
    path = globdata.plotpath + 'llik_' + globdata.rtag + '.png'
    tdpy.util.plot_trac(listllik.flatten(), '$P(D|y)$', path=path)

    # log-prior
    path = globdata.plotpath + 'lpri_' + globdata.rtag + '.png'
    tdpy.util.plot_trac(listlpri.flatten(), '$P(y)$', path=path)
    

    # number, expected number of PS and flux conditional prior power law index 
    for l in range(globdata.numbpopl):
        
        # number of point sources
        path = globdata.plotpath + 'numbpntsdist_popl%d_' % l + globdata.rtag + '.png'
        if globdata.trueinfo and globdata.truenumbpnts != None:
            truepara = globdata.truenumbpnts[l]
        else:
            truepara = None
        tdpy.util.plot_trac(globdata.listnumbpnts[:, l], '$N$', truepara=truepara, path=path)

        # mean number of point sources
        path = globdata.plotpath + 'fdfnnorm_popl%d_' % l + globdata.rtag + '.png'
        truepara = None
        tdpy.util.plot_trac(globdata.listfdfnnorm[:, l], '$\mu$', truepara=truepara, path=path)

        # flux distribution power law index
        for i in globdata.indxenerfdfn:
            path = globdata.plotpath + 'fdfnslopdist_popl%d%d_' % (l, i) + globdata.rtag + '.png'
            if globdata.trueinfo and globdata.truefdfnslop != None:
                truepara = globdata.truefdfnslop[l, i]
            else:
                truepara = None
            titl = globdata.binsenerstrg[i]
            labl =  r'$\alpha_{%d}$' % i
            tdpy.util.plot_trac(globdata.listfdfnslop[:, l, i], labl, truepara=truepara, path=path, titl=titl)
        
        
    # isotropic background normalization
    for i in globdata.indxener:
        path = globdata.plotpath + 'nisodist%d_' % i + globdata.rtag + '.png'
        if globdata.trueinfo:
            if globdata.datatype == 'mock':
                truepara = globdata.truenormback[0, i]
            else:
                truepara = None
        titl = globdata.binsenerstrg[i]
        labl = r'$\mathcal{I}_{%d}$' % i
        tdpy.util.plot_trac(globdata.listnormback[:, 0, i], labl, truepara=truepara, path=path, titl=titl)
       
    if globdata.exprtype == 'ferm':
        # diffuse model normalization
        for i in globdata.indxener:
            path = globdata.plotpath + 'nfdmdist%d_' % i + globdata.rtag + '.png'
            if globdata.trueinfo:
                if globdata.datatype == 'mock':
                    truepara = globdata.truenormback[1, i]
                else:
                    truepara = None
            else:
                truepara = None
            titl = globdata.binsenerstrg[i]
            labl = r'$\mathcal{D}_{%d}$' % i
            tdpy.util.plot_trac(globdata.listnormback[:, 1, i], labl, truepara=truepara, path=path, titl=titl)

    # plot log-likelihood
    figr, axrw = plt.subplots(2, 1, figsize=(7, 12))
    figr.suptitle(r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (globdata.info, globdata.levi))
    for k, axis in enumerate(axrw):
        if k == 0:
            if amin(listllik) != amax(listllik):
                axis.hist(listllik.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(D|x)$')
        else:
            if amin(listlpri) != amax(listllik):
                axis.hist(listlpri.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(x)$')
    figr.savefig(globdata.plotpath + 'leviinfo_' + globdata.rtag + '.png')
    plt.close(figr)

    #make_anim()

    tim1 = time.time()
    print 'Plots are produced in %.3g seconds.' % (tim1 - tim0)

def plot_compfrac(globdata, postpntsfluxmean=None, postnormback=None):
    
    if postpntsfluxmean != None:
        post = True
    else:
        post = False
        
        
    listlinestyl = ['-', '--', '-.', ':']
    listcolr = ['black', 'b', 'b', 'b']
    listlabl = ['Data', 'PS', 'Iso', 'FDM']

    figr, axis = plt.subplots(figsize=(7 * globdata.numbener, 7))
    
    listydat = empty((globdata.numbback + 2, globdata.numbener))
    listyerr = zeros((2, globdata.numbback + 2, globdata.numbener))
    
    listydat[0, :] = globdata.datafluxmean
    if post:
        listydat[1, :] = postpntsfluxmean[0, :]
        listyerr[:, 1, :] = retr_errrvarb(postpntsfluxmean)
        for c in globdata.indxback:
            listydat[c+2, :] = postnormback[0, c, :] * globdata.backfluxmean[c]
            listyerr[:, c+2, :] = retr_errrvarb(postnormback[:, c, :]) * globdata.backfluxmean[c]
    else:
        listydat[1, :] = mean(sum(globdata.thispntsflux * globdata.expo, 2) / sum(globdata.expo, 2), 1)
        for c in globdata.indxback:
            listydat[c+2, :] = globdata.thissampvarb[globdata.indxsampnormback[c, :]] * globdata.backfluxmean[c]

    
    xdat = globdata.meanener
    for k in range(globdata.numbback + 2):
        ydat = globdata.meanener**2 * listydat[k, :]
        yerr = globdata.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5,                     ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    # Fermi-LAT results
    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            
            pass
            
        else:
            
            if globdata.exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PCAT_DATA_PATH"] + '/fermspec' + name + '.csv'
                    data = loadtxt(path)
                    enertemp = data[:, 0] # [GeV]
                    fluxtemp = data[:, 1] * 1e-3 # [GeV/cm^2/s/sr]
                    fluxtemp = interp(globdata.meanener, enertemp, fluxtemp)
                    #fluxtemp = interpolate.interp1d(enertemp, fluxtemp)(globdata.meanener)
                    axis.plot(globdata.meanener, fluxtemp, marker=listmrkr[k], color=listcolr[k], label=listlabl[k])


    axis.set_xlim([amin(globdata.binsener), amax(globdata.binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [GeV]')
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [GeV/cm$^2$/s/sr]')
    axis.legend()

    if post:
        path = globdata.plotpath + 'compfracspec_' + globdata.rtag + '.png'
    else:
        path = globdata.plotpath + 'compfracspec_' + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    
    listlabl = ['PS', 'Iso', 'FDM']
    listexpl = [0.1, 0, 0]

    listsize = zeros(globdata.numbback + 1)
    for k in range(globdata.numbback + 1):
        if globdata.numbener == 1:
            listsize[k] = diffener * listydat[k+1, :]
        else:
            listsize[k] = trapz(listydat[k+1, :], globdata.meanener)
    listsize *= 100. / sum(listsize)
        
    figr, axis = plt.subplots()

    axis.pie(listsize, explode=listexpl, labels=listlabl, autopct='%1.1f%%')
    axis.axis('equal')

    if post:
        path = globdata.plotpath + 'compfrac_' + globdata.rtag + '.png'
    else:
        path = globdata.plotpath + 'compfrac_' + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
     

def plot_histsind(globdata, l, postsindhist=None):
    
    if postsindhist == None:
        post = False
    else:
        post = True
        
    figr, axis = plt.subplots()
    if post:
        xdat = globdata.meansind[i, :]
        ydat = postsindhist[0, :]
        yerr = retr_errrvarb(postsindhist)
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    else:
        axis.hist(globdata.thissampvarb[thisindxsampsind[l]], globdata.binssind, alpha=0.5, color='b', log=True, label='Sample')
    if globdata.trueinfo:
        axis.hist(globdata.truesind[l], globdata.binssind, alpha=0.5, color='g', log=True, label=globdata.truelabl)
        if globdata.datatype == 'mock':
            axis.hist(fgl3sind, globdata.binssind, alpha=0.1, color='red', log=True, label='3FGL')
    axis.set_yscale('log')
    axis.set_xlabel('$s$')
    axis.set_xlim([globdata.minmsind, globdata.maxmsind])
    axis.set_ylabel('$N$')
    if post:
        path = globdata.plotpath + 'histsind_popl%d' % l + globdata.rtag + '.png'
    else:
        path = globdata.plotpath + 'histsind_popl%d' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    

def plot_histspec(globdata, l, listspechist=None):
    
    if listspechist == None:
        post = False
    else:
        post = True
        
    figr, axcl = plt.subplots(1, globdata.numbener, figsize=(7 * globdata.numbener, 7))
    if globdata.numbener == 1:
        axcl = [axcl]
    for i, axis in enumerate(axcl):
        if post:
            xdat = globdata.meanspec[i, :]
            yerr = empty((2, globdata.numbspec))
            ydat = percentile(listspechist[:, i, :], 50., axis=0)
            yerr[0, :] = percentile(listspechist[:, i, :], 16., axis=0)
            yerr[1, :] = percentile(listspechist[:, i, :], 84., axis=0)
            yerr = abs(yerr - ydat)
            axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        else:
            axis.hist(globdata.thissampvarb[globdata.thisindxsampspec[l]][i, :], globdata.binsspec[i, :], alpha=0.5, color='b',                     log=True, label='Sample')
            
            if not globdata.colrprio or i == globdata.indxenerfdfn:
                fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[l, i]]  
                fluxhistmodl = retr_fdfn(globdata, globdata.thissampvarb[globdata.indxsampfdfnnorm[l]], fdfnslop, i)
                axis.plot(globdata.meanspec[i, :], fluxhistmodl, ls='--', alpha=0.5, color='b')
            
        if globdata.trueinfo:
            truehist = axis.hist(globdata.truespec[l][0, i, :], globdata.binsspec[i, :],                                alpha=0.5, color='g', log=True, label=globdata.truelabl)
            if globdata.datatype == 'mock':
                axis.hist(globdata.fgl3spec[0, i, :], globdata.binsspec[i, :], color='red', alpha=0.1, log=True, label='3FGL')

        axis.set_yscale('log')
        axis.set_xlabel('$f$ ' + globdata.strgfluxunit)
        axis.set_xscale('log')
        axis.set_title(globdata.enerstrg[i])
        if globdata.trueinfo:
            axis.set_ylim([0.1, 1e3])
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        if i == 0:
            axis.set_ylabel('$N$')
        if i == globdata.numbener / 2:
            axis.legend()
        
    figr.subplots_adjust(wspace=0.3, bottom=0.2)
    
    if post:
        path = globdata.plotpath + 'histspec%d_' % l + globdata.rtag + '.png'
    else:
        path = globdata.plotpath + 'histspec%d_' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    

def plot_scatspec(globdata, l, postspecmtch=None, thisspecmtch=None):
    
    figr, axrw = plt.subplots(1, globdata.numbener, figsize=(7 * globdata.numbener, 6))
    if globdata.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        xdat = globdata.truespec[l][0, i, :]
        xerr = retr_errrvarb(globdata.truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))

 
        if postspecmtch != None or thisspecmtch != None:
            labl = '$f_{samp}$ ' + globdata.strgfluxunit
            axis.plot(globdata.meanspec[i, :], globdata.meanspec[i, :], ls='--', color='black', alpha=0.2)
            if postspecmtch != None:
                yerr[0, :] = postspecmtch[0, i, :]
                ydat = postspecmtch[1, i, :]
                yerr[1, :] = postspecmtch[2, i, :]
                yerr = abs(yerr - ydat)
            else:
                ydat = thisspecmtch[i, :]

        indx = where(ydat > 0.)[0]
        if indx.size > 0:
            axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], lw=1, marker='o', markersize=5, color='black')
    
        if globdata.indxtruepntstimevari[l].size > 0:
            axis.errorbar(xdat[globdata.indxtruepntstimevari[l]], ydat[globdata.indxtruepntstimevari[l]], ls='', yerr=yerr[:, globdata.indxtruepntstimevari[l]], lw=1, marker='o', markersize=5, color='red')
    
        axis.set_xlabel('$f_{true}$ ' + globdata.strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_title(globdata.enerstrg[i])

        ylim = [globdata.minmspec[i], globdata.maxmspec[i]]
        
        axis.set_ylim(ylim)
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        axis.set_title(globdata.binsenerstrg[i])

    figr.subplots_adjust(wspace=0.4, bottom=0.2)

    if postspecmtch != None:
        path = globdata.plotpath + 'scatspec%d_' % l + globdata.rtag + '.png'
    elif thisspecmtch != None:
        path = globdata.plotpath + 'scatspec%d_' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep

    plt.savefig(path)
    plt.close(figr)


def plot_scatpixl(globdata, l):
    
    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(globdata.numbener * 7, globdata.numbevtt * 7))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            slop, intc, coef, pval, stdv = sp.stats.linregress(globdata.datacnts[i, :, m], globdata.thismodlcnts[i, :, m])
            axis.scatter(globdata.datacnts[i, :, m], globdata.thismodlcnts[i, :, m], alpha=0.5)

            axislimt = [0., amax(globdata.datacnts[i, :, m]) * 1.5]
            axis.set_xlim(axislimt)
            axis.set_ylim(axislimt)
            
            axis.text(0.5, 0.9, '$m = %.3g, b = %.3g, r = %.3g$' % (slop, intc, coef),                     va='center', ha='center', transform=axis.transAxes, fontsize=16)
            if m == globdata.numbevtt - 1:
                axis.set_xlabel('Data Counts')
            if i == 0:
                labl = 'Model Counts'
                if globdata.exprtype == 'ferm':
                    labl += ', ' + globdata.evttstrg[m]
                axis.set_ylabel(labl)
            if m == 0:
                axis.set_title(globdata.enerstrg[i])


            
    figr.subplots_adjust(hspace=0.4, wspace=0.4, top=0.8)
    plt.savefig(globdata.plotpath + 'scatpixl%d_' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep)
    plt.close(figr)
    
    
def plot_discspec(globdata, l, discspecmtch=None):
    
    figr, axrw = plt.subplots(1, globdata.numbener, figsize=(7 * globdata.numbener, 6))
    if globdata.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):
        xdat = globdata.truespec[l][0, i, :]
        xerr = retr_errrvarb(globdata.truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
        if discspecmtch != None:
            ydat = discspecmtch[i, :]
            labl = '$f_{hit}$'
            ylim = [0., 1.]
        if postspecmtch != None or thisspecmtch != None:
            labl = '$f_{samp}$ ' + globdata.strgfluxunit
            axis.set_yscale('log')
            axis.plot(globdata.meanspec[i, :], globdata.meanspec[i, :], ls='--', color='black', alpha=0.2)
            if postspecmtch != None:
                yerr[0, :] = postspecmtch[0, i, :]
                ydat = postspecmtch[1, i, :]
                yerr[1, :] = postspecmtch[2, i, :]
                yerr = abs(yerr - ydat)
            else:
                ydat = thisspecmtch[i, :]
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        axis.set_xlabel('$f_{true}$ ' + globdata.strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_title(globdata.enerstrg[i])
        if globdata.colrprio:
            if i == globdata.indxenerfdfn[0]:
                ylim = [globdata.minmspec[i], globdata.maxmspec[i]]
        else:
            ylim = [globdata.minmspec[i], globdata.maxmspec[i]]      
        axis.set_ylim(ylim)
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        axis.set_title(globdata.binsenerstrg[i])
    figr.subplots_adjust(wspace=0.4, bottom=0.2)
    path = globdata.plotpath + 'discspec%d_' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)

    
def plot_look(globdata):

    indxpixlproxsize = zeros(globdata.numbpixl)
    figr, axis = plt.subplots(figsize=(10, 6))
    for h in range(globdata.numbspecprox):
        for j in globdata.indxpixl:
            indxpixlproxsize[j] = globdata.indxpixlprox[h][j].size
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
    figr.savefig(os.environ["PCAT_DATA_PATH"] + '/png/minmspecinfo.png')
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

    plt.savefig(globdata.plotpath + 'evidtest_' + globdata.rtag + '.png')
    plt.close(figr)
    
    
def plot_pntsprob(globdata, pntsprobcart, ptag, full=False, cumu=False):
    
    if cumu:
        numbcols = 1
    else:
        numbcols = 2
        
    if cumu:
        numbrows = 1
    elif full:
        numbrows = globdata.numbspec / 2
    else:
        numbrows = 2
        
    if globdata.exprtype == 'ferm':
        strgvarb = '$f$'
        strgunit = ' [1/cm$^2$/s/GeV]'
    if globdata.exprtype == 'sdss':
        strgvarb = '$C$'
        strgunit = ' [counts]'
    titl = strgvarb + strgunit

    for l in globdata.indxpopl:
        for i in globdata.indxenerfdfn:
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
                    if globdata.trueinfo:
                        if h < 3 or full:
                            indxpnts = where((globdata.binsspec[i, h] < globdata.truespec[l][0, i, :]) &                                           (globdata.truespec[l][0, i, :] < globdata.binsspec[i, h+1]))[0]
                        else:
                            indxpnts = where(globdata.binsspec[i, 3] < globdata.truespec[l][0, i, :])[0]

                        mar1 = axis.scatter(globdata.truelgal[l][indxpnts],                                           globdata.truebgal[l][indxpnts],                                           s=100, alpha=0.5, marker='x', lw=2, color='g')
                        
                    axis.set_xlabel(globdata.longlabl)
                    axis.set_ylabel(globdata.latilabl)
                    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
                    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
                    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
                    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
                    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
                    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
                    if not cumu:
                        if h < 3 or full:
                            axis.set_title(tdpy.util.mexp(globdata.binsspec[i, h]) + ' $<$ ' + strgvarb + ' $<$ ' + tdpy.util.mexp(globdata.binsspec[i, h+1]))
                        else:
                            axis.set_title(tdpy.util.mexp(globdata.binsspec[i, h]) + ' $<$ ' + strgvarb)
            figr.savefig(globdata.plotpath + 'pntsbind' + ptag + '%d%d' % (l, globdata.indxenerincl[i]) + '_' + globdata.rtag + '.png')
            plt.close(figr)
       
    
def plot_king(globdata):

    angl = rad2deg(globdata.angldisp)

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
                axis.plot(angl, retr_singking(angl, sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        axis.legend(loc=lloc)
        axis.set_yscale('log')
        axis.set_xlabel(r'$\theta$ ' + globdata.strganglunit)
        axis.set_xlabel(r'$\mathcal{K}(\theta)$')
        plt.figtext(0.7, 0.7, '$\mathcal{K}(\theta) = \frac{1}{2\pi\sigma^2}(1-\frac{1}{\gamma}(\frac{x^2}{2\sigma^2\gamma})^{-\gamma})$')
        
    figr.subplots_adjust()
    plt.savefig(globdata.plotpath + 'king.png')
    plt.close(figr)
    
    
def plot_psfn(globdata, thispsfn):
    
    if globdata.exprtype == 'sdss':
        globdata.angldisptemp = rad2deg(globdata.angldisp) * 3600.
    if globdata.exprtype == 'ferm':
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
                if globdata.trueinfo:
                    axis.plot(globdata.angldisptemp, globdata.truepsfn[i, :, m], label='Mock', color='g', ls='--')
                axis.set_yscale('log')
                if m == globdata.numbevtt - 1:
                    axis.set_xlabel(r'$\theta$ ' + globdata.strganglunit)
                if i == 0 and globdata.exprtype == 'ferm':
                    axis.set_ylabel(globdata.evttstrg[m])
                if m == 0:
                    axis.set_title(globdata.enerstrg[i])  
                if i == globdata.numbener - 1 and m == globdata.numbevtt - 1:
                    axis.legend(loc=2)
                indxsamp = globdata.indxsamppsfipara[i*globdata.numbformpara+m*globdata.numbener*globdata.numbformpara]
                if globdata.psfntype == 'singgaus':
                    strg = r'$\sigma = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp]) + globdata.strganglunit
                elif globdata.psfntype == 'singking':
                    strg = r'$\sigma = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp]) + globdata.strganglunit + '\n'
                    strg += r'$\gamma = %.3g$' % globdata.thissampvarb[indxsamp+1]
                elif globdata.psfntype == 'doubgaus':
                    strg = r'$f = %.3g$' % globdata.thissampvarb[indxsamp] + '\n'
                    if globdata.exprtype == 'sdss':
                        paratemp = rad2deg(globdata.thissampvarb[indxsamp+1]) * 3600.
                    if globdata.exprtype == 'ferm':
                        paratemp = rad2deg(globdata.thissampvarb[indxsamp+1])
                    strg += r'$\sigma = %.3g$ ' % paratemp + globdata.strganglunit + '\n'
                    if globdata.exprtype == 'sdss':
                        paratemp = rad2deg(globdata.thissampvarb[indxsamp+2]) * 3600.
                    if globdata.exprtype == 'ferm':
                        paratemp = rad2deg(globdata.thissampvarb[indxsamp+2])
                    strg += r'$\sigma = %.3g$ ' % paratemp + globdata.strganglunit
                elif globdata.psfntype == 'gausking':
                    strg = r'$f_G = %.3g$' % globdata.thissampvarb[indxsamp] + '\n'
                    strg += r'$\sigma_G = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp+1]) + globdata.strganglunit + '\n'
                    strg += r'$\sigma_K = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp+2]) + globdata.strganglunit + '\n'
                    strg += r'$\gamma = %.3g$' % globdata.thissampvarb[indxsamp+3]
                elif globdata.psfntype == 'doubking':
                    strg = r'$f_c = %.3g$' % globdata.thissampvarb[indxsamp] + '\n'
                    strg += r'$\sigma_c = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp+1]) + globdata.strganglunit + '\n'
                    strg += r'$\gamma_c = %.3g$' % globdata.thissampvarb[indxsamp+2] + '\n'
                    strg += r'$\sigma_t = %.3g$ ' % rad2deg(globdata.thissampvarb[indxsamp+3]) + globdata.strganglunit + '\n'
                    strg += r'$\gamma_t = %.3g$' % globdata.thissampvarb[indxsamp+4]
                axis.text(0.75, 0.75, strg, va='center', ha='center', transform=axis.transAxes, fontsize=18)
                
                if globdata.exprtype == 'ferm':
                    axis.set_ylim([1e0, 1e6])
                if globdata.exprtype == 'sdss':
                    axis.set_ylim([1e7, 1e11])

        plt.savefig(globdata.plotpath + 'psfnprof_' + globdata.rtag + '_%09d.png' % globdata.cntrswep)
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
    axis.set_xticklabels(['%.2g' % binsener[i] for i in globdata.indxener])
    axis.set_title('PSF FWHM')
    for i in globdata.indxener:
        for m in indxevtt:
            axis.text(meanener[i], indxevttincl[m]+0.5, r'$%.3g^\circ$' % rad2deg(tranfwhm[m, i]), ha='center', va='center', fontsize=14)

    figr.subplots_adjust(bottom=0.2)
    plt.savefig(globdata.plotpath + 'fwhmcnts_' + globdata.rtag + '_%09d.png' % globdata.cntrswep)
    plt.close(figr)
    
    
def plot_backcntsmean(globdata, backcntsmean):
    
    figr, axis = plt.subplots()

    traglobdata.numbbackcntsrofimean = transpose(backcntsmean)
    
    imag = axis.imshow(traglobdata.numbbackcntsrofimean, origin='lower', extent=[binsener[0], binsener[-1], 0, 4],                      cmap='BuPu', interpolation='none', alpha=0.4)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.set_xscale('log')
    axis.set_ylabel('PSF Class')
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_yticks([0.5, 1.5, 2.5, 3.5])
    axis.set_yticklabels(['0', '1', '2', '3'])
    axis.set_title('Mean FDM counts inumbside a PSF FWHM')
    for i in globdata.indxener:
        for m in indxevtt:
            axis.text(meanener[i], indxevttincl[m]+0.5, '%.3g' % traglobdata.numbbackcntsrofimean[m, i], ha='center', va='center')
            
    figr.subplots_adjust(bottom=0.2)
    plt.savefig(globdata.plotpath + 'backcnts_' + globdata.rtag + '_%09d.png' % globdata.cntrswep)
    plt.close(figr)
    
    
def plot_datacntshist(globdata):

    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(7 * globdata.numbener, 7 * globdata.numbevtt))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            datacntstemp = globdata.datacnts[i, :, m]
            
            maxmdatacntstemp = amax(datacntstemp)
            minmdatacntstemp = amin(datacntstemp[where(datacntstemp > 0.)])
            globdata.binscntstemp = linspace(minmdatacntstemp, maxmdatacntstemp, 20)
            meancntstemp = (globdata.binscntstemp[1:] + globdata.binscntstemp[0:-1]) / 2.
            diffcntstemp = globdata.binscntstemp[1:] - globdata.binscntstemp[0:-1]
            
            datacntshist = axis.hist(datacntstemp, globdata.binscntstemp, color='b')[0]

            init = [meancntstemp[argmax(datacntshist)], 1.]
            
            axis.set_xlim([amin(globdata.binscntstemp), amax(globdata.binscntstemp)])
            if m == globdata.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            #axis.set_xscale('log')
            if m == 0:
                axis.set_title(globdata.binsenerstrg[i])
            if i == 0 and globdata.exprtype == 'ferm':
                axis.set_ylabel(globdata.evttstrg[m])
        
    figr.subplots_adjust(wspace=0.3, hspace=0.2)
    plt.savefig(globdata.plotpath + 'datacntshist' + globdata.rtag + '.png')
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
        plt.savefig(os.environ["PCAT_DATA_PATH"] + '/png/talkintr.png', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_heal(globdata, heal, rofi=True, titl=''):
    
    if rofi:
        healtemp = copy(heal)
        heal = zeros(globdata.numbpixlheal)
        heal[globdata.indxpixlrofi] = healtemp

    cart = tdpy.util.retr_cart(heal, minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    imag = plt.imshow(cart, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.title(titl)
    plt.show()
    
    
def plot_3fgl_thrs(globdata):

    path = os.environ["PCAT_DATA_PATH"] + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = pf.getdata(path, 0)

    bgalfgl3 = linspace(-90., 90., 481)
    lgalfgl3 = linspace(-180., 180., 960)

    bgalexpo = linspace(-90., 90., 400)
    lgalexpo = linspace(-180., 180., 800)

    #fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)
    fluxthrs = griddata([lgalfgl3, bgalfgl3], fluxthrs, [globdata.lgalheal, globdata.bgalheal])

    cntsthrs = fluxthrs * globdata.expo

    jbgal = where(abs(bgalexpo) < 10.)[0]
    jlgal = where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)

    axis.set_title('3FGL Detection Flux Threshold [1/cm$^2$/s], 1.0 GeV - 10. GeV')
    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.savefig(globdata.plotpath + 'thrs_' + globdata.rtag + '.png')
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

        
def plot_histcnts(globdata, l, thiscnts):

    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(7 * globdata.numbener, 7 * globdata.numbevtt))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            if globdata.trueinfo:
                truehist = axis.hist(globdata.truecnts[l][i, :, m], globdata.binscnts[i, :], alpha=0.5, color='g', log=True, label=globdata.truelabl)
                if globdata.datatype == 'mock':
                    axis.hist(globdata.fgl3cnts[i, :, m], globdata.binscnts[i, :], alpha=0.5, color='red', log=True, label='3FGL')
            axis.hist(thiscnts[l][i, :, m], globdata.binscnts[i, :], alpha=0.5, color='b', log=True, label='Sample')
            if m == globdata.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            axis.set_xscale('log')
            if globdata.trueinfo:
                axis.set_ylim([0.1, 1e3])
            if m == 0:
                axis.set_title(globdata.binsenerstrg[i])
            if i == 0 and globdata.exprtype == 'ferm':
                axis.set_ylabel(globdata.evttstrg[m])
            if m == globdata.numbevtt / 2 and i == globdata.numbener / 2:
                axis.legend()
        
    figr.subplots_adjust(wspace=0.3)
    plt.savefig(globdata.plotpath + 'histcnts%d_' % l + globdata.rtag + '_%09d.png' % globdata.cntrswep)
    plt.close(figr)
    
def plot_datacnts(globdata, pener, pevtt, nextstat=False):

    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
    if pevtt == None:
        if pener == None:
            axis.set_title('')
        else:
            axis.set_title(globdata.binsenerstrg[pener])
    else:
        titl = globdata.binsenerstrg[pener]
        if globdata.exprtype == 'ferm':
            titl += ', ' + globdata.evttstrg[pevtt]
        axis.set_title(titl)

    # plot the model count map
    if pevtt == None:
        if pener == None:
            imag = axis.imshow(sum(globdata.datacntscart, axis=3), origin='lower', extent=globdata.exttrofi, interpolation='none')
        else:
            imag = axis.imshow(sum(globdata.datacntscart[:, :, pener, :], axis=2), origin='lower',                              interpolation='none', cmap='Reds', extent=globdata.exttrofi)
    else:
        imag = axis.imshow(globdata.datacntscart[:, :, pener, pevtt], origin='lower', interpolation='none',                          cmap='Reds', extent=globdata.exttrofi)
    
    if pevtt != None or pener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
    
    # superimpose catalogs
    for l in globdata.indxpopl:

        # true catalog
        if globdata.trueinfo:
            mrkrsize = retr_mrkrsize(globdata, globdata.truespec[l][0, pener, :], pener)
            lgal = copy(globdata.truelgal[l])
            bgal = copy(globdata.truebgal[l])
            if globdata.exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal[globdata.trueindxpntsmiss], bgal[globdata.trueindxpntsmiss], s=mrkrsize[globdata.trueindxpntsmiss],                        alpha=globdata.mrkralph, label=globdata.truelabl + ', missed', marker='x', linewidth=2, color='g')
            axis.scatter(lgal[globdata.trueindxpntsbias], bgal[globdata.trueindxpntsbias], s=mrkrsize[globdata.trueindxpntsbias],                        alpha=globdata.mrkralph, label=globdata.truelabl + ', biased', marker='o', linewidth=2, color='g')
            indxpnts = setdiff1d(arange(globdata.truenumbpnts, dtype=int), concatenate((globdata.trueindxpntsbias, globdata.trueindxpntsmiss)))
            axis.scatter(lgal[indxpnts], bgal[indxpnts], s=mrkrsize[indxpnts],                        alpha=globdata.mrkralph, label=globdata.truelabl + ', hit', marker='D', linewidth=2, color='g')
            for l in globdata.indxpopl:
                if globdata.indxtruepntstimevari[l].size > 0:
                    axis.scatter(lgal[globdata.indxtruepntstimevari[l]], bgal[globdata.indxtruepntstimevari[l]], s=100,                                label=globdata.truelabl + ', variable', marker='*', linewidth=2, color='y')

        # model catalog
        mrkrsize = retr_mrkrsize(globdata, globdata.thissampvarb[globdata.thisindxsampspec[l]][pener, :], pener)
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[l]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[l]]
        if globdata.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
            
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=globdata.mrkralph, label='Sample', marker='+', linewidth=2, color='b')

    if nextstat:
        
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(globdata, abs(modispec[pener, k]), pener)
            
            if globdata.exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=globdata.mrkralph,                        marker='o', linewidth=2, color=colr)
            
            if False:
                print 'modilgal[k]'
                print modilgal[k]
                print 'modibgal[k]'
                print modibgal[k]
                print
                
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == globdata.indxpropsplt or indxprop == globdata.indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        if pener == None:
            path = globdata.plotpath + 'datacntsAA_' + globdata.rtag + '_%09d.png' % globdata.cntrswep
        else:
            path = globdata.plotpath + 'datacnts%dA_' % pener + globdata.rtag + '_%09d.png' % globdata.cntrswep
    else:
        path = globdata.plotpath + 'datacnts%d%d_' % (pener, indxevttincl[pevtt]) + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_modlcnts(globdata, pener, pevtt):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
    if pevtt == None:
        axis.set_title(globdata.binsenerstrg[pener])
    else:
        titl = globdata.binsenerstrg[pener]
        if globdata.exprtype == 'ferm':
            titl += ', ' + globdata.evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        modlcntstemp = sum(globdata.thismodlcnts[pener, :, :], axis=1)
    else:
        modlcntstemp = globdata.thismodlcnts[pener, :, pevtt]
    if globdata.pixltype == 'heal':
        modlcntstemp = tdpy.util.retr_cart(modlcntstemp, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal, minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, \
            minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        modlcntstemp = modlcntstemp.reshape((globdata.numbsidecart, globdata.numbsidecart)).T
    modlcntstemp[where(modlcntstemp > globdata.datacntssatu[pener])] = globdata.datacntssatu[pener]
    
    imag = plt.imshow(modlcntstemp, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    
    # superimpose catalogs
    for l in globdata.indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(globdata, globdata.thissampvarb[globdata.thisindxsampspec[l]][pener, :], pener)
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[l]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[l]]
        if globdata.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=globdata.mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if globdata.trueinfo:
            mrkrsize = retr_mrkrsize(globdata, globdata.truespec[l][0, pener, :], pener)
            lgal = globdata.truelgal[l]
            bgal = globdata.truebgal[l]
            if globdata.exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=globdata.mrkralph, label=globdata.truelabl, marker='x', linewidth=2, color='g')

    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = globdata.plotpath + 'modlcnts%dA_' % pener + globdata.rtag + '_%09d.png' % globdata.cntrswep
    else:
        path = globdata.plotpath + 'modlcnts%d%d_' % (pener, indxevttincl[pevtt]) + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_resicnts(globdata, pener, pevtt, resicnts, nextstat=False):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
    if pevtt == None:
        axis.set_title(globdata.binsenerstrg[pener])
    else:
        titl = globdata.binsenerstrg[pener]
        if globdata.exprtype == 'ferm':
            titl += ', ' + globdata.evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        resicntstemp = sum(resicnts[pener, :, :], axis=1)
    else:
        resicntstemp = resicnts[pener, :, pevtt]
    if globdata.pixltype == 'heal':
        resicntstemp = tdpy.util.retr_cart(resicntstemp, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal, minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, \
            minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        resicntstemp = resicntstemp.reshape((globdata.numbsidecart, globdata.numbsidecart))
    resicntstemp[where(resicntstemp > globdata.resicntssatu[pener])] = globdata.resicntssatu[pener]
    resicntstemp[where(resicntstemp < -globdata.resicntssatu[pener])] = -globdata.resicntssatu[pener]
    
    imag = axis.imshow(resicntstemp, origin='lower', cmap='RdGy', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    # superimpose catalogs
    for l in globdata.indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(globdata, globdata.thissampvarb[globdata.thisindxsampspec[l]][pener, :], pener)
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[l]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[l]]
        if globdata.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=globdata.mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if globdata.trueinfo:
            mrkrsize = retr_mrkrsize(globdata, globdata.truespec[l][0, pener, :], pener)
            lgal = copy(globdata.truelgal[l])
            bgal = copy(globdata.truebgal[l])
            if globdata.exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=globdata.mrkralph, label=globdata.truelabl, marker='x', linewidth=2, color='g')

        
    if nextstat:
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(globdata, abs(modispec[pener, k]), pener)
            
            
            if globdata.exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=globdata.mrkralph,                        marker='o', linewidth=2, color=colr)

            
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == globdata.indxpropsplt or indxprop == globdata.indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

            
    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = globdata.plotpath + 'resicnts%dA_' % pener + globdata.rtag + '_%09d.png' % globdata.cntrswep
    else:
        path = globdata.plotpath + 'resicnts%d%d_' % (pener, indxevttincl[pevtt]) + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    
    plt.close(figr)
    
    
def plot_errrcnts(globdata, pener, pevtt, errrcntsrofi):

    if pevtt == None:
        errrcntstemp = sum(errrcntsrofi[pener, :, :], axis=1)
    else:
        errrcntstemp = errrcntsrofi[pener, :, pevtt]
    
    if globdata.pixltype == 'heal':
        errrcntstemp = tdpy.util.retr_cart(errrcntstemp, indxpixlrofi=globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal, \
            minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    else:
        errrcntstemp = errrcntstemp.reshape((globdata.numbsidecart, globdata.numbsidecart))
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
    if pevtt == None:
        axis.set_title(globdata.binsenerstrg[pener])
    else:
        titl = globdata.binsenerstrg[pener]
        if globdata.exprtype == 'ferm':
            titl += ', ' + globdata.evttstrg[pevtt]
        axis.set_title(titl)

    # plot the error count map
    imag = axis.imshow(errrcntstemp, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    
    if pevtt == None:
        path = globdata.plotpath + 'errrcnts%dA_' % pener + globdata.rtag + '_%09d.png' % globdata.cntrswep
    else:
        path = globdata.plotpath + 'errrcnts%d%d_' % (pener, indxevttincl[pevtt]) + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_catl(globdata, pener, pevtt, thiscnts):
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    axis.set_xlim([globdata.frambndrmarg, -globdata.frambndrmarg])
    axis.set_ylim([-globdata.frambndrmarg, globdata.frambndrmarg])
    if pevtt == None:
        axis.set_title(globdata.binsenerstrg[pener])
    else:
        titl = globdata.binsenerstrg[pener]
        if globdata.exprtype == 'ferm':
            titl += ', ' + globdata.evttstrg[pevtt]
        axis.set_title(titl)
        

    # superimpose catalogs
    for l in globdata.indxpopl:
        
        # model catalog
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[l]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[l]]
        if globdata.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=300, alpha=globdata.mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if globdata.trueinfo:
            lgal = copy(globdata.truelgal[l])
            bgal = copy(globdata.truebgal[l])
            if globdata.exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=300, alpha=globdata.mrkralph, label=globdata.truelabl, marker='x', linewidth=2, color='g')

    
    
    if globdata.trueinfo:
        for l in globdata.indxpopl:
            numbpnts = int(globdata.truenumbpnts[l])
            for a in range(numbpnts):
                if pevtt == None:
                    cnts = sum(globdata.truecnts[l][pener, a, :])
                    sigm = sqrt(sum(truesigm[l][pener, a, :]**2))
                else:
                    cnts = globdata.truecnts[l][pener, a, pevtt]
                    sigm = truesigm[l][pener, a, pevtt]
                axis.text(globdata.truelgal[l][a] + 0.7, globdata.truelgal[l][a] - 0.7,                         '%d/%.2f' % (cnts, sigm), color='g', fontsize=13)

    for l in globdata.indxpopl:
        numbpnts = int(globdata.thissampvarb[globdata.indxsampnumbpnts[l]])
        for a in range(numbpnts):
            if pevtt == None:
                cnts = sum(thiscnts[l][pener, a, :])
            else:
                cnts = thiscnts[l][pener, a, pevtt]
            axis.text(globdata.thissampvarb[globdata.thisindxsamplgal[l][a]] - 0.5, globdata.thissampvarb[globdata.thisindxsampbgal[l][a]] + 0.3,                     '%d' % cnts, color='b', fontsize=13)
    
    axis.axvline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(globdata.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-globdata.frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % mean(globdata.thisbackcntsmean[pener, :]), fontsize=18)
    else:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % globdata.thisbackcntsmean[pener, pevtt], fontsize=18)
        
    if pevtt == None:
        path = globdata.plotpath + 'catlcnts%dA_' % pener + globdata.rtag + '_%09d.png' % globdata.cntrswep
    else:
        path = globdata.plotpath + 'catlcnts%d%d_' % (pener, indxevttincl[pevtt]) + globdata.rtag + '_%09d.png' % globdata.cntrswep
    plt.savefig(path)
    plt.close(figr)

    
def plot_pntsdiff():

    tempindxpntsfull = []
    tempnumbpnts = 10
    tempindxpntsfull.append(range(tempnumbpnts))
    tempsamp = rand(maxmsampsize)
    tempsamp[globdata.indxsampnumbpnts] = array([tempnumbpnts])
    tempsamp[globdata.indxsampfdfnnorm] = cdfn_logt(array([tempnumbpnts]), minmfdfnnorm, factfdfnnorm)
    tempsamp[globdata.indxsampfdfnslop] = cdfn_atan(array([1.5]), minmfdfnslop, factfdfnslop)
    tempsamp[globdata.indxsamppsfipara] = 0.5
    for c in globdata.indxback:
        tempsamp[globdata.indxsampnormback[c, 0]] = cdfn_logt(array([1.]), minmnormback[c], factnormback[c])
    tempsampvarb, tempppixl, tempcnts,         temppntsflux, tempflux, tempcntsrofi = pars_samp(tempindxpntsfull, tempsamp)


    pntscnts = tempcnts[0, :, 0] * globdata.expo[0, :, 0] * globdata.apix * globdata.diffener[0]
    isotcnts = isotflux[0, :, 0] * globdata.expo[0, :, 0] * globdata.apix * globdata.diffener[0]
    
    totlcnts0 = isotcnts
    totlcnts1 = pntscnts / 1e6 + isotcnts
    
    nrept = 100
    totlcntsheal0 = zeros((nrept, globdata.numbpixlheal))
    totlcntsheal1 = zeros((nrept, globdata.numbpixlheal))
    for k in range(nrept):
        for n in range(globdata.indxpixlrofi.size):
            totlcntsheal0[k, globdata.indxpixlrofi[n]] = poisson(totlcnts0[n])
            totlcntsheal1[k, globdata.indxpixlrofi[n]] = poisson(totlcnts1[n])
          
    maxmcnts = max(amax(totlcntsheal0), amax(totlcntsheal1))
    globdata.binscnts = linspace(0., maxmcnts, maxmcnts + 2)
    diffcnts = globdata.binscnts[1:] - globdata.binscnts[:-1]
    meancnts = (globdata.binscnts[1:] + globdata.binscnts[:-1]) / 2.
    
    hist0 = empty((nrept, maxmcnts + 1))
    hist1 = empty((nrept, maxmcnts + 1))
    for k in range(nrept):
        hist0[k, :] = histogram(totlcntsheal0[k, globdata.indxpixlrofi], globdata.binscnts)[0].astype(float)
        hist0[k, :] *= 1. / sum(hist0[k, :]) / diffcnts
        hist1[k, :] = histogram(totlcntsheal1[k, globdata.indxpixlrofi], globdata.binscnts)[0].astype(float)
        hist1[k, :] *= 1. / sum(hist1[k, :]) / diffcnts
        
    
    totlcntscart0 = tdpy.util.retr_cart(totlcntsheal0[0, :], minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    totlcntscart1 = tdpy.util.retr_cart(totlcntsheal1[0, :], minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
    
    fig = plt.figure(figsize=(12, 12))
    axis = figr.add_subplot(221)
    imag = axis.imshow(totlcntscart0, origin='lower', cmap='Reds', extent=globdata.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    
    axis.set_xlabel(globdata.longlabl)
    axis.set_ylabel(globdata.latilabl)
    if globdata.exprtype == 'ferm':
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


    tdpy.util.plot_braz(ax, meancnts, hist0,  lcol='lightgreen',                         alpha=0.3, dcol='green', mcol='darkgreen', labl='Isotropic')
    tdpy.util.plot_braz(ax, meancnts, hist1, lcol='lightblue',                         alpha=0.3, dcol='blue', mcol='darkblue', labl='Isotropic + Unresolved PS')

    axis.set_title('Count PDF')
    axis.set_xlabel('$k$')
    axis.set_ylabel('d$P$/d$k$')
    axis.set_yscale('log')
    axis.legend()

    plt.savefig(globdata.plotpath + 'pntsdiff.png')
    plt.close(figr)
    

def pair_catl(globdata, thisindxpopl, modllgal, modlbgal, modlspec):

    indxmodl = zeros_like(globdata.truelgal[thisindxpopl], dtype=int) - 1
    dir2 = array([modllgal, modlbgal])
    for k in range(globdata.truelgal[thisindxpopl].size):
        dir1 = array([globdata.truelgal[thisindxpopl][k], globdata.truebgal[thisindxpopl][k]])
        dist = angdist(dir1, dir2, lonlat=True)
        indxdist = argmin(dist) 
        if dist[indxdist] < deg2rad(0.5):
            indxmodl[k] = indxdist

    indxtruepntsbias = where(amax(abs(modlspec[:, indxmodl] - globdata.truespec[thisindxpopl][0, :, :]) / globdata.truespec[thisindxpopl][0, :, :], axis=0) > 1.2)[0]
    indxtruepntsmiss = where(indxmodl == -1)[0]
    
    return indxmodl, indxtruepntsbias, indxtruepntsmiss
     

