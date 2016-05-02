
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
    globdata.modlpsfntype = hdun[0].header['modlpsfntype']
    if globdata.datatype == 'mock':
        globdata.mockpsfntype = hdun[0].header['mockpsfntype']
    globdata.exprtype = hdun[0].header['exprtype']
    globdata.pixltype = hdun[0].header['pixltype']

    if globdata.pixltype == 'heal':
        globdata.numbsideheal = hdun[0].header['numbsideheal']
    else:
        globdata.numbsidecart = hdun[0].header['numbsidecart']

    globdata.trueglobdata.info = hdun[0].header['trueglobdata.info']
    globdata.colrprio = hdun[0].header['colrprio']
    
    globdata.margsize = hdun[0].header['margsize']
    globdata.datetimestrg = hdun[0].header['datetimestrg']

    globdata.levi = hdun[0].header['levi']
    globdata.info = hdun[0].header['info']
    
    rtag = retr_rtag(None)

    if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
        plotfold = '/n/pan/www/tansu/png/pnts_tran/'
    else:
        plotfold = os.environ["PNTS_TRAN_DATA_PATH"] + '/png/'
    plotpath = plotfold + globdata.datetimestrg + '_' + rtag + '/'
    cmnd = 'mkdir -p ' + plotpath
    os.system(cmnd)


    lgalheal, bgalheal, globdata.numbsideheal, numbpixlheal, apix = tdpy_util.util.retr_heal(globdata.numbsideheal)
    globdata.indxpixlrofi = where((abs(lgalheal) < globdata.maxmgang) & (abs(bgalheal) < globdata.maxmgang))[0]

    globdata.indxenerincl = hdun['indxenerincl'].data
    globdata.indxevttincl = hdun['indxevttincl'].data
    
    globdata.minmmodlpsfipara, globdata.maxmmodlpsfipara,         globdata.factmodlpsfipara, globdata.strgmodlpsfipara,         globdata.scalmodlpsfipara, globdata.indxmodlpsfipara = retr_psfimodl(globdata, globdata.modlpsfntype, 'modl')

    if globdata.datatype == 'mock':
        minmmockpsfipara, maxmmockpsfipara, factmockpsfipara, strgmockpsfipara,             scalmockpsfipara, indxmockpsfipara = retr_psfimodl(globdata, globdata.mockpsfntype, 'mock')


        
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
        
    listspechist = hdun['spechist'].data
    pntsprob = hdun['pntsprob'].data
    
    listllik = hdun['llik'].data
    listlpri = hdun['lpri'].data
    
    gmrbstat = hdun['gmrbstat'].data
    listmodlcnts = hdun['modlcnts'].data
    
    # truth globdata.information
    if globdata.trueglobdata.info:
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
        globdata.truefdfnnorm = hdun['truefdfnnorm'].data
        globdata.truefdfnslop = hdun['truefdfnslop'].data
        globdata.truenormback = hdun['truenormback'].data
        globdata.truepsfipara = hdun['truepsfipara'].data
        
    # prior boundaries
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
    globdata.meanener = hdun['meanener'].data
    
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
    
    numbprop = len(globdata.probprop)

    globdata.maxmgangmarg = globdata.maxmgang + globdata.margsize
     
    globdata.exttrofi = [globdata.minmlgal, globdata.maxmlgal, globdata.minmbgal, globdata.maxmbgal]
    if globdata.exprtype == 'sdss':
        globdata.exttrofi *= 3600.
        globdata.frambndr = globdata.maxmgang * 3600.
        globdata.frambndrmarg = globdata.maxmgangmarg * 3600.
    else:
        globdata.frambndr = globdata.maxmgang
        globdata.frambndrmarg = globdata.maxmgangmarg
        
    globdata.strgfluxunit = retr_strgfluxunit(globdata.exprtype)

    # energy bin string
    globdata.enerstrg, globdata.globdata.binsenerstrg = retr_enerstrg(globdata.exprtype)
    
    if globdata.regitype == 'igal':
        globdata.longlabl = '$l$'
        globdata.latilabl = '$b$'
    else:
        globdata.longlabl = r'$\nu$'
        globdata.latilabl = r'$\mu$'
        
    if globdata.exprtype == 'ferm':
        longlabl += r' [$^\circ$]'
        latilabl += r' [$^\circ$]'
    if globdata.exprtype == 'sdss':
        longlabl += ' [arcsec]'
        latilabl += ' [arcsec]'
        
    if globdata.trueglobdata.info:
        if globdata.datatype == 'mock':
            globdata.truelabl = 'Mock data'
        else:
            globdata.truelabl = '3FGL'

        
    # Gelman-Rubin test
    if globdata.numbproc > 1:
        print 'Making the Gelman-Rubin TS plot...'
        tim0 = time.time()
        
        figr, axis = plt.subplots()
        axis.hist(gmrbstat, bins=linspace(1., amax(gmrbstat), 40))
        axis.set_title('Gelman-Rubin Convergence Test')
        axis.set_xlabel('PSRF')
        axis.set_ylabel('$N_{pix}$')
        figr.savefig(plotpath + 'gmrbdist_' + rtag + '.png')
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
    figr.savefig(plotpath + 'propeffi_' + rtag + '.png')
    figr.subplots_adjust(hspace=0.)
    plt.close(figr)
     

    indxsampsplt = where(globdata.listindxprop == indxpropsplt)[0]
    indxsampmerg = where(globdata.listindxprop == indxpropmerg)[0]
            
    listname = ['laccfrac', 'numbpair', 'combfact', 'jcbnfact']
    listvarb = [globdata.listlaccfrac, globdata.listnumbpair, globdata.listcombfact, globdata.listjcbnfact]
    for k in range(4):
        figr, axis = plt.subplots()
        axis.hist(listvarb[k][indxsampsplt])
        axis.hist(listvarb[k][indxsampmerg])
        axis.set_ylabel('$N_{samp}$')
        axis.set_xlabel(listname[k])
        figr.subplots_adjust(bottom=0.2)
        figr.savefig(plotpath + listname[k] + rtag + '.png')
        plt.close(figr)
    
    
    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Calculating proposal execution times...'
    tim0 = time.time()

    binstime = logspace(log10(amin(globdata.listchro[where(globdata.listchro > 0)] * 1e3)), log10(amax(globdata.listchro * 1e3)), 50)
    with sns.color_palette("nipy_spectral", numbprop):
        figr, axis = plt.subplots(figsize=(14, 12))
        
        axis.hist(globdata.listchro[where(globdata.listchro[:, 0] > 0)[0], 0] * 1e3, binstime, facecolor='none', log=True, lw=1, ls='--', edgecolor='black')
        for g in range(numbprop):
            indxglobdata.listchro = where((globdata.listindxprop == g) & (globdata.listchro[:, 0] > 0))[0]
            # temp
            axis.hist(globdata.listchro[indxglobdata.listchro, 0] * 1e3, binstime, edgecolor='none', log=True, alpha=0.3)#, label=propstrg[g])
        axis.set_xlabel('$t$ [ms]')
        axis.set_xscale('log')
        axis.set_xlim([amin(binstime), amax(binstime)])
        axis.set_ylim([0.5, None])
        axis.legend(loc=2)
        figr.savefig(plotpath + 'chroprop_' + rtag + '.png')
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
    figr.savefig(plotpath + 'chrototl_' + rtag + '.png')
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
                path = plotpath + 'postdist_%d_' % k + rtag
                if globdata.trueglobdata.info:
                    truepara = truesampvarb[indxsamppsfipara[ipsfipara]]
                else:
                    truepara = None
                tdpy_util.util.plot_mcmc(listpostdist, parastrgpsfipara[ipsfipara],                                     truepara=truepara, path=path, numbbins=numbbins, quan=True)
                
                
            
    # flux match with the true catalog
    if globdata.trueglobdata.info:
        for l in globdata.indxpopl:
            
            discspecmtch = zeros(globdata.truenumbpnts) + globdata.numbsamp
            listindxmodl = []
            for k in range(globdata.numbsamp):
                indxmodl, jtruepntsbias, jtruepntsmiss = pair_catl(globdata.truelgal[l],                                      globdata.truebgal[l],                                      globdata.truespec[l][0, :, :],                                      listlgal[l][k, :],                                      listbgal[l][k, :],                                      listspec[l][k, :, :])
                listindxmodl.append(indxmodl)
                discspecmtch[jtruepntsmiss] -= 1.
            discspecmtch /= globdata.numbsamp

            
            postspecmtch = zeros((3, globdata.numbener, globdata.truenumbpnts))
            for i in globdata.indxener:
                listspecmtch = zeros((globdata.numbsamp, globdata.truenumbpnts))
                for k in range(globdata.numbsamp):
                    indxpntstrue = where(listindxmodl[k] >= 0)[0]
                    listspecmtch[k, indxpntstrue] = listspec[l][k][i, listindxmodl[k][indxpntstrue]]
                postspecmtch[0, i, :] = percentile(listspecmtch, 16., axis=0)
                postspecmtch[1, i, :] = percentile(listspecmtch, 50., axis=0)
                postspecmtch[2, i, :] = percentile(listspecmtch, 84., axis=0)
            
            plot_scatspec(l, postspecmtch=postspecmtch)

            # store the comparison
            path = os.environ["PNTS_TRAN_DATA_PATH"] + '/pcatcomp_popl%d_' % l + rtag + '.fits'
            compbund = stack((globdata.truespec[l], postspecmtch))

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Making posterior distribution plots...'

    if globdata.pixltype == 'heal':
        
        reso = 60. / globdata.numbsideheal
        numbbinslgcr = (globdata.maxmlgal - globdata.minmlgal) / reso
        numbbinsbgcr = (globdata.maxmbgal - globdata.minmbgal) / reso

        pntsprobcart = zeros((numbbinslgcr, numbbinsbgcr, globdata.numbpopl, globdata.numbener, globdata.numbspec))
        for l in globdata.indxpopl:
            for i in globdata.indxener:
                for h in range(globdata.numbspec):
                    pntsprobcart[:, :, l, i, h] = tdpy_util.util.retr_cart(pntsprob[l, i, :, h], 
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
    plot_pntsprob(pntsprobcart, ptag='quad')
    plot_pntsprob(pntsprobcart, ptag='full', full=True)
    plot_pntsprob(pntsprobcart, ptag='cumu', cumu=True)
   
    # flux distribution
    for l in globdata.indxpopl:
        plot_histspec(l, listspechist=listspechist[:, l, :, :])

    # fraction of emission components
    postpntsfluxmean = retr_postvarb(globdata.listpntsfluxmean)
    postnormback = retr_postvarb(globdata.listnormback)
    plot_compfrac(postpntsfluxmean=postpntsfluxmean, postnormback=postnormback)


    # PSF parameters
    path = plotpath + 'psfipara_' + rtag
    if globdata.modlpsfntype == 'singgaus' or globdata.modlpsfntype == 'singking':
        globdata.listpsfipara[:, jpsfiparainit] = rad2deg(globdata.listpsfipara[:, jpsfiparainit])
        if globdata.trueglobdata.info:
            globdata.truepsfipara[jpsfiparainit] = rad2deg(globdata.truepsfipara[jpsfiparainit])
    elif globdata.modlpsfntype == 'doubgaus' or globdata.modlpsfntype == 'gausking':
        globdata.listpsfipara[:, jpsfiparainit+1] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+1])
        globdata.listpsfipara[:, jpsfiparainit+2] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+2])
        if globdata.trueglobdata.info:
            globdata.truepsfipara[jpsfiparainit+1] = rad2deg(globdata.truepsfipara[jpsfiparainit+1])
            globdata.truepsfipara[jpsfiparainit+2] = rad2deg(globdata.truepsfipara[jpsfiparainit+2])
    elif globdata.modlpsfntype == 'doubking':
        globdata.listpsfipara[:, jpsfiparainit+1] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+1])
        globdata.listpsfipara[:, jpsfiparainit+3] = rad2deg(globdata.listpsfipara[:, jpsfiparainit+3])
        if globdata.trueglobdata.info:
            globdata.truepsfipara[jpsfiparainit+1] = rad2deg(globdata.truepsfipara[jpsfiparainit+1])
            globdata.truepsfipara[jpsfiparainit+3] = rad2deg(globdata.truepsfipara[jpsfiparainit+3])
    if globdata.trueglobdata.info and psfntype == 'doubking':
        truepara = globdata.truepsfipara
    else:
        truepara = array([None] * globdata.numbpsfipara)
    tdpy_util.util.plot_mcmc(globdata.listpsfipara, strgpsfipara, truepara=truepara,                         nplot=globdata.numbformpara, path=path, numbbins=numbbins, quan=True, ntickbins=3)
    
    for k in range(globdata.numbpsfipara):
        path = plotpath + 'psfipara%d_' % k + rtag + '.png'
        tdpy_util.util.plot_trac(globdata.listpsfipara[:, k], strgpsfipara[k], path=path, quan=True)
    
    
    # log-likelihood
    path = plotpath + 'llik_' + rtag + '.png'
    tdpy_util.util.plot_trac(listllik.flatten(), '$P(D|y)$', path=path)

    # log-prior
    path = plotpath + 'lpri_' + rtag + '.png'
    tdpy_util.util.plot_trac(listlpri.flatten(), '$P(y)$', path=path)
    

    # number, expected number of PS and flux conditional prior power law index 
    for l in range(globdata.numbpopl):
        
        # number of point sources
        path = plotpath + 'numbpntsdist_popl%d_' % l + rtag + '.png'
        if globdata.trueglobdata.info and globdata.truenumbpnts != None:
            truepara = globdata.truenumbpnts[l]
        else:
            truepara = None
        tdpy_util.util.plot_trac(globdata.listnumbpnts[:, l], '$N$', truepara=truepara, path=path)

        # mean number of point sources
        path = plotpath + 'fdfnnorm_popl%d_' % l + rtag + '.png'
        if globdata.trueglobdata.info and globdata.truefdfnnorm != None:
            truepara = globdata.truefdfnnorm[l]
        else:
            truepara = None
        tdpy_util.util.plot_trac(globdata.listfdfnnorm[:, l], '$\mu$', truepara=truepara, path=path)

        # flux distribution power law index
        for i in globdata.indxenerfdfn:
            path = plotpath + 'fdfnslopdist_popl%d%d_' % (l, i) + rtag + '.png'
            if globdata.trueglobdata.info and globdata.truefdfnslop != None:
                truepara = globdata.truefdfnslop[l, i]
            else:
                truepara = None
            titl = globdata.binsenerstrg[i]
            labl =  r'$\alpha_{%d}$' % i
            tdpy_util.util.plot_trac(globdata.listfdfnslop[:, l, i], labl, truepara=truepara, path=path, titl=titl)
        
        
    # isotropic background normalization
    for i in globdata.indxener:
        path = plotpath + 'nisodist%d_' % i + rtag + '.png'
        if globdata.trueglobdata.info:
            if globdata.datatype == 'mock':
                truepara = globdata.truenormback[0, i]
            else:
                truepara = None
        titl = globdata.binsenerstrg[i]
        labl = r'$\mathcal{I}_{%d}$' % i
        tdpy_util.util.plot_trac(globdata.listnormback[:, 0, i], labl, truepara=truepara, path=path, titl=titl)
       
    if globdata.exprtype == 'ferm':
        # diffuse model normalization
        for i in globdata.indxener:
            path = plotpath + 'nfdmdist%d_' % i + rtag + '.png'
            if globdata.trueglobdata.info:
                if globdata.datatype == 'mock':
                    truepara = globdata.truenormback[1, i]
                else:
                    truepara = None
            else:
                truepara = None
            titl = globdata.binsenerstrg[i]
            labl = r'$\mathcal{D}_{%d}$' % i
            tdpy_util.util.plot_trac(globdata.listnormback[:, 1, i], labl, truepara=truepara, path=path, titl=titl)

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
    figr.savefig(plotpath + 'globdata.leviglobdata.info_' + rtag + '.png')
    plt.close(figr)

    #make_anim()

    tim1 = time.time()
    print 'Plots are produced in %.3g seconds.' % (tim1 - tim0)

def plot_compfrac(postpntsfluxmean=None, postnormback=None):
    
    if postpntsfluxmean != None:
        post = True
    else:
        post = False
        
        
    listlinestyl = ['-', '--', '-.', ':']
    listcolr = ['black', 'b', 'b', 'b']
    listlabl = ['Data', 'PS', 'Iso', 'FDM']

    figr, axis = plt.subplots(figsize=(7 * globdata.numbener, 7))
    
    listydat = empty((numbback + 2, globdata.numbener))
    listyerr = zeros((2, numbback + 2, globdata.numbener))
    
    listydat[0, :] = datafluxmean
    if post:
        listydat[1, :] = postpntsfluxmean[0, :]
        listyerr[:, 1, :] = retr_errrvarb(postpntsfluxmean)
        for c in indxback:
            listydat[c+2, :] = postnormback[0, c, :] * backfluxmean[c]
            listyerr[:, c+2, :] = retr_errrvarb(postnormback[:, c, :]) * backfluxmean[c]
    else:
        listydat[1, :] = mean(sum(thispntsflux * expo, 2) / sum(expo, 2), 1)
        for c in indxback:
            listydat[c+2, :] = thissampvarb[indxsampnormback[c, :]] * backfluxmean[c]

    
    xdat = globdata.meanener
    for k in range(numbback + 2):
        ydat = globdata.meanener**2 * listydat[k, :]
        yerr = globdata.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5,                     ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    # Fermi-LAT results
    if globdata.trueglobdata.info:
        if globdata.datatype == 'mock':
            
            pass
            
        else:
            
            if globdata.exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/fermspec' + name + '.csv'
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
        path = plotpath + 'compfracspec_' + rtag + '.png'
    else:
        path = plotpath + 'compfracspec_' + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    
    listlabl = ['PS', 'Iso', 'FDM']
    listexpl = [0.1, 0, 0]

    listsize = zeros(numbback + 1)
    for k in range(numbback + 1):
        if globdata.numbener == 1:
            listsize[k] = diffener * listydat[k+1, :]
        else:
            listsize[k] = trapz(listydat[k+1, :], globdata.meanener)
    listsize *= 100. / sum(listsize)
        
    figr, axis = plt.subplots()

    axis.pie(listsize, explode=listexpl, labels=listlabl, autopct='%1.1f%%')
    axis.axis('equal')

    if post:
        path = plotpath + 'compfrac_' + rtag + '.png'
    else:
        path = plotpath + 'compfrac_' + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
     

def plot_histsind(l, postsindhist=None):
    
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
        axis.hist(thissampvarb[thisindxsampsind[l]], globdata.binssind, alpha=0.5, color='b', log=True, label='Sample')
    if globdata.trueglobdata.info:
        axis.hist(globdata.truesind[l], globdata.binssind, alpha=0.5, color='g', log=True, label=truelabl)
        if globdata.datatype == 'mock':
            axis.hist(fgl3sind, globdata.binssind, alpha=0.1, color='red', log=True, label='3FGL')
    axis.set_yscale('log')
    axis.set_xlabel('$s$')
    axis.set_xlim([globdata.minmsind, globdata.maxmsind])
    axis.set_ylabel('$N$')
    if post:
        path = plotpath + 'histsind_popl%d' % l + rtag + '.png'
    else:
        path = plotpath + 'histsind_popl%d' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    

def plot_histspec(l, listspechist=None):
    
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
            axis.hist(thissampvarb[thisindxsampspec[l]][i, :], globdata.binsspec[i, :], alpha=0.5, color='b',                     log=True, label='Sample')
            
            if not globdata.colrprio or i == globdata.indxenerfdfn:
                fdfnslop = thissampvarb[indxsampfdfnslop[l, i]]  
                fluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[l]], fdfnslop, i)
                axis.plot(globdata.meanspec[i, :], fluxhistmodl, ls='--', alpha=0.5, color='b')
            
        if globdata.trueglobdata.info:
            truehist = axis.hist(globdata.truespec[l][0, i, :], globdata.binsspec[i, :],                                alpha=0.5, color='g', log=True, label=truelabl)
            if globdata.datatype == 'mock':
                axis.hist(fgl3spec[0, i, :], globdata.binsspec[i, :], color='red', alpha=0.1, log=True, label='3FGL')

        axis.set_yscale('log')
        axis.set_xlabel('$f$ ' + strgfluxunit)
        axis.set_xscale('log')
        axis.set_title(enerstrg[i])
        if globdata.trueglobdata.info:
            axis.set_ylim([0.1, 1e3])
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        if i == 0:
            axis.set_ylabel('$N$')
        if i == globdata.numbener / 2:
            axis.legend()
        
    figr.subplots_adjust(wspace=0.3, bottom=0.2)
    
    if post:
        path = plotpath + 'histspec%d_' % l + rtag + '.png'
    else:
        path = plotpath + 'histspec%d_' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
 

def plot_scatspec(l, postspecmtch=None, thisspecmtch=None):
    
    figr, axrw = plt.subplots(1, globdata.numbener, figsize=(7 * globdata.numbener, 6))
    if globdata.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        xdat = globdata.truespec[l][0, i, :]
        xerr = retr_errrvarb(globdata.truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
 
        if postspecmtch != None or thisspecmtch != None:
            labl = '$f_{samp}$ ' + strgfluxunit
            axis.plot(globdata.meanspec[i, :], globdata.meanspec[i, :], ls='--', color='black', alpha=0.2)
            if postspecmtch != None:
                yerr[0, :] = postspecmtch[0, i, :]
                ydat = postspecmtch[1, i, :]
                yerr[1, :] = postspecmtch[2, i, :]
                yerr = abs(yerr - ydat)
            else:
                ydat = thisspecmtch[i, :]

        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    
        if jtruepntstimevari[l].size > 0:
            axis.errorbar(xdat[jtruepntstimevari[l]], ydat[jtruepntstimevari[l]],                           ls='', yerr=yerr[:, jtruepntstimevari[l]],                           lw=1, marker='o', markersize=5, color='red')
            
    
        axis.set_xlabel('$f_{true}$ ' + strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_title(enerstrg[i])

        ylim = [globdata.minmspec[i], globdata.maxmspec[i]]
        
        axis.set_ylim(ylim)
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        axis.set_title(globdata.binsenerstrg[i])

    figr.subplots_adjust(wspace=0.4, bottom=0.2)

    if postspecmtch != None:
        path = plotpath + 'scatspec%d_' % l + rtag + '.png'
    elif thisspecmtch != None:
        path = plotpath + 'scatspec%d_' % l + rtag + '_%09d.png' % j

    plt.savefig(path)
    plt.close(figr)


def plot_scatpixl(l):
    
    figr, axgr = plt.subplots(globdata.numbevtt, globdata.numbener, figsize=(globdata.numbener * 7, globdata.numbevtt * 7))
    if globdata.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if globdata.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            slop, intc, coef, pval, stdv = sp.stats.linregress(datacnts[i, :, m], thismodlcnts[i, :, m])
            axis.scatter(datacnts[i, :, m], thismodlcnts[i, :, m], alpha=0.5)

            axislimt = [0., amax(datacnts[i, :, m]) * 1.5]
            axis.set_xlim(axislimt)
            axis.set_ylim(axislimt)
            
            axis.text(0.5, 0.9, '$m = %.3g, b = %.3g, r = %.3g$' % (slop, intc, coef),                     va='center', ha='center', transform=axis.transAxes, fontsize=16)
            if m == globdata.numbevtt - 1:
                axis.set_xlabel('Data Counts')
            if i == 0:
                labl = 'Model Counts'
                if globdata.exprtype == 'ferm':
                    labl += ', ' + evttstrg[m]
                axis.set_ylabel(labl)
            if m == 0:
                axis.set_title(enerstrg[i])


            
    figr.subplots_adjust(hspace=0.4, wspace=0.4, top=0.8)
    plt.savefig(plotpath + 'scatpixl%d_' % l + rtag + '_%09d.png' % j)
    plt.close(figr)
    
    
def plot_discspec(l, discspecmtch=None):
    
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
            labl = '$f_{samp}$ ' + strgfluxunit
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
        axis.set_xlabel('$f_{true}$ ' + strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_title(enerstrg[i])
        if globdata.colrprio:
            if i == globdata.indxenerfdfn[0]:
                ylim = [globdata.minmspec[i], globdata.maxmspec[i]]
        else:
            ylim = [globdata.minmspec[i], globdata.maxmspec[i]]      
        axis.set_ylim(ylim)
        axis.set_xlim([globdata.minmspec[i], globdata.maxmspec[i]])
        axis.set_title(globdata.binsenerstrg[i])
    figr.subplots_adjust(wspace=0.4, bottom=0.2)
    path = plotpath + 'discspec%d_' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)

     

