# common imports
#from __init__ import *

# internal functions
#from main import init, initarry
import pcat
import sys
import numpy as np
#from astropy.coordinates import SkyCoord
#from pcat.util import retr_refrchaninit


def pcat_lionwrap():
    
    dictglob = pcat.main.init( \
                             pixltype='cart', \
                             listnameback='isot', \
                             backtype=[[1.]], \
                             datatype='inpt', \
                             exprtype='tess', \
                             strgexpo=1., \
                             lionmode=True, \
                             )
    


def pcat_fittminmflux_fittparsnone(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = ['fittvlow', 'fittloww', 'fittnomi', 'fitthigh', 'fittvhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittvlow']['fittminmflux'] = 3e-10
    dictargsvari['fittloww']['fittminmflux'] = 1e-9
    dictargsvari['fittnomi']['fittminmflux'] = 3e-9
    dictargsvari['fitthigh']['fittminmflux'] = 1e-8
    dictargsvari['fittvhig']['fittminmflux'] = 3e-8
    
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['fittminmflux']) for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['truelliktotl']
    listscalvarbcomp = ['self']
    listlablvarbcomp = ['$\ln P(D|M_{true})$']
    listtypevarbcomp = ['']
    listpdfnvarbcomp = ['']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                  
                                  strgpara='$PCAT_PATH/cnfg.py', \
                                  
                                  #forcprev=True, \
                                  #execpara=True, \
                                  
                                  namexaxi='fittminmflux', \
                                  lablxaxi=lablxaxi, \
                                  scalxaxi=scalxaxi, \
                                  listtickxaxi=listtickxaxi, \
                                  
                                  listnamevarbcomp=listnamevarbcomp, \
                                  listscalvarbcomp=listscalvarbcomp, \
                                  listlablvarbcomp=listlablvarbcomp, \
                                  listtypevarbcomp=listtypevarbcomp, \
                                  listpdfnvarbcomp=listpdfnvarbcomp, \
                                 )


def pcat_fittminmflux_fittparsnomi(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = ['fittvlow', 'fittloww', 'fittnomi', 'fitthigh', 'fittvhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittvlow']['fittminmflux'] = 3e-10
    dictargsvari['fittloww']['fittminmflux'] = 1e-9
    dictargsvari['fittnomi']['fittminmflux'] = 3e-9
    dictargsvari['fitthigh']['fittminmflux'] = 1e-8
    dictargsvari['fittvhig']['fittminmflux'] = 3e-8
    
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['fittminmflux']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \
                        
                        namexaxi='fittminmflux', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_fittpars_trueback(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truemaxmnumbelempop0reg0'] = 0
    dictargs['truenumbelempop0reg0'] = 0
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = [ \
                        'parsnega', 'parsnone', 'parsloww', 'parsnomi', 'parshigh', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['parsnega']['priofactdoff'] = -0.5
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    dictargsvari['parsnomi']['priofactdoff'] = 1.
    dictargsvari['parshigh']['priofactdoff'] = 1.5
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \

                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_fittpars_truepnts(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['minmflux'] = 1e-8
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = [ \
                        'parsnega', 'parsnone', 'parsloww', 'parsnomi', 'parshigh', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['parsnega']['priofactdoff'] = -0.5
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    dictargsvari['parsnomi']['priofactdoff'] = 1.
    dictargsvari['parshigh']['priofactdoff'] = 1.5
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \

                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_fittpars_truepntsfittminmfluxloww(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['fittminmflux'] = 1e-9
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = [ \
                        'parsnega', 'parsnone', 'parsloww', 'parsnomi', 'parshigh', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['parsnega']['priofactdoff'] = -0.5
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    dictargsvari['parsnomi']['priofactdoff'] = 1.
    dictargsvari['parshigh']['priofactdoff'] = 1.5
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \

                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_fittpars_truepntsfittminmfluxhigh(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['fittminmflux'] = 1e-8
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = [ \
                        'parsnega', 'parsnone', 'parsloww', 'parsnomi', 'parshigh', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['parsnega']['priofactdoff'] = -0.5
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    dictargsvari['parsnomi']['priofactdoff'] = 1.
    dictargsvari['parshigh']['priofactdoff'] = 1.5
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \

                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_truenumbelem(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.
    dictargs['fittmaxmnumbelempop0reg0'] = 1000
    dictargs['truemaxmnumbelempop0reg0'] = 1000
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = ['numbvlow', 'numbloww', 'numbnomi', 'numbhigh', 'numbvhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['numbvlow']['truenumbelempop0reg0'] = 10
    dictargsvari['numbloww']['truenumbelempop0reg0'] = 30
    dictargsvari['numbnomi']['truenumbelempop0reg0'] = 100
    dictargsvari['numbhigh']['truenumbelempop0reg0'] = 300
    dictargsvari['numbvhig']['truenumbelempop0reg0'] = 1000
    
    lablxaxi = r'$N_{pts}$'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['truenumbelempop0reg0']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \
                        
                        namexaxi='truenumbelempop0reg0', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_trueminmflux(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['indxenerincl'] = array([2])
    
    listnamecnfgextn = [ \
                        'truevlow', 'trueloww', 'truenomi', 'truehigh', 'truevhig', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['truevlow']['trueminmflux'] = 3e-10
    dictargsvari['trueloww']['trueminmflux'] = 1e-9
    dictargsvari['truenomi']['trueminmflux'] = 3e-9
    dictargsvari['truehigh']['trueminmflux'] = 1e-8
    dictargsvari['truevhig']['trueminmflux'] = 3e-8
    
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['trueminmflux']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        execpara=True, \

                        namexaxi='trueminmflux', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )


def pcat_perf(strgcnfgextnexec=None):
   
    ## pixelization type
    #[100,   'cart', 1e2, 1,          numbswepcomm, 1, 'Cartesian'], \
    ## sigc
    #[11.44, 'heal', 5e1, 1,          numbswepcomm, 1, '2X Max PS, 1/2X $f_{min}$'], \
    ## numbener
    #[11.44, 'heal', 1e2, 3,          numbswepcomm, 1, '3 energy bins'], \
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.2
    dictargs['truenumbelempop0reg0'] = 100
    
    dictargs['numbswep'] = 1000000
    dictargs['numbsamp'] = 1000
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['sigcloww', 'sigcnomi', 'sigchigh']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['sigcloww']['truesigc'] = 0.5 / anglfact
    dictargsvari['sigcnomi']['truesigc'] = 2. / anglfact
    dictargsvari['sigchigh']['truesigc'] = 4. / anglfact

    lablxaxi = r'$\sigma_c$ [$^\circ$]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(anglfact * dictargsvari[namecnfgextn]['truesigc']) for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['timereal', 'timeproctotl', 'timeproctotlswep', 'timeatcrcntpmaxm', 'timeprocnorm', 'meanmemoresi', 'timerealtotl']
    listnamevarbcomp += ['numbpixlprox%04d' % k for k in range(3)]
    listscalvarbcomp = ['self' for namevarbcomp in listnamevarbcomp]
    listlablvarbcomp = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{MC}$', r'$t_{CPU}^{\prime\prime}$ [s]', r'$\bar{M}$', r'$\partial_t\bar{M}$', '$t$ [s]']
    listlablvarbcomp += ['$N_{pxp,%d}$' % k for k in range(3)]
    listtypevarbcomp = ['' for namevarbcomp in listnamevarbcomp]
    listpdfnvarbcomp = ['' for namevarbcomp in listnamevarbcomp]
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        
                        namexaxi='truesigc', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listscalvarbcomp=listscalvarbcomp, \
                        listlablvarbcomp=listlablvarbcomp, \
                        listtypevarbcomp=listtypevarbcomp, \
                        listpdfnvarbcomp=listpdfnvarbcomp, \
                       )

    
def pcat_psfn(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 2e8
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['numbsidecart'] = 10
    dictargs['maxmgangdata'] = 0.492 / anglfact * 10 / 2.
    dictargs['minmflux'] = 1e-7
    dictargs['priofactdoff'] = 0.2
    dictargs['truenumbelempop0reg0'] = 10
    
    listnamecnfgextn = ['nomi', 'psfntfix', 'psfnwfix']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['psfntfix']['proppsfp'] = False
    
    dictargsvari['psfnwfix']['proppsfp'] = False
    dictargsvari['psfnwfix']['sigcen00evt0'] = 0.5 / anglfact
   
    lablxaxi = 'PSF'
    scalxaxi = 'self'
    listtickxaxi = ['Float', 'Fixed/Wrong', 'Fixed/True'] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        
                        namexaxi='truesigc', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                       )

    
def pcat_anglassc(strgcnfgextnexec=None):
    
    anglfact = 180. / pi
    
    dictargs = {}
    dictargs['truemaxmnumbelempop0reg0'] = 400
    dictargs['truenumbelempop0reg0'] = 400
    
    dictargs['listnameback'] = ['isot']
    dictargs['backtype'] = [[10.]]
    dictargs['truenumbpopl'] = 1
    dictargs['refrlegdpopl'] = ['PS']
    dictargs['trueelemtype'] = ['lghtpnts']
    dictargs['maxmgangdata'] = 10. / anglfact
    dictargs['truespatdisttype'] = ['self']
    dictargs['spectype'] = ['powr']
    dictargs['psfnevaltype'] = 'kern'
    dictargs['trueelemregitype'] = [True]
    dictargs['proppsfp'] = False
    
    dictargs['fittnumbpopl'] = 1
    dictargs['fittelemtype'] = ['lghtpnts']
    dictargs['fittspatdisttype'] = ['self']
    #dictargs['fittspectype'] = ['colr']
    dictargs['fittmaxmnumbelempop0reg0'] = 1000
    
    #dictargs['strgexprsbrt'] = 'sbrtfermrec8pntsigal0256.fits'
    #dictargs['spectype'] = ['colr']
    #dictargs['listnameback'] = ['isot', 'fdfm', 'dark']
    #dictargs['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'sbrtdarkpntssmthrec8.fits']]
    #dictargs['psfnevaltype'] = 'kern'
    
    dictargs['forccart'] = True
    dictargs['pixltype'] = 'cart'
    dictargs['numbsidecart'] = 100
    
    #dictargs['forccart'] = True
    #dictargs['pixltype'] = 'cart'
    #dictargs['numbsidecart'] = 100
    
    dictargs['numbswep'] = 100000
    dictargs['inittype'] = 'refr'
    #dictargs['makeplotfinlpost'] = False
    dictargs['numbsamp'] = 1000
    
    listnamecnfgextn = ['vlow', 'loww', 'nomi', 'high', 'vhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['vlow']['anglassc'] = 4. / anglfact
    dictargsvari['loww']['anglassc'] = 2. / anglfact
    dictargsvari['nomi']['anglassc'] = 1. / anglfact
    dictargsvari['high']['anglassc'] = 0.5 / anglfact
    dictargsvari['vhig']['anglassc'] = 0.2 / anglfact
    
    lablxaxi = r'$\theta_{asc}$ [deg]'
    scalxaxi = 'self'
    listtickxaxi = [tdpy.util.mexp(anglfact * dictargsvari[namecnfgextn]['anglassc']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                        
                                  namexaxi='anglassc', \
                                  lablxaxi=lablxaxi, \
                                  scalxaxi=scalxaxi, \
                                  listtickxaxi=listtickxaxi, \
                                 )


def pcat_errr(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'sdyn'
    dictargs['backtype'] = [['tgasback.fits']]
    dictargs['numbsidecart'] = 200
    dictargs['strgexpo'] = 1.
    dictargs['elemtype'] = ['clus']
    dictargs['psfnevaltype'] = 'kern'
    
    # temp
    dictargs['numbswep'] = 1000000
    dictargs['numbsamp'] = 1000
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['nomi', 'loww', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['nomi']['specfraceval'] = 0.1
    dictargsvari['loww']['specfraceval'] = 0.01
    dictargsvari['high']['specfraceval'] = 1.
    
    lablxaxi = r'$\Delta_f$ '
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['specfraceval']) for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['cntperrrreg0pop0maxm', 'cntperrrreg0pop0mean']
    #listlablvarbcomp = [r'$\epsilon_m$ [%]', r'$\epsilon_d$ [%]']
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        
                        namexaxi='specfraceval', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        
                       )

    #['logt', 100], \
    #numbangl=numbangl, \
    #binsangltype=binsangltype, \


def pcat_plot(strgcnfgextnexec=None):
  
    dictargs = {}
    
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    
    dictargs['mockonly'] = True
    dictargs['makeplotinit'] = True
    dictargs['makeplotintr'] = True
    
    listnamecnfgextn = ['nomi']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                       )


def pcat_elemspatevaltype(strgcnfgextnexec=None):
  
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    
    dictargs['numbswep'] = 1000000
    dictargs['numbsamp'] = 1000
    dictargs['mockonly'] = True
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['full', 'locllarg', 'loclsmal']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['full']['elemspatevaltype'] = ['full']
    
    dictargsvari['locllarg']['elemspatevaltype'] = ['loclhash']
    dictargsvari['locllarg']['maxmangleval'] = array([5., 7., 10.]) / anglfact
    
    dictargsvari['loclsmal']['elemspatevaltype'] = ['loclhash']
    dictargsvari['loclsmal']['maxmangleval'] = array([2., 3., 4.]) / anglfact
    
    scalxaxi = 'self'
    lablxaxi = ''
    listtickxaxi = ['Full', 'Local, LK', 'Local, SK']
    
    listnamevarbcomp = ['truelliktotl']
    listscalvarbcomp = ['self']
    listlablvarbcomp = ['$\ln P(D|M_{true})$']
    listtypevarbcomp = ['']
    listpdfnvarbcomp = ['']
   
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='elemspatevaltype', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listscalvarbcomp=listscalvarbcomp, \
                        listlablvarbcomp=listlablvarbcomp, \
                        listtypevarbcomp=listtypevarbcomp, \
                        listpdfnvarbcomp=listpdfnvarbcomp, \
                       )


def pcat_spmr(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['probtran'] = 1.
    dictargs['probspmr'] = 1.
    
    dictargs['numbswep'] = 1000000
    dictargs['numbsamp'] = 1000
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['radiloww', 'radihigh']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['radiloww']['radispmr'] = 0.5 / anglfact
    dictargsvari['radihigh']['radispmr'] = 2. / anglfact
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    
    listnamevarbcomp = ['accpsplt', 'accpmerg', 'timeatcrcntpmaxm']
    listscalvarbcomp = ['self' for namevarbcomp in listnamevarbcomp]
    listlablvarbcomp = [r'$\alpha_{splt}$', r'$\alpha_{merg}$', r'$\tau_{ac}$']
    listtypevarbcomp = ['', '', '']
    listpdfnvarbcomp = ['', '', '']

    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['radispmr']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        #forcprev=True, \
                        
                        namexaxi='radispmr', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listscalvarbcomp=listscalvarbcomp, \
                        listlablvarbcomp=listlablvarbcomp, \
                        listtypevarbcomp=listtypevarbcomp, \
                        listpdfnvarbcomp=listpdfnvarbcomp, \
                       )


def cnfg_fermi():
    '''
    The output folder will contain initial, indermediate, and output image and data files. To define this folder, tou should either
        1) set $PCAT_DATA_PATH environment variable, as in
            os.environ['PCAT_DATA_PATH'] = '/path/to/your/pcat/data/folder'
            inside your .bashrc or .zshrc, or
        2) input pathpcat as an argument.
    '''
    
    
    numbener = 5
    numbside = 256
    numbpixl = 12 * numbside**2
    numbevtt = 4
    fluxisot = 1e-6 * np.ones((numbener, numbpixl, numbevtt))
    
    # write the isotropic template to the PCAT input folder
    # Note that the default expected name is fermisotflux.fits
    # but could be different if strgback argument was provided 
    path = os.environ['PCAT_DATA_PATH'] + '/data/inpt/isottuto.fits'
    pf.writeto(path, fluxisot, clobber=True)
    
    # For this run, we will assume a Gaussian prior on the PSF parameters, centered on the best-fit Fermi-LAT PSF 
    # Therefore, we need to download the Fermi-LAT IRF files.
    cmnd = 'wget https://faun.rc.fas.harvard.edu/tansu/pcat/tuto/psf_P7REP_SOURCE_V15_back.fits $PCAT_DATA_PATH/data/inpt/psf_P7REP_SOURCE_V15_back.fits'
    #os.system(cmnd)
    
    # Now we can run PCAT
    pcat.init( \
                   forccart=True, \
                   pixltype='cart', \
                   diagmode=False, \
                   backtype=[1.], \
                   numbswep=2000000, \
                   strgexpo=1e11, \
                   probbrde=0.5, \
                   #lablback=['Isotropic'], \
                   )
                

def cnfg_GaussianMix():
    '''
    Gaussian Mixture
    '''
    
    pcat.init( \
             )


def cnfg_GaussianMix_unbinned():
    '''
    Unbinned Gaussian Mixture
    '''
    
    pcat.init( \
              boolbins=False, \
             )


def cnfg_lens_simu():

    anglfact = 3600. * 180. / np.pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21

    # half-size of the image
    numbside = 400
    maxmgangdata = numbside * 0.5 * sizepixl
    
    # name of the data file
    strgexprsbrt = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0en00 = 1.1e-7
        fittmeanbacpbac0en00 = 1.1e-7
        fittstdvbacpbac0en00 = fittmeanbacpbac0en00 * 1e-3
        fittscalbacpbac0en00 = 'gaus'
    else:
        initbacpbac0en00 = None
        fittmeanbacpbac0en00 = None
        fittstdvbacpbac0en00 = None
        fittscalbacpbac0en00 = None
    
    maxmnumbelem = np.array([0])

    dictfitt = dict()
    dictfitt['init'] = dict()
    dictfitt['fitt'] = dict()
    dictfitt['fitt']['maxm'] = dict()
    
    dictfitt['fitt']['maxm']['maxmnumbelem'] = maxmnumbelem
    dictfitt['init']['initbacpbac0en00'] = initbacpbac0en00
    dictfitt['init']['fittmeanbacpbac0en00'] = fittmeanbacpbac0en00
    dictfitt['init']['fittstdvbacpbac0en00'] = fittstdvbacpbac0en00
    dictfitt['init']['fittscalbacpbac0en00'] = fittscalbacpbac0en00
    pcat.main.init_image( \
                         typeexpr='HST_WFC3_IR', \
                         makeplotinit=False, \
                         typedata='simu', \
                         dictfitt=dictfitt, \

                         strgexpo=strgexpo, \
                         inittype='reco', \
                         namerecostat='pcat_lens_inpt', \
                         maxmgangdata=maxmgangdata, \
                         strgexprsbrt=strgexprsbrt, \
                        )
    

def intr_lens_evalcntpresi():
   
    pcat.main.init( \
                   typeexpr='HST_WFC3_IR', \
                   intrevalcntpresi=True, \
                   numbelempop0=1, \
                   maxmnumbelempop0=1, \
                   numbelempop1=1, \
                   maxmnumbelempop1=1, \
                   numbelempop2=1, \
                   maxmnumbelempop2=1, \
                  )
    

def intr_lens_evalcntpmodl():
   
    pcat.main.init( \
                   typeexpr='HST_WFC3_IR', \
                   intrevalcntpmodl=True, \
                   numbelempop0=1, \
                   maxmnumbelempop0=1, \
                   numbelempop1=1, \
                   maxmnumbelempop1=1, \
                   numbelempop2=1, \
                   maxmnumbelempop2=1, \
                  )
    

def pcat_lens_mock_truesgnl(strgcnfgextnexec=None):
   
    dictargs = {}
    
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['truenumbelempop0'] = 25
    dictargs['fittmaxmnumbelempop0'] = 25
    dictargs['typeelem'] = ['lens']
   
    listnamecnfgextn = ['nomi', 'truenone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['truenone']['truenumbelempop0'] = 0
    
    lablxaxi = ''
    scalxaxi = 'self'
    listtickxaxi = ['Nominal signal', 'No signal']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='truesgnl', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  tickxaxivari=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_trueminmdefs(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / np.pi
    
    dictpcatinpt = pcat.retr_dictpcatinpt()
    
    listlablcnfg = ['', '', '', '', '']
    listnamecnfgextn = ['truevlow', 'trueloww', 'nomi', 'truehigh', 'truevhig']
    dictpcatinptvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictpcatinptvari[namecnfgextn] = pcat.retr_dictpcatinpt()
    
    #dictpcatinpt['verbtype'] = 3
    #dictpcatinpt['probtran'] = 0.
    #dictpcatinpt['probspmr'] = 0.
    #dictpcatinpt['probspmr'] = 0.
    #dictpcatinpt['makeplot'] = False
    dictpcatinpt['typeexpr'] = 'HST_WFC3_IR'
    dictpcatinpt['boolbinsener'] = False
    
    #dictpcatinpt['sqzeexpl'] = True
    #dictpcatinpt['optitype'] = 'none'
    
    #dictpcatinpt['inittype'] = 'refr'
    dictpcatinpt['numbswep'] = 20000
    dictpcatinpt['numbswepplot'] = 5000
    dictpcatinpt['numbburn'] = 1000
    dictpcatinpt['numbsamp'] = 1000
    
    dictpcatinpt['plot_func'] = pcat.plot_lens

    #dictpcatinpt['typeelem'] = []
   
    #dictpcatinptvari['truevlow']['dicttrue']['minmdefs'] = 3e-4 / anglfact
    #dictpcatinptvari['trueloww']['dicttrue']['minmdefs'] = 1e-3 / anglfact
    #
    ##dictpcatinptvari['nomi']['dictfitt']['minmdefs'] = 3e-20 / anglfact
    ##dictpcatinptvari['nomi']['dicttrue']['minmdefs'] = 3e-20 / anglfact
    #
    #dictpcatinptvari['nomi']['dicttrue']['minmdefs'] = 3e-2 / anglfact
    #dictpcatinptvari['truehigh']['dicttrue']['minmdefs'] = 1e-2 / anglfact
    #dictpcatinptvari['truevhig']['dicttrue']['minmdefs'] = 3e-2 / anglfact

    lablxaxi = r'$\alpha_{min}$ [^{\prime\prime}]'
    scalxaxi = 'logt'
    #listtickxaxi = [tdpy.retr_lablmexp(anglfact * dictpcatinptvari[namecnfgextn]['dicttrue']['minmdefs']) for namecnfgextn in listnamecnfgextn] 
    
    dictpcatinpt['typeexpr'] = 'HST_WFC3_IR'
    
    #dictpcatinpt['dicttrue']['numbelempop0'] = 4
    #dictpcatinpt['dictfitt']['maxmpara'].numbelempop0 = 4
    dictpcatinpt['dictboth'] = dict()
    dictpcatinpt['dictboth']['typeelem'] = ['lens']

    dictglob = pcat.main.initarry( \
                                  dictpcatinptvari, \
                                  listlablcnfg, \
                                  
                                  dictpcatinpt=dictpcatinpt, \
                        
                                  boolexecpara=False, \
                                  
                                  namexaxivari='trueminmdefs', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  #tickxaxivari=listtickxaxi, \

                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_pars(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    
    dictargs = {}
    
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['truenumbelempop0'] = 25
    dictargs['typeelem'] = ['lens']
   
    listnamecnfgextn = ['none', 'loww', 'nomi', 'high', 'vhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['none']['priofactdoff'] = 0.
    dictargsvari['loww']['priofactdoff'] = 0.25
    dictargsvari['nomi']['boolsampprio'] = True
    dictargsvari['nomi']['priofactdoff'] = 0.5
    dictargsvari['high']['priofactdoff'] = 0.75
    dictargsvari['vhig']['priofactdoff'] = 1.
    
    #dictargs['numbswep'] = 1000
    #dictargs['numbsamp'] = 10
    #dictargs['verbtype'] = 2
    
    lablxaxi = r'$\alpha_{p}$'
    scalxaxi = 'self'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                 
                                  namexaxivari='priofactdoff', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  tickxaxivari=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_fittnumb(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / np.pi
    
    dictargs = {}
    
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['truenumbelempop0'] = 25
    dictargs['typeelem'] = ['lens']
   
    listnamecnfgextn = ['fittnone', 'fittsing', 'fittdoub', 'fittmany', 'nomi']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittnone']['fittminmnumbelempop0'] = 0
    dictargsvari['fittnone']['fittmaxmnumbelempop0'] = 0

    dictargsvari['fittsing']['fittminmnumbelempop0'] = 1
    dictargsvari['fittsing']['fittmaxmnumbelempop0'] = 1

    dictargsvari['fittdoub']['fittminmnumbelempop0'] = 2
    dictargsvari['fittdoub']['fittmaxmnumbelempop0'] = 2

    dictargsvari['fittmany']['fittminmnumbelempop0'] = 10
    dictargsvari['fittmany']['fittmaxmnumbelempop0'] = 10

    listtickxaxi = ['0', '1', '2', '10', 'Tran.']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='fittnumbelempop0', \
                                  listtickxaxi=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_trueback(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['truenumbelempop0'] = 25
    dictargs['typeelem'] = ['lens']
    
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    
    anglfact = 3600. * 180. / np.pi
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['loww']['truebacpbac0en00'] = 1e-7
    dictargsvari['nomi']['truebacpbac0en00'] = 2e-7
    dictargsvari['high']['truebacpbac0en00'] = 4e-7

    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['truebacpbac0en00']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='truebacp', \
                                  tickxaxivari=listtickxaxi, \

                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_lens_mock_sour(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / np.pi
    
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    
    dictargs['typeelem'] = ['lens', 'lghtgausbgrd']
    dictargs['numbelempop0'] = 20
    dictargs['maxmnumbelempop0'] = 100
    dictargs['numbelempop1'] = 10
    dictargs['maxmnumbelempop1'] = 100
    dictargs['spatdisttype'] = ['unif', 'dsrcexpo']
    
    numbelem = int(25. * 10.**0.9)

    listnamecnfgextn = ['nomi', 'datanone', 'subhsing', 'truelowr', 'parsnone', 'truevlow', 's2nrhigh', 's2nrvhig', 'amplhigh', 'bgrdunif']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['datanone']['killexpo'] = True
    
    dictargsvari['bgrdunif']['spatdisttype'] = ['unif', 'unif']
    
    dictargsvari['subhsing']['fittminmnumbelempop0'] = 1
    dictargsvari['subhsing']['fittmaxmnumbelempop0'] = 1

    dictargsvari['truevlow']['truenumbelempop0'] = 30
    
    dictargsvari['truelowr']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['truevlow']['truenumbelempop0'] = int(25. * 10.**0.9)
    dictargsvari['truevlow']['trueminmdefs'] = 3e-4 / anglfact
    dictargsvari['truevlow']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['parsnone']['priofactdoff'] = 0.
    
    dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    
    dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    dictargsvari['amplhigh']['minmdefs'] = 1e-1 / anglfact
        
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def test_lens_mock_sele():
    
    numbitermacr = 30
    numbiterelem = 10
    
    anglfact = 3600. * 180. / np.pi
    dictargs = {}
    dictargs['mockonly'] = True
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['minmdefs'] = 5e-4 / anglfact
    dictargs['numbelempop0'] = 1000
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['allwfixdtrue'] = False
    dictargs['maxmnumbelempop0'] = 1000
    
    listnamesele = ['pars', 'nrel']
    numbsele = len(listnamesele)
    listnamefeatsele = ['defs', 'mcut', 'rele']
    numbfeatsele = len(listnamefeatsele)

    dictargsvari = {}
    
    matrcutf = empty((numbitermacr, numbiterelem, numbfeatsele))
    
    liststrgvarboutp = []
    for strgvarbelem in listnamefeatsele:
        liststrgvarboutp += ['truehist' + strgvarbelem]
        for namesele in listnamesele:
            liststrgvarboutp += ['true' + strgvarbelem + namesele]
            liststrgvarboutp += ['truehist' + strgvarbelem + namesele]
    
    for k in range(numbitermacr):
        listgdat, dictglob = pcat.main.initarry( \
                                      dictargsvari, \
                                      dictargs, \
                                      listnamecnfgextn, \
                                      seedtypeelem='rand', \
                                      liststrgvarboutp=liststrgvarboutp, \
                                    
                                      strgcnfgextnexec=strgcnfgextnexec, \
                                     )
        
        gdat = listgdat[0]
        
        if k == 0:
            hist = zeros((numbiterelem, gdat.numbbinsplot))
            histsele = zeros((numbiterelem, gdat.numbbinsplot))
            histfitt = zeros((numbiterelem, gdat.numbbinsplot))
            factsele = zeros((numbiterelem, gdat.numbbinsplot))
            for namesele in listnamesele:
                pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
                os.system('mkdir -p %s' % pathimagsele)
            truefixp = empty((numbitermacr, gdat.truenumbfixp))
            corrfixpcutf = empty(gdat.truenumbfixp)
            pvalfixpcutf = empty(gdat.truenumbfixp)
        
        for namefixp in gdat.truenamefixp:
            truefixp[k, :] = gdat.truefixp

        for b, namefeat in enumerate(listnamefeatsele):
            lablvarbtotl = getattr(gdat, 'labl' + namefeat + 'totl')
            meanvarb = getattr(gdat, 'mean' + namefeat)
            deltvarb = getattr(gdat, 'delt' + namefeat)
            factvarbplot = getattr(gdat, 'fact' + namefeat + 'plot')
            limtvarb = [factvarbplot * getattr(gdat, 'minm' + namefeat), factvarbplot * getattr(gdat, 'maxm' + namefeat)]
            for m in range(numbiterelem):
                hist[m, :] = dictglob['truehist' + namefeat][m][0, :]
            for namesele in listnamesele:
                pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
                for m in range(numbiterelem):
                    histsele[m, :] = dictglob['truehist' + namefeat + namesele][m][0, :]
                    factsele[m, :] = histsele[m, :] / hist[m, :]
                    alph, loca, cutf = sp.stats.invgamma.fit(dictglob['true' + namefeat + namesele][m][0], floc=0., f0=1.9)
                    histfitt[m, :] = sum(histsele[m, :]) * sp.stats.invgamma.pdf(meanvarb, alph, loc=loca, scale=cutf) * deltvarb
                    matrcutf[k, m, b] = cutf
            
                meanfactsele = mean(factsele, 0)
                meanmatrcutf = mean(matrcutf, axis=1)

                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                for m in range(numbiterelem):
                    axis.loglog(meanvarb * factvarbplot, hist[m, :], color='g', alpha=0.1, ls='-.')
                    axis.loglog(meanvarb * factvarbplot, histsele[m, :], color='g', alpha=0.1, ls='--')
                    axis.loglog(meanvarb * factvarbplot, histfitt[m, :], color='g', alpha=0.1, ls='-')
                    axis.axvline(matrcutf[k, m, b] * factvarbplot, color='m', alpha=0.1)
                axis.axvline(np.power(prod(matrcutf[k, m, b]), array([1. / numbiterelem])) * factvarbplot, color='m')
    
                axis.loglog(meanvarb * factvarbplot, mean(hist, 0), color='g', ls='-.')
                axis.loglog(meanvarb * factvarbplot, mean(histsele, 0), color='g', ls='--')
                axis.loglog(meanvarb * factvarbplot, mean(histfitt, 0), color='g', ls='-')
    
                axis.set_xlabel(lablvarbtotl)
                axis.set_ylabel('$N$')
                axis.set_ylim([0.5, None])
                axis.set_xlim(limtvarb)
                plt.tight_layout()
                figr.savefig(pathimagsele + 'hist' + namefeat + '%04d.pdf' % k)
                plt.close(figr)
    
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.loglog(meanvarb * factvarbplot, meanfactsele, color='black')
                for m in range(numbiterelem):
                    axis.loglog(meanvarb * factvarbplot, factsele[m, :], alpha=0.1, color='g')
                    axis.axvline(matrcutf[k, m, b] * factvarbplot, color='m', alpha=0.1)
                axis.axvline(np.power(prod(matrcutf[k, :, b]), array([1. / numbiterelem])) * factvarbplot, color='m')
                axis.set_xlabel(lablvarbtotl)
                axis.set_ylabel('$f$')
                axis.set_xlim(limtvarb)
                plt.tight_layout()
                figr.savefig(pathimagsele + 'fact' + namefeat + '%04d.pdf' % k)
                plt.close(figr)
    
    for namesele in listnamesele:
        pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
        for b, namefeat in enumerate(listnamefeatsele):
            for a, namefixp in enumerate(gdat.truenamefixp):
                corrfixpcutf[a], pvalfixpcutf[a] = sp.stats.stats.pearsonr(truefixp[:, a], meanmatrcutf[:, b])
            indx = where(isfinite(corrfixpcutf) & (pvalfixpcutf < 0.1))[0]
            numb = indx.size
            figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
            for k in range(numb):
                size = 100. * (0.1 - pvalfixpcutf[k]) + 5.
                axis.plot(k + 0.5, corrfixpcutf[indx][k], ls='', marker='o', markersize=size, color='black')
            axis.set_xticks(arange(numb) + 0.5)
            axis.set_xticklabels(gdat.truelablfixp[indx])
            axis.set_xlim([0., numb])
            axis.set_ylabel(r'$\xi$')
            plt.tight_layout()
            figr.savefig(pathimagsele + 'corr' + namefeat + namesele + '.pdf')
            plt.close(figr)
    

def test_lens_mock_tmpr():
    
    anglfact = 3600. * 180. / np.pi
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['burntmpr'] = True
    dictargs['maxmnumbelempop0'] = 0
    dictargs['numbelempop0'] = 0

    listnamecnfgextn = ['nomi', 'refr', 'pert']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_many(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / np.pi
    
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['burntmpr'] = True
    for k in range(5):
        dictargs['maxmnumbelempop0reg%d' % k] = 0
        dictargs['numbelempop0reg%d' % k] = 0
    
    dictargs['inittype'] = 'pert'
    dictargs['typeelem'] = ['lens']
    dictargs['numbregi'] = 3
    dictargs['backtype'] = [1., 1., 1.]
    
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    
    listnamecnfgextn = ['nomi', 'regising', 'regimany']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['regising']['numbregi'] = 1 
    dictargsvari['regising']['backtype'] = [1.]
    
    dictargsvari['regimany']['numbregi'] = 5
    dictargsvari['regimany']['backtype'] = [1.] * 5
    
    for k in indxregi:
        dictglob = pcat.main.initarry( \
                                      dictargsvari, \
                                      dictargs, \
                                      listnamecnfgextn, \
                                      
                                      strgcnfgextnexec=strgcnfgextnexec, \
                                     )
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_spmr(strgcnfgextnexec=None):
  
    anglfact = 3600. * 180. / np.pi
    
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['inittype'] = 'refr'
    dictargs['numbelempop0'] = 1
    dictargs['truelgalpop00000'] = 1. / anglfact
    dictargs['truebgalpop00000'] = 0.5 / anglfact
    dictargs['truedefspop00000'] = 1e-2 / anglfact
    dictargs['probtran'] = 1.
    dictargs['typeelem'] = ['lens']
    dictargs['probspmr'] = 1.
    dictargs['indxenerincl'] = array([0])
    
    listnamecnfgextn = ['nomi', 'tranboth', 'parshigh', 'masshigh', 'massloww', 'trannone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['tranboth']['inittype'] = 'pert'
    dictargsvari['tranboth']['probtran'] = 0.4
    dictargsvari['tranboth']['probspmr'] = 0.3
    
    dictargsvari['parshigh']['priofactdoff'] = 1.
    
    dictargsvari['masshigh']['truedefspop00000'] = 3e-2 / anglfact
    
    dictargsvari['massloww']['truedefspop00000'] = 3e-3 / anglfact
    
    dictargsvari['trannone']['probtran'] = 0.
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def writ_data():

    # RA/DEC lists
    liststrgrade = []
    listrade = [[], []]
    
    pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_inpt/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    
    # read SLACS tables
    print('Reading SLACS tables...')
    pathslacpara = pathbase + 'data/slacpara.fits'
    pathslacfull = pathbase + 'data/slacfull.fits'
    hdun = pf.open(pathslacfull)
    numbhead = len(hdun)
    print('%s extensions found.' % numbhead)
    for k in range(numbhead):
        print('Extension %d' % k)
        head = hdun[k].header
        data = hdun[k].data
        
        if data == None:
            print('Data is None, skipping...')
            continue
        else:
            pass
    
        arry = array(stack((head.keys(), head.values()), 1))
        listtype = []
    
        for n in range(arry.shape[0]):
            if arry[n, 0].startswith('TTYPE'):
                listtype.append(arry[n, 1])
        
        if len(listtype) != len(data[0]):
            raise Exception('Number of types does not match the number of fields.')
        
    # find the RA/DEC of relevant SLACS systems 
    indxgold = where((data['Mph'] == 'E') & (data['Mul'] == 'S') & (data['Lens'] == 'A'))[0]
    numbslac = indxgold.size
    
    path = pathdata + 'slacdownlist.txt'
    fileobjt = open(path, 'w')
    for k in indxgold:
        
        # construct the delimited RA/DEC string
        strgrade = '%s %s %s %s %s %s' % (data['SDSS'][k][:2], data['SDSS'][k][2:4], data['SDSS'][k][4:9], data['SDSS'][k][9:12], data['SDSS'][k][12:14], data['SDSS'][k][14:])
        
        ## fill the RA/DEC lists
        liststrgrade.append(strgrade)
        listrade[0].append(data['_RA'][k])
        listrade[1].append(data['_DE'][k])
        
        ## write the RA/DEC list of relevant SLACS systems to disc
        strgline = strgrade + ' \n'
        fileobjt.write(strgline)
    
    for k in range(len(indxgold)):
        print('%20s %20s %20g %20g' % (data['SDSS'][indxgold[k]], data['Name'][indxgold][k], data['_RA'][indxgold][k], data['_DE'][indxgold][k]))
    
    # cutout properties
    numbside = 400
    numbsidehalf = numbside / 2
    
    # data path
    pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    
    ## RA/DEC string of the reference star
    #strgradestar = '00 29 12.65 -00 53 59.7'
    strgradestar = '00 29 06.79 -00 54 07.5'
    liststrgrade.append(strgradestar)
    coorstar = ap.coordinates.SkyCoord(strgradestar, unit=(ap.units.hourangle, ap.units.deg))
    listrade[0].append(coorstar.ra.degree)
    listrade[1].append(coorstar.dec.degree)
    
    numbrade = len(listrade[0])
    print('%d coordinates found.' % numbrade)
    
    # list of files to be read
    listnamefile = ['hst_10886_02_acs_wfc_f814w_drz.fits']
    numbfile = len(listnamefile)
    print('%d files found.' % numbfile)
    
    for k, namefile in enumerate(listnamefile):
        
        print('File number %d' % k)
            
        # read the data fields
        pathfile = pathdata + namefile
        listdata = tdpy.util.read_fits(pathfile, verbtype=0)
        
        # read the WCS header
        listhdun = ap.io.fits.open(pathfile)
        wcso = ap.wcs.WCS(listhdun[2].header)
        
        # RA/DEC string
        strgrade = liststrgrade[k]
    
        # iterate over the RA/DEC list    
        for n in range(numbrade):
            
            # RA/DEC
            strgrade = liststrgrade[n]
            indxyaxi, indxxaxi = wcso.wcs_world2pix(listrade[0][n], listrade[1][n], 0)
            # check if the coordinate is inside the image
            if not isfinite(indxyaxi) or not isfinite(indxxaxi) or indxxaxi - numbsidehalf < 0 or indxyaxi - numbsidehalf < 0 or \
                                                                        indxxaxi + numbsidehalf > listdata[1].shape[1] or indxyaxi + numbsidehalf > listdata[1].shape[0]:
                continue
                #raise Exception('')
    
            path = pathdatapcat + 'lens%s%s%s%s_%04d.fits' % (liststrgrade[n][3:5], liststrgrade[n][6:8], liststrgrade[n][16:18], liststrgrade[n][19:21], numbside)
           
            indxxaxi = int(indxxaxi)
            indxyaxi = int(indxyaxi)
           
            # cut out the image
            rate = listdata[1][indxxaxi-numbsidehalf:indxxaxi+numbsidehalf, indxyaxi-numbsidehalf:indxyaxi+numbsidehalf] # s^-1
    
            # gather different bands
            rate = rate[None, :, :, None]
            
            # find the number of photons per area per time per A per solid angle
            effa = 1. / listdata[4]['PHOTFLAM'][0] # erg^-1 cm^2 A
            timeobsv = listdata[4]['EXPTIME'][0] # s
            apix = (0.05 * np.pi / 3600. / 180.)**2 # sr^2
            expo = effa * timeobsv # erg^-1 cm^2 s A 
            sbrt = rate / effa / apix
            cntp = sbrt * expo * apix
            
            print('Writing to %s...' % path)
            pf.writeto(path, sbrt, clobber=True)
            
        

def pcat_lens_inpt(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / np.pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21
    
    # half-size of the image in pixels
    maxmgangdata = 100 * 0.5 * sizepixl
    maxmgangdatalarg = 400 * 0.5 * sizepixl

    # name of the data file
    strgexprsbrt = namedatasets + '_0100.fits'
    strgexprsbrtlarg = namedatasets + '_0400.fits'
    
    if namedatasets == 'lens29075550':
        initlgalsour = -0.1 / anglfact
        initbgalsour = 0.1 / anglfact
        initbacpbac0en00 = 1e-7
        fittmeanbacpbac0en00 = 1.115e-7
        fittstdvbacpbac0en00 = fittmeanbacpbac0en00 * 1e-3
        fittscalbacpbac0en00 = 'gaus'
    else:
        initlgalsour = None
        initbgalsour = None
        initbacpbac0en00 = None
        fittmeanbacpbac0en00 = None
        fittstdvbacpbac0en00 = None
        fittscalbacpbac0en00 = None

    listmask = [
                ['sqre', -0.3, 0.1, -0.1, 0.2] , \
                ['circ', -9, 8, 1] , \
               ]
    for k, mask in enumerate(listmask):
        for n, valu in enumerate(mask):
            if not isinstance(valu, str):
                listmask[k][n] = valu / anglfact
    
    dictargs = {}
    dictargs['typeelem'] = ['lens']
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['strgexpo'] = strgexpo
    dictargs['indxenerincl'] = array([0])
    dictargs['savestat'] = True
    
    # temp
    dictargs['sqzeprop'] = True
    dictargs['numbswep'] = 10000
    dictargs['numbswepplot'] = 1000
    
    dictargs['serstype'] = 'intp'
    #dictargs['verbtype'] = 2
    dictargs['numbsamp'] = 100
    dictargs['inittype'] = 'reco'
    
    #dictargs['inittype'] = 'rand'
    
    listnamecnfgextn = ['largrofi', 'largrofimask', 'nomi', 'mask', 'sour', 'dsrcexpo', 'sourmask', 'hostmult']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargs['initlgalhostisf0'] = -0.1 / anglfact
    dictargs['initbgalhostisf0'] = 0.
    dictargs['initlgalhostisf1'] = -0.1 / anglfact
    dictargs['initbgalhostisf1'] = 0.
    dictargs['initlgalsour'] = -0.2 / anglfact
    dictargs['initbgalsour'] = 0.2 / anglfact
    
    dictargsvari['largrofi']['maxmnumbelempop0'] = 0
    dictargsvari['largrofi']['maxmnumbelempop1'] = 0
    dictargsvari['largrofi']['maxmnumbelempop2'] = 0
    dictargsvari['largrofi']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofi']['maxmgangdata'] = maxmgangdatalarg
    #dictargsvari['largrofi']['thindata'] = True
    
    dictargsvari['largrofimask']['maxmnumbelempop0'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop1'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop2'] = 0
    dictargsvari['largrofimask']['listmask'] = listmask
    dictargsvari['largrofimask']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofimask']['maxmgangdata'] = maxmgangdatalarg
    
    dictargsvari['dsrcexpo']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['dsrcexpo']['maxmnumbelempop0'] = 2
    dictargsvari['dsrcexpo']['typeelem'] = ['lghtgausbgrd']
    #dictargsvari['dsrcexpo']['spatdisttype'] = ['unif', 'dsrcexpo']
    dictargsvari['dsrcexpo']['spatdisttype'] = ['dsrcexpo']
    dictargsvari['dsrcexpo']['dsrcdisttype'] = ['expo']
    dictargsvari['dsrcexpo']['dsrcdistsexppop0'] = 1e-2 / anglfact
    #dictargsvari['dsrcexpo']['numbburn'] = 0
    #dictargsvari['dsrcexpo']['forcsavestat'] = True
    
    dictargsvari['sour']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['sour']['maxmnumbelempop0'] = 1
    dictargsvari['sour']['typeelem'] = ['lghtgausbgrd']
    #dictargsvari['sour']['forcsavestat'] = True
    
    dictargsvari['hostmult']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['hostmult']['numbsersfgrd'] = array([2])
    #dictargsvari['hostmult']['proppsfp'] = False
    dictargsvari['hostmult']['maxmnumbelempop0'] = 0
    dictargsvari['hostmult']['maxmnumbelempop1'] = 0
    dictargsvari['hostmult']['typeelem'] = ['lens', 'lghtgausbgrd']
    #dictargsvari['hostmult']['forcsavestat'] = True
    #dictargsvari['hostmult']['listmask'] = listmask
    
    dictargsvari['sourmask']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['sourmask']['maxmnumbelempop0'] = 0
    dictargsvari['sourmask']['maxmnumbelempop1'] = 0
    dictargsvari['sourmask']['typeelem'] = ['lens', 'lghtgausbgrd']
    dictargsvari['sourmask']['listmask'] = listmask
    
    dictargsvari['nomi']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['nomi']['maxmgangdata'] = maxmgangdata
    dictargsvari['nomi']['maxmnumbelempop0'] = 0
    
    dictargsvari['mask']['listmask'] = listmask
    dictargsvari['mask']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['mask']['maxmgangdata'] = maxmgangdata
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def test_lens_mock_psfn():
   
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['truenumbelempop0'] = 25
    dictargs['typeelem'] = ['lens']
    dictargs['truesigcen00evt0'] = 4e-7
    
    anglfact = 3600. * 180. / np.pi
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['psfnfrwr']['truesigcen00evt0'] = 1e-7
    dictargsvari['psfnfrwr']['fittsigcen00evt0'] = 2e-7
    
    dictargsvari['psfnfxwr']['proppsfp'] = False
    dictargsvari['psfnfxwr']['initsigcen00evt0'] = 2e-7
    
    dictargsvari['psfnfrtr']['truesigcen00evt0'] = 4e-7

    dictargsvari['psfnfxtr']['proppsfp'] = False
    dictargsvari['psfnfxtr']['initsigcen00evt0'] = 4e-7

    listtickxaxi = ['Free, Wrong PSF', 'Fixed, Wrong PSF', 'Free, True PSF', 'Fixed, True PSF']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='psfnvari', \
                                  tickxaxivari=listtickxaxi, \

                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_papr(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['typeexpr'] = 'HST_WFC3_IR'
    dictargs['typeelem'] = ['lens']
    dictargs['seedtype'] = 4
    dictargs['priofactdoff'] = 0.
   
    # temp
    dictargs['limtydathistfeat'] = [0.5, 10.]
    
    dictargs['truemaxmnumbelempop0'] = 25
    dictargs['truenumbelempop0'] = 25
    #dictargs['maxmnumbelempop0'] = 0
    #dictargs['numbelempop0'] = 0
    
    anglfact = 3600. * 180. / np.pi
    
    numbelem = int(25. * 10.**0.9)

    listnamecnfgextn = ['fittlhig', 'fitthigh', 'fittvhig', 'truenone']
    #listnamecnfgextn = ['nomi', 'fittlhig', 'fitthigh', 'fittvhig', 'truenone', 'truenoneparsnomi', 'parsnomi', \
    #                                                                            'subhsing', 'trueloww', 's2nrhigh', 's2nrvhig', 'boolsampprio']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittlhig']['fittminmdefs'] = 0.005 / anglfact
    
    dictargsvari['fitthigh']['fittminmdefs'] = 0.01 / anglfact

    dictargsvari['fittvhig']['fittminmdefs'] = 0.02 / anglfact
    
    dictargsvari['truenone']['fittminmdefs'] = 0.01 / anglfact
    dictargsvari['truenone']['truenumbelempop0'] = 0
    
    #dictargsvari['parsnomi']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['parsnomi']['priofactdoff'] = 1.
    
    #dictargsvari['subhsing']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['subhsing']['probtran'] = 0.
    #dictargsvari['subhsing']['initnumbelempop0'] = 1
    #dictargsvari['subhsing']['fittminmnumbelempop0'] = 0
    #dictargsvari['subhsing']['fittmaxmnumbelempop0'] = 1
    #
    #dictargsvari['truenoneparsnomi']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['truenoneparsnomi']['truenumbelempop0'] = 0
    #dictargsvari['truenoneparsnomi']['priofactdoff'] = 1.
    #
    #dictargsvari['trueloww']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['trueloww']['truenumbelempop0'] = int(25. * 10.**0.9)
    #dictargsvari['trueloww']['trueminmdefs'] = 3e-4 / anglfact
    #
    #dictargsvari['s2nrhigh']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    #
    #dictargsvari['s2nrvhig']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    #dictargsvari['prio']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['prio']['boolsampprio'] = True
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )

    

globals().get(sys.argv[1])(*sys.argv[2:])
