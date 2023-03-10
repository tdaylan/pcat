# common imports
#from __init__ import *

# internal functions
#from main import init, initarry

from __init__ import *
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


globals().get(sys.argv[1])(*sys.argv[2:])
