# common imports
from __init__ import *

# internal functions
from main import init, initarry

def test_info(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.2
    dictargs['truenumbelempop0reg0'] = 100
    #dictargs['verbtype'] = 0
    
    dictargs['numbswep'] = 100000
    dictargs['numbsamp'] = 100
    
    listnamecnfgextn = ['fittlowr', 'fittnomi', 'fitthigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittlowr']['fittminmflux'] = 1e-9
    dictargsvari['fittnomi']['fittminmflux'] = 3e-9
    dictargsvari['fitthigr']['fittminmflux'] = 1e-8
    
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['fittminmflux']) for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['lliktotl', 'maxmlliktotl', 'levi', 'infoharm', 'bcom']
    listscalvarbcomp = ['self', 'self', 'self', 'self', 'self']
    listlablvarbcomp = ['$\ln P(D|M)$', '$\ln P(D|M_{max})$', '$\ln P(D)$', '$D_{KL}$', '$\eta_B$']
    listtypevarbcomp = ['errr', '', '', '', '']
    listpdfnvarbcomp = ['post', '', 'post', 'post', '']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        #forcprev=True, \
                        
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


def test_perf(strgcnfgextnexec=None):
   
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
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['sigclowr', 'sigcnomi', 'sigchigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['sigclowr']['truesigc'] = 0.5 / anglfact
    dictargsvari['sigcnomi']['truesigc'] = 2. / anglfact
    dictargsvari['sigchigr']['truesigc'] = 4. / anglfact

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
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        forcprev=True, \
                        
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

    
def test_psfn(strgcnfgextnexec=None):
    
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
    dictargs['optitype'] = 'none'
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    #oaxifree
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
    
    listnamevarbcomp = ['numbelempop0reg0']
    
    listnamevarbcomp = ['numbelempop0reg0']
    listscalvarbcomp = ['self']
    listlablvarbcomp = ['$N_{pts}$']
    listtypevarbcomp = ['errr']
    listpdfnvarbcomp = ['post']

    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        forcprev=True, \
                        
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

    
def test_anglassc(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'sdyn'
    dictargs['backtype'] = [['tgasback.fits']]
    dictargs['numbsidecart'] = 200
    dictargs['strgexpo'] = 1.
    dictargs['elemtype'] = ['clus']
    dictargs['psfnevaltype'] = 'kern'
    
    # temp
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['nomi', 'vhig', 'high', 'vlow', 'loww']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
   
    dictargsvari['vlow']['anglassc'] = 0.001
    
    dictargsvari['loww']['anglassc'] = 0.01
   
    dictargsvari['nomi']['anglassc'] = 0.1
    
    dictargsvari['high']['anglassc'] = 0.5
    
    dictargsvari['vhig']['anglassc'] = 1.
    
    lablxaxi = r'$\theta_{asc}$ [deg]'
    scalxaxi = 'self'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['anglassc']) for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['cmplpop0pop0reg0', 'fdispop0pop0reg0']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='anglassc', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                       )


def test_errr(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'sdyn'
    dictargs['backtype'] = [['tgasback.fits']]
    dictargs['numbsidecart'] = 200
    dictargs['strgexpo'] = 1.
    dictargs['elemtype'] = ['clus']
    dictargs['psfnevaltype'] = 'kern'
    
    # temp
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
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
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        forcprev=True, \
                        
                        namexaxi='specfraceval', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        
                       )

    #['logt', 100], \
    #numbangl=numbangl, \
    #binsangltype=binsangltype, \


def test_pars(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['optitype'] = 'none'
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['parsnega', 'parsnone', 'parsloww', 'parsnomi', 'parshigh']
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
    listnamevarbcomp = ['numbelempop0reg0', 'cmplpop0pop0reg0', 'fdispop0pop0reg0', 'fluxdistsloppop0']
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                       )


def test_spmr(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['optitype'] = 'none'
    dictargs['numbproc'] = 2
    dictargs['probtran'] = 1.
    dictargs['probspmr'] = 1.
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['radilowr', 'radihigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['radilowr']['radispmr'] = 0.5 / anglfact
    dictargsvari['radihigr']['radispmr'] = 2. / anglfact
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listnamevarbcomp = ['accpsplt', 'accpmerg', 'timeatcrcntpmaxm']
    listscalvarbcomp = ['self' for namevarbcomp in listnamevarbcomp]
    listlablvarbcomp = [r'$\alpha_{splt}$', r'$\alpha_{merg}$', r'$\tau_{ac}$']
    listtypevarbcomp = ['', '', '']
    listpdfnvarbcomp = ['', '', '']

    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['radispmr']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = initarry( \
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


globals().get(sys.argv[1])()
