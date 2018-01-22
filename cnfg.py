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
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 10
    dictargs['makeplot'] = False
    
    listnamecnfgextn = ['fittlowr', 'fittnomi', 'fitthigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittlowr']['minmflux'] = 1e-9
    dictargsvari['fittnomi']['minmflux'] = 3e-9
    dictargsvari['fitthigr']['minmflux'] = 1e-8
    
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listvarbxaxi = [dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn]
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['levi', 'infoharm', 'bcom']
    listscalvarbcomp = ['self', 'self', 'self']
    listtypevarbcomp = ['', '', '']
    listpdfnvarbcomp = ['post', 'post', '']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='minmflux', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listscalvarbcomp=listscalvarbcomp, \
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
    
    dictargsvari['sigclowr']['truesigc'] = 1e-9
    dictargsvari['sigcnomi']['truesigc'] = 3e-9
    dictargsvari['sigchigr']['truesigc'] = 1e-8

    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    scalxaxi = 'logt'
    listvarbxaxi = [dictargsvari[namecnfgextn]['truesigc'] for namecnfgextn in listnamecnfgextn]
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['truesigc'] for namecnfgextn in listnamecnfgextn] 
    
    #liststrgvarb = ['timereal', 'timeproctotl', 'timeproctotlswep', 'timeatcr', 'timeprocnorm', 'meanmemoresi', 'derimemoresi']
    #listlablvarb = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{MC}$', '$t_{CPU}^{\prime\prime}$ [s]', '$\bar{M}$', '$\partial_t\bar{M}$']
    listnamevarbcomp = ['timereal', 'timeproctotl', 'timeproctotlswep', 'timeatcr', 'timeprocnorm', 'meanmemoresi', 'derimemoresi', 'timerealtotl']
    listtypevarbcomp = ['' for namevarbcomp in listnamevarbcomp]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='truesigc', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listtypevarbcomp=listtypevarbcomp, \
                       )

    
def test_psfn(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.2
    dictargs['truenumbelempop0reg0'] = 100
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
    dictargsvari['psfnwfix']['initsigc'] = 0.5 / anglfact
   
    lablxaxi = 'PSF'
    scalxaxi = 'self'
    listvarbxaxi = [0., 1., 2.]
    listtickxaxi = ['Float', 'Fixed/Wrong', 'Fixed/True'] 
    
    listnamevarbcomp = ['timereal']
    listtypevarbcomp = ['' for namevarbcomp in listnamevarbcomp]

    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
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
   
    dictargsvari['vhig']['anglassc'] = 1.
    
    dictargsvari['high']['anglassc'] = 0.1
    
    dictargsvari['nomi']['anglassc'] = 0.5
    
    dictargsvari['loww']['anglassc'] = 0.01
   
    dictargsvari['vlow']['anglassc'] = 0.001
    
    lablxaxi = r'$theta_{asc}$ [deg]'
    scalxaxi = 'self'
    listvarbxaxi = [dictargsvari[namecnfgextn]['anglassc'] for namecnfgextn in listnamecnfgextn]
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['anglassc'] for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['postcmplpop0pop0reg0', 'postfdispop0pop0reg0']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='anglassc', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
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
    listvarbxaxi = [dictargsvari[namecnfgextn]['specfraceval'] for namecnfgextn in listnamecnfgextn]
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['specfraceval'] for namecnfgextn in listnamecnfgextn] 
    
    listnamevarbcomp = ['posterrrfracdimm', 'posterrrfrac']
    #listlablvarbcomp = [r'$\epsilon_m$ [%]', r'$\epsilon_d$ [%]']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='specfraceval', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
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
    listnamevarbcomp = ['postnumbelempop0reg0', 'postcmplpop0pop0reg0', 'postfdispop0pop0reg0', 'postfluxdistsloppop0']
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['priofactdoff'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['priofactdoff'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='priofactdoff', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
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
    listnamevarbcomp = ['accp', 'gmrb']
    listtypevarbcomp = ['', '']
    listpdfnvarbcomp = ['', '']
    
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['radispmr'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['radispmr'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        
                        namexaxi='radispmr', \
                        lablxaxi=lablxaxi, \
                        scalxaxi=scalxaxi, \
                        listvarbxaxi=listvarbxaxi, \
                        listtickxaxi=listtickxaxi, \
                        
                        listnamevarbcomp=listnamevarbcomp, \
                        listtypevarbcomp=listtypevarbcomp, \
                        listpdfnvarbcomp=listpdfnvarbcomp, \
                       )


globals().get(sys.argv[1])()
