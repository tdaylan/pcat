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
    
    listnamecnfgextn = ['fittlowr', 'fittnomi', 'fitthigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittlowr']['minmflux'] = 1e-9
    dictargsvari['fittnomi']['minmflux'] = 3e-9
    dictargsvari['fitthigr']['minmflux'] = 1e-8
    
    scalxaxi = 'logt'
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    listnamevarbcomp = ['levi', 'infoharm', 'bcom']
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        listnamevarbcomp=listnamevarbcomp, \
                        listvarbxaxi=listvarbxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        lablxaxi=lablxaxi, \
                        namexaxi='minmflux', \
                       )


def test_time():
   
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
    
    listnamecnfgextn = ['sigclowr', 'sigcnomi', 'sigchigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['sigclowr']['truesigc'] = 1e-9
    dictargsvari['sigcnomi']['truesigc'] = 3e-9
    dictargsvari['sigchigr']['truesigc'] = 1e-8
   
    #liststrgvarb = ['timereal', 'timeproctotl', 'timeproctotlswep', 'timeatcr', 'timeprocnorm', 'meanmemoresi', 'derimemoresi']
    #listlablvarb = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{MC}$', '$t_{CPU}^{\prime\prime}$ [s]', '$\bar{M}$', '$\partial_t\bar{M}$']

    scalxaxi = 'logt'
    lablxaxi = r'$f_{min}$ [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'
    listnamevarbcomp = ['levi', 'infoharm', 'bcom']
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['minmflux'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        listnamevarbcomp=listnamevarbcomp, \
                        listvarbxaxi=listvarbxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        lablxaxi=lablxaxi, \
                        namexaxi='minmflux', \
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
    dictargs['numbswep'] = 10000
    
    #oaxifree
    listnamecnfgextn = ['nomi', 'psfntfix', 'psfnwfix', 'psfngaus', 'elemfeww']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['psfntfix']['proppsfp'] = False
    
    dictargsvari['psfnwfix']['proppsfp'] = False
    dictargsvari['psfnwfix']['initsigc'] = 0.5 / anglfact
    
    dictargsvari['psfngaus']['psfntype'] = 'singgaus'
    
    dictargsvari['elemfeww']['fittmaxmnumbelem'] = 3
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                       )

    
def pcat_anglassc(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'sdyn'
    dictargs['backtype'] = [['tgasback.fits']]
    dictargs['numbsidecart'] = 200
    dictargs['strgexpo'] = 1.
    dictargs['elemtype'] = ['clus']
    dictargs['psfnevaltype'] = 'kern'
    # temp
    dictargs['numbswep'] = 10000
    
    listnamecnfgextn = ['nomi', 'vhig', 'high', 'vlow', 'loww']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
   
    dictargsvari['vhig']['anglassc'] = 1.
    
    dictargsvari['high']['anglassc'] = 0.1
    
    dictargsvari['loww']['anglassc'] = 0.01
   
    dictargsvari['vlow']['anglassc'] = 0.001
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        listnamecomp=['anglassc'], \
                       )


def test_errr():
      
    tupl = [ \
            [0.1, 'logt', 100], \
            [0.1, 'linr', 100], \
            [0.1, 'logt', 200, '0.1, log, 200'], \
             
           ]
    numbtupl = len(tupl)
    indxtupl = arange(numbtupl)
    stdverrrfracdimm = empty((2, numbtupl))
    stdverrrfrac = empty((2, numbtupl))
    strgtupl = empty(numbtupl, dtype=object)
    liststrgvarb = ['stdverrrfracdimm', 'stdverrrfrac']
    for k in range(numbtupl):
        
        specfraceval = tupl[k][0]
        binsangltype = tupl[k][1]
        numbangl = tupl[k][2]
        strgtupl[k] = '%3.1f, %s, %d' % (specfraceval, binsangltype, numbangl)
        gdat = init( \
                    diagmode=True, \
                    makeanim=True, \
                    proppsfp=False, \
                    numbangl=numbangl, \
                    binsangltype=binsangltype, \
                    indxenerincl=arange(2, 3), \
                    indxevttincl=arange(3, 4), \
                    back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                    strgexpo='fermexpo_cmp0_ngal.fits', \
                    psfntype='doubking', \
                    maxmnumbpnts=array([100]), \
                    maxmgangdata=deg2rad(10.), \
                    specfraceval=specfraceval, \
                    minmflux=1e-9, \
                    maxmflux=1e-5, \
                    truenumbpnts=array([50]), \
                   )
        for strg in liststrgvarb:
            varbdict[strg] = getattr(gdat, strg)

    size = 0.5
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_errr/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    liststrg = ['stdverrrfracdimm', 'stdverrrfrac']
    listlabl = [r'$\epsilon_m$ [%]', r'$\epsilon_d$ [%]']
    listvarb = [stdverrrfracdimm, stdverrrfrac]
    numbplot = len(liststrg)
    for k in range(numbplot):
        figr, axis = plt.subplots()
        print 'listvarb[k]'
        print listvarb[k]
        print 'listvarb[k]'
        print listvarb[k]

        axis.errorbar(indxtupl, listvarb[k], yerr=listvarb[k], marker='o')
        axis.set_ylabel(listlabl[k])
        axis.set_xticks(indxtupl + size)
        axis.set_xticklabels(strgtupl, rotation=45)
        plt.tight_layout()
        figr.savefig(path + '%s_%s.pdf' % (liststrg[k], strgtimestmp))
        plt.close(figr)


def test_tuto():

	init( \
	     maxmgang=deg2rad(20.), \
	     indxenerincl=arange(1, 4), \
	     indxevttincl=arange(2, 4), \
	     bgalcntr=pi/2., \
	     strgback=['fermisotflux.fits', 'fdfmflux_ngal.fits'], \
	     strgexpo='fermexpo_cmp0_ngal.fits', \
	     strgexprflux='fermflux_cmp0_ngal.fits', \
	    )


def test_pars(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['optitype'] = 'none'
    dictargs['numbswep'] = 100000
    dictargs['numbsamp'] = 100
    
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
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['priofactdoff'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['priofactdoff'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        listnamevarbcomp=listnamevarbcomp, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        listvarbxaxi=listvarbxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        lablxaxi=lablxaxi, \
                        namexaxi='priofactdoff', \
                       )


def test_atcr():
    
    #listnumbpntsmodi = array([1, 2, 3, 5, 10])
    listnumbpntsmodi = array([3, 5, 10])
    numbiter = listnumbpntsmodi.size
    timereal = empty(numbiter)
    timeatcr = empty(numbiter)
    timeproc = empty(numbiter)
    # temp
    #timeatcr = array([5670., 3420., 3042., 1023., 403.])
    #timeproc = array([103., 114., 105., 134., 140.])
    for k in range(numbiter):
        gdat = init( \
                    factthin=1, \
                    #verbtype=2, \
                    makeplot=False, \
                    numbpntsmodi=listnumbpntsmodi[k], \
                    indxenerincl=arange(2, 3), \
                    indxevttincl=arange(3, 4), \
                    back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                    strgexpo='fermexpo_cmp0_ngal.fits', \
                    psfntype='doubking', \
                    maxmnumbpnts=array([5]), \
                    maxmgangdata=deg2rad(10.), \
                    minmflux=3e-11, \
                    maxmflux=3e-7, \
                    truenumbpnts=array([2]), \
                   )
        timeatcr[k] = pcat.timeatcr
        timereal[k] = pcat.timerealtotl
        timeproc[k] = pcat.timeproctotl
        break

    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_atcr/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    figr, axis = plt.subplots()
    axis.plot(listnumbpntsmodi, timeatcr, ls='', lw=1, marker='o')
    axis.set_xlabel('$\delta N$')
    axis.set_ylabel(r'$\tau$')
    axis.set_xlim([0, amax(listnumbpntsmodi) + 1])
    plt.tight_layout()
    figr.savefig(path + 'timeatcr_%s.pdf' % strgtimestmp)
    plt.close(figr)

    figr, axis = plt.subplots()
    axis.plot(listnumbpntsmodi, numbswep / timeatcr / timeproc, ls='', lw=1, marker='o')
    axis.set_xlabel('$\delta N$')
    axis.set_ylabel(r'$\mathcal{P}$')
    axis.set_xlim([0, amax(listnumbpntsmodi) + 1])
    plt.tight_layout()
    figr.savefig(path + 'perfmetr_%s.pdf' % strgtimestmp)
    plt.close(figr)
        

def test_spmr(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['optitype'] = 'none'
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    dictargs['probtran'] = 1.
    dictargs['probspmr'] = 1.
    
    listnamecnfgextn = ['radilowr', 'radihigr']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['radilowr']['radispmr'] = 0.5 / anglfact
    dictargsvari['radihigr']['radispmr'] = 2. / anglfact
    
    scalxaxi = 'self'
    lablxaxi = r'$\alpha_p$'
    listnamevarbcomp = ['accp', 'gmrb']
    listtickxaxi = ['%.2g' % dictargsvari[namecnfgextn]['radispmr'] for namecnfgextn in listnamecnfgextn] 
    listvarbxaxi = [dictargsvari[namecnfgextn]['radispmr'] for namecnfgextn in listnamecnfgextn]
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamecnfgextn, \
                        listnamevarbcomp=listnamevarbcomp, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                        listvarbxaxi=listvarbxaxi, \
                        scalxaxi=scalxaxi, \
                        listtickxaxi=listtickxaxi, \
                        lablxaxi=lablxaxi, \
                        namexaxi='radispmr', \
                       )


globals().get(sys.argv[1])()
