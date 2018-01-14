# common imports
from __init__ import *

# internal functions
from main import init, initarry

def test_info():
    
    listminmflux = logspace(-9., -8., 2)
    numbiter = listminmflux.size
    liststrgvarb = ['levi', 'info']
    varbdict = {}
    for k in range(numbiter):
        seedstat = get_state()
        gdat = init( \
                    seedstat=seedstat, \
                    numbswep=50000, \
                    factthin=500, \
                    indxenerincl=arange(1, 3), \
                    indxevttincl=arange(3, 4), \
                    minmflux=listminmflux[k], \
                    strgexpo=1e-20, \
                    maxmflux=1e-7, \
                    trueminmflux=1e-8, \
                    truenumbpnts=array([50]), \
                    back=['fermisotflux.fits'], \
                    #strgexpo='fermexpo_cmp0_ngal.fits', \
                   )
        for strg in liststrgvarb:
            if k == 0:
                varb = getattr(gdat, strg)
                if isinstance(varb, ndarray):
                    shap = [numbiter] + list(varb.shap)
                elif isinstance(varb, float):
                    shap = [numbiter, 1]
                varbdict[strg] = empty(shap)
            else:
                varbdict[strg][k, :] = getattr(gdat, strg)

    strgtimestmp = tdpy.util.retr_strgtimestmp()
    figr, axis = plt.subplots()
    axistwin = axis.twinx()
    axis.plot(listminmflux, varbdict['info'][:, 0], label='Relative entropy')
    axis.legend(bbox_to_anchor=[0.2, 1.08], loc=2)
    axistwin.plot(listminmflux, varbdict['levi'][:, 0], label='Log-evidence', color='g')
    axistwin.legend(bbox_to_anchor=[0.8, 1.08])
    axis.set_ylabel('$D_{KL}$ [nats]')
    axistwin.set_ylabel(r'$\log P(D)$ [nats]')
    axis.set_xlabel('$f_{min}$ [1/cm$^2$/s/GeV]')
    axis.set_xscale('log')
    plt.tight_layout()
    pathfold = os.environ["PCAT_DATA_PATH"] + '/imag/test_info/'
    os.system('mkdir -p ' + pathfold)
    figr.savefig(pathfold + 'infolevi_%s.pdf' % strgtimestmp)
    plt.close(figr)


def test_time():
   
    print 'Time-test suite for PCAT'

    numbswepcomm = 100

    tupl = [ \
            # reference
            [11.44, 'heal', 1e2, 1,          numbswepcomm, 1, 'Reference'], \
            
            # numbpixl
            [200,   'cart', 1e2, 1,          numbswepcomm, 1, '2X pixels'], \
            
            # pixelization type
            [100,   'cart', 1e2, 1,          numbswepcomm, 1, 'Cartesian'], \

            # minmflux
            [11.44, 'heal', 5e1, 1,          numbswepcomm, 1, '2X Max PS, 1/2X $f_{min}$'], \
            
            # numbener
            [11.44, 'heal', 1e2, 3,          numbswepcomm, 1, '3 energy bins'], \
            
            # numbproc
            [100,   'heal', 1e2, 1, int(numbswepcomm / 2), 2, '2 Processes'], \

           ]
    
    truenumbpnts = array([100])

    numbiter = len(tupl)
    indxiter = np.arange(numbiter)
    
    liststrgvarb = ['timereal', 'timeproctotl', 'timeproctotlswep', 'timeatcr', 'timeprocnorm', 'meanmemoresi', 'derimemoresi']
    listlablvarb = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{MC}$', '$t_{CPU}^{\prime\prime}$ [s]', '$\bar{M}$', '$\partial_t\bar{M}$']
    
    listticklabl = []
    for k in indxiter:
        temp, pixltype, minmflux, numbener, numbswep, numbproc, ticklabl = tupl[k]
        listticklabl.append(ticklabl)
        if pixltype == 'heal':
            maxmgangdata = deg2rad(temp)
            numbsidecart = None
        else:
            maxmgangdata = None
            numbsidecart = temp

        if numbener == 1:
            indxenerincl = arange(2, 3)
        if numbener == 3:
            indxenerincl = arange(1, 4)
     
        print 'tupl[k]'
        print tupl[k]

        gdat = init( \
                    numbswep=numbswep, \
                    numbproc=numbproc, \
                    back=[1.], \
                    makeplot=False, \
                    verbtype=0, \
                    strgexpo=1., \
                    indxenerincl=indxenerincl, \
                    pixltype=pixltype, \
                    indxevttincl=arange(3, 4), \
                    psfntype='doubking', \
                    maxmgangdata=maxmgangdata, \
                    minmflux=minmflux, \
                    maxmflux=1e6, \
                    numbsidecart=numbsidecart, \
                    truenumbpnts=truenumbpnts, \
                   )
        
        print
        
        if k == 0:
            varbdict = dict() 
            for strg in liststrgvarb:
                varb = getattr(gdat, strg)
                if isinstance(varb, float):
                    shap = [1]
                else:
                    shap = list(getattr(gdat, strg).shape)
                varbdict[strg] = empty([numbiter] + shap)
        
        for strg in liststrgvarb:
            print 'varbdict'
            print varbdict
            varbdict[strg][k, :] = getattr(gdat, strg)

    size = 0.5
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_time/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()

    numbplot = len(liststrgvarb)
    for k in range(numbplot):
        figr, axis = plt.subplots()
        axis.bar(indxiter, dictvarb[k, :], 2 * size)
        axis.set_ylabel(listlablvarb[k])
        axis.set_xticks(indxtupl + size)
        axis.set_xticklabels(listtticklabl, rotation=45)
        plt.tight_layout()
        figr.savefig(path + '%s_%04d_%s.pdf' % (liststrgvarb[k], log10(numbswep), strgtimestmp))
        plt.close(figr)

    
def test_psfn(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['priofactdoff'] = 0.2
    dictargs['truemaxmnumbelempop0reg0'] = 400
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
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )

    
def test_nomi():
      
    init( \
         proppsfp=False, \
         indxenerincl=arange(1, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([600]), \
         maxmgangdata=deg2rad(20.), \
         minmflux=5e-11, \
         maxmflux=1e-7, \
         truenumbpnts=array([300]), \
        )


def pcat_anglassc(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'sdyn'
    dictargs['backtype'] = [['tgasback.fits']]
    
    dictargs['numbsidecart'] = 200
    
    #dictargs['minmdatacnts'] = 0.
    dictargs['strgexpo'] = 1.
    dictargs['elemtype'] = ['clus']
    dictargs['psfnevaltype'] = 'kern'
    
    # temp
    dictargs['numbswep'] = 10000
    
    listnamecnfgextn = ['nomi', 'vhig', 'high', 'vlow', 'loww']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
   
    # to test splits and merges
    dictargsvari['vhig']['anglassc'] = 1.
    dictargsvari['high']['anglassc'] = 0.1
    dictargsvari['loww']['anglassc'] = 0.01
    dictargsvari['vlow']['anglassc'] = 0.001
    
    dictglob = initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                  listnamecomp=['anglassc'], \
                                 )


def plot_main_lens():
        init( \
             numbswep=1000, \
             makeplotintr=True, \
             pntstype='lens', \
             exprtype='hubb', \
             truenumbpnts=array([40]), \
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


def test_uppr():
      
    init( \
         #factthin=100, \
         proppsfp=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([800]), \
         maxmgangdata=deg2rad(10.), \
         minmflux=1e-9, \
         maxmflux=1e-5, \
         truenumbpnts=array([400]), \
        )


def test_spatprio():
      
    init( \
         indxenerincl=arange(1, 4), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         spatdisttype=['gaus'], \
        )


def test_pars(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['truemaxmnumbelempop0reg0'] = 400
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['optitype'] = 'none'
    dictargs['makeplot'] = False
    dictargs['numbswep'] = 100
    dictargs['numbsamp'] = 10
    
    listnamecnfgextn = ['parsnone', 'parsloww']
    #listnamecnfgextn = ['parsnone', 'parsloww', 'parsnomi', 'parshigh']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    #dictargsvari['parsnomi']['priofactdoff'] = 1.
    #dictargsvari['parshigh']['priofactdoff'] = 1.5
    
    listnamevarbcomp = ['numbelemreg0pop0']
    
    dictglob = initarry( \
                        dictargsvari, \
                        dictargs, \
                        listnamevarbcomp=listnamevarbcomp, \
                        strgcnfgextnexec=strgcnfgextnexec, \
                       )


def test_lowr():
      
    init( \
         numbswep=10000, \
         verbtype=1, \
         proppsfp=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexprflux='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         maxmnumbpnts=array([3]), \
         maxmgangdata=deg2rad(10.), \
         minmflux=1e-24, \
         maxmflux=1e-20, \
         truenumbpnts=array([3]), \
        )


def test_post():
     
    indxenerincl = arange(2, 4)
    indxevttincl = arange(3, 4)
    #tupl = [[0.3, 0.001, 0.001, 0.01, 0.01, 3e-8, 1e-7], \
    #        [0.3, 0.01, 0.01, 0.1, 0.1, 3e-9, 1e-8]]
    tupl = [[0.3, 0.001, 0.001, 0.01, 0.01, 3e-10, 1e-9], \
            [0.3, 0.01, 0.01, 0.1, 0.1, 3e-11, 1e-10]]
    numbiter = len(tupl)
    for k in range(numbiter):
        stdvback, stdvlbhlminm, stdvlbhlmaxm, stdvflux, stdvsind, minmflux, maxmflux = tupl[k]
        init( \
             numbburn=0, \
    		 factthin=1, \
             indxenerincl=indxenerincl, \
             indxevttincl=indxevttincl, \
             proptran=False, \
             prophypr=False, \
             back=['fermisotflux.fits'], \
             lablback=[r'$\mathcal{I}$'], \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             stdvback=stdvback, \
             stdvlbhlminm=stdvlbhlminm, \
             stdvlbhlmaxm=stdvlbhlmaxm, \
             stdvflux=stdvflux, \
             stdvsind=stdvsind, \
             psfntype='doubking', \
             maxmnumbpnts=array([3]), \
             maxmgangdata=deg2rad(1.), \
             minmflux=minmflux, \
             maxmflux=maxmflux, \
             truenumbpnts=array([3]), \
             truefluxdistslop=array([1.9]), \
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
        

def test_spmr():
    
    listminmflux = array([5e-21, 5e-11])
    listmaxmflux = array([1e-17, 1e-7])
    numbiter = listminmflux.size
    for k in range(numbiter):
        init( \
             #verbtype=2, \
             factthin=1, \
             indxenerincl=arange(2, 3), \
             indxevttincl=arange(3, 4), \
             back=['fermisotflux.fits'], \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             prophypr=False, \
             proppsfp=False, \
             propbacp=False, \
             propcomp=False, \
             probbrde=0., \
             psfntype='doubking', \
             maxmgangdata=deg2rad(1.), \
             minmflux=listminmflux[k], \
             maxmflux=listmaxmflux[k], \
             maxmnumbpnts=array([3]), \
             truenumbpnts=array([2]), \
             #maxmnumbpnts=array([3]), \
             #truenumbpnts=array([2]), \
            )
        

def test_cova():
     
    init( \
         numbswep=100, \
         verbtype=2, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         exprtype='ferm', \
         minmflux=3e-8, \
         maxmflux=3e-7, \
         trueminmflux=3e-11, \
         maxmnumbpnts=array([3]), \
         truenumbpnts=array([2]), \
        )


def test_leak():
     
    init( \
         numbswep=10000, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         exprtype='ferm', \
         minmflux=3e-8, \
         maxmflux=3e-7, \
         trueminmflux=3e-11, \
         truenumbpnts=array([100]), \
        )


def test_popl():
     
    init( \
         factthin=2, \
         #verbtype=2, \
         indxenerincl=arange(2, 4), \
         indxevttincl=arange(3, 4), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         maxmnumbpnts=array([3, 3]), \
         maxmgangdata=deg2rad(20.), \
         psfntype='gausking', \
         minmflux=3e-11, \
         maxmflux=3e-7, \
         sinddiststdv=array([.5, .5]), \
         sinddistmean=array([2., 2.]), \
         truefluxdisttype=['powr', 'powr'], \
         truenumbpnts=array([2, 3]), \
         truefluxdistslop=array([1.9, 1.1]), \
        )


def test_unbd():
    
    init( \
         pixltype='unbd', \
         exprtype='chem', \
         #verbtype=2, \
         numbswep=10000, \
         numbswepplot=10000, \
         factthin=100, \
         numbburn=0, \
         maxmnumbpnts=array([20]), \
        )


def test_empt():
    
    init( \
         numbswep=10000000, \
         emptsamp=True, \
        )


def test():
    
    init( \
         numbswep=1000, \
         factthin=900, \
         maxmnumbpnts=array([10]), \
         makeplot=False, \
         diagmode=False, \
        )

def defa():
    
    init()


globals().get(sys.argv[1])()
