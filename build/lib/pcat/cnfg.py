# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def test_info():
    
    listminmflux = logspace(-12., -8., 10)
    numbiter = listminmflux.size
    liststrgvarb = ['levi', 'info']
    for k in range(numbiter):
        seedstat = get_state()
        gdat = init( \
                    seedstat=seedstat, \
                    numbswep=10000, \
                    indxenerincl=arange(2, 3), \
                    randinit=False, \
                    proppsfp=False, \
                    indxevttincl=arange(3, 4), \
                    minmflux=listminmflux[k], \
                    maxmflux=1e-7, \
                    back=['fermisotflux.fits'], \
                    strgexpo='fermexpo_cmp0_ngal.fits', \
                   )
        for strg in liststrgvarb:
            if k == 0:
                varb = getattr(gdat, strg)
                shap = [numbiter] + list(varb.shap)
                varbdict[strg] = empty(shap)
            else:
                varbdict[k, :] = getattr(gdat, strg)

    strgtimestmp = tdpy.util.retr_strgtimestmp()
    figr, axis = plt.subplots()
    axistwin = axis.twinx()
    axis.plot(listminmflux, listinfo, label='Relative entropy')
    axis.legend(bbox_to_anchor=[0.2, 1.08], loc=2)
    axistwin.plot(listminmflux, listlevi, label='Log-evidence', color='g')
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

    numbswepcomm = 10000

    tupl = [ \
            # reference
            [50, 100, 11.44, 'heal', 1e2, 1,          numbswepcomm, 1, 'Reference'], \
            
            # numbproc
            # temp
            #[50, 100,   100, 'heal', 1e2, 1, int(numbswepcomm / 2), 2, '2 Processes'], \

            # cart
            [50, 100,   100, 'cart', 1e2, 1,          numbswepcomm, 1, 'Cartesian'], \

            # maxmnumbpnts
            [50, 200, 11.44, 'heal', 1e2, 1,          numbswepcomm, 1, '2X Max PS'], \
            
            # minmflux
            [50, 400, 11.44, 'heal', 5e1, 1,          numbswepcomm, 1, '2X Max PS, 1/2X $f_{min}$'], \
            
            # truenumbpnts
            [100, 100, 11.44, 'heal', 1e2, 1,          numbswepcomm, 1, '2X Mock PS'], \

            # numbpixl
            [50, 100, 22.88, 'heal', 1e2, 1,          numbswepcomm, 1, '2X pixels'], \
            
            # numbener
            [50, 100, 11.44, 'heal', 1e2, 3,          numbswepcomm, 1, '3 energy bins'], \
           ]
    numbtupl = len(tupl)
    indxtupl = np.arange(numbtupl)
    strgtupl = []
    liststrgvarb = ['timereal', 'timeproc', 'timeatcr', 'meanmemoresi', 'derimemoresi']
    varbdict = dict() 
    for k in range(numbtupl):
        truenumbpnts, maxmnumbpnts, temp, pixltype, minmflux, numbener, numbswep, numbproc, strg = tupl[k]
        strgtupl.append(strg)
        if pixltype == 'heal':
            maxmgang = deg2rad(temp)
            numbsidecart = None
        else:
            maxmgang = None
            numbsidecart = temp

        if numbener == 1:
            indxenerincl = arange(2, 3)
        if numbener == 3:
            indxenerincl = arange(1, 4)
        binsenerfull = linspace(1., 1. + numbener, numbener + 1)
        
        gdat = init( \
                    numbswep=numbswep, \
                    numbproc=numbproc, \
                    back=[1.], \
                    strgexpo=1., \
                    exprinfo=False, \
                    indxenerincl=indxenerincl, \
                    pixltype=pixltype, \
                    indxevttincl=arange(3, 4), \
                    psfntype='doubking', \
                    maxmnumbpnts=array([maxmnumbpnts]), \
                    maxmgang=maxmgang, \
                    minmflux=minmflux, \
                    maxmflux=1e6, \
                    numbsidecart=numbsidecart, \
                    truenumbpnts=array([mocknumbpnts]), \
                   )
        for strg in liststrgvarb:
            varbdict[strg] = getattr(gdat, strg)

    size = 0.5
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_time/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    liststrg = ['timereal', 'timeproc', 'timeprocsamp', 'timeatcr', 'timeprocnorm', 'meanmemoresi', 'derimemoresi']
    listlabl = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{MC}$', '$t_{CPU}^{\prime\prime}$ [s]', '$\bar{M}$', '$\partial_t\bar{M}$']
    listvarb = [timereal, timeproc, timeproc / numbswep, timeatcr, timeproc / (numbswep / timeatcr), meanmemoresi, derimemoresi]
    numbplot = len(liststrg)
    for k in range(numbplot):
        figr, axis = plt.subplots()
        axis.bar(indxtupl, listvarb[k], 2 * size)
        axis.set_ylabel(listlabl[k])
        axis.set_xticks(indxtupl + size)
        axis.set_xticklabels(strgtupl, rotation=45)
        plt.tight_layout()
        figr.savefig(path + '%s_%04d_%s.pdf' % (liststrg[k], log10(numbswep), strgtimestmp))
        plt.close(figr)

    
def test_psfn():
    
    tupl = [ \
            [                         None, 'singgaus',  True, 'chan'], \
            [                         None, 'singgaus', False, 'chan'], \
            ['chanfluxback_0200_4msc.fits', 'singgaus', False, 'chan'], \
            ['chanfluxback_0200_4msc.fits', 'singgaus',  True, 'chan'], \
            [                         None, 'doubking', False, 'ferm'], \
            [                         None, 'gausking', False, 'ferm'], \
            [                         None, 'doubgaus', False, 'ferm'], \
            [                         None, 'singking', False, 'ferm'], \
            [                         None, 'singgaus', False, 'ferm'], \
            [    'fermflux_cmp0_ngal.fits', 'doubking', False, 'ferm'] \
           ]
    numbtupl = len(tupl)
    for k in range(numbtupl):
        
        strgexprflux = tupl[k][0]
        psfntype = tupl[k][1]
        varioaxi = tupl[k][2]
        exprtype = tupl[k][3]
    
        if exprtype == 'ferm':
            indxenerincl = arange(2, 3)
            indxevttincl = arange(3, 4)
            back = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
            strgexpo = 'fermexpo_cmp0_ngal.fits'
            bgalcntr = deg2rad(90.)
            minmflux = 1e-7
            maxmflux = 2e-7
        else:
            indxenerincl = arange(2)
            indxevttincl = arange(1)
            back = [1.]
            strgexpo = 'chanexpo_0200_4msc.fits'
            bgalcntr = 0.
            minmflux = 1e-4
            maxmflux = 2e-4

        init( \
             numbswep=10000, \
             factthin=1, \
             numbswepplot=2000, \
             exprinfo=False, \
             indxenerincl=indxenerincl, \
             indxevttincl=indxevttincl, \
             exprtype=exprtype, \
             back=back, \
             strgexpo=strgexpo, \
             propbacp=False, \
             propcomp=False, \
             probtran=0., \
             prophypr=False, \
             psfntype=psfntype, \
             varioaxi=varioaxi, \
             bgalcntr=bgalcntr, \
             maxmnumbpnts=array([3]), \
             minmflux=minmflux, \
             maxmflux=maxmflux, \
             strgexprflux=strgexprflux, \
             truenumbpnts=array([3]), \
            )
                
    
def test_nomi():
      
    init( \
         exprinfo=False, \
         proppsfp=False, \
         indxenerincl=arange(1, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([600]), \
         maxmgang=deg2rad(20.), \
         minmflux=5e-11, \
         maxmflux=1e-7, \
         truenumbpnts=array([300]), \
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
                    exprinfo=False, \
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
                    maxmgang=deg2rad(10.), \
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


def test_uppr():
      
    init( \
         #factthin=100, \
         exprinfo=False, \
         proppsfp=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([800]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e-9, \
         maxmflux=1e-5, \
         truenumbpnts=array([400]), \
        )


def test_prio():
    
    truenumbpnts = array([200])
    priofactdoff = array([-1., -0.5, 0., 0.5, 1., 1.5, 2.])
    arry = array([1., 0.8, 1.2])
    numbiter = priofactdoff.size
    liststrgvarb = ['postnumbpnts', 'postmeanpnts']
    for k in range(numbiter):
        gdat = init( \
                    exprinfo=False, \
                    proppsfp=False, \
                    indxenerincl=arange(2, 3), \
                    indxevttincl=arange(3, 4), \
                    priofactdoff=priofactdoff[k], \
                    back=['fermisotflux.fits'], \
                    strgexpo='fermexpo_cmp0_ngal.fits', \
                    psfntype='doubking', \
                    maxmnumbpnts=array([1000]), \
                    maxmgang=deg2rad(30.), \
                    minmflux=5e-11, \
                    maxmflux=1e-7, \
                    truenumbpnts=mocknumbpnts, \
                   )
        for strg in liststrgvarb:
            if k == 0:
                varb = getattr(gdat, strg)
                shap = [numbiter] + list(varb.shap)
                varbdict[strg] = empty(shap)
            else:
                varbdict[strg] = getattr(gdat, strg)
    
    postnumbpnts = 100. * postnumbpnts / truenumbpnts[0] - 100.
    postmeanpnts = 100. * postmeanpnts / truenumbpnts[0] - 100.
    errrnumbpnts = tdpy.util.retr_errrvarb(postnumbpnts)
    errrmeanpnts = tdpy.util.retr_errrvarb(postmeanpnts)
    figr, axis = plt.subplots()
    axis.errorbar(priofactdoff, postnumbpnts[0, :], ls='', yerr=errrnumbpnts, lw=1, marker='o', label='$N$')
    axis.errorbar(priofactdoff, postmeanpnts[0, :], ls='', yerr=errrmeanpnts, lw=1, marker='o', label='$\mu$')
    axis.set_xlabel('IC Factor')
    axis.set_ylabel('$\Delta N$ [%]')
    axis.set_xlim([amin(priofactdoff) - 1., amax(priofactdoff) + 1.])
    axis.legend()
    plt.tight_layout()
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_prio/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    figr.savefig(path + 'percbias_%s.pdf' % strgtimestmp)
    plt.close(figr)


def test_lowr():
      
    init( \
         numbswep=10000, \
         verbtype=1, \
         exprinfo=False, \
         proppsfp=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexprflux='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         maxmnumbpnts=array([3]), \
         maxmgang=deg2rad(10.), \
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
             exprinfo=False, \
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
             maxmgang=deg2rad(1.), \
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
                    exprinfo=False, \
                    indxenerincl=arange(2, 3), \
                    indxevttincl=arange(3, 4), \
                    back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                    strgexpo='fermexpo_cmp0_ngal.fits', \
                    psfntype='doubking', \
                    maxmnumbpnts=array([5]), \
                    maxmgang=deg2rad(10.), \
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
             exprinfo=False, \
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
             maxmgang=deg2rad(1.), \
             minmflux=listminmflux[k], \
             maxmflux=listmaxmflux[k], \
             maxmnumbpnts=array([3]), \
             truenumbpnts=array([2]), \
             #maxmnumbpnts=array([3]), \
             #mocknumbpnts=array([2]), \
            )
        

def test_cova():
     
    init( \
         numbswep=100, \
         exprinfo=False, \
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
         exprinfo=False, \
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
         maxmgang=deg2rad(20.), \
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
