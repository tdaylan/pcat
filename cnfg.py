# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def test_info():
    
    listminmflux = logspace(-12., -8., 10)
    numbiter = listminmflux.size
    maxmnumbpnts = zeros(numbiter, dtype=int) + 3000
    numbswep = zeros(numbiter, dtype=int) + 5000000
    numbburn = 4 * numbswep / 5
    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    for k in range(numbiter):
        gridchan, dictpcat = init( \
                                  modlpsfntype='doubking', \
                                  numbswep=numbswep[k], \
                                  numbburn=numbburn[k], \
                                  indxenerincl=arange(2, 3), \
                                  randinit=False, \
                                  boolproppsfn=False, \
                                  boolpropsind=False, \
                                  lgalcntr=deg2rad(0.), \
                                  bgalcntr=deg2rad(90.), \
                                  indxevttincl=arange(3, 4), \
                                  maxmnumbpnts=array([maxmnumbpnts[k]]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=listminmflux[k], \
                                  maxmflux=1e-7, \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  datatype='inpt', \
                                  strgexpr='fermflux_cmp0_ngal.fits', \
                                 )
        listlevi[k] = dictpcat['levi']
        listinfo[k] = dictpcat['info']

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

    numbswepcomm = 500

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
            
            # mocknumbpnts
            [100, 100, 11.44, 'heal', 1e2, 1,          numbswepcomm, 1, '2X Mock PS'], \

            # numbpixl
            [50, 100, 22.88, 'heal', 1e2, 1,          numbswepcomm, 1, '2X pixels'], \
            
            # numbener
            [50, 100, 11.44, 'heal', 1e2, 3,          numbswepcomm, 1, '3 energy bins'], \
           ]
    numbtupl = len(tupl)
    indxtupl = np.arange(numbtupl)
    timereal = empty(numbtupl)
    timeproc = empty(numbtupl)
    timeatcr = empty(numbtupl)
    meanmemoresi = empty(numbtupl)
    derimemoresi = empty(numbtupl)
    strgtupl = []
    for k in range(numbtupl):
        mocknumbpnts, maxmnumbpnts, temp, pixltype, minmflux, numbener, numbswep, numbproc, strg = tupl[k]
        strgtupl.append(strg)
        if pixltype == 'heal':
            maxmgang = deg2rad(temp)
            numbsidecart = None
        else:
            maxmgang = None
            numbsidecart = temp

        if numbener == 1:
            indxenerincl = arange(2, 3)
            boolpropsind = False
        if numbener == 3:
            indxenerincl = arange(1, 4)
            boolpropsind = False
        binsenerfull = linspace(1., 1. + numbener, numbener + 1)
        
        print 'k'
        print k
        
        gridchan, dictpcat = init( \
                                  # temp
                                  verbtype=0, \
                                  numbswep=numbswep, \
                                  numbproc=numbproc, \
                                  #makeplot=False, \
                                  strgback=['unit'], \
                                  strgexpo='unit', \
                                  boolpropsind=False, \
                                  exprinfo=False, \
                                  indxenerincl=indxenerincl, \
                                  pixltype=pixltype, \
                                  indxevttincl=arange(3, 4), \
                                  modlpsfntype='doubking', \
                                  maxmnumbpnts=array([maxmnumbpnts]), \
                                  maxmgang=maxmgang, \
                                  minmflux=minmflux, \
                                  maxmflux=1e6, \
                                  datatype='mock', \
                                  numbsidecart=numbsidecart, \
                                  mocknumbpnts=array([mocknumbpnts]), \
                                 )
        timereal[k] = dictpcat['timerealtotl']
        timeproc[k] = dictpcat['timeproctotl']
        timeatcr[k] = dictpcat['timeatcr']
        listmemoresi = dictpcat['listmemoresi']
        meanmemoresi[k] = mean(listmemoresi)
        derimemoresi[k] = (listmemoresi[-1] - listmemoresi[0]) / numbswep
        print 'timeatcr'
        print timeatcr[k]
        print

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
    
    tupl = [['mock', None, 'doubking', 'doubking'], \
            ['mock', None, 'gausking', 'gausking'], \
            ['mock', None, 'gausking', 'doubking'], \
            ['mock', None, 'doubking', 'gausking'], \
            ['inpt', 'fermflux_cmp0_ngal.fits', 'gausking', None], \
            ['inpt', 'fermflux_cmp0_ngal.fits', 'doubking', None]]
    numbtupl = len(tupl)
    for k in range(numbtupl):
        
        datatype = tupl[k][0]
        strgexpr = tupl[k][1]
        modlpsfntype = tupl[k][2]
        mockpsfntype = tupl[k][3]

        init( \
             numbswep=10000, \
             numbburn=0, \
             factthin=10, \
             randinit=False, \
             boolpropsind=False, \
             indxenerincl=arange(2, 3), \
             indxevttincl=arange(3, 4), \
             strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             probprop=array([0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
             modlpsfntype=modlpsfntype, \
             lgalcntr=deg2rad(0.), \
             bgalcntr=deg2rad(90.), \
             mockpsfntype=mockpsfntype, \
             maxmnumbpnts=array([100]), \
             maxmgang=deg2rad(10.), \
             minmflux=3e-11, \
             maxmflux=1e-7, \
             datatype=datatype, \
             strgexpr=strgexpr, \
             mocknumbpnts=array([100]), \
            )
                
    
def test_nomi():
      
    init( \
         numbswep=2000000, \
         exprinfo=False, \
         boolproppsfn=False, \
         boolpropsind=False, \
         indxenerincl=arange(1, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         modlpsfntype='doubking', \
         maxmnumbpnts=array([600]), \
         maxmgang=deg2rad(20.), \
         minmflux=5e-11, \
         maxmflux=1e-7, \
         datatype='mock', \
         mocknumbpnts=array([300]), \
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
    for k in range(numbtupl):
        
        specfraceval = tupl[k][0]
        binsangltype = tupl[k][1]
        numbangl = tupl[k][2]
        strgtupl[k] = '%3.1f, %s, %d' % (specfraceval, binsangltype, numbangl)
        gridchan, dictpcat = init( \
                                  numbswep=10000, \
                                  diagmode=True, \
                                  exprinfo=False, \
                                  makeanim=True, \
                                  boolproppsfn=False, \
                                  boolpropsind=False, \
                                  numbangl=numbangl, \
                                  binsangltype=binsangltype, \
                                  indxenerincl=arange(2, 3), \
                                  indxevttincl=arange(3, 4), \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  modlpsfntype='doubking', \
                                  maxmnumbpnts=array([100]), \
                                  maxmgang=deg2rad(10.), \
                                  specfraceval=specfraceval, \
                                  minmflux=1e-9, \
                                  maxmflux=1e-5, \
                                  datatype='mock', \
                                  mocknumbpnts=array([50]), \
                                 )
        stdverrrfracdimm[:, k] = dictpcat['stdverrrfracdimm'].flatten()
        stdverrrfrac[:, k] = dictpcat['stdverrrfrac'].flatten()

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
         numbswep=10000, \
         #factthin=100, \
         exprinfo=False, \
         boolproppsfn=False, \
         boolpropsind=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         modlpsfntype='doubking', \
         maxmnumbpnts=array([800]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e-9, \
         maxmflux=1e-5, \
         datatype='mock', \
         mocknumbpnts=array([400]), \
        )


def test_prio():
    
    mocknumbpnts = array([200])
    priofactdoff = array([-1., -0.5, 0., 0.5, 1., 1.5, 2.])
    numbiter = priofactdoff.size
    postnumbpnts = empty((3, numbiter))
    postmeanpnts = empty((3, numbiter))

    arry = array([1., 0.8, 1.2])
    for k in range(numbiter):
        gridchan, dictpcat = init( \
                                  # temp
                                  numbproc=1, \
                                  numbswep=100000, \
                                  exprinfo=False, \
                                  boolproppsfn=False, \
                                  boolpropsind=False, \
                                  indxenerincl=arange(2, 3), \
                                  indxevttincl=arange(3, 4), \
                                  priofactdoff=priofactdoff[k], \
                                  strgback=['fermisotflux.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  modlpsfntype='doubking', \
                                  maxmnumbpnts=array([1000]), \
                                  maxmgang=deg2rad(30.), \
                                  minmflux=5e-11, \
                                  maxmflux=1e-7, \
                                  datatype='mock', \
                                  mocknumbpnts=mocknumbpnts, \
                                 )
        postnumbpnts[:, k] = dictpcat['postnumbpnts']
        postmeanpnts[:, k] = dictpcat['postmeanpnts']
    postnumbpnts = 100. * postnumbpnts / mocknumbpnts[0] - 100.
    postmeanpnts = 100. * postmeanpnts / mocknumbpnts[0] - 100.
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
         numbswep=500000, \
         exprinfo=False, \
         boolproppsfn=False, \
         boolpropsind=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpr='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         modlpsfntype='doubking', \
         maxmnumbpnts=array([300]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e-24, \
         maxmflux=1e-20, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
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
    		 numbswep=200000, \
    		 numbproc=1, \
             numbburn=0, \
    		 factthin=1, \
             indxenerincl=indxenerincl, \
             indxevttincl=indxevttincl, \
             probprop=array([0., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 1., 1., 1., 1.], dtype=float), \
             exprinfo=False, \
             strgback=['fermisotflux.fits'], \
             lablback=[r'$\mathcal{I}$'], \
             nameback=['normisot'], \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             stdvback=stdvback, \
             stdvlbhlminm=stdvlbhlminm, \
             stdvlbhlmaxm=stdvlbhlmaxm, \
             stdvflux=stdvflux, \
             stdvsind=stdvsind, \
             modlpsfntype='doubking', \
             maxmnumbpnts=array([3]), \
             maxmgang=deg2rad(1.), \
             minmflux=minmflux, \
             maxmflux=maxmflux, \
             datatype='mock', \
             mocknumbpnts=array([3]), \
             mockfluxdistslop=array([1.9]), \
            )


def test_atcr():
    
    #listnumbpntsmodi = array([1, 2, 3, 5, 10])
    listnumbpntsmodi = array([3, 5, 10])
    numbiter = listnumbpntsmodi.size
    timereal = empty(numbiter)
    timeatcr = empty(numbiter)
    timeproc = empty(numbiter)
    numbswep = 100
    # temp
    #timeatcr = array([5670., 3420., 3042., 1023., 403.])
    #timeproc = array([103., 114., 105., 134., 140.])
    for k in range(numbiter):
        gridchan, dictpcat = init( \
	                              numbswep=numbswep, \
                                  factthin=1, \
                                  verbtype=2, \
                                  makeplot=False, \
                                  probprop=array([0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.], dtype=float), \
                                  numbpntsmodi=listnumbpntsmodi[k], \
                                  exprinfo=False, \
                                  indxenerincl=arange(2, 3), \
                                  boolpropsind=False, \
                                  indxevttincl=arange(3, 4), \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  modlpsfntype='doubking', \
                                  maxmnumbpnts=array([5]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=3e-11, \
                                  maxmflux=3e-7, \
                                  datatype='mock', \
                                  mocknumbpnts=array([2]), \
                                 )
        timeatcr[k] = dictpcat['timeatcr']
        timereal[k] = dictpcat['timerealtotl']
        timeproc[k] = dictpcat['timeproctotl']
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
	    	 numbswep=100, \
             verbtype=2, \
             factthin=1, \
             exprinfo=False, \
             indxenerincl=arange(2, 3), \
             indxevttincl=arange(3, 4), \
             strgback=['fermisotflux.fits'], \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             probprop=array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.], dtype=float), \
             modlpsfntype='doubking', \
             maxmgang=deg2rad(1.), \
             minmflux=listminmflux[k], \
             maxmflux=listmaxmflux[k], \
             datatype='mock', \
             maxmnumbpnts=array([3]), \
             mocknumbpnts=array([2]), \
             #maxmnumbpnts=array([3]), \
             #mocknumbpnts=array([2]), \
            )
        

def test_popl():
     
    init( \
		 numbswep=500000, \
         indxenerincl=arange(2, 4), \
         indxevttincl=arange(3, 4), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         maxmnumbpnts=array([500, 500]), \
         maxmgang=deg2rad(20.), \
         minmfluxdistslop=array([1., 1.]), \
         maxmfluxdistslop=array([3., 3.]), \
         sinddiststdv=array([.5, .5]), \
         sinddistmean=array([2., 2.]), \
         modlpsfntype='gausking', \
         fluxdisttype='powr', \
         minmflux=3e-11, \
         maxmflux=3e-7, \
         datatype='mock', \
         mockfluxdisttype='powr', \
         mocknumbpnts=array([300, 200]), \
         mockfluxdistslop=array([1.9, 1.1]), \
        )


def test():
    
    init( \
         numbswep=100, \
         verbtype=2, \
         numbburn=0, \
         factthin=1, \
         exprinfo=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', ], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         probprop=array([0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0.]), \
         modlpsfntype='doubking', \
         maxmnumbpnts=array([4]), \
         maxmgang=deg2rad(10.), \
         minmflux=3e-11, \
         maxmflux=1e-7, \
         datatype='mock', \
         mocknumbpnts=array([3]), \
        )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

