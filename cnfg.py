# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def test_info():
    
    minmflux = logspace(-11., -9., 5)
    numbruns = minmflux.size
    maxmnumbpnts = zeros(numbruns, dtype=int) + 1000
    numbswep = zeros(numbruns, dtype=int) + 200000
    numbburn = numbswep / 2
    
    numbiter = minmflux.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    
    for k in range(numbiter):
        gridchan, dictpcat = init( \
                                  psfntype='doubking', \
                                  numbswep=numbswep[k], \
                                  numbburn=numbburn[k], \
                                  exprinfo=False, \
                                  randinit=True, \
                                  indxenerincl=arange(1, 4), \
                                  indxevttincl=arange(3, 4), \
                                  maxmnumbpnts=array([maxmnumbpnts[k]]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=minmflux[k], \
                                  maxmflux=1e-7, \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  datatype='inpt', \
                                  strgexpr='fermflux_cmp0_ngal.fits', \
                                 )
        listlevi[k] = dictpcat['levi']
        listinfo[k] = dictpcat['info']

    plot_minmfluxinfo(minmflux, listinfo, listlevi)


def test_time():
   
    print 'Time-test suite for PCAT'

    numbswep = 200000

    tupl = [ \
            # reference
            [100, 400, 11.44, 'heal', 1, 1, 'Reference'], \
            
            # numbproc
            [100, 400, 100, 'cart', 1, 3, 'Cartesian'], \

            # cart
            [100, 400, 100, 'cart', 1, 1, 'Cartesian'], \

            # maxmnumbpnts
            [100, 800, 10., 'heal', 1, 1, '2X Max PS'], \
            
            # mocknumbpnts
            [200, 400, 10., 'heal', 1, 1, '2X Mock PS'], \

            # numbpixl
            [100, 400, 20., 'heal', 1, 1, '2X pixels'], \
            
            # numbener
            [100, 400, 10., 'heal', 3, 1, '3 energy bins'], \
           ]
    numbtupl = len(tupl)
    indxtupl = np.arange(numbtupl)
    timereal = empty(numbtupl)
    timeproc = empty(numbtupl)
    timeatcr = empty(numbtupl)
    strgtupl = []
    for k in range(numbtupl):
        mocknumbpnts, maxmnumbpnts, temp, pixltype, numbener, numbproc, strg = tupl[k]
        strgtupl.append(strg)
        if tupl[k][3] == 'heal':
            maxmgang = deg2rad(temp)
            numbsidecart = None
        else:
            maxmgang = None
            numbsidecart = temp

        binsenerfull = linspace(1., 1. + numbener, numbener + 1)
        gridchan, dictpcat = init( \
                                  numbswep=numbswep, \
                                  numbproc=numbproc, \
                                  numbburn=0, \
                                  factthin=1, \
                                  #makeplot=False, \
                                  strgback=['unit'], \
                                  strgexpo='unit', \
                                  randinit=False, \
                                  exprinfo=False, \
                                  binsenerfull=binsenerfull, \
                                  indxenerincl=arange(numbener), \
                                  pixltype=pixltype, \
                                  indxevttincl=arange(3, 4), \
                                  psfntype='doubking', \
                                  maxmnumbpnts=array([maxmnumbpnts]), \
                                  maxmgang=maxmgang, \
                                  minmflux=3e-1, \
                                  maxmflux=1e3, \
                                  datatype='mock', \
                                  numbsidecart=numbsidecart, \
                                  mocknumbpnts=array([mocknumbpnts]), \
                                 )
        timereal[k] = dictpcat['timerealtotl']
        timeproc[k] = dictpcat['timeproctotl']
        timeatcr[k] = dictpcat['timeatcr']
    
    size = 0.5
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'test_time/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    liststrg = ['timereal', 'timeproc', 'timeprocsamp', 'timeprocnorm']
    listlabl = ['$t$ [s]', '$t_{CPU}$ [s]', '$t_{CPU}^\prime$ [s]', '$t_{CPU}^{\prime\prime}$ [s]']
    listvarb = [timereal, timeproc, timeproc / numbswep, timeatcr, timeproc / numbswep / timeatcr]
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
     
    init( \
         numbswep=200000, \
         factthin=1, \
         exprinfo=False, \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         probprop=array([0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
         psfntype='doubking', \
         maxmnumbpnts=array([3]), \
         maxmgang=deg2rad(10.), \
         minmflux=3e-11, \
         maxmflux=1e-7, \
         datatype='mock', \
         mocknumbpnts=array([3]), \
        )
                
    
def test_uppr():
      
    init( \
         numbswep=200000, \
         randinit=False, \
         exprinfo=False, \
         boolproppsfn=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([300]), \
         maxmgang=deg2rad(5.), \
         minmflux=1e0, \
         maxmflux=1e4, \
         datatype='mock', \
         mocknumbpnts=array([400]), \
        )


def test_prio():
     
    mocknumbpnts = array([30])
    priofactdoff = array([-2., 0., 1., 2., 5.])
    numbiter = priofactdoff.size
    medinumbpnts = empty(numbiter)
    medimeanpnts = empty(numbiter)
    for k in range(numbiter):
        gridchan, dictpcat = init( \
		                          numbproc=1, \
                                  numbswep=200000, \
                                  numbswepplot=30000, \
                                  numbburn=0, \
                                  randinit=False, \
                                  exprinfo=False, \
                                  boolproppsfn=False, \
                                  indxenerincl=arange(1, 4), \
                                  indxevttincl=arange(2, 4), \
                                  priofactdoff=priofactdoff[k], \
                                  strgback=['fermisotflux.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  psfntype='doubking', \
                                  maxmnumbpnts=array([1000]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=5e-11, \
                                  maxmflux=1e-7, \
                                  datatype='mock', \
                                  mocknumbpnts=array([200]), \
                                 )
        medinumbpnts[k] = dictpcat['medinumbpnts']
        medimeanpnts[k] = dictpcat['medimeanpnts']
    figr, axis = plt.subplots()
    axis.plot(priofactdoff, 100. * medinumbpnts / mocknumbpnts[0] - 100., 'o', label='$N$')
    axis.plot(priofactdoff, 100. * medimeanpnts / mocknumbpnts[0] - 100., 'o', label='$\mu$')
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
         numbswep=200000, \
         randinit=False, \
         exprinfo=False, \
         boolproppsfn=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpr='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([300]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e-24, \
         maxmflux=1e-20, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
        )


def test_post():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    
    init( \
		 numbswep=200000, \
		 numbproc=1, \
         numbburn=0, \
		 factthin=1, \
         randinit=False, \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         probprop=array([0., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 1., 1., 1., 1.], dtype=float), \
         exprinfo=False, \
         strgback=['fermisotflux.fits'], \
         lablback=[r'$\mathcal{I}$'], \
         nameback=['normisot'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         stdvback=0.01, \
         stdvlbhlminm=0.01, \
         stdvlbhlmaxm=0.01, \
         stdvflux=0.05, \
         stdvsind=0.05, \
         psfntype='doubking', \
         maxmnumbpnts=array([3]), \
         maxmgang=deg2rad(1.5), \
         minmflux=3e-8, \
         maxmflux=1e-7, \
         datatype='mock', \
         mocknumbpnts=array([3]), \
         mockfluxdistslop=array([1.9]), \
        )


def test_atcr():
    
    timeatcr = empty(numbiter)
    timereal = empty(numbiter)
    timeproc = empty(numbiter)
    numbpntsmodi = array([1, 2, 3, 5, 10])
    numbiter = numbpntsmodi.size
    for k in range(numbpntsmodi.size):
        gridchan, dictpcat = init( \
	                              numbswep=200000, \
                                  numbburn=0, \
                                  factthin=1, \
                                  makeplot=False, \
                                  numbpntsmodi=numbpntsmodi[k], \
                                  randinit=False, \
                                  exprinfo=False, \
                                  indxenerincl=arange(1, 3), \
                                  indxevttincl=arange(3, 4), \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo='fermexpo_cmp0_ngal.fits', \
                                  psfntype='doubking', \
                                  maxmnumbpnts=array([200]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=3e-11, \
                                  maxmflux=3e-7, \
                                  datatype='mock', \
                                  mocknumbpnts=array([100]), \
                                 )
        timeatcr[k] = dictpcat['timeatcr']
        timereal[k] = dictpcat['timerealtotl']
        timeproc[k] = dictpcat['timeproctotl']

    plot_numbpntsmodi(numbpntsmodi, timeatcr, timereal, timeproc)
        

def test_spmr():
     
    init( \
		 numbswep=200000, \
         randinit=False, \
         exprinfo=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(2, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         probprop=array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.], dtype=float), \
         psfntype='doubking', \
         maxmnumbpnts=array([100]), \
         maxmgang=deg2rad(10.), \
         minmflux=3e-15, \
         maxmflux=3e-12, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
         minmfluxdistslop=array([-1.5]), \
         mockfluxdistslop=array([-1.]), \
        )
    

def test_popl():
     
    init( \
		 numbswep=200000, \
         randinit=False, \
         indxenerincl=arange(1, 4), \
         indxevttincl=arange(2, 4), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         maxmnumbpnts=array([500, 500]), \
         maxmgang=deg2rad(20.), \
         minmfluxdistslop=array([1., 1.]), \
         maxmfluxdistslop=array([3., 3.]), \
         sinddiststdv=array([.5, .5]), \
         sinddistmean=array([2., 2.]), \
         psfntype='gausking', \
         fluxdisttype='powr', \
         minmflux=3e-11, \
         maxmflux=3e-7, \
         datatype='mock', \
         mockfluxdisttype='powr', \
         mocknumbpnts=array([300, 200]), \
         mockfluxdistslop=array([1.9, 1.1]), \
        )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

