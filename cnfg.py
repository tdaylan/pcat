# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def test_info():
    
    minmflux = array([3e-10, 1e-10, 3e-11, 1e-11])
    numbruns = minmflux.size
    maxmnumbpnts = zeros(numbruns, dtype=int) + 1000
    numbswep = zeros(numbruns, dtype=int) + 20000
    numbburn = numbswep / 2
    
    numbiter = minmflux.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    
    for k in range(numbiter):
        gridchan, dictpcat = init( \
                                  psfntype='doubking', \
                                  numbswep=numbswep[k], \
                                  numbburn=numbburn[k], \
                                  probprop=array([0.01, 0.01, 0., 0., 0., 0., 0., 1., 1., 0., 0., 1., 1., 1., 1.], dtype=float), \
                                  exprinfo=False, \
                                  randinit=True, \
                                  indxenerincl=arange(2, 3), \
                                  indxevttincl=arange(3, 4), \
                                  maxmnumbpnts=array([maxmnumbpnts[k]]), \
                                  maxmgang=deg2rad(10.), \
                                  minmflux=minmflux[k], \
                                  maxmflux=1e-7, \
                                  pathdata=os.environ["FERM_NGAL_DATA_PATH"], \
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
            [100, 400, 11.44, 'heal', 1, 'Reference'], \
            
            # cart
            [100, 400, 100, 'cart', 1, 'Cartesian'], \

            # maxmnumbpnts
            [100, 800, 10., 'heal', 1, '2X Max PS'], \
            
            # mocknumbpnts
            [200, 400, 10., 'heal', 1, '2X Mock PS'], \

            # numbpixl
            [100, 400, 20., 'heal', 1, '2X pixels'], \
            
            # numbener
            [100, 400, 10., 'heal', 3, '3 energy bins'], \
           ]
    numbtupl = len(tupl)
    timereal = empty(numbtupl)
    strgtupl = []
    for k in range(numbtupl):
        mocknumbpnts, maxmnumbpnts, temp, pixltype, numbener, strg = tupl[k]
        strgtupl.append(strg)
        #strgtupl.append(['%d, %d, %d, %s, %d' % (tupl[k][0], tupl[k][1], tupl[k][2], tupl[k][3], tupl[k][4])])
        if tupl[k][3] == 'heal':
            maxmgang = deg2rad(temp)
            numbsidecart = None
        else:
            maxmgang = None
            numbsidecart = temp

        if k == 0:
            continue
    
        binsenerfull = linspace(1., 1. + numbener, numbener + 1)
        gridchan, dictpcat = init( \
                                  pathdata=os.environ["PCAT_DATA_PATH"], \
                                  numbswep=numbswep, \
                                  #makeplot=False, \
                                  strgback=['unit'], \
                                  strgexpo='unit', \
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
                
    indxtupl = np.arange(numbtupl)
    size = 0.5
    figr, axis = plt.subplots()
    axis.bar(indxtupl, timereal, 2. * size)
    axis.set_ylabel('$t$ [s]')
    axis.set_xticks(indxtupl + size)
    axis.set_xticklabels(strgtupl, rotation=45)
    plt.tight_layout()
    path = tdpy.util.retr_path('pcat', onlyimag=True) + 'timereal/'
    os.system('mkdir -p %s' % path)
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    figr.savefig(path + 'timereal%04d_%s.pdf' % (log10(numbswep), strgtimestmp))
    plt.close(figr)


    
def test_psfn():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         verbtype=2, \
         numbswep=10, \
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
         pathdata=os.environ["PCAT_DATA_PATH"], \
         numbswep=2000000, \
         verbtype=1, \
         randinit=False, \
         exprinfo=False, \
         boolproppsfn=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(3, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([300]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e0, \
         maxmflux=1e4, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
        )


def test_prio():
     
    priofactdoff = array([0., 1., 2., 10.])
    for k in range(priofactdoff.size):
        init( \
             pathdata=os.environ["PCAT_DATA_PATH"], \
             numbswep=2000, \
             numbswepplot=10000, \
             numbburn=0, \
             verbtype=1, \
             randinit=False, \
             #exprinfo=False, \
             boolproppsfn=False, \
             indxenerincl=arange(1, 4), \
             indxevttincl=arange(2, 4), \
             priofactdoff=priofactdoff[k], \
             #probprop=array([1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
             strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
             strgexpr='fermflux_cmp0_ngal.fits', \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             psfntype='doubking', \
             maxmnumbpnts=array([30]), \
             maxmgang=deg2rad(10.), \
             minmflux=5e-11, \
             maxmflux=1e-7, \
             datatype='mock', \
             mocknumbpnts=array([30]), \
            )
    

def test_lowr():
      
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         numbswep=1000000, \
         verbtype=1, \
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
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=5000, \
		 numbburn=0, \
		 factthin=1, \
         randinit=False, \
         #verbtype=2, \
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


def test_numbpntsmodi():
    
    numbiter = 5
    timeatcr = empty(numbiter)
    timereal = empty(numbiter)
    timeproc = empty(numbiter)
    numbpntsmodi = arange(numbiter)
    for k in range(numbpntsmodi.size):
        gridchan, dictpcat = init( \
                                  pathdata=os.environ["PCAT_DATA_PATH"], \
	                              numbswep=1000, \
                                  numbburn=0, \
                                  factthin=1, \
                                  makeplot=False, \
	                              numbswepplot=1000, \
                                  verbtype=1, \
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
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=100000, \
		 numbswepplot=1000, \
         #verbtype=2, \
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
         mockfluxdistslop=array([-1.]), \
        )
    

def test_popl():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=100000, \
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

