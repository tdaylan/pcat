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
                                  regitype='ngal', \
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


def test_psfn():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         verbtype=2, \
         numbswep=10, \
         factthin=1, \
         exprinfo=False, \
         regitype='ngal', \
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
         regitype='ngal', \
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
             numbswep=500000, \
             numbburn=0, \
             verbtype=1, \
             randinit=False, \
             #exprinfo=False, \
             boolproppsfn=False, \
             indxenerincl=arange(1, 4), \
             indxevttincl=arange(2, 4), \
             priofactdoff=priofactdoff[k], \
             regitype='ngal', \
             #probprop=array([1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
             strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
             strgexpr='fermflux_cmp0_ngal.fits', \
             strgexpo='fermexpo_cmp0_ngal.fits', \
             psfntype='doubking', \
             maxmnumbpnts=array([1000]), \
             maxmgang=deg2rad(10.), \
             minmflux=3e-11, \
             maxmflux=1e-7, \
             datatype='mock', \
             mocknumbpnts=array([300]), \
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
         regitype='ngal', \
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
     
    indxenerincl = arange(1, 3)
    indxevttincl = arange(2, 4)
    
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=500000, \
         randinit=False, \
         margsize=0., \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         probprop=array([0., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 1., 1., 1., 1.], dtype=float), \
         regitype='ngal', \
         strgback=['fermisotflux.fits'], \
         lablback=[r'$\mathcal{I}$'], \
         nameback=['normisot'], \
         strgexpo='fermexpo_ngal_cmp0.fits', \
         stdvback=0.3, \
         stdvlbhl=0.01, \
         stdvflux=0.05, \
         stdvsind=0.05, \
         psfntype='gausking', \
         maxmnumbpnts=array([3]), \
         maxmgang=deg2rad(1.5), \
         minmflux=array([3e-8]), \
         maxmflux=array([1e-7]), \
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
                                  regitype='ngal', \
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
		 numbswep=10, \
		 numbswepplot=1000, \
         verbtype=2, \
         randinit=False, \
         exprinfo=False, \
         indxenerincl=arange(2, 3), \
         indxevttincl=arange(2, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         regitype='ngal', \
         probprop=array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.], dtype=float), \
         psfntype='doubking', \
         maxmnumbpnts=array([3]), \
         maxmgang=deg2rad(10.), \
         minmflux=3e-15, \
         maxmflux=3e-12, \
         datatype='mock', \
         mocknumbpnts=array([2]), \
         mockfluxdistslop=array([-1.]), \
        )
    

def test_popl():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=100000, \
         randinit=False, \
         indxenerincl=arange(1, 4), \
         indxevttincl=arange(2, 4), \
         regitype='ngal', \
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

