# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def cnfg_test_psfn():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         numbswep=100000, \
         factthin=100, \
         exprinfo=False, \
         regitype='ngal', \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         probprop=array([0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
         psfntype='doubking', \
         maxmgang=deg2rad(10.), \
         minmflux=3e-11, \
         maxmflux=1e-7, \
         datatype='mock', \
         mocknumbpnts=array([10]), \
        )
                
    
def cnfg_test_uppr():
      
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
         strgexpo='fermexpo_cmp0_ngal.fits', \
         psfntype='doubking', \
         maxmnumbpnts=array([300]), \
         maxmgang=deg2rad(10.), \
         minmflux=1e0, \
         maxmflux=1e4, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
        )


def cnfg_test_lowr():
      
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


def cnfg_ferm_post():
     
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


def cnfg_test_spmr():
     
    init( \
         pathdata=os.environ["PCAT_DATA_PATH"], \
		 numbswep=10000, \
		 numbswepplot=6000, \
         verbtype=1, \
         randinit=False, \
         exprinfo=False, \
         indxenerincl=arange(1, 4), \
         indxevttincl=arange(2, 4), \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         regitype='ngal', \
         probprop=array([0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.], dtype=float), \
         maxmgang=deg2rad(10.), \
         maxmnumbpnts=array([300]), \
         minmflux=3e-15, \
         maxmflux=3e-12, \
         datatype='mock', \
         mocknumbpnts=array([100]), \
         mockfluxdistslop=array([0.]), \
        )
    

def cnfg_test_popl():
     
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

