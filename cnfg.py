# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def cnfg_ferm_psfn_expr(psfntype):
     
    init( \
         numbswep=100000, \
         factthin=1, \
         trueinfo=True, \
         datatype='inpt', \
         psfntype=psfntype, \
         maxmgang=10., \
         minmflux=array([1e-8]), \
         maxmflux=array([1e-7]), \
         regitype='ngal', \
         strgexpr='fermflux_ngal.fits', \
         strgexpo='fermexpo_ngal.fits', \
         probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
        )
                
    
def cnfg_test_prio():
      
    indxenerincl = arange(3, 4)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    init(psfntype='doubking', \
         verbtype=1, \
         numbswepplot=100000, \
         numbswep=1000000, \
         numbburn=0, \
         factthin=10000, \
         randinit=False, \
         boolproppsfn=False, \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         regitype='ngal', \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         probprop=array([0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 1., 1., 1., 1.], dtype=float), \
         strgexpr='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         
         #fluxdisttype=['powr'], \
         maxmnumbpnts=array([400]), \
         maxmgang=20., \
         minmflux=3e-31, \
         maxmflux=1e-27, \
         sinddistmean=array([2.2]), \
         sinddiststdv=array([0.3]), \
         
         datatype='mock', \
         numbsideheal=256, \
         mocknumbpnts=array([200]), \
         #mockspatdisttype=['unif'], \
         #mockfluxdisttype=['powr'], \
         #mocksinddisttype=['gaus'], \
         #mockspectype=['powr'], \
         #mockfluxdistbrek=array([1e-29]), \
         #mockfluxdistsloplowr=array([-1.]), \
         #mockfluxdistslopuppr=array([2.5]), \
         mocksinddistmean=array([2.2]), \
         mocksinddiststdv=array([0.3]), \
        )


def cnfg_test_prob():
      
    indxenerincl = arange(3, 4)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    init(psfntype='doubking', \
         verbtype=1, \
         numbswepplot=50000, \
         numbswep=100000, \
         numbburn=0, \
         factthin=1000, \
         randinit=False, \
         boolproppsfn=False, \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         regitype='ngal', \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         probprop=array([1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1.], dtype=float), \
         strgexpr='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         
         #fluxdisttype=['powr'], \
         maxmnumbpnts=array([400]), \
         maxmgang=20., \
         minmflux=3e-31, \
         maxmflux=1e-27, \
         sinddistmean=array([2.2]), \
         sinddiststdv=array([0.3]), \
         
         datatype='mock', \
         numbsideheal=256, \
         mocknumbpnts=array([200]), \
        )

    
def cnfg_test():
      
    indxenerincl = arange(3, 4)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size
    minmflux = 3e-8
    maxmflux = 1e-7
    mockfluxdistslop = array([1.9])
        
    init(psfntype='doubking', \
         numbswep=10, \
         numbswepplot=1, \
         numbburn=0, \
         verbtype=2, \
         randinit=False, \
         maxmgang=20., \
         mocknumbpnts=array([1]), \
         maxmnumbpnts=array([1]), \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         probprop=array([0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], dtype=float), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         makeplot=True, \
         pathdata=os.environ["PCAT_DATA_PATH"], \
         strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
         strgexpr='fermflux_cmp0_ngal.fits', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         datatype='mock', \
         numbsideheal=256, \
         mockfluxdistslop=mockfluxdistslop, \
         mocknormback=ones((2, numbener)), \
        )
    
    
def cnfg_ferm_post():
     
    indxenerincl = arange(1, 3)
    indxevttincl = arange(2, 4)
    
    numbener = indxenerincl.size
    init(psfntype='gausking', \
		 numbswep=500000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=1.5, \
         margsize=0., \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         maxmnumbpnts=array([3]), \
         mocknumbpnts=array([3]), \
         probprop=array([0., 0., 0., 0.1, 0., 0., 0, 0, 1., 1., 1., 1.], dtype=float), \
         minmflux=array([3e-8]), \
         maxmflux=array([1e-7]), \
         regitype='ngal', \
         strgback=['fermisotflux.fits'], \
         lablback=[r'$\mathcal{I}$'], \
         nameback=['normisot'], \
         strgexpo='fermexpo_ngal_cmp0.fits', \
         stdvback=0.3, \
         stdvlbhl=0.01, \
         stdvflux=0.05, \
         stdvsind=0.05, \
         datatype='mock', \
         numbsideheal=256, \
         mockfluxdistslop=array([1.9]), \
         mocknormback=ones((1, numbener)), \
        )


def cnfg_test_spmr():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 1e-7
    mockfluxdistslop = array([1.9, 1.])
      
    init(psfntype='gausking', \
		 numbproc=1, \
		 numbswep=100, \
         makeplot=True, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=2., \
         fluxdisttype='powr', \
         verbtype=2, \
         indxenerincl=indxenerincl, \
         indxevttincl=indxevttincl, \
         maxmnumbpnts=array([3]), \
         #maxmnumbpnts=array([1000]), \
         probprop=array([0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.], dtype=float), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         datatype='mock', \
         mockfluxdisttype='powr', \
         mocknumbpnts=array([2]), \
         #mocknumbpnts=array([500]), \
         numbsideheal=256, \
         mockfluxdistslop=mockfluxdistslop, \
         mocknormback=ones((2, numbener)), \
        )


def cnfg_test_popl():
     
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 1e-7
    mockfluxdistslop = array([1.9, 1.1])
      
    init(psfntype='gausking', \
		 numbswep=100000, \
         makeplot=True, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         fluxdisttype='powr', \
         #verbtype=2, \
         indxenerincl=indxenerincl, \
         indxevttincl=arange(2, 4), \
         maxmnumbpnts=array([500, 500]), \
         minmfluxdistslop=array([1., 1.]), \
         maxmfluxdistslop=array([3., 3.]), \
         sinddiststdv=array([.5, .5]), \
         sinddistmean=array([2., 2.]), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         datatype='mock', \
         mockfluxdisttype='powr', \
         mocknumbpnts=array([300, 200]), \
         numbsideheal=256, \
         mockfluxdistslop=mockfluxdistslop, \
         mocknormback=ones((2, numbener)), \
        )


def cnfg_sdss_mock():

    indxenerincl = arange(3)
    indxevttincl = arange(1)
    numbener = indxenerincl.size
    numbener = indxenerincl.size

    init(psfntype='doubgaus', \
         numbswep=100000, \
         minmflux=array([1e3]), \
         maxmflux=array([1e5]), \
         initnumbpnts=array([100]), \
         exprtype='sdss', \
         pixltype='cart', \
         regitype='mes5', \
         stdvlbhl=2./3600., \
         lgalcntr=202., \
         bgalcntr=2., \
         radispmrlbhl=5./3600., \
         maxmgang=30./3600., \
         margsize=2./3600., \
         strgback=['unit'], \
         strgexpo='unit', \
         indxevttincl=indxevttincl, \
         indxenerincl=indxenerincl, \
         datatype='mock', \
         numbsidecart=100, \
         mockfluxdistslop=array([1.9]), \
         mocknormback=ones((1, numbener)), \
        )
    
    
def cnfg_sdss_expr():

    init(psfntype='doubgaus', \
         trueinfo=False, \
         numbswep=1000000, \
         minmflux=ones(3) * 1e3, \
         maxmflux=ones(3) * 1e5, \
         initnumbpnts=array([10]), \
         exprtype='sdss', \
         datatype='inpt', \
         pixltype='cart', \
         regitype='mes5', \
         stdvlbhl=2./3600., \
         lgalcntr=202., \
         bgalcntr=2., \
         radispmrlbhl=0.5/3600., \
         stdvflux=0.05, \
         margsize=2./3600., \
         maxmgang=30./3600., \
         strgexpr='sdssflux.fits', \
         strgexpo='sdssexpo.fits', \
         stdvback=1e-4, \
         indxevttincl=arange(1), \
         indxenerincl=arange(1), \
        )
    

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:

    pass
    #cnfg_ferm_info()
    
    #cnfg_ferm_psfn_mock('gausking')
    #cnfg_ferm_psfn_mock('doubking')

    #cnfg_ferm_psfn_expr('gausking')
    #cnfg_ferm_psfn_expr('doubking')
    
    #cnfg_ferm_post()
    #cnfg_test()
    
    #cnfg_sdss_mock()
    #cnfg_sdss_expr()

