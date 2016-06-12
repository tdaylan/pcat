# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def cnfg_ferm_psfn_expr(psfntype):
     
    init( \
         numbswep=100000, \
         factthin=1, \
         numbswepplot=50000, \
         trueinfo=True, \
         datatype='inpt', \
         psfntype=psfntype, \
         maxmgang=10., \
         minmflux=array([1e-8]), \
         maxmflux=array([1e-7]), \
         regitype='ngal', \
         strgexpr='fermflux_ngal.fits', \
         strgexpo='fermexpo_ngal.fits', \
         maxmnormback=array([5., 5.]), \
         minmnormback=array([0.2, 0.2]), \
         probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
        )
                
    
def cnfg_ferm_info():
    
    minmflux = array([1e-10, 3e-11, 1e-11, 3e-12, 1e-12])
    maxmnumbpnts = zeros(5, dtype=int) + 500
    numbswep = 20000 * array([1, 3, 10, 30, 100], dtype=int) + 10000
    numbburn = numbswep / 2
    
    numbiter = minmflux.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    
    strgexpo='fermexpo_comp_ngal.fits'
    strgexpr='fermflux_comp_ngal.fits'

    indxenerincl = arange(2, 3)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    numbproc = 1
    for k in range(numbiter):
        
        gridchan = init( \
                        psfntype='doubking', \
                        numbswep=numbswep[k], \
                        numbburn=numbburn[k], \
                        numbswepplot=50000, \
                        probprop=array([0.1, 0.1, 0., 0.1, 0., 0., 0, 0, 1., 1., 1., 1.], dtype=float), \
                        trueinfo=True, \
                        randinit=False, \
                        factthin=1000, \
                        makeplot=False, \
                        maxmgang=10., \
                        maxmnumbpnts=array([maxmnumbpnts[k]]), \
                        indxenerincl=indxenerincl, \
                        indxevttincl=indxevttincl, \
                        minmflux=array([minmflux[k]]), \
                        maxmflux=array([3e-7]), \
                        regitype='ngal', \
                        maxmnormback=array([5., 5.]), \
                        minmnormback=array([0.2, 0.2]), \
                        strgexpo=strgexpo, \
                        datatype='inpt', \
                        strgexpr=strgexpr, \
                       )
        listlevi[k] = gridchan[numbproc]
        listinfo[k] = gridchan[numbproc+1]

    plot_minmfluxinfo(minmflux, listinfo, listlevi)


def cnfg_ferm_expr_igal(strgexpr, strgexpo):
      
    init( \
         psfntype='gausking', \
         numbswep=3000000, \
         numbburn=1500000, \
         numbswepplot=50000, \
         initnumbpnts=array([100]), \
         trueinfo=True, \
         maxmgang=20., \
         indxenerincl=arange(1), \
         #indxevttincl=arange(3, 4), \
         minmflux=array([3e-11]), \
         maxmflux=array([3e-6]), \
         regitype='igal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo=strgexpo, \
         datatype='inpt', \
         strgexpr=strgexpr, \
        )
    
    
def cnfg_ferm_mock_igal():
     
    indxevttincl = arange(3, 4)
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = array([1e-11])
    maxmflux = array([1e-7])
    mockfdfnslop = array([1.9])
      
    init( \
         psfntype='singking', \
         numbswep=100000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         indxevttincl=indxevttincl, \
         indxenerincl=indxenerincl, \
         numbsideheal=256, \
         mocknumbpnts=array([800]), \
         maxmnumbpnts=array([1200]), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         mocknormback=ones((2, numbener)), \
         maxmnormback=array([2., 2.]), \
         mockfdfnslop=mockfdfnslop, \
         minmnormback=array([0.5, 0.5]), \
         strgexpo='fermexpo_cmp0_igal.fits', \
         regitype='igal', \
         datatype='mock', \
        )

    
def cnfg_ferm_expr_ngal(strgexpr='fermflux_comp_ngal.fits', strgexpo='fermexpo_comp_ngal.fits'):
    
    indxenerincl = arange(1, 4)

    minmflux = array([1e-11])
    maxmflux = array([1e-7])
        
    init(psfntype='doubking', \
         numbswep=1000000, \
         numbburn=600000, \
         proppsfn=False, \
         numbswepplot=50000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         indxenerincl=indxenerincl, \
         indxevttincl=arange(2, 4), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo=strgexpo, \
         datatype='inpt', \
         strgexpr=strgexpr, \
        )
    
    
def cnfg_ferm_post():
     
    indxenerincl = arange(1, 3)
    numbener = indxenerincl.size
    init(psfntype='doubking', \
		 numbswep=500000, \
         numbswepplot=50000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=1.5, \
         margsize=0., \
         indxevttincl=arange(3, 4), \
         indxenerincl=indxenerincl, \
         maxmnumbpnts=array([3]), \
         mocknumbpnts=array([3]), \
         probprop=array([0., 0., 0., 0.1, 0., 0., 0, 0, 1., 1., 1., 1.], dtype=float), \
         minmflux=array([3e-8]), \
         maxmflux=array([1e-7]), \
         regitype='ngal', \
         maxmnormback=array([2.]), \
         minmnormback=array([0.5]), \
         strgback=['fermisotflux.fits'], \
         lablback=[r'$\mathcal{I}$'], \
         nameback=['normisot'], \
         strgexpo='fermexpo_ngal_comp.fits', \
         stdvback=0.3, \
         stdvlbhl=0.01, \
         stdvflux=0.05, \
         stdvpropsind=0.05, \
         datatype='mock', \
         numbsideheal=256, \
         mockfdfnslop=array([1.9]), \
         mocknormback=ones((1, numbener)), \
        )


def cnfg_test_fdfn():
     
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = array([1e-11])
    maxmflux = array([1e-7])
    mockfdfnslop = array([1.9])
    mockfdfnsloplowr = array([2.9])
    mockfdfnslopuppr = array([1.9])
    mockfdfnbrek = array([3e-9])
      
    init(psfntype='doubking', \
		 numbproc=1, \
         #verbtype=3, \
		 numbswep=100000, \
         factthin=100, \
         numbswepplot=1000000, \
         makeplot=True, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         #fdfntype='brok', \
         fdfntype='powr', \
         indxevttincl=arange(3, 4), \
         indxenerincl=indxenerincl, \
         maxmnumbpnts=array([600]), \
         probprop=array([1., 1., 0., 0., 0., 0., 0, 0, 0., 0., 0., 0.], dtype=float), \
         #probprop=array([1., 1., 1., 1., 0., 0., 0., 0., 0, 0, 0., 0., 0., 0.], dtype=float), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         datatype='mock', \
         #mockfdfntype='brok', \
         mockfdfntype='powr', \
         mocknumbpnts=array([600]), \
         numbsideheal=256, \
         mockfdfnslop=mockfdfnslop, \
         mockfdfnsloplowr=mockfdfnsloplowr, \
         mockfdfnslopuppr=mockfdfnslopuppr, \
         mockfdfnbrek=mockfdfnbrek, \
         mocknormback=ones((2, numbener)), \
        )


def cnfg_test():
     
    #indxenerincl = arange(2, 3)
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = array([1e-11])
    maxmflux = array([1e-7])
    mockfdfnslop = array([1.9])
    mockfdfnsloplowr = array([2.9])
    mockfdfnslopuppr = array([1.9])
    mockfdfnbrek = array([3e-9])
      
    init(psfntype='doubking', \
		 numbproc=1, \
		 numbswep=10000, \
         #factthin=1, \
         numbswepplot=5000, \
         makeplot=True, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         fdfntype='brok', \
         #fdfntype='powr', \
         #verbtype=3, \
         indxevttincl=arange(3, 4), \
         indxenerincl=indxenerincl, \
         #maxmnumbpnts=array([3]), \
         maxmnumbpnts=array([700]), \
         probprop=array([1., 0., 0., 0., 0., 0., 1., 1., 0, 0, 1., 1., 1., 1.], dtype=float), \
         #probprop=array([1., 1., 0., 0., 0., 0., 0, 0, 0., 0., 0., 0.], dtype=float), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         #initnumbpnts=array([2]), \
         #datatype='inpt', \
         #strgexpr='fermflux_cmp0_ngal.fits', \
         datatype='mock', \
         mockfdfntype='brok', \
         #mockfdfntype='powr', \
         mocknumbpnts=array([600]), \
         #mocknumbpnts=array([100]), \
         numbsideheal=256, \
         mockfdfnslop=mockfdfnslop, \
         mockfdfnsloplowr=mockfdfnsloplowr, \
         mockfdfnslopuppr=mockfdfnslopuppr, \
         mockfdfnbrek=mockfdfnbrek, \
         mocknormback=ones((2, numbener)), \
        )


def cnfg_ferm_mock_ngal():
     
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = array([1e-11])
    maxmflux = array([1e-7])
    mockfdfnslop = array([1.9])
      
    init(psfntype='doubking', \
         numbswep=1000000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         indxevttincl=arange(2, 4), \
         indxenerincl=indxenerincl, \
         mocknumbpnts=array([500]), \
         minmflux=minmflux, \
         maxmflux=maxmflux, \
         regitype='ngal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo='fermexpo_comp_ngal.fits', \
         datatype='mock', \
         numbsideheal=256, \
         mockfdfnslop=mockfdfnslop, \
         mocknormback=ones((2, numbener)), \
        )

    
def cnfg_sdss_mock():

    indxenerincl = arange(3)
    indxevttincl = arange(1)
    numbener = indxenerincl.size
    numbener = indxenerincl.size

    init(psfntype='doubgaus', \
         numbswep=100000, \
         numbswepplot=20000, \
         minmflux=array([1e3]), \
         maxmflux=array([1e5]), \
         initnumbpnts=array([100]), \
         exprtype='sdss', \
         pixltype='cart', \
         regitype='mes5', \
         stdvlbhl=2./3600., \
         lgalcntr=202., \
         bgalcntr=2., \
         spmrlbhl=5./3600., \
         maxmnormback=array([1e3]), \
         minmnormback=array([1e2]), \
         maxmgang=30./3600., \
         margsize=2./3600., \
         strgback=['unit'], \
         strgexpo='unit', \
         indxevttincl=indxevttincl, \
         indxenerincl=indxenerincl, \
         datatype='mock', \
         numbsidecart=100, \
         mockfdfnslop=array([1.9]), \
         mocknormback=ones((1, numbener)), \
        )
    
    
def cnfg_sdss_expr():

    init(psfntype='doubgaus', \
         trueinfo=False, \
         numbswep=1000000, \
         numbswepplot=20000, \
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
         spmrlbhl=0.5/3600., \
         stdvflux=0.05, \
         maxmnormback=array([1e3]), \
         minmnormback=array([1e2]), \
         margsize=2./3600., \
         maxmgang=30./3600., \
         strgexpr='sdssflux.fits', \
         strgexpo='sdssexpo.fits', \
         stdvback=1e-4, \
         indxevttincl=arange(1), \
         indxenerincl=arange(1), \
        )
    
    
if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        name = globals().copy()
        name.update(locals())
        numbargs = len(sys.argv) - 2
        if numbargs == 0:
            print 'Running PCAT configuration %s...' % sys.argv[1]
            name.get(sys.argv[1])()
        else:
            listargs = []
            for k in range(numbargs):
                listargs.append(sys.argv[k+2])
            print 'listargs'
            print listargs

            name.get(sys.argv[1])(listargs)
    else:

        pass
        #cnfg_ferm_info()
        
        #cnfg_ferm_psfn_mock('gausking')
        #cnfg_ferm_psfn_mock('doubking')
    
        #cnfg_ferm_psfn_expr('gausking')
        #cnfg_ferm_psfn_expr('doubking')
        
        #cnfg_ferm_expr_igal('fermflux_igal_comp_time0.fits', 'fermexpo_igal_comp_time0.fits')
        #cnfg_ferm_mock_igal()
        
        #cnfg_ferm_expr_ngal('fermflux_comp_ngal.fits', 'fermexpo_comp_ngal.fits')
        #cnfg_ferm_expr_ngal('fermsflxneww_ngal.fits', 'fermsexpneww_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_cmp1_ngal.fits', 'fermexpo_cmp1_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_cmp2_ngal.fits', 'fermexpo_cmp2_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_cmp3_ngal.fits', 'fermexpo_cmp3_ngal.fits')
        
        #cnfg_ferm_post()
        #cnfg_ferm_mock_ngal()
        #cnfg_test()
        
        #cnfg_sdss_mock()
        #cnfg_sdss_expr()

