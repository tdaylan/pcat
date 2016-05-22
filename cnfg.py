# common imports
from __init__ import *

# internal functions
from main import init
from visu import *

def cnfg_ferm_psfn_expr(psfntype):
     
    init( \
         numbswep=100000, \
         factthin=1, \
         plotperd=20000, \
         trueinfo=True, \
         datatype='inpt', \
         psfntype=psfntype, \
         maxmgang=10., \
         minmspec=array([3e-10, 3e-11, 3e-12]), \
         maxmspec=array([1e-6, 1e-7, 1e-8]), \
         regitype='ngal', \
         strgexpr='fermflux_ngal.fits', \
         strgexpo='fermexpo_ngal.fits', \
         maxmnormback=array([5., 5.]), \
         minmnormback=array([0.2, 0.2]), \
         probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), \
        )
                
    
def cnfg_ferm_info():
    
    minmspec = array([3e-12, 1e-11, 3e-11, 1e-10, 3e-10])
    mocknumbpnts = 20 * array([100, 30, 10, 3, 1], dtype=int)
    numbswep = mocknumbpnts * 100 + 10000
    maxmnumbpnts = 2 * mocknumbpnts
    
    numbiter = minmspec.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)

    indxenerincl = arange(2, 3)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    for k in range(numbiter):
        
        gridchan = init( \
                        psfntype='doubking', \
                        numbswep=numbswep[k], \
                        plotperd=50000, \
                        probprop=array([0.1, 0.1, 0., 0., 1., 1., 0, 0, 1., 1., 1., 0.], dtype=float), \
                        trueinfo=True, \
                        randinit=False, \
                        factthin=1, \
                        maxmgang=20., \
                        maxmnumbpnts=array([maxmnumbpnts[k]]), \
                        indxenerincl=indxenerincl, \
                        indxevttincl=indxevttincl, \
                        minmspec=array([minmspec[k]]), \
                        maxmspec=array([3e-7]), \
                        regitype='ngal', \
                        maxmnormback=array([5., 5.]), \
                        minmnormback=array([0.2, 0.2]), \
                        strgexpo='fermexpo_cmp0_ngal.fits', \
                        datatype='mock', \
                        mocknumbpnts=array([mocknumbpnts[k]]), \
                        numbsideheal=256, \
                        mockfdfnslop=array([[1.9]]), \
                        mocknormback=ones((2, numbener)), \
                       )
        numbproc = len(gridchan)
        for l in range(numbproc):
            listchan = gridchan[l]
            listlevi[k] = listchan[14]
            listinfo[k] = listchan[15]
    plot_minmspecinfo(minmspec, listinfo, listlevi)


def cnfg_ferm_expr_igal(strgexpr, strgexpo):
      
    init( \
         psfntype='gausking', \
         numbswep=3000000, \
         numbburn=1500000, \
         plotperd=50000, \
         initnumbpnts=array([100]), \
         trueinfo=True, \
         maxmgang=20., \
         indxenerincl=arange(1), \
         #indxevttincl=arange(3, 4), \
         #minmspec=array([1e-8]), \
         #maxmspec=array([3e-6]), \
         minmspec=array([1e-8, 1e-9, 1e-10, 1e-11, 1e-12]), \
         maxmspec=array([3e-6, 3e-7, 3e-8, 3e-9, 3e-10]), \
         regitype='igal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo=strgexpo, \
         datatype='inpt', \
         strgexpr=strgexpr, \
        )
    
    
def cnfg_ferm_mock_igal():
     
    init( \
         psfntype='singking', \
         numbswep=1000000, \
         plotperd=50000, \
         numbsideheal=256, \
         minmspec=array([1e-9, 1e-10, 1e-11]), \
         maxmspec=array([1e-6, 1e-7, 1e-8]), \
         maxmnormback=array([5., 5.]), \
         minmnormback=array([0.2, 0.2]), \
         mocknormback=ones((2, 3)), \
         regitype='igal', \
         strgexpo='fermexpo_igal.fits', \
         datatype='mock', \
        )
    
    
def cnfg_ferm_expr_ngal(strgexpr='fermflux_cmp0_ngal.fits', strgexpo='fermexpo_cmp0_ngal.fits'):
    
    indxenerincl = arange(1, 4)

    minmspec = array([1e-11])
    maxmspec = array([1e-7])
        
    init(psfntype='singking', \
         numbswep=2000000, \
         numbburn=50000, \
         plotperd=50000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         indxenerincl=indxenerincl, \
         indxevttincl=arange(2, 4), \
         minmspec=minmspec, \
         maxmspec=maxmspec, \
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
		 numbswep=100000, \
         plotperd=50000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=0.5, \
         margsize=0., \
         indxevttincl=arange(3, 4), \
         indxenerincl=indxenerincl, \
         maxmnumbpnts=array([3]), \
         mocknumbpnts=array([3]), \
         probprop=array([0., 0., 0., 0.1, 0., 0., 0, 0, 1., 1., 1., 1.], dtype=float), \
         minmspec=array([3e-8]), \
         maxmspec=array([1e-7]), \
         regitype='ngal', \
         maxmnormback=array([2.]), \
         minmnormback=array([0.5]), \
         strgback=['fermisotflux.fits'], \
         lablback=[r'$\mathcal{I}$'], \
         nameback=['normisot'], \
         strgexpo='fermexpo_ngal_comp.fits', \
         stdvback=0.15, \
         stdvlbhl=0.01, \
         stdvspec=0.05, \
         stdvpropsind=0.05, \
         datatype='mock', \
         numbsideheal=256, \
         mockfdfnslop=array([[1.9]]), \
         mocknormback=ones((1, numbener)), \
        )


def cnfg_test():
     
    colrprio = True
    
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmspec = array([3e-11])
    maxmspec = array([1e-7])
    mockfdfnsloplowr = array([[1.]])
    mockfdfnslopuppr = array([[1.9]])
    mockfdfnbrek = array([[1e-10]])
      
    init(psfntype='doubking', \
		 numbproc=1, \
		 numbswep=200000, \
         plotperd=50000, \
         makeplot=True, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         fdfntype='brok', \
         mockfdfntype='brok', \
         colrprio=colrprio, \
         verbtype=1, \
         indxevttincl=arange(2, 4), \
         indxenerincl=indxenerincl, \
         maxmnumbpnts=array([300]), \
         mocknumbpnts=array([200]), \
         probprop=array([0., 0., 0.1, 0.1, 1., 1., 0, 0, 1., 1., 1., 1.], dtype=float), \
         minmspec=minmspec, \
         maxmspec=maxmspec, \
         regitype='ngal', \
         maxmnormback=array([2., 2.]), \
         minmnormback=array([0.5, 0.5]), \
         strgexpo='fermexpo_cmp0_ngal.fits', \
         #datatype='inpt', \
         #strgexpr='fermflux_cmp0_ngal.fits', \
         datatype='mock', \
         numbsideheal=256, \
         #mockfdfnslop=mockfdfnslop, \
         mockfdfnsloplowr=mockfdfnsloplowr, \
         mockfdfnslopuppr=mockfdfnslopuppr, \
         mockfdfnbrek=mockfdfnbrek, \
         mocknormback=ones((2, numbener)), \
        )


def cnfg_ferm_mock_ngal():
     
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmspec = array([3e-11])
    maxmspec = array([1e-7])
    mockfdfnslop = array([[1.9]])
      
    init(psfntype='doubking', \
         numbswep=300000, \
         plotperd=50000, \
         randinit=False, \
         trueinfo=True, \
         maxmgang=20., \
         indxevttincl=arange(2, 4), \
         indxenerincl=indxenerincl, \
         mocknumbpnts=array([300]), \
         minmspec=minmspec, \
         maxmspec=maxmspec, \
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
         plotperd=20000, \
         minmspec=array([1e3]), \
         maxmspec=array([1e5]), \
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
         mockfdfnslop=array([[1.9]]), \
         mocknormback=ones((1, numbener)), \
        )
    
    
def cnfg_sdss_expr():

    init(psfntype='doubgaus', \
         trueinfo=False, \
         numbswep=1000000, \
         plotperd=20000, \
         minmspec=ones(3) * 1e3, \
         maxmspec=ones(3) * 1e5, \
         initnumbpnts=array([10]), \
         exprtype='sdss', \
         datatype='inpt', \
         pixltype='cart', \
         regitype='mes5', \
         stdvlbhl=2./3600., \
         lgalcntr=202., \
         bgalcntr=2., \
         spmrlbhl=0.5/3600., \
         stdvspec=0.05, \
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
            print sys.argv[1]
            name.get(sys.argv[1])()
        else:
            listargs = []
            for k in range(numbargs):
                listargs.append(sys.argv[k+2])
            name.get(sys.argv[1])(listargs)
    else:

        #cnfg_ferm_info()
        
        #cnfg_ferm_psfn_mock('gausking')
        #cnfg_ferm_psfn_mock('doubking')
    
        #cnfg_ferm_psfn_expr('gausking')
        #cnfg_ferm_psfn_expr('doubking')
        
        #cnfg_ferm_expr_igal('fermflux_igal_comp_time0.fits', 'fermexpo_igal_comp_time0.fits')
        #cnfg_ferm_mock_igal()
        
        #cnfg_ferm_expr_ngal()
        #cnfg_ferm_expr_ngal('fermflux_cmp1_ngal.fits', 'fermexpo_cmp1_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_cmp2_ngal.fits', 'fermexpo_cmp2_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_cmp3_ngal.fits', 'fermexpo_cmp3_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_full_ngal.fits', 'fermexpo_full_ngal.fits')
        
        #cnfg_ferm_post()
        #cnfg_ferm_mock_ngal()
        cnfg_test()
        
        #cnfg_sdss_mock()
        #cnfg_sdss_expr()

