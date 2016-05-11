# common imports
from __init__ import *

# internal functions
from main import wrap

def retr_cnfg(               verbtype=1,               
              plotperd=50000, \
              
              makeplot=True, \
              
              diagsamp=False, \
              
              numbswep=1000000, \
              numbburn=None, \
              factthin=None, \
              
              regitype='ngal', \
              datatype='inpt', \
              randinit=True, \
              
              maxmgang=None, \
              
              minmspec=None, \
              maxmspec=None, \
              minmsind=None, \
              maxmsind=None, \
              
              minmfdfnnorm=None, \
              maxmfdfnnorm=None, \
              minmfdfnslop=None, \
              maxmfdfnslop=None, \

              fdfntype='singpowr', \
              
              psfntype='doubking', \
              
              proppsfn=True, \
              colrprio=False, \
              
              numbpopl=1, \
              
              indxevttincl=arange(2, 4), \
              indxenerincl=arange(5), \
              
              maxmnumbpnts=array([200]), \

              initnumbpnts=None, \
              trueinfo=False, \
              pntscntr=False, \
              
              numbproc=None, \
              
              liketype='pois', \
              pixltype='heal', \
              exprtype='ferm', \
              
              lgalcntr=0., \
              bgalcntr=0., \
              
              margsize=None, \

              maxmnormback=None, \
              minmnormback=None, \
              
              numbsidecart=None, \
              numbsideheal=None, \
              
              maxmangleval=None, \
              
              spmrlbhl=2., \
            
              stdvlbhl=0.1, \
              stdvback=0.04, \
              stdvspec=0.15, \
              
              stdvfdfnnorm=0.05, \
              stdvfdfnslop=0.001, \
              stdvpsfipara=0.1, \
              
              fracrand=0.05, \
              
              mocknumbpnts=None, \
              mockfdfnslop=None, \
              mocknormback=None, \
              mockpsfntype='doubking', \
              
              strgexpr=None, \
              
              strgback=None, \
              lablback=None, \
              nameback=None, \
              
              strgexpo=None, \
              
              probprop=None, \
              
             ):
    
    
    # experiment-specific defaults
    if exprtype == 'ferm':
        
        # prior boundaries
        if maxmgang == None:
            maxmgang = 20.
            
        if minmsind == None:
            minmsind = array([1.])
        if maxmsind == None:
            maxmsind = array([3.])
            
            
        if minmfdfnnorm == None:
            minmfdfnnorm = 1e0
        if maxmfdfnnorm == None:
            maxmfdfnnorm = 1e2
        if minmfdfnslop == None:
            minmfdfnslop = 1.5
        if maxmfdfnslop == None:
            maxmfdfnslop = 2.5
            
        if margsize == None:
            margsize = 1.
        
    
    cnfg = dict()

    # verbosity level
    ## 0 - no output
    ## 1 - progress
    ## 2 - sampler diagnostics
    cnfg['verbtype'] = verbtype
    
    # plot settings
    ## MCMC time period over which a frame is produced
    cnfg['plotperd'] = plotperd
    ## flag to control generation of plots
    cnfg['makeplot'] = makeplot

    # flag to diagnose the sampler
    cnfg['diagsamp'] = diagsamp
    
    # MCMC setup
    ## number of sweeps, i.e., number of samples before thinning and including burn-in
    cnfg['numbswep'] = numbswep
    ## number of burn-in samples
    cnfg['numbburn'] = numbburn
    ## number of processes
    cnfg['numbproc'] = numbproc
    ## the factor by which to thin the chain
    cnfg['factthin'] = factthin
    
    # region type
    #- ngal - NGPR (North Galactic Polar Region)
    #- igal - IG (inner Galaxy)
    cnfg['regitype'] = regitype

    # data type
    #- mock - mock data
    #- inpt - input data
    cnfg['datatype'] = datatype
    
    # likelihood function type
    cnfg['liketype'] = liketype
    
    # experiment type
    cnfg['exprtype'] = exprtype
    
    # pixelization type
    cnfg['pixltype'] = pixltype
    
    # PSF model type
    cnfg['psfntype'] = psfntype
    
    if datatype == 'mock':
        cnfg['mockpsfntype'] = mockpsfntype
        
    # color priors
    cnfg['colrprio'] = colrprio
   
    # flag to turn off PSF parameter updates
    cnfg['proppsfn'] = proppsfn

    # input data
    ## measured data
    cnfg['strgexpr'] = strgexpr
    ## background
    cnfg['strgback'] = strgback
    cnfg['lablback'] = lablback
    cnfg['nameback'] = nameback
    ## exposure
    cnfg['strgexpo'] = strgexpo
    
    # flag to use truth information
    cnfg['trueinfo'] = trueinfo
    
    # Flux distribution function model
    cnfg['fdfntype'] = fdfntype
    
    # initial state setup
    ## number of point sources
    cnfg['initnumbpnts'] = initnumbpnts
    ## flag to draw the initial state from the prior
    cnfg['randinit'] = randinit

    # energy bins to be included
    cnfg['indxenerincl'] = indxenerincl
    
    # PSF bins to be included
    cnfg['indxevttincl'] = indxevttincl
    
    # number of Point Source (PS) populations
    cnfg['numbpopl'] = numbpopl
    
    # maximum angle from the PSs to evaluate the likelihood
    cnfg['maxmangleval'] = maxmangleval
    
    # parameter limits
    ## flux distribution function normalization
    cnfg['minmfdfnnorm'] = minmfdfnnorm
    cnfg['maxmfdfnnorm'] = maxmfdfnnorm
    ## flux distribution function power law index
    cnfg['minmfdfnslop'] = minmfdfnslop
    cnfg['maxmfdfnslop'] = maxmfdfnslop
    ## flux
    cnfg['minmsind'] = minmsind
    cnfg['maxmsind'] = maxmsind
    ## spectral power-law index
    cnfg['minmspec'] = minmspec
    cnfg['maxmspec'] = maxmspec
    ## background normalizations
    cnfg['minmnormback'] = minmnormback
    cnfg['maxmnormback'] = maxmnormback

    
    # model indicator limits
    cnfg['maxmnumbpnts'] = maxmnumbpnts
    
    # Region of interest
    ## image center
    cnfg['lgalcntr'] = lgalcntr
    cnfg['bgalcntr'] = bgalcntr
    ## half of the image size
    cnfg['maxmgang'] = maxmgang
    ## image margin size
    cnfg['margsize'] = margsize

    # proposals
    ## proposal scales
    cnfg['stdvfdfnnorm'] = stdvfdfnnorm
    cnfg['stdvfdfnslop'] = stdvfdfnslop
    cnfg['stdvpsfipara'] = stdvpsfipara
    cnfg['stdvback'] = stdvback
    cnfg['stdvlbhl'] = stdvlbhl
    cnfg['stdvspec'] = stdvspec

    ## fraction of heavy-tailed proposals
    cnfg['fracrand'] = fracrand
    
    ## maximum angle over which splits and merges are proposed
    cnfg['spmrlbhl'] = spmrlbhl
    
    # mock data setup
    ## mock parameters
    cnfg['mocknumbpnts'] = mocknumbpnts
    cnfg['mockfdfnslop'] = mockfdfnslop
    cnfg['mocknormback'] = mocknormback
    ## flag to position mock point sources at the image center
    cnfg['pntscntr'] = pntscntr
    ## mock image resolution
    cnfg['numbsidecart'] = numbsidecart
    cnfg['numbsideheal'] = numbsideheal

    # proposal frequencies
    cnfg['probprop'] = probprop

    return cnfg


def cnfg_topo():
    
    cnfg = retr_cnfg( numbswep=10000, factthin=1, plotperd=20000, trueinfo=True, datatype='inpt', psfntype=psfntype, maxmgang=3., minmspec=array([3e-10, 3e-11, 3e-12]), maxmspec=array([1e-6, 1e-7, 1e-8]), regitype='ngal', strgexpr='fermflux_ngal.fits', strgexpo='fermexpo_ngal.fits', maxmnormback=array([5., 5.]), minmnormback=array([0.2, 0.2]), stdvback=0.05, probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),)
                
    wrap(cnfg)
    
    labl = ['A', '$\psi_1$', '$\phi_1$', '$f_1$', '$\psi_2$', '$\phi_2$', '$f_2$', '$\psi_3$', '$\phi_3$', '$f_3$']
    labllist = list(itertools.combinations(labl, 2))
    indxlist = list(itertools.combinations(jpara, 2))
    ncomb = len(labllist)
    
    
def cnfg_ferm_psfn_expr(psfntype):
     

    cnfg = retr_cnfg( numbswep=100000, factthin=1, plotperd=20000, trueinfo=True, datatype='inpt', psfntype=psfntype, maxmgang=10., minmspec=array([3e-10, 3e-11, 3e-12]), maxmspec=array([1e-6, 1e-7, 1e-8]), regitype='ngal', strgexpr='fermflux_ngal.fits', strgexpo='fermexpo_ngal.fits', maxmnormback=array([5., 5.]), minmnormback=array([0.2, 0.2]), stdvback=0.05, probprop=array([0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),)
                
    wrap(cnfg)
    
    
def cnfg_ferm_info():
    
    nruns = 2
    listlevi = zeros(nruns)
    listinfo = zeros(nruns)
    minmspec = logspace(-12., -7., nruns)

    for k in range(nruns):
        
        cnfg = retr_cnfg(     psfntype='gausking',     numbswep=50000,     plotperd=50000,     trueinfo=True,     maxmgang=10.,     maxmnumbpnts=array([3000]),     colrprio=True,     indxenerincl=arange(1),     indxevttincl=arange(3, 4),     minmspec=array([minmspec[k]]),     maxmspec=array([3e-7]),     regitype='ngal',     maxmnormback=array([5., 5.]),     minmnormback=array([0.2, 0.2]), strgexpo='fermexpo_ngal.fits',     stdvback=0.1,     datatype='mock',     mocknumbpnts=array([100]),     numbsideheal=256,     makeplot=False,     mocknormback=ones((2, 3)),    )
        
        gridchan = wrap(cnfg)
        numbproc = len(gridchan)
        for l in range(numbproc):
            listchan = gridchan[l]
            listlevi[k] = listchan[14]
            listinfo[k] = listchan[15]

    plot_minmspecinfo(minmspec, listinfo, listlevi)


def cnfg_ferm_expr_igal(strgexpr, strgexpo):
      
    cnfg = retr_cnfg( psfntype='gausking', numbswep=3000000, numbburn=1500000, verbtype=1, makeplot=True, plotperd=50000, initnumbpnts=array([100]), maxmnumbpnts=array([600]), trueinfo=True, maxmgang=20., colrprio=False, #indxenerincl=arange(1), \
#indxevttincl=arange(3, 4), \
#minmspec=array([1e-8]), \
#maxmspec=array([3e-6]), \
minmspec=array([1e-8, 1e-9, 1e-10, 1e-11, 1e-12]), \
maxmspec=array([3e-6, 3e-7, 3e-8, 3e-9, 3e-10]), \
regitype='igal', \
maxmnormback=array([2., 2.]), \
minmnormback=array([0.5, 0.5]), \
strgexpo=strgexpo, \
stdvback=0.1, \
datatype='inpt', \
strgexpr=strgexpr, \
                    )
        
    wrap(cnfg)
    
    
def cnfg_ferm_mock_igal():
     
    cnfg = retr_cnfg( psfntype='singking', numbswep=1000000, plotperd=50000, numbsideheal=256, maxmgang=10., minmspec=array([1e-9, 1e-10, 1e-11]), maxmspec=array([1e-6, 1e-7, 1e-8]), maxmnormback=array([5., 5.]), minmnormback=array([0.2, 0.2]), mocknormback=ones((2, 3)), regitype='igal', strgexpo='fermexpo_igal.fits', stdvback=0.05, trueinfo=True, datatype='mock')

    wrap(cnfg)
    
    
def cnfg_ferm_expr_ngal(strgexpr, strgexpo):
     
    colrprio = True
    
    indxenerincl = arange(1, 4)

    if colrprio:
        minmspec = array([1e-11])
        maxmspec = array([1e-7])
    else:
        minmspec = array([3e-9, 3e-10, 3e-11, 3e-12, 3e-13])[indxenerincl]
        maxmspec = array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9])[indxenerincl]
        
    cnfg = retr_cnfg(psfntype='doubking', \
                    numbswep=60000, \
                    numbburn=50000, 
                    verbtype=1, \
                    makeplot=True, \
                    plotperd=50000, \
                    randinit=False, \
                    proppsfn=True, \
                    #initnumbpnts=array([100]), \
                    maxmnumbpnts=array([200]), \
                    trueinfo=True, \
                    maxmgang=20., \
                    colrprio=colrprio, \
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
    wrap(cnfg)
    
    
def cnfg_ferm_post():
     
    indxenerincl = arange(1, 3)
    numbener = indxenerincl.size
    cnfg = retr_cnfg(psfntype='doubking', \
					 numbproc=1, \
					 numbswep=100000, \
                     plotperd=50000, \
                     makeplot=True, \
                     randinit=False, \
                     trueinfo=True, \
                     maxmgang=2., \
                     margsize=0., \
                     colrprio=True, \
                     verbtype=1, \
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
                     datatype='mock', \
                     numbsideheal=256, \
                     mockfdfnslop=array([[1.8]]), \
                     mocknormback=ones((1, numbener)), \
                                         )
    wrap(cnfg)


def cnfg_test():
     
    colrprio = False
    
    indxenerincl = arange(1, 3)
    numbener = indxenerincl.size

    if colrprio:
        minmspec = array([3e-11])
        maxmspec = array([1e-7])
        mockfdfnslop = array([[1.8]])
    else:
        minmspec = array([3e-9, 3e-10, 3e-11, 3e-12, 3e-13])[indxenerincl]
        maxmspec = array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9])[indxenerincl]
        mockfdfnslop = tile(array([1.8]), (1, numbener))
      
    cnfg = retr_cnfg(psfntype='doubking', \
					 numbproc=1, \
					 numbswep=10, \
                     plotperd=9, \
                     makeplot=True, \
                     randinit=False, \
                     trueinfo=True, \
                     maxmgang=20., \
                     colrprio=colrprio, \
                     verbtype=3, \
                     indxevttincl=arange(3, 4), \
                     indxenerincl=indxenerincl, \
                     maxmnumbpnts=array([3]), \
                     mocknumbpnts=array([2]), \
                     probprop=array([0., 0., 1., 0., 0., 0., 0, 0, 0., 0., 0., 0], dtype=float), \
                     minmspec=minmspec, \
                     maxmspec=maxmspec, \
                     regitype='ngal', \
                     maxmnormback=array([2., 2.]), \
                     minmnormback=array([0.5, 0.5]), \
                     strgexpo='fermexpo_ngal_comp.fits', \
                     stdvback=0.1, \
                     datatype='mock', \
                     numbsideheal=256, \
                     mockfdfnslop=mockfdfnslop, \
                     mocknormback=ones((2, numbener)), \
                                         )
    wrap(cnfg)


def cnfg_ferm_mock_ngal():
     
    colrprio = True
    
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    if colrprio:
        minmspec = array([3e-11])
        maxmspec = array([1e-7])
        mockfdfnslop = array([[1.8]])
    else:
        minmspec = array([3e-9, 3e-10, 3e-11, 3e-12, 3e-13])[indxenerincl]
        maxmspec = array([1e-5, 1e-6, 1e-7, 1e-8, 1e-9])[indxenerincl]
        mockfdfnslop = tile(array([1.8]), (1, numbener))
      
    cnfg = retr_cnfg(psfntype='doubking', \
					 numbswep=60000, \
                     plotperd=50000, \
                     makeplot=True, \
                     randinit=False, \
                     trueinfo=True, \
                     maxmgang=20., \
                     colrprio=colrprio, \
                     verbtype=1, \
                     indxevttincl=arange(2, 4), \
                     indxenerincl=indxenerincl, \
                     maxmnumbpnts=array([200]), \
                     mocknumbpnts=array([100]), \
                     minmspec=minmspec, \
                     maxmspec=maxmspec, \
                     regitype='ngal', \
                     maxmnormback=array([2., 2.]), \
                     minmnormback=array([0.5, 0.5]), \
                     strgexpo='fermexpo_ngal_comp.fits', \
                     datatype='mock', \
                     numbsideheal=256, \
                     mockfdfnslop=mockfdfnslop, \
                     mocknormback=ones((2, numbener)), \
                    )

    wrap(cnfg)
    
    
def cnfg_sdss_mock():

    
    cnfg = retr_cnfg(psfntype='doubgaus', trueinfo=False, numbswep=100000, plotperd=20000, verbtype=1, minmspec=ones(3) * 1e3, maxmspec=ones(3) * 1e5, initnumbpnts=array([100]), exprtype='sdss', datatype='mock', pixltype='cart', regitype='mes5', stdvlbhl=2./3600., lgalcntr=202., bgalcntr=2., mocknormback=ones((1, 3)), spmrlbhl=5./3600., maxmnormback=array([1e3]), minmnormback=array([1e2]), maxmgang=30./3600., numbsidecart=100, margsize=2./3600., strgexpo='sdssexpo.fits', stdvback=0.01, indxevttincl=arange(1), indxenerincl=arange(1))

    wrap(cnfg)
    
    
def cnfg_sdss_expr():

    cnfg = retr_cnfg(psfntype='doubgaus', trueinfo=False, numbswep=1000000, plotperd=20000, verbtype=1, minmspec=ones(3) * 1e3, maxmspec=ones(3) * 1e5, initnumbpnts=array([10]), maxmnumbpnts=20, exprtype='sdss', datatype='inpt', pixltype='cart', regitype='mes5', stdvlbhl=2./3600., lgalcntr=202., bgalcntr=2., spmrlbhl=0.5/3600., stdvspec=0.05, maxmnormback=array([1e3]), minmnormback=array([1e2]), margsize=2./3600., maxmgang=30./3600., strgexpr='sdssflux.fits', strgexpo='sdssexpo.fits', stdvback=1e-4, indxevttincl=arange(1), indxenerincl=arange(1))

    wrap(cnfg)
    
    
if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        eval(sys.argv[1])()
    else:

        #cnfg_ferm_info()
        
        #cnfg_ferm_psfn_mock('gausking')
        #cnfg_ferm_psfn_mock('doubking')
    
        #cnfg_ferm_psfn_expr('gausking')
        #cnfg_ferm_psfn_expr('doubking')
        
        #cnfg_ferm_expr_igal('fermflux_igal_comp_time0.fits', 'fermexpo_igal_comp_time0.fits')
        #cnfg_ferm_mock_igal()
        
        cnfg_ferm_expr_ngal('fermflux_comp_ngal.fits', 'fermexpo_comp_ngal.fits')
        #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time0.fits', 'fermexpo_ngal_comp_time0.fits')
        #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time1.fits', 'fermexpo_ngal_comp_time1.fits')
        #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time2.fits', 'fermexpo_ngal_comp_time2.fits')
        #cnfg_ferm_expr_ngal('fermflux_ngal_comp_time3.fits', 'fermexpo_ngal_comp_time3.fits')
        #cnfg_ferm_expr_ngal('fermflux_ngal_full.fits', 'fermexpo_ngal_full.fits')
        
        #cnfg_ferm_post()
        #cnfg_ferm_mock_ngal()
        #cnfg_test()
        
        #cnfg_sdss_mock()
        #cnfg_sdss_expr()

