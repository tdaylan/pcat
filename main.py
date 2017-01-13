# common imports
from __init__ import *

# internal functions
from util import *
from visu import *

def init( \
         # user interaction
         verbtype=1, \
         pathbase=os.environ["PCAT_DATA_PATH"], \

         # diagnostics
         diagmode=False, \

         # sampler
         numbswep=None, \
         numbburn=None, \
         factthin=None, \

         seedstat=None, \
         indxevttincl=None, \
         indxenerincl=None, \
         
         # empty run
         emptsamp=False, \

         # comparison with the reference catalog
         anglassc=None, \
         margfactcomp=0.9, \
         nameexpr=None, \
        
         numbspatdims=2, \
         pntstype='lght', \

         randinit=None, \
         loadvaripara=False, \
         optiprop=None, \
         regulevi=False, \
         strgexprflux=None, \
         strgcatl=None, \
         
         back=None, \
         strgback=None, \
         nameback=None, \
         lablback=None, \
         
         strgexpo=1., \
         
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         exprinfo=None, \
         pixltype=None, \
         
         # plotting
         numbswepplot=50000, \
         makeplot=True, \
         scalmaps='asnh', \
         satumaps=None, \
         makeanim=True, \
         anotcatl=False, \
         strgbinsener=None, \
         strgexprname=None, \
         strganglunit=None, \
         strganglunittext=None, \
         anglfact=None, \
         fluxfactplot=None, \
         enerfact=None, \
         
         # misc
         strgfunctime='clck', \
         strgxaxi=None, \
         strgyaxi=None, \

         # model
         ## PSF
         specfraceval=0.1, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=400, \
         strgfluxunit=None, \
         strgflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxevttfull=None, \
         binsenerfull=None, \
         maxmnumbpnts=array([1000]), \
         asymfluxprop=False, \
         psfninfoprio=True, \
         ## spectral

         # prior
         priotype='logt', \
         priofactdoff=0., \
         margfactmodl=0.9, \
         bindprio=False, \
         maxmbacp=None, \
         minmbacp=None, \
         maxmgang=None, \
         minmmeanpnts=None, \
         maxmmeanpnts=None, \
    
         spatdisttype=None, \
         spatdistslop=None, \
         
         fluxdisttype=None, \
         minmfluxdistslop=None, \
         maxmfluxdistslop=None, \
         minmfluxbrek=None, \
         maxmfluxbrek=None, \
         minmfluxdistbrek=None, \
         maxmfluxdistbrek=None, \
         minmfluxdistsloplowr=None, \
         maxmfluxdistsloplowr=None, \
         minmfluxdistslopuppr=None, \
         maxmfluxdistslopuppr=None, \
         minmsinddistmean=None, \
         maxmsinddistmean=None, \
         minmsinddiststdv=None, \
         maxmsinddiststdv=None, \

         psfntype=None, \
         varioaxi=None, \
         meansigm=None, \
         stdvsigm=None, \
         
         minmsigm=None, \
         maxmsigm=None, \
         minmgamm=None, \
         maxmgamm=None, \
         
         meangamm=None, \
         stdvgamm=None, \
         meanpsff=None, \
         stdvpsff=None, \

         minmspecsour=None, \
         maxmspecsour=None, \
         minmsizesour=None, \
         maxmsizesour=None, \
         minmellpsour=None, \
         maxmellpsour=None, \
         minmbeinhost=None, \
         maxmbeinhost=None, \
         minmspechost=None, \
         maxmspechost=None, \
         minmsizehost=None, \
         maxmsizehost=None, \
         minmellphost=None, \
         maxmellphost=None, \
         minmsherhost=None, \
         maxmsherhost=None, \
    
         spectype=None, \
         
         curvdistmean=None, \
         curvdiststdv=None, \
         
         minmflux=None, \
         maxmflux=None, \
        
         # proposals
         numbpntsmodi=1, \
         stdvprophypr=0.01, \
         stdvproppsfp=0.1, \
         stdvpropbacp=0.01, \
         stdvproplenp=1e-4, \
         varistdvlbhl=True, \
         stdvlgal=0.001, \
         stdvbgal=0.001, \
         stdvflux=0.001, \
         stdvspep=0.001, \
         stdvspmrsind=0.2, \
         probrand=0.0, \
         propfluxdist=True, \
         propsinddist=True, \
         propfluxdistbrek=True, \
         propnumbpnts=True, \
         prophypr=True, \
         proppsfp=True, \
         propbacp=True, \
         proplenp=True, \
         propcomp=True, \
         probtran=None, \
         probbrde=1., \
         radispmr=None, \

         exprvarioaxi=None, \
         exprpsfntype=None, \

         # true data
         truespatdisttype=None, \
         truespatdistslop=None, \
         truefluxdisttype=None, \
         trueminmflux=None, \
         truemaxmflux=None, \
         truefluxdistslop=None, \
         truefluxdistbrek=None, \
         truefluxdistsloplowr=None, \
         truefluxdistslopuppr=None, \
         truespectype=None, \
         truesinddistmean=None, \
         truesinddiststdv=None, \
         truevarioaxi=None, \
         truepsfntype=None, \
         trueback=None, \
         truebacp=None, \
         truelgalsour=None, \
         truebgalsour=None, \
         truespecsour=None, \
         truesizesour=None, \
         trueellpsour=None, \
         trueanglsour=None, \
         truelgalhost=None, \
         truebgalhost=None, \
         truespechost=None, \
         truesizehost=None, \
         truebeinhost=None, \
         trueellphost=None, \
         trueanglhost=None, \
         truesherhost=None, \
         truesanghost=None, \
         
         truenumbpnts=None, \
         numbsidecart=200, \
         numbsideheal=256, \
         numbdatasamp=100, \
        ):

    # construct the global object 
    gdat = tdpy.util.gdatstrt()
    for attr, valu in locals().iteritems():
        if '__' not in attr:
            setattr(gdat, attr, valu)

    # defaults
    if gdat.back == None:
        if gdat.exprtype == 'ferm':
            gdat.back = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
        elif gdat.exprtype == 'chan':
            gdat.back = ['chanfluxisot.fits']
        else:
            gdat.back = [1.]

    if gdat.lablback == None:
        gdat.lablback = [r'$\mathcal{I}$']
        if gdat.exprtype == 'ferm':
            gdat.lablback.append(r'$\mathcal{D}$')
    
    if gdat.strgback == None:
        gdat.strgback = ['fluxisot']
        if gdat.exprtype == 'ferm':
            gdat.strgback.append('fluxfdfm')
    
    if gdat.nameback == None:
        gdat.nameback = ['Isotropic']
        if gdat.exprtype == 'ferm':
            gdat.nameback.append(r'Fermi Diffuse Model')
    
    if gdat.strgexprflux == None:
        gdat.datatype = 'mock'
    else:
        gdat.datatype = 'inpt'
   
    if gdat.datatype == 'mock':
        if gdat.truenumbpnts == None:
            gdat.truenumbpopl = 1
        else:
            gdat.truenumbpopl = gdat.truenumbpnts.size
        if gdat.trueback == None:
            gdat.trueback = gdat.back
        gdat.truenumbback = len(gdat.trueback)
        gdat.trueindxpopl = arange(gdat.truenumbpopl, dtype=int)
    
    ### number of populations
    gdat.numbpopl = gdat.maxmnumbpnts.size
    gdat.numbback = len(gdat.back)
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)
    
    if gdat.indxevttincl == None:
        if gdat.exprtype == 'ferm':
            gdat.indxevttincl = arange(2, 4)
        if gdat.exprtype == 'chan':
            gdat.indxevttincl = arange(1)
    
    if gdat.indxenerincl == None:
        if gdat.exprtype == 'ferm':
            gdat.indxenerincl = arange(1, 3)
        if gdat.exprtype == 'chan':
            gdat.indxenerincl = arange(2)
    
    ### number of energy bins
    if gdat.indxenerincl != None:
        gdat.numbener = gdat.indxenerincl.size
    else:
        gdat.numbener = 1

    if gdat.randinit == None:
        if gdat.datatype == 'mock':
            gdat.randinit = False
        else:
            gdat.randinit = True

    if gdat.optiprop == None:
        if gdat.datatype == 'mock':
            gdat.optiprop = True
        else:
            gdat.optiprop = False

    # if the images are arcsinh scaled, do not saturate them
    if gdat.satumaps == None:
        if gdat.scalmaps == 'asnh':
            gdat.satumaps = False
        else:
            gdat.satumaps = True

    if gdat.strgflux == None:
        if gdat.pntstype == 'lens':
            gdat.strgflux = r'\theta_E'
        else:
            if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
                gdat.strgflux = 'f'
            if gdat.exprtype == 'sdyn':
                gdat.strgflux = 'p'
    
    if gdat.strgfluxunit == None:
        if gdat.pntstype == 'lens':
            gdat.strgfluxunit = 'arcsec'
        else:
            if gdat.exprtype == 'sdss' or gdat.exprtype == 'hubb':
                gdat.strgfluxunit = 'mag'
            if gdat.exprtype == 'ferm':
                gdat.strgfluxunit = '1/cm$^2$/s/GeV'
            if gdat.exprtype == 'chan':
                gdat.strgfluxunit = '1/cm$^2$/s/KeV'

    if gdat.indxevttfull == None:
        if gdat.exprtype == 'ferm':
            gdat.indxevttfull = arange(4)
        else:
            gdat.indxevttfull = arange(1)

    if gdat.strgenerunit == None:
        if gdat.exprtype == 'ferm':
            gdat.strgenerunit = 'GeV'
        if gdat.exprtype == 'chan':
            gdat.strgenerunit = 'KeV'
        if gdat.exprtype == 'sdyn':
            gdat.strgenerunit = ''

    if gdat.exprtype == 'ferm':
        if gdat.anglassc == None:
            gdat.anglassc = deg2rad(0.5)

    if gdat.pixltype == None:
        if gdat.exprtype == 'ferm':
            gdat.pixltype = 'heal'
        else:
            gdat.pixltype = 'cart'

    if gdat.strgexprname == None:
        if gdat.exprtype == 'chan':
            gdat.strgcatl = 'Chandra'
        if gdat.exprtype == 'ferm':
            gdat.strgexprname = 'Fermi-LAT'
        if gdat.exprtype == 'sche':
            gdat.strgexprname = 'XXXXX'
        if gdat.exprtype == 'sdyn':
            gdat.strgexprname = 'TGAS-RAVE'
    
    if gdat.strganglunit == None:
        if gdat.exprtype == 'ferm':
            gdat.strganglunit = '$^o$'
        if gdat.exprtype == 'sdyn':
            gdat.strganglunit = ''
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan' or gdat.exprtype == 'hubb':
            gdat.strganglunit = '$^{\prime\prime}$'

    if gdat.strganglunittext == None:
        if gdat.exprtype == 'ferm':
            gdat.strganglunittext = 'degree'
        if gdat.exprtype == 'sdyn':
            gdat.strganglunittext = ''
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan' or gdat.exprtype == 'hubb':
            gdat.strganglunittext = 'arcsec'
    
    if gdat.anglfact == None:
        if gdat.exprtype == 'ferm':
            gdat.anglfact = 180. / pi
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan' or gdat.exprtype == 'hubb':
            gdat.anglfact = 3600 * 180. / pi
        if gdat.exprtype == 'sche' or gdat.exprtype == 'sdyn':
            gdat.anglfact = 1.

    if gdat.fluxfactplot == None:
        if gdat.pntstype == 'lens':
            gdat.fluxfactplot = gdat.anglfact
        else:
            gdat.fluxfactplot = 1.

    if gdat.enerfact == None:
        if gdat.exprtype == 'ferm':
            gdat.enerfact = 1.
        if gdat.exprtype == 'chan':
            gdat.enerfact = 1e3

    if gdat.strgxaxi == None:
        if gdat.exprtype == 'sdyn':
            gdat.strgxaxi = r'$L_z^{\prime}$'
        else:
            if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
                gdat.strgxaxi = '$x$'
            else:
                gdat.strgxaxi = '$l$'

    if gdat.strgyaxi == None:
        if gdat.exprtype == 'sdyn':
            gdat.strgyaxi = r'$E_k^{\prime}$'
        else:
            if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
                gdat.strgyaxi = '$y$'
            else:
                gdat.strgyaxi = '$b$'

    ## experiment defaults
    if gdat.binsenerfull == None:
        if gdat.exprtype == 'ferm':
            gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
        elif gdat.exprtype == 'chan':
            gdat.binsenerfull = array([5e-4, 2e-3, 8e-3])
        else:
            gdat.binsenerfull = array([0.])
   
    # energy bin indices
    if gdat.indxenerfull == None:
        if gdat.exprtype == 'ferm':
            gdat.indxenerfull = arange(5)
        elif gdat.exprtype == 'chan':
            gdat.indxenerfull = arange(2)
        else:   
            gdat.indxenerfull = arange(gdat.binsenerfull.size - 1)
   
    # energy band string
    if gdat.strgbinsener == None and gdat.binsenerfull != None:
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'hubb':
            if gdat.exprtype == 'sdss':
                gdat.strgbinsener = ['z-band', 'i-band', 'r-band', 'g-band', 'u-band']
            if gdat.exprtype == 'hubb':
                gdat.strgbinsener = ['F606W']
        else: 
            gdat.strgbinsener = []
            for i in gdat.indxenerfull:
                gdat.strgbinsener.append('%.3g %s - %.3g %s' % (gdat.enerfact * gdat.binsenerfull[i], gdat.strgenerunit, \
                                                                                                                   gdat.enerfact * gdat.binsenerfull[i+1], gdat.strgenerunit))
   
    # flux bin for likelihood evaluation inside circles
    if gdat.exprtype == 'ferm':
        gdat.numbfluxprox = 3
    else:
        gdat.numbfluxprox = 1
    
    ## Lensing
    if gdat.exprtype == 'hubb':
        if gdat.anglassc == None:
            gdat.anglassc = 0.15 / gdat.anglfact

    ## Chandra and SDSS
    if gdat.exprtype == 'chan' or gdat.exprtype == 'sdss':
        if gdat.anglassc == None:
            gdat.anglassc = 0.5 / gdat.anglfact
    
    if gdat.exprinfo == None:
        if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
            gdat.exprinfo = True
        else:
            gdat.exprinfo = False

    if gdat.nameexpr == None:
        if gdat.exprtype == 'ferm':
            gdat.nameexpr = 'Fermi-LAT'
        if gdat.exprtype == 'sdss':
            gdat.nameexpr = 'SDSS'
        if gdat.exprtype == 'chan':
            gdat.nameexpr = 'Chandra'
        if gdat.exprtype == 'hubb':
            gdat.nameexpr = 'HST'
        # temp
        if gdat.exprtype == 'sdyn':
            gdat.nameexpr = 'Gaia'
    
    if gdat.radispmr == None:
        gdat.radispmr = 2. / gdat.anglfact
   
    if gdat.pathbase[-1] != '/':
        gdat.pathbase += '/'
    
    # paths
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    
    ## PSF class
    if gdat.indxevttincl != None:
        gdat.evttbins = True
    else:
        gdat.evttbins = False
    if gdat.evttbins:
        gdat.numbevtt = gdat.indxevttincl.size
        gdat.numbevttfull = gdat.indxevttfull.size
    else:
        gdat.numbevtt = 1
        gdat.numbevttfull = 1
        gdat.indxevttincl = array([0])
    gdat.indxevtt = arange(gdat.numbevtt)

    ## energy
    if gdat.binsenerfull.size > 1:
        gdat.enerbins = True
    else:
        gdat.enerbins = False
    if gdat.enerbins:
        gdat.numbener = gdat.indxenerincl.size
        gdat.numbenerfull = gdat.binsenerfull.size - 1
        gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
        gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
        gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
        gdat.binsener = gdat.binsenerfull[gdat.indxenerinclbins]
        gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]
        gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
        gdat.minmener = gdat.binsener[0]
        gdat.maxmener = gdat.binsener[-1]
        gdat.indxenerfull = gdat.binsenerfull.size - 1
    else:
        gdat.indxenerfull = []
        gdat.numbener = 1
        gdat.numbenerfull = 1
        gdat.indxenerincl = array([0])
        gdat.indxenerfluxdist = array([0])
        gdat.factspecener = array([1.])
    gdat.indxener = arange(gdat.numbener, dtype=int)
    gdat.indxenerfluxdist = ceil(array([gdat.numbener]) / 2.).astype(int) - 1
       
    # construct the experimental PSF
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
        gdat.exprpsfntype = 'doubking'
    if gdat.exprtype == 'chan':
        retr_chanpsfn(gdat)
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'hubb':
        retr_hubbpsfn(gdat)
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'sdyn':
        gdat.exprvarioaxi = False
        gdat.exprpsfntype = 'singgaus'
        gdat.exprpsfp = array([0.1 / gdat.anglfact])
 
    # model
    ## hyperparameters
    setp_varbfull(gdat, 'meanpnts', [1., 1e3], numbpopl=gdat.numbpopl)
    
    ### spatial
    if gdat.exprtype == 'chan':
        gdat.maxmgang = 0.492 / gdat.anglfact * gdat.numbsidecart / 2.
    if gdat.exprtype == 'ferm':
        gdat.maxmgang = 20. / gdat.anglfact
    if gdat.exprtype == 'sdyn':
        gdat.maxmgang = 1.
    if gdat.exprtype == 'hubb':
        gdat.maxmgang = 2. / gdat.anglfact
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
    
    if gdat.spatdisttype == None:
        gdat.spatdisttype = array(['unif' for l in gdat.indxpopl])
    
    ### flux
    if gdat.fluxdisttype == None:
        gdat.fluxdisttype = 'powr'

    if gdat.minmflux == None:
        if gdat.exprtype == 'ferm':
            gdat.minmflux = 3e-11
        if gdat.exprtype == 'chan':
            gdat.minmflux = 1e-7
        if gdat.exprtype == 'sdyn':
            gdat.minmflux = 1e0
        if gdat.pntstype == 'lens':
            gdat.minmflux = 2e-4
    
    if gdat.maxmflux == None:
        if gdat.exprtype == 'ferm':
            gdat.maxmflux = 1e-7
        if gdat.exprtype == 'chan':
            gdat.maxmflux = 1e-3
        if gdat.exprtype == 'sdyn':
            gdat.maxmflux = 1e4
        if gdat.pntstype == 'lens':
            gdat.maxmflux = 0.5
   
    # PS spectral model
    if gdat.spectype == None:
        if gdat.truespectype != None:
            gdat.spectype = gdat.truespectype
        else:
            gdat.spectype = array(['powr' for l in gdat.indxpopl])

    if gdat.datatype == 'mock':
        if gdat.truespectype == None:
            gdat.truespectype = ['powr' for l in gdat.trueindxpopl]
        if gdat.truevarioaxi == None:
            gdat.truevarioaxi = gdat.exprvarioaxi
        if gdat.truepsfntype == None:
            gdat.truepsfntype = gdat.exprpsfntype

    ## PSF
    if gdat.psfntype == None:
        gdat.psfntype = gdat.exprpsfntype

    if gdat.varioaxi == None:
        if gdat.exprtype == 'chan':
            gdat.varioaxi = True
        else:
            gdat.varioaxi = False

    # parameter defaults
    ## distribution
    ### flux
    setp_varbfull(gdat, 'fluxdistslop', [1.5, 3.5], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistbrek', [gdat.minmflux, gdat.maxmflux], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistsloplowr', [-1.5, 3.5], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistslopuppr', [1.5, 3.5], numbpopl=gdat.numbpopl)
    
    ### spectral index
    if gdat.numbener > 1:
        setp_varbfull(gdat, 'sinddistmean', [0.5, 3.], numbpopl=gdat.numbpopl)
        setp_varbfull(gdat, 'sinddiststdv', [0.1, 1.], numbpopl=gdat.numbpopl)
   
    # PSF
    if gdat.psfninfoprio:
        gdat.meanpsfp = gdat.exprpsfp
        gdat.stdvpsfp = 0.1 * gdat.exprpsfp
    else:
        setp_varbfull(gdat, 'sigm', [1e-2 * gdat.maxmgang, 5e-1 * gdat.maxmgang])
        setp_varbfull(gdat, 'gamm', [1.5, 20.])
        gdat.minmpsff = 0.
        gdat.maxmpsff = 1.
 
    # background parameters
    setp_varbfull(gdat, 'bacp', [0.5, 2.])
    
    # lensing
    gdat.minmlgalsour = gdat.minmlgal
    gdat.maxmlgalsour = gdat.maxmlgal
    gdat.minmbgalsour = gdat.minmbgal
    gdat.maxmbgalsour = gdat.maxmbgal
    
    
    gdat.minmlgalhost = gdat.minmlgal
    gdat.maxmlgalhost = gdat.maxmlgal
    gdat.minmbgalhost = gdat.minmbgal
    gdat.maxmbgalhost = gdat.maxmbgal
    setp_varbfull(gdat, 'specsour', array([1e-21, 1e-17]) )
    setp_varbfull(gdat, 'sizesour', [0.1 / gdat.anglfact, 1. / gdat.anglfact])
    setp_varbfull(gdat, 'ellpsour', [0., 0.3])
    setp_varbfull(gdat, 'spechost', array([1e-21, 1e-17]) )
    setp_varbfull(gdat, 'sizehost', [0.2 / gdat.anglfact, 1. / gdat.anglfact])
    setp_varbfull(gdat, 'beinhost', [0.5 / gdat.anglfact, 1. / gdat.anglfact])
    setp_varbfull(gdat, 'ellphost', [0., 0.5])
    setp_varbfull(gdat, 'sherhost', [0., 0.3])
    gdat.minmanglsour = 0.
    gdat.maxmanglsour = pi
    gdat.minmanglhost = 0.
    gdat.maxmanglhost = pi
    gdat.minmsanghost = 0.
    gdat.maxmsanghost = pi

    if gdat.maxmangl == None:
        # temp
        if gdat.exprtype == 'ferm':
            gdat.maxmangl = 20. / gdat.anglfact
        elif gdat.exprtype == 'chan':
            gdat.maxmangl = 10. / gdat.anglfact
        else:
            gdat.maxmangl = 3. * gdat.maxmgang

    # number of total sweeps
    if gdat.numbswep == None:
        gdat.numbswep = 100000

    # number of burned sweeps
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 10

    # factor by which to thin the sweeps to get samples
    if gdat.factthin == None:
        gdat.factthin = int(ceil(1e-3 * (gdat.numbswep - gdat.numbburn)))

    if gdat.strgcatl == None:
        if gdat.datatype == 'mock':
            gdat.strgcatl = 'Mock'
        else:
            if gdat.exprtype == 'ferm':
                gdat.strgcatl = '3FGL'
            else:
                gdat.strgcatl = gdat.strgexprname

    # number of processes
    if gdat.numbproc == None:
        gdat.strgproc = os.uname()[1]
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu':
            gdat.numbproc = 20
        else:
            gdat.numbproc = 1

    # conditional imports
    ## import the lensing solver by Francis-Yan if the PSs are lenses

    # get the time stamp
    gdat.strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # check the call stack for the name of the configuration function
    gdat.strgcnfg = inspect.stack()[1][3]
   
    if gdat.pntstype == 'lens':
        gdat.evalpsfnpnts = False
    else:
        gdat.evalpsfnpnts = True

    # check inputs
    if gdat.numbburn > gdat.numbswep:
        raise Exception('Bad number of burn-in sweeps.')
    if gdat.factthin > gdat.numbswep - gdat.numbburn:
        raise Exception('Bad thinning factor.')
    
    if not gdat.randinit and not gdat.exprinfo and gdat.datatype == 'inpt':
        raise Exception('If the data is provided by the user and no experimental information is given, initial state must be random.')
        
    if gdat.pixltype == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')

    if gdat.pntstype == 'lens' and (gdat.numbback > 1 or not isinstance(gdat.back[0], float)):
        print 'gdat.numbback'
        print gdat.numbback
        print 'gdat.back'
        print gdat.back
        print
        raise Exception('In a lensing problem, the background can only be uniform.')
    
    if gdat.minmflux >= gdat.maxmflux:
        raise Exception('Minimum flux is greater than maximum flux.')

    if gdat.verbtype > 0:
        print 'PCAT v0.1 started at %s' % gdat.strgtimestmp
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)
    
    # initial setup
    setpinit(gdat, True) 
   
    # create the PCAT folders
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathoutpthis = gdat.pathoutp + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    os.system('mkdir -p %s' % gdat.pathoutpthis)
    pathcatllite = gdat.pathoutpthis + 'catllite.fits'  
    pathcatl = gdat.pathoutpthis + 'catl.fits'  
    
    # generate true data
    if gdat.datatype == 'mock':
        
        if gdat.verbtype > 0:
            print 'Generating mock data...'

        if gdat.seedstat != None:
            if gdat.verbtype > 0:
                print 'Setting the seed for the RNG...'
            set_state(gdat.seedstat)

        if gdat.trueminmflux == None:
            gdat.trueminmflux = gdat.minmflux
        
        ## unit sample vector
        gdat.truesamp = zeros(gdat.truenumbpara)
   
        gdat.truefixp = zeros(gdat.truenumbfixp) + nan
        if gdat.truenumbpnts == None:
            if gdat.numbtrap > 0:
                for l in gdat.trueindxpopl:
                    gdat.truefixp[gdat.trueindxfixpnumbpnts[l]] = random_integers(0, gdat.maxmnumbpnts[l])
        else:
            gdat.truefixp[gdat.trueindxfixpnumbpnts] = gdat.truenumbpnts
    
        gdat.trueindxpntsfull = []
        if gdat.numbtrap > 0:
            for l in gdat.trueindxpopl:
                gdat.trueindxpntsfull.append(range(gdat.truenumbpnts[l]))
        else:
            gdat.trueindxpntsfull = []
        gdat.trueindxsamplgal, gdat.trueindxsampbgal, gdat.trueindxsampflux, gdat.trueindxsampspec, gdat.trueindxsampsind, gdat.trueindxsampcurv, gdat.trueindxsampexpo, \
                                                                                      gdat.trueindxsampcompcolr = retr_indx(gdat, gdat.trueindxpntsfull, gdat.truespectype)
        
        gdat.truefixp[gdat.trueindxfixpmeanpnts] = gdat.truefixp[gdat.trueindxfixpnumbpnts]

        if gdat.truespatdisttype == None:
            gdat.truespatdisttype = ['unif' for l in gdat.trueindxpopl]
        
        if gdat.truefluxdisttype == None:
            gdat.truefluxdisttype = 'powr'
     
        print 'gdat.truefluxdisttype'
        print gdat.truefluxdisttype
        # fix some true parameters
        for l in gdat.trueindxpopl:
            if gdat.truefluxdisttype == 'powr':
                defn_defa(gdat, 2., 'fluxdistslop', 'true')
            else:
                defn_defa(gdat, sqrt(gdat.trueminmflux * gdat.maxmflux), 'fluxdistbrek', 'true')
                defn_defa(gdat, 1., 'fluxdistsloplowr', 'true')
                defn_defa(gdat, 2., 'fluxdistslopuppr', 'true')
            defn_defa(gdat, 2., 'sinddistmean', 'true')
            defn_defa(gdat, 0.5, 'sinddiststdv', 'true')
       
        for k in gdat.trueindxfixp:

            # assume the true PSF
            if k in gdat.trueindxfixppsfp:
                gdat.truefixp[k] = gdat.exprpsfp[k-gdat.trueindxfixppsfp[0]]

            else:
                
                ## read input mock model parameters
                try:
                    if getattr(gdat, 'true' + gdat.namefixp[k]) != None:
                        gdat.truefixp[k] = getattr(gdat, 'true' + gdat.namefixp[k])
                except:
                    pass
                
                # randomly sample the rest of the mock model parameters
                if not isfinite(gdat.truefixp[k]):
                    gdat.truefixp[k] = icdf_fixp(gdat, 'true', rand(), k)
        
        # temp
        if gdat.pntstype == 'lens':
            gdat.truefixp[gdat.trueindxfixplgalsour] = 0.2 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalsour] = 0.2 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixplgalhost] = 0.05 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalhost] = 0.05 * gdat.maxmgang * randn()
        
        gdat.truesampvarb = empty(gdat.truenumbpara)
        gdat.truesampvarb[gdat.trueindxfixp] = gdat.truefixp
        
        gdat.truenumbpntstotl = sum(gdat.truefixp[gdat.trueindxfixpnumbpnts])
        gdat.trueindxpntstotl = arange(gdat.truenumbpntstotl)
    
        gdat.truecnts = [[] for l in gdat.trueindxpopl]
        gdat.truelgal = [[] for l in gdat.trueindxpopl]
        gdat.truebgal = [[] for l in gdat.trueindxpopl]
        gdat.truegang = [[] for l in gdat.trueindxpopl]
        gdat.trueaang = [[] for l in gdat.trueindxpopl]
        gdat.truespec = [[] for l in gdat.trueindxpopl]
        gdat.truenumbspep, gdat.trueliststrgspep, gdat.trueliststrgfluxspep = retr_numbspep(gdat.truespectype)
        if gdat.truenumbtrap > 0:
            gdat.truesind = [empty(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]) for l in gdat.trueindxpopl]
            gdat.truecurv = [empty(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]) for l in gdat.trueindxpopl]
            gdat.trueexpo = [empty(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]) for l in gdat.trueindxpopl]
        
            for l in gdat.trueindxpopl:
                if gdat.truespatdisttype[l] == 'unif':
                    gdat.truelgal[l] = icdf_self(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
                    gdat.truebgal[l] = icdf_self(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl) 

                if gdat.truespatdisttype[l] == 'disc':
                    gdat.truebgal[l] = icdf_logt(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.minmgang, gdat.factgang) * \
                                                                                                choice(array([1., -1.]), size=gdat.truefixp[gdat.trueindxfixpnumbpnts[l]])
                    gdat.truelgal[l] = icdf_self(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl) 
                if gdat.truespatdisttype[l] == 'gang':
                    gdat.truegang[l] = icdf_logt(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.minmgang, gdat.factgang)
                    gdat.trueaang[l] = icdf_self(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), 0., 2. * pi)
                    gdat.truelgal[l], gdat.truebgal[l] = retr_lgalbgal(gdat.truegang[l], gdat.trueaang[l])
                
                gdat.truespec[l] = empty((3, gdat.numbener, gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]))
                if gdat.truefluxdisttype == 'powr':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_powr(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                                                    gdat.truefixp[gdat.trueindxfixpfluxdistslop[l]])
                if gdat.truefluxdisttype == 'brok':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_brok(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), \
                                                gdat.trueminmflux, gdat.maxmflux, gdat.truefluxdistbrek[l], gdat.truefluxdistsloplowr[l], gdat.truefluxdistslopuppr[l])
                
                # temp -- make sure this reordering does not mess up other things
                gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = sort(gdat.truespec[l][:, gdat.indxenerfluxdist[0], :], axis=1)[::-1]
    
                if gdat.numbener > 1:
                    # spectral parameters
                    gdat.truesind[l] = icdf_gaus(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.truefixp[gdat.trueindxfixpsinddistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpsinddiststdv[l]])
                    if gdat.truespectype[l] == 'curv':
                        gdat.truecurv[l] = icdf_gaus(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.truecurvdistmean[l], gdat.truecurvdiststdv[l])
                    
                    if gdat.truespectype[l] == 'expo':
                        gdat.trueexpo[l] = icdf_logt(rand(gdat.truefixp[gdat.trueindxfixpnumbpnts[l]]), gdat.minmener, gdat.factener)
                
                    # spectra
                    gdat.truespec[l][:] = retr_spec(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truesind[l], \
                                                                                                             gdat.truecurv[l], gdat.trueexpo[l], gdat.truespectype[l])[None, :, :]
                
                gdat.truesampvarb[gdat.trueindxsamplgal[l]] = gdat.truelgal[l]
                gdat.truesampvarb[gdat.trueindxsampbgal[l]] = gdat.truebgal[l]
                gdat.truesampvarb[gdat.trueindxsampspec[l]] = gdat.truespec[l]
                if gdat.numbener > 1:
                    gdat.truesampvarb[gdat.trueindxsampsind[l]] = gdat.truesind[l]
                    if gdat.truespectype[l] == 'curv':
                        gdat.truesampvarb[gdat.trueindxsampcurv[l]] = gdat.truecurv[l]
                    if gdat.truespectype[l] == 'expo':
                        gdat.truesampvarb[gdat.trueindxsampexpo[l]] = gdat.trueexpo[l]
                if gdat.verbtype > 1:
                    print 'l'
                    print l
                    print 'truelgal[l]'
                    print gdat.truelgal[l]
                    print 'truebgal[l]'
                    print gdat.truebgal[l]
                    if gdat.numbener > 1:
                        print 'truesind[l]'
                        print gdat.truesind[l]
                        print 'truecurv[l]'
                        print gdat.truecurv[l]
                        print 'trueexpo[l]'
                        print gdat.trueexpo[l]
                        print 'truespectype[l]'
                        print gdat.truespectype[l]
                    print 'truespec[l]'
                    print gdat.truespec[l]
               
                if gdat.pixltype != 'unbd':
                    indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
                    gdat.truecnts[l] = gdat.truespec[l][0, :, :, None] * gdat.expo[:, indxpixltemp, :]
                    if gdat.enerbins:
                        gdat.truecnts[l] *= gdat.diffener[:, None, None]
        
        if gdat.pntstype == 'lens':
            gdat.truesourtype = 'gaus'
            gdat.truelenstype = 'SIE'
        
        proc_samp(gdat, None, 'true', raww=True)
        proc_samp(gdat, None, 'true')
            
        if gdat.makeplot:
            plot_samp(gdat, None, 'true')
        
        # temp
        if gdat.pixltype == 'unbd':
            gdat.numbdims = 2
            gdat.truedatacnts = zeros((gdat.numbener, gdat.numbdatasamp, gdat.numbevtt, gdat.numbdims))
            truelgaltemp = concatenate(gdat.truelgal)
            truebgaltemp = concatenate(gdat.truebgal)
            truefluxtemp = concatenate(gdat.truespec)[gdat.indxenerfluxdist[0], :]
            probpntsflux = truefluxtemp / sum(mockfluxtemp)

            gdat.truepsfncdfn = roll(cumsum(gdat.truepsfn, axis=1), 1)[0, :, 0]
            gdat.truepsfncdfn[0] = 0.
            gdat.truepsfncdfn /= amax(gdat.truepsfncdfn)
            truepsfnicdfintp = interp1d(gdat.truepsfncdfn, gdat.binsangl)

            if gdat.trueindxpntstotl.size > 0:
                for k in gdat.indxdatasamp:
                    indxpntsthis = choice(gdat.trueindxpntstotl, p=probpntsflux) 
                    radi = truepsfnicdfintp(rand())
                    aang = rand() * 2. * pi
                    gdat.truedatacnts[0, k, 0, 0] = truelgaltemp[indxpntsthis] + radi * sin(aang)
                    gdat.truedatacnts[0, k, 0, 1] = truebgaltemp[indxpntsthis] + radi * cos(aang)
        
        if gdat.seedstat != None:
            seed()

    # final setup
    setpfinl(gdat, True) 

    if gdat.makeplot:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.pixltype == 'cart' and gdat.pntstype == 'lght':
                    figr, axis, path = init_figr(gdat, None, 'datacntspeak', '', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, gdat.datacnts, i, m, vmin=gdat.minmdatacnts, vmax=gdat.maxmdatacnts)
                    make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                    axis.scatter(gdat.anglfact * gdat.lgalcart[gdat.indxxaximaxm], gdat.anglfact * gdat.bgalcart[gdat.indxyaximaxm], alpha=0.6, s=20, edgecolor='none')
                    
                    plt.tight_layout()
                    plt.savefig(path)
                    plt.close(figr)
    
                if gdat.correxpo:
                    figr, axis, path = init_figr(gdat, None, 'expo', '', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, gdat.expo, i, m)
                    make_cbar(gdat, axis, imag, i)
                    plt.tight_layout()
                    plt.savefig(path)
                    plt.close(figr)
            
                    for c in gdat.indxback:
                        figr, axis, path = init_figr(gdat, None, 'backcnts', '', indxenerplot=i, indxevttplot=m)
                        imag = retr_imag(gdat, axis, gdat.backcnts[c], i, m, vmin=gdat.minmdatacnts, vmax=gdat.maxmdatacnts)
                        make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                        plt.tight_layout()
                        plt.savefig(path)
                        plt.close(figr)
        
                    if gdat.numbback > 1:
                        figr, axis, path = init_figr(gdat, None, 'backcntstotl', '', indxenerplot=i, indxevttplot=m)
                        imag = retr_imag(gdat, axis, gdat.backcntstotl, i, m, vmin=gdat.minmdatacnts, vmax=gdat.maxmdatacnts)
                        make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                        plt.tight_layout()
                        plt.savefig(path)
                        plt.close(figr)
            
                    figr, axis, path = init_figr(gdat, None, 'diffcntstotl', '', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, gdat.datacnts - gdat.backcntstotl, i, m, vmin=gdat.minmdatacnts, vmax=gdat.maxmdatacnts)
                    make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                    plt.tight_layout()
                    plt.savefig(path)
                    plt.close(figr)
    
                if gdat.pntstype == 'lens':
                    
                    figr, axis, path = init_figr(gdat, None, 'modlcntsraww', 'true', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, gdat.truemodlcntsraww, 0, 0, vmin=gdat.minmdatacnts, vmax=gdat.maxmdatacnts)
                    make_cbar(gdat, axis, imag, 0, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                    plt.tight_layout()
                    plt.savefig(path)
                    plt.close(figr)
    
    # write the list of arguments to file
    # temp
    #fram = inspect.currentframe()
    #listargs, temp, temp, listargsvals = inspect.getargvalues(fram)
    #fileargs = open(gdat.pathoutpthis + 'args.txt', 'w')
    #fileargs.write('PCAT call arguments')
    #for args in listargs:
    #    fileargs.write('%s = %s\n' % (args, listargsvals[args]))
    #fileargs.close()

    # start the timer
    gdat.timerealtotl = time.time()
    gdat.timeproctotl = time.clock()
   
    if gdat.verbtype > 1:
        print 'minmflux'
        print gdat.minmflux
        print 'maxmflux'
        print gdat.maxmflux
        if gdat.pixltype != 'unbd':
            print 'minmcnts'
            print gdat.minmcnts
            print 'maxmcnts'
            print gdat.maxmcnts
        print 'indxsampnumbpnts: ', gdat.indxfixpnumbpnts
        print 'indxsampmeanpnts: ', gdat.indxfixpmeanpnts
        print 'indxsampfluxdistslop: ', gdat.indxfixpfluxdistslop
        print 'indxsampfluxdistbrek: ', gdat.indxfixpfluxdistbrek
        print 'indxsampfluxdistsloplowr: ', gdat.indxfixpfluxdistsloplowr
        print 'indxsampfluxdistslopuppr: ', gdat.indxfixpfluxdistslopuppr
        print 'indxsamppsfp: ', gdat.indxfixppsfp
        print 'indxsampbacp: '
        print gdat.indxfixpbacp
        if gdat.evalcirc:
            print 'maxmangleval'
            print gdat.anglfact * gdat.maxmangleval, ' ', gdat.strganglunit
        print
        if gdat.trueinfo:
            print 'Truth information'
            print 'truelgal'
            for l in gdat.indxpopl:
                print gdat.truelgal[l]
            print 'truebgal'
            for l in gdat.indxpopl:
                print gdat.truebgal[l]
            print 'truespec'
            for l in gdat.indxpopl:
                print gdat.truespec[l]
            print 'truepsfp'
            print gdat.truepsfp
            print
            
    # make initial plots
    if gdat.makeplot:
        #plot_3fgl_thrs(gdat)
        if gdat.pixltype != 'unbd':
            plot_datacntshist(gdat)
            if gdat.pntstype == 'lght':
                plot_indxprox(gdat)
        #if gdat.exprtype == 'ferm':
        #    plot_fgl3(gdat)
        # temp
        #plot_intr()
        #plot_plot()
        #plot_king(gdat)
        if gdat.evalcirc:
            plot_eval(gdat)
        #if gdat.datatype == 'mock':
        #    plot_pntsdiff()

    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdat, 'gdat')

    # lock the global object againts any future modifications
    gdat.lockmodi()

    gdat.timereal = zeros(gdat.numbproc)
    gdat.timeproc = zeros(gdat.numbproc)
    if gdat.numbproc == 1:
        listgdatmodi = [work(gdat, 0)]
    else:
        if gdat.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        gdat.lock = mp.Manager().Lock()
        
        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(work, gdat)
        listgdatmodi = pool.map(workpart, gdat.indxproc)

        pool.close()
        pool.join()

    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        timeinit = gdat.functime()

    # unlock the global object
    gdat.unlkmodi()
    
    # aggregate samples from the chains
    ## list of parameters to be gathered
    gdat.liststrgchan = []
    for attr, valu in listgdatmodi[0].__dict__.iteritems():
        if attr.startswith('list') and not isinstance(valu, list):
            gdat.liststrgchan.append(attr[4:])
    
    for strg in gdat.liststrgchan:
        for k in gdat.indxproc:
            if k == 0:
                shap = getattr(listgdatmodi[k], 'list' + strg).shape
                shap = [shap[0], gdat.numbproc] + list(shap[1:])
                temp = zeros(shap) - 1
            if len(shap) > 2:
                temp[:, k, :] = getattr(listgdatmodi[k], 'list' + strg)
            else:
                temp[:, k] = getattr(listgdatmodi[k], 'list' + strg)
        setattr(gdat, 'list' + strg, temp)
    
    gdat.maxmllikswep = empty(gdat.numbproc)
    gdat.indxswepmaxmllik = empty(gdat.numbproc, dtype=int)
    gdat.sampvarbmaxmllik = empty((gdat.numbproc, gdat.numbpara))
    gdat.maxmlposswep = empty(gdat.numbproc)
    gdat.indxswepmaxmlpos = empty(gdat.numbproc, dtype=int)
    gdat.sampvarbmaxmlpos = empty((gdat.numbproc, gdat.numbpara))
    for k in gdat.indxproc:
        gdat.maxmllikswep[k] = listgdatmodi[k].maxmllikswep
        gdat.indxswepmaxmllik[k] = listgdatmodi[k].indxswepmaxmllik
        gdat.sampvarbmaxmllik[k] = listgdatmodi[k].sampvarbmaxmllik
        gdat.maxmlposswep[k] = listgdatmodi[k].maxmlposswep
        gdat.indxswepmaxmlpos[k] = listgdatmodi[k].indxswepmaxmlpos
        gdat.sampvarbmaxmlpos[k] = listgdatmodi[k].sampvarbmaxmlpos

    # Gelman-Rubin test
    if gdat.numbproc > 1:
        gdat.gmrbfixp = zeros(gdat.numbfixp)
        gdat.gmrbstat = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        if gdat.verbtype > 0:
            print 'Computing the Gelman-Rubin TS...'
            timeinit = gdat.functime()
        for k in gdat.indxfixp:
            gdat.gmrbfixp[k] = tdpy.mcmc.gmrb_test(gdat.listsampvarb[:, :, gdat.indxfixp[k]])
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxevtt:
                    gdat.gmrbstat[i, j, m] = tdpy.mcmc.gmrb_test(gdat.listmodlcnts[:, :, i, j, m])
        if gdat.verbtype > 0:
            timefinl = gdat.functime()
            print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # flatten the chain output
    ## PS indices
    gdat.listindxpntsfull = []
    for j in gdat.indxsamp:      
        for k in gdat.indxproc:
            gdat.listindxpntsfull.append(listgdatmodi[k].listindxpntsfull[j])

    ## list of other parameters to be flattened
    gdat.liststrgchanflat = deepcopy(gdat.liststrgchan)
    for strg in ['deltlpri', 'deltllik', 'memoresi']:
        gdat.liststrgchanflat.remove(strg)
    
    ## other parameters
    for strg in gdat.liststrgchanflat:
        inpt = getattr(gdat, 'list' + strg)
        shap = [inpt.shape[0] * inpt.shape[1]] + list(inpt.shape[2:])
        setattr(gdat, 'list' + strg, inpt.reshape(shap))
    
    # add execution times to the chain output
    for k in gdat.indxproc:
        gdat.timereal[k] = listgdatmodi[k].timereal
        gdat.timeproc[k] = listgdatmodi[k].timeproc
    
    # correct the likelihoods for the constant data dependent factorial
    llikoffs = sum(sp.special.gammaln(gdat.datacnts + 1))
    gdat.listllik -= llikoffs
    gdat.maxmllikswep -= llikoffs
    gdat.maxmlposswep -= llikoffs
    
    # find the maximum likelihood and posterior over the chains
    gdat.indxprocmaxmllik = argmax(gdat.maxmllikswep)
    gdat.maxmllikswep = gdat.maxmllikswep[gdat.indxprocmaxmllik]
    gdat.indxswepmaxmllik = gdat.indxprocmaxmllik * gdat.numbsamp + gdat.indxswepmaxmllik[gdat.indxprocmaxmllik]
    gdat.sampvarbmaxmllik = gdat.sampvarbmaxmllik[gdat.indxprocmaxmllik, :]
    
    gdat.indxprocmaxmlpos = argmax(gdat.maxmlposswep)
    gdat.maxmlposswep = gdat.maxmlposswep[gdat.indxprocmaxmlpos]
    gdat.indxswepmaxmlpos = gdat.indxswepmaxmlpos[gdat.indxprocmaxmlpos]
    gdat.sampvarbmaxmlpos = gdat.sampvarbmaxmlpos[gdat.indxprocmaxmlpos, :]

    # calculate log-evidence using the harmonic mean estimator
    if gdat.verbtype > 0:
        print 'Estimating the Bayesian evidence...'
        timeinit = gdat.functime()
    
    if gdat.regulevi:
        # regularize the harmonic mean estimator
        ## get an ellipse around the median posterior 
        gdat.elpscntr = percentile(listsamp, 50., axis=0)
        thissamp = rand(gdat.numbpara) * 1e-6
        stdvpara = ones(gdat.numbpara) * 1e-6
        limtpara = zeros((2, gdat.numbpara))
        limtpara[1, :] = 1.
        ## find the samples that lie inside the ellipse
        elpsaxis, minmfunc = tdpy.util.minm(thissamp, retr_elpsfrac, stdvpara=stdvpara, limtpara=limtpara, tolrfunc=1e-6, verbtype=gdat.verbtype, optiprop=True)
        lnorregu = -0.5 * gdat.numbpara * log(pi) + sp.special.gammaln(0.5 * gdat.numbpara + 1.) - sum(elpsaxis)
        indxsampregu = 0
        listlliktemp = listllik[indxsampregu]
    else:
        listlliktemp = gdat.listllik
    gdat.levi = retr_levi(listlliktemp)
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # calculate relative entropy
    gdat.info = retr_info(gdat.listllik, gdat.levi)

    # parse the sample vector
    gdat.listfixp = gdat.listsampvarb[:, gdat.indxfixp]

    # post process samples
    
    if gdat.numbtrap > 0:
        # collect PS parameters from the chains
        gdat.listlgal = [[] for l in gdat.indxpopl]
        gdat.listbgal = [[] for l in gdat.indxpopl]
        gdat.listspec = [[] for l in gdat.indxpopl]
        gdat.listflux = [[] for l in gdat.indxpopl]
        if gdat.numbener > 1:
            gdat.listsind = [[] for l in gdat.indxpopl]
            gdat.listcurv = [[] for l in gdat.indxpopl]
            gdat.listexpo = [[] for l in gdat.indxpopl]
        for n in gdat.indxsamptotl: 
            for l in gdat.indxpopl:
                indxsamplgal, indxsampbgal, indxsampflux, indxsampspec, indxsampsind, indxsampcurv, \
                                                                    indxsampexpo, indxsampcompcolr = retr_indx(gdat, gdat.listindxpntsfull[n], gdat.spectype)
                gdat.listlgal[l].append(gdat.listsampvarb[n, indxsamplgal[l]])
                gdat.listbgal[l].append(gdat.listsampvarb[n, indxsampbgal[l]])
                gdat.listspec[l].append(gdat.listsampvarb[n, indxsampspec[l]])
                gdat.listflux[l].append(gdat.listsampvarb[n, indxsampspec[l]][gdat.indxenerfluxdist[0], :])
                if gdat.numbener > 1:
                    gdat.listsind[l].append(gdat.listsampvarb[n, indxsampsind[l]])
                    if gdat.spectype[l] == 'curv':
                        gdat.listcurv[l].append(gdat.listsampvarb[n, indxsampcurv[l]])
                    if gdat.spectype[l] == 'expo':
                        gdat.listexpo[l].append(gdat.listsampvarb[n, indxsampexpo[l]])
        
        # calculate the PS counts, radial and azimuthal positions 
        if gdat.pntstype == 'lght':
            gdat.listcnts = [[] for l in gdat.indxpopl]
        gdat.listgang = [[] for l in gdat.indxpopl]
        gdat.listaang = [[] for l in gdat.indxpopl]
        for l in gdat.indxpopl:
            for n in gdat.indxsamptotl:
                gdat.listgang[l].append(retr_gang(gdat.listlgal[l][n], gdat.listbgal[l][n]))
                gdat.listaang[l].append(retr_aang(gdat.listlgal[l][n], gdat.listbgal[l][n]))
                if gdat.pntstype == 'lght':
                    gdat.listcnts[l].append(retr_pntscnts(gdat, gdat.listlgal[l][n], gdat.listbgal[l][n], gdat.listspec[l][n]))

    ## bin PS parameters
    if gdat.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        timeinit = gdat.functime()
    
    gdat.listlgalhist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numblgal))
    gdat.listbgalhist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbgal))
    gdat.listspechist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbinsplot, gdat.numbener))
    if gdat.pntstype == 'lght':
        gdat.listcntshist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbinsplot, gdat.numbener, gdat.numbevtt))
    if gdat.numbener > 1:
        gdat.listspephist = []
        gdat.listfluxspephist = []
        for l in gdat.indxpopl:
            gdat.listspephist.append(empty((gdat.numbsamptotl, gdat.numbspep[l], gdat.numbspepbins)))
            gdat.listfluxspephist.append(empty((gdat.numbsamptotl, gdat.numbspep[l], gdat.numbbinsplot, gdat.numbspepbins)))
    gdat.listganghist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbgang))
    gdat.listaanghist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbaang))
    if gdat.numbtrap > 0:
        for n in gdat.indxsamptotl: 
            for l in gdat.indxpopl:
                # find the indices of the model PSs that are in the comparison area
                indxmodlpntscomp = retr_indxpntscomp(gdat, gdat.listlgal[l][n], gdat.listbgal[l][n])
        
                gdat.listlgalhist[n, l, :] = histogram(gdat.listlgal[l][n][indxmodlpntscomp], gdat.binslgal)[0]
                gdat.listbgalhist[n, l, :] = histogram(gdat.listbgal[l][n][indxmodlpntscomp], gdat.binsbgal)[0]
                for i in gdat.indxener:
                    gdat.listspechist[n, l, :, i] = histogram(gdat.listspec[l][n][i, indxmodlpntscomp], gdat.binsspecplot[i, :])[0]
                    if gdat.pntstype == 'lght':
                        for m in gdat.indxevtt:
                            gdat.listcntshist[n, l, :, i, m] = histogram(gdat.listcnts[l][n][i, indxmodlpntscomp, m], gdat.binscnts[i, :])[0]
                if gdat.numbener > 1:
                    for p in gdat.indxspep[l]:
                        gdat.listspephist[l][n, p, :] = histogram(gdat.listspep[l][n][indxmodlpntscomp, p], gdat.binsspep[:, p])[0]
                        gdat.listfluxspephist[l][n, p, :, :] = histogram2d(gdat.listspec[l][n][gdat.indxenerfluxdist[0], indxmodlpntscomp], \
                                                                                   gdat.listspep[l][n][indxmodlpntscomp, p], [gdat.binsfluxplot, gdat.binsspep[:, p]])[0]
                gdat.listganghist[n, l, :] = histogram(gdat.listgang[l][n][indxmodlpntscomp], gdat.binsgang)[0]
                gdat.listaanghist[n, l, :] = histogram(gdat.listaang[l][n][indxmodlpntscomp], gdat.binsaang)[0]
        
        gdat.pntsprob = zeros((gdat.numbpopl, gdat.numbbgalpntsprob, gdat.numblgalpntsprob, gdat.numbbinsplot))
        for l in gdat.indxpopl:
            temparry = concatenate([gdat.listlgal[l][n] for n in gdat.indxsamptotl])
            temp = empty((len(temparry), 3))
            temp[:, 0] = temparry
            temp[:, 1] = concatenate([gdat.listbgal[l][n] for n in gdat.indxsamptotl])
            temp[:, 2] = concatenate([gdat.listspec[l][n][gdat.indxenerfluxdist[0], :] for n in gdat.indxsamptotl])
            gdat.pntsprob[l, :, :, :] = histogramdd(temp, bins=(gdat.binslgalpntsprob, gdat.binsbgalpntsprob, gdat.binsfluxplot))[0]

    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # calculate the autocorrelation of the chains
    if gdat.verbtype > 0:
        print 'Computing the autocorrelation of the chains...'
        timeinit = gdat.functime()
   
    gdat.atcr, gdat.timeatcr = tdpy.mcmc.retr_timeatcr(gdat.listmodlcnts, verbtype=gdat.verbtype)
    if gdat.timeatcr == 0.:
        print 'Autocorrelation time estimation failed.'

    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    ## construct a deterministic catalog
    # temp
    if gdat.strgcnfg == 'test':
        retr_detrcatl(gdat)
    
    # construct lists of samples for each proposal type
    gdat.listindxsamptotlprop = []
    gdat.listindxsamptotlpropaccp = []
    gdat.listindxsamptotlpropreje = []
    for n in gdat.indxprop:
        gdat.listindxsamptotlprop.append(where(gdat.listindxprop == gdat.indxprop[n])[0])
        gdat.listindxsamptotlpropaccp.append(intersect1d(where(gdat.listindxprop == gdat.indxprop[n])[0], where(gdat.listaccp)[0]))
        gdat.listindxsamptotlpropreje.append(intersect1d(where(gdat.listindxprop == gdat.indxprop[n])[0], where(logical_not(gdat.listaccp))[0]))

    ## cross correlation with the reference catalog
    if gdat.numbtrap > 0 and gdat.trueinfo:
        
        if gdat.verbtype > 0:
            print 'Cross-correlating the sample catalogs with the reference catalog...'
        timeinit = gdat.functime()
        
        gdat.listspecassc = []
        for l in gdat.indxpopl:
            gdat.listspecassc.append(zeros((gdat.numbsamptotl, gdat.numbener, gdat.truefixp[gdat.indxfixpnumbpnts[l]])))
            for n in gdat.indxsamptotl:
                indxmodl, trueindxpntsassc = corr_catl(gdat, None, l, gdat.listlgal[l][n], gdat.listbgal[l][n], gdat.listspec[l][n], None)
                indxpntstrue = where(indxmodl >= 0)[0]
                for i in gdat.indxener:
                    gdat.listspecassc[l][n, i, indxpntstrue] = gdat.listspec[l][n][i, indxmodl[indxpntstrue]]
        
        timefinl = gdat.functime()
        if gdat.verbtype > 0:
            print 'Done in %.3g seconds.' % (timefinl - timeinit)
   
    
    ## spatial mean of maps
    gdat.liststrgmean = ['modlcnts']
    if gdat.pntstype == 'lght':
        gdat.liststrgmean.append('pntsflux')
    for strg in gdat.liststrgmean:
        tempinpt = getattr(gdat, 'list' + strg)
        shap = tempinpt.shape
        shap = [gdat.numbsamptotl, gdat.numbener]
        temp = empty(shap)
        for n in gdat.indxsamptotl:
            temp[n, :] = sum(sum(tempinpt[n, :, :, :] * gdat.expo, 2), 1) / sum(sum(gdat.expo, 2), 1)
        setattr(gdat, 'list' + strg + 'mean', temp)

    # compute credible intervals
    # temp
    gdat.liststrgpostsaveodim = ['modlcnts', 'resicnts']
    if gdat.numbtrap > 0:
        if gdat.pntstype == 'lght':
            gdat.liststrgpostsaveodim += ['pntsflux']
    
    # temp
    gdat.liststrgchan += ['specassc', 'lgalhist', 'bgalhist', 'spechist', 'fixp']
    if gdat.numbener > 1:
        gdat.liststrgchan += ['spephist', 'fluxspephist']
    liststrgexcl = []#['boolreje', 'indxprop', 'lacc']
    for strg in gdat.liststrgchan:
        if not (strg in liststrgexcl):
            listtemp = getattr(gdat, 'list' + strg)
            if isinstance(listtemp, list):
                posttemp = []
                meditemp = []
                errrtemp = []
                for k in range(len(listtemp)):
                    posttempsing = tdpy.util.retr_postvarb(listtemp[k])
                    meditempsing = posttempsing[0, :]
                    errrtempsing = tdpy.util.retr_errrvarb(posttempsing)
                    posttemp.append(posttempsing)
                    meditemp.append(meditempsing)
                    errrtemp.append(errrtempsing)
            else:
                posttemp = tdpy.util.retr_postvarb(listtemp)
                meditemp = posttemp[0, ...]
                errrtemp = tdpy.util.retr_errrvarb(posttemp)
            setattr(gdat, 'post' + strg, posttemp)
            setattr(gdat, 'medi' + strg, meditemp)
            setattr(gdat, 'errr' + strg, errrtemp)
    
    # temp
    if gdat.evalcirc and gdat.correxpo:
        gdat.medicntsbackfwhm = retr_cntsbackfwhm(gdat, gdat.postfixp[0, gdat.indxfixpbacp], gdat.postfwhm[0, :])
        gdat.medibinssigm = retr_sigm(gdat, gdat.binscnts, gdat.medicntsbackfwhm)
   
    # memory usage
    gdat.meanmemoresi = mean(gdat.listmemoresi, 1)
    gdat.derimemoresi = (gdat.meanmemoresi[-1] - gdat.meanmemoresi[0]) / gdat.numbswep
    
    # temp
    if False:
        # write the PCAT output to disc
        if gdat.verbtype > 0:
            print 'Writing the PCAT output to %s...' % pathcatl
  
        shel = shelve.open(pathcatl, 'n')
        for attr, valu in locals().iteritems():
            if attr == 'gdat':
                shel[attr] = valu
        shel.close()

        if gdat.makeplot:
            plot_post(pathcatl=pathcatl, verbtype=gdat.verbtype, makeanim=gdat.makeanim)
    else:
        if gdat.makeplot:
            plot_post(gdat=gdat, writ=False)

    gdat.timerealtotl = time.time() - gdat.timerealtotl
    gdat.timeproctotl = time.clock() - gdat.timeproctotl
     
    if gdat.verbtype > 0:
        for k in gdat.indxproc:
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdat.timereal[k], gdat.timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (gdat.timerealtotl, gdat.timeproctotl)
        if gdat.makeplot:
            print 'The plots are in ' + gdat.pathplot
        print 'The ensemble of catalogs is at ' + pathcatl
    return gdat
    
    
def work(gdat, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for this chain
    seed()
    
    # empty object to hold chain-specific variables that will be modified by the chain
    gdatmodi = tdpy.util.gdatstrt()
  
    # temporary data structures
    gdatmodi.tempfluxdistslop = empty(gdat.numbpopl)
    gdatmodi.tempfluxdistbrek = empty(gdat.numbpopl)
    gdatmodi.tempfluxdistsloplowr = empty(gdat.numbpopl)
    gdatmodi.tempfluxdistslopuppr = empty(gdat.numbpopl)
    gdatmodi.tempsinddistmean = empty(gdat.numbpopl)
    gdatmodi.tempsinddiststdv = empty(gdat.numbpopl)
    gdatmodi.templpri = zeros(gdat.numblpri)
    gdatmodi.templlik = zeros_like(gdat.datacnts)
    gdatmodi.tempspep = empty((gdat.numbpara, 2))

    # data structure to hold the indices of model PS to be compared to the reference catalog 
    gdatmodi.indxmodlpntscomp = [[] for l in gdat.indxpopl]
    
    # construct the initial state
    if gdat.verbtype > 0:
        print 'Initializing the sampler state...'
   
    ## unit sample vector
    gdatmodi.drmcsamp = zeros((gdat.numbpara, 2))
   
    ## Fixed-dimensional parameters
    for k in gdat.indxfixp:
        if gdat.randinit or not isfinite(gdat.truefixp[k]) or k in gdat.indxfixpnumbpnts and gdat.truefixp[k] > gdat.maxmnumbpntstotl:
            if k in gdat.indxfixpnumbpnts:
                gdatmodi.drmcsamp[gdat.indxfixp[k], 0] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[k] + 1))
            else:
                gdatmodi.drmcsamp[gdat.indxfixp[k], 0] = rand()
        else:
            gdatmodi.drmcsamp[gdat.indxfixp[k], 0] = cdfn_fixp(gdat, gdat.truefixp[k], k)

    ## lists of occupied and empty transdimensional parameters
    thisnumbpnts = gdatmodi.drmcsamp[gdat.indxfixpnumbpnts, 0].astype(int)
    gdatmodi.thisindxpntsfull = []
    gdatmodi.thisindxpntsempt = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            gdatmodi.thisindxpntsfull.append(range(thisnumbpnts[l]))
            gdatmodi.thisindxpntsempt.append(range(thisnumbpnts[l], gdat.maxmnumbpnts[l]))
    else:
        gdatmodi.thisindxpntsfull = []
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampflux, gdatmodi.thisindxsampspec, gdatmodi.thisindxsampsind, \
                                    gdatmodi.thisindxsampcurv, gdatmodi.thisindxsampexpo, gdatmodi.thisindxsampcompcolr = retr_indx(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
        
    if gdat.numbtrap > 0:
        ## PS components
        if gdat.randinit:
            randinittemp = True
        else:
            try:
                for l in gdat.indxpopl:
                    gdatmodi.drmcsamp[gdatmodi.thisindxsamplgal[l], 0] = copy(cdfn_self(gdat.truelgal[l], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl))
                    gdatmodi.drmcsamp[gdatmodi.thisindxsampbgal[l], 0] = copy(cdfn_self(gdat.truebgal[l], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl))
                    if gdat.fluxdisttype == 'powr':
                        fluxdistslop = icdf_atan(gdatmodi.drmcsamp[gdat.indxfixpfluxdistslop[l], 0], gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
                        fluxunit = cdfn_flux_powr(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.minmflux, gdat.maxmflux, fluxdistslop)
                    if gdat.fluxdisttype == 'brok':
                        flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                        fluxdistbrek = icdf_logt(gdatmodi.drmcsamp[gdat.indxfixpfluxdistbrek[l], 0], gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
                        fluxdistsloplowr = icdf_atan(gdatmodi.drmcsamp[gdat.indxfixpfluxdistsloplowr[l], 0], gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
                        fluxdistslopuppr = icdf_atan(gdatmodi.drmcsamp[gdat.indxfixpfluxdistslopuppr[l], 0], gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])
                        fluxunit = cdfn_flux_brok(flux, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
                    gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :], 0] = copy(fluxunit)

                    if gdat.numbener > 1:
                        # color parameters
                        gdatmodi.drmcsamp[gdatmodi.thisindxsampsind[l], 0] = cdfn_gaus(gdat.truesind[l], gdat.truefixp[gdat.indxfixpsinddistmean[l]], \
                                                                                                                                gdat.truefixp[gdat.indxfixpsinddiststdv[l]])
                        if gdat.spectype[l] == 'curv':
                            gdatmodi.drmcsamp[gdatmodi.thisindxsampcurv[l], 0] = cdfn_gaus(gdat.truecurv[l], gdat.curvdistmean[l], gdat.curvdiststdv[l])
                        if gdat.spectype[l] == 'expo':
                            gdatmodi.drmcsamp[gdatmodi.thisindxsampexpo[l], 0] = cdfn_logt(gdat.trueexpo[l], gdat.minmener, gdat.factener)
                
                randinittemp = False
            except:
                
                randinittemp = True
                print 'Reference catalog is inappropriate for deterministic initial state. Seeding the initial state randomly...'
       
        if randinittemp:
            for l in gdat.indxpopl:
                gdatmodi.drmcsamp[gdatmodi.thisindxsampcompcolr[l], 0] = rand(gdatmodi.thisindxsampcompcolr[l].size)

    if gdat.verbtype > 1:
        print 'drmcsamp'
        for k in gdat.indxpara:
            print gdatmodi.drmcsamp[k, :]
    
    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.)[0] + gdat.numbpopl
    indxsampbadduppr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.)[0] + gdat.numbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'Initial unit sample vector went outside the unit interval...'
        gdatmodi.drmcsamp[indxsampbaddlowr, 0] = 0.
        gdatmodi.drmcsamp[indxsampbadduppr, 0] = 1.

    ## sample vector
    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.drmcsamp[:, 0], 'this')
    
    ## sample index
    gdatmodi.cntrswep = 0
   
    proc_samp(gdat, gdatmodi, 'this')
   
    ## initial predicted count maps
    if gdat.pntstype == 'lght':
        # temp
        if gdat.boolintpanglcosi:
            binsangltemp = gdat.binsanglcosi
        else:
            binsangltemp = gdat.binsangl
        
    # log-prior
    gdatmodi.thislpri = zeros(gdat.numblpri)

    # allocate memory for variables to hold the proposed state
    ## sample vector
    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
        
    gdatmodi.nextmodlflux = empty_like(gdat.datacnts)
    gdatmodi.nextmodlcnts = empty_like(gdat.datacnts)

    ## likelihood
    gdatmodi.nextllik = zeros_like(gdat.datacnts)
    gdatmodi.deltllik = 0.
        
    ## modification catalog
    gdatmodi.modilgal = zeros(gdat.maxmnumbpntstotl + 2)
    gdatmodi.modibgal = zeros(gdat.maxmnumbpntstotl + 2)
    gdatmodi.modispec = zeros((gdat.numbener, gdat.maxmnumbpntstotl + 2))
    gdatmodi.modispep = zeros((gdat.maxmnumbpntstotl + 2, 2))
    
    # log-prior
    gdatmodi.nextlpri = empty_like(gdatmodi.thislpri)

    # log the initial state
    if gdat.verbtype > 1 and gdat.numbtrap > 0:
        show_samp(gdat, gdatmodi)
        tdpy.util.show_memo(gdatmodi, 'gdatmodi')

    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    gdatmodi.listindxpntsfull = []
    gdatmodi.listsamp = zeros((gdat.numbsamp, gdat.numbpara)) + -1.
    gdatmodi.listsampvarb = zeros((gdat.numbsamp, gdat.numbpara)) - 1
   
    # secondary variables
    gdatmodi.listmodlcnts = zeros((gdat.numbsamp, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    if gdat.pntstype == 'lght':
        gdatmodi.listfwhm = empty((gdat.numbsamp, gdat.numbener, gdat.numbevtt))
        if gdat.varioaxi:
            gdatmodi.listpsfn = zeros((gdat.numbsamp, gdat.numbener, gdat.numbangl + 1, gdat.numbevtt, gdat.numboaxi + 1))
            gdatmodi.listfactoaxi = empty((gdat.numbsamptotl, gdat.numbener, gdat.numbevtt, gdat.numboaxi + 1))
        else:
            gdatmodi.listpsfn = zeros((gdat.numbsamp, gdat.numbener, gdat.numbangl + 1, gdat.numbevtt))
    if gdat.pntstype == 'lens':
        gdatmodi.listdefl = zeros((gdat.numbsamp, gdat.numbsidecart, gdat.numbsidecart, 2))
    
    ## tertiary variables
    gdatmodi.listresicnts = zeros((gdat.numbsamp, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    if gdat.pntstype == 'lght':
        gdatmodi.listpntsflux = zeros((gdat.numbsamp, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        gdatmodi.listpntsfluxmean = zeros((gdat.numbsamp, gdat.numbener))
    if gdat.pntstype == 'lens':
        gdatmodi.listhistdefl = zeros((gdat.numbsamp, gdat.numbdefl))
        gdatmodi.listhostcntsmaps = zeros((gdat.numbsamp, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
        gdatmodi.listlenscnts = zeros((gdat.numbsamp, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
        if gdat.trueinfo:
            gdatmodi.listdeflresi = zeros((gdat.numbsamp, gdat.numbsidecart, gdat.numbsidecart, 2))
            gdatmodi.listdeflcomp = zeros((gdat.numbsamp, gdat.numbsidecart, gdat.numbsidecart))
        gdatmodi.listconv = empty((gdat.numbsamptotl, gdat.numbsidecart, gdat.numbsidecart))
        gdatmodi.listconvpsec = empty((gdat.numbsamptotl, gdat.numbsidewvec, gdat.numbsidewvec))
        gdatmodi.listconvpsecodim = empty((gdat.numbsamptotl, gdat.numbwvecodim))

    gdatmodi.listfluxhistprio = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbinsplot))
    gdatmodi.listsindhistprio = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbinsplot))

    gdatmodi.listllik = zeros((gdat.numbsamp, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    gdatmodi.listlpri = zeros((gdat.numbsamp, gdat.numblpri))
    gdatmodi.listdeltllik = zeros(gdat.numbswep)
    gdatmodi.listdeltlpri = zeros(gdat.numbswep)
    
    gdatmodi.listchrototl = zeros((gdat.numbswep, gdat.numbchrototl))
    gdatmodi.listchrollik = zeros((gdat.numbswep, gdat.numbchrollik))
    gdatmodi.listlprinorm = zeros(gdat.numbsamp)
    gdatmodi.listaccp = zeros(gdat.numbswep, dtype=bool)
    gdatmodi.listboolreje = empty(gdat.numbswep, dtype=bool)
    gdatmodi.listindxfixpmodi = zeros(gdat.numbswep, dtype=int)
    gdatmodi.listindxprop = zeros(gdat.numbswep)
    if gdat.calcerrr:
        gdatmodi.listerrrcnts = empty((gdat.numbsamp, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    gdatmodi.listmemoresi = empty(gdat.numbsamp)
    gdatmodi.listauxipara = zeros((gdat.numbswep, gdat.maxmnumbcompcolr))
    gdatmodi.listlaccfact = zeros(gdat.numbswep)
    gdatmodi.listnumbpair = zeros((gdat.numbswep, 2))
    gdatmodi.listjcbnfact = zeros(gdat.numbswep)
    gdatmodi.listcombfact = zeros(gdat.numbswep)

    ## saved state of the sample index used for logging progress status
    gdatmodi.cntrswepsave = -1
   
    # store the initial sample as the best fit sample
    gdatmodi.maxmllikswep = sum(gdatmodi.thisllik)
    gdatmodi.maxmlposswep = sum(gdatmodi.thisllik) + sum(gdatmodi.thislpri)
    gdatmodi.indxswepmaxmllik = -1 
    gdatmodi.indxswepmaxmlpos = -1 
    gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
    gdatmodi.sampvarbmaxmlpos = copy(gdatmodi.thissampvarb)
   
    # proposal scale optimization
    gdatmodi.stdvstdp = copy(gdat.stdvstdp)
    if gdat.optiprop:
        
        # perform a maximum likelihood fit to bring the sampler to a high likelihood region 
        if False:
            gdat.probtrantemp = gdat.probtran
            gdat.probtran = 0.
            listgdatmodi = [work(gdat, 0)]
            gdat.probtran = gdat.probtrantemp
            gdat.optiproptemp = False
        
        pathstdvprop = gdat.pathopti + '%s.fits' % gdat.rtag
        if os.path.isfile(pathstdvprop) and gdat.loadvaripara:
            if gdat.verbtype > 0 and indxprocwork == 0:
                print 'Reading the previously computed proposal scale from %s...' % pathstdvprop
            gdatmodi.optidone = True
            varipara = pf.getdata(pathstdvprop)
        else:
            if gdat.verbtype > 0 and indxprocwork == 0:
                print 'Optimizing proposal scale...'
            targpropeffi = 0.25
            gdat.factpropeffi = 4.
            minmpropeffi = targpropeffi / gdat.factpropeffi
            maxmpropeffi = targpropeffi * gdat.factpropeffi
            perdpropeffi = 20
            cntrprop = zeros(gdat.numbstdp)
            cntrproptotl = zeros(gdat.numbstdp)
            gdatmodi.optidone = False
            gdatmodi.cntrswepopti = 0
            gdatmodi.cntrswepoptistep = 0
    
        gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
        gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
        gdatmodi.nextindxsamplgal, gdatmodi.nextindxsampbgal, gdatmodi.nextindxsampflux, gdatmodi.nextindxsampspec, gdatmodi.nextindxsampsind, \
                    gdatmodi.nextindxsampcurv, gdatmodi.nextindxsampexpo, gdatmodi.nextindxsampcompcolr = retr_indx(gdat, gdatmodi.nextindxpntsfull, gdat.spectype)
        lposcntr = retr_negalpos(gdat, gdatmodi)
        
        gdatmodi.stdvstdp[gdat.indxstdpcomp] = 0.
        
        gdatmodi.indxstdppara = zeros(gdat.numbpara, dtype=int) - 1
        cntr = 0
        gdat.indxsampactvprop = zeros(gdat.numbfixpactvprop, dtype=int)
        for k in gdat.indxpara:
            if k in gdat.indxfixpactvprop:
                gdatmodi.indxstdppara[k] = cntr
                gdat.indxsampactvprop[cntr] = k
                cntr += 1
            if k in concatenate(gdatmodi.thisindxsamplgal):
                gdatmodi.indxstdppara[k] = gdat.indxstdplgal
            if k in concatenate(gdatmodi.thisindxsampbgal):
                gdatmodi.indxstdppara[k] = gdat.indxstdpbgal
            if k in concatenate(gdatmodi.thisindxsampflux):
                gdatmodi.indxstdppara[k] = gdat.indxstdpflux
            
            if gdat.numbener > 1:
                if k in concatenate(gdatmodi.thisindxsampsind):
                    gdatmodi.indxstdppara[k] = gdat.indxstdpsind
                if k in concatenate(gdatmodi.thisindxsampcurv):
                    gdatmodi.indxstdppara[k] = gdat.indxstdpcurv
                if k in concatenate(gdatmodi.thisindxsampexpo):
                    gdatmodi.indxstdppara[k] = gdat.indxstdpexpo
       
        # temp
        deltparastep = 1e-3
        maxmstdv = 1e-2
        fudgstdv = 1.
        diffpara = deltparastep * array([-1., 0., 1])
        lposdelt = zeros(3)
        gdatmodi.stdvlgaltemp = []
        
        lliktemp = empty(gdat.numbstdp)
        numbiter = diffpara.size
        indxcntr = (numbiter - 1) / 2
        for k in gdat.indxpara:
            if k in gdat.indxfixpactvprop or k in concatenate(gdatmodi.thisindxsampcompcolr):
                for n in range(numbiter):
                    if n == indxcntr:
                        lposdelt[n] = lposcntr
                    else:
                        lposdelt[n] = pert_llik(gdat, gdatmodi, array([k]), array([diffpara[n]]))

                stdv = deltparastep * sqrt(0.5 / amax(fabs(lposcntr - lposdelt))) / sqrt(gdat.numbpara) / fudgstdv
                
                if k in concatenate(gdatmodi.thisindxsampcompcolr):
                    if k in concatenate(gdatmodi.thisindxsamplgal):
                        gdatmodi.stdvlgaltemp.append(stdv)
                        gdatmodi.stdvstdp[gdatmodi.indxstdppara[k]] += stdv * (gdatmodi.thissampvarb[k+2] / gdat.minmflux)
                    else:
                        gdatmodi.stdvstdp[gdatmodi.indxstdppara[k]] += stdv
                else:
                    gdatmodi.stdvstdp[gdatmodi.indxstdppara[k]] = stdv
    
                if False:
                    print 'lposcntr - lposdelt'
                    print lposcntr - lposdelt
                    print 'stdv'
                    print stdv

        gdatmodi.stdvlgaltemp = array(gdatmodi.stdvlgaltemp)
   
        if False:
            print vstack((gdat.namestdp, gdatmodi.stdvstdp)).T
            print
            print 'gdatmodi.stdvlgaltemp'
            print gdatmodi.stdvlgaltemp
            print 'gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampflux)]'
            print gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampflux)]
            print
            raise

        for k in gdat.indxstdp:
            if k in gdat.indxstdpcomp:
                gdatmodi.stdvstdp[k] /= sum(gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts])
            if gdatmodi.stdvstdp[k] > maxmstdv or not isfinite(gdatmodi.stdvstdp[k]):
                gdatmodi.stdvstdp[k] = maxmstdv
   
        gdatmodi.optidone = True
   
        if gdat.makeplot:
            try:
                xdat = gdat.indxstdp
                ydat = gdatmodi.stdvstdp
                path = gdat.pathinit + 'stdv%d.pdf' % indxprocwork
                tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', hist=True)
            
                xdat = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampflux)]
                ydat = gdatmodi.stdvlgaltemp
                path = gdat.pathinit + 'stdvlgalflux.pdf'
                tdpy.util.plot_gene(path, xdat * gdat.anglfact, ydat, scalxdat='logt', lablxdat='$%s$' % gdat.strgflux, lablydat=r'$\sigma_l$ [%s]' % gdat.strganglunit, scat=True)
            except:
                pass
    else:
        gdatmodi.optidone = True
        if gdat.verbtype > 0 and indxprocwork == 0:
            print 'Skipping proposal scale optimization...'

    if gdat.verbtype > 1:
        print 'gdatmodi.stdvstdp'
        print gdatmodi.stdvstdp
        print gdatmodi.stdvstdp[:, None]
        print

    while gdatmodi.cntrswep < gdat.numbswep:
        
        if gdat.emptsamp:
            continue

        timetotlinit = gdat.functime()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                                                                                                        and gdat.makeplot and gdatmodi.optidone
        gdatmodi.boolreje = False
    
        # choose a proposal type
        retr_thisindxprop(gdat, gdatmodi)
            
        if gdat.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            if gdatmodi.cntrswep > 0 and gdat.diagmode:
                print 'gdatmodi.thislpostotlprev'
                print gdatmodi.thislpostotlprev
            print 'thislliktotl'
            print sum(gdatmodi.thisllik)
            print 'thislpritotl'
            print sum(gdatmodi.thislpri)
            print 'thislpostotl'
            print sum(gdatmodi.thisllik) + sum(gdatmodi.thislpri)
            print

        # propose the next sample
        timeinit = gdat.functime()
        
        retr_prop(gdat, gdatmodi)
        
        timefinl = gdat.functime()
        gdatmodi.listchrototl[gdatmodi.cntrswep, 1] = timefinl - timeinit

        # diagnostics
        if gdat.diagmode:
            # sanity checks
            indxsampbadd = where((gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.) | (gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.))[0] + 1
            if indxsampbadd.size > 0:
                print 'cntrswep'
                print gdatmodi.cntrswep
                print 'thisindxprop'
                print gdat.strgprop[gdatmodi.thisindxprop]
                print 'indxsampbadd'
                print indxsampbadd
                print 'drmcsamp'
                print gdatmodi.drmcsamp[indxsampbadd, :]
                raise Exception('Unit sample vector went outside [0,1].')
            
            if gdat.pntstype == 'lght':
                if gdatmodi.propllik:
                    if gdatmodi.indxenermodi.ndim != 1:
                        raise Exception('indxenermodi is inappropriate.')
            if not (gdatmodi.propmerg and gdatmodi.boolreje):
                if not isinstance(gdatmodi.indxsampmodi, ndarray):
                    raise Exception('indxsampmodi is inappropriate.')
           
            gdatmodi.thislpostotl = sum(gdatmodi.thisllik) + sum(gdatmodi.thislpri)
            if False and gdatmodi.thislpostotl - gdatmodi.thislpostotlprev < -30.:
                print 'gdatmodi.thislpostotl'
                print gdatmodi.thislpostotl
                print 'gdatmodi.thislpostotlprev'
                print gdatmodi.thislpostotlprev
                raise Exception('loglikelihood drop is very unlikely!')
            gdatmodi.thislpostotlprev = gdatmodi.thislpostotl
       
            # temp
            if False and gdat.pntstype == 'lght':
                
                if amin(gdatmodi.thispntsflux) < 0.:
                    raise Exception('thispntsflux went negative.')
                
                if amin(gdatmodi.nextpntsflux) < 0.:
                    raise Exception('nextpntsflux went negative.')

            # check what has been changed
            if gdatmodi.cntrswep != 0:
                # temp
                pass
            
            # temp -- only works for single population
            # temp
            #if gdatmodi.cntrswep > 10000 and std(gdatmodi.listsampvarb[where(gdatmodi.listsampvarb[:, gdat.indxfixpnumbpnts[0]] >= 0)[0], gdat.indxfixpnumbpnts[0]]) == 0.:
            #    raise Exception('Number of PS is not changing.')

            # check the population index
            try:
                if gdatmodi.indxpoplmodi < 0:
                    raise
            except:
                pass
            
            if (gdatmodi.thismodlcnts <= 0.).any() or not (isfinite(gdatmodi.thismodlcnts)).all():
                raise Exception('Current flux model is not positive')
            if gdatmodi.cntrswep > 0:
                if not isfinite(gdatmodi.deltlpri):
                    raise Exception('Delta log-prior is not finite.')
                if not isfinite(gdatmodi.deltllik):
                    raise Exception('Delta log-likelihood is not finite.')

            if (gdatmodi.drmcsamp[gdat.indxsampunsd, 1] != 0.).any():
                show_samp(gdat, gdatmodi)
                print 'gdatmodi.drmcsamp[gdat.indxsampunsd, 1]'
                print gdatmodi.drmcsamp[gdat.indxsampunsd, 1]
                raise Exception('Unused vector elements are nonzero.')

            for l in gdat.indxpopl:
                if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]] != len(gdatmodi.thisindxpntsfull[l]):
                    print 'gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]'
                    print gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]
                    print 'gdatmodi.thisindxpntsfull'
                    print gdatmodi.thisindxpntsfull
                    raise Exception('Number of PS is inconsistent with the PS index list.')

                if gdat.numbtrap > 0:
                    if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]] != len(gdatmodi.thisindxpntsfull[l]):
                        raise Exception('Number of PS is inconsistent across data structures.')
                
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]]
                indxtemp = where((flux < gdat.minmflux) | (flux > gdat.maxmflux))[0]
                if indxtemp.size > 0:
                    print
                    print 'Spectrum of a PS went outside the prior range.'
                    print 'l'
                    print l
                    print 'minmflux'
                    print gdat.minmflux
                    print 'maxmflux'
                    print gdat.maxmflux
                    print 'indxtemp'
                    print indxtemp
                    print 'flux'
                    print flux
                    raise Exception('')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            timeinit = gdat.functime()
        
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this')
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]

            # fill the sample lists
            gdatmodi.listsamp[indxsampsave, :] = gdatmodi.drmcsamp[:, 0]
            gdatmodi.listsampvarb[indxsampsave, :] = gdatmodi.thissampvarb
            gdatmodi.listindxpntsfull.append(deepcopy(gdatmodi.thisindxpntsfull))
            
            gdatmodi.listllik[indxsampsave, :] = gdatmodi.thisllik
            gdatmodi.listlpri[indxsampsave, :] = gdatmodi.thislpri
            gdatmodi.listlprinorm[indxsampsave] = gdatmodi.thislprinorm
    
            gdatmodi.listfluxhistprio[indxsampsave, :] = gdatmodi.thisfluxhistprio
            
            if gdat.correxpo:
                gdatmodi.listmodlcnts[indxsampsave, :] = gdatmodi.thismodlcnts
                gdatmodi.listresicnts[indxsampsave, :] = gdatmodi.thisresicnts
                if gdat.pntstype == 'lght':
                    gdatmodi.listpntsflux[indxsampsave, :, :, :] = gdatmodi.thispntsflux
                    gdatmodi.listpntsfluxmean[indxsampsave, :] = gdatmodi.thispntsfluxmean
                if gdat.calcerrr:
                    if gdat.pntstype == 'lght':
                        gdatmodi.listerrrcnts[indxsampsave, :, :] = gdatmodi.thiserrrcnts
            if gdat.evalcirc:
                if gdat.varioaxi:
                    gdatmodi.listfactoaxi[indxsampsave, ...] = gdatmodi.thisfactoaxi
                gdatmodi.listpsfn[indxsampsave, :] = gdatmodi.thispsfn
                gdatmodi.listfwhm[indxsampsave, ...] = gdatmodi.thisfwhm[..., 0]
            
            if gdat.pntstype == 'lens':
                gdatmodi.listdefl[indxsampsave, :, :] = gdatmodi.thisdefl
                gdatmodi.listdeflresi[indxsampsave, :, :, :] = gdatmodi.thisdeflresi
                gdatmodi.listdeflcomp[indxsampsave, :, :] = gdatmodi.thisdeflcomp
                gdatmodi.listhistdefl[indxsampsave, :] = gdatmodi.thishistdefl
                gdatmodi.listhostcntsmaps[indxsampsave, :, :, :] = gdatmodi.thishostcntsmaps
                gdatmodi.listlenscnts[indxsampsave, :, :, :] = gdatmodi.thislenscnts
                gdatmodi.listconv[indxsampsave, :, :] = gdatmodi.thisconv
                gdatmodi.listconvpsec[indxsampsave, :, :] = gdatmodi.thisconvpsec
                gdatmodi.listconvpsecodim[indxsampsave, :] = gdatmodi.thisconvpsecodim
            
            gdatmodi.listmemoresi[indxsampsave] = tdpy.util.retr_memoresi()[0]
            
            timefinl = gdat.functime()
            gdatmodi.listchrototl[gdatmodi.cntrswep, 2] = timefinl - timeinit

        # plot the current sample
        if thismakefram:
            if gdat.verbtype > 0:
                print 'Process %d is in queue for making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                gdat.lock.acquire()
            if gdat.verbtype > 0:
                print 'Process %d started making a frame.' % indxprocwork
            
            timeinit = gdat.functime()
            
            proc_samp(gdat, gdatmodi, 'this')
            plot_samp(gdat, gdatmodi, 'this')
            if gdat.verbtype > 0:
                print 'Process %d finished making a frame.' % indxprocwork
        
            timefinl = gdat.functime()
            gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timeinit
    
            if gdat.numbproc > 1:
                gdat.lock.release()
            
        # reject the sample if proposal is outside the prior
        if where((gdatmodi.drmcsamp[gdatmodi.indxsampchec, 1] < 0.) | (gdatmodi.drmcsamp[gdatmodi.indxsampchec, 1] > 1.))[0].size > 0:
            gdatmodi.boolreje = True

        if gdat.verbtype > 1:
            print 'gdatmodi.boolreje'
            print gdatmodi.boolreje
            show_samp(gdat, gdatmodi)

        if gdat.pntstype == 'lght':
            if gdat.psfntype == 'doubking':
                if gdatmodi.nextsampvarb[gdat.indxfixppsfp[1]] >= gdatmodi.nextsampvarb[gdat.indxfixppsfp[3]]:
                    print 'Rejecting because of the PSF'
                    gdatmodi.boolreje = True
            elif gdat.psfntype == 'doubgaus':
                if gdatmodi.nextsampvarb[gdat.indxfixppsfp[1]] >= gdatmodi.nextsampvarb[gdat.indxfixppsfp[2]]:
                    gdatmodi.boolreje = True
                    print 'Rejecting because of the PSF'
            
        if not gdatmodi.boolreje:
            
            # evaluate the log-prior
            timeinit = gdat.functime()
            
            proc_samp(gdat, gdatmodi, 'next')
            
            timefinl = gdat.functime()
            gdatmodi.listchrototl[gdatmodi.cntrswep, 3] = timefinl - timeinit
   
            if gdat.diagmode:
                if not isfinite(gdatmodi.deltllik):
                    raise Exception('deltlpri is not finite.')
            if gdat.diagmode:
                if not isfinite(gdatmodi.deltllik):
                    raise Exception('deltllik is not finite.')
            
            # evaluate the acceptance probability
            accpprob = exp(gdatmodi.deltllik + gdatmodi.deltlpri + gdatmodi.laccfact)
            
            if gdat.verbtype > 1:
                print 'deltlpri'
                print gdatmodi.deltlpri
                print 'deltllik'
                print gdatmodi.deltllik
                print 'laccfact'
                print gdatmodi.laccfact
                print
        else:
            accpprob = 0.
            gdatmodi.deltllik = 0.
            gdatmodi.deltlpri = 0.
    
        # accept the sample
        if accpprob >= rand():

            if gdat.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(gdat, gdatmodi)

            # check if the accepted sample has
            ## maximal likelihood
            llikswep = gdatmodi.deltllik + sum(gdatmodi.thisllik)
            if llikswep > gdatmodi.maxmllikswep:
                gdatmodi.maxmllikswep = llikswep
                gdatmodi.indxswepmaxmllik = gdatmodi.cntrswep
                gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
            ## maximal posterior
            lposswep = llikswep + gdatmodi.deltlpri + sum(gdatmodi.thislpri)
            if llikswep > gdatmodi.maxmlposswep:
                gdatmodi.maxmlposswep = lposswep
                gdatmodi.indxswepmaxmlpos = gdatmodi.cntrswep
                gdatmodi.sampvarbmaxmlpos = copy(gdatmodi.thissampvarb)
            
            # register the sample as accepted
            gdatmodi.listaccp[gdatmodi.cntrswep] = True

        # reject the sample
        else:

            if gdat.verbtype > 1:
                print 'Rejected.'

            gdatmodi.listaccp[gdatmodi.cntrswep] = False
        
        ## proposal type
        gdatmodi.listindxprop[gdatmodi.cntrswep] = gdatmodi.thisindxprop
        
        ## others
        gdatmodi.listboolreje[gdatmodi.cntrswep] = gdatmodi.boolreje
        gdatmodi.listdeltllik[gdatmodi.cntrswep] = gdatmodi.deltllik
        gdatmodi.listdeltlpri[gdatmodi.cntrswep] = gdatmodi.deltlpri
        if gdatmodi.propsplt or gdatmodi.propmerg:
            gdatmodi.listauxipara[gdatmodi.cntrswep, :gdat.numbcompcolr[gdatmodi.indxpoplmodi], gdatmodi.indxpoplmodi] = gdatmodi.auxipara
            if gdatmodi.thisnumbpair != 0:
                gdatmodi.listnumbpair[gdatmodi.cntrswep, 0] = gdatmodi.thisnumbpair
                if not gdatmodi.boolreje:
                    gdatmodi.listnumbpair[gdatmodi.cntrswep, 1] = gdatmodi.nextnumbpair
                    gdatmodi.listcombfact[gdatmodi.cntrswep] = gdatmodi.combfact
                    gdatmodi.listjcbnfact[gdatmodi.cntrswep] = gdatmodi.jcbnfact
                    gdatmodi.listlaccfact[gdatmodi.cntrswep] = gdatmodi.laccfact
        
        # save the execution time for the sweep
        if not thismakefram:
            timetotlfinl = gdat.functime()
            gdatmodi.listchrototl[gdatmodi.cntrswep, 0] = timetotlfinl - timetotlinit

        # log the progress
        if gdat.verbtype > 0:
            minm = max(0, gdatmodi.cntrswep - 50)
            maxm = gdatmodi.cntrswep + 1
            accp = 100. * where(gdatmodi.listaccp[minm:maxm])[0].size / float(maxm - minm)
            if maxm - minm < 2:
                accp = None
            gdatmodi.cntrswepsave = tdpy.util.show_prog(gdatmodi.cntrswep, gdat.numbswep, gdatmodi.cntrswepsave, indxprocwork=indxprocwork, showmemo=True, accp=accp)

        if gdat.verbtype > 1:
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
            print
        
        # update the sweep counter
        if gdatmodi.optidone:
            gdatmodi.cntrswep += 1
        else:
        
            set_printoptions(precision=2)
            if False:
                print 'perdpropeffi'
                print perdpropeffi
                print 'cntrproptotl'
                print cntrproptotl
                print 'cntrprop'
                print cntrprop
                print 'cntrprop / cntrproptotl'
                print cntrprop / cntrproptotl
                print 'perdpropeffi'
                print perdpropeffi
                print 'gdatmodi.cntrswepopti'
                print gdatmodi.cntrswepopti
                print 'gdatmodi.thisindxprop'
                print gdatmodi.thisindxprop
                print 'gdat.indxfixpconvprop'
                print gdat.indxfixpconvprop
                print 'gdat.numbfixp'
                print gdat.numbfixp
                print 'gdat.strgstdp'
                print gdat.strgstdp

            cntrproptotl[0] += 1.
            if gdatmodi.listaccp[gdatmodi.cntrswep]:
                cntrprop[0] += 1.
            
            if gdatmodi.cntrswepopti % perdpropeffi == 0 and (cntrproptotl[gdat.indxfixpactvprop] > 0).all():
            
                thispropeffi = cntrprop / cntrproptotl
                if gdat.verbtype > 0:
                    print 'Proposal scale optimization step %d' % gdatmodi.cntrswepoptistep
                    print 'Current proposal efficiency'
                    print thispropeffi[gdat.indxstdpactv]
                    print '%.3g +- %.3g' % (mean(thispropeffi[gdat.indxstdpactv]), std(thispropeffi[gdat.indxstdpactv])) 
                if (thispropeffi > minmpropeffi).all() and (thispropeffi < maxmpropeffi).all():
                    if gdat.verbtype > 0:
                        print 'Optimized proposal scale: '
                        print gdat.stdvstdp[gdat.indxstdpactv]
                        print 'Writing the optimized proposal scale to %s...' % pathstdvprop
                    gdatmodi.optidone = True
                    pf.writeto(pathstdvprop, gdat.stdvstdp, clobber=True)

                    plot_opti(gdat, gdatmodi)
                else:
                    factcorr = 2.**(thispropeffi / targpropeffi - 1.)
                    gdat.stdvstdp *= factcorr
                    cntrprop[:] = 0.
                    cntrproptotl[:] = 0.
                    if gdat.verbtype > 0:
                        print 'Correction factor'
                        print factcorr[gdat.indxstdpactv]
                        print 'Current proposal scale: '
                        print gdat.stdvstdp[gdat.indxstdpactv]
                        print
                gdatmodi.cntrswepoptistep += 1
            gdatmodi.cntrswepopti += 1

    gdatmodi.listchrollik = array(gdatmodi.listchrollik)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    return gdatmodi
