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
         diagmode=True, \

         # sampler
         numbswep=None, \
         numbburn=None, \
         factthin=None, \

         seedstat=None, \
         indxevttincl=None, \
         indxenerincl=None, \
        
         evalcirc=None, \

         procrtag=None, \

         # empty run
         emptsamp=False, \

         # comparison with the reference catalog
         anglassc=None, \
         nameexpr=None, \
        
         numbspatdims=2, \
         pntstype='lght', \

         inittype=None, \
         loadvaripara=False, \
         optiprop=None, \
         regulevi=False, \
         strgexprflux=None, \
         strgcatl=None, \
         
         back=None, \
         strgback=None, \
         nameback=None, \
         lablback=None, \
         
         strgexpo=None, \
         
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         exprinfo=None, \
         pixltype=None, \
        
         asscmetrtype='dist', \

         # plotting
         numbswepplot=10000, \
         makeplot=True, \
         makeplotintr=False, \
         scalmaps='asnh', \
         satumaps=None, \
         makeanim=True, \
         anotcatl=False, \
         strgenerfull=None, \
         strgexprname=None, \
         strganglunit=None, \
         strganglunittext=None, \
         anglfact=None, \
         fluxfactplot=None, \
         
         # misc
         lablgangunit=None, \
         labllgal=None, \
         lablbgal=None, \

         # model
         ## PSF
         specfraceval=0.1, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=400, \
         lablfluxunit=None, \
         lablflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxevttfull=None, \
         binsenerfull=None, \
         minmnumbpnts=None, \
         maxmnumbpnts=array([1000]), \
         asymfluxprop=False, \
         psfninfoprio=True, \
         ## spectral

         # prior
         priotype='logt', \
         priofactdoff=0., \
         margfactmodl=1.1, \
         bindprio=False, \
         maxmbacp=None, \
         minmbacp=None, \
         maxmgangdata=None, \
         minmmeanpnts=None, \
         maxmmeanpnts=None, \
    
         spatdistslop=None, \
         
         minmgangdistscal=None, \
         maxmgangdistscal=None, \
         minmbgaldistscal=None, \
         maxmbgaldistscal=None, \
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
         minmcurvdistmean=None, \
         maxmcurvdistmean=None, \
         minmcurvdiststdv=None, \
         maxmcurvdiststdv=None, \
         minmexpodistmean=None, \
         maxmexpodistmean=None, \
         minmexpodiststdv=None, \
         maxmexpodiststdv=None, \

         psfntype=None, \
         oaxitype=None, \
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
    
         curvdistmean=None, \
         curvdiststdv=None, \
         
         minmgang=None, \
         minmflux=None, \
         maxmflux=None, \
        
         # proposals
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
         proppsfp=None, \
         propbacp=True, \
         proplenp=True, \
         propcomp=True, \
         probtran=None, \
         probbrde=1., \
         radispmr=None, \

         exproaxitype=None, \
         exprpsfntype=None, \

         # true data
         truespatdistslop=None, \
         truegangdistscalpop2=None, \
         trueminmflux=None, \
         truemaxmflux=None, \
         truefluxdistslop=None, \
         truefluxdistbrek=None, \
         truefluxdistsloplowr=None, \
         truefluxdistslopuppr=None, \
         truesinddistmean=None, \
         truesinddiststdv=None, \
         truecurvdistmean=None, \
         truecurvdiststdv=None, \
         trueexpodistmean=None, \
         trueexpodiststdv=None, \
         trueoaxitype=None, \
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

         defa=False, \
         **args \
        ):

    
    # construct the global object 
    gdat = tdpy.util.gdatstrt()
    for attr, valu in locals().iteritems():
        if '__' not in attr:
            setattr(gdat, attr, valu)

    if gdat.procrtag != None:
        
        path = gdat.pathdata + gdat.procrtag + '/outp.fits'
        gdat = pf.getdata(path, 1)

        
        pf.writeto(gdat, gdat.stdvstdp, clobber=True)

        
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
            gdat.nameback.append(r'FDM')
    
    if gdat.strgexprflux == None:
        gdat.datatype = 'mock'
    else:
        gdat.datatype = 'inpt'
   
    if gdat.datatype == 'mock':
        if gdat.trueback == None:
            gdat.trueback = gdat.back
        gdat.truenumbback = len(gdat.trueback)
    
    ### number of populations
    gdat.numbpopl = gdat.maxmnumbpnts.size
    gdat.numbback = len(gdat.back)
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)
   
    if gdat.minmnumbpnts == None:
        gdat.minmnumbpnts = zeros(gdat.numbpopl, dtype=int)
    
    if gdat.indxenerincl == None:
        if gdat.exprtype == 'ferm':
            gdat.indxenerincl = arange(1, 4)
        elif gdat.exprtype == 'chan':
            gdat.indxenerincl = arange(2)
        else:
            gdat.indxenerincl = arange(1)
    
    if gdat.indxevttincl == None:
        if gdat.exprtype == 'ferm':
            gdat.indxevttincl = arange(2, 4)
        else:
            gdat.indxevttincl = arange(1)
    
    if gdat.inittype == None:
        if gdat.datatype == 'mock':
            gdat.inittype = 'refr'
        else:
            gdat.inittype = 'blob'

    if gdat.pntstype == 'lens':
        gdat.hubbexpofact = 1.63050e-19
    
    if gdat.strgexpo == None:
        if gdat.pntstype == 'lens':
            gdat.strgexpo = 1e3 / gdat.hubbexpofact
        else:
            gdat.strgexpo = 1.

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

    if gdat.lablflux == None:
        if gdat.pntstype == 'lens':
            gdat.lablflux = r'\theta_E'
        else:
            if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
                gdat.lablflux = 'f'
            if gdat.exprtype == 'sdyn':
                gdat.lablflux = 'p'
    
    if gdat.lablfluxunit == None:
        if gdat.pntstype == 'lens':
            gdat.lablfluxunit = 'arcsec'
        else:
            if gdat.exprtype == 'sdss' or gdat.exprtype == 'hubb':
                gdat.lablfluxunit = 'mag'
            if gdat.exprtype == 'ferm':
                gdat.lablfluxunit = '1/cm$^2$/s/GeV'
            if gdat.exprtype == 'chan':
                gdat.lablfluxunit = '1/cm$^2$/s/KeV'

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
    
    if gdat.proppsfp == None:
        if gdat.pntstype == 'lens':
            gdat.proppsfp = False
        else:
            gdat.proppsfp = True
    
    if gdat.lablgangunit == None:
        if gdat.exprtype == 'ferm':
            gdat.lablgangunit = '$^o$'
        if gdat.exprtype == 'sdyn':
            gdat.lablgangunit = ''
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan' or gdat.exprtype == 'hubb':
            gdat.lablgangunit = '$^{\prime\prime}$'

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

    if gdat.labllgal == None:
        if gdat.exprtype == 'sdyn':
            gdat.labllgal = r'L_z^{\prime}'
        else:
            gdat.labllgal = r'\theta_1'

    if gdat.lablbgal == None:
        if gdat.exprtype == 'sdyn':
            gdat.lablbgal = r'E_k^{\prime}'
        else:
            gdat.lablbgal = r'\theta_2'

    ## experiment defaults
    if gdat.binsenerfull == None:
        if gdat.exprtype == 'ferm':
            gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
        elif gdat.exprtype == 'chan':
            gdat.binsenerfull = array([0.5, 2., 8.])
   
    # energy bin indices
    if gdat.indxenerfull == None:
        if gdat.exprtype == 'ferm':
            gdat.indxenerfull = arange(5)
        elif gdat.exprtype == 'chan':
            gdat.indxenerfull = arange(2)
        else:   
            gdat.indxenerfull = arange(1)
   
    # energy band string
    if gdat.strgenerfull == None:
        if gdat.exprtype == 'sdss':
            gdat.strgenerfull = ['z-band', 'i-band', 'r-band', 'g-band', 'u-band']
        elif gdat.exprtype == 'hubb':
            gdat.strgenerfull = ['F606W']
        elif gdat.exprtype == 'ferm' or gdat.exprtype == 'chan': 
            gdat.strgenerfull = []
            for i in gdat.indxenerfull:
                gdat.strgenerfull.append('%.3g %s - %.3g %s' % (gdat.binsenerfull[i], gdat.strgenerunit, gdat.binsenerfull[i+1], gdat.strgenerunit))
    
    if gdat.evalcirc == None:
        if gdat.pntstype == 'lens':
            gdat.evalcirc = 'bein'
        else:
            gdat.evalcirc = 'psfn'

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
    gdat.pathdataopti = gdat.pathdata + 'opti/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    

    if gdat.exprtype == 'sdyn':
        gdat.enerbins = False
    else:
        gdat.enerbins = True
    
    if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
        gdat.enerdiff = True
    if gdat.exprtype == 'hubb':
        gdat.enerdiff = False
    
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
    gdat.numbenerplot = 100
    gdat.numbener = gdat.indxenerincl.size
    gdat.numbenerfull = len(gdat.strgenerfull)
    gdat.strgener = [gdat.strgenerfull[k] for k in gdat.indxenerincl]
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    if gdat.enerdiff:
        gdat.binsener = gdat.binsenerfull[gdat.indxenerinclbins]
        gdat.minmener = gdat.binsener[0]
        gdat.maxmener = gdat.binsener[-1]
        for strg in ['', 'plot']:
            if strg == '':
                numbbins = gdat.numbener
            else:
                numbbins = gdat.numbenerplot
            retr_axis(gdat, 'ener' + strg, gdat.minmener, gdat.maxmener, numbbins)
        
        gdat.indxenerfull = gdat.binsenerfull.size - 1
    gdat.indxener = arange(gdat.numbener, dtype=int)
    gdat.indxenerfluxdist = ceil(array([gdat.numbener]) / 2.).astype(int) - 1
    # temp
    #gdat.indxenerfluxdist = array([0])
       
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
        gdat.exproaxitype = False
        gdat.exprpsfntype = 'singgaus'
        gdat.exprpsfp = array([0.1 / gdat.anglfact])
 
    # model
    ## hyperparameters
    setp_varbfull(gdat, 'meanpnts', [1., 1e3], numbpopl=gdat.numbpopl)
    
    ### spatial
    if gdat.maxmgangdata == None:
        if gdat.exprtype == 'chan':
            gdat.maxmgangdata = 0.492 / gdat.anglfact * gdat.numbsidecart / 2.
        if gdat.exprtype == 'ferm':
            gdat.maxmgangdata = 20. / gdat.anglfact
        if gdat.exprtype == 'sdyn':
            gdat.maxmgangdata = 1.
        if gdat.exprtype == 'hubb':
            gdat.maxmgangdata = 2. / gdat.anglfact
    
    for strg, valu in args.iteritems():
        setattr(gdat, strg, valu)

    gdat.liststrgtype = ['spatdisttype', 'fluxdisttype', 'spectype']
    for strg in gdat.liststrgtype:
   
        # default the true model types
        try:
            truetype = getattr(gdat, 'true' + strg)
        except:
            if strg == 'spatdisttype':
                strgtemp = 'unif'
            if strg == 'fluxdisttype' or strg == 'spectype':
                strgtemp = 'powr'

            truetype = [strgtemp for l in gdat.indxpopl]
            setattr(gdat, 'true' + strg, truetype)
        
        # default the model types to the true model
        try:
            getattr(gdat, strg)
        except:
            setattr(gdat, strg, truetype)
    
    if gdat.minmgang == None:
        if gdat.exprtype == 'ferm':
            gdat.minmgang = 1e-1 / gdat.anglfact
        else:
            gdat.minmgang = 1e-2 / gdat.anglfact
    
    if gdat.minmflux == None:
        if gdat.exprtype == 'ferm':
            gdat.minmflux = 3e-11
        if gdat.exprtype == 'chan':
            gdat.minmflux = 1e-7
        if gdat.exprtype == 'sdyn':
            gdat.minmflux = 1e0
        if gdat.pntstype == 'lens':
            gdat.minmflux = 5e-3 / gdat.anglfact
    
    if gdat.maxmflux == None:
        if gdat.exprtype == 'ferm':
            gdat.maxmflux = 1e-7
        if gdat.exprtype == 'chan':
            gdat.maxmflux = 1e-3
        if gdat.exprtype == 'sdyn':
            gdat.maxmflux = 1e4
        if gdat.pntstype == 'lens':
            gdat.maxmflux = 0.1 / gdat.anglfact
   
    # PS spectral model
    if gdat.spectype == None:
        if gdat.truespectype != None:
            gdat.spectype = gdat.truespectype
        else:
            gdat.spectype = array(['powr' for l in gdat.indxpopl])

    ## PSF
    if gdat.psfntype == None:
        gdat.psfntype = gdat.exprpsfntype

    if gdat.oaxitype == None:
        if gdat.exprtype == 'chan':
            gdat.oaxitype = True
        else:
            gdat.oaxitype = False

    # maximum horizontal/vertical distance of the elements from the image center
    gdat.maxmgang = gdat.maxmgangdata * gdat.margfactmodl

    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
    
    # parameter defaults
    ## distribution
    ### flux
    setp_varbfull(gdat, 'gangdistscal', [1. / gdat.anglfact, 10. / gdat.anglfact], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'bgaldistscal', [0.5 / gdat.anglfact, 5. / gdat.anglfact], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistslop', [0.5, 3.5], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistbrek', [gdat.minmflux, gdat.maxmflux], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistsloplowr', [-1.5, 3.5], numbpopl=gdat.numbpopl)
    setp_varbfull(gdat, 'fluxdistslopuppr', [1.5, 3.5], numbpopl=gdat.numbpopl)
    
    ### spectral index
    if gdat.numbener > 1:
        setp_varbfull(gdat, 'sinddistmean', [1., 3.])
        setp_varbfull(gdat, 'sinddiststdv', [0.01, 1.])
        setp_varbfull(gdat, 'curvdistmean', [-1., 1.])
        setp_varbfull(gdat, 'curvdiststdv', [0.1, 1.])
        setp_varbfull(gdat, 'expodistmean', [1., 8.])
        setp_varbfull(gdat, 'expodiststdv', [0.01 * gdat.maxmener, gdat.maxmener])
   
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
    
    ## lensing
    gdat.minmlgalsour = gdat.minmlgal
    gdat.maxmlgalsour = gdat.maxmlgal
    gdat.minmbgalsour = gdat.minmbgal
    gdat.maxmbgalsour = gdat.maxmbgal
    gdat.minmlgalhost = gdat.minmlgal
    gdat.maxmlgalhost = gdat.maxmlgal
    gdat.minmbgalhost = gdat.minmbgal
    gdat.maxmbgalhost = gdat.maxmbgal
    setp_varbfull(gdat, 'specsour', array([1e-21, 1e-17]) )
    setp_varbfull(gdat, 'sizesour', [0.01 / gdat.anglfact, 1. / gdat.anglfact])
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
        gdat.numbburn = gdat.numbswep / 5

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
    if gdat.pixltype == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')
    if gdat.pntstype == 'lens' and (gdat.numbback > 1 or not isinstance(gdat.back[0], float)):
        raise Exception('In a lensing problem, the background can only be uniform.')
    if gdat.minmflux >= gdat.maxmflux:
        raise Exception('Minimum flux is greater than maximum flux.')

    if gdat.verbtype > 0:
        print 'PCAT v0.1 started at %s' % gdat.strgtimestmp
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)

    # true model defaults
    # temp
    if gdat.trueminmflux == None:
        gdat.trueminmflux = gdat.minmflux
    if gdat.truemaxmflux == None:
        gdat.truemaxmflux = gdat.maxmflux
        
    defn_truedefa(gdat, 4. / gdat.anglfact, 'gangdistscalpop2')
    defn_truedefa(gdat, 2. / gdat.anglfact, 'bgaldistscalpop1')
    defn_truedefa(gdat, 2.6, 'fluxdistsloppop0')
    defn_truedefa(gdat, 2.2, 'fluxdistsloppop1')
    defn_truedefa(gdat, 3.2, 'fluxdistsloppop2')
    defn_truedefa(gdat, sqrt(gdat.trueminmflux * gdat.maxmflux), 'fluxdistbrek')
    defn_truedefa(gdat, 1., 'fluxdistsloplowr')
    defn_truedefa(gdat, 2., 'fluxdistslopuppr')
    
    defn_truedefa(gdat, 2.5, 'sinddistmeanpop0')
    defn_truedefa(gdat, 0.2, 'sinddiststdvpop0')
    defn_truedefa(gdat, 1.5, 'sinddistmeanpop1')
    defn_truedefa(gdat, 0.2, 'sinddiststdvpop1')
    defn_truedefa(gdat, 1.5, 'sinddistmeanpop2')
    defn_truedefa(gdat, 0.2, 'sinddiststdvpop2')
    
    defn_truedefa(gdat, 2., 'curvdistmean')
    defn_truedefa(gdat, 0.5, 'curvdiststdv')
    
    defn_truedefa(gdat, 2., 'expodistmeanpop1')
    defn_truedefa(gdat, 2., 'expodistmeanpop2')
    defn_truedefa(gdat, 0.2, 'expodiststdvpop1')
    defn_truedefa(gdat, 0.2, 'expodiststdvpop2')
    
    defn_truedefa(gdat, 1., 'bacp')
    
    defn_truedefa(gdat, 0.05 / gdat.anglfact, 'sizesour')
    
    if gdat.defa:
        return gdat

    # initial setup
    setpinit(gdat, True) 
   
    # create the PCAT folders
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathoutpthis = gdat.pathoutp + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    os.system('mkdir -p %s' % gdat.pathoutpthis)
    pathcatllite = gdat.pathoutpthis + 'catllite.fits'  
    pathcatl = gdat.pathoutpthis + 'catl.fits'  
    
    if gdat.exprtype == 'ferm':
        retr_fermdata(gdat)
    if gdat.exprtype == 'chan':
        retr_chandata(gdat)

    # rotate PS coordinates to the ROI center
    if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
        rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='ZYX')
        gdat.exprbgal, gdat.exprlgal = rttr(pi / 2. - gdat.exprbgal, gdat.exprlgal)
        gdat.exprbgal = pi / 2. - gdat.exprbgal

    # external reference catalog
    if gdat.exprinfo:
        gdat.indxexprpntsrofi = where((fabs(gdat.exprlgal) < gdat.maxmgangdata) & (fabs(gdat.exprbgal) < gdat.maxmgangdata))[0]
        for strgfeat in gdat.liststrgfeat:
            try:
                feat = getattr(gdat, 'expr' + strgfeat)
            except:
                feat = None
            setattr(gdat, 'expr' + strgfeat + 'totl', feat)
            if feat == None:
                setattr(gdat, 'expr' + strgfeat, None)
            else:
                setattr(gdat, 'expr' + strgfeat, feat[..., gdat.indxexprpntsrofi])
        gdat.exprnumbpnts = gdat.exprlgal.size
        
        # reorder PS with respect to flux
        indxpnts = argsort(gdat.exprspec[0, gdat.indxenerfluxdist[0], :])[::-1]
        gdat.exprlgal = gdat.exprlgal[indxpnts]
        gdat.exprbgal = gdat.exprbgal[indxpnts]
        gdat.exprspec[0, :, :] = gdat.exprspec[0, :, indxpnts].T
        if gdat.exprcnts != None:
            gdat.exprcnts = gdat.exprcnts[:, indxpnts, :]

        # compute the catalog counts based on the exposure
        gdat.exprcntscalc = empty((gdat.numbener, gdat.exprnumbpnts, gdat.numbevtt))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                indxpixltemp = retr_indxpixl(gdat, gdat.exprbgal, gdat.exprlgal)
                gdat.exprcntscalc[i, :, m] = gdat.exprspec[0, i, :] * gdat.expo[i, indxpixltemp, m] * gdat.deltener[i]
        
        if gdat.strgcnfg == 'pcat_chan_inpt':
            print 'gdat.exprcnts'
            print gdat.exprcnts[:, :, 0].T
            print 'gdat.exprcntscalc'
            print gdat.exprcntscalc[:, :, 0].T

        if gdat.exprcnts != None and gdat.exprlgal.size > 0 and gdat.verbtype > 0:
            if amax(fabs((gdat.exprcnts - gdat.exprcntscalc) / gdat.exprcnts)) > 0.01:
                print 'Experimental information on PS counts is inconsistent.'
        
        gdat.exprcnts = gdat.exprcntscalc

        gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
        gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)
        
        if not isfinite(gdat.exprspec).all():
            raise Exception('exprspec is not finite.')
        
        if gdat.exprnumbpnts > 0:
            gdat.exprfluxbrgt, gdat.exprfluxbrgtassc = retr_fluxbrgt(gdat, gdat.exprlgal, gdat.exprbgal, gdat.exprspec[0, gdat.indxenerfluxdist[0], :])

    # generate true data
    if gdat.datatype == 'mock':
        
        if gdat.verbtype > 0:
            print 'Generating mock data...'

        if gdat.seedstat != None:
            if gdat.verbtype > 0:
                print 'Setting the seed for the RNG...'
            set_state(gdat.seedstat)

        if gdat.truenumbpnts == None:
            gdat.truenumbpnts = empty(gdat.numbpopl, dtype=int)
            if gdat.numbtrap > 0:
                for l in gdat.indxpopl:
                    gdat.truenumbpnts[l] = random_integers(0, gdat.maxmnumbpnts[l])
    
    else:
        

        gdat.truelgal = [gdat.exprlgal]
        gdat.truebgal = [gdat.exprbgal]
        gdat.truenumbpnts = array([gdat.exprlgal.size])
        gdat.truespec = [gdat.exprspec]
        gdat.truecnts = [gdat.exprcnts]
        gdat.truesind = [gdat.exprsind]
        #gdat.truestrg = [gdat.exprstrg]
        #gdat.truestrgclss = [gdat.exprstrgclss]
        #gdat.truestrgassc = [gdat.exprstrgassc]
        
        gdat.trueminmflux = amin(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
        gdat.truemaxmflux = amax(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
        for l in gdat.indxpopl: 
            gdat.trueminmflux = min(gdat.trueminmflux, amin(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
            gdat.truemaxmflux = max(gdat.truemaxmflux, amax(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
            
    gdat.truenumbpopl = gdat.truenumbpnts.size
    gdat.trueindxpopl = arange(gdat.truenumbpopl, dtype=int)
    
    if gdat.datatype == 'mock':
        
        if gdat.truespectype == None:
            gdat.truespectype = ['powr' for l in gdat.trueindxpopl]
        if gdat.trueoaxitype == None:
            gdat.trueoaxitype = gdat.exproaxitype
        if gdat.truepsfntype == None:
            gdat.truepsfntype = gdat.exprpsfntype

        # set mock sample vector indices
        retr_indxsamp(gdat, strgpara='true')
    
        ## unit sample vector
        gdat.truesamp = zeros(gdat.truenumbpara)
        gdat.truefixp = zeros(gdat.truenumbfixp) + nan
        gdat.truefixp[gdat.trueindxfixpnumbpnts] = gdat.truenumbpnts
   
        gdat.trueindxpntsfull = []
        if gdat.numbtrap > 0:
            for l in gdat.trueindxpopl:
                gdat.trueindxpntsfull.append(range(gdat.truenumbpnts[l]))
        else:
            gdat.trueindxpntsfull = []
        gdat.trueindxsamplgal, gdat.trueindxsampbgal, gdat.trueindxsampflux, gdat.trueindxsampsind, gdat.trueindxsampcurv, gdat.trueindxsampexpo, \
                                                                                      gdat.trueindxsampcomp = retr_indx(gdat, gdat.trueindxpntsfull, gdat.truespectype)
        
        gdat.truefixp[gdat.trueindxfixpmeanpnts] = gdat.truefixp[gdat.trueindxfixpnumbpnts]


        for k in gdat.trueindxfixp:
            
            if k in gdat.trueindxfixpnumbpnts or k in gdat.trueindxfixpmeanpnts:
                continue

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
        gdat.truespecplot = [[] for l in gdat.trueindxpopl]
        if gdat.truenumbtrap > 0:
            gdat.truesind = [empty(gdat.truenumbpnts[l]) for l in gdat.trueindxpopl]
            gdat.truecurv = [empty(gdat.truenumbpnts[l]) for l in gdat.trueindxpopl]
            gdat.trueexpo = [empty(gdat.truenumbpnts[l]) for l in gdat.trueindxpopl]
        
            for l in gdat.trueindxpopl:
                if gdat.truespatdisttype[l] == 'unif':
                    gdat.truelgal[l] = icdf_self(rand(gdat.truenumbpnts[l]), -gdat.maxmgangdata, 2. * gdat.maxmgangdata)
                    gdat.truebgal[l] = icdf_self(rand(gdat.truenumbpnts[l]), -gdat.maxmgangdata, 2. * gdat.maxmgangdata) 

                if gdat.truespatdisttype[l] == 'gang':
                    gdat.truegang[l] = icdf_expo(rand(gdat.truenumbpnts[l]), gdat.maxmgang, gdat.truefixp[gdat.trueindxfixpgangdistscal[l]])
                    gdat.trueaang[l] = icdf_self(rand(gdat.truenumbpnts[l]), 0., 2. * pi)
                    gdat.truelgal[l], gdat.truebgal[l] = retr_lgalbgal(gdat.truegang[l], gdat.trueaang[l])
                
                if gdat.truespatdisttype[l] == 'disc':
                    gdat.truelgal[l] = icdf_self(rand(gdat.truenumbpnts[l]), -gdat.maxmgangdata, 2. * gdat.maxmgangdata)
                    gdat.truebgal[l] = icdf_expo(rand(gdat.truenumbpnts[l]), gdat.maxmgang, gdat.truefixp[gdat.trueindxfixpbgaldistscal[l]]) * \
                                                                                                choice(array([1., -1.]), size=gdat.truenumbpnts[l])
                
                gdat.truespec[l] = empty((3, gdat.numbener, gdat.truenumbpnts[l]))
                gdat.truespecplot[l] = empty((gdat.numbenerplot, gdat.truenumbpnts[l]))
                if gdat.truefluxdisttype[l] == 'powr':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_powr(rand(gdat.truenumbpnts[l]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                                                    gdat.truefixp[gdat.trueindxfixpfluxdistslop[l]])
                if gdat.truefluxdisttype[l] == 'brok':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_brok(rand(gdat.truenumbpnts[l]), \
                                                gdat.trueminmflux, gdat.maxmflux, gdat.truefluxdistbrek[l], gdat.truefluxdistsloplowr[l], gdat.truefluxdistslopuppr[l])
                
                # temp -- make sure this reordering does not mess up other things
                gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = sort(gdat.truespec[l][:, gdat.indxenerfluxdist[0], :], axis=1)[::-1]
    
                if gdat.numbener > 1:
                    # spectral parameters
                    gdat.truesind[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpsinddistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpsinddiststdv[l]])
                    if gdat.truespectype[l] == 'curv':
                        gdat.truecurv[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpcurvdistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpcurvdiststdv[l]])
                    
                    if gdat.truespectype[l] == 'expo':
                        gdat.trueexpo[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpexpodistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpexpodiststdv[l]])
                
                    # spectra
                    gdat.truespec[l][:] = retr_spec(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truesind[l], \
                                                                                                             gdat.truecurv[l], gdat.trueexpo[l], gdat.truespectype[l])[None, :, :]
                    
                    gdat.truespecplot[l] = retr_spec(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truesind[l], \
                                                                                              gdat.truecurv[l], gdat.trueexpo[l], gdat.truespectype[l], plot=True)
                
                gdat.truesampvarb[gdat.trueindxsamplgal[l]] = gdat.truelgal[l]
                gdat.truesampvarb[gdat.trueindxsampbgal[l]] = gdat.truebgal[l]
                gdat.truesampvarb[gdat.trueindxsampflux[l]] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                if gdat.numbener > 1:
                    gdat.truesampvarb[gdat.trueindxsampsind[l]] = gdat.truesind[l]
                    if gdat.truespectype[l] == 'curv':
                        gdat.truesampvarb[gdat.trueindxsampcurv[l]] = gdat.truecurv[l]
                    if gdat.truespectype[l] == 'expo':
                        gdat.truesampvarb[gdat.trueindxsampexpo[l]] = gdat.trueexpo[l]
                if gdat.verbtype > 1:
                    print 'l'
                    print l
                    print 'truebgal[l]'
                    print gdat.truebgal[l] * gdat.anglfact
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
                    if gdat.enerdiff:
                        gdat.truecnts[l] *= gdat.deltener[:, None, None]
        
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
            truepsfnicdfintp = interp1d_pick(gdat.truepsfncdfn, gdat.binsangl)

            if gdat.trueindxpntstotl.size > 0:
                for k in gdat.indxdatasamp:
                    indxpntsthis = choice(gdat.trueindxpntstotl, p=probpntsflux) 
                    radi = truepsfnicdfintp(rand())
                    aang = rand() * 2. * pi
                    gdat.truedatacnts[0, k, 0, 0] = truelgaltemp[indxpntsthis] + radi * sin(aang)
                    gdat.truedatacnts[0, k, 0, 1] = truebgaltemp[indxpntsthis] + radi * cos(aang)
        
        if gdat.seedstat != None:
            seed()
    
    gdat.trueflux = [[] for l in gdat.trueindxpopl]
    for l in gdat.trueindxpopl:
        gdat.trueflux[l] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]

    # final setup
    setpfinl(gdat, True) 

    if gdat.makeplot:
        plot_init(gdat)

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
        if gdat.evalcirc != 'full':
            print 'maxmangleval'
            print gdat.anglfact * gdat.maxmangleval, ' [%s]' % gdat.strganglunit
        print
            
    if gdat.verbtype > 0:
        tdpy.util.show_memo(gdat, 'gdat')

    # lock the global object againts any future modifications
    gdat.lockmodi()

    gdat.timereal = zeros(gdat.numbproc)
    gdat.timeproc = zeros(gdat.numbproc)
    if gdat.numbproc == 1:
        listgdatmodi = [worktrac(gdat, 0)]
    else:
        if gdat.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        gdat.lock = mp.Manager().Lock()
        
        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(worktrac, gdat)
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
    gdat.liststrgchanarry = []
    for strg, valu in listgdatmodi[0].__dict__.iteritems():
        if strg.startswith('list') and isinstance(valu, ndarray):
            gdat.liststrgchanarry.append(strg[4:])
    
    gdat.liststrgchan = gdat.liststrgchanarry + ['specassc']

    for strg in gdat.liststrgchanarry:
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
    for k in gdat.indxproc:
        gdat.maxmllikswep[k] = listgdatmodi[k].maxmllikswep
        gdat.indxswepmaxmllik[k] = listgdatmodi[k].indxswepmaxmllik
        gdat.sampvarbmaxmllik[k] = listgdatmodi[k].sampvarbmaxmllik
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

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
    gdat.listspecassc = []
    for j in gdat.indxsamp:      
        for k in gdat.indxproc:
            gdat.listindxpntsfull.append(listgdatmodi[k].listindxpntsfull[j])
            gdat.listspecassc.append(listgdatmodi[k].listspecassc[j])

    ## list of other parameters to be flattened
    gdat.liststrgchanflat = deepcopy(gdat.liststrgchanarry)
    for strg in ['deltllik', 'memoresi']:
        gdat.liststrgchanflat.remove(strg)
   
    gdat.listsampvarbproc = copy(gdat.listsampvarb)

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
    
    # find the maximum likelihood and posterior over the chains
    gdat.indxprocmaxmllik = argmax(gdat.maxmllikswep)
    gdat.maxmllikswep = gdat.maxmllikswep[gdat.indxprocmaxmllik]
    gdat.indxswepmaxmllik = gdat.indxprocmaxmllik * gdat.numbsamp + gdat.indxswepmaxmllik[gdat.indxprocmaxmllik]
    gdat.sampvarbmaxmllik = gdat.sampvarbmaxmllik[gdat.indxprocmaxmllik, :]
    
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
        gdat.listflux = [[] for l in gdat.indxpopl]
        if gdat.numbener > 1:
            gdat.listsind = [[] for l in gdat.indxpopl]
            gdat.listcurv = [[] for l in gdat.indxpopl]
            gdat.listexpo = [[] for l in gdat.indxpopl]
        for n in gdat.indxsamptotl: 
            for l in gdat.indxpopl:
                indxsamplgal, indxsampbgal, indxsampflux, indxsampsind, indxsampcurv, indxsampexpo, indxsampcomp = retr_indx(gdat, gdat.listindxpntsfull[n], gdat.spectype)
                gdat.listlgal[l].append(gdat.listsampvarb[n, indxsamplgal[l]])
                gdat.listbgal[l].append(gdat.listsampvarb[n, indxsampbgal[l]])
                gdat.listflux[l].append(gdat.listsampvarb[n, indxsampflux[l]])
                if gdat.numbener > 1:
                    gdat.listsind[l].append(gdat.listsampvarb[n, indxsampsind[l]])
                    if gdat.spectype[l] == 'curv':
                        gdat.listcurv[l].append(gdat.listsampvarb[n, indxsampcurv[l]])
                    if gdat.spectype[l] == 'expo':
                        gdat.listexpo[l].append(gdat.listsampvarb[n, indxsampexpo[l]])
        
    ## bin PS parameters
    if gdat.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        timeinit = gdat.functime()
        
    gdat.pntsprob = zeros((gdat.numbpopl, gdat.numbbgalpntsprob, gdat.numblgalpntsprob, gdat.numbbinsplot))
    for l in gdat.indxpopl:
        temparry = concatenate([gdat.listlgal[l][n] for n in gdat.indxsamptotl])
        temp = empty((len(temparry), 3))
        temp[:, 0] = temparry
        temp[:, 1] = concatenate([gdat.listbgal[l][n] for n in gdat.indxsamptotl])
        temp[:, 2] = concatenate([gdat.listflux[l][n] for n in gdat.indxsamptotl])
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
    gdat.listindxsamptotlaccpprop = []
    gdat.listindxsamptotl = []
    gdat.listindxsamptotlaccp = []
    gdat.listindxsamptotlreje = []
    for n in gdat.indxproptype:
        gdat.listindxsamptotl.append(where(gdat.listindxproptype == gdat.indxproptype[n])[0])
        gdat.listindxsamptotlaccp.append(intersect1d(where(gdat.listindxproptype == gdat.indxproptype[n])[0], where(gdat.listaccp)[0]))
        gdat.listindxsamptotlaccpprop.append(intersect1d(where(gdat.listindxproptype == gdat.indxproptype[n])[0], where(gdat.listaccpprop)[0]))
        gdat.listindxsamptotlreje.append(intersect1d(where(gdat.listindxproptype == gdat.indxproptype[n])[0], where(logical_not(gdat.listaccp))[0]))
    
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
    gdat.liststrgchan += ['fixp']
    for strg in gdat.liststrgchan:
        listtemp = getattr(gdat, 'list' + strg)
        if isinstance(listtemp, list):
            posttemp = []
            meditemp = []
            errrtemp = []
            numbelem = len(listtemp[0])
            for k in range(numbelem):
                shap = [gdat.numbsamptotl] + list(listtemp[0][k].shape)
                temp = zeros(shap)
                for n in gdat.indxsamptotl:
                    temp[n, ...] = listtemp[n][k]
                posttempsing = tdpy.util.retr_postvarb(temp)
                meditempsing = posttempsing[0, ...]
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
    if gdat.pntstype == 'lght':
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
    gdat.timeproctotlswep = gdat.timeproctotl / gdat.numbswep
    gdat.timeprocnorm = gdat.timeproctotlswep / gdat.timeatcr

    if gdat.verbtype > 0:
        for k in gdat.indxproc:
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdat.timereal[k], gdat.timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (gdat.timerealtotl, gdat.timeproctotl)
        if gdat.makeplot:
            print 'The plots are in ' + gdat.pathplot
        print 'The ensemble of catalogs is at ' + pathcatl
    return gdat
    

def worktrac(gdat, indxprocwork):
	
    try:
        return work(gdat, indxprocwork)
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def work(gdat, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for this chain
    seed()
    
    # empty object to hold chain-specific variables that will be modified by the chain
    gdatmodi = tdpy.util.gdatstrt()
  
    # construct the initial state
    if gdat.verbtype > 0:
        print 'Initializing the sampler state...'
   
    ## unit sample vector
    gdatmodi.thissamp = zeros(gdat.numbpara)
    gdatmodi.thissampvarb = zeros(gdat.numbpara)
   
    ## Fixed-dimensional parameters
    for k in gdat.indxfixp:
        if gdat.inittype == 'rand':
            if k in gdat.indxfixpnumbpnts:
                gdatmodi.thissamp[gdat.indxfixp[k]] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[k] + 1))
            else:
                gdatmodi.thissamp[gdat.indxfixp[k]] = rand()
        else:
            if gdat.datatype == 'mock':
                gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, gdat.truefixp[k], k)
            else:
                if k == 0:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, gdat.truenumbpnts[0], k)
                else:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = rand()
        
        gdatmodi.thissampvarb[k] = icdf_fixp(gdat, '', gdatmodi.thissamp[k], k)
                
    ## lists of occupied and empty transdimensional parameters
    gdatmodi.thisindxpntsfull = []
    gdatmodi.thisindxpntsempt = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            gdatmodi.thisindxpntsfull.append(range(gdatmodi.thissamp[gdat.indxfixpnumbpnts[l]].astype(int)))
            gdatmodi.thisindxpntsempt.append(range(gdatmodi.thissamp[gdat.indxfixpnumbpnts[l]].astype(int), gdat.maxmnumbpnts[l]))
    else:
        gdatmodi.thisindxpntsfull = []
    
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampflux, gdatmodi.thisindxsampsind, \
                             gdatmodi.thisindxsampcurv, gdatmodi.thisindxsampexpo, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
    
    if gdat.numbtrap > 0:
        
        for l in gdat.indxpopl:
            gdatmodi.thissamp[gdatmodi.thisindxsampcomp[l]] = rand(gdatmodi.thisindxsampcomp[l].size)
        
        ## PS components
        if gdat.inittype == 'refr':
            for l in gdat.indxpopl:
                gdatmodi.thissamp[gdatmodi.thisindxsamplgal[l]] = cdfn_self(gdat.truelgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                gdatmodi.thissamp[gdatmodi.thisindxsampbgal[l]] = cdfn_self(gdat.truebgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                
                indxtruepntsgood = where(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :] > gdat.minmflux)[0]
                if gdat.fluxdisttype[l] == 'powr':
                    fluxunit = cdfn_powr(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.minmflux, gdat.maxmflux, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]])
                if gdat.fluxdisttype[l] == 'brok':
                    flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                    fluxunit = cdfn_flux_brok(flux, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
                gdatmodi.thissamp[gdatmodi.thisindxsampflux[l][indxtruepntsgood]] = fluxunit[indxtruepntsgood]

                if gdat.numbener > 1:
                    gdatmodi.thissamp[gdatmodi.thisindxsampsind[l]] = cdfn_gaus(gdat.truesind[l], gdatmodi.thissampvarb[gdat.indxfixpsinddistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpsinddiststdv[l]])
                    if gdat.spectype[l] == 'curv':
                        try:
                            if gdat.verbtype > 0:
                                print 'Initializing the spectral curvatures randomly for population %d' % l
                            gdatmodi.thissamp[gdatmodi.thisindxsampcurv[l]] = cdfn_gaus(gdat.truecurv[l], gdatmodi.thissampvarb[gdat.indxfixpcurvdistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpcurvdiststdv[l]])
                        except:
                            gdatmodi.thissamp[gdatmodi.thisindxsampcurv[l]] = rand(gdat.truenumbpnts[l])
                            
                    if gdat.spectype[l] == 'expo':
                        try:
                            gdatmodi.thissamp[gdatmodi.thisindxsampexpo[l]] = cdfn_gaus(gdat.trueexpo[l], gdatmodi.thissampvarb[gdat.indxfixpexpodistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpexpodiststdv[l]])
                        except:
                            if gdat.verbtype > 0:
                                print 'Initializing the spectral cutoff energies randomly for population %d' % l
                            gdatmodi.thissamp[gdatmodi.thisindxsampexpo[l]] = rand(gdat.truenumbpnts[l])
            
            randinittemp = False
       
    if gdat.verbtype > 1:
        print 'thissamp'
        for k in gdat.indxpara:
            if k in concatenate(gdatmodi.thisindxsamplgal):
                print
            print '%15f' % gdatmodi.thissamp[k]

    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.thissamp[gdat.numbpopl:] < 0.)[0] + gdat.numbpopl
    indxsampbadduppr = where(gdatmodi.thissamp[gdat.numbpopl:] > 1.)[0] + gdat.numbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'indxsampbadd'
        print indxsampbadd
        print gdat.namepara[indxsampbadd]
        print 'gdatmodi.thissamp'
        print gdatmodi.thissamp[:, None]
        raise Exception('Initial unit sample vector went outside the unit interval...')
        #gdatmodi.thissamp[indxsampbaddlowr] = 0.
        #gdatmodi.thissamp[indxsampbadduppr] = 1.
        #print 'Initial unit sample vector went outside the unit interval...'
    
    ## sample vector
    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.thissamp, 'this')
    
    if gdat.verbtype > 1:
        print 'thissamp, thissampvarb'
        for k in gdat.indxpara:
            if k in concatenate(gdatmodi.thisindxsamplgal):
                print
            print '%15f %15f' % (gdatmodi.thissamp[k], gdatmodi.thissampvarb[k])
    
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
        
    # allocate memory for variables to hold the proposed state
    ## sample vector
    gdatmodi.nextsamp = copy(gdatmodi.thissamp)
    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
        
    gdatmodi.nextmodlflux = empty_like(gdat.datacnts)
    gdatmodi.nextmodlcnts = empty_like(gdat.datacnts)

    ## likelihood
    gdatmodi.nextllik = zeros_like(gdat.datacnts)
        
    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    gdatmodi.listindxpntsfull = []
    gdatmodi.listspecassc = []
    
    gdatmodi.liststrgvarbswep = ['memoresi', 'lpri', 'lfctprop', 'lpriprop', 'lpau', 'deltllik', 'chrototl', 'chrollik', 'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype']
    if gdat.probbrde < 1.:
        gdatmodi.liststrgvarbswep += ['auxipara', 'numbpair', 'jcbnfact', 'combfact']
   
    gdatmodi.liststrgvarbsamp = []
    for strg, valu in gdatmodi.__dict__.iteritems():
        if strg.startswith('this') and isinstance(valu, ndarray) and not strg[4:] in gdatmodi.liststrgvarbswep:
            gdatmodi.liststrgvarbsamp.append(strg[4:])
    
    # temp
    gdatmodi.thismemoresi = zeros(1)
    gdatmodi.thisdeltllik = zeros(1)
    gdatmodi.thischrollik = zeros(gdat.numbchrollik)
    gdatmodi.thischrototl = zeros(gdat.numbchrototl)
    gdatmodi.thisaccp = zeros(1, dtype=bool)
    gdatmodi.thisaccppsfn = zeros(1, dtype=bool)
    gdatmodi.thisaccpprio = zeros(1, dtype=bool)
    gdatmodi.thisaccpprop = zeros(1, dtype=bool)
    gdatmodi.thisindxproptype = zeros(1, dtype=int)
    gdatmodi.thisauxipara = zeros(gdat.maxmnumbcomp)
    gdatmodi.thisnumbpair = zeros(1, dtype=int)
    gdatmodi.thisjcbnfact = zeros(1)
    gdatmodi.thiscombfact = zeros(1)
    gdatmodi.thislpau = zeros(gdat.numblpau)
    gdatmodi.thislfctprop = zeros(1)
    gdatmodi.thislpriprop = zeros(gdat.numblpri)
    
    workdict = {}
    for strg in gdatmodi.liststrgvarbswep + gdatmodi.liststrgvarbsamp:
        valu = getattr(gdatmodi, 'this' + strg)
        if strg in gdatmodi.liststrgvarbswep:
            shap = [gdat.numbswep] + list(valu.shape)
        else:
            shap = [gdat.numbsamp] + list(valu.shape)
        workdict['list' + strg] = zeros(shap)
    
    ## saved state of the sample index used for logging progress status
    gdatmodi.cntrswepsave = -1
   
    # store the initial sample as the best fit sample
    gdatmodi.maxmllikswep = sum(gdatmodi.thisllik)
    gdatmodi.indxswepmaxmllik = -1 
    gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
   
    gdatmodi.indxstdppara = zeros(gdat.numbpara, dtype=int) - 1
    cntr = 0
    gdat.indxparaprop = zeros(gdat.numbfixpprop, dtype=int)
    for k in gdat.indxpara:
        if k in gdat.indxfixpprop:
            gdatmodi.indxstdppara[k] = cntr
            gdat.indxparaprop[cntr] = k
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
       
    # proposal scale optimization
    gdatmodi.stdvstdp = copy(gdat.stdvstdp)
    if gdat.optiprop:
        
        gdatmodi.propbrth = False      
        gdatmodi.propdeth = False      
        # perform a maximum likelihood fit to bring the sampler to a high likelihood region 
        if False:
            gdat.probtrantemp = gdat.probtran
            gdat.probtran = 0.
            listgdatmodi = [work(gdat, 0)]
            gdat.probtran = gdat.probtrantemp
            gdat.optiproptemp = False
        
        pathstdvprop = gdat.pathdataopti + '%s.fits' % gdat.rtag
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
    
        gdatmodi.stdvstdp[gdat.indxstdpcomp] = 0.
        
        gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
        gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
        gdatmodi.nextindxsamplgal, gdatmodi.nextindxsampbgal, gdatmodi.nextindxsampflux, gdatmodi.nextindxsampsind, \
                    gdatmodi.nextindxsampcurv, gdatmodi.nextindxsampexpo, gdatmodi.nextindxsampcomp = retr_indx(gdat, gdatmodi.nextindxpntsfull, gdat.spectype)
        lposcntr = retr_negalpos(gdat, gdatmodi)
        
        # temp
        deltparastep = 1e-3
        maxmstdv = 0.01
        fudgstdv = 100.
        #diffpara = zeros((2, 2, 2))
        #diffpara[0, 0, :] = deltparastep * array([-1., -1.])
        #diffpara[0, 1, :] = deltparastep * array([-1., 1.])
        #diffpara[1, 0, :] = deltparastep * array([1., -1.])
        #diffpara[1, 1, :] = deltparastep * array([1., 1.])
        diffpara = deltparastep * array([-1., 0., 1])
        lposdelt = zeros(3)
        gdatmodi.dictmodi = {}
        for strg in gdat.liststrgcomptotl:
            gdatmodi.dictmodi['stdv' + strg + 'indv'] = []
            gdatmodi.dictmodi['stdv' + strg + 'indvflux'] = []
        gdatmodi.cntrparasave = 0
        lliktemp = empty(gdat.numbstdp)
        numbiter = diffpara.size
        indxcntr = (numbiter - 1) / 2
        cntr = zeros(gdat.maxmnumbcomp)
        for k in gdat.indxpara:
            if k in gdat.indxfixpprop or k in concatenate(gdatmodi.thisindxsampcomp):
                                
                for n in range(numbiter):
                    if n == indxcntr:
                        lposdelt[n] = lposcntr
                    else:
                        lposdelt[n] = pert_llik(gdat, gdatmodi, array([k]), array([diffpara[n]]))
                
                #stdv = deltparastep * sqrt(0.5 / amax(fabs(lposcntr - lposdelt)))# * fudgstdv # / sqrt(gdat.numbpara)
                hess = 4. * fabs(lposdelt[0] + lposdelt[2] - 2. * lposdelt[1]) / deltparastep**2 
                stdv = 1. / sqrt(hess)
                
                #for n in gdat.indxpara:
                #    if n in gdat.indxfixpprop or n in concatenate(gdatmodi.thisindxsampcomp):
                #        if k == n:
                #            gdat.hess[k, n] = 1. / 4. / deltparastep * (lposdelt[1, 1] - lposdelt[1, 0] - lposdelt[0, 1] + lposdelt[0, 0])
                #            for a in range(2):
                #                
                #        else:
                #            f
                #            for a in range(2):
                #                for b in range(2):
                #                    lposdelt[a, b] = pert_llik(gdat, gdatmodi, array([k, n]), diffpara[a, b, :])
                #                    print 'a'
                #                    print a
                #                    print 'b'
                #                    print b
                #                    print 'lposdelt[a, b]'
                #                    print lposdelt[a, b]
                #                    print

                print 'gdat.namepara[k]'
                print gdat.namepara[k]
                print 'stdv'
                print stdv
                print

                if stdv > maxmstdv or not isfinite(stdv):
                    stdv = maxmstdv
                
                if k in concatenate(gdatmodi.thisindxsampcomp):
                    cntr = 0
                    indxpnts = (k - gdat.indxsampcomp[0])
                    for strg in gdat.liststrgcomptotl:
                        if k in concatenate(getattr(gdatmodi, 'thisindxsamp' + strg)):
                            indxsampflux = k + 2 - cntr
                            gdatmodi.dictmodi['stdv' + strg + 'indv'].append(stdv)
                            gdatmodi.stdvstdp[gdatmodi.indxstdppara[k]] += stdv * (gdatmodi.thissampvarb[indxsampflux] / gdat.minmflux)**0.5
                            gdatmodi.dictmodi['stdv' + strg + 'indvflux'].append(gdatmodi.thissampvarb[indxsampflux])
                        cntr += 1
                else:
                    gdatmodi.stdvstdp[gdatmodi.indxstdppara[k]] = stdv
            
            if gdat.verbtype > 0:
                gdatmodi.cntrparasave = tdpy.util.show_prog(k, gdat.numbpara, gdatmodi.cntrparasave, indxprocwork=indxprocwork, showmemo=True)

        for strg in gdat.liststrgcomptotl:
            gdatmodi.dictmodi['stdv' + strg + 'indv'] = array(gdatmodi.dictmodi['stdv' + strg + 'indv'])
            gdatmodi.dictmodi['stdv' + strg + 'indvflux'] = array(gdatmodi.dictmodi['stdv' + strg + 'indvflux'])
        
        for k in gdat.indxstdp:
            if k in gdat.indxstdpcomp:
                gdatmodi.stdvstdp[k] /= sum(gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts])
            # temp
            if k == 0 or k == 1:
                gdatmodi.stdvstdp[k] = 0.05
                
        gdatmodi.optidone = True
   
        if gdat.makeplot:
            
            if gdat.numbproc > 1:
                gdat.lock.acquire()
        
            xdat = gdat.indxstdp
            ydat = gdatmodi.stdvstdp
            path = gdat.pathopti + 'stdv%d.pdf' % indxprocwork
            tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist')
            
            for strg in gdat.liststrgcomptotl:
                path = gdat.pathopti + 'stdv' + strg + '.pdf'
                xdat = [gdatmodi.dictmodi['stdv' + strg + 'indvflux'] * gdat.fluxfactplot, gdat.meanfluxplot * gdat.fluxfactplot]
                ydat = [gdatmodi.dictmodi['stdv' + strg + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strg)] / (gdat.meanfluxplot / gdat.minmflux)**0.5]
                lablxdat = gdat.lablfeattotl['flux']
                scalxdat = gdat.dictglob['scalfluxplot']
                limtxdat = array(gdat.dictglob['limtfluxplot']) * gdat.fluxfactplot
                tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                 lablydat=r'$\sigma_{%s}$' % gdat.lablfeat[strg], plottype=['scat', 'line'])
                #tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                #                                 lablydat=r'$\sigma_{%s}$%s' % (gdat.lablfeat[strg], gdat.lablfeatunit[strg]), plottype=['scat', 'line'])
           
            if gdat.numbproc > 1:
                gdat.lock.release()
            
    else:
        gdatmodi.optidone = True
        if gdat.verbtype > 0 and indxprocwork == 0:
            print 'Skipping proposal scale optimization...'
   
    # log the initial state
    if gdat.verbtype > 0 and gdat.numbtrap > 0:
        tdpy.util.show_memo(gdatmodi, 'gdatmodi')

    if gdat.verbtype > 0:
        print 'Sampling...'

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
        gdatmodi.thisaccppsfn = True
        gdatmodi.thisaccpprio = True
    
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
            print gdatmodi.thislliktotl
            print 'thislpritotl'
            print gdatmodi.thislpritotl
            print 'thislpostotl'
            print gdatmodi.thislliktotl + gdatmodi.thislpritotl
            print

        # propose the next sample
        timeinit = gdat.functime()
        
        retr_prop(gdat, gdatmodi)
        
        timefinl = gdat.functime()
        gdatmodi.thischrototl[1] = timefinl - timeinit
       
        # diagnostics
        if gdat.diagmode:
            # sanity checks
            indxsampbadd = where((gdatmodi.thissamp[gdat.numbpopl:] > 1.) | (gdatmodi.thissamp[gdat.numbpopl:] < 0.))[0] + 1
            if indxsampbadd.size > 0:
                print 'cntrswep'
                print gdatmodi.cntrswep
                print 'thisnameproptype'
                print gdat.nameproptype[gdatmodi.thisindxproptype]
                print 'indxsampbadd'
                print indxsampbadd
                print 'thissamp'
                print gdatmodi.thissamp[indxsampbadd]
                raise Exception('Unit sample vector went outside [0,1].')
            
            gdatmodi.thislpostotl = gdatmodi.thislliktotl + gdatmodi.thislpritotl
            if False and gdatmodi.thislpostotl - gdatmodi.thislpostotlprev < -30.:
                print 'gdatmodi.thislpostotl'
                print gdatmodi.thislpostotl
                print 'gdatmodi.thislpostotlprev'
                print gdatmodi.thislpostotlprev
                raise Exception('loglikelihood drop is very unlikely!')
            gdatmodi.thislpostotlprev = gdatmodi.thislpostotl
       
            if not isfinite(gdatmodi.nextsampvarb).all() or not isfinite(gdatmodi.thissampvarb).all() or not isfinite(gdatmodi.thissamp).all():
                raise Exception('Sample vector is not finite.')
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
                    raise Exception('Bad population index')
            except:
                pass
            
            if (gdatmodi.thismodlcnts <= 0.).any() or not (isfinite(gdatmodi.thismodlcnts)).all():
                raise Exception('Current flux model is not positive')
            if gdatmodi.cntrswep > 0:
                if isnan(gdatmodi.thislpri).any():
                    raise Exception('Delta log-prior is not finite.')
                if not isfinite(gdatmodi.thisdeltllik):
                    raise Exception('Delta log-likelihood is not finite.')

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
                
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]]
                indxtemp = where((flux < gdat.minmflux) | (flux > gdat.maxmflux))[0]
                if indxtemp.size > 0:
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
                    print
                    raise Exception('Spectrum of a PS went outside the prior range.')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            timeinit = gdat.functime()
        
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this')
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]
            
            # fill the sample lists
            for strg in gdatmodi.liststrgvarbsamp:
                valu = getattr(gdatmodi, 'this' + strg)
                workdict['list' + strg][indxsampsave, ...] = valu
            gdatmodi.listindxpntsfull.append(deepcopy(gdatmodi.thisindxpntsfull))
            gdatmodi.listspecassc.append(deepcopy(gdatmodi.thisspecassc))

            timefinl = gdat.functime()
            gdatmodi.thischrototl[2] = timefinl - timeinit

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
            gdatmodi.thischrototl[3] = timefinl - timeinit
    
            if gdat.numbproc > 1:
                gdat.lock.release()
            
        if gdat.pntstype == 'lght':
            if gdat.psfntype == 'doubking':
                if gdatmodi.nextsampvarb[gdat.indxfixppsfp[1]] >= gdatmodi.nextsampvarb[gdat.indxfixppsfp[3]]:
                    for k in range(20):
                        print 'Proposal rejected due to PSF'
                    gdatmodi.thisaccppsfn = False
                    print 'gdatmodi.nextsampvarb'
                    print gdatmodi.nextsampvarb
                    print 'gdatmodi.nextsampvarb[gdat.indxfixppsfp]'
                    print gdatmodi.nextsampvarb[gdat.indxfixppsfp]
                    print 'gdatmodi.propbrth'
                    print gdatmodi.propbrth
            elif gdat.psfntype == 'doubgaus':
                if gdatmodi.nextsampvarb[gdat.indxfixppsfp[1]] >= gdatmodi.nextsampvarb[gdat.indxfixppsfp[2]]:
                    gdatmodi.thisaccppsfn = False
                    for k in range(20):
                        print 'Proposal rejected due to PSF'
                    print 'gdatmodi.nextsampvarb[gdat.indxfixppsfp]'
                    print gdatmodi.nextsampvarb[gdat.indxfixppsfp]
                    print 'gdatmodi.propbrth'
                    print gdatmodi.propbrth
        
        gdatmodi.thisaccpprop = gdatmodi.thisaccpprio and gdatmodi.thisaccppsfn
        if gdatmodi.thisaccpprop:
            
            # evaluate the log-prior
            timeinit = gdat.functime()
            
            proc_samp(gdat, gdatmodi, 'next')
            
            timefinl = gdat.functime()
            gdatmodi.thischrototl[4] = timefinl - timeinit
   
            if gdat.diagmode:
                if not isfinite(gdatmodi.thisdeltllik):
                    raise Exception('deltllik is not finite.')
            
            # evaluate the acceptance probability
            accpprob = exp(gdatmodi.thisdeltllik + gdatmodi.nextlpritotl - gdatmodi.thislpritotl + gdatmodi.thislpautotl + gdatmodi.thislfctprop + \
                                                                                                                                gdatmodi.thisjcbnfact + gdatmodi.thiscombfact)
            
            if gdat.verbtype > 1:
                print 'deltllik'
                print gdatmodi.thisdeltllik
                print
        else:
            accpprob = 0.
            gdatmodi.thisdeltllik = 0.
    
        # accept the sample
        if accpprob >= rand():

            if gdat.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(gdat, gdatmodi)

            # check if the accepted sample has maximal likelihood
            if gdatmodi.thislliktotl > gdatmodi.maxmllikswep:
                gdatmodi.maxmllikswep = gdatmodi.thislliktotl
                gdatmodi.indxswepmaxmllik = gdatmodi.cntrswep
                gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
            
            # register the sample as accepted
            gdatmodi.thisaccp = True

        # reject the sample
        else:

            if gdat.verbtype > 1:
                print 'Rejected.'

            gdatmodi.thisaccp = False
        
        # save the execution time for the sweep
        if not thismakefram:
            timetotlfinl = gdat.functime()
            gdatmodi.thischrototl[0] = timetotlfinl - timetotlinit
       
        ## variables to be saved for each sweep
        for strg in gdatmodi.liststrgvarbswep:
            workdict['list' + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi, 'this' + strg)
        
        # log the progress
        if gdat.verbtype > 0:
            minm = max(0, gdatmodi.cntrswep - 100)
            maxm = gdatmodi.cntrswep + 1
            accp = 100. * where(workdict['listaccp'][minm:maxm])[0].size / float(maxm - minm)
            accpprio = 100. * where(workdict['listaccpprio'][minm:maxm])[0].size / float(maxm - minm)
            if maxm - minm < 2:
                accp = None
            gdatmodi.cntrswepsave = tdpy.util.show_prog(gdatmodi.cntrswep, gdat.numbswep, gdatmodi.cntrswepsave, indxprocwork=indxprocwork, showmemo=True, \
                                                                                                                                            accp=accp, accpprio=accpprio)
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
        
            cntrproptotl[0] += 1.
            if gdatmodi.listaccp[gdatmodi.cntrswep]:
                cntrprop[0] += 1.
            
            if gdatmodi.cntrswepopti % perdpropeffi == 0 and (cntrproptotl[gdat.indxfixpprop] > 0).all():
            
                thispropeffi = cntrprop / cntrproptotl
                if gdat.verbtype > 0:
                    print 'Proposal scale optimization step %d' % gdatmodi.cntrswepoptistep
                    print 'Current proposal efficiency'
                    print thispropeffi[gdat.indxstdp]
                    print '%.3g +- %.3g' % (mean(thispropeffi[gdat.indxstdp]), std(thispropeffi[gdat.indxstdp])) 
                if (thispropeffi > minmpropeffi).all() and (thispropeffi < maxmpropeffi).all():
                    if gdat.verbtype > 0:
                        print 'Optimized proposal scale: '
                        print gdat.stdvstdp[gdat.indxstdp]
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
                        print factcorr[gdat.indxstdp]
                        print 'Current proposal scale: '
                        print gdat.stdvstdp[gdat.indxstdp]
                        print
                gdatmodi.cntrswepoptistep += 1
            gdatmodi.cntrswepopti += 1
  
    for strg in gdatmodi.liststrgvarbsamp + gdatmodi.liststrgvarbswep:
        valu = workdict['list' + strg]
        setattr(gdatmodi, 'list' + strg, valu)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    return gdatmodi
