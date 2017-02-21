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
         emptsamp=False, \

         # sampler
         numbswep=None, \
         numbburn=None, \
         factthin=None, \

         seedstat=None, \
         indxevttincl=None, \
         indxenerincl=None, \
       
         # evaluate the likelihood inside circles around elements
         evalcirc=None, \
        
         numbspatdims=2, \
         
         # model type
         pntstype='lght', \
        
         # initialization type
         inittype=None, \
         
         procrtag=None, \
         loadvaripara=False, \
         
         propcova=True, \
         optillik=False, \
         optiprop=False, \
         optipropsimp=True, \
         regulevi=False, \
         
         strgexprflux=None, \
         strgcatl=None, \
         anglassc=None, \
         nameexpr=None, \
         
         lgalprio=None, \
         bgalprio=None, \

         strgexpo=None, \
         
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         pixltype=None, \
         
         asscmetrtype='dist', \

         # plotting
         numbswepplot=10000, \
         makeplot=True, \
         makeplotfram=True, \
         numbframpost=None, \
         makeplotintr=False, \
         scalmaps='asnh', \
         makeanim=True, \
         strgenerfull=None, \
         strgexprname=None, \
         strganglunit=None, \
         strganglunittext=None, \
         anglfact=None, \
         fluxfactplot=None, \
         
         # model
         ## PSF
         specfraceval=0.1, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=400, \

         lablgangunit=None, \
         labllgal=None, \
         lablbgal=None, \
         lablfluxunit=None, \
         lablflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxevttfull=None, \
         binsenerfull=None, \
         asymfluxprop=False, \
         psfninfoprio=None, \
         ## spectral

         # prior
         priotype='logt', \
         priofactdoff=0., \
         # temp
         margfactmodl=1., \
         bindprio=False, \
         maxmgangdata=None, \
    
         # proposals
         stdvprophypr=0.01, \
         stdvproppsfp=0.1, \
         stdvpropbacp=0.01, \
         stdvproplenp=1e-4, \
         stdvlgal=0.001, \
         stdvbgal=0.001, \
         stdvflux=0.001, \
         stdvspep=0.001, \
         stdvspmrsind=0.2, \
         varistdvlbhl=True, \
         
         propfluxdist=True, \
         propsinddist=True, \
         propnumbpnts=True, \
         prophypr=True, \
         proppsfp=True, \
         propbacp=True, \
         proplenp=True, \
         propcomp=True, \
         probrand=0.0, \
         probtran=None, \
         probbrde=1., \
         
         radispmr=None, \

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
    
    # copy all provided inputs to the global object
    for strg, valu in args.iteritems():
        setattr(gdat, strg, valu)

    if gdat.procrtag != None:
        path = gdat.pathdata + gdat.procrtag + '/outp.fits'
        gdat = pf.getdata(path, 1)
        pf.writeto(gdat, gdat.stdvstdp, clobber=True)

    # preliminary setup
    gdat.numbfluxdistnorm = 4

    # defaults
    if gdat.strgexprflux == None:
        gdat.datatype = 'mock'
    else:
        gdat.datatype = 'inpt'
  
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

    if gdat.lablflux == None:
        if gdat.pntstype == 'lens':
            gdat.lablflux = r'\alpha_s'
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
    if gdat.anglassc == None:
        if gdat.exprtype == 'ferm':
            gdat.anglassc = 0.5 / gdat.anglfact
        if gdat.exprtype == 'hubb':
            gdat.anglassc = 0.15 / gdat.anglfact
        if gdat.exprtype == 'chan' or gdat.exprtype == 'sdss':
            gdat.anglassc = 0.5 / gdat.anglfact
    
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
        gdat.meanener = sqrt(gdat.binsener[1:] * gdat.binsener[:-1])
        gdat.deltener = gdat.binsener[1:] - gdat.binsener[:-1]
        gdat.minmener = gdat.binsener[0]
        gdat.maxmener = gdat.binsener[-1]
        for strg in ['plot']:
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
       
    ### spatial extent of the data
    if gdat.maxmgangdata == None:
        if gdat.exprtype == 'chan':
            gdat.maxmgangdata = 0.492 / gdat.anglfact * gdat.numbsidecart / 2.
        if gdat.exprtype == 'ferm':
            gdat.maxmgangdata = 20. / gdat.anglfact
        if gdat.exprtype == 'sdyn':
            gdat.maxmgangdata = 1.
        if gdat.exprtype == 'hubb':
            gdat.maxmgangdata = 2. / gdat.anglfact
    
    ### experimental PSFs
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
    if gdat.exprtype == 'chan':
        retr_chanpsfn(gdat)
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)
    if gdat.exprtype == 'hubb':
        retr_hubbpsfn(gdat)
    if gdat.exprtype == 'sdyn':
        psfp = array([0.1 / gdat.anglfact])
   
    # number of processes
    if gdat.numbproc == None:
        gdat.strgproc = os.uname()[1]
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu':
            gdat.numbproc = 20
        else:
            gdat.numbproc = 1
    
    # number of total sweeps
    if gdat.numbswep == None:
        gdat.numbswep = 100000

    # number of burned sweeps
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 5

    # factor by which to thin the sweeps to get samples
    if gdat.factthin == None:
        gdat.factthin = int(ceil(1e-3 * (gdat.numbswep - gdat.numbburn) * gdat.numbproc))

    if gdat.strgcatl == None:
        if gdat.datatype == 'mock':
            gdat.strgcatl = 'Mock'
        else:
            if gdat.exprtype == 'ferm':
                gdat.strgcatl = '3FGL'
            else:
                gdat.strgcatl = gdat.strgexprname
    
    # evaluation of the PSF
    if gdat.pntstype == 'lens':
        gdat.evalpsfnpnts = False
    else:
        gdat.evalpsfnpnts = True

    ## generative model
    setp_true(gdat, 'minmnumbpnts', array([0]))
    setp_true(gdat, 'maxmnumbpnts', array([1000]))
    
    # set mock sample vector indices
    setp_true(gdat, 'numbpnts', array([50]))
    gdat.truenumbpopl = gdat.truenumbpnts.size
    gdat.trueindxpopl = arange(gdat.truenumbpopl)
    
    setp_true(gdat, 'minmnumbpnts', zeros(gdat.truenumbpopl))
    setp_true(gdat, 'maxmnumbpnts', zeros(gdat.truenumbpopl) + 500)

    ### element parameter distributions
    setp_true(gdat, 'spatdisttype', ['unif' for l in gdat.trueindxpopl])
    setp_true(gdat, 'fluxdisttype', ['powr' for l in gdat.trueindxpopl])
    setp_true(gdat, 'sinddisttype', ['gaus' for l in gdat.trueindxpopl])
    setp_true(gdat, 'spectype', ['powr' for l in gdat.trueindxpopl])
    
    ### PSF model
    #### angular profile
    if gdat.exprtype == 'ferm':
        gdat.exprpsfntype = 'doubking'
    if gdat.exprtype == 'chan':
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'sdss':
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'hubb':
        gdat.exprpsfntype = 'singgaus'
    if gdat.exprtype == 'sdyn':
        gdat.exprpsfntype = 'singgaus'
    
    if psfninfoprio == None:
        if gdat.exprtype == 'ferm':
            psfninfoprio = True
        else:
            psfninfoprio = False

    psfntype = gdat.exprpsfntype
    setp_true(gdat, 'psfntype', psfntype)
    
    #### off-axis profile
    if gdat.exprtype == 'chan':
        oaxitype = True
    else:
        oaxitype = False
    setp_true(gdat, 'oaxitype', oaxitype)

    ### background
    #### template
    if gdat.pntstype == 'lght':
        if gdat.exprtype == 'ferm':
            back = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
        elif gdat.exprtype == 'chan':
            back = ['chanfluxisot.fits']
        else:
            back = [1.]
    if gdat.pntstype == 'lens':
        back = [1e-7]
    setp_true(gdat, 'back', back)
    
    #### label
    lablback = [r'$\mathcal{I}$']
    if gdat.exprtype == 'ferm':
        lablback.append(r'$\mathcal{D}$')
    setp_true(gdat, 'lablback', lablback)
    
    #### name
    strgback = ['fluxisot']
    if gdat.exprtype == 'ferm':
        strgback.append('fluxfdfm')
    setp_true(gdat, 'strgback', strgback)
    
    #### legend
    nameback = ['Isotropic']
    if gdat.exprtype == 'ferm':
        nameback.append(r'FDM')
    setp_true(gdat, 'nameback', nameback)
    
    retr_indxsamp(gdat, strgpara='true')
    
    ### maximum horizontal/vertical distance of the elements from the image center
    gdat.maxmgang = gdat.maxmgangdata * gdat.margfactmodl

    ### PSF parameters
    if gdat.psfninfoprio:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.truepsfntype == 'doubking' or gdat.truepsfntype == 'doubgaus' or gdat.truepsfntype == 'gausking':
                    meanpsff = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener]
                    stdvpsff = meanpsff * 0.1
                    setp_truedefa(gdat, 'psffene%devt%d' % (i, m), [meanpsff, stdvpsff], typelimt='meanstdv')
                meansigc = gdat.exprpsfp[i * gdat.truenumbpsfptotl + m * gdat.truenumbpsfptotl * gdat.numbener + 1]
                stdvsigc = meansigc * 0.1
                setp_truedefa(gdat, 'sigcene%devt%d' % (i, m), [meansigc, stdvsigc], typelimt='meanstdv')
                if gdat.truepsfntype == 'doubking' or gdat.truepsfntype == 'singking':
                    meangamc = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 2]
                    stdvgamc = meangamc * 0.1
                    setp_truedefa(gdat, 'gamcene%devt%d' % (i, m), [meangamc, stdvgamc], typelimt='meanstdv')
                if gdat.truepsfntype == 'doubking' or gdat.truepsfntype == 'doubgaus' or gdat.truepsfntype == 'gausking':
                    meansigt = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 3]
                    stdvsigt = meansigt * 0.1
                    setp_truedefa(gdat, 'sigtene%devt%d' % (i, m), [meansigt, stdvsigt], typelimt='meanstdv')
                if gdat.truepsfntype == 'doubking' or gdat.truepsfntype == 'gausking':
                    meangamt = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 4]
                    stdvgamt = meangamt * 0.1
                    setp_truedefa(gdat, 'gamtene%devt%d' % (i, m), [meangamt, stdvgamt], typelimt='meanstdv')
    else:
        if gdat.exprtype == 'ferm':
            minmsigm = 0.1
            maxmsigm = 10.
        if gdat.exprtype == 'hubb':
            minmsigm = 0.01 / gdat.anglfact
            maxmsigm = 0.1 / gdat.anglfact
        if gdat.exprtype == 'chan':
            minmsigm = 0.1 / gdat.anglfact
            maxmsigm = 2. / gdat.anglfact
        minmgamm = 1.5
        maxmgamm = 20.
        minmonor = 0.01
        maxmonor = 1.
        minmoind = 1.5
        maxmoind = 2.5
        setp_truedefa(gdat, 'sigc', [minmsigm, maxmsigm], ener=True, evtt=True)
        setp_truedefa(gdat, 'sigt', [minmsigm, maxmsigm], ener=True, evtt=True)
        setp_truedefa(gdat, 'gamc', [minmgamm, maxmgamm], ener=True, evtt=True)
        setp_truedefa(gdat, 'gamt', [minmgamm, maxmgamm], ener=True, evtt=True)
        setp_truedefa(gdat, 'onor', [minmonor, maxmonor], ener=True, evtt=True)
        setp_truedefa(gdat, 'oind', [minmoind, maxmoind], ener=True, evtt=True)
    setp_truedefa(gdat, 'psff', [0., 1.], ener=True, evtt=True)
 
    ### normalization
    setp_truedefa(gdat, 'bacp', [0.5, 2.], ener=True, back=True)

    ## hyperparameters
    setp_truedefa(gdat, 'meanpnts', [1., 1e3], popl=True)
    
    ### element parameter boundaries
    #### spatial
    if gdat.exprtype == 'ferm':
        minmgang = 1e-1 / gdat.anglfact
    else:
        minmgang = 1e-2 / gdat.anglfact
    setp_true(gdat, 'minmgang', minmgang)
    
    if gdat.exprtype == 'ferm':
        minmflux = 5e-11
    if gdat.exprtype == 'chan':
        minmflux = 5e-10
    if gdat.exprtype == 'sdyn':
        minmflux = 1e0
    if gdat.pntstype == 'lens':
        minmflux = 0.04 / gdat.anglfact
    setp_true(gdat, 'minmflux', minmflux)
    
    if gdat.exprtype == 'ferm':
        maxmflux = 1e-7
    if gdat.exprtype == 'chan':
        maxmflux = 1e-7
    if gdat.exprtype == 'sdyn':
        maxmflux = 1e4
    if gdat.pntstype == 'lens':
        maxmflux = 0.5 / gdat.anglfact
    setp_true(gdat, 'maxmflux', maxmflux)
   
    # parameter defaults
    ## distribution
    ### flux
    setp_truedefa(gdat, 'gangdistscal', [1. / gdat.anglfact, 10. / gdat.anglfact], popl=True)
    setp_truedefa(gdat, 'spatdistcons', [1e-2, 1e0], popl=True)
    setp_truedefa(gdat, 'bgaldistscal', [0.5 / gdat.anglfact, 5. / gdat.anglfact], popl=True)
    setp_truedefa(gdat, 'fluxdistslop', [1., 4.], popl=True)
    
    for k in range(gdat.numbfluxdistnorm):
        setp_truedefa(gdat, 'fluxdistnormbin%d' % k, [1e-3, 1e3], popl=True)

    ### spectral index
    if gdat.numbener > 1:
        sind = [1., 3.]
        setp_truedefa(gdat, 'sind', sind, popl=True)
        setp_truedefa(gdat, 'sinddistmean', sind, popl=True)
        #### standard deviations should not be too small
        setp_truedefa(gdat, 'sinddiststdv', [0.5, 2.], popl=True)
        setp_truedefa(gdat, 'curvdistmean', [-1., 1.], popl=True)
        setp_truedefa(gdat, 'curvdiststdv', [0.1, 1.], popl=True)
        setp_truedefa(gdat, 'expodistmean', [1., 8.], popl=True)
        setp_truedefa(gdat, 'expodiststdv', [0.01 * gdat.maxmener, gdat.maxmener], popl=True)
   
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
    
    ## lensing
    setp_truedefa(gdat, 'lgalsour', [gdat.minmlgal, gdat.maxmlgal])
    setp_truedefa(gdat, 'bgalsour', [gdat.minmbgal, gdat.maxmbgal])
    setp_truedefa(gdat, 'lgalhost', [gdat.minmlgal, gdat.maxmlgal])
    setp_truedefa(gdat, 'bgalhost', [gdat.minmbgal, gdat.maxmbgal])
    for i in gdat.indxener:
        setp_truedefa(gdat, 'specsourene%d' % i, array([1e-21, 1e-17]))
    setp_truedefa(gdat, 'sizesour', [0.1 / gdat.anglfact, 2. / gdat.anglfact])
    setp_truedefa(gdat, 'ellpsour', [0., 0.3])
    for i in gdat.indxener:
        setp_truedefa(gdat, 'spechostene%d' % i, array([1e-21, 1e-17]))
    setp_truedefa(gdat, 'sizehost', [0.1 / gdat.anglfact, 2. / gdat.anglfact])
    setp_truedefa(gdat, 'beinhost', [0.5 / gdat.anglfact, 2. / gdat.anglfact])
    setp_truedefa(gdat, 'ellphost', [0., 0.5])
    setp_truedefa(gdat, 'sherhost', [0., 0.3])
    setp_truedefa(gdat, 'anglsour', [0., pi])
    setp_truedefa(gdat, 'anglhost', [0., pi])
    setp_truedefa(gdat, 'sanghost', [0., pi])
    
    # copy the true model to the inference model if the inference model parameter has not been specified
    temp = deepcopy(gdat.__dict__)
    for strg, valu in temp.iteritems():
        if strg.startswith('true'):
            try:
                valumodl = getattr(gdat, strg[4:])
                if valumodl == None:
                    raise
            except:
                setattr(gdat, strg[4:], getattr(gdat, strg))
   
    # get the time stamp
    gdat.strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # check the call stack for the name of the configuration function
    gdat.strgcnfg = inspect.stack()[1][3]
   
    # check inputs
    if gdat.numbburn > gdat.numbswep:
        raise Exception('Bad number of burn-in sweeps.')
    if gdat.factthin > gdat.numbswep - gdat.numbburn:
        raise Exception('Bad thinning factor.')
    if gdat.pixltype == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')
    if gdat.minmflux >= gdat.maxmflux:
        raise Exception('Minimum flux is greater than maximum flux.')

    if gdat.verbtype > 0:
        print 'PCAT v0.1 started at %s' % gdat.strgtimestmp
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)

    setp_true(gdat, 'spatdistcons', 1e-1, popl=True)
    setp_true(gdat, 'gangdistscal', 4. / gdat.anglfact, popl=True)
    setp_true(gdat, 'bgaldistscal', 2. / gdat.anglfact, popl=True)
    if gdat.pntstype == 'lens':
        fluxdistslop = 3.
    else:
        fluxdistslop = 2.6
    setp_true(gdat, 'fluxdistslop', fluxdistslop, popl=True)

    retr_axis(gdat, 'flux', gdat.trueminmflux, gdat.maxmflux, gdat.numbfluxdistnorm - 1, scal='logt', strg='true')
    fluxdistnorm = gdat.truebinsflux**(-fluxdistslop)
    fluxdistnorm *= 1e-2 / amin(fluxdistnorm)
    for k in range(gdat.numbfluxdistnorm):
        setp_true(gdat, 'fluxdistnormbin%d' % k, fluxdistnorm[k], popl=True)
    
    setp_true(gdat, 'sinddistmean', 2.15, popl=True)
    setp_true(gdat, 'sinddiststdv', 1., popl=True)
    
    setp_true(gdat, 'curvdistmean', 2., popl=True)
    setp_true(gdat, 'curvdiststdv', 0.2, popl=True)
    
    setp_true(gdat, 'expodistmeanpop1', 2., popl=True)
    setp_true(gdat, 'expodiststdvpop1', 0.2, popl=True)
    
    setp_true(gdat, 'bacp', 1.)
    
    if gdat.pntstype == 'lens':
        setp_true(gdat, 'beinhost', 1.5 / gdat.anglfact)
        setp_true(gdat, 'sizesour', 0.5 / gdat.anglfact)
        setp_true(gdat, 'sizehost', 1. / gdat.anglfact)
        for i in gdat.indxener:
            setp_true(gdat, 'specsourene%d' % i, 10. * gdat.hubbexpofact)
            setp_true(gdat, 'spechostene%d' % i, 20. * gdat.hubbexpofact)
   
    if gdat.defa:
        return gdat
    
    # initial setup
    setpinit(gdat, True) 
    setp_fixp(gdat, strgpara='true')
    
    # intermediate setup
    if gdat.numbener > 1:
        gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
        if gdat.enerfluxdist == 0.:
            raise Exception('Pivot energy cannot be zero.')
        gdat.factspecener = (gdat.meanener / gdat.enerfluxdist)**(-sqrt(amin(gdat.minmsinddistmeanpop0) * amax(gdat.maxmsinddistmeanpop0)))
        gdat.enerexpofact = gdat.enerfluxdist - gdat.meanener
    else:
        gdat.factspecener = array([1.])
  
    # color bars
    gdat.minmconv = 1e-2
    gdat.maxmconv = 1e1
    gdat.minmdeflcomp = 0.
    gdat.maxmdeflcomp = 1e-4
    gdat.minmlpdfspatpriointp = -8.#2. * log(1. / 2. / gdat.maxmgang) - 2.
    gdat.maxmlpdfspatpriointp = -15.#2. * log(1. / 2. / gdat.maxmgang) + 2.
    gdat.maxmllik = 0.
    gdat.minmllik = -1e2
    
    gdat.cmapdatacnts = 'Reds'
    gdat.cmapresicnts = 'RdBu'
    gdat.cmapdeflcomp = 'Oranges'
    gdat.cmapconv = 'Purples'
    gdat.cmapllik = 'YlGn'
    gdat.cmaplpdfspatpriointp = 'PuBu'
    
    liststrgcbar = ['conv', 'deflcomp', 'llik', 'lpdfspatpriointp']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)
    
    # element features
    if gdat.numbener > 1:
        # temp
        gdat.minmexpo = gdat.minmexpodistmeanpop0
        gdat.maxmexpo = gdat.maxmexpodistmeanpop0
    gdat.minmaang = -pi
    gdat.maxmaang = pi
    gdat.minmgangplot = 0.
    gdat.maxmgangplot = gdat.maxmlgal
    gdat.minmdeltllik = -1.
    gdat.maxmdeltllik = 4.
    gdat.minmdistsour = 0.
    gdat.maxmdistsour = 3. * gdat.maxmgang
    gdat.minmdotpsour = 0.
    gdat.maxmdotpsour = 5e-4

    print 'gdat.psfninfoprio'
    print gdat.psfninfoprio
    print

    ## plot limits for element parameters
    for strgfeat in gdat.liststrgfeat:
        for strglimt in ['minm', 'maxm']:
            
            # temp
            if strgfeat == 'spec' or strgfeat == 'cnts':
                continue

            try:
                truelimt = getattr(gdat, 'true' + strglimt + strgfeat)
            except:
                truelimt = None
             
            if strgfeat in ['sind', 'curv']:
                if strglimt == 'minm':
                    limt = getattr(gdat, 'minm' + strgfeat + 'distmeanpop0') - getattr(gdat, 'maxm' + strgfeat + 'diststdvpop0')
                else:
                    limt = getattr(gdat, 'maxm' + strgfeat + 'distmeanpop0') + getattr(gdat, 'maxm' + strgfeat + 'diststdvpop0')
            else:
                limt = getattr(gdat, strglimt + strgfeat)
            
            setattr(gdat, strglimt + strgfeat, limt)

            if truelimt == None:
                setattr(gdat, strglimt + strgfeat + 'plot', limt)
            else:
                limt = min(limt, truelimt)
                setattr(gdat, strglimt + strgfeat + 'plot', limt)
    
    gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
    gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
    gdat.minmcntsplot = 0.1 * gdat.minmflux * mean(mean(gdat.expo, 1), 1)
    gdat.maxmcntsplot = gdat.maxmflux * mean(mean(gdat.expo, 1), 1)
    if gdat.enerdiff:
        gdat.minmcntsplot *= gdat.deltener * gdat.factspecener
        gdat.maxmcntsplot *= gdat.deltener * gdat.factspecener

    # temp
    gdat.minmspec = gdat.minmflux * gdat.factspecener
    gdat.maxmspec = gdat.maxmflux * gdat.factspecener

    for strgfeat in gdat.liststrgfeat:
        if strgfeat == 'spec':
            gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
            gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
            gdat.binsspecplot = gdat.binsfluxplot[None, :] * gdat.factspecener[:, None]
            gdat.meanspecplot = empty((gdat.numbener, gdat.numbbinsplot))
            for i in gdat.indxener:
                gdat.meanspecplot[i, :] = sqrt(gdat.binsspecplot[i, 1:] * gdat.binsspecplot[i, :-1])
        elif strgfeat == 'cnts':
            gdat.minmcnts = 0.1 * gdat.minmflux * mean(mean(gdat.expo, 1), 1)
            gdat.maxmcnts = gdat.maxmflux * mean(mean(gdat.expo, 1), 1)
            if gdat.enerdiff:
                gdat.minmcnts *= gdat.deltener * gdat.factspecener
                gdat.maxmcnts *= gdat.deltener * gdat.factspecener
            retr_axis(gdat, 'cntsplot', gdat.minmcnts[gdat.indxenerfluxdist[0]], gdat.maxmcnts[gdat.indxenerfluxdist[0]], gdat.numbbinsplot, scal='logt')
            gdat.binscnts = empty((gdat.numbener, gdat.numbbinsplot + 1))
            for i in gdat.indxener:
                gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbbinsplot + 1) # [1]
            gdat.meancnts = sqrt(gdat.binscnts[:, :-1] * gdat.binscnts[:, 1:]) 
        else:
            if strgfeat in ['flux', 'expo']:
                scal = 'logt'
            else:
                scal = 'self'
            maxm = getattr(gdat, 'maxm' + strgfeat + 'plot')
            minm = getattr(gdat, 'minm' + strgfeat + 'plot')
            retr_axis(gdat, strgfeat + 'plot', minm, maxm, gdat.numbbinsplot, scal=scal)
            retr_axis(gdat, strgfeat + 'plotprio', minm, maxm, gdat.numbbinsplotprio, scal=scal)
    
        gdat.dictglob['limt' + strgfeat + 'plot'] = [getattr(gdat, 'minm' + strgfeat), getattr(gdat, 'maxm' + strgfeat)]






    # create the PCAT folders
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathoutpthis = gdat.pathoutp + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    os.system('mkdir -p %s' % gdat.pathoutpthis)
    pathcatllite = gdat.pathoutpthis + 'catllite.fits'  
    pathcatl = gdat.pathoutpthis + 'catl.fits'  
    
    if gdat.verbtype > 0:
        sizetotl = 0.
        for root, dirs, listfile in os.walk(gdat.pathoutp):
            for thisfile in listfile:
                sizetotl += os.path.getsize(root + '/' + thisfile) / 2**30
        if sizetotl > 10.:
            print 'Warning: PCAT data path size is %d GB' % sizetotl

    if gdat.exprtype == 'ferm':
        retr_fermdata(gdat)
        gdat.exprinfo = True
    if gdat.exprtype == 'chan':
        retr_chandata(gdat)
        gdat.exprinfo = True
    
    gdat.trueinfo = gdat.exprinfo or gdat.datatype == 'mock'

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
        for k in range(3):
            gdat.exprspec[k, :, :] = gdat.exprspec[k, :, indxpnts].T
        gdat.exprsind = gdat.exprsind[indxpnts]
        if gdat.exprcnts != None:
            gdat.exprcnts = gdat.exprcnts[:, indxpnts, :]

        # compute the catalog counts based on the exposure
        gdat.exprcntscalc = empty((gdat.numbener, gdat.exprnumbpnts, gdat.numbevtt))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                indxpixltemp = retr_indxpixl(gdat, gdat.exprbgal, gdat.exprlgal)
                gdat.exprcntscalc[i, :, m] = gdat.exprspec[0, i, :] * gdat.expo[i, indxpixltemp, m] * gdat.deltener[i]
        
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
    
    ## unit sample vector
    gdat.truesamp = zeros(gdat.truenumbpara)
    gdat.truefixp = zeros(gdat.truenumbfixp) + nan
    gdat.truefixp[gdat.trueindxfixpnumbpnts] = gdat.truenumbpnts
    
    if gdat.datatype == 'mock':
        
        gdat.trueindxpntsfull = []
        if gdat.numbtrap > 0:
            for l in gdat.trueindxpopl:
                gdat.trueindxpntsfull.append(range(gdat.truenumbpnts[l]))
        else:
            gdat.trueindxpntsfull = []

        gdat.trueindxsamplgal, gdat.trueindxsampbgal, gdat.trueindxsampflux, gdat.trueindxsampsind, gdat.trueindxsampcurv, gdat.trueindxsampexpo, \
                                                                                      gdat.trueindxsampcomp = retr_indxsampcomp(gdat, gdat.trueindxpntsfull, gdat.truespectype)
        
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
            gdat.truefixp[gdat.trueindxfixplgalsour] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalsour] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixplgalhost] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalhost] = 0.04 * gdat.maxmgang * randn()
        
        gdat.truesampvarb = empty(gdat.truenumbpara)
        gdat.truesampvarb[gdat.trueindxfixp] = gdat.truefixp
        
        gdat.truenumbpntstotl = sum(gdat.truefixp[gdat.trueindxfixpnumbpnts])
        gdat.trueindxpntstotl = arange(gdat.truenumbpntstotl)
        
        # temp
        #for l in gdat.trueindxpopl:
        #    if gdat.truespatdisttype[l] == 'gaus':
        #        gdat.truelpdfspatprio, gdat.truelpdfspatprioobjt = retr_spatprio(gdat, gdat.truesampvarb[gdat.trueindxfixpspatdistcons[l]], gdat.truepdfnspatpriotemp)
               
        gdat.truecnts = [[] for l in gdat.trueindxpopl]
        gdat.truelgal = [[] for l in gdat.trueindxpopl]
        gdat.truebgal = [[] for l in gdat.trueindxpopl]
        gdat.truegang = [[] for l in gdat.trueindxpopl]
        gdat.trueaang = [[] for l in gdat.trueindxpopl]
        gdat.truespec = [[] for l in gdat.trueindxpopl]
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
                
                #if gdat.truespatdisttype[l] == 'gaus':
                #    numbrejesamp = gdat.truenumbpnts[l] * 1
                #    lgaltemp = icdf_self(rand(numbrejesamp), -gdat.maxmgangdata, 2. * gdat.maxmgangdata)
                #    bgaltemp = icdf_self(rand(numbrejesamp), -gdat.maxmgangdata, 2. * gdat.maxmgangdata) 
                #    lgaltemp * gdat.anglfact
                #    bgaltemp * gdat.anglfact
                #    lpdftemp = gdat.truelpdfspatprioobjt(bgaltemp, lgaltemp, grid=False)
                #    probtemp = exp(lpdftemp)
                #    probtemp /= sum(probtemp)
                #    indxtemp = choice(arange(numbrejesamp), p=probtemp, size=gdat.truenumbpnts[l])
                #    gdat.truelgal[l] = lgaltemp[indxtemp]
                #    gdat.truebgal[l] = bgaltemp[indxtemp]
                
                gdat.truespec[l] = empty((3, gdat.numbener, gdat.truenumbpnts[l]))

                if gdat.truefluxdisttype[l] == 'powr':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_powr(rand(gdat.truenumbpnts[l]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                                                    gdat.truefixp[gdat.trueindxfixpfluxdistslop[l]])
                if gdat.truefluxdisttype[l] == 'bind':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_bind(rand(gdat.truenumbpnts[l]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                           gdat.truebinsflux, gdat.truefixp[gdat.trueindxfixpfluxdistnorm[l, :]])
                
                # temp -- make sure this reordering does not mess up other things
                gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = sort(gdat.truespec[l][:, gdat.indxenerfluxdist[0], :], axis=1)[::-1]
    
                if gdat.numbener > 1:
                    # spectral parameters
                    if gdat.truesinddisttype[l] == 'gaus':
                        gdat.truesind[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpsinddistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpsinddiststdv[l]])
                    if gdat.truesinddisttype[l] == 'atan':
                        gdat.truesind[l] = icdf_atan(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpsinddistmean[l]], \
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
                    
                gdat.truesampvarb[gdat.trueindxsamplgal[l]] = gdat.truelgal[l]
                gdat.truesampvarb[gdat.trueindxsampbgal[l]] = gdat.truebgal[l]
                gdat.truesampvarb[gdat.trueindxsampflux[l]] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                if gdat.numbener > 1:
                    gdat.truesampvarb[gdat.trueindxsampsind[l]] = gdat.truesind[l]
                    if gdat.truespectype[l] == 'curv':
                        gdat.truesampvarb[gdat.trueindxsampcurv[l]] = gdat.truecurv[l]
                    if gdat.truespectype[l] == 'expo':
                        gdat.truesampvarb[gdat.trueindxsampexpo[l]] = gdat.trueexpo[l]
                
                indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
                gdat.truecnts[l] = gdat.truespec[l][0, :, :, None] * gdat.expo[:, indxpixltemp, :]
                if gdat.enerdiff:
                    gdat.truecnts[l] *= gdat.deltener[:, None, None]
        
        if gdat.verbtype > 1:
            print
            print '%20s %20s %15s' % ('truenamepara', 'truesampvarb', 'truescalpara')
            for k in gdat.trueindxpara:
                if k in concatenate(gdat.trueindxsamplgal):
                    print
                print '%20s %20f %15s' % (gdat.truenamepara[k], gdat.truesampvarb[k], gdat.truescalpara[k])
    
        if gdat.lgalprio == None or gdat.bgalprio == None:
            gdat.lgalprio = concatenate((gdat.truelgal))
            gdat.bgalprio = concatenate((gdat.truebgal))
            gdat.numbspatprio = gdat.lgalprio.size

    if gdat.datatype == 'inpt':
        if gdat.lgalprio == None or gdat.bgalprio == None:
            gdat.numbspatprio = 20
            gdat.lgalprio = icdf_self(rand(gdat.numbspatprio), -gdat.maxmgangdata, 2. * gdat.maxmgangdata)
            gdat.bgalprio = icdf_self(rand(gdat.numbspatprio), -gdat.maxmgangdata, 2. * gdat.maxmgangdata)
    
    # spatial template for the catalog prior
    # temp -- this should move outside the if
    gdat.pdfnspatpriotemp = zeros((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
    for k in range(gdat.numbspatprio):
        gdat.pdfnspatpriotemp[:] += 1. / sqrt(2. * pi) / gdat.stdvspatprio * exp(-0.5 * (gdat.binslgalcartmesh - gdat.lgalprio[k])**2 / gdat.stdvspatprio**2) * \
                                                                                    exp(-0.5 * (gdat.binsbgalcartmesh - gdat.bgalprio[k])**2 / gdat.stdvspatprio**2)
    gdat.pdfnspatpriotemp /= amax(gdat.pdfnspatpriotemp)
    
    if gdat.datatype == 'mock':
        
        proc_samp(gdat, None, 'true', raww=True)
        proc_samp(gdat, None, 'true')
            
        if gdat.makeplot:
            plot_samp(gdat, None, 'true')
        
        if gdat.seedstat != None:
            seed()
    
    if gdat.datatype == 'inpt':
        retr_datatick(gdat)
    
    gdat.trueflux = [[] for l in gdat.trueindxpopl]
    for l in gdat.trueindxpopl:
        gdat.trueflux[l] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
    
    if gdat.makeplot:
        plot_init(gdat)

    # final setup
    setpfinl(gdat, True) 

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
        print 'minmcntsplot'
        print gdat.minmcntsplot
        print 'maxmcntsplot'
        print gdat.maxmcntsplot
        if gdat.evalcirc != 'full':
            print 'maxmangleval'
            print gdat.anglfact * gdat.maxmangleval, ' [%s]' % gdat.strganglunit
        print
            
    # process lock for simultaneous plotting
    gdat.lock = mp.Manager().Lock()
        
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdat, 'gdat')
    
    if gdat.verbtype > 0:
        print 'Writing the global state to the disc before spawning workers...'
    path = gdat.pathoutpthis + 'gdatinit'
    writfile(gdat, path) 
    #writoutp(gdat, path, catl=False) 
    
    if gdat.numbproc == 1:
        worktrac(gdat.pathoutpthis, gdat.lock, 0)
    else:
        if gdat.verbtype > 0:
            print 'Forking the sampler...'

        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(worktrac, gdat.pathoutpthis, gdat.lock)
        pool.map(workpart, gdat.indxproc)

        pool.close()
        pool.join()

    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        timeinit = gdat.functime()
    
    listgdatmodi = []
    for k in gdat.indxproc:
        path = gdat.pathoutpthis + 'gdatmodi%04d' % k
        listgdatmodi.append(readfile(path))
   
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
    for strg in ['deltlliktotl', 'memoresi']:
        gdat.liststrgchanflat.remove(strg)
   
    gdat.listsampvarbproc = copy(gdat.listsampvarb)

    ## other parameters
    for strg in gdat.liststrgchanflat:
        inpt = getattr(gdat, 'list' + strg)
        shap = [inpt.shape[0] * inpt.shape[1]] + list(inpt.shape[2:])
        setattr(gdat, 'list' + strg, inpt.reshape(shap))
    
    # add execution times to the chain output
    gdat.timereal = zeros(gdat.numbproc)
    gdat.timeproc = zeros(gdat.numbproc)
    for k in gdat.indxproc:
        gdat.timereal[k] = listgdatmodi[k].timereal
        gdat.timeproc[k] = listgdatmodi[k].timeproc
    
    # correct the likelihoods for the constant data dependent factorial
    llikoffs = sum(sp.special.gammaln(gdat.datacnts + 1))
    gdat.listlliktotl -= llikoffs
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
        listlliktemp = listlliktotl[indxsampregu]
    else:
        listlliktemp = gdat.listlliktotl
    gdat.levi = retr_levi(listlliktemp)
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # calculate relative entropy
    gdat.info = retr_info(gdat.listlliktotl, gdat.levi)

    # parse the sample vector
    gdat.listfixp = gdat.listsampvarb[:, gdat.indxfixp]

    # post process samples
    if gdat.numbtrap > 0:
        # collect PS parameters from the chains
        dicttemp = {}
        for strg in gdat.liststrgcomptotl:
            gdat.dictglob['list' + strg] = [[] for l in gdat.indxpopl]
        for n in gdat.indxsamptotl: 
            dicttemp['indxsamplgal'], dicttemp['indxsampbgal'], dicttemp['indxsampflux'], dicttemp['indxsampsind'], dicttemp['indxsampcurv'], \
															dicttemp['indxsampexpo'], dicttemp['indxsampcomp'] = retr_indxsampcomp(gdat, gdat.listindxpntsfull[n], gdat.spectype)
            for l in gdat.indxpopl:
                for strg in gdat.liststrgcomp[l]:
                    gdat.dictglob['list' + strg][l].append(gdat.listsampvarb[n, dicttemp['indxsamp' + strg][l]])
        for strg in gdat.liststrgcomptotl:
            setattr(gdat, 'list' + strg, gdat.dictglob['list' + strg])

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

    if gdat.verbtype > 0:
        if gdat.timeatcr == 0.:
            print 'Autocorrelation time estimation failed.'
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    ## construct a deterministic catalog
    if gdat.verbtype > 0:
        print 'Constructing a labeled catalog...'
        timeinit = gdat.functime()
    
    #retr_detrcatl(gdat)
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

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
            stdvtemp = []
            numbelem = len(listtemp[0])
            for k in range(numbelem):
                shap = [gdat.numbsamptotl] + list(listtemp[0][k].shape)
                temp = zeros(shap)
                for n in gdat.indxsamptotl:
                    temp[n, ...] = listtemp[n][k]
                posttempsing = tdpy.util.retr_postvarb(temp)
                meditempsing = posttempsing[0, ...]
                errrtempsing = tdpy.util.retr_errrvarb(posttempsing)
                stdvtempsing = std(temp)
                posttemp.append(posttempsing)
                meditemp.append(meditempsing)
                errrtemp.append(errrtempsing)
                stdvtemp.append(stdvtempsing)
        else:
            posttemp = tdpy.util.retr_postvarb(listtemp)
            meditemp = posttemp[0, ...]
            errrtemp = tdpy.util.retr_errrvarb(posttemp)
            stdvtemp = std(posttemp)
        setattr(gdat, 'post' + strg, posttemp)
        setattr(gdat, 'medi' + strg, meditemp)
        setattr(gdat, 'errr' + strg, errrtemp)
        setattr(gdat, 'stdv' + strg, stdvtemp)
    
    # temp
    if gdat.pntstype == 'lght':
        gdat.medicntsbackfwhm = retr_cntsbackfwhm(gdat, gdat.postfixp[0, gdat.indxfixpbacp], gdat.postfwhm[0, :])
        gdat.medibinssigm = retr_sigm(gdat, gdat.binscnts, gdat.medicntsbackfwhm)
   
    # memory usage
    gdat.meanmemoresi = mean(gdat.listmemoresi, 1)
    gdat.derimemoresi = (gdat.meanmemoresi[-1] - gdat.meanmemoresi[0]) / gdat.numbswep
    
    path = gdat.pathoutpthis + 'pcat.h5'
    writoutp(gdat, path)
    os.system('rm -rf %sgdat*' % gdat.pathoutpthis) 
    # temp
    if False:
        if gdat.makeplot:
            plot_post(pathcatl=pathcatl, verbtype=gdat.verbtype, makeanim=gdat.makeanim)
    else:
        if gdat.makeplot:
            plot_post(gdat=gdat, writ=False)

    gdat.timerealtotl = time.time() - gdat.timerealtotl
    gdat.timeproctotl = time.clock() - gdat.timeproctotl
    gdat.timeproctotlswep = gdat.timeproctotl / gdat.numbswep
    if gdat.timeatcr == 0.:
        gdat.timeprocnorm = 0.
    else:
        gdat.timeprocnorm = gdat.timeproctotlswep / gdat.timeatcr

    if gdat.verbtype > 0:
        for k in gdat.indxproc:
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdat.timereal[k], gdat.timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (gdat.timerealtotl, gdat.timeproctotl)
        if gdat.makeplot:
            print 'The plots are in ' + gdat.pathplot
        print 'The ensemble of catalogs is at ' + pathcatl
    return gdat
    

def retr_deltlpos(gdat, gdatmodi, indxparapert, stdvparapert):

    numbpert = indxparapert.size 
    
    gdatmodi.thissamp = copy(gdatmodi.thissamptemp)
    for k in range(numbpert):
        gdatmodi.thissamp[indxparapert[k]] += stdvparapert[k]
    
    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.thissamp, 'this')
    
    proc_samp(gdat, gdatmodi, 'this', fast=True)
    if False and indxparapert[0] in concatenate(gdatmodi.thisindxsampbgal):
        plot_samp(gdat, gdatmodi, 'this')
   
    deltlpos = gdatmodi.thislpostotl - gdatmodi.thislpostotltemp
    deltlpri = gdatmodi.thislpritotl - gdatmodi.thislpritotltemp
    deltllik = gdatmodi.thislliktotl - gdatmodi.thislliktotltemp
    
    if True:
        print 'retr_deltlpos'
        print 'gdatmodi.thislpostotltemp'
        print gdatmodi.thislpostotltemp
        print 'gdatmodi.thislpritotltemp'
        print gdatmodi.thislpritotltemp
        print 'gdatmodi.thislliktotltemp'
        print gdatmodi.thislliktotltemp
        print 'gdatmodi.thislpostotl'
        print gdatmodi.thislpostotl
        print 'gdatmodi.thislpritotl'
        print gdatmodi.thislpritotl
        print 'gdatmodi.thislliktotl'
        print gdatmodi.thislliktotl
        print 'deltlpos'
        print deltlpos
        print 'deltllik'
        print deltllik
        print 'deltlpri'
        print deltlpri
        print
    
    return deltlpos


#def worktrac(gdat, indxprocwork):
def worktrac(pathoutpthis, lock, indxprocwork):
	
    try:
        #return work(gdat, indxprocwork)
        return work(pathoutpthis, lock, indxprocwork)
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def optiprop(gdat, gdatmodi, indxprocwork):
    
    gdatmodi.propbrth = False      
    gdatmodi.propdeth = False      
    
    pathstdvprop = gdat.pathdataopti + '%s.fits' % gdat.rtag
    if os.path.isfile(pathstdvprop) and gdat.loadvaripara:
        if gdat.verbtype > 0 and indxprocwork == 0:
            print 'Reading the previously computed proposal scale from %s...' % pathstdvprop
        gdatmodi.optipropdone = True
        varipara = pf.getdata(pathstdvprop)
    else:
        if gdat.verbtype > 0:
            print 'Optimizing proposal scale...'
        gdatmodi.optipropdone = False
    
    gdatmodi.stdvstdp[gdat.indxstdpcomp] = 0.
    gdatmodi.thislpritotltemp = copy(gdatmodi.thislpritotl)
    gdatmodi.thislliktotltemp = copy(gdatmodi.thislliktotl)
    gdatmodi.thislpostotltemp = copy(gdatmodi.thislpostotl)
    gdatmodi.thissampvarbtemp = copy(gdatmodi.thissampvarb)
    gdatmodi.thissamptemp = copy(gdatmodi.thissamp)
    gdatmodi.thisindxpntsfulltemp = deepcopy(gdatmodi.thisindxpntsfull)
    
    # temp
    deltparastep = 1e-5
    maxmstdv = 1.
    fudgstdv = 1.
    #diffpara = zeros((2, 2, 2))
    #diffpara[0, 0, :] = deltparastep * array([-1., -1.])
    #diffpara[0, 1, :] = deltparastep * array([-1., 1.])
    #diffpara[1, 0, :] = deltparastep * array([1., -1.])
    #diffpara[1, 1, :] = deltparastep * array([1., 1.])
    diffpara = deltparastep * array([-1., 0., 1])
    deltlpos = zeros(3)
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
            
            print 'gdat.namepara[k]'
            print gdat.namepara[k]
            
            for n in range(numbiter):
                if n != indxcntr:
                    deltlpos[n] = retr_deltlpos(gdat, gdatmodi, array([k]), array([diffpara[n]]))
                    gdatmodi.cntrswep += 1
            
            hess = 4. * fabs(deltlpos[0] + deltlpos[2] - 2. * deltlpos[1]) / deltparastep**2
            stdv = 1. / sqrt(hess)# / sqrt(gdat.numbpara) / fudgstdv
            
            #for n in gdat.indxpara:
            #    if n in gdat.indxfixpprop or n in concatenate(gdatmodi.thisindxsampcomp):
            #        if k == n:
            #            gdat.hess[k, n] = 1. / 4. / deltparastep * (deltlpos[1, 1] - deltlpos[1, 0] - deltlpos[0, 1] + deltlpos[0, 0])
            #            for a in range(2):
            #                
            #        else:
            #            f
            #            for a in range(2):
            #                for b in range(2):
            #                    deltlpos[a, b] = retr_deltlpos(gdat, gdatmodi, array([k, n]), diffpara[a, b, :])
            #                    print 'a'
            #                    print a
            #                    print 'b'
            #                    print b
            #                    print 'deltlpos[a, b]'
            #                    print deltlpos[a, b]
            #                    print

            if False or gdat.verbtype > 1:
                print 'deltlpos'
                print deltlpos
                print 'stdv'
                print stdv
                print
                print
                print
                print
                print
                print
                print
            
            if gdat.verbtype > 0 and not isfinite(stdv):
                print 'Proposal scale estimate went infinite.'
                print gdat.namepara[k]
                raise Exception('')
                print

            #if stdv > maxmstdv or not isfinite(stdv):
            #    stdv = maxmstdv
            
            if k in concatenate(gdatmodi.thisindxsampcomp):
                cntr = 0
                indxpnts = (k - gdat.indxsampcompinit)
                for strg in gdat.liststrgcomptotl:
                    if k in concatenate(getattr(gdatmodi, 'thisindxsamp' + strg)):
                        indxsampflux = k + 2 - cntr
                        gdatmodi.dictmodi['stdv' + strg + 'indv'].append(stdv)
                        if strg == 'flux':
                            gdatmodi.stdvstdp[gdat.indxstdppara[k]] += stdv * (gdatmodi.thissampvarb[indxsampflux] / gdat.minmflux)**2.
                        else:
                            gdatmodi.stdvstdp[gdat.indxstdppara[k]] += stdv * (gdatmodi.thissampvarb[indxsampflux] / gdat.minmflux)**0.5
                        
                        gdatmodi.dictmodi['stdv' + strg + 'indvflux'].append(gdatmodi.thissampvarb[indxsampflux])
                    cntr += 1
            else:
                gdatmodi.stdvstdp[gdat.indxstdppara[k]] = stdv
        
        if gdat.verbtype > 0:
            gdatmodi.cntrparasave = tdpy.util.show_prog(k, gdat.numbpara, gdatmodi.cntrparasave, indxprocwork=indxprocwork)

    for strg in gdat.liststrgcomptotl:
        gdatmodi.dictmodi['stdv' + strg + 'indv'] = array(gdatmodi.dictmodi['stdv' + strg + 'indv'])
        gdatmodi.dictmodi['stdv' + strg + 'indvflux'] = array(gdatmodi.dictmodi['stdv' + strg + 'indvflux'])
    
    for k in gdat.indxstdp:
        if k in gdat.indxstdpcomp:
            gdatmodi.stdvstdp[k] /= sum(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts])
            
            # temp 
            #gdatmodi.stdvstdp[k] = 1.
    
    gdatmodi.thissamp = copy(gdatmodi.thissamptemp)
    gdatmodi.thissampvarb = copy(gdatmodi.thissampvarbtemp)

    proc_samp(gdat, gdatmodi, 'this')
    
    if gdat.makeplot:
        
        if gdat.numbproc > 1:
            lock.acquire()
    
        xdat = gdat.indxstdp
        ydat = gdatmodi.stdvstdp
        path = gdat.pathopti + 'stdv%d.pdf' % indxprocwork
        tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist')
        
        for strg in gdat.liststrgcomptotl:
            path = gdat.pathopti + 'stdv' + strg + '.pdf'
            xdat = [gdatmodi.dictmodi['stdv' + strg + 'indvflux'] * gdat.fluxfactplot, gdat.meanfluxplot * gdat.fluxfactplot]
            if strg == 'flux':
                ydat = [gdatmodi.dictmodi['stdv' + strg + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strg)] / (gdat.meanfluxplot / gdat.minmflux)**2.]
            else:
                ydat = [gdatmodi.dictmodi['stdv' + strg + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strg)] / (gdat.meanfluxplot / gdat.minmflux)**0.5]
            lablxdat = gdat.lablfeattotl['flux']
            scalxdat = gdat.dictglob['scalfluxplot']
            limtxdat = array(gdat.dictglob['limtfluxplot']) * gdat.fluxfactplot
            tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                             lablydat=r'$\sigma_{%s}$' % gdat.lablfeat[strg], plottype=['scat', 'line'])
            #tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
            #                                 lablydat=r'$\sigma_{%s}$%s' % (gdat.lablfeat[strg], gdat.lablfeatunit[strg]), plottype=['scat', 'line'])
            
            tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                             lablydat=r'$\sigma_{%s}$' % gdat.lablfeat[strg], plottype=['scat', 'line'])

        if gdat.numbproc > 1:
            lock.release()
    
    gdatmodi.cntrswep = 0
    gdatmodi.optipropdone = True


def work(pathoutpthis, lock, indxprocwork):

    path = pathoutpthis + 'gdatinit'
    gdat = readfile(path) 
    
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
    gdatmodi.thissamp = rand(gdat.numbpara)
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
                if k in gdat.indxfixpnumbpnts and (gdat.truefixp[k] < gdat.minmnumbpnts[k] or gdat.truefixp[k] > gdat.maxmnumbpnts[k]):
                    gdatmodi.thissamp[gdat.indxfixp[k]] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[k] + 1))
                else:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, gdat.truefixp[k], k)
            if gdat.datatype == 'inpt':
                if k in gdat.indxfixpnumbpnts:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, gdat.truenumbpnts[k], k)
                elif k in gdat.indxfixppsfp:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, gdat.exprpsfp[k-gdat.indxfixppsfpinit], k)
                else:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = rand()
        gdatmodi.thissampvarb[k] = icdf_fixp(gdat, '', gdatmodi.thissamp[k], k)
                
    ## lists of occupied and empty transdimensional parameters
    gdatmodi.thisindxpntsfull = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            gdatmodi.thisindxpntsfull.append(range(gdatmodi.thissamp[gdat.indxfixpnumbpnts[l]].astype(int)))
    else:
        gdatmodi.thisindxpntsfull = []
    
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampflux, gdatmodi.thisindxsampsind, \
                             gdatmodi.thisindxsampcurv, gdatmodi.thisindxsampexpo, gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
    
    if gdat.numbtrap > 0:
        
        for l in gdat.indxpopl:
            gdatmodi.thissamp[gdatmodi.thisindxsampcomp[l]] = rand(gdatmodi.thisindxsampcomp[l].size)
        
        if gdat.strgcnfg == 'pcat_ferm_inpt_ngal':
            print 'gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]'
            print gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]
            print 'len(gdatmodi.thisindxpntsfull[0])'
            print len(gdatmodi.thisindxpntsfull[0])
            print 'gdat.truenumbpnts[0]'
            print gdat.truenumbpnts[0]
            print 'gdat.truesind[0]'
            print gdat.truesind[0].size
            print 'gdatmodi.thisindxsampsind[0]'
            print gdatmodi.thisindxsampsind[l].size
            print

        ## PS components
        if gdat.inittype == 'refr':
            for l in gdat.indxpopl:
                gdatmodi.thissamp[gdatmodi.thisindxsamplgal[l]] = cdfn_self(gdat.truelgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                gdatmodi.thissamp[gdatmodi.thisindxsampbgal[l]] = cdfn_self(gdat.truebgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                
                try:
                    indxtruepntsgood = where(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :] > gdat.minmflux)[0]
                    flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                    if gdat.fluxdisttype[l] == 'powr':
                        fluxunit = cdfn_powr(flux, gdat.minmflux, gdat.maxmflux, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]])
                    if gdat.fluxdisttype[l] == 'bind':
                        fluxdistnorm = gdatmodi.thissampvarb[gdat.indxfixpfluxdistnorm[l, :]]
                        fluxunit = cdfn_bind(flux, gdat.minmflux, gdat.maxmflux, gdat.binsflux, fluxdistnorm)
                    gdatmodi.thissamp[gdatmodi.thisindxsampflux[l][indxtruepntsgood]] = fluxunit[indxtruepntsgood]
                except:
                    gdatmodi.thissamp[gdatmodi.thisindxsampflux[l]] = rand(gdat.truenumbpnts[l])
                    
                if gdat.numbener > 1:
                    if gdat.sinddisttype[l] == 'gaus':
                        gdatmodi.thissamp[gdatmodi.thisindxsampsind[l]] = cdfn_gaus(gdat.truesind[l], gdatmodi.thissampvarb[gdat.indxfixpsinddistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpsinddiststdv[l]])
                    if gdat.sinddisttype[l] == 'atan':
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
            
    if gdat.verbtype > 1:
        print 'thissamp'
        for k in gdat.indxpara:
            if k in concatenate(gdatmodi.thisindxsamplgal):
                print
            print '%15f' % gdatmodi.thissamp[k]

    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.thissamp[gdat.numbpopl:] <= 0.)[0] + gdat.numbpopl
    indxsampbadduppr = where(gdatmodi.thissamp[gdat.numbpopl:] >= 1.)[0] + gdat.numbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'indxsampbadd'
        print indxsampbadd
        print gdat.namepara[indxsampbadd]
        print 'gdatmodi.thissamp'
        print gdatmodi.thissamp[:, None]
        raise Exception('Initial unit sample vector went outside the unit interval...')
    
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
   
    ## initial predicted count maps
    if gdat.pntstype == 'lght':
        # temp
        if gdat.boolintpanglcosi:
            binsangltemp = gdat.binsanglcosi
        else:
            binsangltemp = gdat.binsanglplot
        
    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    gdatmodi.listindxpntsfull = []
    gdatmodi.listspecassc = []
    
    gdatmodi.liststrgvarbswep = ['memoresi', 'lpri', 'lfctprop', 'lpriprop', 'lpau', 'deltlliktotl', 'lliktotl', 'chrototl', \
                                                                    'chrollik', 'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype']
    if gdat.probbrde < 1.:
        gdatmodi.liststrgvarbswep += ['auxipara', 'numbpair', 'jcbnfact', 'combfact']
   
    # process the initial sample, define the variables to be processed in each sample
    proc_samp(gdat, gdatmodi, 'this')

    # temp
    gdatmodi.thismemoresi = zeros(1)
    gdatmodi.thisdeltlliktotl = zeros(1)
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
    
    gdatmodi.liststrgvarbsamp = []
    for strg, valu in gdatmodi.__dict__.iteritems():
        if strg.startswith('this') and not strg[4:] in gdatmodi.liststrgvarbswep and isinstance(valu, ndarray):
            gdatmodi.liststrgvarbsamp.append(strg[4:])
    
    workdict = {}
    for strg in gdatmodi.liststrgvarbswep + gdatmodi.liststrgvarbsamp:
        valu = getattr(gdatmodi, 'this' + strg)
        if strg in gdatmodi.liststrgvarbswep:
            shap = [gdat.numbswep] + list(valu.shape)
        else:
            shap = [gdat.numbsamp] + list(valu.shape)
        workdict['list' + strg] = zeros(shap)
    
    ## saved state of the sample index used for logging progress status
    gdatmodi.percswepsave = -1
   
    # store the initial sample as the best fit sample
    gdatmodi.maxmllikswep = sum(gdatmodi.thisllik)
    gdatmodi.indxswepmaxmllik = -1 
    gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
   
    # proposal scale optimization
    gdatmodi.stdvstdp = copy(gdat.stdvstdp)
    gdatmodi.nextstdvstdp = copy(gdat.stdvstdp)
    gdatmodi.listllikopti = []
    gdatmodi.optillikdone = not gdat.optillik
    gdatmodi.optipropdone = not gdat.optiprop

    # log the initial state
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdatmodi, 'gdatmodi')

    if gdat.verbtype > 0:
        print 'Sampling...'

    if gdat.verbtype > 1:
        print 'gdatmodi.stdvstdp'
        print gdatmodi.stdvstdp
        print gdatmodi.stdvstdp[:, None]
        print

    while gdatmodi.cntrswep < gdat.numbswep:
        
        gdatmodi.thischrototl[:] = 0.
        #print 'gdatmodi.listllikopti'
        #print gdatmodi.listllikopti
        #print 'gdatmodi.listllikopti'
        #print gdatmodi.listllikopti
        #print 'gdatmodi.stdvstdp'
        #print gdatmodi.stdvstdp
        #print 
        if not gdatmodi.optipropdone:
            if not gdatmodi.optillikdone:
                booltemp = len(gdatmodi.listllikopti) > 10 and amax(diff(array(gdatmodi.listllikopti))[-10:]) < 1e-3
            else:
                booltemp = True
            if booltemp:
                optiprop(gdat, gdatmodi, indxprocwork)
            
            # temp
            for k in gdat.indxfixp:
                if not k in gdat.indxfixpnumbpnts:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = rand()
            for l in gdat.indxpopl:
                print 'hey'
                gdatmodi.thissamp[gdatmodi.thisindxsampcomp[l]] = rand(gdatmodi.thisindxsampcomp[l].size)

        if gdat.emptsamp:
            continue

        timetotlinit = gdat.functime()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                                                                          and gdat.makeplotfram and gdat.makeplot and gdatmodi.optipropdone
        # choose a proposal type
        timeinit = gdat.functime()
        
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
        
        timefinl = gdat.functime()
        gdatmodi.thischrototl[0] = timefinl - timeinit

        # propose the next sample
        timeinit = gdat.functime()
        
        retr_prop(gdat, gdatmodi)
        
        timefinl = gdat.functime()
        gdatmodi.thischrototl[1] = timefinl - timeinit
       
        # diagnostics
        if gdat.diagmode:
            
            timeinit = gdat.functime()
            
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
            
            if gdatmodi.propwith and gdatmodi.thislpautotl != 0.:
                raise Exception('Auxiliary variable PDF should is not zero during a within-model proposal.')

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
                if not isfinite(gdatmodi.thisdeltlliktotl):
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
        
            timefinl = gdat.functime()
            gdatmodi.thischrototl[2] = timefinl - timeinit
       
    
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
            gdatmodi.thischrototl[3] = timefinl - timeinit

        # plot the current sample
        if thismakefram:
            
            timeinit = gdat.functime()
            
            if gdat.verbtype > 0:
                print 'Process %d is in queue for making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                lock.acquire()
            if gdat.verbtype > 0:
                print 'Process %d started making a frame.' % indxprocwork
            
            proc_samp(gdat, gdatmodi, 'this')
            plot_samp(gdat, gdatmodi, 'this')
            if gdat.verbtype > 0:
                print 'Process %d finished making a frame.' % indxprocwork
        
            if gdat.numbproc > 1:
                lock.release()
        
            timefinl = gdat.functime()
            gdatmodi.thischrototl[4] = timefinl - timeinit
    
        # temp
        if False and gdat.pntstype == 'lght':
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
       
        # determine the acceptance probability
        timeinit = gdat.functime()

        gdatmodi.thisaccpprop = gdatmodi.thisaccpprio and gdatmodi.thisaccppsfn
        if gdatmodi.thisaccpprop:
            
            proc_samp(gdat, gdatmodi, 'next')
            gdatmodi.thisdeltlliktotl = gdatmodi.nextlliktotl - gdatmodi.thislliktotl
            
            # evaluate the acceptance probability
            accpprob = exp(gdatmodi.thisdeltlliktotl + gdatmodi.nextlpritotl - gdatmodi.thislpritotl + gdatmodi.thislpautotl + gdatmodi.thislfctprop + \
                                                                                                                                gdatmodi.thisjcbnfact + gdatmodi.thiscombfact)
            
        else:
            accpprob = 0.
    
        timefinl = gdat.functime()
        gdatmodi.thischrototl[5] = timefinl - timeinit
   
        # accept or reject the proposal
        timeinit = gdat.functime()

        if gdatmodi.optillikdone:
            booltemp = accpprob >= rand()
        else:
            booltemp = gdatmodi.thisaccpprop and gdatmodi.thisdeltlliktotl > 0.
        if booltemp:
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
        
        timefinl = gdat.functime()
        gdatmodi.thischrototl[6] = timefinl - timeinit
   
        # save the execution time for the sweep
        if not thismakefram:
            timetotlfinl = gdat.functime()
            gdatmodi.thischrototl[7] = timetotlfinl - timetotlinit
       
        ## variables to be saved for each sweep
        for strg in gdatmodi.liststrgvarbswep:
            workdict['list' + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi, 'this' + strg)
        
        # log the progress
        if gdat.verbtype > 0:
            gdatmodi.nextpercswep = 10 * int(10. * gdatmodi.cntrswep / gdat.numbswep) 
            if gdatmodi.nextpercswep > gdatmodi.percswepsave:
                gdatmodi.percswepsave = gdatmodi.nextpercswep
                
                minm = max(0, gdatmodi.cntrswep - 1000)
                maxm = gdatmodi.cntrswep + 1
                if maxm > minm:
                    fact = 100. / float(maxm - minm)
                    accp = fact * where(workdict['listaccp'][minm:maxm])[0].size
                    print '%3d%% completed.' % gdatmodi.nextpercswep
                    print 'Acceptance rate: %.3g%%' % accp
                    for k in gdat.indxproptype:
                        numb = where(workdict['listindxproptype'][minm:maxm] == k)[0].size
                        if numb > 0:
                            fact =  100. / float(where(workdict['listindxproptype'][minm:maxm] == k)[0].size)
                            accp = fact * where(logical_and(workdict['listaccp'][minm:maxm], workdict['listindxproptype'][minm:maxm] == k))[0].size
                            print '%s acceptance rate: %.3g%%' % (gdat.legdproptype[k], accp)
                        
                        if gdat.optipropsimp and k == gdat.indxproptypewith and accp < 0.1:
                            gdatmodi.cntrswep = 0
                            gdatmodi.percswepsave = -1.
                            gdatmodi.stdvstdp /= 2.
                            if gdat.verbtype > 0:
                                print 'Acceptance ratio went below 0.1%%.'
                                print 'Restarting the chain with a smaller proposal scale...'
                                
                    print 'Chronometers: '
                    for k in range(gdatmodi.thischrototl.size):
                        print '%.3g msec' % (gdatmodi.thischrototl[k] * 1e3)
                    print



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
        if gdatmodi.optillikdone:
            gdatmodi.cntrswep += 1
        else:
            if gdatmodi.thisaccp:
                gdatmodi.nextstdvstdp = gdatmodi.stdvstdp[k] * 2**(0.5 * randn(gdat.numbpara))
                gdatmodi.listllikopti.append(gdatmodi.nextlliktotl)
                gdatmodi.stdvstdp = copy(gdatmodi.nextstdvstdp)
            else:
                gdatmodi.nextstdvstdp = gdatmodi.stdvstdp[k] / 1.5**(0.5 * randn(gdat.numbpara))
  
    for strg in gdatmodi.liststrgvarbsamp + gdatmodi.liststrgvarbswep:
        valu = workdict['list' + strg]
        setattr(gdatmodi, 'list' + strg, valu)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    path = gdat.pathoutpthis + 'gdatmodi%04d' % indxprocwork
    writfile(gdatmodi, path) 
    #writoutp(gdatmodi, path, catl=False) 
    
