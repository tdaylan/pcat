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
         # temp
         diagmode=True, \
         emptsamp=False, \

         # sampler
         numbswep=None, \
         numbburn=None, \
         factthin=None, \
        
         condcatl=True, \
         seedstat=None, \
         indxevttincl=None, \
         indxenerincl=None, \
       
         # evaluate the likelihood inside circles around elements
         evalcirc=None, \
        
         numbspatdims=2, \
         
         # model type
         elemtype='lght', \
        
         # initialization type
         inittype=None, \
         
         procrtag=None, \
         loadvaripara=False, \
         
         propcova=True, \
         propwithsing=False, \
         optihess=False, \
         optillik=False, \
         optiprop=False, \
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
         numbswepplot=2000, \
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
         specfraceval=None, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=400, \
    
         listspecrefrplot=None, \
         listenerreftplot=None, \
         listlablreftplot=None, \

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
        
         checprio=False, \

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
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
         
         radispmr=None, \

         numbsidecart=100, \
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
            if gdat.elemtype == 'lens':
                gdat.inittype = 'pert'
            if gdat.elemtype == 'lght':
                gdat.inittype = 'refr'
        else:
            gdat.inittype = 'rand'

    if gdat.elemtype == 'lens':
        gdat.hubbexpofact = 1.63050e-19
    
    if gdat.strgexpo == None:
        if gdat.elemtype == 'lens':
            gdat.strgexpo = 1e3 / gdat.hubbexpofact
        else:
            gdat.strgexpo = 1.

    if gdat.lablflux == None:
        if gdat.elemtype == 'lens':
            gdat.lablflux = r'\alpha_s'
        else:
            if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
                gdat.lablflux = 'f'
            if gdat.exprtype == 'sdyn':
                gdat.lablflux = 'p'
    
    if gdat.lablfluxunit == None:
        if gdat.elemtype == 'lens':
            gdat.lablfluxunit = u'$^{\prime\prime}$'
        else:
            if gdat.exprtype == 'sdss' or gdat.exprtype == 'hubb':
                gdat.lablfluxunit = 'mag'
            if gdat.exprtype == 'ferm':
                gdat.lablfluxunit = '1/cm$^2$/s/GeV'
            if gdat.exprtype == 'chan':
                gdat.lablfluxunit = '1/cm$^2$/s/KeV'
            if gdat.exprtype == 'sdyn':
                gdat.lablfluxunit = ''

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
        if gdat.elemtype == 'lens':
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
        if gdat.elemtype == 'lens':
            gdat.evalcirc = 'full'
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
    if gdat.exprtype == 'hubb' or gdat.exprtype == 'sdyn':
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
    gdat.numbener = gdat.indxenerincl.size
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.indxenerfluxdist = ceil(array([gdat.numbener]) / 2.).astype(int) - 1
    if gdat.enerbins:
        gdat.numbenerplot = 100
        gdat.numbenerfull = len(gdat.strgenerfull)
        gdat.strgener = [gdat.strgenerfull[k] for k in gdat.indxenerincl]
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
        retr_sdynpsfn(gdat)
   
    # number of processes
    gdat.strgproc = os.uname()[1]
    if gdat.numbproc == None:
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu' or gdat.strgproc == 'wise':
            gdat.numbproc = 10
        else:
            gdat.numbproc = 1
    
    # number of total sweeps
    if gdat.numbswep == None:
        gdat.numbswep = 100000

    # number of burned sweeps
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 10

    # factor by which to thin the sweeps to get samples
    if gdat.factthin == None:
        gdat.factthin = int(ceil(1e-3 * (gdat.numbswep - gdat.numbburn) * gdat.numbproc))
    
    # model evaluation approximation error tolerance in units of the fraction of the lowest PS flux
    if gdat.specfraceval == None:
        if gdat.exprtype == 'ferm':
            gdat.specfraceval = 0.5
        else:
            gdat.specfraceval = 0.1

    if gdat.strgcatl == None:
        if gdat.exprtype == 'ferm':
            gdat.strgcatl = '3FGL'
        else:
            gdat.strgcatl = gdat.strgexprname
    
    # evaluation of the PSF
    if gdat.elemtype == 'lens':
        gdat.evalpsfnpnts = False
    else:
        gdat.evalpsfnpnts = True

    ## generative model
    # set mock sample vector indices
    setp_true(gdat, 'numbpnts', array([50]))
    gdat.truenumbpopl = gdat.truenumbpnts.size
    gdat.trueindxpopl = arange(gdat.truenumbpopl)
    
    setp_true(gdat, 'minmnumbpnts', zeros(gdat.truenumbpopl, dtype=int) + 1)
    setp_true(gdat, 'maxmnumbpnts', zeros(gdat.truenumbpopl, dtype=int) + 1000)
    
    ### element parameter distributions
    setp_true(gdat, 'spatdisttype', ['unif' for l in gdat.trueindxpopl])
    setp_true(gdat, 'fluxdisttype', ['powr' for l in gdat.trueindxpopl])
    setp_true(gdat, 'spectype', ['powr' for l in gdat.trueindxpopl])
    
    gdat.calcllik = True
    
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
    
    if gdat.psfninfoprio == None:
        if gdat.exprtype == 'ferm':
            gdat.psfninfoprio = True
        else:
            gdat.psfninfoprio = False

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
    if gdat.elemtype == 'lght':
        if gdat.exprtype == 'ferm':
            back = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
        elif gdat.exprtype == 'chan':
            back = [array([1.3e1, 4.]), 1e1]
        else:
            back = [1.]
    if gdat.elemtype == 'lens':
        back = [1e-7]
    setp_true(gdat, 'back', back)
    
    #### spectra for the background
    if gdat.exprtype == 'chan':
        specback = [None, array([1., 1.3])]
    else:
        specback = [None, None]
    setp_true(gdat, 'specback', specback)
    
    #### label
    lablback = [r'$\mathcal{I}$']
    if gdat.exprtype == 'ferm':
        lablback.append(r'$\mathcal{D}$')
    if gdat.exprtype == 'chan':
        lablback.append(r'$\mathcal{I}_p$')
    setp_true(gdat, 'lablback', lablback)
    
    #### name
    strgback = ['fluxisot']
    if gdat.exprtype == 'ferm':
        strgback.append('fluxfdfm')
    if gdat.exprtype == 'chan':
        strgback.append('fluxpart')
    setp_true(gdat, 'strgback', strgback)
    
    #### legend
    nameback = ['Isotropic']
    if gdat.exprtype == 'ferm':
        nameback.append(r'FDM')
    if gdat.exprtype == 'chan':
        nameback.append(r'Particle')
    setp_true(gdat, 'nameback', nameback)
    
    retr_indxsamp(gdat, strgpara='true')

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
        if gdat.exprtype == 'sdyn':
            minmsigm = 0.1
            maxmsigm = 0.2
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
        minmflux = 3e-11
    if gdat.exprtype == 'chan':
        minmflux = 3e-9
    if gdat.exprtype == 'sdyn':
        minmflux = 1e0
    if gdat.elemtype == 'lens':
        minmflux = 2e-3 / gdat.anglfact
    setp_true(gdat, 'minmflux', minmflux)
    
    if gdat.exprtype == 'ferm':
        maxmflux = 1e-7
    if gdat.exprtype == 'chan':
        maxmflux = 1e-7
    if gdat.exprtype == 'sdyn':
        maxmflux = 1e4
    if gdat.elemtype == 'lens':
        maxmflux = 5e-2 / gdat.anglfact
    setp_true(gdat, 'maxmflux', maxmflux)
   
    # parameter defaults
    ## distribution
    ### flux
    setp_truedefa(gdat, 'gangdistscal', [1. / gdat.anglfact, 10. / gdat.anglfact], popl=True)
    setp_truedefa(gdat, 'spatdistcons', [1e-4, 1e-2], popl=True)
    setp_truedefa(gdat, 'bgaldistscal', [0.5 / gdat.anglfact, 5. / gdat.anglfact], popl=True)
    setp_truedefa(gdat, 'fluxdistslop', [1., 4.], popl=True)
    
    for k in range(gdat.numbfluxdistnorm):
        setp_truedefa(gdat, 'fluxdistnormbin%d' % k, [1e-3, 1e3], popl=True)

    ### spectral index
    if gdat.elemtype == 'lght' and gdat.numbener > 1:
        if gdat.exprtype == 'ferm':
            sind = [1., 3.]
        if gdat.exprtype == 'chan':
            sind = [-1., 3.]
        setp_truedefa(gdat, 'sind', sind, popl=True)
        setp_truedefa(gdat, 'sinddistmean', sind, popl=True)
        #### standard deviations should not be too small
        setp_truedefa(gdat, 'sinddiststdv', [0.3, 2.], popl=True)
        setp_truedefa(gdat, 'curvdistmean', [-1., 1.], popl=True)
        setp_truedefa(gdat, 'curvdiststdv', [0.1, 1.], popl=True)
        setp_truedefa(gdat, 'expodistmean', [1., 8.], popl=True)
        setp_truedefa(gdat, 'expodiststdv', [0.01 * gdat.maxmener, gdat.maxmener], popl=True)
    
    ### maximum horizontal/vertical distance of the elements from the image center
    gdat.maxmgang = gdat.maxmgangdata * gdat.margfactmodl
    gdat.truemaxmgang = gdat.maxmgang
    
    if gdat.elemtype == 'lens':
        ### scale angular distance
        setp_truedefa(gdat, 'ascadistmean', [gdat.maxmgang * 1e-2, gdat.maxmgang], popl=True)
        setp_truedefa(gdat, 'ascadiststdv', [gdat.maxmgang * 1e-3, gdat.maxmgang], popl=True)
        ### cutoff angular distance
        setp_truedefa(gdat, 'acutdistmean', [gdat.maxmgang * 1e-2, gdat.maxmgang], popl=True)
        setp_truedefa(gdat, 'acutdiststdv', [gdat.maxmgang * 1e-3, gdat.maxmgang], popl=True)
    
    ## lensing
    setp_truedefa(gdat, 'lgalsour', [-gdat.maxmgang, gdat.maxmgang])
    setp_truedefa(gdat, 'bgalsour', [-gdat.maxmgang, gdat.maxmgang])
    setp_truedefa(gdat, 'lgalhost', [-gdat.maxmgang, gdat.maxmgang])
    setp_truedefa(gdat, 'bgalhost', [-gdat.maxmgang, gdat.maxmgang])
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
    
    gdat.trueminmlgal = -gdat.maxmgang
    gdat.truemaxmlgal = gdat.maxmgang
    gdat.trueminmbgal = -gdat.maxmgang
    gdat.truemaxmbgal = gdat.maxmgang
    
    gdat.trueminmaang = -pi
    gdat.truemaxmaang = pi
    
    # copy the true model to the inference model if the inference model parameter has not been specified
    temp = deepcopy(gdat.__dict__)
    for strg, valu in temp.iteritems():
        if strg.startswith('true') and not strg[4:].startswith('indx'):
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
    if gdat.datatype == 'inpt':
        if 'gaus' in gdat.spatdisttype and (gdat.lgalprio == None or gdat.bgalprio == None):
            raise Exception('If spatdisttype is "gaus", spatial coordinates of the prior catalog should be provided via lgalprio and bgalprio.')

    if gdat.verbtype > 0:
        print 'PCAT v0.1 started at %s' % gdat.strgtimestmp
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)

    setp_true(gdat, 'spatdistcons', 1e-3, popl=True)
    setp_true(gdat, 'gangdistscal', 4. / gdat.anglfact, popl=True)
    setp_true(gdat, 'bgaldistscal', 2. / gdat.anglfact, popl=True)
    if gdat.elemtype == 'lens':
        fluxdistslop = 3.
    else:
        fluxdistslop = 2.2
    setp_true(gdat, 'fluxdistslop', fluxdistslop, popl=True)

    retr_axis(gdat, 'flux', gdat.trueminmflux, gdat.maxmflux, gdat.numbfluxdistnorm - 1, scal='logt', strgtype='true')
    fluxdistnorm = gdat.truebinsflux**(-fluxdistslop)
    fluxdistnorm *= 1e-2 / amin(fluxdistnorm)
    for k in range(gdat.numbfluxdistnorm):
        setp_true(gdat, 'fluxdistnormbin%d' % k, fluxdistnorm[k], popl=True)
    
    setp_true(gdat, 'ascadistmean', 0.2 / gdat.anglfact, popl=True)
    setp_true(gdat, 'ascadiststdv', 0.04 / gdat.anglfact, popl=True)
    setp_true(gdat, 'acutdistmean', 0.6 / gdat.anglfact, popl=True)
    setp_true(gdat, 'acutdiststdv', 0.04 / gdat.anglfact, popl=True)
        
    if gdat.numbener > 1:
        if gdat.exprtype == 'ferm':
            sinddistmean = 2.15
        if gdat.exprtype == 'chan':
            sinddistmean = 1.
        setp_true(gdat, 'sinddistmean', sinddistmean, popl=True)
        setp_true(gdat, 'sinddiststdv', 0.5, popl=True)
        
        setp_true(gdat, 'curvdistmean', 2., popl=True)
        setp_true(gdat, 'curvdiststdv', 0.2, popl=True)
        
        setp_true(gdat, 'expodistmeanpop1', 2., popl=True)
        setp_true(gdat, 'expodiststdvpop1', 0.2, popl=True)
    
    setp_true(gdat, 'bacp', 1., ener=True, back=True)
    
    if gdat.elemtype == 'lens':
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
    
    if gdat.listspecrefrplot == None:
        if gdat.exprtype == 'chan':
            gdat.listenerrefrplot = []
            gdat.listspecrefrplot = []
            gdat.listlablrefrplot = ['CDFS', 'HEAO', 'XMM', 'ROSAT']
            for strgfile in ['cdfs', 'heao', 'xmmm', 'rost']:
                path = gdat.pathinpt + '%s.csv' % strgfile
                gdat.listenerrefrplot.append(loadtxt(path, delimiter=',')[:, 0])
                gdat.listspecrefrplot.append(loadtxt(path, delimiter=',')[:, 1])

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
    gdat.maxmllik = 0.
    gdat.minmllik = -1e2
    gdat.scalllik = 'linr'
    gdat.cmapllik = 'YlGn'
    
    gdat.scaldatacnts = 'logt'
    gdat.cmapdatacnts = 'Greys'
    
    gdat.scalresicnts = 'asnh'
    gdat.cmapresicnts = 'RdBu'
    
    gdat.scalexpo = 'logt'
    gdat.cmapexpo = 'OrRd'
    
    gdat.minmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) - 2.
    gdat.maxmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) + 2.
    gdat.scallpdfspatpriointp = 'linr'
    gdat.cmaplpdfspatpriointp = 'PuBu'
    
    gdat.minmconv = 1e-4
    gdat.maxmconv = 10.
    gdat.scalconv = 'logt'
    gdat.cmapconv = 'Purples'
    
    gdat.minminvm = -1e2
    gdat.maxminvm = 1e2
    gdat.scalinvm = 'asnh'
    gdat.cmapinvm = 'BrBG'
    
    gdat.minmdeflcomp = 1e-7
    gdat.maxmdeflcomp = 1e-4
    gdat.scaldeflcomp = 'logt'
    gdat.cmapdeflcomp = 'Oranges'
    
    liststrgcbar = ['llik', 'lpdfspatpriointp', 'conv', 'invm', 'deflcomp']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)
    
    # element features
    #gdat.minmgangplot = 0.
    #gdat.maxmgangplot = gdat.maxmlgal
    
    gdat.minmdeltllik = -1.
    gdat.maxmdeltllik = 4.
    gdat.minmdiss = 0.
    gdat.maxmdiss = 3. * gdat.maxmgang
    gdat.minmdots = 2e-7
    gdat.maxmdots = 1e-5

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

    for strgtype in ['', 'true']:
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
                scal = gdat.dictglob['scal' + strgfeat + 'plot']
                maxm = getattr(gdat, 'maxm' + strgfeat + 'plot')
                minm = getattr(gdat, 'minm' + strgfeat + 'plot')
                retr_axis(gdat, strgfeat + 'plot', minm, maxm, gdat.numbbinsplot, scal=scal, strgtype=strgtype)
                if strgfeat in gdat.liststrgfeatpriototl:
                    maxm = getattr(gdat, strgtype + 'maxm' + strgfeat)
                    minm = getattr(gdat, strgtype + 'minm' + strgfeat)
                    retr_axis(gdat, strgfeat, minm, maxm, gdat.numbbinsplotprio, scal=scal, strgtype=strgtype)
        
            gdat.dictglob['limt' + strgfeat + 'plot'] = array([getattr(gdat, 'minm' + strgfeat + 'plot'), getattr(gdat, 'maxm' + strgfeat + 'plot')])

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

    # experimental information
    gdat.exprinfo = False
    if gdat.exprtype == 'ferm':
        retr_fermdata(gdat)
        gdat.exprinfo = True
    if gdat.exprtype == 'chan':
        retr_chandata(gdat)
        gdat.exprinfo = True
    
    gdat.trueinfo = gdat.exprinfo or gdat.datatype == 'mock'

    # rotate element coordinates to the ROI center
    if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
        rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='ZYX')
        gdat.exprbgal, gdat.exprlgal = rttr(pi / 2. - gdat.exprbgal, gdat.exprlgal)
        gdat.exprbgal = pi / 2. - gdat.exprbgal

    # external reference catalog
    if gdat.exprinfo:
        
        # find entries inside the ROI
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
        
        # reorder elements with respect to flux
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
                print 'Experimental information on element counts is inconsistent.'
       
        gdat.exprflux = gdat.exprspec[0, gdat.indxenerfluxdist[0], :]
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
        if gdat.exprinfo:
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
   
    for strgfeat in gdat.liststrgfeat:
        try:
            hist = histogram(getattr(gdat, 'expr' + strgfeat), getattr(gdat, 'bins' + strgfeat + 'plot'))[0]
            setattr(gdat, 'expr' + strgfeat + 'hist', hist)
        except:
            pass

    if gdat.datatype == 'mock':
        
        gdat.trueindxpntsfull = []
        if gdat.numbtrap > 0:
            for l in gdat.trueindxpopl:
                gdat.trueindxpntsfull.append(range(gdat.truenumbpnts[l]))
        else:
            gdat.trueindxpntsfull = []

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
                    gdat.truefixp[k] = getattr(gdat, 'true' + gdat.truenamefixp[k])
                except:
                    pass
                
                # randomly sample the rest of the mock model parameters
                if not isfinite(gdat.truefixp[k]):
                    gdat.truefixp[k] = icdf_fixp(gdat, 'true', rand(), k)

        # temp
        if gdat.elemtype == 'lens':
            gdat.truefixp[gdat.trueindxfixplgalsour] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalsour] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixplgalhost] = 0.04 * gdat.maxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalhost] = 0.04 * gdat.maxmgang * randn()
        
        gdat.truesampvarb = empty(gdat.truenumbpara)
        gdat.truesampvarb[gdat.trueindxfixp] = gdat.truefixp
        
        gdat.truenumbpntstotl = sum(gdat.truefixp[gdat.trueindxfixpnumbpnts])
        gdat.trueindxpntstotl = arange(gdat.truenumbpntstotl)
        
        for strgfeat in gdat.liststrgfeat:
            setattr(gdat, 'true' + strgfeat, [[] for l in gdat.trueindxpopl])
        
        # temp
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
                
                gdat.truespec[l] = empty((3, gdat.numbener, gdat.truenumbpnts[l]))

                if gdat.truefluxdisttype[l] == 'powr':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_flux_powr(rand(gdat.truenumbpnts[l]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                                                    gdat.truefixp[gdat.trueindxfixpfluxdistslop[l]])
                if gdat.truefluxdisttype[l] == 'bind':
                    gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = icdf_bind(rand(gdat.truenumbpnts[l]), gdat.trueminmflux, gdat.maxmflux, \
                                                                                           gdat.truebinsflux, gdat.truefixp[gdat.trueindxfixpfluxdistnorm[l, :]])
                
                # temp -- make sure this reordering does not mess up other things
                gdat.truespec[l][:, gdat.indxenerfluxdist[0], :] = sort(gdat.truespec[l][:, gdat.indxenerfluxdist[0], :], axis=1)[::-1]
   
        
                if gdat.elemtype == 'lens':
                    
                    gdat.trueasca[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpascadistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpascadiststdv[l]])
                    
                    gdat.trueacut[l] = icdf_gaus(rand(gdat.truenumbpnts[l]), gdat.truefixp[gdat.trueindxfixpacutdistmean[l]], \
                                                                                                                            gdat.truefixp[gdat.trueindxfixpacutdiststdv[l]])

                if gdat.elemtype == 'lght' and gdat.numbener > 1:
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
                    
                gdat.truesampvarb[gdat.trueindxsampcomp['lgal'][l]] = gdat.truelgal[l]
                gdat.truesampvarb[gdat.trueindxsampcomp['bgal'][l]] = gdat.truebgal[l]
                gdat.truesampvarb[gdat.trueindxsampcomp['flux'][l]] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                if gdat.elemtype == 'lens':
                    gdat.truesampvarb[gdat.trueindxsampcomp['asca'][l]] = gdat.trueasca[l]
                    gdat.truesampvarb[gdat.trueindxsampcomp['acut'][l]] = gdat.trueacut[l]
                if gdat.elemtype == 'lght' and gdat.numbener > 1:
                    gdat.truesampvarb[gdat.trueindxsampcomp['sind'][l]] = gdat.truesind[l]
                    if gdat.truespectype[l] == 'curv':
                        gdat.truesampvarb[gdat.trueindxsampcomp['curv'][l]] = gdat.truecurv[l]
                    if gdat.truespectype[l] == 'expo':
                        gdat.truesampvarb[gdat.trueindxsampcomp['expo'][l]] = gdat.trueexpo[l]
                
                    indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
                    gdat.truecnts[l] = gdat.truespec[l][0, :, :, None] * gdat.expo[:, indxpixltemp, :]
                    if gdat.enerdiff:
                        gdat.truecnts[l] *= gdat.deltener[:, None, None]
        
        if gdat.verbtype > 1:
            print
            print '%20s %20s %15s' % ('truenamepara', 'truesampvarb', 'truescalpara')
            for k in gdat.trueindxpara:
                if k in concatenate(gdat.trueindxsampcomp['lgal']):
                    print
                print '%20s %20f %15s' % (gdat.truenamepara[k], gdat.truesampvarb[k], gdat.truescalpara[k])
    
    if 'gaus' in gdat.spatdisttype:
        if gdat.lgalprio == None or gdat.bgalprio == None:
            gdat.lgalprio = concatenate((gdat.truelgal))
            gdat.bgalprio = concatenate((gdat.truebgal))
        gdat.numbspatprio = gdat.lgalprio.size
    
        # spatial template for the catalog prior
        # temp -- this should move outside the if
        gdat.apixmodl = (gdat.maxmgang / gdat.numbsidecart)**2
        gdat.pdfnspatpriotemp = zeros((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
        for k in range(gdat.numbspatprio):
            gdat.pdfnspatpriotemp[:] += 1. / sqrt(2. * pi) / gdat.stdvspatprio * exp(-0.5 * (gdat.binslgalcartmesh - gdat.lgalprio[k])**2 / gdat.stdvspatprio**2) * \
                                                                                        exp(-0.5 * (gdat.binsbgalcartmesh - gdat.bgalprio[k])**2 / gdat.stdvspatprio**2)
        gdat.pdfnspatpriotemp /= amax(gdat.pdfnspatpriotemp)
   
    # temp
    for l in gdat.trueindxpopl:
        setattr(gdat, 'truenumbpntspop%d' % l, gdat.truenumbpnts[l])
        setattr(gdat, 'truemeanpntspop%d' % l, gdat.truenumbpnts[l])
    
    if gdat.exprinfo:
        # true flux
        gdat.trueflux = [[] for l in gdat.trueindxpopl]
        for l in gdat.trueindxpopl:
            gdat.trueflux[l] = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
        
    if gdat.datatype == 'mock':
        
        proc_samp(gdat, None, 'true', raww=True)
        proc_samp(gdat, None, 'true')
            
        if gdat.makeplot:
            plot_samp(gdat, None, 'true')
        
        if gdat.seedstat != None:
            seed()
    
    if gdat.datatype == 'inpt':
        retr_datatick(gdat)
    
    if gdat.makeplot:
        plot_init(gdat)

    # final setup
    setpfinl(gdat, True) 
    
    # write the list of arguments to file
    fram = inspect.currentframe()
    listargs, temp, temp, listargsvals = inspect.getargvalues(fram)
    fileargs = open(gdat.pathoutpthis + 'args.txt', 'w')
    fileargs.write('PCAT call arguments\n')
    for args in listargs:
        fileargs.write('%s = %s\n' % (args, listargsvals[args]))
    fileargs.close()
    
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
    lock = mp.Manager().Lock()
        
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdat, 'gdat')
    
    # list of variables for which the posterior is calculated at each sweep
    gdat.liststrgvarbarryswep = ['memoresi', 'lpri', 'lfctprop', 'lpriprop', 'lpau', 'deltlliktotl', 'lliktotl', 'chrototl', \
                                                                    'chroproc', 'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype']
    if gdat.probbrde < 1.:
        gdat.liststrgvarbarryswep += ['auxipara', 'numbpair', 'jcbnfact', 'combfact']
    
    # perform a fudicial processing of a sample vector in order to find the list of variables for which the posterior will be calculated
    if gdat.verbtype > 0:
        print 'Processing a fudicial sample...'
    gdatmodifudi = tdpy.util.gdatstrt()
    gdatmodifudi.thischroproc = zeros(gdat.numbchroproc)
    gdatmodifudi.thissamp = rand(gdat.numbpara)
    gdatmodifudi.thisindxpntsfull = [[] for l in gdat.indxpopl]
    gdatmodifudi.thissamp[gdat.indxfixpnumbpnts[0]] = 1
    gdatmodifudi.thisindxpntsfull[0].append(0)
    gdatmodifudi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodifudi.thisindxpntsfull, gdat.spectype)
    gdatmodifudi.thissampvarb = retr_sampvarb(gdat, gdatmodifudi.thisindxsampcomp, gdatmodifudi.thissamp, 'this')
    proc_samp(gdat, gdatmodifudi, 'this')
    gdat.liststrgvarbarrysamp = []
    gdat.liststrgvarblistsamp = []
    for strg, valu in gdatmodifudi.__dict__.iteritems():
        if strg.startswith('this') and not strg[4:] in gdat.liststrgvarbarryswep and isinstance(valu, ndarray):
            gdat.liststrgvarbarrysamp.append(strg[4:])
        if strg.startswith('this') and isinstance(valu, list) and strg != 'thisindxsampcomp' and strg != 'thispsfnkern':
            gdat.liststrgvarblistsamp.append(strg[4:])
    
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
   
    setp_indxswepsave(gdat)
    
    gdat.opti = gdat.optiprop or gdat.optillik or gdat.optihess

    # estimate the covariance
    if gdat.opti:
        gdat.numbprocsave = gdat.numbproc
        gdat.numbproc = 1
    
        if gdat.verbtype > 0:
            print 'Optimizing proposal scale...'
    
        worksamp(gdat, lock)

        if gdat.verbtype > 0:
            print 'Writing the estimated covariance matrix to the disc...'
        
        thisfile = h5py.File(gdat.pathoutpthis + 'opti.h5', 'r')
        gdat.stdvstdp = thisfile['stdvstdp'][()]
        thisfile.close()
        gdat.numbproc = gdat.numbprocsave
        
    gdat.optiprop = False
    gdat.optillik = False
    gdat.optihess = False
    gdat.opti = False

    # perform an initial run, sampling from the prior
    if gdat.checprio:
        
        if gdat.verbtype > 0:
            print 'Sampling from the prior...'
        
        gdat.calcllik = False
        
        # save some variables that will be changed for the prior-only run
        liststrgvarbsave = ['legdsampdist', 'makeplot', 'pathfram', 'pathpost', 'pathdiag', 'pathanim', 'stdvstdp']
        
        for strgvarbsave in liststrgvarbsave:
            varb = getattr(gdat, strgvarbsave)
            setattr(gdat, strgvarbsave + 'saveprio', varb)
       
        ## change the variables
        gdat.stdvstdp = 1e-2
        gdat.legdsampdist = 'Prior'
        gdat.pathpost += 'chec/'
        gdat.pathfram += 'chec/'
        gdat.pathdiag += 'chec/'
        gdat.pathanim += 'chec/'
        os.system('mkdir -p %s %s %s %s' % (gdat.pathpost, gdat.pathfram, gdat.pathdiag, gdat.pathanim))

        ## perform sampling
        worksamp(gdat, lock)
        
        ## post process the samples
        proc_post(gdat, prio=True)
        
        # save the prior median of variables
        for strgvarb in gdat.liststrgchan:
            setattr(gdat, 'prio' + strgvarb, getattr(gdat, 'medi' + strgvarb))
        
        gdat.calcllik = True
        
        ## retrieve saved variables
        for strgvarbsave in liststrgvarbsave:
            varb = getattr(gdat, strgvarbsave + 'saveprio')
            setattr(gdat, strgvarbsave, varb)

    else:
        gdat.calcllik = True
    
    # run the sampler
    worksamp(gdat, lock)
    
    # post process the samples
    proc_post(gdat)

    if gdat.verbtype > 0:
        print 'The ensemble of catalogs is at ' + pathcatl
        if gdat.makeplot:
            print 'The plots are at ' + gdat.pathplot
        print 'PCAT has run successfully. Returning to the OS...'
    
    return gdat
    

def proc_post(gdat, prio=False):

    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        timeinit = gdat.functime()
   
    # read the chains
    listgdatmodi = []
    for k in gdat.indxproc:
        path = gdat.pathoutpthis + 'gdatmodi%04d' % k
        listgdatmodi.append(readfile(path))
        os.system('rm -rf %s*' % path)
    os.system('rm -rf %s' % gdat.pathoutpthis + 'gdatinit*')

    # aggregate samples from the chains
    ## list of parameters to be gathered
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
    gdat.liststrgchan = gdat.liststrgvarbarrysamp + ['fixp']
    # temp
    #gdat.liststrgvarblistsamp + ['fixp']
    if gdat.trueinfo:
        gdat.liststrgchan += ['specassc']

    for strg in gdat.liststrgvarbarry:
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
    ## element indices
    for strgvarb in gdat.liststrgvarblistsamp:
        setattr(gdat, 'list' + strgvarb, [])

    # temp
    #gdat.listindxpntsfull = []
    #gdat.listspecassc = []
    for strgvarb in gdat.liststrgvarblistsamp:
        listtemp = []
        for j in gdat.indxsamp:      
            for k in gdat.indxproc:
                listtemp.append(getattr(listgdatmodi[k], 'list' + strgvarb)[j])
                #gdat.listindxpntsfull.append(listgdatmodi[k].listindxpntsfull[j])
        setattr(gdat, 'list' + strgvarb, listtemp)
                #if gdat.trueinfo:
                #    gdat.listspecassc.append(listgdatmodi[k].listspecassc[j])

    ## list of other parameters to be flattened
    gdat.liststrgvarbarryflat = deepcopy(gdat.liststrgvarbarry)
    for strg in ['deltlliktotl', 'memoresi']:
        gdat.liststrgvarbarryflat.remove(strg)
   
    gdat.listsampvarbproc = copy(gdat.listsampvarb)

    ## other parameters
    for strg in gdat.liststrgvarbarryflat:
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
    gdat.info = retr_infofromlevi(gdat.listlliktotl, gdat.levi)
    
    # parse the sample vector
    gdat.listfixp = gdat.listsampvarb[:, gdat.indxfixp]
    for k, namefixp in enumerate(gdat.namefixp):
        setattr(gdat, 'list' + namefixp, gdat.listfixp[:, k])
    
    # temp
    if gdat.elemtype == 'lens':
        gdat.listfracsubh = gdat.listfracsubh.flatten()
    
    if gdat.checprio:
        for namevarbscal in gdat.listnamevarbscal:
            binsvarbscal = getattr(gdat, 'bins' + namevarbscal)
            meanvarbscal = getattr(gdat, 'mean' + namevarbscal)
            deltvarbscal = getattr(gdat, 'delt' + namevarbscal)
            listvarbscal = getattr(gdat, 'list' + namevarbscal)
            pdfn = histogram(listvarbscal, bins=binsvarbscal)[0].astype(float)
            #try:
            #    pdfn = sp.stats.gaussian_kde(listvarbscal)(meanvarbscal)
            #except: 
            #    pdfn = zeros_like(meanvarbscal)
            
            pdfn[where(pdfn < 1e-50)[0]] = 1e-50

            if prio:
                strgtemp = 'prio'
            else:
                strgtemp = 'post'
            
            setattr(gdat, 'pdfn' + strgtemp + namevarbscal, pdfn)
        
        if not prio:
            for namevarbscal in gdat.listnamevarbscal:
                pdfnpost = getattr(gdat, 'pdfnpost' + namevarbscal)
                pdfnprio = getattr(gdat, 'pdfnprio' + namevarbscal)
                info = retr_info(pdfnpost, pdfnprio)
                deltvarbscal = getattr(gdat, 'delt' + namevarbscal)
                infototl = sum(info * deltvarbscal)
                setattr(gdat, 'info' + namevarbscal, info)
                setattr(gdat, 'infototl' + namevarbscal, infototl)
                
    # post process samples
    ## bin element features
    if gdat.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        timeinit = gdat.functime()
       
    gdat.pntsprob = [[] for l in gdat.indxpopl]
    for l in gdat.indxpopl:
        numb = len(gdat.liststrgfeatsign[l])
        gdat.pntsprob[l] = zeros((gdat.numbbgalpntsprob, gdat.numblgalpntsprob, gdat.numbbinsplot, numb))
        temparry = concatenate([gdat.listlgal[n][l] for n in gdat.indxsamptotl])
        temp = empty((len(temparry), 3))
        temp[:, 0] = temparry
        temp[:, 1] = concatenate([gdat.listbgal[n][l] for n in gdat.indxsamptotl])
        for k, strgfeat in enumerate(gdat.liststrgfeatsign[l]):
            temp[:, 2] = concatenate([getattr(gdat, 'list' + strgfeat)[n][l] for n in gdat.indxsamptotl])
            bins = getattr(gdat, 'bins' + strgfeat + 'plot')
            
            print 'strgfeat'
            print strgfeat
            print 'temp[:, 0]'
            summgene(temp[:, 0])
            print 'temp[:, 1]'
            summgene(temp[:, 1])
            print 'temp[:, 2]'
            summgene(temp[:, 2])
            print 'bins'
            print bins
            print 'gdat.pntsprob[l, :, :, :, k]'
            summgene(gdat.pntsprob[l][:, :, :, k])
            gdat.pntsprob[l][:, :, :, k] = histogramdd(temp, bins=(gdat.binslgalpntsprob, gdat.binsbgalpntsprob, bins))[0]
            print 'gdat.pntsprob[l, :, :, :, k]'
            summgene(gdat.pntsprob[l][:, :, :, k])
            print 

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

    ## construct a condensed catalog of elements
    if gdat.condcatl:
        if gdat.verbtype > 0:
            print 'Constructing a condensed catalog...'
            timeinit = gdat.functime()
        retr_condcatl(gdat)
    
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
    if gdat.elemtype == 'lght':
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
    for strgchan in gdat.liststrgchan:
        listtemp = getattr(gdat, 'list' + strgchan)
        if isinstance(listtemp, list):

            # temp
            if strgchan in gdat.liststrgfeat or strgchan == 'indxpntsfull':
                continue
            
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
            stdvtemp = std(posttemp, axis=0)

        setattr(gdat, 'post' + strgchan, posttemp)
        setattr(gdat, 'medi' + strgchan, meditemp)
        setattr(gdat, 'errr' + strgchan, errrtemp)
        setattr(gdat, 'stdv' + strgchan, stdvtemp)
    
    # memory usage
    gdat.meanmemoresi = mean(gdat.listmemoresi, 1)
    gdat.derimemoresi = (gdat.meanmemoresi[-1] - gdat.meanmemoresi[0]) / gdat.numbswep
   
    # temp
    if gdat.numbtrap > 0:
        liststrgbins = ['quad', 'full']
        for l in gdat.indxpopl:
            plot_postbindmaps(gdat, l, 'cumu')
            for strgbins in liststrgbins:
                for strgfeatsign in gdat.liststrgfeatsign[l]:
                    plot_postbindmaps(gdat, l, strgbins, strgfeatsign)

    path = gdat.pathoutpthis + 'pcat.h5'
    writoutp(gdat, path)
    os.system('rm -rf %sgdat*' % gdat.pathoutpthis) 
    # temp
    if False:
        if gdat.makeplot:
            plot_post(pathcatl=pathcatl, verbtype=gdat.verbtype, makeanim=gdat.makeanim, prio=prio)
    else:
        if gdat.makeplot:
            plot_post(gdat=gdat, writ=False, prio=prio)

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
        print 'Parent process has run in %d real seconds, %d CPU seconds.' % (gdat.timerealtotl, gdat.timeproctotl)


def retr_deltlpos(gdat, gdatmodi, indxparapert, stdvparapert):

    numbpert = indxparapert.size 
    
    gdatmodi.thissamp = copy(gdatmodi.thissamptemp)
    
    for k in range(numbpert):
        gdatmodi.thissamp[indxparapert[k]] += stdvparapert[k]
    
    # rescale
    gdatmodi.nextsamp = copy(gdatmodi.thissamp)
    gdatmodi.nextsampvarb = zeros_like(gdatmodi.thissamp) 
    for k in gdat.indxfixpdist:
        gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, '', gdatmodi.thissamp[k], k)
    for k in range(numbpert):
        rscl_elem(gdat, gdatmodi, indxparapert[k])
    gdatmodi.thissamp = copy(gdatmodi.nextsamp)

    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxsampcomp, gdatmodi.thissamp, 'this')

    proc_samp(gdat, gdatmodi, 'this', fast=True)
   
    deltlpos = gdatmodi.thislpostotl - gdatmodi.thislpostotltemp
    deltlpri = gdatmodi.thislpritotl - gdatmodi.thislpritotltemp
    deltllik = gdatmodi.thislliktotl - gdatmodi.thislliktotltemp
    
    if False:
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


def worktrac(pathoutpthis, lock, indxprocwork):
	
    try:
        return work(pathoutpthis, lock, indxprocwork)
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def optihess(gdat, gdatmodi, indxprocwork):
    
    if gdat.verbtype > 0:
        print 'Calculating the Fisher information...'
    
    gdatmodi.propbrth = False      
    gdatmodi.propdeth = False      
   
    pathstdvprop = gdat.pathdataopti + '%s.fits' % gdat.rtag
    
    gdatmodi.thislpritotltemp = copy(gdatmodi.thislpritotl)
    gdatmodi.thislliktotltemp = copy(gdatmodi.thislliktotl)
    gdatmodi.thislpostotltemp = copy(gdatmodi.thislpostotl)
    gdatmodi.thissampvarbtemp = copy(gdatmodi.thissampvarb)
    gdatmodi.thissamptemp = copy(gdatmodi.thissamp)
    gdatmodi.thisindxpntsfulltemp = deepcopy(gdatmodi.thisindxpntsfull)
    
    # temp
    deltparastep = 1e-5

    maxmstdv = 0.1
    fudgstdv = 1e0
    diffparaodim = zeros(3)
    diffparaodim[0] = -deltparastep
    diffparaodim[2] = deltparastep
    diffpara = zeros((3, 3, 2))
    diffpara[0, 0, :] = deltparastep * array([-1., -1.])
    diffpara[0, 2, :] = deltparastep * array([-1., 1.])
    diffpara[2, 0, :] = deltparastep * array([1., -1.])
    diffpara[2, 2, :] = deltparastep * array([1., 1.])
    gdatmodi.dictmodi = {}
    for strg in gdat.liststrgcomptotl:
        gdatmodi.dictmodi['stdv' + strg + 'indv'] = []
        gdatmodi.dictmodi['stdv' + strg + 'indvflux'] = []
    gdatmodi.cntrparasave = 0
    lliktemp = empty(gdat.numbstdp)
    cntr = zeros(gdat.maxmnumbcomp)
    gdatmodi.stdvstdpmatr = zeros((gdat.numbstdp, gdat.numbstdp)) 
    gdatmodi.hess = zeros((gdat.numbstdp, gdat.numbstdp)) 
    deltlpos = zeros((3, 3))

    deltlpos[1, 1] = retr_deltlpos(gdat, gdatmodi, array([0]), array([0.]))
    
    for k in gdat.indxpara:
        if k in gdat.indxfixpprop or k in concatenate(gdatmodi.thisindxsampcomp['comp']):
            indxstdpfrst = gdat.indxstdppara[k]
            for n in gdat.indxpara:
                if n in gdat.indxfixpprop or n in concatenate(gdatmodi.thisindxsampcomp['comp']):
                    indxstdpseco = gdat.indxstdppara[n]
                    if k == n:
                       
                        if False:
                            print 'k'
                            print k
                            print 'gdat.namepara[k]'
                            print gdat.namepara[k]
                            print 'n'
                            print n
                            print 'indxstdpfrst'
                            print indxstdpfrst
                            print 'indxstdpseco'
                            print indxstdpseco
                            print 
                        
                        for a in [0, 2]:
                            # evaluate the posterior
                            deltlpos[a, 1] = retr_deltlpos(gdat, gdatmodi, array([k]), array([diffparaodim[a]]))
                            
                            # increase sample counter for plots
                            gdatmodi.cntrswep += 1
                        
                        gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / deltparastep**2 * fabs(deltlpos[0, 1] + deltlpos[2, 1] - 2. * deltlpos[1, 1])
                        
                        if False:
                            print 'deltlpos'
                            print deltlpos
                            if gdatmodi.hess[indxstdpfrst, indxstdpseco] == 0.:
                                raise Exception('')

                        if k in concatenate(gdatmodi.thisindxsampcomp['comp']):
                            stdv = 1. / sqrt(gdatmodi.hess[indxstdpfrst, indxstdpseco])
                            
                            cntr = 0
                            indxpnts = (k - gdat.indxsampcompinit)
                            for strg in gdat.liststrgcomptotl:
                                if k in concatenate(gdatmodi.thisindxsampcomp[strg]):
                                    indxsampflux = k + 2 - cntr
                                    fluxfact = gdatmodi.thissampvarb[indxsampflux] / gdat.minmflux
                                    if strg == 'flux':
                                        gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * fluxfact**2. / sum(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts])
                                    else:
                                        gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * fluxfact**0.5 / sum(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts])
                                    gdatmodi.dictmodi['stdv' + strg + 'indv'].append(stdv)
                                    gdatmodi.dictmodi['stdv' + strg + 'indvflux'].append(gdatmodi.thissampvarb[indxsampflux])

                                cntr += 1
        
                    else:
                        continue
                        for a in [0, 2]:
                            for b in [0, 2]:
                                # evaluate the posterior
                                deltlpos[a, b] = retr_deltlpos(gdat, gdatmodi, array([k, n]), diffpara[a, b, :])
                            
                                # increase sample counter for plots
                                gdatmodi.cntrswep += 1
                        
                        gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / deltparastep**2 * (deltlpos[1, 1] - deltlpos[1, 0] - deltlpos[0, 1] + deltlpos[0, 0])
            
                    if gdat.verbtype > 0 and not isfinite(gdatmodi.hess[indxstdpfrst, indxstdpseco]):
                        print 'Proposal scale estimate went infinite.'
                        print gdat.namepara[k]
                        raise Exception('')
                        print

        if gdat.verbtype > 0:
            gdatmodi.cntrparasave = tdpy.util.show_prog(k, gdat.numbpara, gdatmodi.cntrparasave, indxprocwork=indxprocwork)

    #gdatmodi.stdvstdpmatr[:gdat.numbstdpfixp, :gdat.numbstdpfixp] = linalg.inv(gdatmodi.hess[:gdat.numbstdpfixp, :gdat.numbstdpfixp])
    gdatmodi.stdvstdpmatr[:gdat.numbstdpfixp, :gdat.numbstdpfixp] = 1. / sqrt(gdatmodi.hess[:gdat.numbstdpfixp, :gdat.numbstdpfixp])

    gdatmodi.stdvstdpmatr *= 2.38 / sqrt(gdat.numbpara) * fudgstdv
    
    indx = where(logical_not(isfinite(gdatmodi.stdvstdpmatr)))
    for k in range(indx[0].size):
        if indx[0][k] == indx[1][k]:
            print 'Bad estimation of the proposal scale'
            print 'gdat.namestdp[indx[k]]'
            print gdat.namestdp[indx[0][k]]
            print
            #raise Exception('')

        gdatmodi.stdvstdpmatr[indx] = maxmstdv
            
    for strg in gdat.liststrgcomptotl:
        gdatmodi.dictmodi['stdv' + strg + 'indv'] = array(gdatmodi.dictmodi['stdv' + strg + 'indv'])
        gdatmodi.dictmodi['stdv' + strg + 'indvflux'] = array(gdatmodi.dictmodi['stdv' + strg + 'indvflux'])
    
    proc_samp(gdat, gdatmodi, 'this')
    
    gdatmodi.stdvstdp = gdatmodi.stdvstdpmatr[gdat.indxstdp, gdat.indxstdp]
    
    if gdat.makeplot:
        
        if gdat.numbproc > 1:
            lock.acquire()
    
        xdat = gdat.indxstdp
        ydat = gdatmodi.stdvstdp
        path = gdat.pathopti + 'stdv%d.pdf' % indxprocwork
        tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist', limtydat=[amin(ydat) / 2., 2. * amax(ydat)])
        
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


def worksamp(gdat, lock):

    if gdat.verbtype > 0:
        print 'Writing the global state to the disc before spawning workers...'
    path = gdat.pathoutpthis + 'gdatinit'
    writfile(gdat, path) 
    
    if gdat.numbproc == 1:
        worktrac(gdat.pathoutpthis, lock, 0)
    else:
        if gdat.verbtype > 0:
            print 'Forking the sampler...'

        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(worktrac, gdat.pathoutpthis, lock)
        pool.map(workpart, gdat.indxproc)

        pool.close()
        pool.join()


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
    for k, namefixp in enumerate(gdat.namefixp):
        if gdat.inittype == 'rand':
            if k in gdat.indxfixpnumbpnts:
                gdatmodi.thissamp[gdat.indxfixp[k]] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[k] + 1))
            else:
                gdatmodi.thissamp[gdat.indxfixp[k]] = rand()
        else:
            if gdat.datatype == 'mock':
                if k in gdat.indxfixpnumbpnts and (gdat.truefixp[k] < gdat.minmnumbpnts[k] or gdat.truefixp[k] > gdat.maxmnumbpnts[k]):
                    gdatmodi.thissamp[k] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[k] + 1))
                else:  
                    indxtrue = where(gdat.truenamefixp == namefixp)[0]
                    if indxtrue.size > 0:
                        gdatmodi.thissamp[k] = cdfn_fixp(gdat, 'true', gdat.truefixp[indxtrue], indxtrue)
                    else:
                        gdatmodi.thissamp[k] = rand()
                if (gdatmodi.thissamp[k] > 1 or gdatmodi.thissamp[k] < 0) and not k in gdat.indxfixpnumbpnts:
                    raise Exception('')

            if gdat.datatype == 'inpt':
                if k in gdat.indxfixppsfp:
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, 'true', gdat.exprpsfp[k-gdat.indxfixppsfpinit], k)
                elif k in gdat.indxfixp:
                    varb = getattr(gdat, 'true' + gdat.namepara[k])
                    gdatmodi.thissamp[gdat.indxfixp[k]] = cdfn_fixp(gdat, 'true', varb, k)
        
        if gdat.inittype == 'pert' and not k in gdat.indxfixpnumbpnts:
            while True:
                gdatmodi.thissamp[k] += 1e-2 * randn()
                if gdatmodi.thissamp[k] > 0. and gdatmodi.thissamp[k] < 1.:
                    break
        
        gdatmodi.thissampvarb[k] = icdf_fixp(gdat, '', gdatmodi.thissamp[k], k)
    
    ## lists of occupied and empty transdimensional parameters
    gdatmodi.thisindxpntsfull = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            gdatmodi.thisindxpntsfull.append(range(gdatmodi.thissamp[gdat.indxfixpnumbpnts[l]].astype(int)))
    else:
        gdatmodi.thisindxpntsfull = []
    
    gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
    
    if gdat.numbtrap > 0:
        
        for l in gdat.indxpopl:
            gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l]] = rand(gdatmodi.thisindxsampcomp['comp'][l].size)
        
        ## element parameters
        if gdat.inittype == 'refr':
            for l in gdat.indxpopl:
                gdatmodi.thissamp[gdatmodi.thisindxsampcomp['lgal'][l]] = cdfn_self(gdat.truelgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                gdatmodi.thissamp[gdatmodi.thisindxsampcomp['bgal'][l]] = cdfn_self(gdat.truebgal[l], -gdat.maxmgang, 2. * gdat.maxmgang)
                
                try:
                    indxtruepntsgood = where(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :] > gdat.minmflux)[0]
                    flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                    if gdat.fluxdisttype[l] == 'powr':
                        fluxunit = cdfn_powr(flux, gdat.minmflux, gdat.maxmflux, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]])
                    if gdat.fluxdisttype[l] == 'bind':
                        fluxdistnorm = gdatmodi.thissampvarb[gdat.indxfixpfluxdistnorm[l, :]]
                        fluxunit = cdfn_bind(flux, gdat.minmflux, gdat.maxmflux, gdat.binsflux, fluxdistnorm)
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['flux'][l][indxtruepntsgood]] = fluxunit[indxtruepntsgood]
                except:
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['flux'][l]] = rand(gdat.truenumbpnts[l])
                    
                if gdat.elemtype == 'lens':
                    
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['asca'][l]] = cdfn_gaus(gdat.trueasca[l], gdatmodi.thissampvarb[gdat.indxfixpascadistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpascadiststdv[l]])
                    
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['acut'][l]] = cdfn_gaus(gdat.trueacut[l], gdatmodi.thissampvarb[gdat.indxfixpacutdistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpacutdiststdv[l]])
                
                if gdat.elemtype == 'lght' and gdat.numbener > 1:
                    
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['sind'][l]] = cdfn_gaus(gdat.truesind[l], gdatmodi.thissampvarb[gdat.indxfixpsinddistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpsinddiststdv[l]])
                    if gdat.spectype[l] == 'curv':
                        try:
                            if gdat.verbtype > 0:
                                print 'Initializing the spectral curvatures randomly for population %d' % l
                            gdatmodi.thissamp[gdatmodi.thisindxsampcomp['curv'][l]] = cdfn_gaus(gdat.truecurv[l], gdatmodi.thissampvarb[gdat.indxfixpcurvdistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpcurvdiststdv[l]])
                        except:
                            gdatmodi.thissamp[gdatmodi.thisindxsampcomp['curv'][l]] = rand(gdat.truenumbpnts[l])
                            
                    if gdat.spectype[l] == 'expo':
                        try:
                            gdatmodi.thissamp[gdatmodi.thisindxsampcomp['expo'][l]] = cdfn_gaus(gdat.trueexpo[l], gdatmodi.thissampvarb[gdat.indxfixpexpodistmean[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdat.indxfixpexpodiststdv[l]])
                        except:
                            if gdat.verbtype > 0:
                                print 'Initializing the spectral cutoff energies randomly for population %d' % l
                            gdatmodi.thissamp[gdatmodi.thisindxsampcomp['expo'][l]] = rand(gdat.truenumbpnts[l])
            
        if gdat.inittype == 'pert':
            for l in gdat.indxpopl:
                print 'gdatmodi.thissamp[gdatmodi.thisindxsampcomp[comp][l]]'
                print gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l]]
                for k in range(gdatmodi.thisindxsampcomp['comp'][l].size):
                    while True:
                        gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l][k]] += 1e-2 * randn()#randn(gdatmodi.thisindxsampcomp[l].size)
                        if gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l][k]] > 0. and gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l][k]] < 1.:
                            break
                print 'gdatmodi.thissamp[gdatmodi.thisindxsampcomp[comp][l]]'
                print gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l]]
                print
        
    if gdat.verbtype > 1:
        print 'thissamp'
        for k in gdat.indxpara:
            if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
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
    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxsampcomp, gdatmodi.thissamp, 'this')
    
    if gdat.verbtype > 1:
        print 'thissamp, thissampvarb'
        for k in gdat.indxpara:
            if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
                print
            print '%15f %15f' % (gdatmodi.thissamp[k], gdatmodi.thissampvarb[k])
    
    ## sample index
    gdatmodi.cntrswep = 0
    gdatmodi.cntroptiprop = 0
   
    ## initial predicted count maps
    if gdat.elemtype == 'lght':
        # temp
        if gdat.boolintpanglcosi:
            binsangltemp = gdat.binsanglcosi
        else:
            binsangltemp = gdat.binsanglplot
        
    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    # temp
    #gdatmodi.listindxpntsfull = []
    #if gdat.trueinfo:
    #    gdatmodi.listspecassc = []
    
    gdatmodi.thismemoresi = zeros(1)
    gdatmodi.thisdeltlliktotl = zeros(1)
    gdatmodi.thischroproc = zeros(gdat.numbchroproc)
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
    
    # process the initial sample, define the variables to be processed in each sample
    proc_samp(gdat, gdatmodi, 'this')

    workdict = {}
    for strgvarb in gdat.liststrgvarbarry:
        valu = getattr(gdatmodi, 'this' + strgvarb)
        if strgvarb in gdat.liststrgvarbarryswep:
            shap = [gdat.numbswep] + list(valu.shape)
        else:
            shap = [gdat.numbsamp] + list(valu.shape)
        workdict['list' + strgvarb] = zeros(shap)
    for strgvarb in gdat.liststrgvarblistsamp:
        workdict['list' + strgvarb] = []
    
    ## saved state of the sample index used for logging progress status
    gdatmodi.percswepsave = -1
   
    # store the initial sample as the best fit sample
    gdatmodi.maxmllikswep = sum(gdatmodi.thisllik)
    gdatmodi.indxswepmaxmllik = -1 
    gdatmodi.sampvarbmaxmllik = copy(gdatmodi.thissampvarb)
   
    # proposal scale optimization
    if gdat.optillik:
        gdatmodi.listllikopti = []
    if gdat.optiprop:
        gdatmodi.cntrstdpmodi = 0

    # log the initial state
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdatmodi, 'gdatmodi')

    if gdat.verbtype > 0:
        print 'Sampling...'

    if gdat.verbtype > 1:
        print 'gdat.stdvstdp'
        print gdat.stdvstdp[:, None]
        print
    
    if gdat.optiprop:
        gdatmodi.thisstdvstdp = copy(gdat.stdvstdp)

    gdatmodi.optidone = False 

    while gdatmodi.cntrswep < gdat.numbswep:
        
        gdatmodi.thischroproc[:] = 0.
        gdatmodi.thischrototl[:] = 0.
        
        if gdat.optihess:
            if gdat.optillik:
                booltemp = len(gdatmodi.listllikopti) > 10 and amax(diff(array(gdatmodi.listllikopti))[-10:]) < 1e-3
            else:
                booltemp = True
            if booltemp:
                optihess(gdat, gdatmodi, indxprocwork)
                gdatmodi.optidone = True
            
        if gdatmodi.optidone:
            thisfile = h5py.File(gdat.pathoutpthis + 'opti.h5', 'w')
            thisfile.create_dataset('stdvstdp', data=gdat.stdvstdp)
            thisfile.close()
            break

        if gdat.emptsamp:
            continue

        timetotlinit = gdat.functime()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                                                                          and gdat.makeplotfram and gdat.makeplot and not gdat.opti
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
            if False and gdat.elemtype == 'lght':
                
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
                    raise Exception('Number of elements is inconsistent with the element index list.')

                if gdat.numbtrap > 0:
                    if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]] != len(gdatmodi.thisindxpntsfull[l]):
                        raise Exception('Number of elements is inconsistent across data structures.')
                
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['flux'][l]]
                indxtemp = where((flux < gdat.minmflux) | (flux > gdat.maxmflux))[0]
                if indxtemp.size > 0:
                    print 'Spectrum of an element went outside the prior range.'
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
                    raise Exception('Spectrum of an element went outside the prior range.')
        
            timefinl = gdat.functime()
            gdatmodi.thischrototl[2] = timefinl - timeinit
       
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            timeinit = gdat.functime()
        
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this')
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]
            
            # fill the sample lists
            for strgvarb in gdat.liststrgvarbarrysamp:
                valu = getattr(gdatmodi, 'this' + strgvarb)
                workdict['list' + strgvarb][indxsampsave, ...] = valu
            for strgvarb in gdat.liststrgvarblistsamp:
                workdict['list' + strgvarb].append(deepcopy(getattr(gdatmodi, 'this' + strgvarb)))
            
            # temp
            #gdatmodi.listindxpntsfull.append(deepcopy(gdatmodi.thisindxpntsfull))
            #if gdat.trueinfo:
            #    gdatmodi.listspecassc.append(deepcopy(gdatmodi.thisspecassc))
            
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
        if False and gdat.elemtype == 'lght':
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
            
            if gdat.calcllik:
                proc_samp(gdat, gdatmodi, 'next')
                gdatmodi.thisdeltlliktotl = gdatmodi.nextlliktotl - gdatmodi.thislliktotl
            else:
                proc_samp(gdat, gdatmodi, 'next', lprionly=True)
                gdatmodi.thisdeltlliktotl = 0.
        
            if gdat.verbtype > 1:
                print 'gdatmodi.thisdeltlliktotl'
                print gdatmodi.thisdeltlliktotl
                print 'gdatmodi.deltlpritotl'
                print gdatmodi.nextlpritotl - gdatmodi.thislpritotl
                print 'gdatmodi.thislpautotl'
                print gdatmodi.thislpautotl
                print 
            
            # evaluate the acceptance probability
            gdatmodi.accpprob = exp(gdatmodi.thisdeltlliktotl + gdatmodi.nextlpritotl - gdatmodi.thislpritotl + gdatmodi.thislpautotl + gdatmodi.thislfctprop + \
                                                                                                                                gdatmodi.thisjcbnfact + gdatmodi.thiscombfact)
            
        else:
            gdatmodi.accpprob = 0.
    
        timefinl = gdat.functime()
        gdatmodi.thischrototl[5] = timefinl - timeinit
   
        # accept or reject the proposal
        timeinit = gdat.functime()

        if gdat.optillik:
            booltemp = gdatmodi.thisaccpprop and gdatmodi.thisdeltlliktotl > 0.
        else:
            booltemp = gdatmodi.accpprob >= rand()
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
        for strg in gdat.liststrgvarbarryswep:
            workdict['list' + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi, 'this' + strg)
        
        # log the progress
        if gdat.verbtype > 0:
            gdatmodi.nextpercswep = 10 * int(10. * gdatmodi.cntrswep / gdat.numbswep) 
            if gdatmodi.nextpercswep > gdatmodi.percswepsave:
                gdatmodi.percswepsave = gdatmodi.nextpercswep
                
                minm = max(0, gdatmodi.cntrswep - 100)
                maxm = gdatmodi.cntrswep + 1
                if maxm > minm:
                    fact = 100. / float(maxm - minm)
                    print 'Sweep number %d' % gdatmodi.cntrswep
                    print '%3d%% completed.' % gdatmodi.nextpercswep
                    for k in gdat.indxproptype:
                        numb = where(workdict['listindxproptype'][minm:maxm] == k)[0].size
                        if numb > 10:
                            fact =  100. / float(where(workdict['listindxproptype'][minm:maxm] == k)[0].size)
                            accp = fact * where(logical_and(workdict['listaccp'][minm:maxm], workdict['listindxproptype'][minm:maxm] == k))[0].size
                            print '%s acceptance rate: %.3g%%' % (gdat.legdproptype[k], accp)
                        
                    print 'Number of elements:'
                    print gdatmodi.thissampvarb[gdat.indxfixpnumbpnts].astype(int)
                    print 'Chronometers: '
                    print 'chrototl'
                    for k in range(gdatmodi.thischrototl.size):
                        print '%.3g msec' % (gdatmodi.thischrototl[k] * 1e3)
                    print 'chroproc'
                    for k in range(gdatmodi.thischroproc.size):
                        print '%.3g msec' % (gdatmodi.thischroproc[k] * 1e3)
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
        if gdat.opti:
            if gdat.verbtype > 0:
                print 'Proposal optimization step %d' % gdatmodi.cntroptiprop
                print 'gdatmodi.thisstdvstdp'
                print gdatmodi.thisstdvstdp
    
            if gdat.optiprop:
                print 'hey'
                print 'gdatmodi.accpprob'
                print gdatmodi.accpprob
                if gdatmodi.accpprob[gdatmodi.cntrstdpmodi] > 0.75:
                    gdatmodi.thisstdvstdp[gdatmodi.cntrstdpmodi] = gdatmodi.nextstdvstdp[gdatmodi.cntrstdpmodi]
                gdatmodi.cntrstdpmodi += 1
                if gdatmodi.cntrstdpmodi == gdat.numbstdp:
                    gdatmodi.cntrstdpmodi = 0

            if gdatmodi.thisaccp and gdat.optillik:
                gdatmodi.listllikopti.append(gdatmodi.nextlliktotl)

            gdatmodi.cntroptiprop += 1
        else:
            gdatmodi.cntrswep += 1
    
    for strgvarb in gdat.liststrgvarbarry + gdat.liststrgvarblistsamp:
        valu = workdict['list' + strgvarb]
        setattr(gdatmodi, 'list' + strgvarb, valu)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    path = gdat.pathoutpthis + 'gdatmodi%04d' % indxprocwork
    writfile(gdatmodi, path) 
    
