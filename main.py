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
         diagmode=False, \
         emptsamp=False, \

         # sampler
         numbswep=1000000, \
         numbburn=None, \
         factthin=None, \
        
         showmoreaccp=False, \

         condcatl=True, \
         seedstat=None, \
         randseedelem=False, \
         indxevttincl=None, \
         indxenerincl=None, \
         
         listnamefeatsele=None, \
         burntmpr=False, \
         #burntmpr=True, \

         # evaluate the likelihood inside circles around elements
         evalcirc=None, \
        
         truelgalimps=None, \
         truebgalimps=None, \
         truefluximps=None, \
         truedefsimps=None, \
         shrtfram=False, \
        
         relnprio=False, \
        
         mockonly=False, \
        
         variasca=True, \
         variacut=True, \

         numbspatdims=2, \
         
         # model type
         elemtype='lght', \
    
         strgcnfg=None, \

         # initialization type
         inittype=None, \
         
         procrtag=None, \
         loadvaripara=False, \
         
         propcova=True, \
         propwithsing=True, \
         pertmodleval=None, \
         optihess=True, \
         optillik=False, \
         optiprop=False, \
         regulevi=False, \
         
         intreval=False, \
        
         truestdvdefsdistslop=0.5, \

         strgexprflux=None, \
         strgcatl=None, \
         anglassc=None, \
         nameexpr=None, \
         
         lgalprio=None, \
         bgalprio=None, \

         strgexpo=None, \
         
         numbproc=None, \
         liketype='pois', \
         exprtype=None, \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         pixltype=None, \
        
         fittampldisttype=None, \
        
         allwfixdtrue=True, \
         asscmetrtype='dist', \

         # plotting
         numbswepplot=None, \
         makeplot=True, \
         makeplotinit=True, \
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
         psfninfoprio=True, \
         ## spectral

         # prior
         priotype='logt', \
         priofactdoff=0., \
         
         # lensing
         fittrelnpowr=0., \

         # temp
         margfactmodl=1., \
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
         
         propmeanpnts=True, \
         propdist=True, \
         proppsfp=True, \
         propbacp=True, \
         proplenp=None, \
         propcomp=None, \
         probtran=None, \
         probbrde=1., \
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
            
         jitt=False, \

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
    ## time stamp
    gdat.strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    ## name of the configuration function
    if gdat.strgcnfg == None:
        gdat.strgcnfg = inspect.stack()[1][3]
   
    ## number of burned sweeps
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 2
    
    ## number of sweeps between frame plots
    if gdat.numbswepplot == None:
        gdat.numbswepplot = max(gdat.numbswep / 10, 100000)

    ## number of processes
    gdat.strgproc = os.uname()[1]
    if gdat.numbproc == None:
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu' or gdat.strgproc == 'wise':
            gdat.numbproc = 10
        else:
            gdat.numbproc = 1
    
    ## factor by which to thin the sweeps to get samples
    if gdat.factthin == None:
        gdat.factthin = int(ceil(1e-3 * (gdat.numbswep - gdat.numbburn) * gdat.numbproc))
    
    if gdat.verbtype > 0:
        print 'PCAT v0.1 started at %s' % gdat.strgtimestmp
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)

    # feature to be used to sort elements
    if gdat.elemtype == 'lght':
        gdat.namefeatsort = 'flux'
    if gdat.elemtype == 'lens':
        gdat.namefeatsort = 'defs'
    if gdat.elemtype == 'clus':
        gdat.namefeatsort = 'nobj'
    
    # feature correlated with the significance of elements
    if gdat.elemtype == 'lght':
        gdat.namefeatsign = 'flux'
    if gdat.elemtype == 'lens':
        gdat.namefeatsign = 'rele'
    if gdat.elemtype == 'clus':
        gdat.namefeatsign = 'nobj'
    
    # feature used to model the amplitude of elements
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.namecompampl = gdat.namefeatsign
    if gdat.elemtype == 'lens':
        gdat.namecompampl = 'defs'
    gdat.fittindxcompampl = 2

    # defaults
    if gdat.strgexprflux == None:
        gdat.datatype = 'mock'
    else:
        gdat.datatype = 'inpt'
    
    if gdat.exprtype == None:
        if gdat.elemtype == 'lght':
            gdat.exprtype = 'ferm'
        if gdat.elemtype == 'lens':
            gdat.exprtype = 'hubb'
    
    if gdat.pertmodleval == None:
        gdat.pertmodleval = gdat.propwithsing

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
        if gdat.datatype == 'inpt':
            gdat.inittype = 'rand'
        else:
            gdat.inittype = 'pert'

    if gdat.elemtype == 'lens':
        gdat.hubbexpofact = 1.63050e-19
    
    if gdat.strgexpo == None:
        if gdat.elemtype == 'lens':
            gdat.strgexpo = 2.2e3 / gdat.hubbexpofact
        else:
            gdat.strgexpo = 1.

    if gdat.indxevttfull == None:
        if gdat.exprtype == 'ferm':
            gdat.indxevttfull = arange(4)
        else:
            gdat.indxevttfull = arange(1)
    
    if gdat.exprtype == 'ferm':
        gdat.lablenerunit = 'GeV'
    if gdat.exprtype == 'chan':
        gdat.lablenerunit = 'KeV'
    if gdat.exprtype == 'sdyn':
        gdat.lablenerunit = ''

    if gdat.pixltype == None:
        if gdat.exprtype == 'ferm':
            gdat.pixltype = 'heal'
        else:
            gdat.pixltype = 'cart'
    
    if gdat.listnamefeatsele == None:
        if gdat.elemtype == 'lght':
            gdat.listnamefeatsele = ['flux']
        if gdat.elemtype == 'lens':
            gdat.listnamefeatsele = ['defs', 'mcut', 'rele']
    
    if gdat.proplenp == None:
        if gdat.elemtype == 'lens':
            gdat.proplenp = True
        else:
            gdat.proplenp = False

    if gdat.strgexprname == None:
        if gdat.exprtype == 'chan':
            gdat.strgexprname = 'Chandra'
        if gdat.exprtype == 'ferm':
            gdat.strgexprname = 'Fermi-LAT'
        if gdat.exprtype == 'sche':
            gdat.strgexprname = 'XXXXX'
        if gdat.exprtype == 'sdyn':
            gdat.strgexprname = 'TGAS-RAVE'
    
    if gdat.strgcatl == None:
        if gdat.datatype == 'mock':
            gdat.strgcatl = 'true'
        else:
            if gdat.exprtype == 'ferm':
                gdat.strgcatl = '3FGL'
            else:
                gdat.strgcatl = gdat.strgexprname
    
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
            gdat.labllgal = r'L_{z}'
        else:
            gdat.labllgal = r'\theta_1'

    if gdat.lablbgal == None:
        if gdat.exprtype == 'sdyn':
            gdat.lablbgal = r'E_k'
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
                gdat.strgenerfull.append('%.3g %s - %.3g %s' % (gdat.binsenerfull[i], gdat.lablenerunit, gdat.binsenerfull[i+1], gdat.lablenerunit))

    if gdat.evalcirc == None:
        if gdat.elemtype == 'lght':
            gdat.evalcirc = 'locl'
        else:
            gdat.evalcirc = 'full'

    if gdat.elemtype == 'lght':
        gdat.listnamesele = ['pars']
    if gdat.elemtype == 'lens':
        gdat.listnamesele = ['pars', 'nrel']

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
        if gdat.elemtype == 'lght':
            if gdat.exprtype == 'ferm':
                gdat.radispmr = 1. / gdat.anglfact
            if gdat.exprtype == 'chan':
                gdat.radispmr = 0.5 / gdat.anglfact
        if gdat.elemtype == 'lens':
            gdat.radispmr = 0.2 / gdat.anglfact
   
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
    else:
        gdat.numbenerfull = 1
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
   
    gdat.factburntmpr = 0.75
    gdat.numbburntmpr = gdat.factburntmpr * gdat.numbburn

    # model evaluation approximation error tolerance in units of the fraction of the lowest PS flux
    if gdat.specfraceval == None:
        if gdat.exprtype == 'ferm':
            gdat.specfraceval = 0.5
        else:
            gdat.specfraceval = 0.1

    ## generative model
    # set mock sample vector indices
    if gdat.elemtype == 'lens':
        numbpnts = array([20])
    if gdat.elemtype == 'lght':
        numbpnts = array([100])
    setp_namevarbvalu(gdat, 'numbpnts', numbpnts)
    gdat.truenumbpopl = gdat.truenumbpnts.size
    gdat.trueindxpopl = arange(gdat.truenumbpopl)
    
    setp_namevarbvalu(gdat, 'minmnumbpnts', zeros(gdat.truenumbpopl, dtype=int) + 1)
    setp_namevarbvalu(gdat, 'maxmnumbpnts', zeros(gdat.truenumbpopl, dtype=int) + 2000)
     
    for l in gdat.trueindxpopl:
        setattr(gdat, 'trueminmnumbpntspop%d' % l, gdat.trueminmnumbpnts[l])
        setattr(gdat, 'truemaxmnumbpntspop%d' % l, gdat.truemaxmnumbpnts[l])
        
    ### element parameter distributions
    setp_namevarbvalu(gdat, 'spectype', ['powr' for l in gdat.trueindxpopl])
    
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
    
    psfntype = gdat.exprpsfntype
    setp_namevarbvalu(gdat, 'psfntype', psfntype)
    
    #### off-axis profile
    if gdat.exprtype == 'chan':
        oaxitype = True
    else:
        oaxitype = False
    setp_namevarbvalu(gdat, 'oaxitype', oaxitype)

    ### background
    #### template
    if gdat.exprtype == 'ferm':
        back = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    elif gdat.exprtype == 'chan':
        back = [array([1.3e1, 4.]), 1e1]
    else:
        back = [1.]
    if gdat.elemtype == 'clus':
        back = [1.]
    setp_namevarbvalu(gdat, 'back', back)
   
    #### spectra for the background
    if gdat.exprtype == 'chan':
        specback = [None, array([1., 1.3])[gdat.indxenerincl]]
    else:
        specback = [None, None]
    setp_namevarbvalu(gdat, 'specback', specback)
    
    #### label
    lablback = [r'$\mathcal{I}$']
    if gdat.exprtype == 'ferm':
        lablback.append(r'$\mathcal{D}$')
    if gdat.exprtype == 'chan':
        lablback.append(r'$\mathcal{I}_p$')
    setp_namevarbvalu(gdat, 'lablback', lablback)
    
    #### name
    strgback = ['fluxisot']
    if gdat.exprtype == 'ferm':
        strgback.append('fluxfdfm')
    if gdat.exprtype == 'chan':
        strgback.append('fluxpart')
    setp_namevarbvalu(gdat, 'strgback', strgback)
    
    #### legend
    nameback = ['Isotropic']
    if gdat.exprtype == 'ferm':
        nameback.append(r'FDM')
    if gdat.exprtype == 'chan':
        nameback.append(r'Particle')
    setp_namevarbvalu(gdat, 'nameback', nameback)
    
    retr_indxsamp(gdat, strgmodl='true')
    
    ### PSF parameters
    if gdat.psfninfoprio:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                meansigc = gdat.exprpsfp[i * gdat.truenumbpsfptotl + m * gdat.truenumbpsfptotl * gdat.numbener]
                stdvsigc = meansigc * 0.1
                setp_namevarblimt(gdat, 'sigcene%devt%d' % (i, m), [meansigc, stdvsigc], typelimt='meanstdv')
                if gdat.truepsfntype == 'doubking':
                    meangamc = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 1]
                    stdvgamc = meangamc * 0.1
                    setp_namevarblimt(gdat, 'gamcene%devt%d' % (i, m), [meangamc, stdvgamc], typelimt='meanstdv')
                    meansigt = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 2]
                    stdvsigt = meansigt * 0.1
                    setp_namevarblimt(gdat, 'sigtene%devt%d' % (i, m), [meansigt, stdvsigt], typelimt='meanstdv')
                    meangamt = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 3]
                    stdvgamt = meangamt * 0.1
                    setp_namevarblimt(gdat, 'gamtene%devt%d' % (i, m), [meangamt, stdvgamt], typelimt='meanstdv')
                    meanpsff = gdat.exprpsfp[i * 5 + m * 5 * gdat.numbener + 4]
                    stdvpsff = meanpsff * 0.1
                    setp_namevarblimt(gdat, 'psffene%devt%d' % (i, m), [meanpsff, stdvpsff], typelimt='meanstdv')
                elif gdat.trueoaxitype:
                    meanonor = gdat.exprpsfp[i * 3 + m * 3 * gdat.numbener + 1]
                    stdvonor = meanonor * 0.1
                    setp_namevarblimt(gdat, 'onorene%devt%d' % (i, m), [meanonor, stdvonor], typelimt='meanstdv')
                    meanoind = gdat.exprpsfp[i * 3 + m * 3 * gdat.numbener + 2]
                    stdvoind = meanoind * 0.1
                    setp_namevarblimt(gdat, 'oindene%devt%d' % (i, m), [meanoind, stdvoind], typelimt='meanstdv')
    else:
        if gdat.exprtype == 'sdyn':
            minmsigm = 0.1 / gdat.anglfact
            maxmsigm = 0.2 / gdat.anglfact
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
        setp_namevarblimt(gdat, 'sigc', [minmsigm, maxmsigm], ener=True, evtt=True)
        setp_namevarblimt(gdat, 'sigt', [minmsigm, maxmsigm], ener=True, evtt=True)
        setp_namevarblimt(gdat, 'gamc', [minmgamm, maxmgamm], ener=True, evtt=True)
        setp_namevarblimt(gdat, 'gamt', [minmgamm, maxmgamm], ener=True, evtt=True)
        setp_namevarblimt(gdat, 'onor', [minmonor, maxmonor], ener=True, evtt=True)
        setp_namevarblimt(gdat, 'oind', [minmoind, maxmoind], ener=True, evtt=True)
    setp_namevarblimt(gdat, 'psff', [0., 1.], ener=True, evtt=True)
 
    ### normalization
    if gdat.exprtype == 'hubb':
        bacp = [1e-8, 1e-6]
    else:
        bacp = [1e-1, 1e1]
    setp_namevarblimt(gdat, 'bacp', bacp, ener=True, back=True)

    ## hyperparameters
    meanpnts = [0.1, 1000.]
    setp_namevarblimt(gdat, 'meanpnts', meanpnts, popl=True)
    
    ### element parameter boundaries
    #### spatial
    if gdat.exprtype == 'ferm':
        minmgang = 1e-1 / gdat.anglfact
    else:
        minmgang = 1e-2 / gdat.anglfact
    setp_namevarbvalu(gdat, 'minmgang', minmgang)
   
    if gdat.elemtype == 'lght':
        if gdat.exprtype == 'ferm':
            minmflux = 3e-11
        if gdat.exprtype == 'chan':
            minmflux = 1e-9
        if gdat.exprtype == 'sdyn':
            minmflux = 0.1
        setp_namevarbvalu(gdat, 'minmflux', minmflux)
    
    minmdefs = 2e-3 / gdat.anglfact
    setp_namevarbvalu(gdat, 'minmdefs', minmdefs)
    
    minmnobj = 5e0
    setp_namevarbvalu(gdat, 'minmnobj', minmnobj)
    
    setp_namevarblimt(gdat, 'curv', [-1., 1.])

    if gdat.elemtype == 'lght':
        if gdat.exprtype == 'ferm':
            maxmflux = 1e-7
        if gdat.exprtype == 'chan':
            maxmflux = 1e-6
        if gdat.exprtype == 'sdyn':
            maxmflux = 100.
        setp_namevarbvalu(gdat, 'maxmflux', maxmflux)
    
    maxmnobj = 5e2
    setp_namevarbvalu(gdat, 'maxmnobj', maxmnobj)
    
    maxmdefs = 2e-2 / gdat.anglfact
    setp_namevarbvalu(gdat, 'maxmdefs', maxmdefs)
   
    # parameter defaults
    ## distribution
    ### flux
    setp_namevarblimt(gdat, 'gangdistscal', [1. / gdat.anglfact, 10. / gdat.anglfact], popl=True)
    setp_namevarblimt(gdat, 'spatdistcons', [1e-4, 1e-2], popl=True)
    setp_namevarblimt(gdat, 'bgaldistscal', [0.5 / gdat.anglfact, 5. / gdat.anglfact], popl=True)
    if gdat.elemtype == 'lght':
        setp_namevarblimt(gdat, 'fluxdistslop', [1., 4.], popl=True)
    if gdat.elemtype == 'lens':
        if gdat.truestdvdefsdistslop != 'none':
            setp_namevarblimt(gdat, 'defsdistslop', [1.9, gdat.truestdvdefsdistslop], popl=True, typelimt='meanstdv')
        else:
            setp_namevarblimt(gdat, 'defsdistslop', [1., 3.], popl=True)
    if gdat.elemtype == 'clus':
        setp_namevarblimt(gdat, 'nobjdistslop', [1., 3.], popl=True)
    
    ### spectral index
    if gdat.elemtype == 'lght' and gdat.numbener > 1:
        if gdat.exprtype == 'ferm':
            sind = [1., 3.]
            gdat.trueminmsind = 1.
            gdat.truemaxmsind = 3.
        if gdat.exprtype == 'chan':
            gdat.trueminmsind = 0.5
            gdat.truemaxmsind = 2.5
            sind = [0.1, 3.]
        setp_namevarblimt(gdat, 'sinddistmean', sind, popl=True)
        #### standard deviations should not be too small
        setp_namevarblimt(gdat, 'sinddiststdv', [0.3, 2.], popl=True)
        setp_namevarblimt(gdat, 'curvdistmean', [-1., 1.], popl=True)
        setp_namevarblimt(gdat, 'curvdiststdv', [0.1, 1.], popl=True)
        setp_namevarblimt(gdat, 'expodistmean', [1., 8.], popl=True)
        setp_namevarblimt(gdat, 'expodiststdv', [0.01 * gdat.maxmener, gdat.maxmener], popl=True)
    
    ### maximum horizontal/vertical distance of the elements from the image center
    gdat.fittmaxmgang = gdat.maxmgangdata * gdat.margfactmodl
    gdat.truemaxmgang = gdat.fittmaxmgang
    
    if gdat.elemtype == 'lens':
        ### projected scale radius
        gdat.trueminmasca = 0.
        gdat.truemaxmasca = 0.1 / gdat.anglfact
        ### projected cutoff radius
        gdat.trueminmacut = 0.
        gdat.truemaxmacut = 2. / gdat.anglfact
    
    ## lensing
    typelimt='meanstdv'
    # temp
    #setp_namevarblimt(gdat, 'lgalsour', [-gdat.truemaxmgang, gdat.truemaxmgang])
    #setp_namevarblimt(gdat, 'bgalsour', [-gdat.truemaxmgang, gdat.truemaxmgang])
    #setp_namevarblimt(gdat, 'lgalhost', [-gdat.truemaxmgang, gdat.truemaxmgang])
    #setp_namevarblimt(gdat, 'bgalhost', [-gdat.truemaxmgang, gdat.truemaxmgang])
    
    gdat.stdvhostsour = 0.1 / gdat.anglfact
    setp_namevarblimt(gdat, 'lgalsour', [0., gdat.stdvhostsour], typelimt='meanstdv')
    setp_namevarblimt(gdat, 'bgalsour', [0., gdat.stdvhostsour], typelimt='meanstdv')
    setp_namevarblimt(gdat, 'lgalhost', [0., gdat.stdvhostsour], typelimt='meanstdv')
    setp_namevarblimt(gdat, 'bgalhost', [0., gdat.stdvhostsour], typelimt='meanstdv')
    
    for i in gdat.indxener:
        setp_namevarblimt(gdat, 'specsourene%d' % i, array([1e-20, 1e-16]))
    setp_namevarblimt(gdat, 'sizesour', [0.1 / gdat.anglfact, 2. / gdat.anglfact])
    setp_namevarblimt(gdat, 'ellpsour', [0., 0.3])
    for i in gdat.indxener:
        setp_namevarblimt(gdat, 'spechostene%d' % i, array([1e-20, 1e-16]))
    setp_namevarblimt(gdat, 'sizehost', [0.1 / gdat.anglfact, 4. / gdat.anglfact])
    setp_namevarblimt(gdat, 'beinhost', [0.5 / gdat.anglfact, 2. / gdat.anglfact])
    setp_namevarblimt(gdat, 'ellphost', [0., 0.5])
    setp_namevarblimt(gdat, 'sherhost', [0., 0.3])
    setp_namevarblimt(gdat, 'anglsour', [0., pi])
    setp_namevarblimt(gdat, 'anglhost', [0., pi])
    setp_namevarblimt(gdat, 'sanghost', [0., pi])
    
    gdat.trueminmlgal = -gdat.fittmaxmgang
    gdat.truemaxmlgal = gdat.fittmaxmgang
    gdat.trueminmbgal = -gdat.fittmaxmgang
    gdat.truemaxmbgal = gdat.fittmaxmgang
    
    gdat.truefactlgal = gdat.truemaxmlgal - gdat.trueminmlgal
    gdat.truefactbgal = gdat.truemaxmbgal - gdat.trueminmbgal
    gdat.trueminmaang = -pi
    gdat.truemaxmaang = pi
   
    # copy the true model to the inference model if the inference model parameter has not been specified
    temp = deepcopy(gdat.__dict__)
    for strg, valu in temp.iteritems():
        if strg.startswith('true') and not strg[4:].startswith('indx'):
            try:
                valumodl = getattr(gdat, 'fitt' + strg[4:])
                if valumodl == None:
                    raise
                if gdat.verbtype > 1:
                    print 'Received custom input for ' + strg[4:]
            except:
                setattr(gdat, 'fitt' + strg[4:], getattr(gdat, strg))
    
    # check inputs
    if gdat.numbburn > gdat.numbswep:
        raise Exception('Bad number of burn-in sweeps.')
    if gdat.factthin > gdat.numbswep - gdat.numbburn:
        raise Exception('Bad thinning factor.')
    if gdat.pixltype == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')

    if gdat.allwfixdtrue:
        setp_namevarbvalu(gdat, 'spatdistcons', 1e-3, popl=True)
        setp_namevarbvalu(gdat, 'gangdistscal', 4. / gdat.anglfact, popl=True)
        setp_namevarbvalu(gdat, 'bgaldistscal', 2. / gdat.anglfact, popl=True)
        if gdat.elemtype == 'lght':
            setp_namevarbvalu(gdat, 'fluxdistslop', 2.2, popl=True)
        if gdat.elemtype == 'lens':
            setp_namevarbvalu(gdat, 'defsdistslop', 1.9, popl=True)
        if gdat.elemtype == 'clus':
            setp_namevarbvalu(gdat, 'nobjdistslop', 1.9, popl=True)

        setp_namevarbvalu(gdat, 'ascadistmean', 0.05 / gdat.anglfact, popl=True)
        setp_namevarbvalu(gdat, 'ascadiststdv', 0.04 / gdat.anglfact, popl=True)
        setp_namevarbvalu(gdat, 'acutdistmean', 1. / gdat.anglfact, popl=True)
        setp_namevarbvalu(gdat, 'acutdiststdv', 0.04 / gdat.anglfact, popl=True)
            
        if gdat.numbener > 1:
            if gdat.exprtype == 'ferm':
                sinddistmean = 2.15
            if gdat.exprtype == 'chan':
                sinddistmean = 1.
            setp_namevarbvalu(gdat, 'sinddistmean', sinddistmean, popl=True)
            setp_namevarbvalu(gdat, 'sinddiststdv', 0.5, popl=True)
            
            setp_namevarbvalu(gdat, 'curvdistmean', 2., popl=True)
            setp_namevarbvalu(gdat, 'curvdiststdv', 0.2, popl=True)
            
            setp_namevarbvalu(gdat, 'expodistmeanpop1', 2., popl=True)
            setp_namevarbvalu(gdat, 'expodiststdvpop1', 0.2, popl=True)
        
        if gdat.exprtype == 'hubb':
            bacp = 1e-7
        else:
            bacp = 1.
        setp_namevarbvalu(gdat, 'bacp', bacp, ener=True, back=True)
        
        if gdat.elemtype == 'lens':
            setp_namevarbvalu(gdat, 'beinhost', 1.5 / gdat.anglfact)
            setp_namevarbvalu(gdat, 'sizesour', 0.3 / gdat.anglfact)
            setp_namevarbvalu(gdat, 'sizehost', 1. / gdat.anglfact)
            for i in gdat.indxener:
                setp_namevarbvalu(gdat, 'specsourene%d' % i, 1e-19)
                setp_namevarbvalu(gdat, 'spechostene%d' % i, 1e-18)
            setp_namevarbvalu(gdat, 'sanghost', pi / 2.)
   
    if gdat.defa:
        return gdat
    
    if gdat.verbtype > 0:
        if gdat.burntmpr:
            print 'Warning: Tempered burn-in.'

    # initial setup
    setpinit(gdat, True) 
    
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
        # temp
        gdat.factspecener = (gdat.meanener / gdat.enerfluxdist)**(-sqrt(amin(gdat.fittminmsinddistmeanpop0) * amax(gdat.fittmaxmsinddistmeanpop0)))
        gdat.enerexpofact = gdat.enerfluxdist - gdat.meanener
    else:
        gdat.factspecener = array([1.])
 
    for strgmodl in gdat.liststrgmodl:
        for namesele in gdat.listnamesele:
            for namefeat in gdat.listnamefeatsele:
                try:
                    setattr(gdat, 'minm' + namefeat + namesele, getattr(gdat, 'minm' + namefeat))
                    setattr(gdat, 'maxm' + namefeat + namesele, getattr(gdat, 'maxm' + namefeat))
                except:
                    setattr(gdat, strgmodl + 'minm' + namefeat + namesele, getattr(gdat, strgmodl + 'minm' + namefeat))
                    setattr(gdat, strgmodl + 'maxm' + namefeat + namesele, getattr(gdat, strgmodl + 'maxm' + namefeat))
    
    # element features
    ## plot limits for element parameters
    # define minima and maxima
    for namevarb in gdat.fittliststrgfeattotl + list(gdat.fittnamefixp):
        for strglimt in ['minm', 'maxm']:
            for strgmodl in gdat.liststrgmodl:
                try:
                    getattr(gdat, strgmodl + strglimt + namevarb)
                except:
                    setattr(gdat, strgmodl + strglimt + namevarb, getattr(gdat, strglimt + namevarb))
            
            try:
                limt = getattr(gdat, strglimt + namevarb)
            except:
                if strglimt == 'minm':
                    limt = minimum(getattr(gdat, 'fittminm' + namevarb), getattr(gdat, 'fittminm' + namevarb))
                else:
                    limt = maximum(getattr(gdat, 'fittmaxm' + namevarb), getattr(gdat, 'fittmaxm' + namevarb))
                setattr(gdat, strglimt + namevarb, limt)

    # color bars
    gdat.minmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) - 10.
    gdat.maxmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) + 10.
    gdat.scallpdfspatpriointp = 'linr'
    gdat.cmaplpdfspatpriointp = 'PuBu'
    
    gdat.minmllikmaps = -10.
    gdat.maxmllikmaps = 0.
    gdat.scalllikmaps = 'asnh'
    gdat.cmapllikmaps = 'YlGn'
    
    gdat.minmperc = 0.
    gdat.maxmperc = 1e2
    gdat.scalperc = 'asnh'
    gdat.cmapperc = 'afmhot'
    
    gdat.minmpercresi = -1e2
    gdat.maxmpercresi = 1e2
    gdat.scalpercresi = 'asnh'
    gdat.cmappercresi = 'coolwarm'
    
    gdat.scalexpomaps = 'logt'
    gdat.cmapexpomaps = 'OrRd'
    
    gdat.scaldatacnts = 'logt'
    gdat.cmapdatacnts = 'Greys'
    
    gdat.scalresicnts = 'asnh'
    gdat.cmapresicnts = 'RdBu'
    
    gdat.minmconv = 1e-4
    gdat.maxmconv = 10.
    gdat.scalconv = 'logt'
    gdat.cmapconv = 'Purples'
    
    gdat.minmmagn = -1e2
    gdat.maxmmagn = 1e2
    gdat.scalmagn = 'asnh'
    gdat.cmapmagn = 'BrBG'
    
    gdat.minmdeflcomp = 0.01
    gdat.maxmdeflcomp = 10.
    gdat.scaldeflcomp = 'logt'
    gdat.cmapdeflcomp = 'Oranges'
    
    gdat.minmconvelemresi = -0.1
    gdat.maxmconvelemresi = 0.1
    gdat.scalconvelemresi = 'asnh'
    gdat.cmapconvelemresi = 'PiYG'
    
    gdat.minmmagnresi = -10.
    gdat.maxmmagnresi = 10.
    gdat.scalmagnresi = 'asnh'
    gdat.cmapmagnresi = 'PRGn'
    
    liststrgcbar = ['llikmaps', 'perc', 'percresi', 'expomaps', 'lpdfspatpriointp', 'conv', 'magn', 'deflcomp', 'convelemresi', 'magnresi']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)

    gdat.liststrglimt = ['minm', 'maxm']
    for strgmodl in gdat.liststrgmodl:
        for namesele in gdat.listnamesele:
            for namefeat in gdat.listnamefeatsele:
                for strglimt in gdat.liststrglimt:
                    try:
                        getattr(gdat, strglimt + namefeat + namesele)
                    except:
                        setattr(gdat, strglimt + namefeat + namesele, getattr(gdat, strglimt + namefeat))

    # construct the bins for element features
    for strgmodl in gdat.liststrgmodl:
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        liststrgfeatpriototl = getattr(gdat, strgmodl + 'liststrgfeatpriototl')
        for strgfeat in liststrgfeattotl + gdat.liststrgfeatplot:
            scal = getattr(gdat, 'scal' + strgfeat + 'plot')
            maxm = getattr(gdat, 'maxm' + strgfeat)
            minm = getattr(gdat, 'minm' + strgfeat)
            retr_axis(gdat, strgfeat, minm, maxm, gdat.numbbinsplot, scal=scal)
            if strgfeat in liststrgfeatpriototl:
                maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                retr_axis(gdat, strgfeat + 'prio', minm, maxm, gdat.numbbinsplotprio, scal=scal, strginit=strgmodl)
            limt = array([getattr(gdat, 'minm' + strgfeat), getattr(gdat, 'maxm' + strgfeat)])
            setattr(gdat, 'limt' + strgfeat + 'plot', limt)

    # construct bins for the scalar variables
    for namevarbscal in gdat.listnamevarbscal:
        minm = getattr(gdat, 'minm' + namevarbscal)
        maxm = getattr(gdat, 'maxm' + namevarbscal)
        if namevarbscal in gdat.fittnamefixp:
            scal = getattr(gdat, 'fittscal' + namevarbscal)
        else:
            scal = getattr(gdat, 'scal' + namevarbscal)
        
        retr_axis(gdat, namevarbscal, minm, maxm, 50, scal=scal)
    
    # create the PCAT folders
    gdat.pathoutp = gdat.pathdata + 'outp/'
    gdat.pathoutpthis = gdat.pathoutp + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    os.system('mkdir -p %s' % gdat.pathoutpthis)
    pathcatllite = gdat.pathoutpthis + 'catllite.fits'  
    pathcatl = gdat.pathoutpthis + 'catl.fits'  
    
    sys.stdout = logg(gdat)

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
        gdat.exprbgal[0, :], gdat.exprlgal[0, :] = rttr(pi / 2. - gdat.exprbgal[0, :], gdat.exprlgal[0, :])
        gdat.exprbgal[0, :] = pi / 2. - gdat.exprbgal[0, :]

    # external reference catalog
    if gdat.exprinfo:
        
        # find entries inside the ROI
        gdat.indxexprpntsrofi = where((fabs(gdat.exprlgal) < gdat.maxmgangdata) & (fabs(gdat.exprbgal) < gdat.maxmgangdata))[0]

        for strgfeat in gdat.fittliststrgfeattotl:
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
        
        gdat.exprflux = gdat.exprspec[:, gdat.indxenerfluxdist[0], :]
        
        gdat.exprlgal = gdat.exprlgal[:, indxpnts]
        gdat.exprbgal = gdat.exprbgal[:, indxpnts]
        for k in range(3):
            gdat.exprspec[k, :, :] = gdat.exprspec[k, :, indxpnts].T
        gdat.exprflux = gdat.exprflux[:, indxpnts]
        if gdat.numbener > 1:
            gdat.exprsind = gdat.exprsind[:, indxpnts]

        gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
        gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)
        
        if not isfinite(gdat.exprspec).all():
            raise Exception('exprspec is not finite.')
        
        # temp
        #if gdat.exprnumbpnts > 0:
        #    gdat.exprfluxbrgt, gdat.exprfluxbrgtassc = retr_fluxbrgt(gdat, gdat.exprlgal, gdat.exprbgal, gdat.exprflux[0, :])

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
            if gdat.truenumbtrap > 0:
                for l in gdat.trueindxpopl:
                    gdat.truenumbpnts[l] = random_integers(0, gdat.maxmnumbpnts[l])
    
    else:
        if gdat.exprinfo:
            gdat.truelgal = [gdat.exprlgal]
            gdat.truebgal = [gdat.exprbgal]
            gdat.truenumbpnts = array([gdat.exprlgal.shape[1]])
            gdat.trueflux = [gdat.exprflux]
            gdat.truespec = [gdat.exprspec]
            gdat.truesind = [gdat.exprsind]
    
            gdat.trueminmflux = amin(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
            gdat.truemaxmflux = amax(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
            for l in gdat.trueindxpopl: 
                gdat.trueminmflux = min(gdat.trueminmflux, amin(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
                gdat.truemaxmflux = max(gdat.truemaxmflux, amax(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
            
    gdat.truenumbpopl = gdat.truenumbpnts.size
    gdat.trueindxpopl = arange(gdat.truenumbpopl, dtype=int)
    
    if gdat.datatype == 'mock':
        gdat.limtpntshist = [0.5, 10**ceil(log10(max(gdat.truemaxmnumbpntstotl, gdat.fittmaxmnumbpntstotl)))]
    
    ## unit sample vector
    gdat.truesamp = zeros(gdat.truenumbpara)
    gdat.truefixp = zeros(gdat.truenumbfixp) + nan
    if gdat.truenumbtrap > 0:
        gdat.truefixp[gdat.trueindxfixpnumbpnts] = gdat.truenumbpnts
   
    for strgfeat in gdat.fittliststrgfeattotl:
        try:
            hist = histogram(getattr(gdat, 'expr' + strgfeat), getattr(gdat, 'bins' + strgfeat + 'plot'))[0]
            setattr(gdat, 'expr' + strgfeat + 'hist', hist)
        except:
            pass

    if gdat.datatype == 'mock':
        
        if gdat.truenumbtrap > 0:
            gdat.trueindxelemfull = []
            if gdat.truenumbtrap > 0:
                for l in gdat.trueindxpopl:
                    gdat.trueindxelemfull.append(range(gdat.truenumbpnts[l]))
            else:
                gdat.trueindxelemfull = []

            gdat.trueindxsampcomp = retr_indxsampcomp(gdat, gdat.trueindxelemfull, 'true')
            
            gdat.truefixp[gdat.trueindxfixpmeanpnts] = gdat.truefixp[gdat.trueindxfixpnumbpnts]

        gdat.truesampvarb = empty(gdat.truenumbpara)
        gdat.truesamp = rand(gdat.truenumbpara)

        for k in gdat.trueindxfixp:
            
            if gdat.truenumbtrap > 0 and (k in gdat.trueindxfixpnumbpnts or k in gdat.trueindxfixpmeanpnts):
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
                    gdat.truefixp[k] = icdf_fixp(gdat, 'true', gdat.truesamp[k], k)

        # temp
        if gdat.elemtype == 'lens':
            gdat.truefixp[gdat.trueindxfixplgalsour] = 0.04 * gdat.truemaxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalsour] = 0.04 * gdat.truemaxmgang * randn()
            gdat.truefixp[gdat.trueindxfixplgalhost] = 0.04 * gdat.truemaxmgang * randn()
            gdat.truefixp[gdat.trueindxfixpbgalhost] = 0.04 * gdat.truemaxmgang * randn()
        
        gdat.truesampvarb[gdat.trueindxfixp] = gdat.truefixp
        if gdat.truenumbtrap > 0:
            gdat.truesamp[gdat.trueindxfixpnumbpnts] = gdat.truesampvarb[gdat.trueindxfixpnumbpnts]
            gdat.truenumbpntstotl = sum(gdat.truefixp[gdat.trueindxfixpnumbpnts])
            gdat.trueindxpntstotl = arange(gdat.truenumbpntstotl)
        
            if gdat.randseedelem:
                seed()
                for l in gdat.trueindxpopl:
                    gdat.truesamp[gdat.trueindxsampcomp['comp'][l]] = rand(gdat.trueindxsampcomp['comp'][l].size)
   
            # sample element components from the true metamodel
            retr_sampvarbcomp(gdat, 'true', gdat.trueindxsampcomp, gdat.trueindxpopl, gdat.trueliststrgcomp, gdat.truelistscalcomp, gdat.truesamp, gdat.truesampvarb)
    
    if gdat.truelgalimps != None:
        gdat.truesampvarb[gdat.trueindxsampcomp['lgal'][0]] = gdat.truelgalimps
    if gdat.truebgalimps != None:
        gdat.truesampvarb[gdat.trueindxsampcomp['bgal'][0]] = gdat.truebgalimps
    if gdat.truefluximps != None:
        gdat.truesampvarb[gdat.trueindxsampcomp['flux'][0]] = gdat.truefluximps
    if gdat.truedefsimps != None:
        gdat.truesampvarb[gdat.trueindxsampcomp['defs'][0]] = gdat.truedefsimps

    gdat.apixmodl = (gdat.fittmaxmgang / gdat.numbsidecart)**2
    
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        for l in indxpopl:
            liststrgfeatmodu = getattr(gdat, strgmodl + 'liststrgfeatmodu')
            liststrgpdfnmodu = getattr(gdat, strgmodl + 'liststrgpdfnmodu')
            for strgfeat, strgpdfn in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]):
                if strgpdfn == 'tmpl':
                    if gdat.lgalprio == None or gdat.bgalprio == None:
                        gdat.lgalprio = concatenate((gdat.truelgal))
                        gdat.bgalprio = concatenate((gdat.truebgal))
                    gdat.numbspatprio = gdat.lgalprio.size
    
                    # spatial template for the catalog prior
                    # temp -- this should move outside the if
                    gdat.pdfnspatpriotemp = zeros((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                    for k in range(gdat.numbspatprio):
                        gdat.pdfnspatpriotemp[:] += 1. / sqrt(2. * pi) / gdat.stdvspatprio * exp(-0.5 * (gdat.binslgalcartmesh - gdat.lgalprio[k])**2 / gdat.stdvspatprio**2) * \
                                                                                                    exp(-0.5 * (gdat.binsbgalcartmesh - gdat.bgalprio[k])**2 / gdat.stdvspatprio**2)
                    gdat.pdfnspatpriotemp /= amax(gdat.pdfnspatpriotemp)
   
    # temp
    for l in gdat.trueindxpopl:
        setattr(gdat, 'truenumbpntspop%d' % l, gdat.truenumbpnts[l])
        setattr(gdat, 'truemeanpntspop%d' % l, gdat.truenumbpnts[l])
    
    if gdat.datatype == 'mock':
        
        if gdat.elemtype == 'lens':
            proc_samp(gdat, None, 'true', raww=True)
        proc_samp(gdat, None, 'true')
            
        if gdat.intreval:
            plot_genemaps(gdat, None, 'true', 'modlcnts', strgcbar='datacnts', thisindxener=0, thisindxevtt=0, intreval=True)

        if gdat.makeplot and gdat.makeplotinit:
            plot_samp(gdat, None, 'true')
        
        if gdat.seedstat != None:
            seed()
   
    if gdat.datatype == 'inpt':
        retr_datatick(gdat)
    
    if gdat.makeplot and gdat.makeplotinit:
        plot_init(gdat)

    # final setup
    setpfinl(gdat, True) 
    
    if gdat.mockonly:
        if gdat.verbtype > 0:
            print 'Mock dataset is generated. Quitting...'
        return gdat

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
        if gdat.elemtype == 'lght':
            print 'minmflux'
            print gdat.minmflux
            print 'maxmflux'
            print gdat.maxmflux
            print 'minmcnts'
            print gdat.minmcnts
            print 'maxmcnts'
            print gdat.maxmcnts
        if gdat.evalcirc != 'full':
            print 'maxmangleval'
            print gdat.anglfact * gdat.maxmangleval, ' [%s]' % gdat.strganglunit
        print
            
    # process lock for simultaneous plotting
    lock = mp.Manager().Lock()
        
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdat, 'gdat')
    
    # list of variables for which the posterior is calculated at each sweep
    gdat.liststrgvarbarryswep = ['memoresi', 'lpri', 'lfctasym', 'lpriprop', 'lpau', 'deltlliktotl', 'lliktotl', 'chro', 'accpprob', \
                                                                    'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype', 'amplpert']
    # temp
    #gdat.liststrgvarbarryswep = []
    
    if gdat.probbrde < 1.:
        gdat.liststrgvarbarryswep += ['auxipara', 'numbpair', 'ljcbfact', 'lcomfact']
    
    # perform a fudicial processing of a sample vector in order to find the list of variables for which the posterior will be calculated
    if gdat.verbtype > 0:
        print 'Processing a fudicial sample...'
    gdatmodifudi = tdpy.util.gdatstrt()
    gdatmodifudi.thischro = zeros(gdat.numbchro)
    gdatmodifudi.thissamp = rand(gdat.fittnumbpara)
    if gdat.fittnumbtrap > 0:
        gdatmodifudi.thissamp[gdat.fittindxfixpnumbpnts] = 1
        gdatmodifudi.thisindxelemfull = [[] for l in gdat.fittindxpopl]
        for l in gdat.fittindxpopl:
            gdatmodifudi.thisindxelemfull[l].append(0)
        gdatmodifudi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodifudi.thisindxelemfull, 'fitt')
    else:
        gdatmodifudi.thisindxsampcomp = None
    gdatmodifudi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodifudi.thissamp, gdatmodifudi.thisindxsampcomp)
    proc_samp(gdat, gdatmodifudi, 'this')
    gdat.liststrgvarbarrysamp = []
    gdat.liststrgvarblistsamp = []
    for strg, valu in gdatmodifudi.__dict__.iteritems():
        if strg.startswith('this') and not strg[4:] in gdat.liststrgvarbarryswep and isinstance(valu, ndarray):
            gdat.liststrgvarbarrysamp.append(strg[4:])
        if strg.startswith('this') and isinstance(valu, list) and strg != 'thisindxsampcomp' and strg != 'thispsfnkern':
            gdat.liststrgvarblistsamp.append(strg[4:])
   
    # temp
    #gdat.liststrgvarbarrysamp += ['memoresi', 'lpri', 'lfctasym', 'lpriprop', 'lpau', 'deltlliktotl', 'lliktotl', 'chro', 'accpprob', 'stdvsamp', \
    #                                                                                                  'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype']
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
   
    setp_indxswepsave(gdat)
   
    # perform an initial run, sampling from the prior
    if gdat.checprio:
        
        if gdat.verbtype > 0:
            print 'Sampling from the prior...'
        
        gdat.namesampdist = 'prio'
        gdat.legdsampdist = 'Prior'
        gdat.calcllik = False
        
        workopti(gdat, lock)

        # save some variables that will be changed for the prior-only run
        #liststrgvarbsave = ['stdvstdp']
        #
        #for strgvarbsave in liststrgvarbsave:
        #    varb = getattr(gdat, strgvarbsave)
        #    if strgvarbsave == 'stdvstdp':
        #        setattr(gdat, strgvarbsave + 'saveprio', copy(varb))
        #    else:
        #        setattr(gdat, strgvarbsave + 'saveprio', varb)
           
        ## change the variables
        #gdat.stdvstdp[:] = 1e-2

        ## perform sampling
        worksamp(gdat, lock)
        
        ## post process the samples
        proc_post(gdat, prio=True)
        
        # save the prior median of variables
        for strgvarb in gdat.liststrgchan:
            
            # temp
            if strgvarb in gdat.fittliststrgfeattotl or strgvarb == 'indxelemfull':
                continue

            setattr(gdat, 'prio' + strgvarb, getattr(gdat, 'medi' + strgvarb))
        
        gdat.calcllik = True
        
        ## retrieve saved variables
        #for strgvarbsave in liststrgvarbsave:
        #    varb = getattr(gdat, strgvarbsave + 'saveprio')
        #    setattr(gdat, strgvarbsave, varb)

    else:
        gdat.calcllik = True
    
    gdat.namesampdist = 'post'
    gdat.legdsampdist = 'Posterior'

    if gdat.verbtype > 0:
        print 'Sampling from the posterior...'
        
    workopti(gdat, lock)

    # run the sampler
    worksamp(gdat, lock)

    proc_post(gdat)
    #if not os.fork():
    #    # post process the samples
    #    proc_post(gdat)

    if gdat.verbtype > 0:
        print 'The ensemble of catalogs is at ' + pathcatl
        if gdat.makeplot:
            print 'The plots are at ' + gdat.pathplot
        print 'PCAT has run successfully. Returning to the OS...'
    
    return gdat
    

def workopti(gdat, lock):

    inittypesave = gdat.inittype 
    # temp
    #gdat.inittype = 'refr'

    # estimate the covariance
    gdat.opti = gdat.optiprop or gdat.optillik or gdat.optihess
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
    
    gdat.inittype = inittypesave
    
    # temp
    gdat.optiprop = False
    gdat.optillik = False
    gdat.optihess = False
    gdat.opti = False


def initarry( \
             dictvarbvari, \
             dictvarb, \
             sameseed=True, \
             randseedelem=False, \
             makeplotarry=False, \
             liststrgvarboutp=None, \
             listlablinpt=None, \
            ):
    
    if sameseed:
        seedstat = get_state()
    else:
        seedstat = None
    
    strgcnfg = inspect.stack()[1][3]

    numbiter = 0
    for strg, valu in dictvarbvari.iteritems():
        numbiter = len(valu)
        break

    if liststrgvarboutp != None:
        numboutp = len(liststrgvarboutp)
        dictoutp = dict()
        for strgvarb in liststrgvarboutp:
            dictoutp[strgvarb] = [[] for k in range(numbiter)]
    
    dictvarb['seedstat'] = seedstat
    dictvarb['randseedelem'] = randseedelem
    dictvarb['strgcnfg'] = strgcnfg
    
    listgdat = []
    for k in range(numbiter):
        for strgvarb, valu in dictvarbvari.iteritems():
            dictvarb[strgvarb] = valu[k]
        dictvarb['strgcnfg'] = inspect.stack()[1][3] + '_%04d' % k
        listgdat.append(init(**dictvarb))

        if liststrgvarboutp != None:
            for strgvarb in liststrgvarboutp:
                dictoutp[strgvarb][k] = getattr(listgdat[k], strgvarb)
        
    if makeplotarry:
        
        strgtimestmp = tdpy.util.retr_strgtimestmp()
    
        path = os.environ["PCAT_DATA_PATH"] + '/imag/%s/' % strgtimestmp
        os.system('mkdir -p %s' % path)
        for strgvarbvari, varbvari in dictvarbvari.iteritems():
            for strgvarboutp, varboutp in dictoutp.iteritems():
                if strgvarbvari.startswith('fitt'):
                    strgvarbvari = strgvarbvari[4:]
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.plot(varbvari, varboutp)
                axis.set_xticklabels(listlablinpt)
                axis.set_ylabel(getattr(gdat, 'labl' + strgvarboutp))
                print 'strgvarbvari'
                print strgvarbvari
                print 'getattr(gdat, scal + strgvarbvari)'
                print getattr(gdat, 'scal' + strgvarbvari)
                print 'strgvarboutp'
                print strgvarboutp
                print 'getattr(gdat, scal + strgvarboutp)'
                print getattr(gdat, 'scal' + strgvarboutp)
                if getattr(gdat, 'scal' + strgvarbvari) == 'logt':
                    axis.set_xscale('log')
                if getattr(gdat, 'scal' + strgvarboutp) == 'logt':
                    axis.set_yscale('log')
                plt.tight_layout()
                plt.savefig('%s/%s%s.pdf' % (path, strgvarbvari, strgvarboutp))
                plt.close(figr)
    
    if liststrgvarboutp != None:
        return listgdat, dictoutp
    else:
        return listgdat


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
    gdat.liststrgchan = gdat.liststrgvarbarry + ['fixp'] + gdat.liststrgvarblistsamp
    
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
   
    gdat.maxmllikproc = empty(gdat.numbproc)
    gdat.indxswepmaxmllikproc = empty(gdat.numbproc, dtype=int)
    gdat.sampvarbmaxmllikproc = empty((gdat.numbproc, gdat.fittnumbpara))
    for k in gdat.indxproc:
        gdat.maxmllikproc[k] = listgdatmodi[k].maxmllikswep
        gdat.indxswepmaxmllikproc[k] = listgdatmodi[k].indxswepmaxmllik
        gdat.sampvarbmaxmllikproc[k] = listgdatmodi[k].sampvarbmaxmllik
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # Gelman-Rubin test
    if gdat.numbproc > 1:
        gdat.gmrbfixp = zeros(gdat.fittnumbfixp)
        gdat.gmrbstat = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        if gdat.verbtype > 0:
            print 'Computing the Gelman-Rubin TS...'
            timeinit = gdat.functime()
        for k in gdat.fittindxfixp:
            gdat.gmrbfixp[k] = tdpy.mcmc.gmrb_test(gdat.listsampvarb[:, :, gdat.fittindxfixp[k]])
            if not isfinite(gdat.gmrbfixp[k]):
                gdat.gmrbfixp[k] = 0.
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
    for strgvarb in gdat.liststrgvarblistsamp:
        listtemp = []
        for j in gdat.indxsamp:      
            for k in gdat.indxproc:
                listtemp.append(getattr(listgdatmodi[k], 'list' + strgvarb)[j])
        setattr(gdat, 'list' + strgvarb, listtemp)
    
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
    gdat.maxmllikproc -= llikoffs
    
    # find the maximum likelihood and posterior over the chains
    gdat.indxprocmaxmllik = argmax(gdat.maxmllikproc)
    gdat.maxmllik = gdat.maxmllikproc[gdat.indxprocmaxmllik]
    gdat.indxswepmaxmllik = gdat.indxprocmaxmllik * gdat.numbsamp + gdat.indxswepmaxmllikproc[gdat.indxprocmaxmllik]
    gdat.sampvarbmaxmllik = gdat.sampvarbmaxmllikproc[gdat.indxprocmaxmllik, :]
    
    # calculate log-evidence using the harmonic mean estimator
    if gdat.verbtype > 0:
        print 'Estimating the Bayesian evidence...'
        timeinit = gdat.functime()
    
    if gdat.regulevi:
        # regularize the harmonic mean estimator
        ## get an ellipse around the median posterior 
        gdat.elpscntr = percentile(listsamp, 50., axis=0)
        thissamp = rand(gdat.fittnumbpara) * 1e-6
        stdvpara = ones(gdat.fittnumbpara) * 1e-6
        limtpara = zeros((2, gdat.fittnumbpara))
        limtpara[1, :] = 1.
        ## find the samples that lie inside the ellipse
        elpsaxis, minmfunc = tdpy.util.minm(thissamp, retr_elpsfrac, stdvpara=stdvpara, limtpara=limtpara, tolrfunc=1e-6, verbtype=gdat.verbtype, optiprop=True)
        lnorregu = -0.5 * gdat.fittnumbpara * log(pi) + sp.special.gammaln(0.5 * gdat.fittnumbpara + 1.) - sum(elpsaxis)
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
    gdat.listfixp = gdat.listsampvarb[:, gdat.fittindxfixp]
    for k, namefixp in enumerate(gdat.fittnamefixp):
        setattr(gdat, 'list' + namefixp, gdat.listfixp[:, k])
    
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
    
    if gdat.fittnumbtrap > 0:
        gdat.pntsprob = [[] for l in gdat.fittindxpopl]
        for l in gdat.fittindxpopl:
            numb = len(gdat.fittliststrgfeatsign[l])
            gdat.pntsprob[l] = zeros((gdat.numbbgalpntsprob, gdat.numblgalpntsprob, gdat.numbbinsplot, numb))
            temparry = concatenate([gdat.listlgal[n][l] for n in gdat.indxsamptotl])
            temp = empty((len(temparry), 3))
            temp[:, 0] = temparry
            temp[:, 1] = concatenate([gdat.listbgal[n][l] for n in gdat.indxsamptotl])
            for k, strgfeat in enumerate(gdat.fittliststrgfeatsign[l]):
                temp[:, 2] = concatenate([getattr(gdat, 'list' + strgfeat)[n][l] for n in gdat.indxsamptotl])
                bins = getattr(gdat, 'bins' + strgfeat)
                gdat.pntsprob[l][:, :, :, k] = histogramdd(temp, bins=(gdat.binslgalpntsprob, gdat.binsbgalpntsprob, bins))[0]

    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # calculate the autocorrelation of the chains
    if gdat.verbtype > 0:
        print 'Computing the autocorrelation of the chains...'
        timeinit = gdat.functime()
   
    gdat.atcr, gdat.timeatcr = tdpy.mcmc.retr_timeatcr(gdat.listmodlcnts)

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
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
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

            # ensure that transdimensional lists are not included
            # temp
            if strgchan in gdat.fittliststrgfeattotl or strgchan == 'indxelemfull':
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
    
    # write an output file to the disc, indicating that the run has been executed successfully
    filecomp = open(gdat.pathoutpthis + 'comp.txt', 'w')
    filecomp.write('PCAT has run successfully.\n')
    filecomp.close()
    
    # start the timer

    if gdat.verbtype > 0:
        for k in gdat.indxproc:
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdat.timereal[k], gdat.timeproc[k])
        print 'Parent process has run in %d real seconds, %d CPU seconds.' % (gdat.timerealtotl, gdat.timeproctotl)


class logg(object):
    
    def __init__(self, gdat):
        self.terminal = sys.stdout
        gdat.pathstdo = gdat.pathoutpthis + 'stdo.txt'
        self.log = open(gdat.pathstdo, 'a')
        pathlink = gdat.pathplot + 'stdo.txt'
        os.system('ln -s %s %s' % (gdat.pathstdo, pathlink))
    
    def write(self, strg):
        self.terminal.write(strg)
        self.log.write(strg)  

    def flush(self):
        pass


def retr_deltlpos(gdat, gdatmodi, indxparapert, stdvparapert):

    numbpert = indxparapert.size 
    
    gdatmodi.thissamp = copy(gdatmodi.thissamptemp)
    
    for k in range(numbpert):
        gdatmodi.thissamp[indxparapert[k]] += stdvparapert[k]
    
    if gdat.fittnumbtrap > 0:
        indx = arange(gdat.fittnumbpopl, gdat.fittnumbpara)
    else:
        indx = arange(gdat.fittnumbpara)
    if where(gdatmodi.thissamp[indx] < 0.)[0].size > 0 or where(gdatmodi.thissamp[indx] > 1.)[0].size > 0:
        print 'Parameter went outside prior bounds when perturbing...'
        deltlpos = 0.
    
    else:
        # rescale
        gdatmodi.nextsamp = copy(gdatmodi.thissamp)
        gdatmodi.nextsampvarb = zeros_like(gdatmodi.thissamp) 
       
        if gdat.fittnumbtrap > 0:
            for k in gdat.fittindxfixpdist:
                gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, 'fitt', gdatmodi.thissamp[k], k)
            for k in range(numbpert):
                rscl_elem(gdat, gdatmodi, indxparapert[k])
        
        gdatmodi.thissamp = copy(gdatmodi.nextsamp)
        gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)

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


def optihess(gdat, gdatmodi):
    
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
    if gdat.fittnumbtrap > 0:
        gdatmodi.thisindxelemfulltemp = deepcopy(gdatmodi.thisindxelemfull)
    
    # temp
    deltparastep = 1e-4

    maxmstdv = 10.
    if gdat.exprtype == 'ferm':
        fudgstdv = 0.5
    else:
        fudgstdv = 0.5
    diffparaodim = zeros(3)
    diffparaodim[0] = -deltparastep
    diffparaodim[2] = deltparastep
    diffpara = zeros((3, 3, 2))
    diffpara[0, 0, :] = deltparastep * array([-1., -1.])
    diffpara[0, 2, :] = deltparastep * array([-1., 1.])
    diffpara[2, 0, :] = deltparastep * array([1., -1.])
    diffpara[2, 2, :] = deltparastep * array([1., 1.])
    gdatmodi.dictmodi = {}
    for strg in gdat.fittliststrgcomptotl:
        gdatmodi.dictmodi['stdv' + strg + 'indv'] = []
        gdatmodi.dictmodi['stdv' + strg + 'indv' + gdat.namecompampl] = []
    gdatmodi.cntrparasave = 0
    lliktemp = empty(gdat.numbstdp)
    cntr = zeros(gdat.fittmaxmnumbcomp)
    gdatmodi.stdvstdpmatr = zeros((gdat.numbstdp, gdat.numbstdp)) 
    gdatmodi.hess = zeros((gdat.numbstdp, gdat.numbstdp)) 
    deltlpos = zeros((3, 3))

    deltlpos[1, 1] = retr_deltlpos(gdat, gdatmodi, array([0]), array([0.]))
    
    if gdat.propcomp:
        indxsamptranprop = concatenate(gdatmodi.thisindxsampcomp['comp'])
    else:
        indxsamptranprop = []

    for k in gdat.fittindxpara:
        if k in gdat.indxfixpprop or k in indxsamptranprop:
            indxstdpfrst = gdat.indxstdppara[k]
            for n in gdat.fittindxpara:
                if n in gdat.indxfixpprop or n in indxsamptranprop:
                    indxstdpseco = gdat.indxstdppara[n]
                    if k == n:
                        
                        if False:
                            print 'k'
                            print k
                            print 'n'
                            print n
                            print 'gdat.fittnamepara[k]'
                            print gdat.fittnamepara[k]
                            print 'indxstdpfrst'
                            print indxstdpfrst
                            print 'indxstdpseco'
                            print indxstdpseco
                       
                        for a in [0, 2]:
                            # evaluate the posterior
                            deltlpos[a, 1] = retr_deltlpos(gdat, gdatmodi, array([k]), array([diffparaodim[a]]))

                        gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / deltparastep**2 * fabs(deltlpos[0, 1] + deltlpos[2, 1] - 2. * deltlpos[1, 1])
                        
                        if False:
                            print 'deltparastep'
                            print deltparastep
                            print 'deltlpos'
                            print deltlpos
                            print 'gdatmodi.hess[indxstdpfrst, indxstdpseco]'
                            print gdatmodi.hess[indxstdpfrst, indxstdpseco]
                            print
                            #if gdatmodi.hess[indxstdpfrst, indxstdpseco] == 0.:
                            #    raise Exception('')

                        if gdat.fittnumbtrap > 0:
                            if k in concatenate(gdatmodi.thisindxsampcomp['comp']):
                                stdv = 1. / sqrt(gdatmodi.hess[indxstdpfrst, indxstdpseco])
                                
                                cntr = 0
                                indxpnts = (k - gdat.fittindxsampcompinit)
                                for strg in gdat.fittliststrgcomptotl:
                                    if k in concatenate(gdatmodi.thisindxsampcomp[strg]):
                                        indxsampampl = k + gdat.fittindxcompampl - cntr
                                        amplfact = gdatmodi.thissampvarb[indxsampampl] / getattr(gdat, 'minm' + gdat.namecompampl)
                                        if strg == gdat.namecompampl:
                                            gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * amplfact**2. / sum(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts])
                                        else:
                                            gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * amplfact**0.5 / sum(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts])
                                        gdatmodi.dictmodi['stdv' + strg + 'indv'].append(stdv)
                                        gdatmodi.dictmodi['stdv' + strg + 'indv' + gdat.namecompampl].append(gdatmodi.thissampvarb[indxsampampl])

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
                    
                    if not isfinite(gdatmodi.hess[indxstdpfrst, indxstdpseco]) or gdatmodi.hess[indxstdpfrst, indxstdpseco] == 0.:
                        print 'Proposal scale estimate went infinite for %s...' % gdat.fittnamepara[k]

        if gdat.verbtype > 0:
            gdatmodi.cntrparasave = tdpy.util.show_prog(k, gdat.fittnumbpara, gdatmodi.cntrparasave, indxprocwork=gdatmodi.indxprocwork)

    #gdatmodi.stdvstdpmatr[:gdat.numbstdpfixp, :gdat.numbstdpfixp] = linalg.inv(gdatmodi.hess[:gdat.numbstdpfixp, :gdat.numbstdpfixp])
    gdatmodi.stdvstdpmatr[:gdat.numbstdpfixp, :gdat.numbstdpfixp] = 1. / sqrt(gdatmodi.hess[:gdat.numbstdpfixp, :gdat.numbstdpfixp])
    
    if False:
        print 'gdatmodi.hess'
        print gdatmodi.hess
        print 'gdatmodi.stdvstdpmatr'
        print gdatmodi.stdvstdpmatr

    gdatmodi.stdvstdpmatr *= fudgstdv * 2.38 
    if not gdat.propwithsing:
        gdatmodi.stdvstdpmatr /= sqrt(gdat.fittnumbpara)
    
    indx = where(logical_not(isfinite(gdatmodi.stdvstdpmatr)))
    for k in range(indx[0].size):
        if indx[0][k] == indx[1][k]:
            print 'Bad estimation of the proposal scale'
            print 'gdat.namestdp[indx[k]]'
            print gdat.namestdp[indx[0][k]]
            print
            #raise Exception('')

        gdatmodi.stdvstdpmatr[indx] = maxmstdv
            
    for strg in gdat.fittliststrgcomptotl:
        gdatmodi.dictmodi['stdv' + strg + 'indv'] = array(gdatmodi.dictmodi['stdv' + strg + 'indv'])
        gdatmodi.dictmodi['stdv' + strg + 'indv' + gdat.namecompampl] = array(gdatmodi.dictmodi['stdv' + strg + 'indv' + gdat.namecompampl])
    
    proc_samp(gdat, gdatmodi, 'this')
    
    gdatmodi.stdvstdp = gdatmodi.stdvstdpmatr[gdat.indxstdp, gdat.indxstdp]
    
    if gdat.makeplot:
        
        xdat = gdat.indxstdp
        ydat = gdatmodi.stdvstdp
        
        pathopti = getattr(gdat, 'path' + gdat.namesampdist + 'opti')
        path = pathopti + 'stdv%d.pdf' % gdatmodi.indxprocwork
        tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist', limtydat=[amin(ydat) / 2., 2. * amax(ydat)])
        
        if gdat.fittnumbtrap > 0:
            for strgcomp in gdat.fittliststrgcomptotl:
                path = pathopti + 'stdv' + strgcomp + '.pdf'
                factplot = getattr(gdat, 'fact' + strgcomp + 'plot')
                meanplot = getattr(gdat, 'mean' + gdat.namecompampl)
                minm = getattr(gdat, 'minm' + gdat.namecompampl)
                factamplplot = getattr(gdat, 'fact' + gdat.namecompampl + 'plot')
                xdat = [gdatmodi.dictmodi['stdv' + strgcomp + 'indv' + gdat.namecompampl] * factamplplot, meanplot * factamplplot]
                
                if strgcomp == gdat.namecompampl:
                    ydat = [gdatmodi.dictmodi['stdv' + strgcomp + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)] / (meanplot / minm)**2.]
                else:
                    ydat = [gdatmodi.dictmodi['stdv' + strgcomp + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)] / (meanplot / minm)**0.5]
                lablxdat = getattr(gdat, 'labl' + gdat.namecompampl + 'totl')
                scalxdat = getattr(gdat, 'scal' + gdat.namecompampl + 'plot')
                limtxdat = array(getattr(gdat, 'limt' + gdat.namecompampl + 'plot')) * factamplplot
                tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                 lablydat=r'$\sigma_{%s}$' % getattr(gdat, 'labl' + strgcomp), plottype=['scat', 'line'])
                #tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                #                                 lablydat=r'$\sigma_{%s}$%s' % (getattr(gdat, 'labl' + strgcomp), getattr(gdat, 'labl' + strgcomp + 'unit')), plottype=['scat', 'line'])
                
                tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                 lablydat=r'$\sigma_{%s}$' % getattr(gdat, 'labl' + strgcomp), plottype=['scat', 'line'])

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


def samp_randnumbpnts(gdat, gdatmodi, l):
    
    gdatmodi.thissamp[l] = poisson(gdatmodi.thissampvarb[gdat.fittindxfixpmeanpnts[l]])
    gdatmodi.thissamp[l] = min(gdat.fittmaxmnumbpnts[l], gdatmodi.thissamp[l])
    gdatmodi.thissamp[l] = max(gdat.fittminmnumbpnts[l], gdatmodi.thissamp[l])


def work(pathoutpthis, lock, indxprocwork):

    path = pathoutpthis + 'gdatinit'
    gdat = readfile(path) 
    
    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for this chain
    seed()
    
    # empty object to hold chain-specific variables that will be modified by the chain
    gdatmodi = tdpy.util.gdatstrt()
    gdatmodi.lock = lock
    gdatmodi.indxprocwork = indxprocwork

    # construct the initial state
    if gdat.verbtype > 0:
        print 'Initializing the sampler state...'
  
    gdatmodi.thisamplpert = zeros(1)

    ## unit sample vector
    gdatmodi.thissamp = rand(gdat.fittnumbpara)
    gdatmodi.thissampvarb = zeros(gdat.fittnumbpara)
    
    if gdat.elemtype == 'lens' and gdat.inittype == 'rand':
        gdatmodi.thissamp[gdat.fittindxfixplgalhost] = 0.5
        gdatmodi.thissamp[gdat.fittindxfixpbgalhost] = 0.5
        gdatmodi.thissamp[gdat.fittindxfixplgalsour] = 0.5
        gdatmodi.thissamp[gdat.fittindxfixpbgalsour] = 0.5
    
    if gdat.datatype == 'inpt' and gdat.inittype != 'rand':
        raise Exception('')

    ## Fixed-dimensional parameters
    if gdat.inittype == 'refr' or gdat.inittype == 'pert':
        for k, namefixp in enumerate(gdat.fittnamefixp):
            gdatmodi.thissamp[k] = cdfn_fixp(gdat, 'true', gdat.truefixp[k], k)
            gdatmodi.thissampvarb[k] = gdat.truesampvarb[k]#icdf_fixp(gdat, 'fitt', gdatmodi.thissamp[k], k)
    
    if gdat.fittnumbtrap > 0:
        if gdat.inittype == 'rand' or gdat.inittype == 'pert':
            for l in gdat.fittindxpopl:
                samp_randnumbpnts(gdat, gdatmodi, l)
    
        ## lists of occupied and empty transdimensional parameters
        gdatmodi.thisindxelemfull = []
        if gdat.fittnumbtrap > 0:
            for l in gdat.fittindxpopl:
                gdatmodi.thisindxelemfull.append(range(gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[l]].astype(int)))
        else:
            gdatmodi.thisindxelemfull = []
        
        gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxelemfull, 'fitt')
        
        ## element parameters
        if gdat.fittnumbtrap > 0:
            if gdat.inittype == 'refr':
                for l in gdat.fittindxpopl:
                    for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                        try:
                            comp = getattr(gdat, 'true' + strgcomp)[l]
                            minm = getattr(gdat, 'fittminm' + strgcomp)
                            maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                            bins = getattr(gdat, 'bins' + strgcomp)
                            if gdat.fittlistscalcomp[l][k] == 'self':
                                fact = getattr(gdat, 'fittfact' + strgcomp)
                                compunit = cdfn_self(comp, minm, fact)
                            if gdat.fittlistscalcomp[l][k] == 'powrslop' or gdat.fittlistscalcomp[l][k] == 'igam':
                                slop = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[l]]
                                if gdat.fittlistscalcomp[l][k] == 'powrslop':
                                    compunit = cdfn_powr(comp, minm, maxm, slop)
                                if gdat.fittlistscalcomp[l][k] == 'igam':
                                    cutf = getattr(gdat, 'cutf' + strgcomp)
                                    compunit = cdfn_igam(comp, slop, cutf)
                            if gdat.fittlistscalcomp[l][k] == 'gaus':
                                distmean = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[l]]
                                diststdv = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[l]]
                                compunit = cdfn_gaus(comp, distmean, diststdv)
                        except:
                            if gdat.verbtype > 0:
                                print 'Initialization from the reference catalog failed for %s. Sampling randomly...' % strgcomp
                            compunit = rand(gdat.truenumbpnts[l])
                        gdatmodi.thissamp[gdatmodi.thisindxsampcomp[strgcomp][l]] = compunit
                    
            else:
                for l in gdat.fittindxpopl:
                    gdatmodi.thissamp[gdatmodi.thisindxsampcomp['comp'][l]] = rand(gdatmodi.thisindxsampcomp['comp'][l].size)
    else:
        gdatmodi.thisindxsampcomp = None

    if gdat.verbtype > 1:
        print 'thissamp'
        for k in gdat.fittindxpara:
            if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
                print
            print '%15f' % gdatmodi.thissamp[k]

    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.thissamp[gdat.fittnumbpopl:] <= 0.)[0] + gdat.fittnumbpopl
    indxsampbadduppr = where(gdatmodi.thissamp[gdat.fittnumbpopl:] >= 1.)[0] + gdat.fittnumbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'indxsampbadd'
        print indxsampbadd
        print gdat.fittnamepara[indxsampbadd]
        print 'gdatmodi.thissamp'
        print gdatmodi.thissamp[:, None]
        raise Exception('Initial unit sample vector went outside the unit interval...')
    
    ## sample vector
    gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)
    
    if gdat.verbtype > 1:
        print 'thissamp, thissampvarb'
        for k in gdat.fittindxpara:
            if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
                print
            print '%15f %15f' % (gdatmodi.thissamp[k], gdatmodi.thissampvarb[k])
    
    ## sample index
    gdatmodi.cntrswep = 0
    gdatmodi.cntroptiprop = 0
   
    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    gdatmodi.thismemoresi = zeros(1)
    gdatmodi.thisdeltlliktotl = zeros(1)
    # temp
    #gdatmodi.thisstdvsamp = zeros(gdat.fittnumbpara)
    gdatmodi.thisaccpprob = zeros(1)
    gdatmodi.thischro = zeros(gdat.numbchro)
    gdatmodi.thisaccp = zeros(1, dtype=bool)
    gdatmodi.thisaccppsfn = zeros(1, dtype=bool)
    gdatmodi.thisaccpprio = zeros(1, dtype=bool)
    gdatmodi.thisaccpprop = zeros(1, dtype=bool)
    gdatmodi.thisindxproptype = zeros(1, dtype=int)
    gdatmodi.thisauxipara = zeros(gdat.fittmaxmnumbcomp)
    gdatmodi.thisnumbpair = zeros(1, dtype=int)
    gdatmodi.thisljcbfact = zeros(1)
    gdatmodi.thislcomfact = zeros(1)
    gdatmodi.thislpau = zeros(gdat.numblpau)
    gdatmodi.thislfctasym = zeros(1)
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
        
        gdatmodi.thischro[:] = 0.
        
        if gdat.optihess:
            if gdat.optillik:
                booltemp = len(gdatmodi.listllikopti) > 10 and amax(diff(array(gdatmodi.listllikopti))[-10:]) < 1e-3
            else:
                booltemp = True
            if booltemp:
                optihess(gdat, gdatmodi)
                gdatmodi.optidone = True
            
        if gdatmodi.optidone:
            thisfile = h5py.File(gdat.pathoutpthis + 'opti.h5', 'w')
            thisfile.create_dataset('stdvstdp', data=gdatmodi.stdvstdp)
            thisfile.close()
            break

        if gdat.emptsamp:
            continue

        initchro(gdat, gdatmodi, 'totl')
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and gdatmodi.indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                                                   and gdat.makeplotfram and gdat.makeplot and not gdat.opti
        # choose a proposal type
        initchro(gdat, gdatmodi, 'type')
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
        
        stopchro(gdat, gdatmodi, 'type')

        if gdat.burntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
            gdatmodi.thisfactmpr = ((gdatmodi.cntrswep + 1.) / gdat.numbburntmpr)**4
            #gdatmodi.thistmprfactdeltllik = gdatmodi.thisfactmpr
            gdatmodi.thistmprfactdeltllik = gdatmodi.thisfactmpr
            gdatmodi.thistmprfactstdv = 1.#gdatmodi.thisfactmpr
            #gdatmodi.thistmprlposelem = -1000. * (1. - gdatmodi.thisfactmpr) * concatenate(gdatmodi.thisindxsampcomp['comp']).size
            gdatmodi.thistmprlposelem = 0.
        else:
            gdatmodi.thistmprfactdeltllik = 1.
            gdatmodi.thistmprfactstdv = 1.
            gdatmodi.thistmprlposelem = 0. 

        # propose the next sample
        initchro(gdat, gdatmodi, 'prop')
        
        retr_prop(gdat, gdatmodi)
        
        stopchro(gdat, gdatmodi, 'prop')
       
        if gdat.verbtype > 1:
            show_samp(gdat, gdatmodi)
    
        # diagnostics
        if gdat.diagmode:
            
            initchro(gdat, gdatmodi, 'diag')
            
            indxsampbadd = where((gdatmodi.thissamp[gdat.fittnumbpopl:] > 1.) | (gdatmodi.thissamp[gdat.fittnumbpopl:] < 0.))[0] + 1
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
       
            for strgstat in ['this', 'next']:
                for strgvarb in ['samp', 'sampvarb']:
                    varb = getattr(gdatmodi, strgstat + strgvarb)
                    if not isfinite(varb).all():
                        print 'gdatmodi' + strgstat + strgvarb
                        for k in gdat.fittindxpara:
                            print varb[k]
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

            for l in gdat.fittindxpopl:
                if gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[l]] != len(gdatmodi.thisindxelemfull[l]):
                    print 'gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts]'
                    print gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts]
                    print 'gdatmodi.thisindxelemfull'
                    print gdatmodi.thisindxelemfull
                    raise Exception('Number of elements is inconsistent with the element index list.')

                if gdat.fittnumbtrap > 0:
                    if gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[l]] != len(gdatmodi.thisindxelemfull[l]):
                        raise Exception('Number of elements is inconsistent across data structures.')
                
                for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                    if gdat.fittlistscalcomp[l][k] == 'gaus' or gdat.fittlistscalcomp[l][k] == 'igam':
                        continue
                    comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l]]
                    minm = getattr(gdat, 'fittminm' + strgcomp)
                    maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                    indxtemp = where((comp < minm) | (comp > maxm))[0]
                    if indxtemp.size > 0:
                        print 'l'
                        print l
                        print 'strgcomp'
                        print strgcomp
                        print 
                        print 'minm'
                        print minm
                        print 'maxm'
                        print maxm
                        print 'indxtemp'
                        print indxtemp
                        print 'comp'
                        print comp
                        print 'comp[indxtemp]'
                        print comp[indxtemp]
                        print
                        raise Exception('A component of an element went outside the prior range.')
        
            stopchro(gdat, gdatmodi, 'diag')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            initchro(gdat, gdatmodi, 'save')
        
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this')
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]
            
            # fill the sample lists
            for strgvarb in gdat.liststrgvarbarrysamp:
                valu = getattr(gdatmodi, 'this' + strgvarb)
                workdict['list' + strgvarb][indxsampsave, ...] = valu
            
            for strgvarb in gdat.liststrgvarblistsamp:
                workdict['list' + strgvarb].append(deepcopy(getattr(gdatmodi, 'this' + strgvarb)))
                 
            stopchro(gdat, gdatmodi, 'save')

        # plot the current sample
        if thismakefram:
            
            initchro(gdat, gdatmodi, 'plot')
            
            if gdat.verbtype > 0:
                print 'Process %d is in queue for making a frame.' % gdatmodi.indxprocwork
            
            proc_samp(gdat, gdatmodi, 'this')
            
            if gdat.numbproc > 1:
                gdatmodi.lock.acquire()
            
            if gdat.verbtype > 0:
                print 'Process %d started making a frame.' % gdatmodi.indxprocwork
            
            plot_samp(gdat, gdatmodi, 'this')
            
            if gdat.verbtype > 0:
                print 'Process %d finished making a frame.' % gdatmodi.indxprocwork
        
            if gdat.numbproc > 1:
                gdatmodi.lock.release()
        
            stopchro(gdat, gdatmodi, 'plot')
    
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
       
        # determine the acceptance probability
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
                print 'gdatmodi.thislfctasym'
                print gdatmodi.thislfctasym
                print 'gdatmodi.thislcomfact'
                print gdatmodi.thislcomfact
                print 'gdatmodi.thisljcbfact'
                print gdatmodi.thisljcbfact
                print 'gdatmodi.thistmprlposelem'
                print gdatmodi.thistmprlposelem
                print
            
            # temp
            #if gdatmodi.propbrth and (gdatmodi.nextlpritotl - gdatmodi.thislpritotl + gdatmodi.thislpautotl) < -5.:
            #    raise Exception('')

            # evaluate the acceptance probability
            gdatmodi.thisaccpprob[0] = exp(gdatmodi.thistmprfactdeltllik * gdatmodi.thisdeltlliktotl + gdatmodi.thistmprlposelem + gdatmodi.nextlpritotl - \
                                                                                gdatmodi.thislpritotl + gdatmodi.thislpautotl + gdatmodi.thislfctasym + \
                                                                                                                                gdatmodi.thisljcbfact + gdatmodi.thislcomfact)
            
        else:
            gdatmodi.thisaccpprob[0] = 0.
    
        # accept or reject the proposal
        if gdat.optillik:
            booltemp = gdatmodi.thisaccpprop and gdatmodi.thisdeltlliktotl > 0.
        else:
            booltemp = gdatmodi.thisaccpprob[0] >= rand()
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
        
        ## variables to be saved for each sweep
        for strg in gdat.liststrgvarbarryswep:
            workdict['list' + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi, 'this' + strg)
        
        # save the execution time for the sweep
        stopchro(gdat, gdatmodi, 'totl')
        
        workdict['listaccpprob'][gdatmodi.cntrswep, 0] = gdatmodi.thisaccpprob[0]
        
        # temp
        #if gdatmodi.propwith:
        #    workdict['liststdvsamp'][gdatmodi.cntrswep, :] = gdatmodi.thisstdvsamp

        # log the progress
        if gdat.verbtype > 0:
            gdatmodi.nextpercswep = 5 * int(20. * gdatmodi.cntrswep / gdat.numbswep) 
            if gdatmodi.nextpercswep > gdatmodi.percswepsave or thismakefram:
                gdatmodi.percswepsave = gdatmodi.nextpercswep
                
                minm = max(0, gdatmodi.cntrswep - 10000)
                maxm = gdatmodi.cntrswep + 1
                if maxm > minm:
                    
                    proc_samp(gdat, gdatmodi, 'this')
                    
                    print
                    print '--------------'
                    print 'Sweep number %d' % gdatmodi.cntrswep
                    print '%3d%% completed.' % gdatmodi.nextpercswep
                    if gdat.burntmpr:
                        print 'factdeltllik'
                        print gdatmodi.thistmprfactdeltllik
                    indxswepintv = arange(minm, maxm)
                    for k in gdat.indxproptype:
                        boolproptype = workdict['listindxproptype'][indxswepintv, 0] == k
                        boolaccp = workdict['listaccp'][indxswepintv, 0] == 1
                        if gdat.showmoreaccp and gdat.indxproptype[k] in gdat.indxproptypecomp:
                            numbtotl = empty(gdat.numbbinsplot, dtype=int)
                            for a in gdat.indxbinsplot: 
                                binsampl = getattr(gdat, 'bins' + gdat.namecompampl)
                                boolbins = (binsampl[a] < workdict['listamplpert'][indxswepintv, 0]) & (workdict['listamplpert'][indxswepintv, 0]< binsampl[a+1])
                                numbtotl[a] = where(boolproptype & boolbins)[0].size
                        else:
                            numbtotl = where(boolproptype)[0].size

                        if gdat.propwithsing and k in gdat.indxstdp:
                            strgstdvstdp = '%.3g' % gdat.stdvstdp[k]
                        else:
                            strgstdvstdp = ''
                        
                        if gdat.showmoreaccp and gdat.indxproptype[k] in gdat.indxproptypecomp:
                            numbaccp = empty(gdat.numbbinsplot, dtype=int)
                            binsampl = getattr(gdat, 'bins' + gdat.namecompampl)
                            for a in gdat.indxbinsplot: 
                                boolbins = (binsampl[a] < workdict['listamplpert'][indxswepintv, 0]) & (workdict['listamplpert'][indxswepintv, 0]< binsampl[a+1])
                                numbaccp[a] = where(boolaccp & boolproptype & boolbins)[0].size
                            indx = where(numbtotl > 0)[0]
                            percaccp = zeros(gdat.numbbinsplot)
                            if indx.size > 0:
                                percaccp[indx] = 100. / numbtotl[indx].astype(float) * numbaccp[indx]
                        else:
                            numbaccp = where(boolaccp & boolproptype)[0].size
                            percaccp = 0.
                            if numbtotl > 0:
                                percaccp = 100. / float(numbtotl) * numbaccp
                        
                        if gdat.showmoreaccp and gdat.indxproptype[k] in gdat.indxproptypecomp:
                            for a in gdat.indxbinsplot:
                                print '%30s %40s %10s' % ('%s-%02d' % (gdat.legdproptype[k], a), 'acceptance rate: %3d%% (%5d out of %5d)' % \
                                                                                                            (percaccp[a], numbaccp[a], numbtotl[a]), strgstdvstdp)
                        else:
                            print '%30s %40s %10s' % (gdat.legdproptype[k], 'acceptance rate: %3d%% (%5d out of %5d)' % (percaccp, numbaccp, numbtotl), strgstdvstdp)
                        
                    if gdat.burntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
                        print 'Tempered burn-in'
                        print 'gdatmodi.thisfactmpr'
                        print gdatmodi.thisfactmpr
                    print 'gdatmodi.thislliktotl'
                    print gdatmodi.thislliktotl
                    print 'gdatmodi.thislpritotl'
                    print gdatmodi.thislpritotl
                    if gdat.fittnumbtrap > 0:
                        print 'Number of elements:'
                        print gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts].astype(int)
                        print 'Mean number of elements:'
                        print gdatmodi.thissampvarb[gdat.fittindxfixpmeanpnts]
                        print 'Log-slope of the amplitude parameter distribution, population 0:'
                        indxfixp = getattr(gdat, 'fittindxfixp' + gdat.namecompampl + 'distsloppop0')
                        print gdatmodi.thissampvarb[indxfixp]
                        print 'Log-prior penalization term: '
                        print gdatmodi.thislpripena
                        if gdat.trueinfo:
                            print 'Completeness'
                            print gdatmodi.thiscmplpop0
                            print 'Completeness binned in significance parameter'
                            print getattr(gdatmodi, 'thiscmpl' + gdat.namefeatsign)
                            print 'False discovery rate'
                            print gdatmodi.thisfdispop0
                            print 'False discovery rate binned in significance parameter'
                            print getattr(gdatmodi, 'thisfdis' + gdat.namefeatsign)

                    print 'Mean residual'
                    print mean(gdatmodi.thisresicnts)

                    print 'Chronometers: '
                    print 'chro'
                    for k in range(gdat.numbchro):
                        for name, valu in gdat.indxchro.iteritems():
                            if valu == k:
                                print '%s: %.3g msec' % (name, gdatmodi.thischro[k] * 1e3)
                   
                    # temp
                    #if not gdat.propwithsing:
                    #    corrstdv = mean((workdict['liststdvsamp'][indxswepintv, :] - mean(workdict['liststdvsamp'][indxswepintv, :], axis=0)) * \
                    #                                                  (workdict['listaccpprob'][indxswepintv, 0] - mean(workdict['listaccpprob'][indxswepintv, 0]))[:, None], axis=0)
                    #    corrstdv /= std(workdict['liststdvsamp'][indxswepintv, :], axis=0) * std(workdict['listaccpprob'][indxswepintv, 0])
                    #    print 'Acceptance correlations: '
                    #    for k in gdat.indxstdp:
                    #        print '%20s: %.3g ' % (gdat.fittnamepara[k], corrstdv[k])
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
            #if gdat.verbtype > 0:
            #    print 'Proposal optimization step %d' % gdatmodi.cntroptiprop
            #    print gdatmodi.thisstdvstdp
            #if gdat.optiprop:
            #    if gdatmodi.thisaccpprob[gdatmodi.cntrstdpmodi] > 0.75:
            #        gdatmodi.thisstdvstdp[gdatmodi.cntrstdpmodi] = gdatmodi.nextstdvstdp[gdatmodi.cntrstdpmodi]
            #    gdatmodi.cntrstdpmodi += 1
            #    if gdatmodi.cntrstdpmodi == gdat.numbstdp:
            #        gdatmodi.cntrstdpmodi = 0

            if gdatmodi.thisaccp and gdat.optillik:
                gdatmodi.listllikopti.append(gdatmodi.nextlliktotl)

            if gdat.optiprop and gdatmodi.cntrswep % gdat.numbswepoptiprop and gdatmodi.cntrswep != 0 and gdatmodi.cntrswep < gdatmodi.numbburn:
                minm = gdatmodi.cntrswep + 1 - gdat.numbswepoptiprop
                maxm = gdatmodi.cntrswep + 1
                indxswepintv = arange(minm, maxm)
                gdatmodi.accpcumu = zeros(gdat.numbstdp)
                gdatmodi.optifact = zeros(gdat.numbstdp) - 1.
                for k in gdat.indxstdp:
                    numb = where(workdict['listindxproptype'][indxswepintv] == k)[0].size
                    if numb > 0:
                        fact =  100. / float(numb)
                        gdatmodi.accpcumu[k] = fact * where(logical_and(workdict['listaccp'][indxswepintv], workdict['listindxproptype'][indxswepintv] == k))[0].size
                        gdatmodi.optifact[k] = 2.**(gdatmodi.accpcumu[k] - 0.25)
                        gdatmodi.thisstdvstdp[k] *= gdatmodi.optifact[k]
                if gdat.verbtype > 0:
                    print 'Proposal scale adaptation step number %d' % (gdatmodi.cntrswep / gdat.numbswepoptiprop)
                    print 'Sweep number %d' % gdatmodi.cntrswep
                    print '%20s %20s %20s' % ('name', 'factor', 'new scale')
                    print '%20s %20.3g %20.3g' % (gdat.namestdp[k], gdatmodi.optifact[k], gdatmodi.thisstdvstdp[k])

            gdatmodi.cntroptiprop += 1
        else:
            gdatmodi.cntrswep += 1
        
    
    for strgvarb in gdat.liststrgvarbarry + gdat.liststrgvarblistsamp:
        valu = workdict['list' + strgvarb]
        setattr(gdatmodi, 'list' + strgvarb, valu)
    
    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    path = gdat.pathoutpthis + 'gdatmodi%04d' % gdatmodi.indxprocwork
    writfile(gdatmodi, path) 
    
