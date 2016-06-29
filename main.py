# common imports
from __init__ import *

# internal functions
from util import *
from visu import *

def work(gdat, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for the process
    seed()
    
    gdatmodi = gdatstrt()
    gdatmodi.modilgal = empty(gdat.numbmodipnts)
    gdatmodi.modibgal = empty(gdat.numbmodipnts)
    gdatmodi.modisind = empty(gdat.numbmodipnts)
    gdatmodi.modispec = empty((gdat.numbener, gdat.numbmodipnts))

    # construct the run tag
    gdat.rtag = retr_rtag(gdat, indxprocwork)
    
    # boolean flag to indicate that the initial state will not be random
    gdat.detrinit = (not gdat.randinit and gdat.datatype == 'mock')
    
    gdatmodi.deltllik = 0.
    
    # initialize the sample vector
    gdatmodi.drmcsamp = zeros((gdat.numbpara, 2))
    
    ## number of PS
    if gdat.initnumbpnts != None or not gdat.randinit:
        if gdat.initnumbpnts != None:
            thisnumbpnts = gdat.initnumbpnts
        else:
            thisnumbpnts = gdat.truenumbpnts
    else:
        thisnumbpnts = empty(gdat.numbpopl, dtype=int)
        for l in gdat.indxpopl:
            thisnumbpnts[l] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[l] + 1))
    gdatmodi.drmcsamp[gdat.indxsampnumbpnts, 0] = thisnumbpnts
       
    gdatmodi.thisindxpntsfull = []
    gdatmodi.thisindxpntsempt = []
    for l in gdat.indxpopl:
        gdatmodi.thisindxpntsfull.append(range(thisnumbpnts[l]))
        gdatmodi.thisindxpntsempt.append(range(thisnumbpnts[l], gdat.maxmnumbpnts[l]))
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampspec, \
        gdatmodi.thisindxsampsind, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    ## mean number of PSs
    if gdat.initfdfnnorm != None or gdat.detrinit:
        if gdat.initfdfnnorm != None:
            fdfnnorm = gdat.initfdfnnorm
        else:
            fdfnnorm = gdat.truefdfnnorm
        for l in gdat.indxpopl:
            gdatmodi.drmcsamp[gdat.indxsampfdfnnorm[l], 0] = cdfn_logt(fdfnnorm[l], gdat.minmfdfnnorm[l], gdat.factfdfnnorm[l])
    else:
        gdatmodi.drmcsamp[gdat.indxsampfdfnnorm, 0] = rand(gdat.numbpopl)
   
    ## FDF shape
    ### single power law
    if gdat.fdfntype == 'powr':
        if gdat.initfdfnslop != None or gdat.detrinit:
            if gdat.initfdfnslop != None:
                fdfnslop = gdat.initfdfnslop
            else:
                fdfnslop = gdat.truefdfnslop
            for l in gdat.indxpopl:
                gdatmodi.drmcsamp[gdat.indxsampfdfnslop[l], 0] = cdfn_atan(fdfnslop, gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
        else:
            gdatmodi.drmcsamp[gdat.indxsampfdfnslop, 0] = rand(gdat.numbpopl)
    ### broken power law
    if gdat.fdfntype == 'brok':
        if gdat.initfdfnbrek != None or gdat.detrinit:
            if gdat.initfdfnbrek != None:
               fdfnbrek = gdat.initfdfnbrek 
            else:
               fdfnbrek = gdat.truefdfnbrek 
            for l in gdat.indxpopl:
                gdatmodi.drmcsamp[gdat.indxsampfdfnbrek[l], 0] = cdfn_logt(fdfnbrek[l], gdat.minmfdfnbrek[l], gdat.factfdfnbrek[l])
        else:
            gdatmodi.drmcsamp[gdat.indxsampfdfnbrek, 0] = rand(gdat.numbpopl)
            
        if gdat.initfdfnsloplowr != None or gdat.detrinit:
            if gdat.initfdfnsloplowr != None:
                fdfnsloplowr = gdat.initfdfnsloplowr
            else:
                fdfnsloplowr = gdat.truefdfnsloplowr
            for l in gdat.indxpopl:
                gdatmodi.drmcsamp[gdat.indxsampfdfnsloplowr[l], 0] = cdfn_atan(fdfnsloplowr, gdat.minmfdfnsloplowr[l], gdat.factfdfnsloplowr[l])
        else:
            gdatmodi.drmcsamp[gdat.indxsampfdfnsloplowr, 0] = rand(gdat.numbpopl)

        if gdat.initfdfnslopuppr != None or gdat.detrinit:
            if gdat.initfdfnslopuppr != None:
                fdfnslopuppr = gdat.initfdfnslopuppr
            else:
                fdfnslopuppr = gdat.truefdfnslopuppr
            for l in gdat.indxpopl:
                gdatmodi.drmcsamp[gdat.indxsampfdfnslopuppr[l], 0] = cdfn_atan(fdfnslopuppr, gdat.minmfdfnslopuppr[l], gdat.factfdfnslopuppr[l])
        else:
            gdatmodi.drmcsamp[gdat.indxsampfdfnslopuppr, 0] = rand(gdat.numbpopl)

    ## PSF parameters
    if gdat.randinit or gdat.truepsfipara == None or gdat.psfntype != gdat.truepsfntype:
        if gdat.verbtype > 1:
            'Randomly seeding the PSF parameters from the prior...'
        gdatmodi.drmcsamp[gdat.indxsamppsfipara, 0] = retr_randunitpsfipara(gdat)
    else:
        for k in gdat.indxpsfipara:
            gdatmodi.drmcsamp[gdat.indxsamppsfipara[k], 0] = cdfn_psfipara(gdat, gdat.truepsfipara[k], k)

    ## background normalization
    for c in gdat.indxback:
        if gdat.randinit or gdat.datatype != 'mock':
            gdatmodi.drmcsamp[gdat.indxsampnormback[c, :], 0] = rand(gdat.numbener)
        else:
            gdatmodi.drmcsamp[gdat.indxsampnormback[c, :], 0] = cdfn_logt(gdat.truenormback[c, :], gdat.minmnormback[c], gdat.factnormback[c])
    
    ## PS components
    for l in gdat.indxpopl:
        if gdat.randinit:
            gdatmodi.drmcsamp[gdatmodi.thisindxsampcomp[l], 0] = rand(gdatmodi.thisindxsampcomp[l].size)
        else:
            gdatmodi.drmcsamp[gdatmodi.thisindxsamplgal[l], 0] = copy(cdfn_self(gdat.truelgal[l], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg))
            gdatmodi.drmcsamp[gdatmodi.thisindxsampbgal[l], 0] = copy(cdfn_self(gdat.truebgal[l], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg))
            if gdat.fdfntype == 'powr':
                fdfnslop = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfdfnslop[l], 0], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
                fluxunit = cdfn_flux_powr(gdat, gdat.truespec[l][0, gdat.indxenerfdfn[0], :], fdfnslop)
            if gdat.fdfntype == 'brok':
                flux = gdat.truespec[l][0, gdat.indxenerfdfn[0], :]
                fdfnbrek = icdf_logt(gdatmodi.drmcsamp[gdat.indxsampfdfnbrek[l], 0], gdat.minmfdfnbrek[l], gdat.factfdfnbrek[l])
                fdfnsloplowr = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfdfnsloplowr[l], 0], gdat.minmfdfnsloplowr[l], gdat.factfdfnsloplowr[l])
                fdfnslopuppr = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfdfnslopuppr[l], 0], gdat.minmfdfnslopuppr[l], gdat.factfdfnslopuppr[l])
                fluxunit = cdfn_flux_brok(gdat, flux, fdfnbrek, fdfnsloplowr, fdfnslopuppr)
            gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn, :], 0] = copy(fluxunit)
            gdatmodi.drmcsamp[gdatmodi.thisindxsampsind[l], 0] = cdfn_eerr(gdat.truesind[l], gdat.meansdfn[l], gdat.stdvsdfn[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
   
    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.)[0] + gdat.numbpopl
    indxsampbadduppr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.)[0] + gdat.numbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'Initial unit sample vector went outside [0, 1]. Correcting it...'
        print 'bad index vector'
        print indxsampbadd
        gdatmodi.drmcsamp[indxsampbaddlowr, 0] = 0.
        gdatmodi.drmcsamp[indxsampbadduppr, 0] = 1.

    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.drmcsamp[:, 0])
    gdatmodi.thispntsflux, gdatmodi.thispntscnts, gdatmodi.thismodlflux, gdatmodi.thismodlcnts = retr_maps(gdat, gdatmodi.thisindxpntsfull, gdatmodi.thissampvarb)
    gdat.temppntsflux, gdat.temppntscnts, gdat.tempmodlflux, gdat.tempmodlcnts = retr_maps(gdat, gdatmodi.thisindxpntsfull, gdatmodi.thissampvarb)
    indxsamplgaltemp, indxsampbgaltemp, indxsampspectemp, indxsampsindtemp, indxsampcomptemp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    if gdat.verbtype > 1:
        print 'thisindxpntsfull'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxpntsfull[l]
        
        print 'thisindxpntsempt'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxpntsempt  
        
        print 'thisindxsamplgal'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxsamplgal[l]
        
        print 'thisindxsampbgal'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxsampbgal[l]
        
        print 'thisindxsampspec'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxsampspec[l]
        
        print 'thisindxsampsind'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxsampsind[l]
        
        print 'thisindxsampcomp'
        for l in gdat.indxpopl:
            print gdatmodi.thisindxsampcomp[l]

    gdat.nextpntsflux = zeros_like(gdatmodi.thispntsflux)
    gdat.nextmodlflux = zeros_like(gdatmodi.thispntsflux)
    gdat.nextmodlcnts = zeros_like(gdatmodi.thispntsflux)
    gdat.nextllik = zeros_like(gdatmodi.thispntsflux)

    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
    
    if gdat.verbtype > 1:
        print 'drmcsamp'
        for k in gdat.indxpara:
            print gdatmodi.drmcsamp[k, :]
        print 'thispsfipara'
        print gdatmodi.thissampvarb[gdat.indxsamppsfipara]
        print 'thissampvarb'
        for k in gdat.indxpara:
            print gdatmodi.thissampvarb[k]
        
    if False:
        print 'hey'
        print 'gdatmodi'
        print [attr for attr in dir(gdatmodi) if not attr.startswith('__')]
        print
 
    # run the sampler
    listchan = rjmc(gdat, gdatmodi, indxprocwork)
    
    timereal = time.time() - timereal
    timeproc = time.clock() - timeproc
    
    listchan.append(timereal)
    listchan.append(timeproc)
    
    return listchan

    
def init( \
         verbtype=1, \
         numbswepplot=50000, \
         makeplot=True, \
         diagsamp=False, \
         numbswep=2000000, \
         numbburn=None, \
         factthin=None, \
         regitype='ngal', \
         datatype='inpt', \
         randinit=True, \
         maxmgang=None, \
         minmflux=None, \
         maxmflux=None, \
         minmsind=None, \
         maxmsind=None, \
         meansdfn=None, \
         stdvsdfn=None, \
         minmfdfnnorm=None, \
         maxmfdfnnorm=None, \
         minmfdfnslop=None, \
         maxmfdfnslop=None, \
         minmfdfnbrek=None, \
         maxmfdfnbrek=None, \
         minmfdfnsloplowr=None, \
         maxmfdfnsloplowr=None, \
         minmfdfnslopuppr=None, \
         maxmfdfnslopuppr=None, \
         fdfntype='powr', \
         psfntype=None, \
         boolproppsfn=True, \
         boolpropfdfn=True, \
         indxevttincl=arange(2, 4), \
         indxenerincl=arange(5), \
         maxmnumbpnts=array([1000]), \
         initnumbpnts=None, \
         initfdfnnorm=None, \
         initfdfnslop=None, \
         initfdfnsloplowr=None, \
         initfdfnslopuppr=None, \
         # temp -- if datatype == 'inpt' trueinfo should depend on whether truexxxx are provided
         trueinfo=True, \
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
         stdvfdfnnorm=0.05, \
         stdvfdfnslop=0.1, \
         stdvpsfipara=0.1, \
         stdvback=0.04, \
         stdvlbhl=0.1, \
         stdvflux=0.15, \
         stdvsind=0.15, \
         fracrand=0.05, \
         radispmrlbhl=None, \
         stdvspmrsind=0.2, \
         mocknumbpnts=None, \
         mockfdfnslop=None, \
         mockfdfnsloplowr=None, \
         mockfdfnslopuppr=None, \
         mockfdfnbrek=None, \
         mocknormback=ones((1, 1)), \
         mockfdfntype='powr', \
         mockpsfntype=None, \
         strgexpr=None, \
         strgback=None, \
         lablback=None, \
         nameback=None, \
         strgexpo=None, \
         probprop=None, \
        ):
    
    # start the timer
    timetotlreal = time.time()
    timetotlproc = time.clock()
   
    # defaults
    ## Fermi-LAT
    if exprtype == 'ferm':
        if maxmgang == None:
            maxmgang = 20.
        if minmsind == None:
            minmsind = 1.2
        if maxmsind == None:
            maxmsind = 3.2
        if meansdfn == None:
            meansdfn = array([2.2])
        if stdvsdfn == None:
            stdvsdfn = array([0.3])
        if minmfdfnnorm == None:
            minmfdfnnorm = array([1.])
        if maxmfdfnnorm == None:
            maxmfdfnnorm = array([1e3])
        if minmfdfnslop == None:
            minmfdfnslop = array([0.5])
        if maxmfdfnslop == None:
            maxmfdfnslop = array([3.5])
        if minmfdfnbrek == None:
            minmfdfnbrek = array([minmflux])
        if maxmfdfnbrek == None:
            maxmfdfnbrek = array([maxmflux])
        if minmfdfnsloplowr == None:
            minmfdfnsloplowr = minmfdfnslop
        if maxmfdfnsloplowr == None:
            maxmfdfnsloplowr = maxmfdfnslop
        if minmfdfnslopuppr == None:
            minmfdfnslopuppr = minmfdfnslop
        if maxmfdfnslopuppr == None:
            maxmfdfnslopuppr = maxmfdfnslop
        if margsize == None:
            margsize = 1.
        if psfntype == None:
            psfntype = 'singking'
        if mockpsfntype == None:
            mockpsfntype = 'doubking'
        if radispmrlbhl == None:
            radispmrlbhl = 2.
    ## SDSS
    if exprtype == 'sdss':
        if maxmgang == None:
            maxmgang = 20. / 3600.
        if minmsind == None:
            minmsind = array([1.])
        if maxmsind == None:
            maxmsind = array([3.])
        if meansdfn == None:
            meansdfn = array([2.2])
        if stdvsdfn == None:
            stdvsdfn = array([0.5])
        if minmfdfnnorm == None:
            minmfdfnnorm = array([1.])
        if maxmfdfnnorm == None:
            maxmfdfnnorm = array([1e3])
        if minmfdfnslop == None:
            minmfdfnslop = array([1.5])
        if maxmfdfnslop == None:
            maxmfdfnslop = array([2.5])
        if margsize == None:
            margsize = 1.
        if psfntype == None:
            psfntype = 'doubgaus'
        if mockpsfntype == None:
            mockpsfntype = 'doubgaus'
        if radispmrlbhl == None:
            radispmrlbhl = 2. / 3600.
    
    ## number of populations
    numbpopl = maxmnumbpnts.size
        
    # initialize the global object 
    gdat = gdatstrt()
   
    # load the global object
    ## verbosity level
    ### 0 - no output
    ### 1 - progress
    ### 2 - sampler diagnostics
    gdat.verbtype = verbtype
    
    ## plot settings
    ### MCMC time period over which a frame is produced
    gdat.numbswepplot = numbswepplot
    ### flag to control generation of plots
    gdat.makeplot = makeplot

    ## flag to diagnose the sampler
    gdat.diagsamp = diagsamp
    
    ## MCMC setup
    ### number of sweeps, i.e., number of samples before thinning and including burn-in
    gdat.numbswep = numbswep
    ### number of burn-in samples
    gdat.numbburn = numbburn
    ### number of processes
    gdat.numbproc = numbproc
    ### the factor by which to thin the chain
    gdat.factthin = factthin
    
    ## region type
    ###- ngal - NGPR (North Galactic Polar Region)
    ###- igal - IG (inner Galaxy)
    gdat.regitype = regitype

    ## data type
    ###- mock - mock data
    ###- inpt - input data
    gdat.datatype = datatype
    
    ## likelihood function type
    gdat.liketype = liketype
    
    ## experiment type
    gdat.exprtype = exprtype
    
    ## pixelization type
    gdat.pixltype = pixltype
    
    ## PSF model type
    gdat.psfntype = psfntype
    
    ## flag to turn off PSF parameter updates
    gdat.boolproppsfn = boolproppsfn
    gdat.boolpropfdfn = boolpropfdfn

    ## input data
    ### measured data
    gdat.strgexpr = strgexpr
    ### background
    gdat.strgback = strgback
    gdat.lablback = lablback
    gdat.nameback = nameback
    ### exposure
    gdat.strgexpo = strgexpo
    
    ## flag to use truth information
    gdat.trueinfo = trueinfo
    
    ## Flux distribution function model
    gdat.fdfntype = fdfntype
    
    ## initial state setup
    ### number of point sources
    gdat.initnumbpnts = initnumbpnts
    gdat.initfdfnnorm = initfdfnnorm
    if gdat.fdfntype == 'powr':
        gdat.initfdfnslop = initfdfnslop
    if gdat.fdfntype == 'brok':
        gdat.initfdfnbrek = initfdfnbrek
        gdat.initfdfnsloplowr = initfdfnsloplowr
        gdat.initfdfnslopuppr = initfdfnslopuppr
    ### flag to draw the initial state from the prior
    gdat.randinit = randinit

    ## energy bins to be included
    gdat.indxenerincl = indxenerincl
    
    ## PSF bins to be included
    gdat.indxevttincl = indxevttincl
    
    ## number of Point Source (PS) populations
    gdat.numbpopl = numbpopl
    
    ## maximum angle from the PSs to evaluate the likelihood
    gdat.maxmangleval = maxmangleval
    
    ## parameter limits
    ### FDF normalization
    gdat.minmfdfnnorm = minmfdfnnorm
    gdat.maxmfdfnnorm = maxmfdfnnorm
    if fdfntype == 'powr':
        ### FDF power law index
        gdat.minmfdfnslop = minmfdfnslop
        gdat.maxmfdfnslop = maxmfdfnslop
    if fdfntype == 'brok':
        ### FDF break flux
        gdat.minmfdfnbrek = minmfdfnbrek
        gdat.maxmfdfnbrek = maxmfdfnbrek
        ### FDF lower power law index
        gdat.minmfdfnsloplowr = minmfdfnsloplowr
        gdat.maxmfdfnsloplowr = maxmfdfnsloplowr
        ### FDF upper power law index
        gdat.minmfdfnslopuppr = minmfdfnslopuppr
        gdat.maxmfdfnslopuppr = maxmfdfnslopuppr
    ### flux
    gdat.minmflux = minmflux
    gdat.maxmflux = maxmflux
    ### spectral power-law index
    gdat.minmsind = minmsind
    gdat.maxmsind = maxmsind
    gdat.meansdfn = meansdfn
    gdat.stdvsdfn = stdvsdfn
    ### background normalizations
    gdat.minmnormback = minmnormback
    gdat.maxmnormback = maxmnormback

    ## model indicator limits
    gdat.maxmnumbpnts = maxmnumbpnts
    
    ## Region of interest
    ### image center
    gdat.lgalcntr = lgalcntr
    gdat.bgalcntr = bgalcntr
    ### half of the image size
    gdat.maxmgang = maxmgang
    ### image margin size
    gdat.margsize = margsize

    ## proposals
    ### proposal scales
    gdat.stdvfdfnnorm = stdvfdfnnorm
    gdat.stdvfdfnslop = stdvfdfnslop
    gdat.stdvpsfipara = stdvpsfipara
    gdat.stdvback = stdvback
    gdat.stdvlbhl = stdvlbhl
    gdat.stdvflux = stdvflux
    gdat.stdvsind = stdvsind

    ### fraction of heavy-tailed proposals
    gdat.fracrand = fracrand
    
    ### radius of the circle in which splits and merges are proposed
    gdat.stdvspmrsind = stdvspmrsind
    gdat.radispmrlbhl = radispmrlbhl
    
    ## mock data setup
    if datatype == 'mock':
        ### mock PSF model type
        gdat.mockpsfntype = mockpsfntype
        ### mock FDF model type
        gdat.mockfdfntype = mockfdfntype
        ### mock number of PS
        gdat.mocknumbpnts = mocknumbpnts
        if gdat.fdfntype == 'powr':
            gdat.mockfdfnslop = mockfdfnslop
        if gdat.fdfntype == 'brok':
            gdat.mockfdfnbrek = mockfdfnbrek
            gdat.mockfdfnsloplowr = mockfdfnsloplowr
            gdat.mockfdfnslopuppr = mockfdfnslopuppr
        gdat.mocknormback = mocknormback
        ### flag to position mock point sources at the image center
        gdat.pntscntr = pntscntr
        ### mock image resolution
        gdat.numbsidecart = numbsidecart
        gdat.numbsideheal = numbsideheal

    ## proposal frequencies
    gdat.probprop = probprop

    # temp
    # check inputs
    
    # date and time
    gdat.strgtime = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    if gdat.verbtype > 0:
        print 'PCAT started at ', gdat.strgtime
        print 'Initializing...'
    
    # check the call stack for the name of the configuring function
    gdat.strgcnfg = inspect.stack()[1][3]

    # setup the sampler
    setp(gdat) 

    if gdat.verbtype > 1:
        print 'probprop: '
        print vstack((arange(gdat.numbprop), gdat.strgprop, gdat.probprop)).T
        print 'indxsampnumbpnts: ', gdat.indxsampnumbpnts
        print 'indxsampfdfnnorm: ', gdat.indxsampfdfnnorm
        if gdat.fdfntype == 'powr':
            print 'indxsampfdfnslop: ', gdat.indxsampfdfnslop
        if gdat.fdfntype == 'brok':
            print 'indxsampfdfnbrek: ', gdat.indxsampfdfnbrek
            print 'indxsampfdfnsloplowr: ', gdat.indxsampfdfnsloplowr
            print 'indxsampfdfnslopuppr: ', gdat.indxsampfdfnslopuppr
        print 'indxsamppsfipara: ', gdat.indxsamppsfipara
        print 'indxsampnormback: '
        print gdat.indxsampnormback
        print 'indxsampcompinit: ', gdat.indxsampcompinit
        print 'maxmangleval'
        print gdat.maxmangleval
        print
        print 'Truth information'
        if gdat.trueinfo:
            print 'truelgal'
            for l in gdat.indxpopl:
                print gdat.truelgal[l]
            print 'truebgal'
            for l in gdat.indxpopl:
                print gdat.truebgal[l]
            print 'truespec'
            for l in gdat.indxpopl:
                print gdat.truespec[l]
            print 'truenumbpnts: ', gdat.truenumbpnts
            print 'truepsfipara'
            print gdat.truepsfipara
            print
            if gdat.datatype == 'mock':
                print 'Mock data'
                print 'mocknumbpnts: ', gdat.mocknumbpnts
                if gdat.mockfdfntype == 'powr':
                    print 'mockfdfnslop: ', gdat.mockfdfnslop
                if gdat.mockfdfntype == 'brok':
                    print 'mockfdfnbrek: ', gdat.mockfdfnbrek
                    print 'mockfdfnsloplowr: ', gdat.mockfdfnsloplowr
                    print 'mockfdfnslopuppr: ', gdat.mockfdfnslopuppr
                print 'mockpsfipara: '
                print gdat.mockpsfipara
                print 'mocknormback: '
                print gdat.mocknormback

    # make initial plots
    if gdat.makeplot:
        #plot_3fgl_thrs(gdat)
        #plot_datacntshist()
        #if gdat.exprtype == 'ferm':
        #    plot_fgl3(gdat)
        #plot_intr()
        #plot_king(gdat)
        #plot_look(gdat)
        plot_eval(gdat)
        #if gdat.datatype == 'mock':
        #    plot_pntsdiff()

    if gdat.verbtype > 0:
        print 'Sampling...'
    
    # temp
    if False:
        print 'hey'
        print 'gdat'
        print [attr for attr in dir(gdat) if not attr.startswith('__')]
        print
        #print 'sys.getsizeof()'
        #print sys.getsizeof(gdatmodi.thismodlcnts)

    timereal = zeros(gdat.numbproc)
    timeproc = zeros(gdat.numbproc)
    if gdat.numbproc == 1:
        gridchan = [work(gdat, 0)]
    else:
        if gdat.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        gdat.lock = mp.Manager().Lock()

        # process pool
        pool = mp.Pool(gdat.numbproc)
        
        # spawn the processes
        workpart = functools.partial(work, gdat)
        indxproc = arange(gdat.numbproc)
        gridchan = pool.map(workpart, indxproc)

        pool.close()
        pool.join()

    for k in gdat.indxproc:
        timereal[k] = gridchan[k][17]
        timeproc[k] = gridchan[k][18]

    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        tim0 = time.time()

    # parse the sample bundle
    listsamp = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbpara))
    listsampvarb = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbpara))
    listindxprop = zeros((gdat.numbswep, gdat.numbproc))
    listchrollik = zeros((gdat.numbswep, gdat.numbproc, gdat.numbchrollik))
    listchrototl = zeros((gdat.numbswep, gdat.numbproc, gdat.numbchrototl))
    listllik = zeros((gdat.numbsamp, gdat.numbproc))
    listlpri = zeros((gdat.numbsamp, gdat.numbproc))
    listaccp = zeros((gdat.numbswep, gdat.numbproc))
    listmodlcnts = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbpixlsave))
    listpntsfluxmean = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbener))
    listindxpntsfull = []
    listindxparamodi = zeros((gdat.numbswep, gdat.numbproc), dtype=int)
    gdat.listauxipara = empty((gdat.numbswep, gdat.numbproc, gdat.numbcompcolr))
    gdat.listlaccfrac = empty((gdat.numbswep, gdat.numbproc))
    gdat.listnumbpair = empty((gdat.numbswep, gdat.numbproc))
    gdat.listjcbnfact = empty((gdat.numbswep, gdat.numbproc))
    gdat.listcombfact = empty((gdat.numbswep, gdat.numbproc))

    for k in gdat.indxproc:
        gdat.rtag = retr_rtag(gdat, k)
        listchan = gridchan[k]
        listsamp[:, k, :] = listchan[0]
        listsampvarb[:, k, :] = listchan[1]
        listindxprop[:, k] = listchan[2]
        listchrototl[:, k, :] = listchan[3]
        listllik[:, k] = listchan[4]
        listlpri[:, k] = listchan[5]
        listaccp[:, k] = listchan[6]
        listmodlcnts[:, k, :] = listchan[7]    
        listindxpntsfull.append(listchan[8])
        listindxparamodi[:, k] = listchan[9]
        gdat.listauxipara[:, k, :] = listchan[10]
        gdat.listlaccfrac[:, k] = listchan[11]
        gdat.listnumbpair[:, k] = listchan[12]
        gdat.listjcbnfact[:, k] = listchan[13]
        gdat.listcombfact[:, k] = listchan[14]
        listpntsfluxmean[:, k, :] = listchan[15]
        listchrollik[:, k, :] = listchan[16]

    listindxprop = listindxprop.flatten()
    gdat.listauxipara = gdat.listauxipara.reshape((gdat.numbsweptotl, gdat.numbcompcolr))
    gdat.listlaccfrac = gdat.listlaccfrac.reshape(gdat.numbsweptotl)
    gdat.listnumbpair = gdat.listnumbpair.reshape(gdat.numbsweptotl)
    gdat.listjcbnfact = gdat.listjcbnfact.reshape(gdat.numbsweptotl)
    gdat.listcombfact = gdat.listcombfact.reshape(gdat.numbsweptotl)
    
    listchrototl = listchrototl.reshape((gdat.numbsweptotl, gdat.numbchrototl)) 
    listchrollik = listchrollik.reshape((gdat.numbsweptotl, gdat.numbchrollik)) 
    listaccp = listaccp.flatten()
    listindxparamodi = listindxparamodi.flatten()
    
    gdat.rtag = retr_rtag(gdat, None)

    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(gdat.datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    if gdat.verbtype > 0:
        print 'Estimating the Bayesian evidence...'
        tim0 = time.time()
    
    listsamp = listsamp.reshape(gdat.numbsamptotl, -1)
    
    ## get an ellipse 
    gdat.elpscntr = percentile(listsamp, 50., axis=0)
    def retr_elpsfrac(elpsaxis):
        distnorm = sum(((listsamp - gdat.elpscntr[None, :]) / elpsaxis[None, :])**2, axis=1)
        indxsampregu = where(distnorm < 1.)[0]
        thissampfrac = indxsampregu.size / gdat.numbsamp
        vari = (thissampfrac / 0.05 - 1.)**2
        return vari
    thissamp = rand(gdat.numbpara) * 1e-6
    stdvpara = ones(gdat.numbpara) * 1e-6
    limtpara = zeros((2, gdat.numbpara))
    limtpara[1, :] = 1.
    elpsaxis, minmfunc = tdpy.util.minm(thissamp, retr_elpsfrac, stdvpara=stdvpara, limtpara=limtpara, tolrfunc=1e-6, verbtype=gdat.verbtype, optiprop=True)
    lnorregu = -0.5 * gdat.numbpara * log(pi) + sp.special.gammaln(0.5 * gdat.numbpara + 1.) - sum(elpsaxis)
    
    # temp
    #levi = lnorregu - log(mean(1. / exp(listllik[indxsampregu] - minmlistllik))) + minmlistllik
    #if not isfinite(levi):
    #    levi = 0.
    levi = retr_levi(listllik)
    gridchan.append(levi)
  
    # relative entropy
    info = retr_info(listllik, levi)
    gridchan.append(info)

    # collect posterior samples from the processes
    ## number of PS
    listnumbpnts = listsampvarb[:, :, gdat.indxsampnumbpnts].astype(int).reshape(gdat.numbsamptotl, -1)
    ## FDF normalization
    listfdfnnorm = listsampvarb[:, :, gdat.indxsampfdfnnorm].reshape(gdat.numbsamptotl, -1)
    ## FDF shape
    if gdat.fdfntype == 'powr':
        listfdfnslop = listsampvarb[:, :, gdat.indxsampfdfnslop].reshape(gdat.numbsamptotl, gdat.numbpopl)
    if gdat.fdfntype == 'brok':
        listfdfnbrek = listsampvarb[:, :, gdat.indxsampfdfnbrek].reshape(gdat.numbsamptotl, gdat.numbpopl)
        listfdfnsloplowr = listsampvarb[:, :, gdat.indxsampfdfnsloplowr].reshape(gdat.numbsamptotl, gdat.numbpopl)
        listfdfnslopuppr = listsampvarb[:, :, gdat.indxsampfdfnslopuppr].reshape(gdat.numbsamptotl, gdat.numbpopl)
    ## PSF parameters
    listpsfipara = listsampvarb[:, :, gdat.indxsamppsfipara].reshape(gdat.numbsamptotl, -1)
    ## Background normalization
    listnormback = listsampvarb[:, :, gdat.indxsampnormback].reshape(gdat.numbsamptotl, gdat.numbback, gdat.numbener)
    
    ## PS parameters
    listlgal = []
    listbgal = []
    listspec = []
    listsind = []
    listgang = []
    listaang = []
    for l in gdat.indxpopl:
        listlgal.append(zeros((gdat.numbsamptotl, gdat.maxmnumbpnts[l])))
        listbgal.append(zeros((gdat.numbsamptotl, gdat.maxmnumbpnts[l])))
        listspec.append(zeros((gdat.numbsamptotl, gdat.numbener, gdat.maxmnumbpnts[l])))
        listsind.append(zeros((gdat.numbsamptotl, gdat.maxmnumbpnts[l])))
        listgang.append(zeros((gdat.numbsamptotl, gdat.maxmnumbpnts[l])))
        listaang.append(zeros((gdat.numbsamptotl, gdat.maxmnumbpnts[l])))
    
    ## binned PS parameters
    listlgalhist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numblgal))
    listbgalhist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbbgal))
    listspechist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbflux, gdat.numbener))
    listsindhist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbsind))
    listganghist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbgang))
    listaanghist = empty((gdat.numbsamptotl, gdat.numbpopl, gdat.numbaang))
    for k in gdat.indxproc:
        for j in gdat.indxsamp:            
            n = k * gdat.numbsamp + j
            indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, listindxpntsfull[k][j])
            for l in gdat.indxpopl:
                numbpnts = indxsamplgal[l].size
                listlgal[l][n, 0:numbpnts] = listsampvarb[j, k, indxsamplgal[l]]
                listbgal[l][n, 0:numbpnts] = listsampvarb[j, k, indxsampbgal[l]]
                listspec[l][n, :, 0:numbpnts] = listsampvarb[j, k, indxsampspec[l]]
                listsind[l][n, 0:numbpnts] = listsampvarb[j, k, indxsampsind[l]]
                listgang[l][n, 0:numbpnts] = retr_gang(listlgal[l][n, 0:numbpnts], listbgal[l][n, 0:numbpnts])
                listaang[l][n, 0:numbpnts] = retr_aang(listlgal[l][n, 0:numbpnts], listbgal[l][n, 0:numbpnts])
                listlgalhist[n, l, :] = histogram(listlgal[l][n, 0:numbpnts], gdat.binslgal)[0]
                listbgalhist[n, l, :] = histogram(listbgal[l][n, 0:numbpnts], gdat.binsbgal)[0]
                for i in gdat.indxener:
                    listspechist[n, l, :, i] = histogram(listspec[l][n][i, :], gdat.binsspec[i, :])[0]
                listsindhist[n, l, :] = histogram(listsind[l][n], gdat.binssind)[0]
                listganghist[n, l, :] = histogram(listgang[l][n], gdat.binsgang)[0]
                listaanghist[n, l, :] = histogram(listaang[l][n], gdat.binsaang)[0]

    # auxiliary variables
    listpntsfluxmean = listpntsfluxmean.reshape(gdat.numbsamptotl, gdat.numbener)
   
    if gdat.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        tim0 = time.time()

    # binned and stacked posterior
    pntsprob = zeros((gdat.numbpopl, gdat.numbener, gdat.numbpixl, gdat.numbflux))
    for k in range(gdat.numbsamp):
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for h in range(gdat.numbflux):
                    indxpnts = where((gdat.binsspec[i, h] < listspec[l][k, i, :]) & (listspec[l][k, i, :] < gdat.binsspec[i, h+1]))[0]
                    hpixl = retr_indxpixl(gdat, listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    # Gelman-Rubin test
    if gdat.verbtype > 0:
        print 'Performing Gelman-Rubin convergence test...'
        tim0 = time.time()
    gmrbstat = zeros(gdat.numbpixlsave)
    for n in range(gdat.numbpixlsave):
        gmrbstat[n] = tdpy.mcmc.gmrb_test(listmodlcnts[:, :, n])

    # write the PCAT output to disc
    pathpcatlite = os.environ["PCAT_DATA_PATH"] + '/pcatlite_' + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.fits'  
    pathpcat = os.environ["PCAT_DATA_PATH"] + '/pcat_' + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.fits'  
    
    head = pf.Header()
    head['numbener'] = (gdat.numbener, 'Number of energy bins')
    head['numbevtt'] = (gdat.numbevtt, 'Number of PSF class bins')
    head['numbpopl'] = (gdat.numbpopl, 'Number of PS population')
    head['numbpsfipara'] = gdat.numbpsfipara
    head['numbformpara'] = gdat.numbformpara
    
    head['numbsamp'] = gdat.numbsamp
    head['numbburn'] = gdat.numbburn
    head['numbswep'] = gdat.numbswep
    head['factthin'] = gdat.factthin
    head['numbpopl'] = gdat.numbpopl
    head['numbproc'] = gdat.numbproc
    
    head['maxmgang'] = gdat.maxmgang
    head['lgalcntr'] = gdat.lgalcntr
    head['bgalcntr'] = gdat.bgalcntr
    head['numbflux'] = gdat.numbflux
    if gdat.pixltype == 'heal':
        head['numbsideheal'] = gdat.numbsideheal
    if gdat.pixltype == 'cart':
        head['numbsidecart'] = gdat.numbsidecart
    
    # axes
    head['minmlgal'] = gdat.minmlgal
    head['maxmlgal'] = gdat.maxmlgal
    head['minmbgal'] = gdat.minmbgal
    head['maxmbgal'] = gdat.maxmbgal
    
    head['minmflux'] = gdat.minmflux
    head['maxmflux'] = gdat.maxmflux
    head['minmsind'] = gdat.minmsind
    head['maxmsind'] = gdat.maxmsind
    
    head['datatype'] = gdat.datatype
    head['regitype'] = gdat.regitype
    head['psfntype'] = gdat.psfntype
    head['exprtype'] = gdat.exprtype
    head['pixltype'] = gdat.pixltype
    head['fdfntype'] = gdat.fdfntype
    
    head['trueinfo'] = gdat.trueinfo
    head['margsize'] = gdat.margsize
    head['strgcnfg'] = gdat.strgcnfg
    head['strgtime'] = gdat.strgtime
   
    ## proposal scales
    ### parameter updates
    head['stdvfdfnnorm'] = gdat.stdvfdfnnorm
    head['stdvfdfnslop'] = gdat.stdvfdfnslop
    head['stdvpsfipara'] = gdat.stdvpsfipara
    head['stdvback'] = gdat.stdvback
    head['stdvlbhl'] = gdat.stdvlbhl
    head['stdvflux'] = gdat.stdvflux
    head['stdvsind'] = gdat.stdvsind
    ### relative size of the tail of the Gaussian proposals
    head['fracrand'] = gdat.fracrand
    ### split and merge
    head['radispmrlbhl'] = gdat.radispmrlbhl
    head['stdvspmrsind'] = gdat.stdvspmrsind

    ## index of the energy bin, where the FDF is defined
    head['indxenerfdfn'] = gdat.indxenerfdfn[0]

    if gdat.datatype == 'mock':
        head['mockfdfntype'] = gdat.mockfdfntype
        head['mockpsfntype'] = gdat.mockpsfntype

    head['strgexpr'] = gdat.strgexpr
    for k in gdat.indxback:
        head['strgback%04d' % k] = gdat.strgback[k]
        head['lablback%04d' % k] = gdat.lablback[k]
        head['nameback%04d' % k] = gdat.nameback[k]
    head['strgexpo'] = gdat.strgexpo
    
    head['levi'] = levi
    head['info'] = info
    
    listhdun = []
    listhdun.append(pf.PrimaryHDU(header=head))

    # posterior
    ## primaries
    for l in gdat.indxpopl:
        listhdun.append(pf.ImageHDU(listlgal[l]))
        listhdun[-1].header['EXTNAME'] = 'lgalpop%d' % l
        listhdun.append(pf.ImageHDU(listbgal[l]))
        listhdun[-1].header['EXTNAME'] = 'bgalpop%d' % l
        listhdun.append(pf.ImageHDU(listspec[l]))
        listhdun[-1].header['EXTNAME'] = 'specpop%d' % l
        listhdun.append(pf.ImageHDU(listsind[l]))
        listhdun[-1].header['EXTNAME'] = 'sindpop%d' % l
        listhdun.append(pf.ImageHDU(listgang[l]))
        listhdun[-1].header['EXTNAME'] = 'gangpop%d' % l
        listhdun.append(pf.ImageHDU(listaang[l]))
        listhdun[-1].header['EXTNAME'] = 'aangpop%d' % l

        ## save the posterior positions as a CSV file
        path = os.environ["PCAT_DATA_PATH"] + '/pcat_lgalpop%d_' % l + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.csv'  
        savetxt(path, listlgal[l], fmt='%7.5g', delimiter=',')
        path = os.environ["PCAT_DATA_PATH"] + '/pcat_bgalpop%d_' % l + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.csv'  
        savetxt(path, listbgal[l], fmt='%7.5g', delimiter=',')


    listhdun.append(pf.ImageHDU(listfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'fdfnnorm'
   
    if gdat.fdfntype == 'powr':
        listhdun.append(pf.ImageHDU(listfdfnslop))
        listhdun[-1].header['EXTNAME'] = 'fdfnslop'
    
    if gdat.fdfntype == 'brok':
        listhdun.append(pf.ImageHDU(listfdfnbrek))
        listhdun[-1].header['EXTNAME'] = 'fdfnbrek'
    
        listhdun.append(pf.ImageHDU(listfdfnsloplowr))
        listhdun[-1].header['EXTNAME'] = 'fdfnsloplowr'
    
        listhdun.append(pf.ImageHDU(listfdfnslopuppr))
        listhdun[-1].header['EXTNAME'] = 'fdfnslopuppr'
    
    listhdun.append(pf.ImageHDU(listpsfipara))
    listhdun[-1].header['EXTNAME'] = 'psfipara'
    
    listhdun.append(pf.ImageHDU(listnormback))
    listhdun[-1].header['EXTNAME'] = 'normback'
   
    ## log-likelihood
    listhdun.append(pf.ImageHDU(listllik))
    listhdun[-1].header['EXTNAME'] = 'llik'
    
    ## log-prior
    listhdun.append(pf.ImageHDU(listlpri))
    listhdun[-1].header['EXTNAME'] = 'lpri'
   
    # store the lite file
    pf.HDUList(listhdun).writeto(pathpcatlite, clobber=True, output_verify='ignore')

    ## unit sample
    listhdun.append(pf.ImageHDU(listsamp))
    listhdun[-1].header['EXTNAME'] = 'listsamp'
    
    ## number of PS
    listhdun.append(pf.ImageHDU(listnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'numbpnts'

    ## mean PS flux
    listhdun.append(pf.ImageHDU(listpntsfluxmean))
    listhdun[-1].header['EXTNAME'] = 'listpntsfluxmean'
    
    ## binned PS parameters
    listhdun.append(pf.ImageHDU(listlgalhist))
    listhdun[-1].header['EXTNAME'] = 'lgalhist'
    
    listhdun.append(pf.ImageHDU(listbgalhist))
    listhdun[-1].header['EXTNAME'] = 'bgalhist'
    
    listhdun.append(pf.ImageHDU(listganghist))
    listhdun[-1].header['EXTNAME'] = 'ganghist'
    
    listhdun.append(pf.ImageHDU(listaanghist))
    listhdun[-1].header['EXTNAME'] = 'aanghist'
    
    listhdun.append(pf.ImageHDU(listspechist))
    listhdun[-1].header['EXTNAME'] = 'spechist'
    
    listhdun.append(pf.ImageHDU(listsindhist))
    listhdun[-1].header['EXTNAME'] = 'sindhist'
    
    ## stored model counts in random pixels
    listhdun.append(pf.ImageHDU(listmodlcnts))
    listhdun[-1].header['EXTNAME'] = 'modlcnts'
   
    ## proposal type
    listhdun.append(pf.ImageHDU(listindxprop))
    listhdun[-1].header['EXTNAME'] = 'indxprop'
   
    ## section timings
    listhdun.append(pf.ImageHDU(listchrototl))
    listhdun[-1].header['EXTNAME'] = 'listchrototl'
    
    ## likelihood timings
    listhdun.append(pf.ImageHDU(listchrollik))
    listhdun[-1].header['EXTNAME'] = 'listchrollik'
    
    ## acceptance
    listhdun.append(pf.ImageHDU(listaccp))
    listhdun[-1].header['EXTNAME'] = 'accp'
    
    ## index of the modified parameter
    listhdun.append(pf.ImageHDU(listindxparamodi))
    listhdun[-1].header['EXTNAME'] = 'indxparamodi'
    
    ## split and merge diagnostics
    listhdun.append(pf.ImageHDU(gdat.listauxipara))
    listhdun[-1].header['EXTNAME'] = 'auxipara'
   
    listhdun.append(pf.ImageHDU(gdat.listlaccfrac))
    listhdun[-1].header['EXTNAME'] = 'laccfrac'
    
    listhdun.append(pf.ImageHDU(gdat.listnumbpair))
    listhdun[-1].header['EXTNAME'] = 'numbpair'
    
    listhdun.append(pf.ImageHDU(gdat.listjcbnfact))
    listhdun[-1].header['EXTNAME'] = 'jcbnfact'
    
    listhdun.append(pf.ImageHDU(gdat.listcombfact))
    listhdun[-1].header['EXTNAME'] = 'combfact'
    
    # priors
    ## prior limits
    listhdun.append(pf.ImageHDU(gdat.minmpsfipara))
    listhdun[-1].header['EXTNAME'] = 'minmpsfipara'

    listhdun.append(pf.ImageHDU(gdat.maxmpsfipara))
    listhdun[-1].header['EXTNAME'] = 'maxmpsfipara'

    listhdun.append(pf.ImageHDU(gdat.minmnormback))
    listhdun[-1].header['EXTNAME'] = 'minmnormback'

    listhdun.append(pf.ImageHDU(gdat.maxmnormback))
    listhdun[-1].header['EXTNAME'] = 'maxmnormback'

    ## hyperprior limits
    
    listhdun.append(pf.ImageHDU(gdat.minmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'minmfdfnnorm'
    
    listhdun.append(pf.ImageHDU(gdat.maxmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'maxmfdfnnorm'
    
    if gdat.fdfntype == 'powr':
        listhdun.append(pf.ImageHDU(gdat.minmfdfnslop))
        listhdun[-1].header['EXTNAME'] = 'minmfdfnslop'
    
        listhdun.append(pf.ImageHDU(gdat.maxmfdfnslop))
        listhdun[-1].header['EXTNAME'] = 'maxmfdfnslop'
    
    if gdat.fdfntype == 'brok':
        listhdun.append(pf.ImageHDU(gdat.minmfdfnbrek))
        listhdun[-1].header['EXTNAME'] = 'minmfdfnbrek'
    
        listhdun.append(pf.ImageHDU(gdat.maxmfdfnbrek))
        listhdun[-1].header['EXTNAME'] = 'maxmfdfnbrek'
    
        listhdun.append(pf.ImageHDU(gdat.minmfdfnsloplowr))
        listhdun[-1].header['EXTNAME'] = 'minmfdfnsloplowr'
    
        listhdun.append(pf.ImageHDU(gdat.maxmfdfnsloplowr))
        listhdun[-1].header['EXTNAME'] = 'maxmfdfnsloplowr'
    
        listhdun.append(pf.ImageHDU(gdat.minmfdfnslopuppr))
        listhdun[-1].header['EXTNAME'] = 'minmfdfnslopuppr'
    
        listhdun.append(pf.ImageHDU(gdat.maxmfdfnslopuppr))
        listhdun[-1].header['EXTNAME'] = 'maxmfdfnslopuppr'
    
    listhdun.append(pf.ImageHDU(gdat.meansdfn))
    listhdun[-1].header['EXTNAME'] = 'meansdfn'
    
    listhdun.append(pf.ImageHDU(gdat.stdvsdfn))
    listhdun[-1].header['EXTNAME'] = 'stdvsdfn'
    
    listhdun.append(pf.ImageHDU(gdat.binsener))
    listhdun[-1].header['EXTNAME'] = 'binsener'

    # sampler setup
    ## maximum number of point sources
    listhdun.append(pf.ImageHDU(gdat.maxmnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'maxmnumbpnts'
    
    ## indices of energy bins included in the analysis
    listhdun.append(pf.ImageHDU(gdat.indxenerincl))
    listhdun[-1].header['EXTNAME'] = 'indxenerincl'
    
    ## indices of PSF classes included in the analysis
    listhdun.append(pf.ImageHDU(gdat.indxevttincl))
    listhdun[-1].header['EXTNAME'] = 'indxevttincl'
    
    ## proposal frequencies
    listhdun.append(pf.ImageHDU(gdat.probprop))
    listhdun[-1].header['EXTNAME'] = 'probprop'
    
    # processed output products
    ## convergence diagnostic
    listhdun.append(pf.ImageHDU(gmrbstat))
    listhdun[-1].header['EXTNAME'] = 'gmrbstat'
    
    ## PDF of PS positions
    listhdun.append(pf.ImageHDU(pntsprob))
    listhdun[-1].header['EXTNAME'] = 'pntsprob'
   
    # mock data
    # truth information
    if gdat.datatype == 'mock':
        listhdun.append(pf.ImageHDU(gdat.truenumbpnts))
        listhdun[-1].header['EXTNAME'] = 'mocknumbpnts'

        for l in gdat.indxpopl:
            listhdun.append(pf.ImageHDU(gdat.truelgal[l]))
            listhdun[-1].header['EXTNAME'] = 'mocklgalpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truebgal[l]))
            listhdun[-1].header['EXTNAME'] = 'mockbgalpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truespec[l]))
            listhdun[-1].header['EXTNAME'] = 'mockspecpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truesind[l]))
            listhdun[-1].header['EXTNAME'] = 'mocksindpop%d' % l

        if gdat.mockfdfntype == 'powr':
            listhdun.append(pf.ImageHDU(gdat.mockfdfnslop))
            listhdun[-1].header['EXTNAME'] = 'mockfdfnslop'

        if gdat.mockfdfntype == 'brok':
            listhdun.append(pf.ImageHDU(gdat.mockfdfnbrek))
            listhdun[-1].header['EXTNAME'] = 'mockfdfnbrek'

            listhdun.append(pf.ImageHDU(gdat.mockfdfnsloplowr))
            listhdun[-1].header['EXTNAME'] = 'mockfdfnsloplowr'

            listhdun.append(pf.ImageHDU(gdat.mockfdfnslopuppr))
            listhdun[-1].header['EXTNAME'] = 'mockfdfnslopuppr'

        listhdun.append(pf.ImageHDU(gdat.truenormback))
        listhdun[-1].header['EXTNAME'] = 'mocknormback'

        listhdun.append(pf.ImageHDU(gdat.truepsfipara))
        listhdun[-1].header['EXTNAME'] = 'mockpsfipara'
        
    pf.HDUList(listhdun).writeto(pathpcat, clobber=True, output_verify='ignore')

    if gdat.makeplot:
        plot_post(pathpcat)

    timetotlreal = time.time() - timetotlreal
    timetotlproc = time.clock() - timetotlproc
    
    if gdat.verbtype > 0:
        for k in gdat.indxproc:
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, timereal[k], timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (timetotlreal, sum(timeproc))
        print 'The ensemble of catalogs is at ' + pathpcat
        if gdat.makeplot:
            print 'The plots are in ' + gdat.pathplot
        
    return gridchan
    
    
def plot_samp(gdat, gdatmodi):

    gdatmodi.thisresicnts = gdat.datacnts - gdatmodi.thismodlcnts
    
    gdatmodi.thispsfn = retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.angldisp, gdat.psfntype)
    if gdat.pixltype == 'cart':
        gdatmodi.thispsfn = gdatmodi.thispsfn.reshape((gdat.numbener, -1, gdat.numbevtt))
    gdatmodi.thisfwhm = 2. * retr_psfnwdth(gdat, gdatmodi.thispsfn, 0.5)

    plot_psfn(gdat, gdatmodi)

    gdatmodi.thisbackcntsmean = empty((gdat.numbener, gdat.numbevtt))
    for c in gdat.indxback:
        gdatmodi.thisbackcntsmean += mean(gdatmodi.thissampvarb[gdat.indxsampnormback[c, :, None, None]] * gdat.backflux[c] * gdat.expo * \
            gdat.diffener[:, None, None] * pi * gdatmodi.thisfwhm[:, None, :]**2 / 4., 1)

    # temp -- list may not be the ultimate solution to copy gdatmodi.thisindxpntsfull
    gdat.temppntsflux, gdat.temppntscnts, gdat.tempmodlflux, gdat.tempmodlcnts = retr_maps(gdat, list(gdatmodi.thisindxpntsfull), copy(gdatmodi.thissampvarb))
    gdatmodi.thispntscnts = gdatmodi.thispntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    gdatmodi.thiserrrpnts = 100. * (gdatmodi.thispntscnts - gdat.temppntscnts) / gdat.temppntscnts

    gdatmodi.thiscnts = []
    gdatmodi.indxtruepntsassc = []
    for l in gdat.indxpopl:
        indxpixltemp = retr_indxpixl(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]], gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]])
        cntstemp = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]][:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]
        gdatmodi.thiscnts.append(cntstemp)
        if gdat.trueinfo:
            lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]]
            bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]]
            spec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]]
            indxmodl, indxtruepntsassc = pair_catl(gdat, l, lgal, bgal, spec)
            gdatmodi.indxtruepntsassc.append(indxtruepntsassc)
            gdatmodi.thisspecmtch = copy(spec[:, indxmodl])
            gdatmodi.thisspecmtch[:, indxtruepntsassc.miss] = 0.
            plot_scatspec(gdat, l, gdatmodi=gdatmodi)
        
        plot_histspec(gdat, l, gdatmodi=gdatmodi)
        plot_histsind(gdat, l, gdatmodi=gdatmodi)
        plot_histcnts(gdat, l, gdatmodi=gdatmodi)
        plot_compfrac(gdat, gdatmodi=gdatmodi)

    for i in gdat.indxener:
        
        plot_datacnts(gdat, l, gdatmodi, i, None)
        plot_modlcnts(gdat, l, gdatmodi, i, None)
        plot_resicnts(gdat, l, gdatmodi, i, None)
        plot_errrpnts(gdat, gdatmodi, i, None)

        # temp
        #plot_scatpixl(gdat, gdatmodi=gdatmodi)
        #plot_catlfram(gdat, l, gdatmodi, i, None)
        #for m in gdat.indxevtt:
        #    plot_datacnts(gdat, gdatmodi, i, m)
        #    plot_catl(gdat, gdatmodi, i, m, thiscnts)
        #    plot_modlcnts(gdat, gdatmodi, i, m)
        #    plot_resicnts(gdat, gdatmodi, i, m)
        

    #if gdat.numbener == 3:
    #    plot_datacnts(gdat, gdatmodi, None, None)
        
    #plot_fwhm(gdat, gdatmodi)
    #plot_backcntsmean(gdat, gdatmodi)
    
    if amax(abs(gdatmodi.thiserrrpnts)) > 0.1 and False:
        print 'Approximation error went above the limit!'
    
    
def rjmc(gdat, gdatmodi, indxprocwork):

    # sweeps to be saved
    boolsave = zeros(gdat.numbswep, dtype=bool)
    indxswepsave = arange(gdat.numbburn, gdat.numbswep, gdat.factthin)
    boolsave[indxswepsave] = True
    
    indxsampsave = zeros(gdat.numbswep, dtype=int)
    indxsampsave[indxswepsave] = arange(gdat.numbsamp)

    listsamp = zeros((gdat.numbsamp, gdat.numbpara)) + -1.
    listsampvarb = zeros((gdat.numbsamp, gdat.numbpara)) + -1.
    listindxprop = zeros(gdat.numbswep)
    listchrototl = zeros((gdat.numbswep, gdat.numbchrototl))
    gdatmodi.listchrollik = zeros((gdat.numbswep, gdat.numbchrollik))
    listllik = zeros(gdat.numbsamp)
    listlpri = zeros(gdat.numbsamp)
    listlprinorm = zeros(gdat.numbsamp)
    listaccp = zeros(gdat.numbswep, dtype=bool)
    listaccpspec = []
    listindxparamodi = zeros(gdat.numbswep, dtype=int)
    listmodlcnts = zeros((gdat.numbsamp, gdat.numbpixlsave))
    listpntsfluxmean = zeros((gdat.numbsamp, gdat.numbener))
    listindxpntsfull = []
    
    gdat.listauxipara = zeros((gdat.numbswep, gdat.numbcompcolr))
    gdat.listlaccfrac = zeros(gdat.numbswep)
    gdat.listnumbpair = zeros(gdat.numbswep)
    gdat.listjcbnfact = zeros(gdat.numbswep)
    gdat.listcombfact = zeros(gdat.numbswep)

    gdat.cntrswep = 0
    
    # initialize the chain
    retr_llik(gdat, gdatmodi, init=True)
    retr_lpri(gdat, gdatmodi, init=True)

    # current sample index
    thiscntr = -1
    
    while gdat.cntrswep < gdat.numbswep:
        
        timeinit = time.time()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdat.cntrswep

        thismakefram = (gdat.cntrswep % gdat.numbswepplot == 0) and indxprocwork == int(float(gdat.cntrswep) / gdat.numbswep * gdat.numbproc) and gdat.makeplot
        gdatmodi.boolreje = False
    
        # choose a proposal type
        retr_thisindxprop(gdat, gdatmodi)
            
        # save the proposal type
        listindxprop[gdat.cntrswep] = gdatmodi.thisindxprop
        if gdat.verbtype > 1:
            print 'thisindxprop: ', gdat.strgprop[gdatmodi.thisindxprop]
        
        if gdat.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop(gdat, gdatmodi)
        timefinl = time.time()
        listchrototl[gdat.cntrswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                gdat.lock.acquire()
            print 'Process %d started making a frame.' % indxprocwork
            plot_samp(gdat, gdatmodi)
            print 'Process %d finished making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                gdat.lock.release()
            
        # reject the sample if proposal is outside the prior
        if gdat.fdfntype == 'powr':
            if gdatmodi.thisindxprop == gdat.indxpropfdfnslop:
                if gdatmodi.drmcsamp[gdat.indxsampvarbmodi, 1] < 0. or gdatmodi.drmcsamp[gdat.indxsampvarbmodi, 1] > 1.:
                    gdatmodi.boolreje = True
        if gdat.fdfntype == 'brok':
            if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek or gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr or gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
                if gdatmodi.drmcsamp[gdat.indxsampvarbmodi, 1] < 0. or gdatmodi.drmcsamp[gdat.indxsampvarbmodi, 1] > 1.:
                    gdatmodi.boolreje = True
        if (gdatmodi.thisindxprop == gdat.indxpropfdfnnorm or \
                    gdatmodi.thisindxprop >= gdat.indxproppsfipara and gdatmodi.thisindxprop != gdat.indxpropbrth and gdatmodi.thisindxprop != gdat.indxpropdeth) and not gdatmodi.boolreje:
            if where((gdatmodi.drmcsamp[gdat.indxsampmodi, 1] < 0.) | (gdatmodi.drmcsamp[gdat.indxsampmodi, 1] > 1.))[0].size > 0:
                gdatmodi.boolreje = True

        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            if gdat.psfntype == 'doubking':
                if gdatmodi.nextsampvarb[gdat.indxsamppsfipara[1]] >= gdatmodi.nextsampvarb[gdat.indxsamppsfipara[3]]:
                    gdatmodi.boolreje = True
            elif gdat.psfntype == 'doubgaus':
                if gdatmodi.nextsampvarb[gdat.indxsamppsfipara[1]] >= gdatmodi.nextsampvarb[gdat.indxsamppsfipara[2]]:
                    gdatmodi.boolreje = True
            
        if not gdatmodi.boolreje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri(gdat, gdatmodi)
            timefinl = time.time()
            listchrototl[gdat.cntrswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(gdat, gdatmodi) 
            timefinl = time.time()
            listchrototl[gdat.cntrswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(gdatmodi.deltllik + gdatmodi.deltlpri + gdatmodi.laccfrac)

            if gdat.verbtype > 1:
                print 'deltlpri'
                print gdatmodi.deltlpri
                print 'deltllik'
                print gdatmodi.deltllik
                print 'laccfrac'
                print gdatmodi.laccfrac
                print
        else:
            accpprob = 0.
    
        # accept the sample
        if accpprob >= rand():

            if gdat.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(gdat, gdatmodi)

            listaccp[gdat.cntrswep] = True

        # reject the sample
        else:

            if gdat.verbtype > 1:
                print 'Rejected.'

            listaccp[gdat.cntrswep] = False
             
        # sanity checks
        indxsampbadd = where((gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.) | (gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.))[0] + 1
        if indxsampbadd.size > 0:
            print 'Unit sample vector went outside [0,1]!'
            print 'cntrswep'
            print gdat.cntrswep
            print 'thisindxprop'
            print gdat.strgprop[gdatmodi.thisindxprop]
            print 'indxsampbadd'
            print indxsampbadd
            print 'drmcsamp'
            print gdatmodi.drmcsamp[indxsampbadd, :]
            return
            
        for l in gdat.indxpopl:
            indxtemp = where(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], :]] < gdat.minmflux)[0]
            if indxtemp.size > 0:
                print 'Spectrum of some PS went below the prior range!'
                print 'indxtemp'
                print indxtemp
                print 'flux'
                print gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], indxtemp]]
                print 
            indxtemp = where(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], :]] > gdat.maxmflux)[0]
            if indxtemp.size > 0:
                print 'Spectrum of some PS went above the prior range!'          
                print 'indxtemp'
                print indxtemp
                print 'flux'
                print gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], indxtemp]]
                print 

        # save the sample
        if boolsave[gdat.cntrswep]:
            listsamp[indxsampsave[gdat.cntrswep], :] = gdatmodi.drmcsamp[:, 0]
            listsampvarb[indxsampsave[gdat.cntrswep], :] = gdatmodi.thissampvarb
            listmodlcnts[indxsampsave[gdat.cntrswep], :] = gdatmodi.thismodlcnts[0, gdat.indxpixlsave, 0]
            listpntsfluxmean[indxsampsave[gdat.cntrswep], :] = mean(sum(gdatmodi.thispntsflux * gdat.expo, 2) / sum(gdat.expo, 2), 1)
            listindxpntsfull.append(copy(gdatmodi.thisindxpntsfull))
            listllik[indxsampsave[gdat.cntrswep]] = sum(gdatmodi.thisllik)
            listlpri[indxsampsave[gdat.cntrswep]] = sum(gdatmodi.thislpri)
            
            lprinorm = 0.
            for l in gdat.indxpopl:
                # temp
                ## brok terms are not complete
                ## lpri calculation is turned off
                break
                numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[l]]
                fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[l]]
                lpri += numbpnts * gdat.priofactlgalbgal + gdat.priofactfdfnslop + gdat.priofactfdfnnorm - log(fdfnnorm)
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], :]]
                if gdat.fdfntype == 'powr':
                    fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[l]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_flux_powr(gdat, flux, fdfnslop)))
                if gdat.fdfntype == 'brok':
                    fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[l]]
                    fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[l]]
                    fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[l]]
                    lpri += sum(log(pdfn_flux_brok(gdat, flux, fdfnbrek, fdfnsloplowr, fdfnslopuppr)))
            listlprinorm[indxsampsave[gdat.cntrswep]] = lprinorm
            
            if gdat.tracsamp:
                
                numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[gdat.cntrswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            timefinl = time.time()
            listchrototl[gdat.cntrswep, 0] = timefinl - timeinit

        # log the progress
        if gdat.verbtype > 0:
            thiscntr = tdpy.util.show_prog(gdat.cntrswep, gdat.numbswep, thiscntr, indxprocwork=indxprocwork)
    
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
        gdat.cntrswep += 1

    if gdat.verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    gdatmodi.listchrollik = array(gdatmodi.listchrollik)
    
    listchan = [listsamp, listsampvarb, listindxprop, listchrototl, listllik, listlpri, listaccp, listmodlcnts, listindxpntsfull, listindxparamodi, \
        gdat.listauxipara, gdat.listlaccfrac, gdat.listnumbpair, gdat.listjcbnfact, gdat.listcombfact, \
        listpntsfluxmean, gdatmodi.listchrollik]
    
    return listchan

