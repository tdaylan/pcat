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
    
    # construct the run tag
    gdat.rtag = retr_rtag(gdat, indxprocwork)
    
    # initialize the sample vector 
    if gdat.randinit or not gdat.trueinfo:
        if gdat.initnumbpnts != None:
            thisnumbpnts = gdat.initnumbpnts
        else:
            thisnumbpnts = empty(gdat.numbpopl, dtype=int)
            for l in gdat.indxpopl:
                thisnumbpnts[l] = choice(arange(gdat.minmnumbpnts, gdat.maxmnumbpnts[l] + 1))
    else:
        thisnumbpnts = gdat.truenumbpnts
        
    gdat.thisindxpntsfull = []
    gdat.thisindxpntsempt = []
    for l in gdat.indxpopl:
        gdat.thisindxpntsfull.append(range(thisnumbpnts[l]))
        gdat.thisindxpntsempt.append(range(thisnumbpnts[l], gdat.maxmnumbpnts[l]))
    gdat.thisindxsamplgal, gdat.thisindxsampbgal, gdat.thisindxsampspec, \
        gdat.thisindxsampsind, gdat.thisindxsampcomp = retr_indx(gdat, gdat.thisindxpntsfull)
   
    gdat.deltllik = 0.

    gdat.drmcsamp = zeros((gdat.maxmsampsize, 2))
    
    gdat.drmcsamp[gdat.indxsampnumbpnts, 0] = thisnumbpnts
    gdat.drmcsamp[gdat.indxsampfdfnnorm, 0] = rand(gdat.numbpopl)
    if gdat.trueinfo and gdat.datatype == 'mock':
        for l in gdat.indxpopl:
            if gdat.fdfntype == 'brok':
                gdat.drmcsamp[gdat.indxsampfdfnsloplowr[l], 0] = cdfn_atan(gdat.truefdfnsloplowr[l], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
                gdat.drmcsamp[gdat.indxsampfdfnslopuppr[l], 0] = cdfn_atan(gdat.truefdfnslopuppr[l], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
                gdat.drmcsamp[gdat.indxsampfdfnbrek[l], 0] = cdfn_logt(gdat.truefdfnbrek[l], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
            else:
                gdat.drmcsamp[gdat.indxsampfdfnslop[l], 0] = cdfn_atan(gdat.truefdfnslop[l], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
    else:
        gdat.drmcsamp[gdat.indxsampfdfnslop, 0] = rand(gdat.numbpopl * gdat.numbener).reshape((gdat.numbpopl, gdat.numbener))
    gdat.drmcsamp[gdat.indxsampnormback, 0] = rand(gdat.numbback * gdat.numbener).reshape((gdat.numbback, gdat.numbener))
    if gdat.randinit or not gdat.trueinfo or gdat.truepsfipara == None:
        gdat.drmcsamp[gdat.indxsamppsfipara, 0] = retr_randunitpsfipara(gdat)
    else:
        if gdat.psfntype == gdat.truepsfntype:
            for k in gdat.indxmodlpsfipara:
                gdat.drmcsamp[gdat.indxsamppsfipara[k], 0] = cdfn_psfipara(gdat, gdat.truepsfipara[k], k)
        else:
            gdat.drmcsamp[gdat.indxsamppsfipara, 0] = retr_randunitpsfipara(gdat)
   
    if False:
        print
        print 'hey'
        for k in range(gdat.truepsfipara.size):
            if k % 5 == 0:
                print
            print gdat.truepsfipara[k]
        print
        print 'hey'
        for k in range(gdat.indxsamppsfipara.size):
            if k % 5 == 0:
                print
            print gdat.drmcsamp[gdat.indxsamppsfipara[k], 0]

    
    for l in gdat.indxpopl:
        if gdat.randinit or not gdat.trueinfo:
            gdat.drmcsamp[gdat.thisindxsampcomp[l], 0] = rand(gdat.thisindxsampcomp[l].size)
        else:
            gdat.drmcsamp[gdat.thisindxsamplgal[l], 0] = copy(cdfn_self(gdat.truelgal[l], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg))
            gdat.drmcsamp[gdat.thisindxsampbgal[l], 0] = copy(cdfn_self(gdat.truebgal[l], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg))
            for i in gdat.indxenerfdfn:
                if gdat.fdfntype == 'brok':
                    fdfnsloplowr = icdf_atan(gdat.drmcsamp[gdat.indxsampfdfnsloplowr[l, i], 0], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
                    fdfnslopuppr = icdf_atan(gdat.drmcsamp[gdat.indxsampfdfnslopuppr[l, i], 0], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
                    fdfnbrek = icdf_logt(gdat.drmcsamp[gdat.indxsampfdfnbrek[l, i], 0], gdat.minmfdfnbrek[l, i], gdat.factfdfnbrek[l, i])
                    spec = icdf_spec_brok(gdat, gdat.truespec[l][0, i, :], fdfnsloplowr, fdfnslopuppr, fdfnbrek, gdat.minmspec[i], gdat.maxmspec[i])
                else:
                    fdfnslop = icdf_atan(gdat.drmcsamp[gdat.indxsampfdfnslop[l, i], 0], gdat.minmfdfnslop, gdat.factfdfnslop)
                    spec = cdfn_spec(gdat, gdat.truespec[l][0, i, :], fdfnslop, gdat.minmspec[i], gdat.maxmspec[i])
                gdat.drmcsamp[gdat.thisindxsampspec[l][i, :], 0] = copy(spec)
            if gdat.colrprio:
                gdat.drmcsamp[gdat.thisindxsampsind[l], 0] = cdfn_eerr(gdat.truesind[l], gdat.meansind[l], gdat.stdvsind[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
    
    gdat.thissampvarb = retr_sampvarb(gdat, gdat.thisindxpntsfull, gdat.drmcsamp[:, 0])
    gdat.thispntsflux, gdat.thispntscnts, gdat.thismodlflux, gdat.thismodlcnts = retr_maps(gdat, gdat.thisindxpntsfull, gdat.thissampvarb)
    
    if gdat.verbtype > 2:
        print 'thisindxpntsfull: ', gdat.thisindxpntsfull
        print 'thisindxpntsempt: ', gdat.thisindxpntsempt  
        print 'thisindxsamplgal: ', gdat.thisindxsamplgal
        print 'thisindxsampbgal: ', gdat.thisindxsampbgal
        print 'thisindxsampspec: '
        print gdat.thisindxsampspec
        if gdat.colrprio:
            print 'thisindxsampsind: ', gdat.thisindxsampsind
        print 'thisindxsampcomp: ', gdat.thisindxsampcomp

    gdat.nextpntsflux = zeros_like(gdat.thispntsflux)
    gdat.nextmodlflux = zeros_like(gdat.thispntsflux)
    gdat.nextmodlcnts = zeros_like(gdat.thispntsflux)
    gdat.nextllik = zeros_like(gdat.thispntsflux)

    gdat.nextsampvarb = copy(gdat.thissampvarb)
    
    if gdat.verbtype > 1:
        print 'thissampvarb: ', gdat.thissampvarb
        
    # sampler setup
    # auxiliary variable standard gdat.deviation for merge/split
    gdat.maxmdistpnts = 2. # [deg]
 
    listchan = rjmc(gdat, indxprocwork)
    
    timereal = time.time() - timereal
    timeproc = time.clock() - timeproc
    
    listchan.append(timereal)
    listchan.append(timeproc)
    
    return listchan

    
def init( \
         verbtype=1, \
         plotperd=50000, \
         makeplot=True, \
         diagsamp=False, \
         numbswep=2000000, \
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
         meansind=None, \
         stdvsind=None, \
         minmfdfnnorm=None, \
         maxmfdfnnorm=None, \
         minmfdfnslop=None, \
         maxmfdfnslop=None, \
         fdfntype='powr', \
         psfntype=None, \
         proppsfn=True, \
         colrprio=True, \
         numbpopl=1, \
         indxevttincl=arange(2, 4), \
         indxenerincl=arange(5), \
         maxmnumbpnts=array([1000]), \
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
         stdvfdfnnorm=0.05, \
         stdvfdfnslop=0.1, \
         stdvpsfipara=0.1, \
         stdvback=0.04, \
         stdvlbhl=0.1, \
         stdvspec=0.15, \
         stdvpropsind=0.15, \
         fracrand=0.05, \
         mocknumbpnts=None, \
         mockfdfnslop=None, \
         mockfdfnsloplowr=None, \
         mockfdfnslopuppr=None, \
         mockfdfnbrek=None, \
         mocknormback=None, \
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
            minmsind = array([1.2])
        if maxmsind == None:
            maxmsind = array([3.2])
        if meansind == None:
            meansind = array([2.2])
        if stdvsind == None:
            stdvsind = array([0.3])
        if minmfdfnnorm == None:
            minmfdfnnorm = array([1e0])
        if maxmfdfnnorm == None:
            maxmfdfnnorm = array([1e2])
        if minmfdfnslop == None:
            minmfdfnslop = array([1.5])
        if maxmfdfnslop == None:
            maxmfdfnslop = array([2.5])
        if margsize == None:
            margsize = 1.
        if psfntype == None:
            psfntype = 'singking'
        if mockpsfntype == None:
            mockpsfntype = 'doubking'
    ## SDSS
    if exprtype == 'sdss':
        if maxmgang == None:
            maxmgang = 20. / 3600.
        if minmsind == None:
            minmsind = array([1.])
        if maxmsind == None:
            maxmsind = array([3.])
        if meansind == None:
            meansind = array([2.2])
        if stdvsind == None:
            stdvsind = array([0.5])
        if minmfdfnnorm == None:
            minmfdfnnorm = 1e0
        if maxmfdfnnorm == None:
            maxmfdfnnorm = 1e2
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
    gdat.plotperd = plotperd
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
    
    ## color priors
    gdat.colrprio = colrprio
   
    ## flag to turn off PSF parameter updates
    gdat.proppsfn = proppsfn

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
    ### flux distribution function normalization
    gdat.minmfdfnnorm = minmfdfnnorm
    gdat.maxmfdfnnorm = maxmfdfnnorm
    ### flux distribution function power law index
    gdat.minmfdfnslop = minmfdfnslop
    gdat.maxmfdfnslop = maxmfdfnslop
    ### flux
    gdat.minmspec = minmspec
    gdat.maxmspec = maxmspec
    ### spectral power-law index
    gdat.minmsind = minmsind
    gdat.maxmsind = maxmsind
    gdat.meansind = meansind
    gdat.stdvsind = stdvsind
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
    gdat.stdvspec = stdvspec
    gdat.stdvpropsind = stdvpropsind

    ### fraction of heavy-tailed proposals
    gdat.fracrand = fracrand
    
    ### maximum angle over which splits and merges are proposed
    gdat.spmrlbhl = spmrlbhl
    
    ## mock data setup
    if datatype == 'mock':
        ### mock PSF model type
        gdat.mockpsfntype = mockpsfntype
        ### mock FDF model type
        gdat.mockfdfntype = mockfdfntype
        ### mock number of PS
        gdat.mocknumbpnts = mocknumbpnts
        if gdat.fdfntype == 'brok':
            gdat.mockfdfnsloplowr = mockfdfnsloplowr
            gdat.mockfdfnslopuppr = mockfdfnslopuppr
            gdat.mockfdfnbrek = mockfdfnbrek
        else:
            gdat.mockfdfnslop = mockfdfnslop
        gdat.mocknormback = mocknormback
        ### flag to position mock point sources at the image center
        gdat.pntscntr = pntscntr
        ### mock image resolution
        gdat.numbsidecart = numbsidecart
        gdat.numbsideheal = numbsideheal

    ## proposal frequencies
    gdat.probprop = probprop


    # check inputs
    if gdat.colrprio:
        if gdat.minmspec.size != 1 or gdat.minmspec.size != 1:
            print 'Spectral limits must be numpy scalars when color priors are used!'
    
    # date and time
    gdat.strgtime = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
    if gdat.verbtype > 0:
        print 'PCAT started at ', gdat.strgtime
        print 'Initializing...'
    
    # setup the sampler
    setp(gdat) 

    if gdat.verbtype > 1:
        print 'probprop: '
        print vstack((arange(gdat.numbprop), gdat.strgprop, gdat.probprop)).T
        print 'indxsampnumbpnts: ', gdat.indxsampnumbpnts
        print 'indxsampfdfnnorm: ', gdat.indxsampfdfnnorm
        print 'indxsampfdfnslop: ', gdat.indxsampfdfnslop
        print 'indxsamppsfipara: ', gdat.indxsamppsfipara
        print 'indxsampnormback: '
        print gdat.indxsampnormback
        print 'indxsampcompinit: ', gdat.indxsampcompinit
        print 'maxmangleval'
        print gdat.maxmangleval
        if gdat.trueinfo:
            print 'truelgal: ', gdat.truelgal
            print 'truebgal: ', gdat.truebgal
            print 'truespec: '
            print gdat.truespec
            print 'truenumbpnts: ', gdat.truenumbpnts
            print 'truefdfnslop: ', gdat.truefdfnslop
            print 'truenormback: '
            print gdat.truenormback
            print
            if gdat.datatype == 'mock':
                print 'mocknumbpnts: ', gdat.mocknumbpnts
                print 'mockfdfnslop: ', gdat.mockfdfnslop
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

    for k in range(gdat.numbproc):
        timereal[k] = gridchan[k][18]
        timeproc[k] = gridchan[k][19]


    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        tim0 = time.time()

    # parse the sample bundle
    listsampvarb = zeros((gdat.numbsamp, gdat.numbproc, gdat.maxmsampsize))
    listindxprop = zeros((gdat.numbswep, gdat.numbproc))
    listchrollik = zeros((gdat.numbswep, gdat.numbproc, gdat.numbchrollik))
    listchrototl = zeros((gdat.numbswep, gdat.numbproc, gdat.numbchrototl))
    listllik = zeros((gdat.numbsamp, gdat.numbproc))
    listlpri = zeros((gdat.numbsamp, gdat.numbproc))
    listaccp = zeros((gdat.numbswep, gdat.numbproc))
    listmodlcnts = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbpixlsave))
    listpntsfluxmean = zeros((gdat.numbsamp, gdat.numbproc, gdat.numbener))
    listindxpntsfull = []
    listindxsampmodi = zeros((gdat.numbswep, gdat.numbproc), dtype=int)
    gdat.listauxipara = empty((gdat.numbswep, gdat.numbproc, gdat.numbcomp))
    gdat.listlaccfrac = empty((gdat.numbswep, gdat.numbproc))
    gdat.listnumbpair = empty((gdat.numbswep, gdat.numbproc))
    gdat.listjcbnfact = empty((gdat.numbswep, gdat.numbproc))
    gdat.listcombfact = empty((gdat.numbswep, gdat.numbproc))

    levi = 0.
    info = 0.
    
    for k in range(gdat.numbproc):
        gdat.rtag = retr_rtag(gdat, k)
        listchan = gridchan[k]
        listsampvarb[:, k, :] = listchan[0]
        listindxprop[:, k] = listchan[1]
        listchrototl[:, k, :] = listchan[2]
        listllik[:, k] = listchan[3]
        listlpri[:, k] = listchan[4]
        listaccp[:, k] = listchan[5]
        listmodlcnts[:, k, :] = listchan[6]    
        listindxpntsfull.append(listchan[7])
        listindxsampmodi[:, k] = listchan[8]
        gdat.listauxipara[:, k, :] = listchan[9]
        gdat.listlaccfrac[:, k] = listchan[10]
        gdat.listnumbpair[:, k] = listchan[11]
        gdat.listjcbnfact[:, k] = listchan[12]
        gdat.listcombfact[:, k] = listchan[13]
        levi += listchan[14]
        info += listchan[15]
        listpntsfluxmean[:, k, :] = listchan[16]
        listchrollik[:, k, :] = listchan[17]

    listindxprop = listindxprop.flatten()
    gdat.listauxipara = gdat.listauxipara.reshape((gdat.numbswep * gdat.numbproc, gdat.numbcomp))
    gdat.listlaccfrac = gdat.listlaccfrac.reshape(gdat.numbswep * gdat.numbproc)
    gdat.listnumbpair = gdat.listnumbpair.reshape(gdat.numbswep * gdat.numbproc)
    gdat.listjcbnfact = gdat.listjcbnfact.reshape(gdat.numbswep * gdat.numbproc)
    gdat.listcombfact = gdat.listcombfact.reshape(gdat.numbswep * gdat.numbproc)
    
    listchrototl = listchrototl.reshape((gdat.numbproc * gdat.numbswep, gdat.numbchrototl)) 
    listchrollik = listchrollik.reshape((gdat.numbproc * gdat.numbswep, gdat.numbchrollik)) 
    listaccp = listaccp.flatten()
    listindxsampmodi = listindxsampmodi.flatten()
    
    gdat.rtag = retr_rtag(gdat, None)

    # collect posterior samples from the processes
    ## number of PS
    listnumbpnts = listsampvarb[:, :, gdat.indxsampnumbpnts].astype(int).reshape(gdat.numbsamp * gdat.numbproc, -1)
    ## FDF normalization
    listfdfnnorm = listsampvarb[:, :, gdat.indxsampfdfnnorm].reshape(gdat.numbsamp * gdat.numbproc, -1)
    ## FDF slope
    listfdfnslop = listsampvarb[:, :, gdat.indxsampfdfnslop].reshape(gdat.numbsamp * gdat.numbproc, gdat.numbpopl, gdat.numbener)
    ## PSF parameters
    listpsfipara = listsampvarb[:, :, gdat.indxsamppsfipara].reshape(gdat.numbsamp * gdat.numbproc, -1)
    ## Background normalization
    listnormback = listsampvarb[:, :, gdat.indxsampnormback].reshape(gdat.numbsamp * gdat.numbproc, gdat.numbback, gdat.numbener)
    ## PS parameters
    listlgal = [[] for l in gdat.indxpopl]
    listbgal = [[] for l in gdat.indxpopl]
    listspec = [[] for l in gdat.indxpopl]
    if gdat.colrprio:
        listsind = [[] for l in gdat.indxpopl]
    listspechist = empty((gdat.numbsamp * gdat.numbproc, gdat.numbpopl, gdat.numbener, gdat.numbspec))
    for k in range(gdat.numbproc):
        for j in range(gdat.numbsamp):            
            n = k * gdat.numbsamp + j
            indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, listindxpntsfull[k][j])
            for l in gdat.indxpopl:
                listlgal[l].append(listsampvarb[j, k, indxsamplgal[l]])
                listbgal[l].append(listsampvarb[j, k, indxsampbgal[l]])
                listspec[l].append(listsampvarb[j, k, indxsampspec[l]])
                if gdat.colrprio:
                    listsind[l].append(listsampvarb[j, k, indxsampsind[l]])
                for i in gdat.indxener:
                    listspechist[n, l, i, :] = histogram(listspec[l][n][i, :], gdat.binsspec[i, :])[0]
    for l in gdat.indxpopl:
        listlgaltemp = zeros((gdat.numbsamp, gdat.maxmnumbpnts[l])) - 1.
        listbgaltemp = zeros((gdat.numbsamp, gdat.maxmnumbpnts[l])) - 1.
        listspectemp = zeros((gdat.numbsamp, gdat.numbener, gdat.maxmnumbpnts[l])) - 1.
        listsindtemp = zeros((gdat.numbsamp, gdat.maxmnumbpnts[l])) - 1.
        for k in range(gdat.numbsamp):
            listlgaltemp[k, 0:listlgal[l][k].size] = listlgal[l][k]
            listbgaltemp[k, 0:listbgal[l][k].size] = listbgal[l][k]
            listspectemp[k, :, 0:listspec[l][k].shape[1]] = listspec[l][k]
            listsindtemp[k, 0:listsind[l][k].size] = listsind[l][k]
        listlgal[l] = listlgaltemp
        listbgal[l] = listbgaltemp 
        listspec[l] = listspectemp    
        listsind[l] = listsindtemp

    # auxiliary variables
    listpntsfluxmean = listpntsfluxmean.reshape(gdat.numbsamp * gdat.numbproc, gdat.numbener)
   
        
        
    if gdat.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        tim0 = time.time()

    # posterior maps
    pntsprob = zeros((gdat.numbpopl, gdat.numbener, gdat.numbpixl, gdat.numbspec))
    for k in range(gdat.numbsamp):
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for h in range(gdat.numbspec):
                    indxpnts = where((gdat.binsspec[i, h] < listspec[l][k, i, :]) & (listspec[l][k, i, :] < gdat.binsspec[i, h+1]))[0]
                    hpixl = retr_indxpixl(gdat, listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    
        
    if gdat.verbtype > 0:
        print 'Performing Gelman-Rubin convergence test...'
        tim0 = time.time()

    gmrbstat = zeros(gdat.numbpixlsave)
    for n in range(gdat.numbpixlsave):
        gmrbstat[n] = tdpy.mcmc.gmrb_test(listmodlcnts[:, :, n])


            
    pathprobcatl = os.environ["PCAT_DATA_PATH"] + '/probcatl_' + gdat.strgtime + '_' + gdat.rtag + '.fits'  
    
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
    head['numbspec'] = gdat.numbspec
    if gdat.pixltype == 'heal':
        head['numbsideheal'] = gdat.numbsideheal
    if gdat.pixltype == 'cart':
        head['numbsidecart'] = gdat.numbsidecart
    
    head['minmlgal'] = gdat.minmlgal
    head['maxmlgal'] = gdat.maxmlgal
    head['minmbgal'] = gdat.minmbgal
    head['maxmbgal'] = gdat.maxmbgal
    
    head['datatype'] = gdat.datatype
    head['regitype'] = gdat.regitype
    head['psfntype'] = gdat.psfntype
    head['exprtype'] = gdat.exprtype
    head['pixltype'] = gdat.pixltype
    head['fdfntype'] = gdat.fdfntype
    
    head['colrprio'] = gdat.colrprio
    head['trueinfo'] = gdat.trueinfo
    head['margsize'] = gdat.margsize
    head['strgtime'] = gdat.strgtime
    
    head['stdvfdfnnorm'] = gdat.stdvfdfnnorm
    head['stdvfdfnslop'] = gdat.stdvfdfnslop
    head['stdvpsfipara'] = gdat.stdvpsfipara
    head['stdvback'] = gdat.stdvback
    head['stdvlbhl'] = gdat.stdvlbhl
    head['stdvspec'] = gdat.stdvspec
    head['spmrlbhl'] = gdat.spmrlbhl
    head['fracrand'] = gdat.fracrand

    if gdat.trueinfo and gdat.datatype == 'mock':
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

    for l in gdat.indxpopl:
        listhdun.append(pf.ImageHDU(listlgal[l]))
        listhdun[-1].header['EXTNAME'] = 'lgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listbgal[l]))
        listhdun[-1].header['EXTNAME'] = 'bgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listspec[l]))
        listhdun[-1].header['EXTNAME'] = 'specpopl%d' % l
        if gdat.colrprio:
            listhdun.append(pf.ImageHDU(listsind[l]))
            listhdun[-1].header['EXTNAME'] = 'sindpopl%d' % l

    
    listhdun.append(pf.ImageHDU(gdat.minmnormback))
    listhdun[-1].header['EXTNAME'] = 'minmnormback'

    listhdun.append(pf.ImageHDU(gdat.maxmnormback))
    listhdun[-1].header['EXTNAME'] = 'maxmnormback'


    listhdun.append(pf.ImageHDU(gdat.maxmnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'maxmnumbpnts'

    listhdun.append(pf.ImageHDU(listnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'numbpnts'

    listhdun.append(pf.ImageHDU(listfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'fdfnnorm'
    
    listhdun.append(pf.ImageHDU(listfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'fdfnslop'
    
    listhdun.append(pf.ImageHDU(listpsfipara))
    listhdun[-1].header['EXTNAME'] = 'psfipara'
    
    listhdun.append(pf.ImageHDU(listnormback))
    listhdun[-1].header['EXTNAME'] = 'normback'
    
    listhdun.append(pf.ImageHDU(listspechist))
    listhdun[-1].header['EXTNAME'] = 'spechist'
    
    
    listhdun.append(pf.ImageHDU(listllik))
    listhdun[-1].header['EXTNAME'] = 'llik'
    listhdun.append(pf.ImageHDU(listlpri))
    listhdun[-1].header['EXTNAME'] = 'lpri'
    
    if gdat.trueinfo and gdat.datatype == 'mock':
        listhdun.append(pf.ImageHDU(gdat.mocknumbpnts))
        listhdun[-1].header['EXTNAME'] = 'mocknumbpnts'
        
        listhdun.append(pf.ImageHDU(gdat.mockfdfnslop))
        listhdun[-1].header['EXTNAME'] = 'mockfdfnslop'
        
        listhdun.append(pf.ImageHDU(gdat.mocknormback))
        listhdun[-1].header['EXTNAME'] = 'mocknormback'

    # convergence diagnostics
    listhdun.append(pf.ImageHDU(gmrbstat))
    listhdun[-1].header['EXTNAME'] = 'gmrbstat'
    listhdun.append(pf.ImageHDU(listmodlcnts))
    listhdun[-1].header['EXTNAME'] = 'modlcnts'
    
    
    listhdun.append(pf.ImageHDU(pntsprob))
    listhdun[-1].header['EXTNAME'] = 'pntsprob'
   
    # truth information
    if gdat.trueinfo:
        listhdun.append(pf.ImageHDU(gdat.truenumbpnts))
        listhdun[-1].header['EXTNAME'] = 'truenumbpnts'

        for l in gdat.indxpopl:
            listhdun.append(pf.ImageHDU(gdat.truelgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truelgalpopl%d' % l

            listhdun.append(pf.ImageHDU(gdat.truebgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truebgalpopl%d' % l

            listhdun.append(pf.ImageHDU(gdat.truespec[l]))
            listhdun[-1].header['EXTNAME'] = 'truespecpopl%d' % l

            if gdat.colrprio:
                listhdun.append(pf.ImageHDU(gdat.truesind[l]))
                listhdun[-1].header['EXTNAME'] = 'truesindpopl%d' % l

        listhdun.append(pf.ImageHDU(gdat.truefdfnslop))
        listhdun[-1].header['EXTNAME'] = 'truefdfnslop'

        listhdun.append(pf.ImageHDU(gdat.truenormback))
        listhdun[-1].header['EXTNAME'] = 'truenormback'

        listhdun.append(pf.ImageHDU(gdat.truepsfipara))
        listhdun[-1].header['EXTNAME'] = 'truepsfipara'

        listhdun.append(pf.ImageHDU(gdat.indxtruepntstimevari))
        listhdun[-1].header['EXTNAME'] = 'indxtruepntstimevari'

    if gdat.colrprio:
        gdat.minmspec = gdat.minmspec[gdat.indxenerfdfn]
        gdat.maxmspec = gdat.maxmspec[gdat.indxenerfdfn]

    # hyperprior limits
    listhdun.append(pf.ImageHDU(gdat.minmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'minmfdfnnorm'
    listhdun.append(pf.ImageHDU(gdat.maxmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'maxmfdfnnorm'
    listhdun.append(pf.ImageHDU(gdat.minmfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'minmfdfnslop'
    listhdun.append(pf.ImageHDU(gdat.maxmfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'maxmfdfnslop'

    listhdun.append(pf.ImageHDU(gdat.minmspec))
    listhdun[-1].header['EXTNAME'] = 'minmspec'
    listhdun.append(pf.ImageHDU(gdat.maxmspec))
    listhdun[-1].header['EXTNAME'] = 'maxmspec'
    if gdat.colrprio:
        listhdun.append(pf.ImageHDU(gdat.minmsind))
        listhdun[-1].header['EXTNAME'] = 'minmsind'
        listhdun.append(pf.ImageHDU(gdat.maxmsind))
        listhdun[-1].header['EXTNAME'] = 'maxmsind'
        listhdun.append(pf.ImageHDU(gdat.meansind))
        listhdun[-1].header['EXTNAME'] = 'meansind'
        listhdun.append(pf.ImageHDU(gdat.stdvsind))
        listhdun[-1].header['EXTNAME'] = 'stdvsind'
    
    listhdun.append(pf.ImageHDU(gdat.binsener))
    listhdun[-1].header['EXTNAME'] = 'binsener'
    listhdun.append(pf.ImageHDU(gdat.indxenerfdfn))
    listhdun[-1].header['EXTNAME'] = 'indxenerfdfn'
    
    listhdun.append(pf.ImageHDU(gdat.indxenerincl))
    listhdun[-1].header['EXTNAME'] = 'indxenerincl'
    
    listhdun.append(pf.ImageHDU(gdat.indxevttincl))
    listhdun[-1].header['EXTNAME'] = 'indxevttincl'
 
    # utilities
    listhdun.append(pf.ImageHDU(listpntsfluxmean))
    listhdun[-1].header['EXTNAME'] = 'listpntsfluxmean'
    
    listhdun.append(pf.ImageHDU(gdat.probprop))
    listhdun[-1].header['EXTNAME'] = 'probprop'

    listhdun.append(pf.ImageHDU(listindxprop))
    listhdun[-1].header['EXTNAME'] = 'indxprop'
    
    listhdun.append(pf.ImageHDU(listchrototl))
    listhdun[-1].header['EXTNAME'] = 'listchrototl'
    
    listhdun.append(pf.ImageHDU(listchrollik))
    listhdun[-1].header['EXTNAME'] = 'listchrollik'
    
    listhdun.append(pf.ImageHDU(listaccp))
    listhdun[-1].header['EXTNAME'] = 'accp'
    
    listhdun.append(pf.ImageHDU(listindxsampmodi))
    listhdun[-1].header['EXTNAME'] = 'sampmodi'
    
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
    
    pf.HDUList(listhdun).writeto(pathprobcatl, clobber=True, output_verify='ignore')

    if gdat.makeplot:
        plot_post(pathprobcatl)

    timetotlreal = time.time() - timetotlreal
    timetotlproc = time.clock() - timetotlproc
     
    if gdat.verbtype > 0:
        for k in range(gdat.numbproc):
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, timereal[k], timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (timetotlreal, sum(timeproc))
        print 'The ensemble of catalogs is at ' + pathprobcatl
        if gdat.makeplot:
            print 'The plots are in ' + gdat.plotpath
        
    return gridchan
    
    
def plot_samp(gdat):

    gdat.thisresicnts = gdat.datacnts - gdat.thismodlcnts

    gdat.thispsfn = retr_psfn(gdat, gdat.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.angldisp, gdat.psfntype)
    if gdat.pixltype == 'cart':
        gdat.thispsfn = gdat.thispsfn.reshape((gdat.numbener, -1, gdat.numbevtt))
    thisfwhm = 2. * retr_psfnwdth(gdat, gdat.thispsfn, 0.5)

    plot_psfn(gdat)

    gdat.thisbackcntsmean = empty((gdat.numbener, gdat.numbevtt))
    for c in gdat.indxback:
        gdat.thisbackcntsmean += mean(gdat.thissampvarb[gdat.indxsampnormback[c, :, None, None]] * gdat.backflux[c] * gdat.expo * \
            gdat.diffener[:, None, None] * pi * thisfwhm[:, None, :]**2 / 4., 1)

    gdat.temppntsflux, gdat.temppntscnts, gdat.tempmodlflux, gdat.tempmodlcnts = retr_maps(gdat, gdat.thisindxpntsfull, gdat.thissampvarb)
    
    gdat.thispntscnts = gdat.thispntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    gdat.errrpnts = 100. * (gdat.thispntscnts - gdat.temppntscnts) / gdat.temppntscnts

    thiscnts = []
    for l in gdat.indxpopl:
    
        indxpixltemp = retr_indxpixl(gdat, gdat.thissampvarb[gdat.thisindxsampbgal[l]], gdat.thissampvarb[gdat.thisindxsamplgal[l]])
        cntstemp = gdat.thissampvarb[gdat.thisindxsampspec[l]][:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]
        thiscnts.append(cntstemp)
        if gdat.colrprio:
            plot_histsind(gdat, l)
        plot_scatpixl(gdat, l)
        if gdat.trueinfo:
            indxmodl, gdat.trueindxpntsbias, gdat.trueindxpntsmiss = pair_catl(gdat, l, gdat.thissampvarb[gdat.thisindxsamplgal[l]], \
                gdat.thissampvarb[gdat.thisindxsampbgal[l]], gdat.thissampvarb[gdat.thisindxsampspec[l]])
            thisspecmtch = gdat.thissampvarb[gdat.thisindxsampspec[l]][:, indxmodl]
            thisspecmtch[:, gdat.trueindxpntsmiss] = 0.
            plot_scatspec(gdat, l, thisspecmtch=thisspecmtch)
        plot_histspec(gdat, l)
        plot_histcnts(gdat, l, thiscnts)
        if gdat.numbback == 2:
            plot_compfrac(gdat)

    for i in gdat.indxener:
        
        plot_datacnts(gdat, i, None)
        plot_modlcnts(gdat, i, None)
        plot_resicnts(gdat, i, None)
        
        plot_thispntscnts(gdat, i, None)
        plot_temppntscnts(gdat, i, None)
        plot_errrpnts(gdat, i, None)
        
        #plot_catlfram(gdat, i, None, thiscnts)

        #for m in gdat.indxevtt:
        #    plot_datacnts(gdat, i, m)
        #    plot_catl(gdat, i, m, thiscnts)
        #    plot_modlcnts(gdat, i, m)
        #    plot_resicnts(gdat, i, m, gdat.thisresicnts)
        

    #if gdat.numbener == 3:
    #    plot_datacnts(gdat, None, None)
        
    #plot_fwhm(gdat, thisfwhm)
    #plot_backcntsmean(gdat, gdat.thisbackcntsmean)
    
    if amax(abs(gdat.errrpnts)) > 0.1 and False:
        print 'Approximation error went above the limit!'
    
    
def rjmc(gdat, indxprocwork):

    # sweeps to be saved
    boolsave = zeros(gdat.numbswep, dtype=bool)
    indxswepsave = arange(gdat.numbburn, gdat.numbswep, gdat.factthin)
    boolsave[indxswepsave] = True
    
    sampindx = zeros(gdat.numbswep, dtype=int)
    sampindx[indxswepsave] = arange(gdat.numbsamp)

    listsampvarb = zeros((gdat.numbsamp, gdat.maxmsampsize)) + -1.
    listindxprop = zeros(gdat.numbswep)
    listchrototl = zeros((gdat.numbswep, gdat.numbchrototl))
    listchrollik = zeros((gdat.numbswep, gdat.numbchrollik))
    listllik = zeros(gdat.numbsamp)
    listlprising = zeros(gdat.numbsamp)
    listlpri = zeros(gdat.numbsamp)
    listaccp = zeros(gdat.numbswep, dtype=bool)
    listaccpspec = []
    listindxsampmodi = zeros(gdat.numbswep, dtype=int)
    listmodlcnts = zeros((gdat.numbsamp, gdat.numbpixlsave))
    listpntsfluxmean = zeros((gdat.numbsamp, gdat.numbener))
    listindxpntsfull = []
    
    gdat.listauxipara = zeros((gdat.numbswep, gdat.numbcomp))
    gdat.listlaccfrac = zeros(gdat.numbswep)
    gdat.listnumbpair = zeros(gdat.numbswep)
    gdat.listjcbnfact = zeros(gdat.numbswep)
    gdat.listcombfact = zeros(gdat.numbswep)

    gdat.cntrswep = 0
    
    # initialize the chain
    retr_llik(gdat, listchrollik, init=True)
    retr_lpri(gdat, init=True)

    # current sample index
    thiscntr = -1
    
    while gdat.cntrswep < gdat.numbswep:
        
        timeinit = time.time()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdat.cntrswep

        thismakefram = (gdat.cntrswep % gdat.plotperd == 0) and indxprocwork == int(float(gdat.cntrswep) / gdat.numbswep * gdat.numbproc) and gdat.makeplot
        gdat.reje = False
    
        # choose a proposal type
        retr_thisindxprop(gdat, gdat.drmcsamp[:, 0])
            
        # save the proposal type
        listindxprop[gdat.cntrswep] = gdat.thisindxprop
        if gdat.verbtype > 1:
            print 'thisindxprop: ', gdat.strgprop[gdat.thisindxprop]
        
        if gdat.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop(gdat)
        timefinl = time.time()
        listchrototl[gdat.cntrswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                gdat.lock.acquire()
            print 'Process %d started making a frame.' % indxprocwork
            plot_samp(gdat)
            print 'Process %d finished making a frame.' % indxprocwork
            if gdat.numbproc > 1:
                gdat.lock.release()
            
        # reject the sample if proposal is outside the prior
        if gdat.thisindxprop == gdat.indxpropfdfnslop:
            if gdat.drmcsamp[gdat.indxsampfdfnslopmodi, 1] < 0. or gdat.drmcsamp[gdat.indxsampfdfnslopmodi, 1] > 1.:
                gdat.reje = True
        elif gdat.thisindxprop != gdat.indxpropbrth and gdat.thisindxprop != gdat.indxpropdeth and not gdat.reje:
            if where((gdat.drmcsamp[gdat.indxsampmodi, 1] < 0.) | (gdat.drmcsamp[gdat.indxsampmodi, 1] > 1.))[0].size > 0:
                gdat.reje = True

        if gdat.thisindxprop == gdat.indxproppsfipara:
            if gdat.psfntype == 'doubking':
                if gdat.nextsampvarb[gdat.indxsamppsfipara[1]] >= gdat.nextsampvarb[gdat.indxsamppsfipara[3]]:
                    gdat.reje = True
            elif gdat.psfntype == 'doubgaus':
                if gdat.nextsampvarb[gdat.indxsamppsfipara[1]] >= gdat.nextsampvarb[gdat.indxsamppsfipara[2]]:
                    gdat.reje = True
            
        if not gdat.reje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri(gdat)
            timefinl = time.time()
            listchrototl[gdat.cntrswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(gdat, listchrollik) 
            timefinl = time.time()
            listchrototl[gdat.cntrswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(gdat.deltllik + gdat.deltlpri + gdat.laccfrac)

            if gdat.verbtype > 1:
                print 'deltlpri'
                print gdat.deltlpri
                print 'deltllik'
                print gdat.deltllik
                print 'laccfrac'
                print gdat.laccfrac
                print
            
            # make diagnostic plots for the next step
            # temp
            if False and thismakefram and gdat.thisindxprop >= gdat.indxproppsfipara:
                gdat.thispntsfluxmodi = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt)) - 1e-6
                gdat.nextpntsfluxmodi = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt)) - 1e-6
                gdat.diffpntsfluxmodi = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt)) - 1e-6
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        gdat.thispntsfluxmodi[i, gdat.indxpixlmodi, m] = copy(gdat.thispntsflux[i, gdat.indxpixlmodi, m])
                        gdat.nextpntsfluxmodi[i, gdat.indxpixlmodi, m] = copy(gdat.nextpntsflux[i, gdat.indxpixlmodi, m])
                        gdat.diffpntsfluxmodi[i, gdat.indxpixlmodi, m] = gdat.nextpntsfluxmodi[i, gdat.indxpixlmodi,  m] - gdat.thispntsfluxmodi[i, gdat.indxpixlmodi, m]
                        plot_nextstat(gdat, i, m)
                        plot_diag(gdat, i, m)
        else:
            accpprob = 0.
    
        # accept the sample
        if accpprob >= rand():

            if gdat.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(gdat)

            listaccp[gdat.cntrswep] = True

        # reject the sample
        else:

            if gdat.verbtype > 1:
                print 'Rejected.'

            listaccp[gdat.cntrswep] = False
             
        # sanity checks
        indxsampbadd = where((gdat.drmcsamp[1:, 0] > 1.) | (gdat.drmcsamp[1:, 0] < 0.))[0] + 1
        if indxsampbadd.size > 0:
            print 'Unit sample vector went outside [0,1]!'
            print 'cntrswep'
            print gdat.cntrswep
            print 'thisindxprop'
            print gdat.thisindxprop
            print 'indxsampbadd'
            print indxsampbadd
            print 'drmcsamp'
            print gdat.drmcsamp[indxsampbadd, :]
            return
            
        for l in gdat.indxpopl:
            for i in gdat.indxenerfdfn:
                if where(gdat.thissampvarb[gdat.thisindxsampspec[l][i, :]] < gdat.minmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went below the prior range!'
                if where(gdat.thissampvarb[gdat.thisindxsampspec[l][i, :]] > gdat.maxmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went above the prior range!'          

        # save the sample
        if boolsave[gdat.cntrswep]:
            listsampvarb[sampindx[gdat.cntrswep], :] = gdat.thissampvarb
            listmodlcnts[sampindx[gdat.cntrswep], :] = gdat.thismodlcnts[0, gdat.indxpixlsave, 0]
            listpntsfluxmean[sampindx[gdat.cntrswep], :] = mean(sum(gdat.thispntsflux * gdat.expo, 2) / sum(gdat.expo, 2), 1)
            listindxpntsfull.append(gdat.thisindxpntsfull)
            listllik[sampindx[gdat.cntrswep]] = sum(gdat.thisllik)
            
            lpri = 0.
            for l in gdat.indxpopl:
                numbpnts = gdat.thissampvarb[gdat.indxsampnumbpnts[l]]
                fdfnnorm = gdat.thissampvarb[gdat.indxsampfdfnnorm[l]]
                lpri += numbpnts * gdat.priofactlgalbgal + gdat.priofactfdfnslop + gdat.priofactfdfnnorm - log(fdfnnorm)
                for i in gdat.indxenerfdfn:
                    flux = gdat.thissampvarb[gdat.thisindxsampspec[l][i, :]]
                    fdfnslop = gdat.thissampvarb[gdat.indxsampfdfnslop[l, i]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_spec(gdat, flux, fdfnslop, gdat.minmspec[i], gdat.maxmspec[i])))
            listlpri[sampindx[gdat.cntrswep]] = lpri
            
            if gdat.tracsamp:
                
                numbpnts = gdat.thissampvarb[gdat.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[gdat.cntrswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            tim1 = time.time()
            listchrototl[gdat.cntrswep, 0] = tim1 - timeinit

        # log the progress
        if gdat.verbtype > 0:
            thiscntr = tdpy.util.show_prog(gdat.cntrswep, gdat.numbswep, thiscntr, indxprocwork=indxprocwork)
            
        if gdat.diagsamp:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_diagfram(i, m)
                    plot_nextfram(i, m)
        
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
    
    listchrollik = array(listchrollik)

    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(gdat.datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    minmlistllik = amin(listllik)
    levi = -log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    info = mean(listllik) - levi

    listchan = [listsampvarb, listindxprop, listchrototl, listllik, listlpri, listaccp, listmodlcnts, listindxpntsfull, listindxsampmodi, \
        gdat.listauxipara, gdat.listlaccfrac, gdat.listnumbpair, gdat.listjcbnfact, gdat.listcombfact, \
        levi, info, listpntsfluxmean, listchrollik]
    
    return listchan

