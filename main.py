# common imports
from __init__ import *

# internal functions
from util import *
from visu import *

def work(gdat, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for this chain
    seed()
    
    # construct the run tag for this chain
    gdat.rtag = retr_rtag(gdat, indxprocwork)
    
    # empty object to hold chain-specific variables that will be modified by the chain
    gdatmodi = gdatstrt()
    
    # construct the initial state
    ## unit sample vector
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
    if gdat.initmeanpnts != None or (not gdat.randinit and gdat.datatype == 'mock'):
        if gdat.initmeanpnts != None:
            meanpnts = gdat.initmeanpnts
        else:
            meanpnts = gdat.truemeanpnts
        for l in gdat.indxpopl:
            gdatmodi.drmcsamp[gdat.indxsampmeanpnts[l], 0] = cdfn_logt(meanpnts[l], gdat.minmmeanpnts[l], gdat.factmeanpnts[l])
    else:
        gdatmodi.drmcsamp[gdat.indxsampmeanpnts, 0] = rand(gdat.numbpopl)
   
    ## hyperparameters
    ### boolean flag to indicate that the initial state will not be random
    detrinittemp = not gdat.randinit and gdat.datatype == 'mock'
    
    ### flux distribution
    for l in gdat.indxpopl:
        detrinit = detrinittemp and gdat.fluxdisttype[l] == gdat.mockfluxdisttype[l]
    
        #### single power law
        if gdat.fluxdisttype[l] == 'powr':
            if gdat.initfluxdistslop[l] != None or detrinit:
                if gdat.initfluxdistslop[l] != None:
                    fluxdistslop = gdat.initfluxdistslop[l]
                else:
                    fluxdistslop = gdat.truefluxdistslop[l]
                gdatmodi.drmcsamp[gdat.indxsampfluxdistslop[l], 0] = cdfn_atan(fluxdistslop, gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
            else:
                gdatmodi.drmcsamp[gdat.indxsampfluxdistslop[l], 0] = rand()
       
        #### broken power law
        if gdat.fluxdisttype[l] == 'brok':
            if gdat.initfluxdistbrek[l] != None or detrinit:
                if gdat.initfluxdistbrek[l] != None:
                   fluxdistbrek = gdat.initfluxdistbrek[l]
                else:
                   fluxdistbrek = gdat.truefluxdistbrek[l]
                gdatmodi.drmcsamp[gdat.indxsampfluxdistbrek[l], 0] = cdfn_logt(fluxdistbrek, gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
            else:
                gdatmodi.drmcsamp[gdat.indxsampfluxdistbrek, 0] = rand()
                
            if gdat.initfluxdistsloplowr[l] != None or detrinit:
                if gdat.initfluxdistsloplowr[l] != None:
                    fluxdistsloplowr = gdat.initfluxdistsloplowr[l]
                else:
                    fluxdistsloplowr = gdat.truefluxdistsloplowr[l]
                gdatmodi.drmcsamp[gdat.indxsampfluxdistsloplowr[l], 0] = cdfn_atan(fluxdistsloplowr, gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
            else:
                gdatmodi.drmcsamp[gdat.indxsampfluxdistsloplowr, 0] = rand()

            if gdat.initfluxdistslopuppr[l] != None or detrinit:
                if gdat.initfluxdistslopuppr[l] != None:
                    fluxdistslopuppr = gdat.initfluxdistslopuppr[l]
                else:
                    fluxdistslopuppr = gdat.truefluxdistslopuppr[l]
                gdatmodi.drmcsamp[gdat.indxsampfluxdistslopuppr[l], 0] = cdfn_atan(fluxdistslopuppr, gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])
            else:
                gdatmodi.drmcsamp[gdat.indxsampfluxdistslopuppr, 0] = rand()

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
            if gdat.fluxdisttype[l] == 'powr':
                fluxdistslop = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfluxdistslop[l], 0], gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
                fluxunit = cdfn_flux_powr(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], fluxdistslop)
            if gdat.fluxdisttype[l] == 'brok':
                flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
                fluxdistbrek = icdf_logt(gdatmodi.drmcsamp[gdat.indxsampfluxdistbrek[l], 0], gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
                fluxdistsloplowr = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfluxdistsloplowr[l], 0], gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
                fluxdistslopuppr = icdf_atan(gdatmodi.drmcsamp[gdat.indxsampfluxdistslopuppr[l], 0], gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])
                fluxunit = cdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
            gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :], 0] = copy(fluxunit)
            gdatmodi.drmcsamp[gdatmodi.thisindxsampsind[l], 0] = cdfn_eerr(gdat.truesind[l], gdat.sinddistmean[l], gdat.sinddiststdv[l], \
                                                                                                        gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])

    ## sample vector
    gdatmodi.thissampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.drmcsamp[:, 0])
    
    ## PS and total flux and count maps
    gdatmodi.thispntsflux, gdatmodi.thispntscnts, gdatmodi.thismodlflux, gdatmodi.thismodlcnts = retr_maps(gdat, gdatmodi.thisindxpntsfull, gdatmodi.thissampvarb)
    
    ## indices of the PS parameters
    indxsamplgaltemp, indxsampbgaltemp, indxsampspectemp, indxsampsindtemp, indxsampcomptemp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    ## PSF
    if gdat.boolintpanglcosi:
        gdatmodi.thispsfnintp = interp1d(gdat.binsanglcosi, retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype), axis=1)
    else:
        gdatmodi.thispsfnintp = interp1d(gdat.binsangl, retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype), axis=1)
    gdatmodi.thispsfn = gdatmodi.thispsfnintp(gdat.binsangl)
    
    # log-prior
    if gdat.bindprio:
        gdatmodi.thislpri = empty((gdat.numbpopl, gdat.numbflux))
    else:
        gdatmodi.thislpri = empty((gdat.numbpopl, 2))

    # check the initial unit sample vector for bad entries
    indxsampbaddlowr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.)[0] + gdat.numbpopl
    indxsampbadduppr = where(gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.)[0] + gdat.numbpopl
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'Initial unit sample vector went outside [0, 1]. Correcting it...'
        print 'bad index vector'
        print indxsampbadd
        print (indxsampbadd - gdat.indxsampcompinit ) % gdat.numbcomp
        print 
        gdatmodi.drmcsamp[indxsampbaddlowr, 0] = 0.
        gdatmodi.drmcsamp[indxsampbadduppr, 0] = 1.

    # allocate memory for variables to hold the proposed state
    ## sample vector
    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
    
    ## flux and count maps
    gdatmodi.nextpntsflux = zeros_like(gdatmodi.thispntsflux)
    gdatmodi.nextmodlflux = zeros_like(gdatmodi.thispntsflux)
    gdatmodi.nextmodlcnts = zeros_like(gdatmodi.thispntsflux)
    
    ## likelihood
    gdatmodi.nextllik = zeros_like(gdatmodi.thispntsflux)
    gdatmodi.deltllik = 0.
    
    ## modification catalog
    gdatmodi.modilgal = empty(gdat.numbmodipnts)
    gdatmodi.modibgal = empty(gdat.numbmodipnts)
    gdatmodi.modisind = empty(gdat.numbmodipnts)
    gdatmodi.modispec = empty((gdat.numbener, gdat.numbmodipnts))
    
    # plotting variables
    gdatmodi.thisbackfwhmcnts = empty((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
   
    # log-prior
    gdatmodi.nextlpri = empty((gdat.numbpopl, gdat.numbflux))

    # log the initial state
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

    if gdat.verbtype > 1:
        print 'drmcsamp'
        for k in gdat.indxpara:
            print gdatmodi.drmcsamp[k, :]
        print 'thispsfipara'
        print gdatmodi.thissampvarb[gdat.indxsamppsfipara]
        print 'thissampvarb'
        for k in gdat.indxpara:
            print gdatmodi.thissampvarb[k]
        
    if gdat.verbtype > 1:
        show_memo(gdatmodi, 'gdatmodi')

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
         boolproppsfn=True, \
         boolpropfluxdist=True, \
         indxevttincl=arange(2, 4), \
         indxenerincl=arange(5), \
         strgexpr=None, \
         strgback=None, \
         lablback=None, \
         nameback=None, \
         strgexpo=None, \
         bindprio=False, \
         probprop=None, \
         pathdata='.', \
    
         strgfluxunit=None, \
         strgenerunit=None, \
    
         binsenerfull=None, \
         evttfull=None, \

         spatdisttype=None, \
         fluxdisttype=None, \
         sinddisttype=None, \
         spectype=None, \
         psfntype=None, \
         
         maxmnormback=None, \
         minmnormback=None, \
         maxmgang=None, \
         maxmgangmarg=None, \
         sinddistmean=None, \
         sinddiststdv=None, \
         minmmeanpnts=None, \
         maxmmeanpnts=None, \
         minmfluxdistslop=None, \
         maxmfluxdistslop=None, \
         minmfluxdistbrek=None, \
         maxmfluxdistbrek=None, \
         minmfluxdistsloplowr=None, \
         maxmfluxdistsloplowr=None, \
         minmfluxdistslopuppr=None, \
         maxmfluxdistslopuppr=None, \
         maxmnumbpnts=array([1000]), \
         minmflux=None, \
         maxmflux=None, \
         minmsind=None, \
         maxmsind=None, \
         
         # temp -- if datatype == 'inpt' trueinfo should depend on whether truexxxx are provided
         pntscntr=False, \
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         margfact=None, \
         maxmangl=None, \
         maxmangleval=None, \
         specfraceval=0.01, \
         anglassc=None, \
         
         stdvmeanpnts=0.05, \
         stdvfluxdistslop=0.1, \
         stdvpsfipara=0.1, \
         stdvback=0.04, \
         stdvlbhl=0.1, \
         stdvflux=0.15, \
         stdvsind=0.15, \
         fracrand=0.05, \
         radispmrlbhl=None, \
         stdvspmrsind=0.2, \
         
         initnumbpnts=None, \
         initmeanpnts=None, \
         initfluxdistslop=None, \
         initfluxdistbrek=None, \
         initfluxdistsloplowr=None, \
         initfluxdistslopuppr=None, \
         
         mockpixltype='heal', \
         mockspatdisttype=None, \
         mockfluxdisttype=None, \
         mocksinddisttype=None, \
         mockspectype=None, \
         mockpsfntype=None, \
         mockfluxdistslop=None, \
         mockfluxdistbrek=None, \
         mockfluxdistsloplowr=None, \
         mockfluxdistslopuppr=None, \
         mocksinddistmean=None, \
         mocksinddiststdv=None, \
         mocknormback=None, \
         mocknumbpnts=None, \
         mocknumbsidecart=None, \
         mocknumbsideheal=None, \
         
        ):
    
    # start the timer
    timetotlreal = time.time()
    timetotlproc = time.clock()
   
    # defaults
    
    ## convenience variables in order to set the defaults
    ### number of backgrounds
    numbback = len(strgback)
    ### number of energy bins
    numbener = indxenerincl.size
    ### number of populations
    numbpopl = maxmnumbpnts.size
    
    # PS parameter distribution models
    if spatdisttype == None:
        spatdisttype = array(['unif' for l in range(numbpopl)])
    if fluxdisttype == None:
        fluxdisttype = array(['powr' for l in range(numbpopl)])
    if sinddisttype == None:
        sinddisttype = array(['gaus' for l in range(numbpopl)])
    
    # PS spectral model
    if spectype == None:
        spectype = array(['powr' for l in range(numbpopl)])
    
    # PS parameter distributions
    if datatype == 'mock':
        mocknumbpopl = mocknumbpnts.size
        mocknormback = ones((numbback, numbener))
        if mockspatdisttype == None:
            mockspatdisttype = array(['unif' for l in range(mocknumbpopl)])
        if mockfluxdisttype == None:
            mockfluxdisttype = array(['powr' for l in range(mocknumbpopl)])
        if mocksinddisttype == None:
            mocksinddisttype = array(['gaus' for l in range(mocknumbpopl)])
        if mockspectype == None:
            mockspectype = array(['powr' for l in range(mocknumbpopl)])
        if mocksinddistmean == None:
            mocksinddistmean = zeros(mocknumbpopl) + 2.2
        if mocksinddiststdv == None:
            mocksinddiststdv = zeros(mocknumbpopl) + 0.3
        if mockfluxdistslop == None:
            mockfluxdistslop = array([2.])
        if mockfluxdistbrek == None:
            mockfluxdistbrek = array([2.])
        if mockfluxdistsloplowr == None:
            mockfluxdistsloplowr = array([1.5])
        if mockfluxdistslopuppr == None:
            mockfluxdistslopuppr = array([2.5])

    ## common
    if minmnormback == None:
        minmnormback = ones(numbback) * 0.5
    if maxmnormback == None:
        maxmnormback = ones(numbback) * 2.
        
    if minmsind == None:
        minmsind = 1.
    if maxmsind == None:
        maxmsind = 3.5
    
    if initfluxdistslop == None:
        initfluxdistslop = array([None for n in range(numbpopl)])
    if initfluxdistbrek == None:
        initfluxdistbrek = array([None for n in range(numbpopl)])
    if initfluxdistsloplowr == None:
        initfluxdistsloplowr = array([None for n in range(numbpopl)])
    if initfluxdistslopuppr == None:
        initfluxdistslopuppr = array([None for n in range(numbpopl)])
        
    if sinddistmean == None:
        sinddistmean = array([2.2])
    if sinddiststdv == None:
        sinddiststdv = array([0.3])

    if margfact == None:
        margfact = 1.1

    ## Fermi-LAT
    if exprtype == 'ferm':
        if maxmgang == None:
            maxmgang = deg2rad(20.)
        if psfntype == None:
            psfntype = 'singking'
        if mockpsfntype == None:
            mockpsfntype = 'doubking'
        if radispmrlbhl == None:
            radispmrlbhl = 2.
        if lablback == None:
            lablback = [r'$\mathcal{I}$', r'$\mathcal{D}$']
        if nameback == None:
            nameback = ['normisot', 'normfdfm']
        if strgenerunit == None:
            strgenerunit = r'GeV'
        if strgfluxunit == None:
            strgfluxunit = r'[1/cm$^2$/s/GeV]'
        if evttfull == None:
            evttfull = arange(4)
        binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
        if lablback == None:
            lablback = [r'$\mathcal{I}$', r'$\mathcal{D}$']
        if nameback == None:
            nameback = ['normisot', 'normfdfm']
        if anglassc == None:
            anglassc = deg2rad(0.5)

    ## Chandra and SDSS
    if exprtype == 'chan' or exprtype == 'sdss':
        if evttfull == None:
            evttfull = arange(4)
        if lablback == None:
            lablback = [r'$\mathcal{I}$']
        if nameback == None:
            nameback = ['normisot']
        if maxmgang == None:
            maxmgang = deg2rad(100. / 3600.)
        if anglassc == None:
            anglassc = deg2rad(0.5 / 3600.)
    
    ## Chandra
    if exprtype == 'chan':
        if psfntype == None:
            psfntype = 'singgaus'
        if strgenerunit == None:
            strgenerunit = r'KeV'
        if strgfluxunit == None:
            strgfluxunit = r'[1/cm$^2$/s/KeV]'
    
    ## SDSS
    if exprtype == 'sdss':
        if minmsind == None:
            minmsind = array([1.])
        if maxmsind == None:
            maxmsind = array([3.])
        if sinddistmean == None:
            sinddistmean = array([2.2])
        if sinddiststdv == None:
            sinddiststdv = array([0.5])
        if minmmeanpnts == None:
            minmmeanpnts = array([1.])
        if maxmmeanpnts == None:
            maxmmeanpnts = array([1e3])
        if psfntype == None:
            psfntype = 'doubgaus'
        if mockpsfntype == None:
            mockpsfntype = 'doubgaus'
        if radispmrlbhl == None:
            radispmrlbhl = 2. / 3600.
        if exprtype == 'sdss':
            strgfluxunit = '[nMgy]'
        if strgfluxunit == None:
            strgfluxunit = r'[mMag]'

    if maxmgangmarg == None:
        maxmgangmarg = margfact * maxmgang
        
    if maxmangl == None:
        maxmangl = 3. * maxmgangmarg
    
    if minmmeanpnts == None:
        minmmeanpnts = zeros(numbpopl) + 1.
    if maxmmeanpnts == None:
        maxmmeanpnts = zeros(numbpopl) + 1e3
        
    if minmfluxdistslop == None:
        minmfluxdistslop = zeros(numbpopl) + 0.5
    if maxmfluxdistslop == None:
        maxmfluxdistslop = zeros(numbpopl) + 3.5
    if minmfluxdistbrek == None:
        minmfluxdistbrek = minmflux
    if maxmfluxdistbrek == None:
        maxmfluxdistbrek = maxmflux
    if minmfluxdistsloplowr == None:
        minmfluxdistsloplowr = minmfluxdistslop
    if maxmfluxdistsloplowr == None:
        maxmfluxdistsloplowr = maxmfluxdistslop
    if minmfluxdistslopuppr == None:
        minmfluxdistslopuppr = minmfluxdistslop
    if maxmfluxdistslopuppr == None:
        maxmfluxdistslopuppr = maxmfluxdistslop
    
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

    ## axes
    ### energy
    gdat.binsenerfull = binsenerfull
    ### event class
    gdat.evttfull = evttfull

    ## boolean flag to evaluate a binned flux prior 
    gdat.bindprio = bindprio

    ## data type
    ###- mock - mock data
    ###- inpt - input data
    gdat.datatype = datatype
    
    ## likelihood function type
    gdat.liketype = liketype
    
    ## experiment type
    gdat.exprtype = exprtype
    
    ## PSF model type
    gdat.psfntype = psfntype
    
    ## flag to turn off PSF parameter updates
    gdat.boolproppsfn = boolproppsfn
    gdat.boolpropfluxdist = boolpropfluxdist

    ## input data
    ### measured data
    gdat.strgexpr = strgexpr
    ### background
    gdat.strgback = strgback
    gdat.lablback = lablback
    gdat.nameback = nameback
    ### exposure
    gdat.strgexpo = strgexpo

    ### plotting strings
    #### flux units
    gdat.maxmgangmarg = maxmgangmarg
    gdat.maxmangl = maxmangl
    gdat.strgenerunit = strgenerunit
    gdat.strgfluxunit = strgfluxunit

    ## PS parameter distribution models
    ### spatial
    gdat.spatdisttype = spatdisttype
    ### flux
    gdat.fluxdisttype = fluxdisttype
    ### color
    gdat.sinddisttype = sinddisttype
    
    ## PS spectral model
    gdat.spectype = spectype
    
    ## initial state setup
    ### number of point sources
    gdat.initnumbpnts = initnumbpnts
    gdat.initmeanpnts = initmeanpnts
    gdat.initfluxdistslop = initfluxdistslop
    gdat.initfluxdistbrek = initfluxdistbrek
    gdat.initfluxdistsloplowr = initfluxdistsloplowr
    gdat.initfluxdistslopuppr = initfluxdistslopuppr
    ### flag to draw the initial state from the prior
    gdat.randinit = randinit

    ## energy bins to be included
    gdat.indxenerincl = indxenerincl
    
    ## PSF bins to be included
    gdat.indxevttincl = indxevttincl
    
    ## number of Point Source (PS) populations
    gdat.numbpopl = numbpopl
    
    ## maximum angle from the PSs to evaluate the likelihood
    gdat.maxmangl = maxmangl
    gdat.maxmangleval = maxmangleval
    
    ## the maximum approximation error tolerated in units of the minimum flux allowed by the model
    gdat.specfraceval = specfraceval

    ## parameter limits
    ### background normalizations
    gdat.minmnormback = minmnormback
    gdat.maxmnormback = maxmnormback

    ### hyperparameters
    #### mean number of PS
    gdat.minmmeanpnts = minmmeanpnts
    gdat.maxmmeanpnts = maxmmeanpnts
    #### FDF power law index
    gdat.minmfluxdistslop = minmfluxdistslop
    gdat.maxmfluxdistslop = maxmfluxdistslop
    #### FDF break flux
    gdat.minmfluxdistbrek = minmfluxdistbrek
    gdat.maxmfluxdistbrek = maxmfluxdistbrek
    #### FDF lower power law index
    gdat.minmfluxdistsloplowr = minmfluxdistsloplowr
    gdat.maxmfluxdistsloplowr = maxmfluxdistsloplowr
    #### FDF upper power law index
    gdat.minmfluxdistslopuppr = minmfluxdistslopuppr
    gdat.maxmfluxdistslopuppr = maxmfluxdistslopuppr
    
    ### number of PS
    gdat.maxmnumbpnts = maxmnumbpnts
    
    ### PS parameters
    #### flux
    gdat.minmflux = minmflux
    gdat.maxmflux = maxmflux
    #### spectral index
    gdat.minmsind = minmsind
    gdat.maxmsind = maxmsind
    gdat.sinddistmean = sinddistmean
    gdat.sinddiststdv = sinddiststdv
    
    ## Region of interest
    ### image center
    gdat.lgalcntr = lgalcntr
    gdat.bgalcntr = bgalcntr
    ### half of the image size
    gdat.maxmgang = maxmgang
    ### the factor by which the model image is larger than the data image
    gdat.margfact = margfact
    
    ## angular radius in which associations can be made
    gdat.anglassc = anglassc

    ## proposals
    ### proposal scales
    gdat.stdvmeanpnts = stdvmeanpnts
    gdat.stdvfluxdistslop = stdvfluxdistslop
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
        
        ## pixelization type for mock data
        gdat.mockpixltype = mockpixltype
    
        gdat.mocknumbpopl = mocknumbpopl
        
        ### mock PSF model type
        gdat.mockpsfntype = mockpsfntype
        
        ### background normalization
        gdat.mocknormback = mocknormback
        
        ### number of PS
        gdat.mocknumbpnts = mocknumbpnts
        
        ### spatial distribution
        gdat.mockspatdisttype = mockspatdisttype

        ### flux distribution
        gdat.mockfluxdisttype = mockfluxdisttype
        gdat.mockfluxdistslop = mockfluxdistslop
        gdat.mockfluxdistbrek = mockfluxdistbrek
        gdat.mockfluxdistsloplowr = mockfluxdistsloplowr
        gdat.mockfluxdistslopuppr = mockfluxdistslopuppr
        
        ### color distribution
        gdat.mocksinddisttype = mocksinddisttype
        gdat.mocksinddistmean = mocksinddistmean
        gdat.mocksinddiststdv = mocksinddiststdv

        gdat.mockspectype = mockspectype

        ### flag to position mock point sources at the image center
        gdat.pntscntr = pntscntr
        ### mock image resolution
        gdat.mocknumbsidecart = mocknumbsidecart
        gdat.mocknumbsideheal = mocknumbsideheal

    ## proposal frequencies
    gdat.probprop = probprop

    ## the path of the input data
    gdat.pathdata = pathdata
    
    # temp
    # check inputs
    
    # date and time
    gdat.strgtime = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    if gdat.verbtype > 0:
        print 'PCAT started at ', gdat.strgtime
        print 'Initializing...'
    
    # check the call stack for the name of the configuring function
    gdat.strgcnfg = inspect.stack()[1][3]
    if gdat.verbtype > 0:
        print 'Configuration %s' % gdat.strgcnfg

    # setup the sampler
    setp(gdat) 

    if gdat.verbtype > 1:
        print 'minmflux'
        print gdat.minmflux
        print 'maxmflux'
        print gdat.maxmflux
        print 'minmcnts'
        print gdat.minmcnts
        print 'maxmcnts'
        print gdat.maxmcnts
        print 'probprop: '
        print vstack((arange(gdat.numbprop), gdat.strgprop, gdat.probprop)).T
        print 'indxsampnumbpnts: ', gdat.indxsampnumbpnts
        print 'indxsampmeanpnts: ', gdat.indxsampmeanpnts
        print 'indxsampfluxdistslop: ', gdat.indxsampfluxdistslop
        print 'indxsampfluxdistbrek: ', gdat.indxsampfluxdistbrek
        print 'indxsampfluxdistsloplowr: ', gdat.indxsampfluxdistsloplowr
        print 'indxsampfluxdistslopuppr: ', gdat.indxsampfluxdistslopuppr
        print 'indxsamppsfipara: ', gdat.indxsamppsfipara
        print 'indxsampnormback: '
        print gdat.indxsampnormback
        print 'indxsampcompinit: ', gdat.indxsampcompinit
        print 'maxmangleval'
        print gdat.maxmangleval
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
            print 'truenumbpnts: ', gdat.truenumbpnts
            print 'truepsfipara'
            print gdat.truepsfipara
            print
            if gdat.datatype == 'mock':
                print 'Mock data'
                print 'mocknumbpnts: ', gdat.mocknumbpnts
                print 'mockspatdisttype'
                print mockspatdisttype
                print 'mockfluxdistslop: ', gdat.mockfluxdistslop
                print 'mockfluxdistbrek: ', gdat.mockfluxdistbrek
                print 'mockfluxdistsloplowr: ', gdat.mockfluxdistsloplowr
                print 'mockfluxdistslopuppr: ', gdat.mockfluxdistslopuppr
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
        if gdat.truepsfn != None:
            plot_eval(gdat)
        #if gdat.datatype == 'mock':
        #    plot_pntsdiff()

    if gdat.verbtype > 0:
        print 'Sampling...'
    
    if gdat.verbtype > 1:
        show_memo(gdat, 'gdat')

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
    listindxparamodi = zeros((gdat.numbswep, gdat.numbproc), dtype=int) - 1
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
    ## PSF parameters
    listpsfipara = listsampvarb[:, :, gdat.indxsamppsfipara].reshape(gdat.numbsamptotl, -1)
    ## Background normalization
    listnormback = listsampvarb[:, :, gdat.indxsampnormback].reshape(gdat.numbsamptotl, gdat.numbback, gdat.numbener)
    ## mean number of PS
    listmeanpnts = listsampvarb[:, :, gdat.indxsampmeanpnts].reshape(gdat.numbsamptotl, -1)
    ## flux distribution
    listfluxdistslop = listsampvarb[:, :, gdat.indxsampfluxdistslop].reshape(gdat.numbsamptotl, gdat.numbpopl)
    listfluxdistbrek = listsampvarb[:, :, gdat.indxsampfluxdistbrek].reshape(gdat.numbsamptotl, gdat.numbpopl)
    listfluxdistsloplowr = listsampvarb[:, :, gdat.indxsampfluxdistsloplowr].reshape(gdat.numbsamptotl, gdat.numbpopl)
    listfluxdistslopuppr = listsampvarb[:, :, gdat.indxsampfluxdistslopuppr].reshape(gdat.numbsamptotl, gdat.numbpopl)
    ## number of PS
    listnumbpnts = listsampvarb[:, :, gdat.indxsampnumbpnts].astype(int).reshape(gdat.numbsamptotl, -1)

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

    # bin the posterior
    pntsprob = zeros((gdat.numbpopl, gdat.numbener, gdat.numbpixl, gdat.numbflux))
    for k in gdat.indxsamp:
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for h in gdat.indxflux:
                    indxpnts = where((gdat.binsspec[i, h] < listspec[l][k, i, :]) & (listspec[l][k, i, :] < gdat.binsspec[i, h+1]))[0]
                    hpixl = retr_indxpixl(gdat, listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    # Gelman-Rubin test
    gmrbstat = zeros(gdat.numbpixlsave)
    if gdat.numbproc > 1:
        
        if gdat.verbtype > 0:
            print 'Computing the Gelman-Rubin TS...'
            tim0 = time.time()
    
        for n in range(gdat.numbpixlsave):
            gmrbstat[n] = tdpy.mcmc.gmrb_test(listmodlcnts[:, :, n])
        
        if gdat.verbtype > 0:
            tim1 = time.time()
            print 'Done in %.3g seconds.' % (tim1 - tim0)


    # calculate the autocorrelation of the chains
    if gdat.verbtype > 0:
        print 'Computing the autocorrelation of the chains...'
        tim0 = time.time()
   
    atcr, timeatcr = tdpy.mcmc.retr_atcr(listmodlcnts)

    if gdat.verbtype > 0:
        tim1 = time.time()
        print 'Done in %.3g seconds.' % (tim1 - tim0)

    # write the PCAT output to disc
    pathpcatlite = gdat.pathdata + '/pcatlite_' + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.fits'  
    pathpcat = gdat.pathdata + '/pcat_' + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.fits'  
    
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
    
    ## axes
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
    
    head['maxmangl'] = gdat.maxmangl
    head['margfact'] = gdat.margfact
    head['strgcnfg'] = gdat.strgcnfg
    head['strgtime'] = gdat.strgtime
    
    head['timeatcr'] = timeatcr
   
    ## proposal scales
    ### parameter updates
    head['stdvmeanpnts'] = gdat.stdvmeanpnts
    head['stdvfluxdistslop'] = gdat.stdvfluxdistslop
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

    ## approximation error in units of the minimum flux
    head['specfraceval'] = gdat.specfraceval

    ## index of the energy bin, where the FDF is defined
    head['indxenerfluxdist'] = gdat.indxenerfluxdist[0]

    # PS parameter distributions
    for l in gdat.indxpopl:
        head['spatdisttypepop%d' % l] = gdat.spatdisttype[l]
        head['fluxdisttypepop%d' % l] = gdat.fluxdisttype[l]
        head['sinddisttypepop%d' % l] = gdat.sinddisttype[l]
        head['spectypepop%d' % l] = gdat.spectype[l]
    
    ## mock data
    if gdat.datatype == 'mock':
        head['mocknumbpopl'] = gdat.mocknumbpopl
        head['mockpsfntype'] = gdat.mockpsfntype
        for l in gdat.indxpopl:
            head['mockspatdisttypepop%d' % l] = gdat.mockspatdisttype[l]
            head['mockfluxdisttypepop%d' % l] = gdat.mockfluxdisttype[l]
            head['mocksinddisttypepop%d' % l] = gdat.mocksinddisttype[l]
            head['mockspectypepop%d' % l] = gdat.mockspectype[l]

    head['strgexpr'] = gdat.strgexpr
    for k in gdat.indxback:
        head['strgback%04d' % k] = gdat.strgback[k]
        head['lablback%04d' % k] = gdat.lablback[k]
        head['nameback%04d' % k] = gdat.nameback[k]
    head['strgexpo'] = gdat.strgexpo
    
    head['pathdata'] = gdat.pathdata
    
    head['levi'] = levi
    head['info'] = info
    
    listhdun = []
    listhdun.append(pf.PrimaryHDU(header=head))

    ## posterior
    ### number of PS
    listhdun.append(pf.ImageHDU(listnumbpnts))
    listhdun[-1].header['EXTNAME'] = 'numbpnts'

    ### hyperparameters
    listhdun.append(pf.ImageHDU(listmeanpnts))
    listhdun[-1].header['EXTNAME'] = 'meanpnts'
   
    listhdun.append(pf.ImageHDU(listfluxdistslop))
    listhdun[-1].header['EXTNAME'] = 'fluxdistslop'
    
    listhdun.append(pf.ImageHDU(listfluxdistbrek))
    listhdun[-1].header['EXTNAME'] = 'fluxdistbrek'
    
    listhdun.append(pf.ImageHDU(listfluxdistsloplowr))
    listhdun[-1].header['EXTNAME'] = 'fluxdistsloplowr'
    
    listhdun.append(pf.ImageHDU(listfluxdistslopuppr))
    listhdun[-1].header['EXTNAME'] = 'fluxdistslopuppr'
    
    ### PSF parameters
    listhdun.append(pf.ImageHDU(listpsfipara))
    listhdun[-1].header['EXTNAME'] = 'psfipara'
    
    ### background normalizations
    listhdun.append(pf.ImageHDU(listnormback))
    listhdun[-1].header['EXTNAME'] = 'normback'
   
    ### PS parameters
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

        #### save the posterior positions as a CSV file
        path = gdat.pathdata + '/pcat_lgalpop%d_' % l + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.csv'  
        savetxt(path, listlgal[l], fmt='%7.5g', delimiter=',')
        path = gdat.pathdata + '/pcat_bgalpop%d_' % l + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '.csv'  
        savetxt(path, listbgal[l], fmt='%7.5g', delimiter=',')


    ### log-likelihood
    listhdun.append(pf.ImageHDU(listllik))
    listhdun[-1].header['EXTNAME'] = 'llik'
    
    ### log-prior
    listhdun.append(pf.ImageHDU(listlpri))
    listhdun[-1].header['EXTNAME'] = 'lpri'
   
    ## store the lite file
    pf.HDUList(listhdun).writeto(pathpcatlite, clobber=True, output_verify='ignore')

    ## unit sample
    listhdun.append(pf.ImageHDU(listsamp))
    listhdun[-1].header['EXTNAME'] = 'listsamp'
    
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
    
    listhdun.append(pf.ImageHDU(gdat.minmmeanpnts))
    listhdun[-1].header['EXTNAME'] = 'minmmeanpnts'
    
    listhdun.append(pf.ImageHDU(gdat.maxmmeanpnts))
    listhdun[-1].header['EXTNAME'] = 'maxmmeanpnts'
    
    listhdun.append(pf.ImageHDU(gdat.minmfluxdistslop))
    listhdun[-1].header['EXTNAME'] = 'minmfluxdistslop'
    
    listhdun.append(pf.ImageHDU(gdat.maxmfluxdistslop))
    listhdun[-1].header['EXTNAME'] = 'maxmfluxdistslop'
    
    listhdun.append(pf.ImageHDU(gdat.minmfluxdistbrek))
    listhdun[-1].header['EXTNAME'] = 'minmfluxdistbrek'
    
    listhdun.append(pf.ImageHDU(gdat.maxmfluxdistbrek))
    listhdun[-1].header['EXTNAME'] = 'maxmfluxdistbrek'
    
    listhdun.append(pf.ImageHDU(gdat.minmfluxdistsloplowr))
    listhdun[-1].header['EXTNAME'] = 'minmfluxdistsloplowr'
    
    listhdun.append(pf.ImageHDU(gdat.maxmfluxdistsloplowr))
    listhdun[-1].header['EXTNAME'] = 'maxmfluxdistsloplowr'
    
    listhdun.append(pf.ImageHDU(gdat.minmfluxdistslopuppr))
    listhdun[-1].header['EXTNAME'] = 'minmfluxdistslopuppr'
    
    listhdun.append(pf.ImageHDU(gdat.maxmfluxdistslopuppr))
    listhdun[-1].header['EXTNAME'] = 'maxmfluxdistslopuppr'
    
    listhdun.append(pf.ImageHDU(gdat.sinddistmean))
    listhdun[-1].header['EXTNAME'] = 'sinddistmean'
    
    listhdun.append(pf.ImageHDU(gdat.sinddiststdv))
    listhdun[-1].header['EXTNAME'] = 'sinddiststdv'
    
    listhdun.append(pf.ImageHDU(gdat.binsenerfull))
    listhdun[-1].header['EXTNAME'] = 'binsenerfull'

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
    ## autocorrelation
    listhdun.append(pf.ImageHDU(atcr))
    listhdun[-1].header['EXTNAME'] = 'atcr'
    
    ## convergence diagnostic
    listhdun.append(pf.ImageHDU(gmrbstat))
    listhdun[-1].header['EXTNAME'] = 'gmrbstat'
    
    ## PDF of PS positions
    listhdun.append(pf.ImageHDU(pntsprob))
    listhdun[-1].header['EXTNAME'] = 'pntsprob'
   
    # mock data
    if gdat.datatype == 'mock':
        listhdun.append(pf.ImageHDU(gdat.truenumbpnts))
        listhdun[-1].header['EXTNAME'] = 'mocknumbpnts'

        # PS parameters
        for l in gdat.indxpopl:
            listhdun.append(pf.ImageHDU(gdat.truelgal[l]))
            listhdun[-1].header['EXTNAME'] = 'mocklgalpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truebgal[l]))
            listhdun[-1].header['EXTNAME'] = 'mockbgalpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truespec[l]))
            listhdun[-1].header['EXTNAME'] = 'mockspecpop%d' % l

            listhdun.append(pf.ImageHDU(gdat.truesind[l]))
            listhdun[-1].header['EXTNAME'] = 'mocksindpop%d' % l

        ## hyperparameters
        listhdun.append(pf.ImageHDU(gdat.mockfluxdistslop))
        listhdun[-1].header['EXTNAME'] = 'mockfluxdistslop'

        listhdun.append(pf.ImageHDU(gdat.mockfluxdistbrek))
        listhdun[-1].header['EXTNAME'] = 'mockfluxdistbrek'

        listhdun.append(pf.ImageHDU(gdat.mockfluxdistsloplowr))
        listhdun[-1].header['EXTNAME'] = 'mockfluxdistsloplowr'

        listhdun.append(pf.ImageHDU(gdat.mockfluxdistslopuppr))
        listhdun[-1].header['EXTNAME'] = 'mockfluxdistslopuppr'

        listhdun.append(pf.ImageHDU(gdat.mocksinddistmean))
        listhdun[-1].header['EXTNAME'] = 'mocksinddistmean'

        listhdun.append(pf.ImageHDU(gdat.mocksinddiststdv))
        listhdun[-1].header['EXTNAME'] = 'mocksinddiststdv'

        ## PSF parameters
        listhdun.append(pf.ImageHDU(gdat.truepsfipara))
        listhdun[-1].header['EXTNAME'] = 'mockpsfipara'
        
        ## Background normalizations
        listhdun.append(pf.ImageHDU(gdat.truenormback))
        listhdun[-1].header['EXTNAME'] = 'mocknormback'

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
    
    # PSF radial profile
    gdatmodi.thispsfn = gdatmodi.thispsfnintp(gdat.binsangl)

    # PSF FWHM
    gdatmodi.thisfwhm = 2. * retr_psfnwdth(gdat, gdatmodi.thispsfn, 0.5)
    
    # number of background counts per PSF
    gdatmodi.thisbackfwhmcnts = retr_backfwhmcnts(gdat, gdatmodi.thissampvarb[gdat.indxsampnormback], gdatmodi.thisfwhm)

    # number of counts and standard deviation of each PS
    gdatmodi.thiscnts = []
    gdatmodi.thissigm = []
    for l in gdat.indxpopl:
        indxpixltemp = retr_indxpixl(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]], gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]])
        cntstemp = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]][:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]
        retr_sigm(gdat, sum(cntstemp, 2), gdatmodi.thisbackfwhmcnts)
        gdatmodi.thiscnts.append(cntstemp)
        gdatmodi.thissigm.append(cntstemp)
    
    # standard deviation axis
    gdatmodi.binssigm = retr_sigm(gdat, gdat.binscnts, gdatmodi.thisbackfwhmcnts)
    
    # plots
    ## PSF radial profile
    # temp
    if gdat.strgcnfg != 'cnfg_test':
        plot_psfn(gdat, gdatmodi)
    
    ## PSF FWHM
    plot_fwhm(gdat, gdatmodi)
    
    # number of background counts per PSF
    for i in gdat.indxener:
        path = gdat.pathplot + 'backfwhmcntsflux%d_%09d.pdf' % (i, gdatmodi.cntrswep)
        tdpy.util.plot_heal(path, sum(gdatmodi.thisbackfwhmcnts, 2)[i, :], pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, \
                                                                              minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)

    # temp -- list may not be the ultimate solution to copy gdatmodi.thisindxpntsfull
    temppntsflux, temppntscnts, tempmodlflux, tempmodlcnts = retr_maps(gdat, list(gdatmodi.thisindxpntsfull), copy(gdatmodi.thissampvarb))
    gdatmodi.thispntscnts = gdatmodi.thispntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    gdatmodi.thiserrrpnts = 100. * (gdatmodi.thispntscnts - temppntscnts) / temppntscnts

    # temp
    if False:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = gdat.pathplot + 'temppntsflux%d%d_%09d.pdf' % (i, m, gdatmodi.cntrswep)
                tdpy.util.plot_heal(path, temppntsflux[i, :, m], pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, \
                                                                                  minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = gdat.pathplot + 'thispntsflux%d%d_%09d.pdf' % (i, m, gdatmodi.cntrswep)
                tdpy.util.plot_heal(path, gdatmodi.thispntsflux[i, :, m], pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, \
                                                                              minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            path = gdat.pathplot + 'errrpntsabsl%d%d_%09d.pdf' % (i, m, gdatmodi.cntrswep)
            tdpy.util.plot_heal(path, gdatmodi.thispntscnts[i, :, m] - temppntscnts[i, :, m], pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, \
                    numbpixl=gdat.numbpixlheal, resi=True, minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
    
    gdatmodi.indxtruepntsassc = []
    for l in gdat.indxpopl:
        if gdat.trueinfo:
            lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]]
            bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]]
            spec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]]
            indxmodl, indxtruepntsassc = corr_catl(gdat, l, lgal, bgal, spec)
            gdatmodi.indxtruepntsassc.append(indxtruepntsassc)
            gdatmodi.thisspecmtch = copy(spec[:, indxmodl])
            gdatmodi.thisspecmtch[:, indxtruepntsassc.miss] = 0.
            # temp
            if gdat.strgcnfg != 'cnfg_test':
                plot_scatspec(gdat, l, gdatmodi=gdatmodi)
            
        # temp
        if gdat.strgcnfg != 'cnfg_test':
            plot_histspec(gdat, l, gdatmodi=gdatmodi)
            plot_histsind(gdat, l, gdatmodi=gdatmodi)
            plot_fluxsind(gdat, l, gdatmodi=gdatmodi)
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
   
    if gdat.verbtype > 1:
        'Variables owned by processes'
        print 'listsamp'
        print sys.getsizeof(listsamp) / 2.**20, 'MB'
        print 'listsampvarb'
        print sys.getsizeof(listsampvarb) / 2.**20, 'MB'
        print 'listmodlcnts'
        print sys.getsizeof(listmodlcnts) / 2.**20, 'MB'
        print 'listchrollik'
        print sys.getsizeof(gdatmodi.listchrollik) / 2.**20, 'MB'
        print 

    gdat.listauxipara = zeros((gdat.numbswep, gdat.numbcompcolr))
    gdat.listlaccfrac = zeros(gdat.numbswep)
    gdat.listnumbpair = zeros(gdat.numbswep)
    gdat.listjcbnfact = zeros(gdat.numbswep)
    gdat.listcombfact = zeros(gdat.numbswep)

    gdatmodi.cntrswep = 0
    
    # initialize the chain
    retr_llik(gdat, gdatmodi, init=True)
    retr_lpri(gdat, gdatmodi, init=True)

    # current sample index
    thiscntr = -1
    
    while gdatmodi.cntrswep < gdat.numbswep:
        
        timeinit = time.time()
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) and gdat.makeplot
        gdatmodi.boolreje = False
    
        # choose a proposal type
        retr_thisindxprop(gdat, gdatmodi)
            
        # save the proposal type
        listindxprop[gdatmodi.cntrswep] = gdatmodi.thisindxprop
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
        listchrototl[gdatmodi.cntrswep, 1] = timefinl - timebegn

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
        if gdatmodi.thisindxprop == gdat.indxpropfluxdistslop or gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or \
                                              gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                if gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, 1] < 0. or gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, 1] > 1.:
                    gdatmodi.boolreje = True
        if (gdatmodi.thisindxprop == gdat.indxpropmeanpnts or gdatmodi.thisindxprop >= gdat.indxproppsfipara and gdatmodi.thisindxprop != gdat.indxpropbrth and \
                                                                                                        gdatmodi.thisindxprop != gdat.indxpropdeth) and not gdatmodi.boolreje:
            if where((gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 1] < 0.) | (gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 1] > 1.))[0].size > 0:
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
            listchrototl[gdatmodi.cntrswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(gdat, gdatmodi) 
            timefinl = time.time()
            listchrototl[gdatmodi.cntrswep, 3] = timefinl - timebegn
   
            # temp
            if gdat.strgcnfg == 'cnfg_test':
                temp = zeros(gdat.numbpixl)
                temp[gdat.indxpixlmodi] = 1.
                path = gdat.pathplot + 'indxpixlmodi_%09d.pdf' % gdatmodi.cntrswep
                tdpy.util.plot_heal(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, pixltype=gdat.pixltype, \
                                                                                  minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
            
                temp = zeros(gdat.numbpixl)
                temp[gdat.indxpixlmodi] = gdatmodi.thispntsflux[gdat.indxcubemodi]
                path = gdat.pathplot + 'thispntsflux_%09d.pdf' % gdatmodi.cntrswep
                tdpy.util.plot_heal(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, pixltype=gdat.pixltype, \
                                                                                  minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
            
                temp = zeros(gdat.numbpixl)
                temp[gdat.indxpixlmodi] = gdatmodi.nextpntsflux[gdat.indxcubemodi]
                path = gdat.pathplot + 'nextpntsflux_%09d.pdf' % gdatmodi.cntrswep
                tdpy.util.plot_heal(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, pixltype=gdat.pixltype, \
                                                                                  minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
            
                temp = zeros(gdat.numbpixl)
                temp[gdat.indxpixlmodi] = gdatmodi.nextpntsflux[gdat.indxcubemodi] - gdatmodi.thispntsflux[gdat.indxcubemodi]
                path = gdat.pathplot + 'diffpntsflux_%09d.pdf' % gdatmodi.cntrswep
                tdpy.util.plot_heal(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlheal, pixltype=gdat.pixltype, \
                                                                                  minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
            
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

            listaccp[gdatmodi.cntrswep] = True

        # reject the sample
        else:

            if gdat.verbtype > 1:
                print 'Rejected.'

            listaccp[gdatmodi.cntrswep] = False
        
        if gdatmodi.thisindxprop < gdat.indxpropbrth:
            listindxparamodi[gdatmodi.cntrswep] = gdatmodi.indxsampvarbmodi

        # sanity checks
        indxsampbadd = where((gdatmodi.drmcsamp[gdat.numbpopl:, 0] > 1.) | (gdatmodi.drmcsamp[gdat.numbpopl:, 0] < 0.))[0] + 1
        if indxsampbadd.size > 0:
            print 'Unit sample vector went outside [0,1]!'
            print 'cntrswep'
            print gdatmodi.cntrswep
            print 'thisindxprop'
            print gdat.strgprop[gdatmodi.thisindxprop]
            print 'indxsampbadd'
            print indxsampbadd
            print 'drmcsamp'
            print gdatmodi.drmcsamp[indxsampbadd, :]
            return
            
        for l in gdat.indxpopl:
            indxtemp = where(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]] < gdat.minmflux)[0]
            if indxtemp.size > 0:
                print 'Spectrum of some PS went below the prior range!'
                print 'indxtemp'
                print indxtemp
                print 'flux'
                print gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], indxtemp]]
                print 
            indxtemp = where(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]] > gdat.maxmflux)[0]
            if indxtemp.size > 0:
                print 'Spectrum of some PS went above the prior range!'          
                print 'indxtemp'
                print indxtemp
                print 'flux'
                print gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], indxtemp]]
                print 

        # save the sample
        if boolsave[gdatmodi.cntrswep]:
            listsamp[indxsampsave[gdatmodi.cntrswep], :] = gdatmodi.drmcsamp[:, 0]
            listsampvarb[indxsampsave[gdatmodi.cntrswep], :] = gdatmodi.thissampvarb
            listmodlcnts[indxsampsave[gdatmodi.cntrswep], :] = gdatmodi.thismodlcnts[0, gdat.indxpixlsave, 0]
            listpntsfluxmean[indxsampsave[gdatmodi.cntrswep], :] = mean(sum(gdatmodi.thispntsflux * gdat.expo, 2) / sum(gdat.expo, 2), 1)
            listindxpntsfull.append(copy(gdatmodi.thisindxpntsfull))
            listllik[indxsampsave[gdatmodi.cntrswep]] = sum(gdatmodi.thisllik)
            listlpri[indxsampsave[gdatmodi.cntrswep]] = sum(gdatmodi.thislpri)
            
            lprinorm = 0.
            for l in gdat.indxpopl:
                # temp
                ## brok terms are not complete
                ## lpri calculation is turned off
                break
                numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[l]]
                meanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]]
                lpri += numbpnts * gdat.priofactlgalbgal + gdat.priofactfluxdistslop + gdat.priofactmeanpnts - log(meanpnts)
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]]
                if gdat.fluxdisttype[l] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]]
                    lpri -= log(1. + fluxdistslop**2)
                    lpri += sum(log(pdfn_flux_powr(gdat, flux, fluxdistslop)))
                if gdat.fluxdisttype[l] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]]
                    lpri += sum(log(pdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))
            listlprinorm[indxsampsave[gdatmodi.cntrswep]] = lprinorm
            
            if gdat.tracsamp:
                
                numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[gdatmodi.cntrswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            timefinl = time.time()
            listchrototl[gdatmodi.cntrswep, 0] = timefinl - timeinit

        # log the progress
        if gdat.verbtype > 0:
            thiscntr = tdpy.util.show_prog(gdatmodi.cntrswep, gdat.numbswep, thiscntr, indxprocwork=indxprocwork)
    
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
        gdatmodi.cntrswep += 1

    if gdat.verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    gdatmodi.listchrollik = array(gdatmodi.listchrollik)
    
    listchan = [listsamp, listsampvarb, listindxprop, listchrototl, listllik, listlpri, listaccp, listmodlcnts, listindxpntsfull, listindxparamodi, \
        gdat.listauxipara, gdat.listlaccfrac, gdat.listnumbpair, gdat.listjcbnfact, gdat.listcombfact, \
        listpntsfluxmean, gdatmodi.listchrollik]
    
    return listchan

