# common imports
from __init__ import *

# internal functions
from util import *
from visu import *

def work(globdata, indxprocwork):

    timereal = time.time()
    timeproc = time.clock()
    
    # re-seed the random number generator for the process
    seed()
    
    # construct the run tag
    globdata.rtag = retr_rtag(globdata, indxprocwork)
    
    # initialize the sample vector 
    if globdata.randinit or not globdata.trueinfo:
        if globdata.initnumbpnts != None:
            thisnumbpnts = globdata.initnumbpnts
        else:
            thisnumbpnts = empty(globdata.numbpopl, dtype=int)
            for l in globdata.indxpopl:
                thisnumbpnts[l] = choice(arange(globdata.minmnumbpnts, globdata.maxmnumbpnts[l] + 1))
    else:
        thisnumbpnts = globdata.truenumbpnts
        
    globdata.thisindxpntsfull = []
    globdata.thisindxpntsempt = []
    for l in globdata.indxpopl:
        globdata.thisindxpntsfull.append(range(thisnumbpnts[l]))
        globdata.thisindxpntsempt.append(range(thisnumbpnts[l], globdata.maxmnumbpnts[l]))
    globdata.thisindxsamplgal, globdata.thisindxsampbgal, globdata.thisindxsampspec, \
        globdata.thisindxsampsind, globdata.thisindxsampcomp = retr_indx(globdata, globdata.thisindxpntsfull)
    
    globdata.drmcsamp = zeros((globdata.maxmsampsize, 2))
    
    globdata.drmcsamp[globdata.indxsampnumbpnts, 0] = thisnumbpnts
    globdata.drmcsamp[globdata.indxsampfdfnnorm, 0] = rand(globdata.numbpopl)
    if globdata.trueinfo and globdata.datatype == 'mock':
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = cdfn_atan(globdata.truefdfnslop, globdata.minmfdfnslop, globdata.factfdfnslop)
    else:
        globdata.drmcsamp[globdata.indxsampfdfnslop, 0] = rand(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener))
    globdata.drmcsamp[globdata.indxsampnormback, 0] = rand(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener))
    if globdata.randinit or not globdata.trueinfo or globdata.truepsfipara == None:
        globdata.drmcsamp[globdata.indxsamppsfipara, 0] = rand(globdata.numbpsfipara)
    else:
        if globdata.psfntype == globdata.truepsfntype:
            for k in globdata.indxmodlpsfipara:
                globdata.drmcsamp[globdata.indxsamppsfipara[k], 0] = cdfn_psfipara(globdata, globdata.truepsfipara[k], k)
        else:
            globdata.drmcsamp[globdata.indxsamppsfipara, 0] = rand(globdata.numbpsfipara)
        
    for l in globdata.indxpopl:
        if globdata.randinit or not globdata.trueinfo:
            globdata.drmcsamp[globdata.thisindxsampcomp[l], 0] = rand(globdata.thisindxsampcomp[l].size)
        else:
            globdata.drmcsamp[globdata.thisindxsamplgal[l], 0] = copy(cdfn_self(globdata.truelgal[l], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg))
            globdata.drmcsamp[globdata.thisindxsampbgal[l], 0] = copy(cdfn_self(globdata.truebgal[l], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg))
            for i in globdata.indxenerfdfn:
                fdfnslop = icdf_atan(globdata.drmcsamp[globdata.indxsampfdfnslop[l, i], 0], globdata.minmfdfnslop, globdata.factfdfnslop)
                spec = cdfn_spec(globdata, globdata.truespec[l][0, i, :], fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])
                globdata.drmcsamp[globdata.thisindxsampspec[l][i, :], 0] = copy(spec)
            if globdata.colrprio:
                #globdata.drmcsamp[globdata.thisindxsampsind[l], 0] = copy(cdfn_atan(globdata.truesind[l], globdata.minmsind, globdata.factsind))
                globdata.drmcsamp[globdata.thisindxsampsind[l], 0] = cdfn_eerr(globdata.truesind[l], globdata.meansind[l], globdata.stdvsind[l], \
                                                                                        globdata.sindcdfnnormminm[l], globdata.sindcdfnnormdiff[l])
    globdata.thissampvarb, thisindxpixlpnts, thiscnts, globdata.thispntsflux, thismodlflux, \
        globdata.thismodlcnts = pars_samp(globdata, globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])
    
    if globdata.verbtype > 2:
        print 'thisindxpntsfull: ', globdata.thisindxpntsfull
        print 'thisindxpntsempt: ', globdata.thisindxpntsempt  
        print 'thisindxsamplgal: ', globdata.thisindxsamplgal
        print 'thisindxsampbgal: ', globdata.thisindxsampbgal
        print 'thisindxsampspec: '
        print globdata.thisindxsampspec
        if globdata.colrprio:
            print 'thisindxsampsind: ', globdata.thisindxsampsind
        print 'thisindxsampcomp: ', globdata.thisindxsampcomp

    globdata.nextpntsflux = zeros_like(globdata.thispntsflux)
    globdata.nextmodlflux = zeros_like(globdata.thispntsflux)
    globdata.nextmodlcnts = zeros_like(globdata.thispntsflux)
    globdata.nextllik = zeros_like(globdata.thispntsflux)

    globdata.nextsampvarb = copy(globdata.thissampvarb)
    
    if globdata.verbtype > 1:
        print 'thissampvarb: ', globdata.thissampvarb
        
    # sampler setup
    # auxiliary variable standard globdata.deviation for merge/split
    globdata.maxmdistpnts = 2. # [deg]
 
    listchan = rjmc(globdata, indxprocwork)
    
    timereal = time.time() - timereal
    timeproc = time.clock() - timeproc
    
    listchan.append(timereal)
    listchan.append(timeproc)
    
    return listchan




def init(cnfg):
    
    timetotlreal = time.time()
    timetotlproc = time.clock()
    
    globdata = globdatastrt()
    
    globdata.verbtype = cnfg['verbtype']
    globdata.plotperd = cnfg['plotperd']
    globdata.makeplot = cnfg['makeplot']
    globdata.diagsamp = cnfg['diagsamp']
    
    globdata.numbproc = cnfg['numbproc']
    globdata.numbswep = cnfg['numbswep']
    globdata.numbburn = cnfg['numbburn']
    globdata.factthin = cnfg['factthin']
    
    globdata.stdvfdfnnorm = cnfg['stdvfdfnnorm']
    globdata.stdvfdfnslop = cnfg['stdvfdfnslop']
    globdata.stdvpsfipara = cnfg['stdvpsfipara']
    globdata.stdvback = cnfg['stdvback']
    globdata.stdvlbhl = cnfg['stdvlbhl']
    globdata.stdvspec = cnfg['stdvspec']
    globdata.stdvpropsind = cnfg['stdvpropsind']
    globdata.spmrlbhl = cnfg['spmrlbhl']
    globdata.fracrand = cnfg['fracrand']
    
    globdata.datatype = cnfg['datatype']
    
    
    globdata.psfntype = cnfg['psfntype']
    
    globdata.liketype = cnfg['liketype']
    globdata.exprtype = cnfg['exprtype']
    globdata.pixltype = cnfg['pixltype']
    
    globdata.regitype = cnfg['regitype']
    globdata.randinit = cnfg['randinit']
    
    globdata.fdfntype = cnfg['fdfntype']

    globdata.proppsfn = cnfg['proppsfn']
    globdata.colrprio = cnfg['colrprio']
    globdata.numbpopl = cnfg['numbpopl']
    globdata.indxenerincl = cnfg['indxenerincl']
    globdata.indxevttincl = cnfg['indxevttincl']
    
    globdata.maxmangleval = cnfg['maxmangleval']

    globdata.minmspec = cnfg['minmspec']
    globdata.maxmspec = cnfg['maxmspec']
    if globdata.colrprio:
        globdata.minmsind = cnfg['minmsind']
        globdata.maxmsind = cnfg['maxmsind']
        globdata.meansind = cnfg['meansind']
        globdata.stdvsind = cnfg['stdvsind']
   

    globdata.minmfdfnnorm = cnfg['minmfdfnnorm']
    globdata.maxmfdfnnorm = cnfg['maxmfdfnnorm']
    globdata.minmfdfnslop = cnfg['minmfdfnslop']
    globdata.maxmfdfnslop = cnfg['maxmfdfnslop']
    

    globdata.maxmnormback = cnfg['maxmnormback']
    globdata.minmnormback = cnfg['minmnormback']
    
    if globdata.datatype == 'mock':
        globdata.mocknumbpnts = cnfg['mocknumbpnts']
        globdata.mockfdfnslop = cnfg['mockfdfnslop']
        globdata.mocknormback = cnfg['mocknormback']
        globdata.mockpsfntype = cnfg['mockpsfntype']

    globdata.maxmnumbpnts = cnfg['maxmnumbpnts']
    
    globdata.probprop = cnfg['probprop']
   
    # input
    ## experimental data
    globdata.strgexpr = cnfg['strgexpr']
    ## background
    globdata.strgback = cnfg['strgback']
    globdata.lablback = cnfg['lablback']
    globdata.nameback = cnfg['nameback']
    ## exposure
    globdata.strgexpo = cnfg['strgexpo']
    
    globdata.initnumbpnts = cnfg['initnumbpnts']
    globdata.trueinfo = cnfg['trueinfo']
    
    globdata.margsize = cnfg['margsize']
    globdata.maxmgang = cnfg['maxmgang']
    globdata.lgalcntr = cnfg['lgalcntr']
    globdata.bgalcntr = cnfg['bgalcntr']
    
    globdata.numbsideheal = cnfg['numbsideheal']
    globdata.numbsidecart = cnfg['numbsidecart']
    
    if globdata.colrprio:
        if globdata.minmspec.size != 1 or globdata.minmspec.size != 1:
            print 'Spectral limits must be numpy scalars when color priors are used!'
    
    # date and time
    globdata.strgtime = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
    if globdata.verbtype > 0:
        print 'PCAT started at ', globdata.strgtime
        print 'Initializing...'
    
    # setup the sampler
    setp(globdata) 

    if globdata.verbtype > 1:
        print 'probprop: '
        print vstack((arange(globdata.numbprop), globdata.strgprop, globdata.probprop)).T
        print 'indxsampnumbpnts: ', globdata.indxsampnumbpnts
        print 'indxsampfdfnnorm: ', globdata.indxsampfdfnnorm
        print 'indxsampfdfnslop: ', globdata.indxsampfdfnslop
        print 'indxsamppsfipara: ', globdata.indxsamppsfipara
        print 'indxsampnormback: '
        print globdata.indxsampnormback
        print 'indxsampcompinit: ', globdata.indxsampcompinit
        print 'maxmangleval'
        print globdata.maxmangleval
        if globdata.trueinfo:
            print 'truelgal: ', globdata.truelgal
            print 'truebgal: ', globdata.truebgal
            print 'truespec: '
            print globdata.truespec
            print 'truenumbpnts: ', globdata.truenumbpnts
            print 'truefdfnslop: ', globdata.truefdfnslop
            print 'truenormback: '
            print globdata.truenormback
            print
            if globdata.datatype == 'mock':
                print 'mocknumbpnts: ', globdata.mocknumbpnts
                print 'mockfdfnslop: ', globdata.mockfdfnslop
                print 'mockpsfipara: '
                print globdata.mockpsfipara
                print 'mocknormback: '
                print globdata.mocknormback


    # make initial plots
    if globdata.makeplot:
        # temp
        #plot_3fgl_thrs(globdata)
        #plot_datacntshist()
        if globdata.exprtype == 'ferm':
            plot_fgl3(globdata)
        #plot_intr()
        plot_king(globdata)
        plot_look(globdata)
        plot_eval(globdata)
        #if globdata.datatype == 'mock':
        #    plot_pntsdiff()

 
    if globdata.verbtype > 0:
        print 'Sampling...'
    
    timereal = zeros(globdata.numbproc)
    timeproc = zeros(globdata.numbproc)
    if globdata.numbproc == 1:
        gridchan = [work(globdata, 0)]
    else:
        if globdata.verbtype > 0:
            print 'Forking the sampler...'

        # process lock for simultaneous plotting
        globdata.lock = mp.Manager().Lock()

        # process pool
        pool = mp.Pool(globdata.numbproc)
        
        # spawn the processes
        workpart = functools.partial(work, globdata)
        indxproc = arange(globdata.numbproc)
        gridchan = pool.map(workpart, indxproc)

        pool.close()
        pool.join()

    for k in range(globdata.numbproc):
        timereal[k] = gridchan[k][18]
        timeproc[k] = gridchan[k][19]


    if globdata.verbtype > 0:
        print 'Accumulating samples from all processes...'
        tim0 = time.time()

    # parse the sample bundle
    listsampvarb = zeros((globdata.numbsamp, globdata.numbproc, globdata.maxmsampsize))
    listindxprop = zeros((globdata.numbswep, globdata.numbproc))
    listchrollik = zeros((globdata.numbswep, globdata.numbproc, globdata.numbchrollik))
    listchrototl = zeros((globdata.numbswep, globdata.numbproc, globdata.numbchrototl))
    listllik = zeros((globdata.numbsamp, globdata.numbproc))
    listlpri = zeros((globdata.numbsamp, globdata.numbproc))
    listaccp = zeros((globdata.numbswep, globdata.numbproc))
    listmodlcnts = zeros((globdata.numbsamp, globdata.numbproc, globdata.numbpixlsave))
    listpntsfluxmean = zeros((globdata.numbsamp, globdata.numbproc, globdata.numbener))
    listindxpntsfull = []
    listindxsampmodi = zeros((globdata.numbswep, globdata.numbproc), dtype=int)
    globdata.listauxipara = empty((globdata.numbswep, globdata.numbproc, globdata.numbcomp))
    globdata.listlaccfrac = empty((globdata.numbswep, globdata.numbproc))
    globdata.listnumbpair = empty((globdata.numbswep, globdata.numbproc))
    globdata.listjcbnfact = empty((globdata.numbswep, globdata.numbproc))
    globdata.listcombfact = empty((globdata.numbswep, globdata.numbproc))

    levi = 0.
    info = 0.
    
    for k in range(globdata.numbproc):
        globdata.rtag = retr_rtag(globdata, k)
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
        globdata.listauxipara[:, k, :] = listchan[9]
        globdata.listlaccfrac[:, k] = listchan[10]
        globdata.listnumbpair[:, k] = listchan[11]
        globdata.listjcbnfact[:, k] = listchan[12]
        globdata.listcombfact[:, k] = listchan[13]
        levi += listchan[14]
        info += listchan[15]
        listpntsfluxmean[:, k, :] = listchan[16]
        listchrollik[:, k, :] = listchan[17]

    listindxprop = listindxprop.flatten()
    globdata.listauxipara = globdata.listauxipara.reshape((globdata.numbswep * globdata.numbproc, globdata.numbcomp))
    globdata.listlaccfrac = globdata.listlaccfrac.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listnumbpair = globdata.listnumbpair.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listjcbnfact = globdata.listjcbnfact.reshape(globdata.numbswep * globdata.numbproc)
    globdata.listcombfact = globdata.listcombfact.reshape(globdata.numbswep * globdata.numbproc)
    
    listchrototl = listchrototl.reshape((globdata.numbproc * globdata.numbswep, globdata.numbchrototl)) 
    listchrollik = listchrollik.reshape((globdata.numbproc * globdata.numbswep, globdata.numbchrollik)) 
    listaccp = listaccp.flatten()
    listindxsampmodi = listindxsampmodi.flatten()
    
    globdata.rtag = retr_rtag(globdata, None)

    # collect posterior samples from the processes
    ## number of PS
    listnumbpnts = listsampvarb[:, :, globdata.indxsampnumbpnts].astype(int).reshape(globdata.numbsamp * globdata.numbproc, -1)
    ## FDF normalization
    listfdfnnorm = listsampvarb[:, :, globdata.indxsampfdfnnorm].reshape(globdata.numbsamp * globdata.numbproc, -1)
    ## FDF slope
    listfdfnslop = listsampvarb[:, :, globdata.indxsampfdfnslop].reshape(globdata.numbsamp * globdata.numbproc, globdata.numbpopl, globdata.numbener)
    ## PSF parameters
    listpsfipara = listsampvarb[:, :, globdata.indxsamppsfipara].reshape(globdata.numbsamp * globdata.numbproc, -1)
    ## Background normalization
    listnormback = listsampvarb[:, :, globdata.indxsampnormback].reshape(globdata.numbsamp * globdata.numbproc, globdata.numbback, globdata.numbener)
    ## PS parameters
    listlgal = [[] for l in globdata.indxpopl]
    listbgal = [[] for l in globdata.indxpopl]
    listspec = [[] for l in globdata.indxpopl]
    if globdata.colrprio:
        listsind = [[] for l in globdata.indxpopl]
    listspechist = empty((globdata.numbsamp * globdata.numbproc, globdata.numbpopl, globdata.numbener, globdata.numbspec))
    for k in range(globdata.numbproc):
        for j in range(globdata.numbsamp):            
            n = k * globdata.numbsamp + j
            indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(globdata, listindxpntsfull[k][j])
            for l in globdata.indxpopl:
                listlgal[l].append(listsampvarb[j, k, indxsamplgal[l]])
                listbgal[l].append(listsampvarb[j, k, indxsampbgal[l]])
                listspec[l].append(listsampvarb[j, k, indxsampspec[l]])
                if globdata.colrprio:
                    listsind[l].append(listsampvarb[j, k, indxsampsind[l]])
                for i in globdata.indxener:
                    listspechist[n, l, i, :] = histogram(listspec[l][n][i, :], globdata.binsspec[i, :])[0]
    for l in globdata.indxpopl:
        listlgaltemp = zeros((globdata.numbsamp, globdata.maxmnumbpnts[l])) - 1.
        listbgaltemp = zeros((globdata.numbsamp, globdata.maxmnumbpnts[l])) - 1.
        listspectemp = zeros((globdata.numbsamp, globdata.numbener, globdata.maxmnumbpnts[l])) - 1.
        listsindtemp = zeros((globdata.numbsamp, globdata.maxmnumbpnts[l])) - 1.
        for k in range(globdata.numbsamp):
            listlgaltemp[k, 0:listlgal[l][k].size] = listlgal[l][k]
            listbgaltemp[k, 0:listbgal[l][k].size] = listbgal[l][k]
            listspectemp[k, :, 0:listspec[l][k].shape[1]] = listspec[l][k]
            listsindtemp[k, 0:listsind[l][k].size] = listsind[l][k]
        listlgal[l] = listlgaltemp
        listbgal[l] = listbgaltemp 
        listspec[l] = listspectemp    
        listsind[l] = listsindtemp

    # auxiliary variables
    listpntsfluxmean = listpntsfluxmean.reshape(globdata.numbsamp * globdata.numbproc, globdata.numbener)
   
        
        
    if globdata.verbtype > 0:
        print 'Binning the probabilistic catalog...'
        tim0 = time.time()

    # posterior maps
    pntsprob = zeros((globdata.numbpopl, globdata.numbener, globdata.numbpixl, globdata.numbspec))
    for k in range(globdata.numbsamp):
        for l in globdata.indxpopl:
            for i in globdata.indxener:
                for h in range(globdata.numbspec):
                    indxpnts = where((globdata.binsspec[i, h] < listspec[l][k, i, :]) & (listspec[l][k, i, :] < globdata.binsspec[i, h+1]))[0]
                    hpixl = retr_indxpixl(globdata, listbgal[l][k, indxpnts], listlgal[l][k, indxpnts])
                    pntsprob[l, i, hpixl, h] += 1.
    
    
        
    if globdata.verbtype > 0:
        print 'Performing Gelman-Rubin convergence test...'
        tim0 = time.time()

    gmrbstat = zeros(globdata.numbpixlsave)
    for n in range(globdata.numbpixlsave):
        gmrbstat[n] = tdpy.mcmc.gmrb_test(listmodlcnts[:, :, n])


            
    pathprobcatl = os.environ["PCAT_DATA_PATH"] + '/probcatl_' + globdata.strgtime + '_' + globdata.rtag + '.fits'  
    
    head = pf.Header()
    head['numbener'] = (globdata.numbener, 'Number of energy bins')
    head['numbevtt'] = (globdata.numbevtt, 'Number of PSF class bins')
    head['numbpopl'] = (globdata.numbpopl, 'Number of PS population')
    head['numbpsfipara'] = globdata.numbpsfipara
    head['numbformpara'] = globdata.numbformpara
    
    head['numbsamp'] = globdata.numbsamp
    head['numbburn'] = globdata.numbburn
    head['numbswep'] = globdata.numbswep
    head['factthin'] = globdata.factthin
    head['numbpopl'] = globdata.numbpopl
    head['numbproc'] = globdata.numbproc
    
    head['maxmgang'] = globdata.maxmgang
    head['lgalcntr'] = globdata.lgalcntr
    head['bgalcntr'] = globdata.bgalcntr
    head['numbspec'] = globdata.numbspec
    if globdata.pixltype == 'heal':
        head['numbsideheal'] = globdata.numbsideheal
    if globdata.pixltype == 'cart':
        head['numbsidecart'] = globdata.numbsidecart
    
    head['minmlgal'] = globdata.minmlgal
    head['maxmlgal'] = globdata.maxmlgal
    head['minmbgal'] = globdata.minmbgal
    head['maxmbgal'] = globdata.maxmbgal
    
    head['datatype'] = globdata.datatype
    head['regitype'] = globdata.regitype
    head['psfntype'] = globdata.psfntype
    head['exprtype'] = globdata.exprtype
    head['pixltype'] = globdata.pixltype
    
    head['colrprio'] = globdata.colrprio
    head['trueinfo'] = globdata.trueinfo
    head['margsize'] = globdata.margsize
    head['strgtime'] = globdata.strgtime
    
    head['stdvfdfnnorm'] = globdata.stdvfdfnnorm
    head['stdvfdfnslop'] = globdata.stdvfdfnslop
    head['stdvpsfipara'] = globdata.stdvpsfipara
    head['stdvback'] = globdata.stdvback
    head['stdvlbhl'] = globdata.stdvlbhl
    head['stdvspec'] = globdata.stdvspec
    head['spmrlbhl'] = globdata.spmrlbhl
    head['fracrand'] = globdata.fracrand

    if globdata.trueinfo and globdata.datatype == 'mock':
        head['mockpsfntype'] = globdata.mockpsfntype

    head['strgexpr'] = globdata.strgexpr
    for k in globdata.indxback:
        head['strgback%04d' % k] = globdata.strgback[k]
        head['lablback%04d' % k] = globdata.lablback[k]
        head['nameback%04d' % k] = globdata.nameback[k]
    head['strgexpo'] = globdata.strgexpo
    
    head['levi'] = levi
    head['info'] = info
    
    listhdun = []
    listhdun.append(pf.PrimaryHDU(header=head))

    for l in globdata.indxpopl:
        listhdun.append(pf.ImageHDU(listlgal[l]))
        listhdun[-1].header['EXTNAME'] = 'lgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listbgal[l]))
        listhdun[-1].header['EXTNAME'] = 'bgalpopl%d' % l
        listhdun.append(pf.ImageHDU(listspec[l]))
        listhdun[-1].header['EXTNAME'] = 'specpopl%d' % l
        if globdata.colrprio:
            listhdun.append(pf.ImageHDU(listsind[l]))
            listhdun[-1].header['EXTNAME'] = 'sindpopl%d' % l

    
    listhdun.append(pf.ImageHDU(globdata.minmnormback))
    listhdun[-1].header['EXTNAME'] = 'minmnormback'

    listhdun.append(pf.ImageHDU(globdata.maxmnormback))
    listhdun[-1].header['EXTNAME'] = 'maxmnormback'


    listhdun.append(pf.ImageHDU(globdata.maxmnumbpnts))
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
    
    if globdata.trueinfo and globdata.datatype == 'mock':
        listhdun.append(pf.ImageHDU(globdata.mocknumbpnts))
        listhdun[-1].header['EXTNAME'] = 'mocknumbpnts'
        
        listhdun.append(pf.ImageHDU(globdata.mockfdfnslop))
        listhdun[-1].header['EXTNAME'] = 'mockfdfnslop'
        
        listhdun.append(pf.ImageHDU(globdata.mocknormback))
        listhdun[-1].header['EXTNAME'] = 'mocknormback'

    # convergence diagnostics
    listhdun.append(pf.ImageHDU(gmrbstat))
    listhdun[-1].header['EXTNAME'] = 'gmrbstat'
    listhdun.append(pf.ImageHDU(listmodlcnts))
    listhdun[-1].header['EXTNAME'] = 'modlcnts'
    
    
    listhdun.append(pf.ImageHDU(pntsprob))
    listhdun[-1].header['EXTNAME'] = 'pntsprob'
   
    # truth information
    if globdata.trueinfo:
        listhdun.append(pf.ImageHDU(globdata.truenumbpnts))
        listhdun[-1].header['EXTNAME'] = 'truenumbpnts'

        for l in globdata.indxpopl:
            listhdun.append(pf.ImageHDU(globdata.truelgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truelgalpopl%d' % l

            listhdun.append(pf.ImageHDU(globdata.truebgal[l]))
            listhdun[-1].header['EXTNAME'] = 'truebgalpopl%d' % l

            listhdun.append(pf.ImageHDU(globdata.truespec[l]))
            listhdun[-1].header['EXTNAME'] = 'truespecpopl%d' % l

            if globdata.colrprio:
                listhdun.append(pf.ImageHDU(globdata.truesind[l]))
                listhdun[-1].header['EXTNAME'] = 'truesindpopl%d' % l

        listhdun.append(pf.ImageHDU(globdata.truefdfnslop))
        listhdun[-1].header['EXTNAME'] = 'truefdfnslop'

        listhdun.append(pf.ImageHDU(globdata.truenormback))
        listhdun[-1].header['EXTNAME'] = 'truenormback'

        listhdun.append(pf.ImageHDU(globdata.truepsfipara))
        listhdun[-1].header['EXTNAME'] = 'truepsfipara'

        listhdun.append(pf.ImageHDU(globdata.indxtruepntstimevari))
        listhdun[-1].header['EXTNAME'] = 'indxtruepntstimevari'

    if globdata.colrprio:
        globdata.minmspec = globdata.minmspec[globdata.indxenerfdfn]
        globdata.maxmspec = globdata.maxmspec[globdata.indxenerfdfn]

    # hyperprior limits
    listhdun.append(pf.ImageHDU(globdata.minmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'minmfdfnnorm'
    listhdun.append(pf.ImageHDU(globdata.maxmfdfnnorm))
    listhdun[-1].header['EXTNAME'] = 'maxmfdfnnorm'
    listhdun.append(pf.ImageHDU(globdata.minmfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'minmfdfnslop'
    listhdun.append(pf.ImageHDU(globdata.maxmfdfnslop))
    listhdun[-1].header['EXTNAME'] = 'maxmfdfnslop'

    listhdun.append(pf.ImageHDU(globdata.minmspec))
    listhdun[-1].header['EXTNAME'] = 'minmspec'
    listhdun.append(pf.ImageHDU(globdata.maxmspec))
    listhdun[-1].header['EXTNAME'] = 'maxmspec'
    if globdata.colrprio:
        listhdun.append(pf.ImageHDU(globdata.minmsind))
        listhdun[-1].header['EXTNAME'] = 'minmsind'
        listhdun.append(pf.ImageHDU(globdata.maxmsind))
        listhdun[-1].header['EXTNAME'] = 'maxmsind'
        listhdun.append(pf.ImageHDU(globdata.meansind))
        listhdun[-1].header['EXTNAME'] = 'meansind'
        listhdun.append(pf.ImageHDU(globdata.stdvsind))
        listhdun[-1].header['EXTNAME'] = 'stdvsind'
    
    listhdun.append(pf.ImageHDU(globdata.binsener))
    listhdun[-1].header['EXTNAME'] = 'binsener'
    listhdun.append(pf.ImageHDU(globdata.indxenerfdfn))
    listhdun[-1].header['EXTNAME'] = 'indxenerfdfn'
    
    listhdun.append(pf.ImageHDU(globdata.indxenerincl))
    listhdun[-1].header['EXTNAME'] = 'indxenerincl'
    
    listhdun.append(pf.ImageHDU(globdata.indxevttincl))
    listhdun[-1].header['EXTNAME'] = 'indxevttincl'
 
    # utilities
    listhdun.append(pf.ImageHDU(listpntsfluxmean))
    listhdun[-1].header['EXTNAME'] = 'listpntsfluxmean'
    
    listhdun.append(pf.ImageHDU(globdata.probprop))
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
    
    listhdun.append(pf.ImageHDU(globdata.listauxipara))
    listhdun[-1].header['EXTNAME'] = 'auxipara'
    
    listhdun.append(pf.ImageHDU(globdata.listlaccfrac))
    listhdun[-1].header['EXTNAME'] = 'laccfrac'
    
    listhdun.append(pf.ImageHDU(globdata.listnumbpair))
    listhdun[-1].header['EXTNAME'] = 'numbpair'
    
    listhdun.append(pf.ImageHDU(globdata.listjcbnfact))
    listhdun[-1].header['EXTNAME'] = 'jcbnfact'
    
    listhdun.append(pf.ImageHDU(globdata.listcombfact))
    listhdun[-1].header['EXTNAME'] = 'combfact'
    
    pf.HDUList(listhdun).writeto(pathprobcatl, clobber=True, output_verify='ignore')

    if globdata.makeplot:
        plot_post(pathprobcatl)

    timetotlreal = time.time() - timetotlreal
    timetotlproc = time.clock() - timetotlproc
     
    if globdata.verbtype > 0:
        for k in range(globdata.numbproc):
            print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, timereal[k], timeproc[k])
        print 'PCAT has run in %d real seconds, %d CPU seconds.' % (timetotlreal, sum(timeproc))
        print 'The ensemble of catalogs is at ' + pathprobcatl
        if globdata.makeplot:
            print 'The plots are in ' + globdata.plotpath
        
    return gridchan
    
    
def plot_samp(globdata):

    globdata.thisresicnts = globdata.datacnts - globdata.thismodlcnts

    globdata.thispsfn = retr_psfn(globdata, globdata.thissampvarb[globdata.indxsamppsfipara], globdata.indxener, globdata.angldisp, globdata.psfntype)
    if globdata.pixltype == 'cart':
        globdata.thispsfn = globdata.thispsfn.reshape((globdata.numbener, -1, globdata.numbevtt))
    thisfwhm = 2. * retr_psfnwdth(globdata, globdata.thispsfn, 0.5)

    plot_psfn(globdata)

    globdata.thisbackcntsmean = empty((globdata.numbener, globdata.numbevtt))
    for c in globdata.indxback:
        globdata.thisbackcntsmean += mean(globdata.thissampvarb[globdata.indxsampnormback[c, :, None, None]] * globdata.backflux[c] * globdata.expo * \
            globdata.diffener[:, None, None] * pi * thisfwhm[:, None, :]**2 / 4., 1)

    thiscnts = []
    for l in globdata.indxpopl:
        indxpixltemp = retr_indxpixl(globdata, globdata.thissampvarb[globdata.thisindxsampbgal[l]], globdata.thissampvarb[globdata.thisindxsamplgal[l]])
        cntstemp = globdata.thissampvarb[globdata.thisindxsampspec[l]][:, :, None] *             globdata.expo[:, indxpixltemp, :] * globdata.diffener[:, None, None]
        thiscnts.append(cntstemp)

        #if globdata.thissampvarb[globdata.indxsampnumbpnts[l]] > 1:
        if globdata.colrprio:
            plot_histsind(globdata, l)
        plot_scatpixl(globdata, l)
       
        if globdata.trueinfo:
            indxmodl, globdata.trueindxpntsbias, globdata.trueindxpntsmiss = pair_catl(globdata, l, globdata.thissampvarb[globdata.thisindxsamplgal[l]], \
                globdata.thissampvarb[globdata.thisindxsampbgal[l]], globdata.thissampvarb[globdata.thisindxsampspec[l]])

            thisspecmtch = globdata.thissampvarb[globdata.thisindxsampspec[l]][:, indxmodl]
            thisspecmtch[:, globdata.trueindxpntsmiss] = 0.
            if globdata.verbtype > 1:
                print 'thisspecmtch: (popl%d) ' % l
                print thisspecmtch
            plot_scatspec(globdata, l, thisspecmtch=thisspecmtch)
        plot_histspec(globdata, l)
        plot_histcnts(globdata, l, thiscnts)
        if globdata.numbback == 2:
            plot_compfrac(globdata)

    for i in globdata.indxener:
        
        plot_datacnts(globdata, i, None)
        #plot_catl(globdata, i, None, thiscnts)
        plot_modlcnts(globdata, i, None)
        plot_resicnts(globdata, i, None, globdata.thisresicnts)

        #for m in indxevtt:
        #    plot_datacnts(globdata, i, m)
        #    plot_catl(globdata, i, m, thiscnts)
        #    plot_modlcnts(globdata, i, m)
        #    plot_resicnts(globdata, i, m, globdata.thisresicnts)
        
    #if globdata.numbener == 3:
    #    plot_datacnts(globdata, None, None)
        
    #plot_fwhm(globdata, thisfwhm)
    #plot_backcntsmean(globdata, globdata.thisbackcntsmean)
    
    tempsampvarb, tempppixl, tempcnts, temppntsflux, tempmodlflux, tempmodlcnts = pars_samp(globdata, globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])
    globdata.errrmodlcnts = 100. * (globdata.thismodlcnts - tempmodlcnts) / tempmodlcnts
    
    for i in globdata.indxener:
        plot_errrcnts(globdata, i, None, globdata.errrmodlcnts)

    if amax(abs(globdata.errrmodlcnts)) > 0.1 and False:
        print 'Approximation error went above the limit!'
    
    
def rjmc(globdata, indxprocwork):

    # sweeps to be saved
    boolsave = zeros(globdata.numbswep, dtype=bool)
    indxswepsave = arange(globdata.numbburn, globdata.numbswep, globdata.factthin)
    boolsave[indxswepsave] = True
    
    sampindx = zeros(globdata.numbswep, dtype=int)
    sampindx[indxswepsave] = arange(globdata.numbsamp)

    listsampvarb = zeros((globdata.numbsamp, globdata.maxmsampsize)) + -1.
    listindxprop = zeros(globdata.numbswep)
    listchrototl = zeros((globdata.numbswep, globdata.numbchrototl))
    listchrollik = zeros((globdata.numbswep, globdata.numbchrollik))
    listllik = zeros(globdata.numbsamp)
    listlprising = zeros(globdata.numbsamp)
    listlpri = zeros(globdata.numbsamp)
    listaccp = zeros(globdata.numbswep, dtype=bool)
    listaccpspec = []
    listindxsampmodi = zeros(globdata.numbswep, dtype=int)
    listmodlcnts = zeros((globdata.numbsamp, globdata.numbpixlsave))
    listpntsfluxmean = zeros((globdata.numbsamp, globdata.numbener))
    listindxpntsfull = []
    
    globdata.listauxipara = zeros((globdata.numbswep, globdata.numbcomp))
    globdata.listlaccfrac = zeros(globdata.numbswep)
    globdata.listnumbpair = zeros(globdata.numbswep)
    globdata.listjcbnfact = zeros(globdata.numbswep)
    globdata.listcombfact = zeros(globdata.numbswep)

    globdata.cntrswep = 0
    
    # initialize the chain
    retr_llik(globdata, listchrollik, init=True)
    retr_lpri(globdata, init=True)

    # current sample index
    thiscntr = -1
    
    while globdata.cntrswep < globdata.numbswep:
        
        timeinit = time.time()
        
        if globdata.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % globdata.cntrswep

        thismakefram = (globdata.cntrswep % globdata.plotperd == 0) and indxprocwork == int(float(globdata.cntrswep) / globdata.numbswep * globdata.numbproc) and globdata.makeplot
        globdata.reje = False
    
        # choose a proposal type
        retr_thisindxprop(globdata, globdata.drmcsamp[:, 0])
            
        # save the proposal type
        listindxprop[globdata.cntrswep] = globdata.thisindxprop
        if globdata.verbtype > 1:
            print 'thisindxprop: ', globdata.strgprop[globdata.thisindxprop]
        
        if globdata.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print

        # propose the next sample
        timebegn = time.time()
        retr_prop(globdata)
        timefinl = time.time()
        listchrototl[globdata.cntrswep, 1] = timefinl - timebegn

        # plot the current sample
        if thismakefram:
            print 'Process %d is in queue for making a frame.' % indxprocwork
            if globdata.numbproc > 1:
                globdata.lock.acquire()
            print 'Process %d started making a frame.' % indxprocwork
            plot_samp(globdata)
            print 'Process %d finished making a frame.' % indxprocwork
            if globdata.numbproc > 1:
                globdata.lock.release()
            
        # reject the sample if proposal is outside the prior
        if globdata.thisindxprop == globdata.indxpropfdfnslop:
            if globdata.drmcsamp[globdata.indxsampfdfnslopmodi, 1] < 0. or globdata.drmcsamp[globdata.indxsampfdfnslopmodi, 1] > 1.:
                globdata.reje = True
        elif globdata.thisindxprop != globdata.indxpropbrth and globdata.thisindxprop != globdata.indxpropdeth and not globdata.reje:
            if where((globdata.drmcsamp[globdata.indxsampmodi, 1] < 0.) | (globdata.drmcsamp[globdata.indxsampmodi, 1] > 1.))[0].size > 0:
                globdata.reje = True

        if globdata.thisindxprop == globdata.indxproppsfipara:
            if globdata.psfntype == 'doubking':
                if globdata.nextsampvarb[globdata.indxsamppsfipara[1]] >= globdata.nextsampvarb[globdata.indxsamppsfipara[3]]:
                    globdata.reje = True
            elif globdata.psfntype == 'doubgaus':
                if globdata.nextsampvarb[globdata.indxsamppsfipara[1]] >= globdata.nextsampvarb[globdata.indxsamppsfipara[2]]:
                    globdata.reje = True
            
        if not globdata.reje:

            # evaluate the log-prior
            timebegn = time.time()
            retr_lpri(globdata)
            timefinl = time.time()
            listchrototl[globdata.cntrswep, 2] = timefinl - timebegn

            # evaluate the log-likelihood
            timebegn = time.time()
            retr_llik(globdata, listchrollik) 
            timefinl = time.time()
            listchrototl[globdata.cntrswep, 3] = timefinl - timebegn
            
            # evaluate the acceptance probability
            accpprob = exp(globdata.deltllik + globdata.deltlpri + globdata.laccfrac)

            if globdata.verbtype > 1:
                print 'deltlpri'
                print globdata.deltlpri
                print 'deltllik'
                print globdata.deltllik
                print 'laccfrac'
                print globdata.laccfrac
                print
                
        else:
            accpprob = 0.
    
    
        # accept the sample
        if accpprob >= rand():

            if globdata.verbtype > 1:
                print 'Accepted.'

            # update the current state
            updt_samp(globdata)

            listaccp[globdata.cntrswep] = True

        # reject the sample
        else:

            if globdata.verbtype > 1:
                print 'Rejected.'

            listaccp[globdata.cntrswep] = False
             
        # sanity checks
        indxsampbadd = where((globdata.drmcsamp[1:, 0] > 1.) | (globdata.drmcsamp[1:, 0] < 0.))[0] + 1
        if indxsampbadd.size > 0:
            print 'Unit sample vector went outside [0,1]!'
            print 'cntrswep'
            print globdata.cntrswep
            print 'thisindxprop'
            print globdata.thisindxprop
            print 'indxsampbadd'
            print indxsampbadd
            print 'drmcsamp'
            print globdata.drmcsamp[indxsampbadd, :]
            return
            
        for l in globdata.indxpopl:
            for i in globdata.indxenerfdfn:
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] < globdata.minmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went below the prior range!'
                if where(globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]] > globdata.maxmspec[i])[0].size > 0:
                    print 'Spectrum of some PS went above the prior range!'          

        # temp
        if False:
            tempsampvarb, tempppixl, tempcnts, temppntsflux, tempmodlflux, tempmodlcnts = pars_samp(globdata, globdata.thisindxpntsfull, globdata.drmcsamp[:, 0])
            globdata.errrmodlcnts = 100. * (globdata.thismodlcnts - tempmodlcnts) / tempmodlcnts
            print 'globdata.errrmodlcnts'
            print amin(globdata.errrmodlcnts), amax(globdata.errrmodlcnts)
            if amax(globdata.errrmodlcnts) > 0.1 or amin(globdata.errrmodlcnts) < -0.1:
                return

        # save the sample
        if boolsave[globdata.cntrswep]:
            listsampvarb[sampindx[globdata.cntrswep], :] = globdata.thissampvarb
            listmodlcnts[sampindx[globdata.cntrswep], :] = globdata.thismodlcnts[0, globdata.indxpixlsave, 0]
            listpntsfluxmean[sampindx[globdata.cntrswep], :] = mean(sum(globdata.thispntsflux * globdata.expo, 2) / sum(globdata.expo, 2), 1)
            listindxpntsfull.append(globdata.thisindxpntsfull)
            listllik[sampindx[globdata.cntrswep]] = sum(globdata.thisllik)
            
            lpri = 0.
            for l in globdata.indxpopl:
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[l]]
                fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[l]]
                lpri += numbpnts * globdata.priofactlgalbgal + globdata.priofactfdfnslop + globdata.priofactfdfnnorm - log(fdfnnorm)
                for i in globdata.indxenerfdfn:
                    flux = globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]]
                    fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[l, i]]
                    lpri -= log(1. + fdfnslop**2)
                    lpri += sum(log(pdfn_spec(globdata, flux, fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])))
            listlpri[sampindx[globdata.cntrswep]] = lpri
            
            if globdata.tracsamp:
                
                numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                listdiffllikdiffpara.append(diffllikdiffpara)

                tranmatr = diffllikdiffpara[:, None] * listdiffllikdiffpara[globdata.cntrswep-1][None, :]
                listtranmatr.append(tranmatr)

        # save the execution time for the sweep
        if not thismakefram:
            tim1 = time.time()
            listchrototl[globdata.cntrswep, 0] = tim1 - timeinit

        # log the progress
        if globdata.verbtype > 0:
            thiscntr = tdpy.util.show_prog(globdata.cntrswep, globdata.numbswep, thiscntr, indxprocwork=indxprocwork)
            
        if globdata.diagsamp:
            plot_datacnts(0, 0, nextstat=True)
            plot_resicnts(0, 0, thisresicnts, nextstat=True)
        
        if globdata.verbtype > 1:
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
        globdata.cntrswep += 1

    
    if globdata.verbtype > 1:
        print 'listsampvarb: '
        print listsampvarb
    
    listchrollik = array(listchrollik)

    # correct the likelihoods for the constant data dependent factorial
    listllik -= sum(sp.special.gammaln(globdata.datacnts + 1))
    
    # calculate the log-evidence and relative entropy using the harmonic mean estimator
    minmlistllik = amin(listllik)
    levi = -log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    info = mean(listllik) - levi

    listchan = [listsampvarb, listindxprop, listchrototl, listllik, listlpri, listaccp, listmodlcnts, listindxpntsfull, listindxsampmodi, \
        globdata.listauxipara, globdata.listlaccfrac, globdata.listnumbpair, globdata.listjcbnfact, globdata.listcombfact, \
        levi, info, listpntsfluxmean, listchrollik]
    
    return listchan

