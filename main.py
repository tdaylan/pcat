# common imports
from __init__ import *

# internal functions
from util import *
from visu import *

def init( \
         # user interaction
         verbtype=1, \
         pathbase=os.environ["PCAT_DATA_PATH"], \
         showmoreaccp=False, \

         # diagnostics
         diagmode=True, \
         emptsamp=False, \

         evoltype='samp', \
        
         killexpo=False, \
         highexpo=False, \
         sqzeprop=False, \
         explprop=False, \
        
         asscrefr=True, \

         # chain setup
         numbswep=2000000, \
         numbsamp=None, \
         numbburn=None, \
         factthin=None, \
        
         # output
         ## condensed catalog
         condcatl=False, \
        
         refrlegd=None, \
         refrlegdpopl=None, \
         fittlegdpopl=None, \

         randseedelem=False, \
         indxevttincl=None, \
         indxenerincl=None, \
        
         listmask=None, \

         listnamefeatsele=None, \
         burntmpr=False, \
         #burntmpr=True, \
        
         numbregi=1, \
         # evaluate the likelihood inside circles around elements
         elemspatevaltype=None, \
        
         namestattrue=None, \

         shrtfram=False, \
        
         suprelem=True, \
         plotelemcorr=True, \
         relnprio=False, \
        
         # elements
         ## vary projected scale radius
         variasca=True, \
         ## vary projected cutoff radius
         variacut=True, \
        
         # name of the configuration
         strgcnfg=None, \
        
         allwrefr=True, \

         # numpy RNG seed
         seedstat=None, \
    
         # metamodel settings
         ## number of spatial dimensions
         numbspatdims=2, \
         
         #hostemistype=None, \

         ## lens model
         #lensmodltype=None, \
        
         penalpridiff=False, \

         ## PSF evaluation type
         trueoaxitype=False, \
         fittoaxitype=False, \
         ## kernel evaluation type
         kernevaltype='ulip', \

         ## lens specific
         ### Sersic type
         serstype='vauc', \

         # initialization
         ## initialization type
         inittype=None, \
        
         loadvaripara=False, \
         
         # save the state of the MCMC
         savestat=False, \
         namesavestat=None, \
         # recover the state from a previous run
         namerecostat=None, \
         forcsavestat=False, \

         # proposals
         propcova=True, \
         propwithsing=True, \
         # Hessian estimation
         optitype='hess', \
         regulevi=False, \
         
         # modes of operation
         ## plot results from a previous chain
         rtagredo=None, \
         ## interactive
         intrevalcntpmodl=False, \
         intrevalcntpresi=False, \
         ## only generate and plot mock data
         mockonly=False, \
         ## perform an additional run sampling from the prior
         checprio=False, \

         strgexprsbrt=None, \
         anglassc=None, \
         nameexpr=None, \
         
         lgalprio=None, \
         bgalprio=None, \
         minmcntpdata=None, \
         strgexpo=None, \
         
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         anlytype=None, \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         pixltype=None, \
         forccart=False, \
         fittampldisttype=None, \
        
         allwfixdtrue=True, \
         asscmetrtype='dist', \

         # plotting
         numbswepplot=None, \
         
         makeplot=True, \
         makeplotinit=True, \
         makeplotfram=True, \
         makeplotpost=True, \
         
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
         priofactdoff=1., \
         
         # lensing
         fittrelnpowr=0., \

         # temp
         margfactmodl=1., \
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
         
         propmeanelem=True, \
         propdist=True, \
         proppsfp=True, \
         propbacp=True, \
         proplenp=None, \
         propcomp=None, \
         probtran=None, \
         probspmr=0.3, \
         # when proposing from the covariance, fracproprand should be very small!
         fracproprand=0., \
            
         jitt=False, \

         radispmr=None, \

         numbsidecart=None, \
         numbsideheal=256, \
         numbdatasamp=100, \

         defa=False, \
         **args \
        ):

    # construct the global object 
    gdat = tdpy.util.gdatstrt()
    for attr, valu in locals().iteritems():
        if '__' not in attr and attr != 'gdat':
            setattr(gdat, attr, valu)
    
    # copy all provided inputs to the global object
    for strg, valu in args.iteritems():
        setattr(gdat, strg, valu)

    # defaults
    if gdat.strgexprsbrt == None:
        gdat.datatype = 'mock'
    else:
        gdat.datatype = 'inpt'
    
    # list of models
    gdat.liststrgmodl = ['fitt']
    gdat.listlegdmodl = ['Fitting']
    if gdat.datatype == 'mock':
        gdat.liststrgmodl += ['true']
        gdat.listlegdmodl += ['True']
    
    # PCAT folders
    if gdat.pathbase[-1] != '/':
        gdat.pathbase += '/'
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathdataopti = gdat.pathdata + 'opti/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    gdat.pathoutp = gdat.pathdata + 'outp/'
    
    # run tag
    gdat.strgswep = '%d' % (gdat.numbswep)
    
    # preliminary setup
    ## time stamp
    gdat.strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    ## name of the configuration function
    if gdat.strgcnfg == None:
        gdat.strgcnfg = inspect.stack()[1][3]
   
    ## number of burned sweeps
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 10
    
    ## number of processes
    gdat.strgproc = os.uname()[1]
    if gdat.numbproc == None:
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu' or gdat.strgproc == 'wise':
            gdat.numbproc = 10
        else:
            gdat.numbproc = 1
    
    # string describing the number of sweeps
    gdat.strgnumbswep = '%d' % gdat.numbswep
    
    # output paths
    gdat.rtag = gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.strgnumbswep
    gdat.pathoutpthis = gdat.pathoutp + gdat.rtag + '/'

    # plot previous run tag
    if gdat.rtagredo != None:
        gdat.pathoutpredo = gdat.pathoutp + gdat.rtagredo + '/'
        path = gdat.pathoutpredo + 'gdatfinl'
        if gdat.verbtype > 0:
            print 'PCAT started in retro-mode.'
            print 'Reading %s...' % path
        gdat = readfile(path) 
        plot_post(gdat)
        return
    
    gdat.strgvers = 'v0.3'
    ## catalog output
    # create output folder for the run
    os.system('mkdir -p %s' % gdat.pathoutpthis)

    ## factor by which to thin the sweeps to get samples
    
    if gdat.factthin != None and gdat.numbsamp != None:
        raise Exception('Both factthin and numbsamp cannot be provided at the same time.')
    elif gdat.factthin == None and gdat.numbsamp == None:
        gdat.factthin = int(ceil(5e-4 * (gdat.numbswep - gdat.numbburn) * gdat.numbproc))
        gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    elif gdat.numbsamp != None:
        gdat.factthin = int((gdat.numbswep - gdat.numbburn) * gdat.numbproc / gdat.numbsamp)
    elif gdat.factthin != None:
        gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    if not isinstance(gdat.numbsamp, int) or not isinstance(gdat.factthin, int) or not isinstance(gdat.numbburn, int) or not isinstance(gdat.numbswep, int):
        raise Exception('Number of samples is not an integer.')

    # samples to be saved
    gdat.indxsamp = arange(gdat.numbsamp)
    
    # samples to be saved from all chains
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc
    
    if gdat.verbtype > 0:
        print 'PCAT %s started at %s' % (gdat.strgvers, gdat.strgtimestmp)
        print 'Configuration %s' % gdat.strgcnfg
        print 'Initializing...'
        print '%d samples will be taken, discarding the first %d. The chain will be thinned by a factor of %d.' % (gdat.numbswep, gdat.numbburn, gdat.factthin)
        print 'The resulting chain will contain %d samples per chain and %d samples in total.' % (gdat.numbsamp, gdat.numbsamptotl)
    
    if gdat.anlytype == None:
        if gdat.exprtype == 'chan':
            gdat.anlytype = 'home'
        elif gdat.exprtype == 'ferm':
            gdat.anlytype = 'manu'
        else:
            gdat.anlytype = 'nomi'
    
    if gdat.exprtype == 'ferm':
        elemtype = ['lghtpnts']
    if gdat.exprtype == 'chan':
        elemtype = ['lghtpnts']
    if gdat.exprtype == 'hubb':
        elemtype = ['lghtpnts', 'lens', 'lghtgausbgrd']
    if gdat.exprtype == 'sdyn':
        elemtype = ['clus']
    setp_varbvalu(gdat, 'elemtype', elemtype)
    
    # feature correlated with the significance of elements
    gdat.namefeatsign = 'deltllik'
    if gdat.datatype == 'mock':
        gdat.namefeatsignrefr = 'deltllik'
    
    if gdat.numbsidecart == None:
        if gdat.datatype == 'mock':
            gdat.numbsidecart = 100

    if gdat.inittype == None:
        if gdat.datatype == 'inpt':
            gdat.inittype = 'rand'
        else:
            gdat.inittype = 'pert'

    gdat.commelemtype = []
    for strgmodl in gdat.liststrgmodl:
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        for elemtypetemp in elemtype:
            if not elemtypetemp in gdat.commelemtype:
                gdat.commelemtype.append(elemtypetemp)
            
    #for strgvarb in ['boolelempsfnanyy']:
    #    varbcomm = [[] for l in indxpopl]
    #    for strgmodl in gdat.liststrgmodl:
    #        varb = getattr(gdat, strgmodl + strgvarb)
    #        varbcomm = varb
    #        for elemtypetemp in :
    #            if not elemtypetemp in gdat.commelemtype:
    #                gdat.commelemtype.append(elemtypetemp)
    
    if 'lens' in gdat.commelemtype:
        gdat.hubbexpofact = 1.63050e-19
    
    gdat.sdenunit = 'degr'

    if gdat.strgexpo == None:
        if gdat.exprtype == 'ferm':
            gdat.strgexpo = 'expofermcmp0igal.fits'
        elif gdat.exprtype == 'hubb':
            gdat.strgexpo = 1000. / gdat.hubbexpofact
        else:
            gdat.strgexpo = 1.
    
    if gdat.indxevttfull == None:
        if gdat.exprtype == 'ferm':
            if gdat.anlytype.startswith('rec7') or gdat.anlytype.startswith('manu'):
                gdat.indxevttfull = arange(2)
            if gdat.anlytype.startswith('rec8'):
                gdat.indxevttfull = arange(4)
        else:
            gdat.indxevttfull = arange(1)
    
    if gdat.indxevttincl == None:
        if gdat.exprtype == 'ferm':
            if gdat.anlytype.startswith('rec7') or gdat.anlytype.startswith('manu'):
                gdat.indxevttincl = array([0, 1])
            if gdat.anlytype.startswith('rec8'):
                gdat.indxevttincl = array([2, 3])
        else:
            gdat.indxevttincl = arange(1)
    
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

    if gdat.exprtype == 'ferm':
        gdat.lablenerunit = 'GeV'
    if gdat.exprtype == 'chan':
        gdat.lablenerunit = 'keV'
    if gdat.exprtype == 'sdyn':
        gdat.lablenerunit = ''
    
    gdat.factergskevv = 1.6e-9
    if gdat.exprtype == 'ferm':
        gdat.listspecconvunit = [['ene2', 'gevv']]
    if gdat.exprtype == 'chan':
        gdat.listspecconvunit = [['ene0', 'kevv'], ['ene2', 'kevv'], ['ene2', 'ergs'], ['ene3', 'ergs', '0520', 0.5,  2.], \
                                                                                       ['ene3', 'ergs', '0210',  2., 10.], \
                                                                                       ['ene3', 'ergs', '0510', 0.5, 10.], \
                                                                                       ['ene3', 'ergs', '0208',  2.,  8.], \
                                                                                       ['ene3', 'ergs', '0508', 0.5,  8.], \
                                                                                       ['ene3', 'ergs', '0207',  2.,  7.], \
                                                                                       ['ene3', 'ergs', '0507', 0.5,  7.]]
    if gdat.exprtype == 'hubb':
        gdat.listspecconvunit = [['ene3', 'ergs']]
    
    if gdat.pixltype == None:
        if gdat.exprtype == 'ferm':
            gdat.pixltype = 'heal'
        else:
            gdat.pixltype = 'cart'
    
    # temp
    #if gdat.exprtype == 'chan' and (gdat.anlytype.startswith('home') or gdat.anlytype.startswith('extr')):
    #    gdat.truelegdpopl = ['AGN', 'Galaxy']

    ## generative model
    if gdat.anglfact == None:
        if gdat.exprtype == 'ferm':
            gdat.anglfact = 180. / pi
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan' or gdat.exprtype == 'hubb':
            gdat.anglfact = 3600 * 180. / pi
        if gdat.exprtype == 'sche' or gdat.exprtype == 'sdyn':
            gdat.anglfact = 1.

    ### spatial extent of the data
    if gdat.maxmgangdata == None:
        if gdat.exprtype == 'chan':
            gdat.maxmgangdata = 0.492 / gdat.anglfact * gdat.numbsidecart / 2.
        if gdat.exprtype == 'ferm':
            gdat.maxmgangdata = 15. / gdat.anglfact
        if gdat.exprtype == 'sdyn':
            gdat.maxmgangdata = 1.
        if gdat.exprtype == 'hubb':
            gdat.maxmgangdata = 2. / gdat.anglfact
    
    ## experiment defaults
    if gdat.binsenerfull == None:
        if gdat.exprtype == 'ferm':
            gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
        if gdat.exprtype == 'chan':
            if gdat.anlytype.startswith('home'):
                gdat.binsenerfull = array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
            if gdat.anlytype.startswith('extr'):
                gdat.binsenerfull = array([0.5, 2., 8.])
            if gdat.anlytype == 'spec':
                gdat.binsenerfull = logspace(log10(0.5), log10(10.), 21)
        if gdat.exprtype == 'hubb':
            # temp
            gdat.binsenerfull = array([500., 750, 1000.])
    
    if gdat.binsenerfull == None:
        gdat.enerbins = False
    else:
        gdat.enerbins = True
    
    # temp
    gdat.enerbinsadje = True

    if gdat.enerbins:
        if gdat.enerbinsadje:
            gdat.meanenerfull = sqrt(gdat.binsenerfull[1:] * gdat.binsenerfull[:-1])
    
    ### background
    #### template
    if gdat.exprtype == 'ferm':
        if gdat.anlytype == 'bfun':
            gdat.ordrexpa = 10
            gdat.numbexpasing = gdat.ordrexpa**2
            gdat.numbexpa = gdat.numbexpasing * 4
            gdat.indxexpa = arange(gdat.numbexpa)
            backtype = ['bfun%04d' % k for k in gdat.indxexpa]
        else:
            backtype = [1., 'sbrtfdfmsmthrec7norm.fits']
    if gdat.exprtype == 'chan':
        
        # particle background
        if gdat.anlytype == 'spec':
            # temp -- this is fake!
            sbrtparttemp = array([70.04, 70.04, 12.12, 15.98, 10.79, 73.59, 73.59])
            binsenerpart = logspace(log10(0.5), log10(10.), 6)
            meanenerpart = sqrt(binsenerpart[:-1] * binsenerpart[1:])
            meanenerparttemp = concatenate((array([0.5]), meanenerpart, array([10.])))
            backtypetemp = interp(gdat.meanenerfull, meanenerparttemp, sbrtparttemp)
        if gdat.anlytype.startswith('home') :
            backtypetemp = array([70.04, 12.12, 15.98, 10.79, 73.59]) / 70.04
        if gdat.anlytype.startswith('extr'):
            backtypetemp = 'sbrtchanback' + gdat.anlytype + '.fits'
        if gdat.anlytype == 'spec':
            backtype = [[1e2, 2.], backtypetemp]
        else:
            backtype = [1., backtypetemp]
    if gdat.exprtype == 'hubb':
        backtype = [1.]
    if gdat.exprtype == 'sdyn':
        backtype = [1.]
    setp_varbvalu(gdat, 'backtype', backtype)
        
    setpprem(gdat)
    
    if gdat.exprtype == 'chan':
        if gdat.anlytype.startswith('extr'):
            meanbacpbac1 = 1.
        else:
            meanbacpbac1 = 70.04
        setp_varbvalu(gdat, 'scalbacp', 'gaus', back=1, regi='full')
        stdvbacpbac1 = 1e-5 * meanbacpbac1
        setp_varblimt(gdat, 'bacp', [meanbacpbac1, stdvbacpbac1], back=1, regi='full', typelimt='meanstdv')

    if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan':
        gdat.enerdiff = True
    if gdat.exprtype == 'hubb' or gdat.exprtype == 'sdyn':
        gdat.enerdiff = False
    
    if gdat.indxenerincl == None:
        if gdat.binsenerfull != None:
            gdat.indxenerincl = arange(gdat.binsenerfull.size - 1)
        if gdat.exprtype == 'ferm':
            gdat.indxenerincl = arange(1, 4)
        if gdat.exprtype == 'chan':
            if gdat.anlytype.startswith('home'):
                gdat.indxenerincl = arange(5)
            if gdat.anlytype.startswith('extr'):
                gdat.indxenerincl = arange(2)
        if gdat.exprtype == 'hubb':
            gdat.indxenerincl = array([1])
            #gdat.indxenerincl = array([0, 1])
        if gdat.exprtype == 'sdyn':
            gdat.indxenerincl = array([0])
    
    # energy band string
    if gdat.strgenerfull == None:
        if gdat.exprtype == 'sdss':
            gdat.strgenerfull = ['z-band', 'i-band', 'r-band', 'g-band', 'u-band']
        if gdat.exprtype == 'hubb':
            gdat.strgenerfull = ['F606W', 'F814W']
        if gdat.exprtype == 'ferm' or gdat.exprtype == 'chan': 
            gdat.strgenerfull = []
            for i in range(len(gdat.binsenerfull) - 1):
                gdat.strgenerfull.append('%.3g %s - %.3g %s' % (gdat.binsenerfull[i], gdat.lablenerunit, gdat.binsenerfull[i+1], gdat.lablenerunit))
        if gdat.exprtype == 'sdyn':
            gdat.strgenerfull = ['']

    gdat.indxenerfull = arange(len(gdat.strgenerfull) - 1)

    ## energy
    gdat.numbener = gdat.indxenerincl.size
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.indxenerpivt = 0
    if gdat.enerbins:
        gdat.numbenerplot = 100
        gdat.numbenerfull = len(gdat.strgenerfull)
        gdat.strgener = [gdat.strgenerfull[k] for k in gdat.indxenerincl]
        if gdat.enerbinsadje:
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

        gdat.limtener = [amin(gdat.binsener), amax(gdat.binsener)] 
    else:
        gdat.numbenerfull = 1
    if gdat.numbener > 1:
        gdat.enerpivt = gdat.meanener[gdat.indxenerpivt]
    gdat.indxener = arange(gdat.numbener, dtype=int)
    
    # set mock sample vector indices
    setp_varbvalu(gdat, 'maxmnumbelem', 400, popl='full', regi='full')
    setp_varbvalu(gdat, 'minmnumbelem', 0, popl='full', regi='full')
    
    # define maximum and minimum number of elements as lists of arrays
    for strgmodl in gdat.liststrgmodl:
        for strglimt in gdat.liststrglimt:
            numbpopl = getattr(gdat, strgmodl + 'numbpopl')
            indxpopl = getattr(gdat, strgmodl + 'indxpopl')
            limtnumbelem = [[] for l in indxpopl]
            for l in indxpopl:
                if not elemtype[l].startswith('hypr'):
                    limtnumbelem[l] = zeros(gdat.numbregi, dtype=int)
                else:
                    limtnumbelem[l] = zeros(1, dtype=int)
                for d in arange(limtnumbelem[l].size):
                    limtnumbelem[l][d] = getattr(gdat, strgmodl + strglimt + 'numbelempop%dreg%d' % (l, d))
            setattr(gdat, strgmodl + strglimt + 'numbelem', limtnumbelem)
    
    ## hyperparameters
    limtmeanelem = [0.1, 1000.]
    setp_varblimt(gdat, 'meanelem', limtmeanelem, popl='full')
    
    #### boolean flag background
    if gdat.exprtype == 'chan':
        if gdat.anlytype == 'spec':
            specback = [True, True]
        else:
            specback = [False, True]
        setp_varbvalu(gdat, 'specback', specback)
    else:
        for strgmodl in gdat.liststrgmodl:
            backtype = getattr(gdat, strgmodl + 'backtype')
            specback = [False for k in range(len(backtype))]
            setp_varbvalu(gdat, 'specback', specback, strgmodl=strgmodl)
    
    if gdat.exprtype == 'hubb':
        hostemistype = 'sers'
    else:
        hostemistype = 'none'
    setp_varbvalu(gdat, 'hostemistype', hostemistype)

    if gdat.exprtype == 'hubb':
        lensmodltype = 'full'
    else:
        lensmodltype = 'none'
    setp_varbvalu(gdat, 'lensmodltype', lensmodltype)
    
    if gdat.strgexprname == None:
        if gdat.exprtype == 'chan':
            gdat.strgexprname = 'Chandra'
        if gdat.exprtype == 'ferm':
            gdat.strgexprname = 'Fermi-LAT'
        if gdat.exprtype == 'hubb':
            gdat.strgexprname = 'HST'
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

    if gdat.strgenerunit == None:
        if gdat.exprtype == 'ferm':
            gdat.strgenerunit = 'GeV'
            gdat.nameenerunit = 'gevv'
        if gdat.exprtype == 'chan':
            gdat.strgenerunit = 'keV'
            gdat.nameenerunit = 'kevv'
        if gdat.exprtype == 'sdyn':
            gdat.strgenerunit = ''
            gdat.nameenerunit = ''
        if gdat.exprtype == 'hubb':
            gdat.strgenerunit = 'erg'
            gdat.nameenerunit = 'ergs'

    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        elemspatevaltype = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lens':
                elemspatevaltype[l] = 'full'
            else:
                elemspatevaltype[l] = 'loclhash'
        setp_varbvalu(gdat, 'elemspatevaltype', elemspatevaltype, strgmodl=strgmodl)

    gdat.commelemspatevaltype = []
    for strgmodl in gdat.liststrgmodl:
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
        for elemspatevaltypetemp in elemspatevaltype:
            if not elemspatevaltypetemp in gdat.commelemspatevaltype:
                gdat.commelemspatevaltype.append(elemspatevaltypetemp)
    
    if gdat.nameexpr == None:
        if gdat.exprtype == 'ferm':
            gdat.nameexpr = 'Fermi-LAT'
        if gdat.exprtype == 'sdss':
            gdat.nameexpr = 'SDSS'
        if gdat.exprtype == 'chan':
            gdat.nameexpr = 'Chandra'
        if gdat.exprtype == 'hubb':
            gdat.nameexpr = 'HST'
        if gdat.exprtype == 'gaia':
            gdat.nameexpr = 'Gaia'
    
    ## Lensing
    if gdat.anglassc == None:
        if gdat.exprtype == 'ferm':
            gdat.anglassc = 0.5 / gdat.anglfact
        if gdat.exprtype == 'hubb':
            gdat.anglassc = 0.15 / gdat.anglfact
        if gdat.exprtype == 'chan' or gdat.exprtype == 'sdss':
            gdat.anglassc = 0.5 / gdat.anglfact
        if gdat.exprtype == 'sdyn':
            gdat.anglassc = 0.2
    
    if gdat.radispmr == None:
        gdat.radispmr = gdat.anglassc
   
    ### experimental PSFs
    if gdat.exprtype == 'ferm':
        retr_psfnferm(gdat)

        #angltemp = pi * linspace(0., 10., 100) / 180.
        #psfn = retr_psfnferm(meanener, angltemp)
        #fwhm = retr_fwhm(psfn, angl) 
    
    if gdat.numbpixlfull > 1:
        if gdat.exprtype == 'chan':
            retr_psfnchan(gdat)
        if gdat.exprtype == 'sdss':
            retr_psfnsdss(gdat)
        if gdat.exprtype == 'hubb':
            retr_psfnhubb(gdat)
        if gdat.exprtype == 'sdyn':
            retr_psfnsdyn(gdat)
   
    gdat.factburntmpr = 0.75
    gdat.numbburntmpr = gdat.factburntmpr * gdat.numbburn

    # model evaluation approximation error tolerance in units of the fraction of the lowest PS flux
    if gdat.specfraceval == None:
        if gdat.exprtype == 'ferm':
            gdat.specfraceval = 0.5
        else:
            gdat.specfraceval = 0.1

    ### element parameter distributions
    gdat.calcllik = True
    
    gdat.binslgalcart = linspace(gdat.minmlgaldata, gdat.maxmlgaldata, gdat.numbsidecart + 1)
    gdat.binsbgalcart = linspace(gdat.minmbgaldata, gdat.maxmbgaldata, gdat.numbsidecart + 1)
    gdat.meanlgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
    gdat.meanbgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
    
    ### PSF model
    #### angular profile
    if gdat.exprtype == 'ferm':
        gdat.psfntypeexpr = 'doubking'
    if gdat.exprtype == 'chan':
        gdat.psfntypeexpr = 'singking'
    if gdat.exprtype == 'sdss':
        gdat.psfntypeexpr = 'singgaus'
    if gdat.exprtype == 'hubb':
        gdat.psfntypeexpr = 'singgaus'
    if gdat.exprtype == 'sdyn':
        gdat.psfntypeexpr = 'singgaus'
    
    psfntype = gdat.psfntypeexpr
    setp_varbvalu(gdat, 'psfntype', psfntype)
    
    #### background names
    listnameback = ['isot']
    if gdat.exprtype == 'ferm':
        listnameback.append('fdfm')
    if gdat.exprtype == 'chan':
        listnameback.append('part')
    setp_varbvalu(gdat, 'listnameback', listnameback)
    
    # reference elements
    gdat.numbrefr = 0
    if gdat.datatype == 'mock':
        gdat.numbrefr = gdat.truenumbpopl
    if gdat.datatype == 'inpt':
        if gdat.exprtype == 'ferm':
            gdat.numbrefr = 2
        if gdat.exprtype == 'chan':
            gdat.numbrefr = 2
    
    for strgmodl in gdat.liststrgmodl:
        maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        maxmnumbelempopl = [[] for l in indxpopl]
        for l in indxpopl:
            maxmnumbelempopl[l] = sum(maxmnumbelem[l])
        numbregipopl = [[] for l in indxpopl]
        indxregipopl = [[] for l in indxpopl]
        for l in indxpopl:
            numbregipopl[l] = maxmnumbelem[l].size
            indxregipopl[l] = arange(numbregipopl[l])
        setattr(gdat, strgmodl + 'numbregipopl', numbregipopl)
        setattr(gdat, strgmodl + 'indxregipopl', indxregipopl)
        setattr(gdat, strgmodl + 'maxmnumbelempopl', maxmnumbelempopl)
        
    gdat.indxrefr = arange(gdat.numbrefr)
    gdat.listnamefeatamplrefr = [[] for q in gdat.indxrefr]
    gdat.listnamerefr = [] 
    gdat.refrliststrgfeat = [[] for q in gdat.indxrefr]
    gdat.refrliststrgfeatodim = [[] for q in gdat.indxrefr]
    gdat.refrinfo = False
    gdat.listpathwcss = [[] for d in gdat.indxregi]
    gdat.numbpixllgalshft = [[] for d in gdat.indxregi]
    gdat.numbpixlbgalshft = [[] for d in gdat.indxregi]
            
    gdat.refrlistcolr = ['darkgreen', 'olivedrab', 'mediumspringgreen']
    if gdat.allwrefr:
        if gdat.datatype == 'mock':
            gdat.refrinfo = True
            gdat.numbrefr = gdat.truenumbpopl
            gdat.listnamerefr = ['mock' for l in gdat.trueindxpopl] 
            gdat.indxrefr = arange(gdat.numbrefr)
            gdat.listcolrrefr = gdat.refrlistcolr[:gdat.truenumbpopl]
        if gdat.datatype == 'inpt':
            if gdat.exprtype == 'ferm':
                gdat.refrinfo = True
                retr_refrferminit(gdat)
            if gdat.exprtype == 'chan':
                gdat.refrinfo = True
                retr_refrchaninit(gdat)
            
            for q in gdat.indxrefr:
                if 'lgal' in gdat.refrliststrgfeat[q] and 'bgal' in gdat.refrliststrgfeat[q]:
                    gdat.refrliststrgfeat[q] += ['gang', 'aang']
                for strgfeat in gdat.refrliststrgfeat[q]:
                    setattr(gdat, 'refr' + strgfeat, [[[] for d in gdat.indxregi] for q in gdat.indxrefr])
            gdat.refrliststrgfeattotl = retr_listconc(gdat.refrliststrgfeat)
            
            if gdat.exprtype == 'ferm':
                retr_refrfermfinl(gdat)
            if gdat.exprtype == 'chan':
                retr_refrchanfinl(gdat)
            
            gdat.listcolrrefr = ['m', 'y', 'orange', 'c'][:gdat.numbrefr]
    
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        spatdisttype = [[] for l in indxpopl]
        for l in indxpopl:
            spatdisttype[l] = 'unif'
        setp_varbvalu(gdat, 'spatdisttype', spatdisttype, strgmodl=strgmodl)

    if gdat.datatype == 'mock':
        for l in gdat.trueindxpopl:
            if gdat.trueelemtype[l] == 'lens':
                numbelem = 25
            else:
                numbelem = 100
            setp_varbvalu(gdat, 'numbelem', numbelem, popl=l, regi='full', strgmodl='true')
    
        retr_indxsamp(gdat, strgmodl='true')
        gdat.refrliststrgfeattotl = gdat.trueliststrgfeattotl
        for l in gdat.trueindxpopl:
            gdat.listnamefeatamplrefr[l] = gdat.truenamefeatampl[l]
            for strgfeat in gdat.trueliststrgfeatodim[l]:
                gdat.refrliststrgfeat[l].append(strgfeat)
    
    ### background template normalizations
    if gdat.exprtype == 'ferm':
        # Fourier basis
        if 'bfun' in gdat.anlytype:
            setp_varblimt(gdat, 'bacp', [1e-10, 1e10], ener='full', back='full', regi='full')
            print 'gdat.fittnumbback'
            print gdat.fittnumbback
            print

        # physical basis
        else:
            # isotropic + unresolved
            setp_varblimt(gdat, 'bacp', [1e-7, 1e-2], ener=0, back=0, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-9, 1e-3], ener=1, back=0, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-10, 1e-4], ener=2, back=0, regi='full')
            # diffuse
            setp_varblimt(gdat, 'bacp', [1e-6, 1e-2], ener=0, back=1, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-7, 1e-3], ener=1, back=1, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-8, 1e-4], ener=2, back=1, regi='full')
            # dark
            setp_varblimt(gdat, 'bacp', [1e-11, 1e-8], ener=0, back=2, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-11, 1e-8], ener=1, back=2, regi='full')
            setp_varblimt(gdat, 'bacp', [1e-11, 1e-8], ener=2, back=2, regi='full')
    else:
        back = 0
        # sky background + unresolved
        if gdat.exprtype == 'chan':
            if gdat.anlytype == 'spec':
                bacp = [1e0, 1e2]
                setp_varblimt(gdat, 'bacp', bacp, back=0, regi='full')
            else:
                bacp = [1e-1, 1e3]
                setp_varblimt(gdat, 'bacp', bacp, ener='full', back=0, regi='full')
        else:
            if gdat.exprtype == 'hubb':
                bacp = [1e-10, 1e-6]
            # background
            if gdat.exprtype == 'sdyn':
                bacp = [1e-1, 1e1]
            setp_varblimt(gdat, 'bacp', bacp, ener='full', back=0, regi='full')
    
        # particle background
        #if gdat.exprtype == 'chan':
        #    if gdat.anlytype == 'spec':
        #        bacp = [1e-8, 1e-6]
        #    else:
        #        bacp = [1e-1, 1e2]
        #    setp_varblimt(gdat, 'bacp', bacp, back=1, regi='full')
        
    ### element parameter boundaries
    #### spatial
    if gdat.exprtype == 'ferm':
        minmgang = 1e-1 / gdat.anglfact
    else:
        minmgang = 1e-2 / gdat.anglfact
    setp_varbvalu(gdat, 'minmgang', minmgang)
    
    if 'lghtline' in gdat.commelemtype:
        setp_varblimt(gdat, 'elin', gdat.limtener)

    minmdefs = 0.003 / gdat.anglfact
    setp_varbvalu(gdat, 'minmdefs', minmdefs)
    minmdefs = 0.01 / gdat.anglfact
    setp_varbvalu(gdat, 'minmdefs', minmdefs, strgmodl='fitt')
    
    minmnobj = 1e0
    setp_varbvalu(gdat, 'minmnobj', minmnobj)
    
    setp_varblimt(gdat, 'curv', [-1., 1.])

    maxmnobj = 1e3
    setp_varbvalu(gdat, 'maxmnobj', maxmnobj)
    
    maxmdefs = 1. / gdat.anglfact
    setp_varbvalu(gdat, 'maxmdefs', maxmdefs)
   
    # parameter defaults
    ## distribution
    ### flux
    
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
        for l in indxpopl:
            if elemtype[l] == 'lens':
                meandistslop = 1.9
                stdvdistslop = 0.5
                scal = 'gaus'
            else:
                minmdistslop = 0.5
                maxmdistslop = 3.
                scal = 'logt'
            if scal == 'gaus':
                typelimt = 'meanstdv'
                limt = [meandistslop, stdvdistslop]
            else:
                typelimt = 'minmmaxm'
                limt = [minmdistslop, maxmdistslop]
            setp_varblimt(gdat, namefeatampl[l] + 'distslop', limt, popl=l, typelimt=typelimt, strgmodl=strgmodl)
            setp_varbvalu(gdat, 'scal' + namefeatampl[l] + 'distslop', scal, popl=l, strgmodl=strgmodl)

            if elemtype[l] == 'lghtgausbgrd':
                setp_varblimt(gdat, 'gwdtdistslop', [0.5, 4.], popl=l, strgmodl=strgmodl)
                setp_varbvalu(gdat, 'scalgwdtdistslop', 'logt', popl=l, strgmodl=strgmodl)

    if 'lens' in gdat.commelemtype:
        ### projected scale radius
        limtasca = array([0., 0.1]) / gdat.anglfact
        setp_varblimt(gdat, 'asca', limtasca)
        ### projected cutoff radius
        limtacut = array([0., 2.]) / gdat.anglfact
        setp_varblimt(gdat, 'acut', limtacut)
   
    # true model parameters
    
    # temp -- add poisson
    #gdat.truenumbelem = empty(gdat.truenumbpopl, dtype=int)
    #if gdat.truenumbtrap > 0:
    #    for l in gdat.trueindxpopl:
    #       for d in gdat.trueindxregipopl[l]:
    #            getattr(gdat, 'truenumbelempop%dreg%d' % (l, d))
    #            gdat.truenumbelem[l][d] = poisson(gdat.truemeannumbelem[l][d])
    #            if gdat.truenumbelem[l][d] > gdat.truemaxmnumbelem[l][d]:
    #                gdat.truenumbelem[l][d] = gdat.truemaxmnumbelem[l][d]
    
    if gdat.datatype == 'mock':
        gdat.truenumbelem = [zeros(gdat.truenumbregipopl[l], dtype=int) for l in gdat.trueindxpopl]
        for l in gdat.trueindxpopl:
            setattr(gdat, 'truemeanelempop%d' % l, getattr(gdat, 'truenumbelempop%dreg%d' % (l, 0)))
            for d in gdat.trueindxregipopl[l]: 
                gdat.truenumbelem[l][d] = getattr(gdat, 'truenumbelempop%dreg%d' % (l, d))
    
                if gdat.truenumbelem[l][d] > gdat.truemaxmnumbelem[l][d]:
                    raise Exception('')

    setp_varbvalu(gdat, 'scalmeanelem', 'logt')
    
    for strgmodl in gdat.liststrgmodl:
        setp_varblimt(gdat, 'lgal', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl=strgmodl)
        setp_varblimt(gdat, 'bgal', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl=strgmodl)
    
    setp_varblimt(gdat, 'lgalsour', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', regi='full')
    setp_varblimt(gdat, 'bgalsour', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', regi='full')
    setp_varblimt(gdat, 'lgalhost', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', regi='full')
    setp_varblimt(gdat, 'bgalhost', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl='fitt', regi='full')
    
    gdat.stdvhostsour = 0.04 / gdat.anglfact
    setp_varblimt(gdat, 'lgalsour', [0., gdat.stdvhostsour], strgmodl='true', typelimt='meanstdv', regi='full')
    setp_varblimt(gdat, 'bgalsour', [0., gdat.stdvhostsour], strgmodl='true', typelimt='meanstdv', regi='full')
    setp_varblimt(gdat, 'lgalhost', [0., gdat.stdvhostsour], strgmodl='true', typelimt='meanstdv', regi='full')
    setp_varblimt(gdat, 'bgalhost', [0., gdat.stdvhostsour], strgmodl='true', typelimt='meanstdv', regi='full')
    
    setp_varblimt(gdat, 'fluxsour', array([1e-22, 1e-17]), regi='full')
    setp_varblimt(gdat, 'sindsour', array([0., 4.]), regi='full')
    setp_varblimt(gdat, 'sizesour', [0.1 / gdat.anglfact, 2. / gdat.anglfact], regi='full')
    setp_varblimt(gdat, 'ellpsour', [0., 0.5], regi='full')
    setp_varblimt(gdat, 'fluxhost', array([1e-20, 1e-15]), regi='full')
    setp_varblimt(gdat, 'sindhost', array([0., 4.]), regi='full')
    setp_varblimt(gdat, 'sizehost', [0.1 / gdat.anglfact, 4. / gdat.anglfact], regi='full')
    setp_varblimt(gdat, 'beinhost', [0.5 / gdat.anglfact, 2. / gdat.anglfact], regi='full')
    setp_varblimt(gdat, 'ellphost', [0., 0.5], regi='full')
    setp_varblimt(gdat, 'sherextr', [0., 0.1], regi='full')
    setp_varblimt(gdat, 'anglsour', [0., pi], regi='full')
    setp_varblimt(gdat, 'anglhost', [0., pi], regi='full')
    setp_varblimt(gdat, 'serihost', [1., 8.], regi='full')
    setp_varblimt(gdat, 'sangextr', [0., pi], regi='full')
    
    # temp -- to be removed
    #gdat.truefactlgal = gdat.truemaxmlgal - gdat.trueminmlgal
    #gdat.truefactbgal = gdat.truemaxmbgal - gdat.trueminmbgal
    #gdat.trueminmaang = -pi
    #gdat.truemaxmaang = pi
    
    setp_varblimt(gdat, 'aang', [-pi, pi])
   
    # copy the true model to the inference model if the inference model parameter has not been specified
    #temp = deepcopy(gdat.__dict__)
    #for strg, valu in temp.iteritems():
    #    if strg.startswith('true') and not strg[4:].startswith('indx'):
    #        try:
    #            valumodl = getattr(gdat, 'fitt' + strg[4:])
    #            if valumodl == None:
    #                raise
    #            if gdat.verbtype > 1:
    #                print 'Received custom input for ' + strg[4:]
    #        except:
    #            setattr(gdat, 'fitt' + strg[4:], getattr(gdat, strg))
    
    # check inputs
    if gdat.numbburn > gdat.numbswep:
        raise Exception('Bad number of burn-in sweeps.')
    if gdat.factthin > gdat.numbswep - gdat.numbburn or gdat.factthin < 1:
        raise Exception('Bad thinning factor.')
    if gdat.pixltype == 'heal' and gdat.numbspatdims > 2:
        raise Exception('More than 2 spatial dimensions require Cartesian binning.')
    
    if gdat.allwfixdtrue and gdat.datatype == 'mock':
        
        for l in gdat.trueindxpopl:
            if gdat.trueboolelemspat[l]:
                setp_varblimt(gdat, 'lgal', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl=strgmodl, popl=l)
                setp_varblimt(gdat, 'bgal', [-gdat.maxmgangdata, gdat.maxmgangdata], strgmodl=strgmodl, popl=l)
                setp_varbvalu(gdat, 'spatdistcons', 1e-3, popl=l)
                setp_varbvalu(gdat, 'gangdistslop', 1.1, popl=l)
                setp_varbvalu(gdat, 'bgaldistscal', 2. / gdat.anglfact, popl=l)
            if gdat.trueelemtype[l] == 'lghtpnts':
                setp_varbvalu(gdat, 'fluxdistslop', 2.2, popl=l)
            if gdat.trueelemtype[l] == 'lghtline':
                setp_varbvalu(gdat, 'fluxdistslop', 2., popl=l)
            if gdat.trueelemtype[l] == 'lens':
                setp_varbvalu(gdat, 'defsdistslop', 1.9, popl=l)
            if gdat.trueelemtype[l] == 'clus':
                setp_varbvalu(gdat, 'nobjdistslop', 2., popl=l)

            if gdat.trueelemtype[l] == 'lens':
                setp_varbvalu(gdat, 'ascadistmean', 0.05 / gdat.anglfact, popl=l)
                setp_varbvalu(gdat, 'ascadiststdv', 0.04 / gdat.anglfact, popl=l)
                setp_varbvalu(gdat, 'acutdistmean', 1. / gdat.anglfact, popl=l)
                setp_varbvalu(gdat, 'acutdiststdv', 0.04 / gdat.anglfact, popl=l)
            
            if gdat.trueelemtype[l] == 'lghtgausbgrd':
                setp_varbvalu(gdat, 'gwdtdistslop', 2., popl=l)
            
            if gdat.trueboolelemspec[l]:
                if gdat.exprtype == 'ferm':
                    sinddistmean = 2.15
                if gdat.exprtype == 'chan':
                    sinddistmean = 1.
                if gdat.exprtype == 'hubb':
                    sinddistmean = 1.
                setp_varbvalu(gdat, 'sinddistmean', sinddistmean, popl=l)
                setp_varbvalu(gdat, 'sinddiststdv', 0.5, popl=l)
                
                setp_varbvalu(gdat, 'curvdistmean', 2., popl=l)
                setp_varbvalu(gdat, 'curvdiststdv', 0.2, popl=l)
                
                setp_varbvalu(gdat, 'expcdistmean', 2., popl=l)
                setp_varbvalu(gdat, 'expcdiststdv', 0.2, popl=l)
        
            if gdat.trueelemtype[l] == 'lghtpntspuls':
                setp_varbvalu(gdat, 'per0distmean', 3e-3, popl=l)
                setp_varbvalu(gdat, 'per0diststdv', 0.3, popl=l)
                setp_varbvalu(gdat, 'magfdistmean', 10**8.5, popl=l)
                setp_varbvalu(gdat, 'magfdiststdv', 0.7, popl=l)
                setp_varbvalu(gdat, 'dglcdistslop', 2., popl=l)

        if gdat.exprtype == 'ferm':

            #setp_varbvalu(gdat, 'bacp', 5e-6, ener=0, back=0, regi='full')
            setp_varbvalu(gdat, 'bacp', 5e-6, ener=0, back=0, regi='full')
            setp_varbvalu(gdat, 'bacp', 2e-8, ener=1, back=0, regi='full')
            setp_varbvalu(gdat, 'bacp', 2e-9, ener=2, back=0, regi='full')
            #setp_varbvalu(gdat, 'bacp', 1e-5, ener=4, back=0, regi='full')
            #setp_varbvalu(gdat, 'bacp', 7e-7, ener=0, back=1, regi='full')
            setp_varbvalu(gdat, 'bacp', 1e-4, ener=0, back=1, regi='full')
            setp_varbvalu(gdat, 'bacp', 1e-5, ener=1, back=1, regi='full')
            setp_varbvalu(gdat, 'bacp', 7e-7, ener=2, back=1, regi='full')
            #setp_varbvalu(gdat, 'bacp', 3e-8, ener=4, back=1, regi='full')

        else:
            # sky background
            if gdat.exprtype == 'chan':
                if gdat.anlytype == 'spec':
                    bacp = 10.
                else:
                    bacp = 1.
            if gdat.exprtype == 'hubb':
                bacp = 2e-7
            if gdat.exprtype == 'sdyn':
                bacp = 1.
            setp_varbvalu(gdat, 'bacp', bacp, ener='full', back=0, regi='full')

            # particle background
            if gdat.exprtype == 'chan':
                bacp = 70.04
                setp_varbvalu(gdat, 'bacp', bacp, back=1, regi='full')
        
        if gdat.truelensmodltype == 'host' or gdat.truelensmodltype == 'full':
            setp_varbvalu(gdat, 'beinhost', 1.5 / gdat.anglfact, regi='full')
            setp_varbvalu(gdat, 'sizesour', 0.3 / gdat.anglfact, regi='full')
            setp_varbvalu(gdat, 'sizehost', 1. / gdat.anglfact, regi='full')
            setp_varbvalu(gdat, 'ellpsour', 0.2, regi='full')
            setp_varbvalu(gdat, 'fluxsour', 1e-18, regi='full')
            setp_varbvalu(gdat, 'sindsour', 1.5, regi='full')
            setp_varbvalu(gdat, 'fluxhost', 1e-16, regi='full')
            setp_varbvalu(gdat, 'sindhost', 2.5, regi='full')
            setp_varbvalu(gdat, 'ellphost', 0.2, regi='full')
            setp_varbvalu(gdat, 'sangextr', pi / 2., regi='full')
            setp_varbvalu(gdat, 'serihost', 4., regi='full')
    #dictargs['truenumbpntsreg0pop0'] = 25
    #dictargs['truemeanpntspop0'] = 25
    #dictargs['truedefsdistsloppop0'] = 1.9
    #dictargs['truesigcene0evt0'] = 4.21788e-07
    #dictargs['truebacp0000reg0ene0'] = 2e-07
    #dictargs['truefluxsourreg0'] = 1e-18
    #dictargs['truesizesourreg0'] = 1.45444e-06
    #dictargs['trueellpsourreg0'] = 0.2
    #dictargs['trueanglsourreg0'] = 2.4485
    #dictargs['truelgalhostreg0'] = 1.10908e-07
    #dictargs['truebgalhostreg0'] = 2.26346e-08
    #dictargs['truefluxhostreg0'] = 1e-16
    #dictargs['truesizehostreg0'] = 4.84814e-06
    #dictargs['truebeinhostreg0'] = 7.27221e-06
    #dictargs['trueellphostreg0'] = 0.2
    #dictargs['trueanglhostreg0'] = 1.21445
    #dictargs['trueserihostreg0'] = 4
    #dictargs['truesherextrreg0'] = 0.0956653
    #dictargs['truesangextrreg0'] = 1.5708
    #dictargs['truelgalreg0pop00000'] = 8.70681e-06
    #dictargs['truebgalreg0pop00000'] = 5.5522e-06
    #dictargs['truedefsreg0pop00000'] = 4.46996e-07
    #dictargs['trueascareg0pop00000'] = 8.3953e-08
    #dictargs['trueacutreg0pop00000'] = 7.26722e-07
    #dictargs['truelgalreg0pop00001'] = 1.95366e-06
    #dictargs['truebgalreg0pop00001'] = -6.43887e-06
    #dictargs['truedefsreg0pop00001'] = 2.0933e-07

    if gdat.defa:
        return gdat
    
    if gdat.verbtype > 0:
        if gdat.burntmpr:
            print 'Warning: Tempered burn-in.'

    # initial setup
    setpinit(gdat, True) 
    
    gdat.minmgang = 1e-3 * sqrt(2.) * gdat.maxmgangdata
    gdat.maxmgang = sqrt(2.) * gdat.maxmgangdata

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
                if gdat.datatype == 'mock':
                    if strglimt == 'minm':
                        limt = minimum(getattr(gdat, 'fittminm' + namevarb), getattr(gdat, 'trueminm' + namevarb))
                    else:
                        limt = maximum(getattr(gdat, 'fittmaxm' + namevarb), getattr(gdat, 'truemaxm' + namevarb))
                else:
                    limt = getattr(gdat, 'fitt' + strglimt + namevarb)
                setattr(gdat, strglimt + namevarb, limt)
    
    # temp
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
        for l in indxpopl:
            for strgfeat, strgpdfn in zip(liststrgfeatprio[l], liststrgpdfnprio[l]):
                if strgpdfn == 'self':
                    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                    fact = maxm - minm
                    setattr(gdat, strgmodl + 'fact' + strgfeat, fact)
    # intermediate setup

    ## reference spectra
    if gdat.listspecrefrplot == None:
        if gdat.exprtype == 'chan':
            gdat.listenerrefrplot = []
            gdat.listsbrtrefrplot = []
            gdat.listlablrefrplot = ['CDFS', 'HEAO', 'XMM', 'ROSAT']
            for strgfile in ['cdfs', 'heao', 'xmmm', 'rost']:
                path = gdat.pathinpt + '%s.csv' % strgfile
                enerrefrplot = loadtxt(path, delimiter=',')[:, 0]
                sbrtrefrplot = loadtxt(path, delimiter=',')[:, 1]
                sbrtrefrplot /= enerrefrplot**2
                gdat.listenerrefrplot.append(enerrefrplot)
                gdat.listsbrtrefrplot.append(sbrtrefrplot)

    if gdat.numbener > 1:
        if gdat.enerpivt == 0.:
            raise Exception('Pivot energy cannot be zero.')
        # temp
        gdat.factspecener = (gdat.meanener / gdat.enerpivt)**(-sqrt(amin(gdat.fittminmsinddistmeanpop0) * amax(gdat.fittmaxmsinddistmeanpop0)))
        gdat.enerexpcfact = gdat.enerpivt - gdat.meanener
    else:
        gdat.factspecener = array([1.])

    # temp -- this assumes square ROI
    gdat.frambndrmodl = gdat.maxmlgal * gdat.anglfact
    
    if 'lens' in gdat.commelemtype:
        
        if gdat.serstype == 'intp':
            # construct pixel-convolved Sersic surface brightness template
            gdat.factsersusam = 10
            maxmlgal = 2. * sqrt(2.) * gdat.maxmlgal
            gdat.numblgalsers = int(ceil(maxmlgal / gdat.sizepixl))
            gdat.numblgalsersusam = (1 + gdat.numblgalsers) * gdat.factsersusam
            retr_axis(gdat, 'lgalsers', 0., maxmlgal, gdat.numblgalsers)
            retr_axis(gdat, 'lgalsersusam', -gdat.sizepixl / 2., maxmlgal + gdat.sizepixl, gdat.numblgalsersusam)
            retr_axis(gdat, 'bgalsersusam', -gdat.sizepixl / 2., gdat.sizepixl / 2., gdat.factsersusam)
            
            gdat.numbhalfsers = 20
            gdat.numbindxsers = 20
                
            retr_axis(gdat, 'halfsers', gdat.sizepixl, 5. / gdat.anglfact, gdat.numbhalfsers)
            retr_axis(gdat, 'indxsers', 0.5, 10., gdat.numbindxsers)
            
            gdat.binslgalsersusammesh, gdat.binsbgalsersusammesh = meshgrid(gdat.binslgalsersusam, gdat.binsbgalsersusam, indexing='ij')
            gdat.binsradisersusam = sqrt(gdat.binslgalsersusammesh**2 + gdat.binsbgalsersusammesh**2)
             
            gdat.sersprofcntr = empty((gdat.numblgalsers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
            gdat.sersprof = empty((gdat.numblgalsers + 1, gdat.numbhalfsers + 1, gdat.numbindxsers + 1))
            
            for n in range(gdat.numbindxsers + 1):
                for k in range(gdat.numbhalfsers + 1):
                    
                    profusam = retr_sbrtsersnorm(gdat.binsradisersusam, gdat.binshalfsers[k], indxsers=gdat.binsindxsers[n])
    
                    ## take the pixel average
                    indxbgallowr = gdat.factsersusam * (gdat.numblgalsers + 1) / 2
                    indxbgaluppr = gdat.factsersusam * (gdat.numblgalsers + 3) / 2
                    for a in range(gdat.numblgalsers):
                        indxlgallowr = gdat.factsersusam * a
                        indxlgaluppr = gdat.factsersusam * (a + 1) + 1
                        gdat.sersprofcntr[a, k, n] = profusam[(indxlgallowr+indxlgaluppr)/2, 0]
                        gdat.sersprof[a, k, n] = mean(profusam[indxlgallowr:indxlgaluppr, :])
            
            temp, indx = unique(gdat.binslgalsers, return_index=True)
            gdat.binslgalsers = gdat.binslgalsers[indx]
            gdat.sersprof = gdat.sersprof[indx, :, :]
            gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]
    
            indx = argsort(gdat.binslgalsers)
            gdat.binslgalsers = gdat.binslgalsers[indx]
            gdat.sersprof = gdat.sersprof[indx, :, :]
            gdat.sersprofcntr = gdat.sersprofcntr[indx, :, :]

    gdatdictcopy = deepcopy(gdat.__dict__)
    for strg, valu in gdatdictcopy.iteritems():
        if strg.startswith('cmap') and strg[4:] != 'cntpdata' and strg[4:] != 'cntpresi' and strg[4:] != 'cntpmodl':
            retr_ticklabl(gdat, strg[4:])
            
    #liststrgcbar = ['llikmaps', 'perc', 'percresi', 'expo', 'lpdfspatpriointp', 'conv', 'magn', 'deflcomp', 'resiconvelem', 'resimagn']
    #for strgcbar in liststrgcbar:
    #    retr_ticklabl(gdat, strgcbar)
    
    # temp
    #for strgmodl in gdat.liststrgmodl:
    #    for namesele in gdat.listnamesele:
    #        for namefeat in gdat.listnamefeatsele:
    #            for strglimt in gdat.liststrglimt:
    #                try:
    #                    getattr(gdat, strglimt + namefeat + namesele)
    #                except:
    #                    setattr(gdat, strglimt + namefeat + namesele, getattr(gdat, strglimt + namefeat))

    # construct the bins for element features
    for strgmodl in gdat.liststrgmodl:
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        liststrgfeatpriototl = getattr(gdat, strgmodl + 'liststrgfeatpriototl')
        for strgfeat in liststrgfeattotl:
            scal = getattr(gdat, 'scal' + strgfeat + 'plot')
            try:
                maxm = getattr(gdat, 'maxm' + strgfeat)
                minm = getattr(gdat, 'minm' + strgfeat)
            except:
                #scal = getattr(gdat, strgmodl + 'scal' + strgfeat + 'plot')
                maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
            retr_axis(gdat, strgfeat, minm, maxm, gdat.numbbinsplot, scal=scal)
            if strgfeat in liststrgfeatpriototl:
                maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                retr_axis(gdat, strgfeat + 'prio', minm, maxm, gdat.numbbinsplotprio, scal=scal, strginit=strgmodl)
            limt = array([minm, maxm])
            setattr(gdat, 'limt' + strgfeat + 'plot', limt)

    if gdat.fittnumbtrap > 0:
        if gdat.allwrefr and gdat.asscrefr:
            for q in gdat.indxrefr:
                for l in gdat.fittindxpopl:
                    for d in gdat.fittindxregipopl[l]:
                        gdat.listnamevarbscal += ['cmplref%dpop%dreg%d' % (q, l, d)]
                        gdat.listnamevarbscal += ['fdisref%dpop%dreg%d' % (q, l, d)]
    
    # construct bins for the scalar variables
    for namevarbscal in gdat.listnamevarbscal:
        minm = getattr(gdat, 'minm' + namevarbscal)
        maxm = getattr(gdat, 'maxm' + namevarbscal)
        if namevarbscal in gdat.fittnamefixp:
            scal = getattr(gdat, 'fittscal' + namevarbscal)
        else:
            scal = getattr(gdat, 'scal' + namevarbscal)
        if gdat.diagmode:
            if not isinstance(scal, str):
                print 'namevarbscal'
                print namevarbscal
                print 'scal'
                print scal
                raise Exception('Parameter scaling is bad.')
        retr_axis(gdat, namevarbscal, minm, maxm, 50, scal=scal)
    
    sys.stdout = logg(gdat)

    if gdat.verbtype > 0:
        sizetotl = 0.
        for root, dirs, listfile in os.walk(gdat.pathoutp):
            for thisfile in listfile:
                sizetotl += os.path.getsize(root + '/' + thisfile) / 2**30
        if sizetotl > 10.:
            print 'Warning: PCAT data path size is %d GB' % sizetotl

    if gdat.allwrefr:
        # external catalog
        if gdat.datatype == 'inpt':
            gdat.numbrefrpnts = zeros((gdat.numbregi, gdat.numbrefr), dtype=int)
            for q in gdat.indxrefr:
                gdat.numbrefrpnts[d, q] = gdat.refrlgal[q][d].size
        
        if gdat.datatype == 'inpt':
        
            ## rotate element coordinates to the ROI center
            if gdat.pixltype == 'heal' and (gdat.lgalcntr != 0. or gdat.bgalcntr != 0.):
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='ZYX')
                            gdat.refrbgal[q][d][0, :], gdat.refrlgal[0, :] = rttr(pi / 2. - gdat.refrbgal[0, :], gdat.refrlgal[0, :])
                            gdat.refrbgal[q][d][0, :] = pi / 2. - gdat.refrbgal[0, :]

            ## assign zero to nonspecified uncertainties for the reference element features
            for q in gdat.indxrefr:
                for strgfeat in gdat.refrliststrgfeat[q]:
                    if strgfeat == 'gang' or strgfeat == 'aang':
                        continue
                    refrfeat = getattr(gdat, 'refr' + strgfeat)
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            if refrfeat[q][d].ndim == 1:
                                refrfeat[q][d] = tile(refrfeat[q][d], (3, 1)) 
        
            ## preprocess reference element features
            for q in gdat.indxrefr:
                for l in gdat.fittindxpopl:
                    for d in gdat.fittindxregipopl[l]:
                        gdat.refrgang[q][d] = empty((3, gdat.refrlgal[q][d].shape[1]))
                        gdat.refraang[q][d] = empty((3, gdat.refrlgal[q][d].shape[1]))
                        gdat.refrgang[q][d][0, :] = retr_gang(gdat.refrlgal[q][d][0, :], gdat.refrbgal[q][d][0, :])
                        gdat.refraang[q][d][0, :] = retr_aang(gdat.refrlgal[q][d][0, :], gdat.refrbgal[q][d][0, :])
        
    # temp
    #if gdat.refrnumbelem > 0:
    #    gdat.refrfluxbrgt, gdat.refrfluxbrgtassc = retr_fluxbrgt(gdat, gdat.refrlgal, gdat.refrbgal, gdat.refrflux[0, :])

    # local kernel evaluation plot
    if gdat.makeplot:   
        if gdat.fittnumbtrap > 0:
            for l in gdat.fittindxpopl:
                if gdat.fittelemspatevaltype[l] != 'full':
                    plot_eval(gdat, l)
    
    if gdat.verbtype > 1 and gdat.datatype == 'mock':
        print 'True state parameters:'
        for name in gdat.truenamepara:
            print name
            try:
                print getattr(gdat, 'true' + name)
            except:
                print 'failed'
        print

    # generate true data
    if gdat.datatype == 'mock':
        
        if gdat.verbtype > 0:
            print 'Generating mock data...'

        if gdat.seedstat != None:
            if gdat.verbtype > 0:
                print 'Setting the seed for the RNG...'
            set_state(gdat.seedstat)
    
    if gdat.datatype == 'mock':
        ## unit sample vector
        gdat.truesamp = zeros(gdat.truenumbpara)
        gdat.truesampvarb = zeros(gdat.truenumbpara)
        if gdat.truenumbtrap > 0:
            for l in gdat.trueindxpopl:
                gdat.truesampvarb[gdat.trueindxfixpnumbelem[l]] = gdat.truenumbelem[l]
            gdat.truesampvarb[gdat.trueindxfixpmeanelem] = mean(gdat.truenumbelem, axis=0)
        
        if gdat.truenumbtrap > 0:
            if gdat.truenumbtrap > 0:
                gdat.trueindxelemfull = [[[] for d in gdat.trueindxregipopl[l]] for l in gdat.trueindxpopl]
                for l in gdat.trueindxpopl:
                    for d in gdat.trueindxregipopl[l]:
                        gdat.trueindxelemfull[l][d] = range(gdat.truenumbelem[l][d])
                gdat.trueindxsampcomp = retr_indxsampcomp(gdat, gdat.trueindxelemfull, 'true')
            else:
                gdat.trueindxelemfull = []

        gdat.truesamp = rand(gdat.truenumbpara)
        if gdat.truenumbtrap > 0 and gdat.randseedelem:
            seed()
            gdat.truesamp[gdat.indxparatrap] = rand(gdat.truenumbtrap)
       
        for k in gdat.trueindxpara:
            
            if gdat.truenumbtrap > 0 and (k in gdat.trueindxfixpnumbelemtotl or k in gdat.trueindxfixpmeanelem):
                continue

            # assume the true PSF
            if gdat.numbpixl > 1 and k in gdat.trueindxfixppsfp:
                gdat.truesampvarb[k] = gdat.psfpexpr[k-gdat.trueindxfixppsfp[0]]
            else:
                ## read input mock model parameters
                try:
                    # impose user-defined true parameter
                    gdat.truesampvarb[k] = getattr(gdat, 'true' + gdat.truenamepara[k])
                    if gdat.verbtype > 1:
                        print 'Imposing true parameter:'
                        print gdat.truenamepara[k]
                        print gdat.truesampvarb[k]
                        print
                except:
                    # randomly sample the rest of the mock model parameters
                    if k < gdat.truenumbfixp:
                        gdat.truesampvarb[k] = icdf_fixp(gdat, 'true', gdat.truesamp[k], k)
                    else:
                        d = int(gdat.truenamepara[k][-5])
                        l = int(gdat.truenamepara[k][-9])
                        g = (k - gdat.truenumbfixp - gdat.truenumbtrapregipoplcuml[l][d]) % gdat.truenumbcomp[l]

                        gdat.truesampvarb[k] = icdf_trap(gdat, 'true', gdat.truesamp[k], gdat.truesampvarb, gdat.truelistscalcomp[l][g], gdat.trueliststrgcomp[l][g], l, d)

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
    
    if gdat.datatype == 'mock':
        
        gdat.truenumbelemtotl = sum(gdat.truenumbelem)

        # temp
        gdatdictcopy = deepcopy(gdat.__dict__)
        for name, valu in gdatdictcopy.iteritems():
            if name.startswith('true'):
                setattr(gdat, 'refr' + name[4:], valu)

    if gdat.allwrefr and gdat.datatype == 'inpt':

        # rotate reference elements to the spatial coordinate system of PCAT
        # temp -- this does not rotate the uncertainties!
        for q in gdat.indxrefr:
            for l in gdat.fittindxpopl:
                for d in gdat.fittindxregipopl[l]:
                    if len(gdat.listpathwcss[d]) > 0:
                        listhdun = ap.io.fits.open(gdat.listpathwcss[d])
                        wcso = ap.wcs.WCS(listhdun[0].header)
                        skycobjt = ap.coordinates.SkyCoord("galactic", l=gdat.refrlgal[q][d][0, :] * 180. / pi, b=gdat.refrbgal[q][d][0, :] * 180. / pi, unit='deg')
                        rasc = skycobjt.fk5.ra.degree
                        decl = skycobjt.fk5.dec.degree
                        lgal, bgal = wcso.wcs_world2pix(rasc, decl, 0)
                        lgal -= gdat.numbpixllgalshft[d] + gdat.numbsidecart / 2
                        bgal -= gdat.numbpixlbgalshft[d] + gdat.numbsidecart / 2
                        lgal *= gdat.sizepixl
                        bgal *= gdat.sizepixl
                        gdat.refrlgal[q][d][0, :] = lgal
                        gdat.refrbgal[q][d][0, :] = bgal

        # save all reference element features
        for strgfeat in gdat.refrliststrgfeattotl:
            refrfeattotl = [[] for q in gdat.indxrefr]
            for q in gdat.indxrefr:
                for strgfeat in gdat.refrliststrgfeat[q]:
                    refrfeat = getattr(gdat, 'refr' + strgfeat)
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            if len(refrfeat[q][d]) > 0:
                                refrfeattotl[q] = refrfeat[q][d]
            setattr(gdat, 'refr' + strgfeat + 'totl', refrfeattotl)
        
        # find the reference elements inside the ROI
        gdat.indxrefrpntsrofi = [[[] for d in gdat.indxregi] for q in gdat.indxrefr]
        for q in gdat.indxrefr:
            for d in gdat.indxregi:
                gdat.indxrefrpntsrofi[q][d] = where((fabs(gdat.refrlgal[q][d][0, :]) < gdat.maxmgangdata) & (fabs(gdat.refrbgal[q][d][0, :]) < gdat.maxmgangdata))[0]
        for strgfeat in gdat.refrliststrgfeattotl:
            refrfeat = getattr(gdat, 'refr' + strgfeat)
            refrfeatrofi = [[[] for d in gdat.indxregi] for q in gdat.indxrefr]
            for q in gdat.indxrefr:
                for d in gdat.indxregi:
                    if len(refrfeat[q][d]) > 0:
                        refrfeatrofi[q][d] = refrfeat[q][d][..., gdat.indxrefrpntsrofi[q][d]]
            setattr(gdat, 'refr' + strgfeat, refrfeatrofi)
        
        gdat.refrnumbelem = [[] for q in gdat.indxrefr]
        gdat.refrnumbelemtotl = 0
        for q in gdat.indxrefr:
            gdat.refrnumbelem[q] = zeros(gdat.numbregi, dtype=int)
            for d in gdat.indxregi:
                gdat.refrnumbelem[q][d] = gdat.refrlgal[q][d].shape[1]
            gdat.refrnumbelemtotl += sum(gdat.refrnumbelem[q]) 
        
        ## check that all reference element features are finite
        for d in gdat.indxregi:
            for q in gdat.indxrefr:
                for strgfeat in gdat.refrliststrgfeat[q]:
                    refrfeat = getattr(gdat, 'refr' + strgfeat)
                    if len(refrfeat[q][d]) > 0:
                        indx = where(logical_not(isfinite(refrfeat[q][d])))[0]
                        if indx.size > 0:
                            refrfeat[q][d][indx] = 0.
                            if gdat.verbtype > 0:
                                print 'Warning: Provided reference element feature is not finite. Defaulting to 0...'
        
        # bin reference element features
        for q in gdat.indxrefr:
            for d in gdat.indxregi:
                for strgfeat in gdat.refrliststrgfeat[q]:
                    refrfeat = getattr(gdat, 'refr' + strgfeat)
                    if len(refrfeat[q][d]) > 0:
                        bins = getattr(gdat, 'bins' + strgfeat)
                        hist = histogram(refrfeat[q][d][0, :], bins)[0]
                        setattr(gdat, 'refrhist' + strgfeat + 'ref%dreg%d' % (q, d), hist)
    
    if gdat.fittnumbtrap > 0:
        # plot settings
        ## upper limit of histograms
        if gdat.allwrefr:
            gdat.limtydathistfeat = [0.5, 10**ceil(log10(max(gdat.fittmaxmnumbelemtotl, gdat.refrnumbelemtotl)))]
        else:
            gdat.limtydathistfeat = [0.5, 10**ceil(log10(gdat.fittmaxmnumbelemtotl))]
            #if gdat.datatype == 'mock':
            #    gdat.limtydathistfeat = [0.5, sum(gdat.truenumbelem) / 2.]
            #if gdat.datatype == 'inpt':

    if gdat.datatype == 'mock':
        if lensmodltype != 'none':
            proc_samp(gdat, None, 'this', 'true', raww=True)
        proc_samp(gdat, None, 'this', 'true')
        
        if gdat.makeplot and gdat.makeplotinit:
            plot_samp(gdat, None, 'this', 'true')
        
    ## element feature indices ordered with respect to the amplitude variable
    refrfeatsort = [[[] for d in gdat.indxregi] for q in gdat.indxrefr]
    if not (gdat.datatype == 'mock' and gdat.truenumbtrap == 0):
        for q in gdat.indxrefr:
            refrfeatampl = getattr(gdat, 'refr' + gdat.listnamefeatamplrefr[q])
            for d in gdat.indxregi:
                if len(refrfeatampl[q][d]) > 0:
                    indxelem = argsort(refrfeatampl[q][d][0, :])[::-1]
                    for strgfeat in gdat.refrliststrgfeat[q]:
                        refrfeat = getattr(gdat, 'refr' + strgfeat)
                        if len(refrfeat[q][d]) > 0:
                            refrfeatsort[q][d] = refrfeat[q][d][..., indxelem]
        setattr(gdat, 'refr' + strgfeat, refrfeatsort)
        
    # initial plots
    if gdat.makeplot and gdat.makeplotinit:
        plot_init(gdat)

    # final setup
    setpfinl(gdat, True) 
    
    # exit before running the sampler
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
    
    # write the numpy RNG state to file
    with open(gdat.pathoutpthis + 'stat.p', 'wb') as thisfile:
    	cPickle.dump(get_state(), thisfile)
    
    # start the timer
    gdat.timerealtotl = time.time()
    gdat.timeproctotl = time.clock()
   
    if gdat.verbtype > 1:
        for strgmodl in gdat.liststrgmodl:
            indxpopl = getattr(gdat, strgmodl + 'indxpopl')
            elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
            if 'locl' in elemspatevaltype:
                print 'maxmangleval'
                print gdat.anglfact * gdat.maxmangleval[l], ' [%s]' % gdat.strganglunit
            
    # process lock for simultaneous plotting
    lock = mp.Manager().Lock()
        
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdat, 'gdat')
    
    # list of variables for which the posterior is calculated at each sweep
    gdat.liststrgvarbarryswep = ['memoresi', 'lpri', 'lfctasym', 'lpriprop', 'lpau', 'deltlliktotl', 'lliktotl', 'chro', 'accpprob', \
                                                                    'accp', 'accppsfn', 'accpprio', 'accpprop', 'indxproptype', 'amplpert']
    if gdat.probspmr > 0.:
        gdat.liststrgvarbarryswep += ['auxipara', 'lcomfact', 'ljcbfact']
    
    # perform a fudicial processing of a sample vector in order to find the list of variables for which the posterior will be calculated
    if gdat.verbtype > 0:
        print 'Processing a fudicial sample...'
    gdatmodifudi = tdpy.util.gdatstrt()
    gdatmodifudi.thischro = zeros(gdat.numbchro)
    gdatmodifudi.thissamp = rand(gdat.fittnumbpara)
    if gdat.fittnumbtrap > 0:
        for l in gdat.fittindxpopl:
            gdatmodifudi.thissamp[gdat.fittindxfixpnumbelem[l]] = 1
    
    retr_elemlist(gdat, gdatmodifudi)
    gdatmodifudi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodifudi.thissamp, gdatmodifudi.thisindxsampcomp)
    proc_samp(gdat, gdatmodifudi, 'this', 'fitt')
    gdat.liststrgvarbarrysamp = []
    gdat.liststrgvarblistsamp = []
    for strg, valu in gdatmodifudi.__dict__.iteritems():
        if strg.startswith('this') and not strg[4:] in gdat.liststrgvarbarryswep:
            if isinstance(valu, ndarray):
                gdat.liststrgvarbarrysamp.append(strg[4:])
            elif isinstance(valu, list) and strg != 'thisindxsampcomp' and strg != 'thispsfnconv' and \
                                                                                            strg != 'thistrueindxelemasscmiss' and strg != 'thistrueindxelemasschits':
                gdat.liststrgvarblistsamp.append(strg[4:])
    
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
    
    setp_indxswepsave(gdat)
    
    gdat.optitypetemp = 'none'
   
    # perform an initial run, sampling from the prior
    if gdat.checprio:
        
        if gdat.verbtype > 0:
            print 'Sampling from the prior...'
        
        gdat.namesampdist = 'prio'
        gdat.legdsampdist = 'Prior'
        gdat.calcllik = False
        
        if gdat.optitype != 'none':
            workopti(gdat, lock)

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
    
    if gdat.optitype != 'none':
        workopti(gdat, lock)
    
    # run the sampler
    worksamp(gdat, lock)
    
    #prid = os.fork()
    #if prid:
    #    print 'Parent process returning. Child will continue...'
    #    return
    #    #os._exit(0)
    
    # post process the samples
    proc_post(gdat)
    
    # make animations
    if gdat.makeanim:
        make_anim(gdat.rtag)

    if gdat.verbtype > 0:
        print 'The ensemble of catalogs is at ' + gdat.pathoutpthis
        if gdat.makeplot:
            print 'The plots are at ' + gdat.pathplot
        print 'PCAT has run successfully. Returning to the OS...'
    
    return gdat


def workopti(gdat, lock):

    gdat.optitypetemp = gdat.optitype
    inittypesave = gdat.inittype 

    # estimate the covariance
    gdat.numbprocsave = gdat.numbproc
    gdat.numbproc = 1

    if gdat.verbtype > 0:
        print 'Optimizing proposal scale...'
    
    worksamp(gdat, lock)
    
    if gdat.verbtype > 0:
        print 'Writing the estimated covariance matrix to the disc...'
    
    gdat.optitypetemp = 'none'
    
    thisfile = h5py.File(gdat.pathoutpthis + 'opti.h5', 'r')
    gdat.stdvstdp = thisfile['stdvstdp'][()]
    thisfile.close()
    gdat.numbproc = gdat.numbprocsave
    
    gdat.inittype = inittypesave
    

def initarry( \
             dictvarbvari, \
             dictvarb, \
             sameseed=True, \
             execpara=False, \
             nameconfexec=None, \
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

    numbiter = len(dictvarbvari)

    if liststrgvarboutp != None:
        numboutp = len(liststrgvarboutp)
        dictoutp = dict()
        for strgvarb in liststrgvarboutp:
            dictoutp[strgvarb] = [[] for k in range(numbiter)]
    
    dictvarb['seedstat'] = seedstat
    dictvarb['randseedelem'] = randseedelem
    dictvarb['strgcnfg'] = strgcnfg
    
    listgdat = []
    listprid = []
    for k, nameconf in enumerate(dictvarbvari):
        
        if nameconfexec != None:
            if nameconf != nameconfexec:
                continue

        dictvarbtemp = deepcopy(dictvarb)
        for strgvarb, valu in dictvarbvari[nameconf].iteritems():
            dictvarbtemp[strgvarb] = valu
        dictvarbtemp['strgcnfg'] = inspect.stack()[1][3] + '_' + nameconf
        
        if execpara:
            prid = os.fork()
            listprid.append(prid)
            if prid == 0:
                #print '%d slept.' % k
                #time.sleep(5)
                #print '%d woke up.' % k
                listgdat.append(init(**dictvarbtemp))
                if liststrgvarboutp != None:
                    for strgvarb in liststrgvarboutp:
                        dictoutp[strgvarb][k] = getattr(listgdat[k], strgvarb)
                os._exit(0)
            if k % 4 == 3 or k == len(dictvarbvari) - 1:
                #print 'parent waiting...'
                for kk in range(len(listprid)):
                    os.waitpid(listprid[kk], 0)
                listprid = []
        else:
            listgdat.append(init(**dictvarbtemp))
            
    #print 'finished'
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
        print 'Post processing'
        timeinit = gdat.functime()
   
    # read the chains
    if gdat.verbtype > 0:
        print 'Reading chains from disk...'
        timeinit = gdat.functime()
    listgdatmodi = []
    for k in gdat.indxproc:
        path = gdat.pathoutpthis + 'gdatmodi%04d' % k
        listgdatmodi.append(readfile(path))
        os.system('rm -rf %s*' % path)
    os.system('rm -rf %s' % gdat.pathoutpthis + 'gdatinit*')
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)
    
    # aggregate samples from the chains
    if gdat.verbtype > 0:
        print 'Accumulating samples from all processes...'
        timeinit = gdat.functime()
    # arrays 
    gdat.liststrgvarbarry = gdat.liststrgvarbarrysamp + gdat.liststrgvarbarryswep
    gdat.liststrgchan = gdat.liststrgvarbarry + ['fixp'] + gdat.liststrgvarblistsamp
    for strgvarb in gdat.liststrgvarbarry:
        for k in gdat.indxproc:
            if k == 0:
                shap = getattr(listgdatmodi[k], 'list' + strgvarb).shape
                shap = [shap[0], gdat.numbproc] + list(shap[1:])
                temp = zeros(shap) - 1
            if len(shap) > 2:
                temp[:, k, :] = getattr(listgdatmodi[k], 'list' + strgvarb)
            else:
                temp[:, k] = getattr(listgdatmodi[k], 'list' + strgvarb)
        setattr(gdat, 'list' + strgvarb, temp)
    
    # lists of lists collected at each sample
    for strgvarb in gdat.liststrgvarblistsamp:
        listtemp = [[[] for k in gdat.indxproc] for j in gdat.indxsamp]
        for j in gdat.indxsamp:      
            for k in gdat.indxproc:
                listtemp[j][k] = getattr(listgdatmodi[k], 'list' + strgvarb)[j]
        setattr(gdat, 'list' + strgvarb, listtemp)

    ## maximum likelihood sample 
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
        if gdat.verbtype > 0:
            print 'Computing the Gelman-Rubin TS...'
            timeinit = gdat.functime()
        gdat.gmrbfixp = zeros(gdat.fittnumbfixp)
        gdat.gmrbstat = [[zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))] for d in gdat.indxregi]
        for k in gdat.fittindxfixp:
            gdat.gmrbfixp[k] = tdpy.mcmc.gmrb_test(gdat.listsampvarb[:, :, k])
            if not isfinite(gdat.gmrbfixp[k]):
                gdat.gmrbfixp[k] = 0.
        for d in gdat.indxregi:
            listcntpmodl = getattr(gdat, 'listcntpmodlreg%d' % d)
            for i in gdat.indxener:
                for j in gdat.indxpixl:
                    for m in gdat.indxevtt:
                        gdat.gmrbstat[d][i, j, m] = tdpy.mcmc.gmrb_test(listcntpmodl[:, :, i, j, m])
        if gdat.verbtype > 0:
            timefinl = gdat.functime()
            print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # calculate the autocorrelation of the chains
    if gdat.verbtype > 0:
        print 'Computing the autocorrelation of the chains...'
        timeinit = gdat.functime()
    gdat.atcrcntp = [[] for d in gdat.indxregi]
    gdat.timeatcrcntp = [[] for d in gdat.indxregi]
    for d in gdat.indxregi:
        gdat.atcrcntp[d] = empty((gdat.numbproc, gdat.numbener, gdat.numbpixl, gdat.numbevtt, gdat.numbsamp / 2))
        gdat.timeatcrcntp[d] = empty((gdat.numbproc, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    gdat.atcrpara = empty((gdat.numbproc, gdat.fittnumbpara, gdat.numbsamp / 2))
    gdat.timeatcrpara = empty((gdat.numbproc, gdat.fittnumbpara))
    for k in gdat.indxproc:
        gdat.atcrpara[k, :, :], gdat.timeatcrpara[k, :] = 0., 0.#tdpy.mcmc.retr_timeatcr(gdat.listsampvarb[:, k, :], verbtype=gdat.verbtype)
        for d in gdat.indxregi:
            listcntpmodl = getattr(gdat, 'listcntpmodlreg%d' % d)
            gdat.atcrcntp[d][k, :], gdat.timeatcrcntp[d][k, :] = 0., 0.#tdpy.mcmc.retr_timeatcr(listcntpmodltemp[:, k, :, :, :], verbtype=gdat.verbtype)
    gdat.timeatcrmaxm = 0.
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # flatten the chains from different walkers
    ## lists collected at each sample
    for strgvarb in gdat.liststrgvarblistsamp:
        listtemp = []
        listinpt = getattr(gdat, 'list' + strgvarb)
        for j in gdat.indxsamp:      
            for k in gdat.indxproc:
                listtemp.append(listinpt[j][k])
        setattr(gdat, 'list' + strgvarb, listtemp)

    ## list of other parameters to be flattened
    gdat.liststrgvarbarryflat = deepcopy(gdat.liststrgvarbarry)
    for strg in ['deltlliktotl', 'memoresi']:
        gdat.liststrgvarbarryflat.remove(strg)
   
    gdat.listsampvarbproc = copy(gdat.listsampvarb)

    ## other parameters
    for strgvarb in gdat.liststrgvarbarryflat:
        inpt = getattr(gdat, 'list' + strgvarb)
        if inpt.shape[2:] == (1,):
            shap = [inpt.shape[0] * inpt.shape[1]]
        else:
            shap = [inpt.shape[0] * inpt.shape[1]] + list(inpt.shape[2:])
        setattr(gdat, 'list' + strgvarb, inpt.reshape(shap))
    
    # maximum likelihood sample
    indxsamptotlmlik = argmax(sum(sum(sum(sum(gdat.listllik, 4), 3), 2), 1))
    
    # copy the maximum likelihood sample
    for strgvarb in gdat.liststrgvarbarrysamp:
        setattr(gdat, 'mlik' + strgvarb, getattr(gdat, 'list' + strgvarb)[indxsamptotlmlik, ...])
    for strgvarb in gdat.liststrgvarblistsamp:
        setattr(gdat, 'mlik' + strgvarb, getattr(gdat, 'list' + strgvarb)[indxsamptotlmlik])

    # temp -- dont gdat.listllik and gdat.listsamp have the same dimensions?
    gdat.mliksamp = gdat.listsamp[indxsamptotlmlik, :]
    gdat.mliksampvarb = gdat.listsampvarb[indxsamptotlmlik, :]
    if gdat.fittnumbtrap > 0:
        gdat.mlikindxelemfull = gdat.listindxelemfull[indxsamptotlmlik]
    gdat.mlikfixp = gdat.mliksampvarb[gdat.fittindxfixp]
    for k, namefixp in enumerate(gdat.fittnamefixp):
        setattr(gdat, 'mlik' + namefixp, gdat.mlikfixp[k])

    # add execution times to the chain output
    gdat.timereal = zeros(gdat.numbproc)
    gdat.timeproc = zeros(gdat.numbproc)
    for k in gdat.indxproc:
        gdat.timereal[k] = listgdatmodi[k].timereal
        gdat.timeproc[k] = listgdatmodi[k].timeproc
    
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
        if gdat.fittboolelemspatanyy:
            gdat.posthistlgalbgalelemstkd = [[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl]
        
        for l in gdat.fittindxpopl:
            if gdat.fittelemtype[l] != 'lghtline':
                numb = len(gdat.fittliststrgfeatsign[l])
                for d in gdat.fittindxregipopl[l]:
                    gdat.posthistlgalbgalelemstkd[l][d] = zeros((gdat.numbbgalpntsprob, gdat.numblgalpntsprob, gdat.numbbinsplot, numb))
                    temparry = concatenate([gdat.listlgal[n][l][d] for n in gdat.indxsamptotl])
                    temp = empty((len(temparry), 3))
                    temp[:, 0] = temparry
                    temp[:, 1] = concatenate([gdat.listbgal[n][l][d] for n in gdat.indxsamptotl])
                    for k, strgfeat in enumerate(gdat.fittliststrgfeatsign[l]):
                        temp[:, 2] = concatenate([getattr(gdat, 'list' + strgfeat)[n][l][d] for n in gdat.indxsamptotl])
                        bins = getattr(gdat, 'bins' + strgfeat)
                        gdat.posthistlgalbgalelemstkd[l][d][:, :, :, k] = histogramdd(temp, bins=(gdat.binslgalpntsprob, gdat.binsbgalpntsprob, bins))[0]

    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    ## construct a condensed catalog of elements
    if gdat.condcatl and gdat.fittnumbtrap > 0:
        
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
    
    # compute credible intervals
    if gdat.verbtype > 0:
        print 'Computing credible intervals...'
        timeinit = gdat.functime()
    
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
            numb = len(listtemp[0])
            
            for k in range(numb):
                if isinstance(listtemp[0][k], list):
                    continue
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
            
            if strgchan in gdat.listnamevarbcpos:
                cpostemp = empty([listtemp.shape[0]] + [3] + list(listtemp.shape[1:]))
                for n in gdat.indxsamptotl:
                    cpostemp[n, ...] = tdpy.util.retr_postvarb(listtemp[:n+1, ...])

        setattr(gdat, 'post' + strgchan, posttemp)
        setattr(gdat, 'medi' + strgchan, meditemp)
        setattr(gdat, 'errr' + strgchan, errrtemp)
        setattr(gdat, 'stdv' + strgchan, stdvtemp)
        if strgchan in gdat.listnamevarbcpos:
            setattr(gdat, 'cpos' + strgchan, cpostemp)
    
    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)
    
    # memory usage
    gdat.meanmemoresi = mean(gdat.listmemoresi, 1)
    gdat.derimemoresi = (gdat.meanmemoresi[-1] - gdat.meanmemoresi[0]) / gdat.numbswep
   
    gdat.pathpcat = gdat.pathoutpthis + 'pcat'
    writoutp(gdat, gdat.pathpcat)
    os.system('rm -rf %sgdat*' % gdat.pathoutpthis) 
   
    # write the final gdat object
    if True:
        path = gdat.pathoutpthis + 'gdatfinl'
        if gdat.verbtype > 0:
            print 'Writing the global object to %s...' % path
        writfile(gdat, path) 
    
    if gdat.makeplot and gdat.makeplotpost:
        plot_post(gdat, prio=prio)

    gdat.timerealtotl = time.time() - gdat.timerealtotl
    gdat.timeproctotl = time.clock() - gdat.timeproctotl
    gdat.timeproctotlswep = gdat.timeproctotl / gdat.numbswep
    if gdat.timeatcrmaxm == 0.:
        gdat.timeprocnorm = 0.
    else:
        gdat.timeprocnorm = gdat.timeproctotlswep / gdat.timeatcrmaxm
    
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
        indx = setdiff1d(gdat.fittindxpara, gdat.fittindxfixpnumbelemtotl)
        indxbadd = where((gdatmodi.thissamp[indx] < 0.) | (gdatmodi.thissamp[indx] > 1.))[0]
        indxbadd = indx[indxbadd]
    else:
        indxbadd = where((gdatmodi.thissamp < 0.) | (gdatmodi.thissamp > 1.))[0]
    if indxbadd.size > 0:
        print '%s went outside prior bounds when perturbing...' % gdat.fittnamepara[indxbadd]
        print 'gdatmodi.thissamp[indx[indxbadd]]'
        print gdatmodi.thissamp[indx[indxbadd]]
        print 'gdatmodi.thissamp[indxparapert[k]]'
        print gdatmodi.thissamp[indxparapert[k]]
        print 'stdvparapert[k]'
        print stdvparapert[k]
        raise Exception('')
        deltlpos = 0.
    
    else:
        # temp
        gdatmodi.nextsamp = copy(gdatmodi.thissamp)
        gdatmodi.nextsampvarb = zeros_like(gdatmodi.thissamp) 
       
        if gdat.fittnumbtrap > 0:
            for k in gdat.fittindxfixpdist:
                gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, 'fitt', gdatmodi.thissamp[k], k)
            # rescale element features
            for k in range(numbpert):
                rscl_elem(gdat, gdatmodi.thissampvarb, gdatmodi.thisindxsampcomp, gdatmodi.nextsampvarb, gdatmodi.nextsamp, indxparapert[k])
        
        gdatmodi.thissamp = copy(gdatmodi.nextsamp)
        gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)

        proc_samp(gdat, gdatmodi, 'this', 'fitt', fast=True)
   
        deltlpos = gdatmodi.thislpostotl - gdatmodi.thislpostotltemp
        deltlpri = gdatmodi.thislpritotl - gdatmodi.thislpritotltemp
        deltllik = gdatmodi.thislliktotl - gdatmodi.thislliktotltemp
        
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
    
    deltparastep = 1e-4

    maxmstdv = 0.1
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
    
    if gdat.fittnumbtrap > 0:
        gdatmodi.dictmodi = [[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl]
        for l in gdat.fittindxpopl:
            for d in gdat.fittindxregipopl[l]:
                gdatmodi.dictmodi[l][d] = dict()
                gdatmodi.dictmodi[l][d][gdat.fittnamefeatampl[l] + 'indv'] = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.fittnamefeatampl[l]][l][d]]
                for strgcomp in gdat.fittliststrgcomp[l]:
                    gdatmodi.dictmodi[l][d]['stdv' + strgcomp + 'indv'] = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l][d]]
            
    gdatmodi.cntrparasave = 0
    lliktemp = empty(gdat.numbstdp)
    cntr = zeros(gdat.fittmaxmnumbcomp)
    gdatmodi.stdvstdpmatr = zeros((gdat.numbstdp, gdat.numbstdp)) 
    gdatmodi.hess = zeros((gdat.numbstdp, gdat.numbstdp)) 
    deltlpos = zeros((3, 3))

    deltlpos[1, 1] = retr_deltlpos(gdat, gdatmodi, array([0]), array([0.]))
    
    if gdat.fittnumbtrap > 0:
        gdatmodi.thisindxsampcompconc = concatenate([concatenate(gdatmodi.thisindxsampcomp['comp'][l]) for l in gdat.fittindxpopl])
    if gdat.propcomp:
        indxsamptranprop = gdatmodi.thisindxsampcompconc
    else:
        indxsamptranprop = []
    
    for k in gdat.fittindxpara:
        if k in gdat.indxfixpprop or k in indxsamptranprop:
            indxstdpfrst = gdat.indxstdppara[k]
            for n in gdat.fittindxpara:
                if n in gdat.indxfixpprop or n in indxsamptranprop:
                    indxstdpseco = gdat.indxstdppara[n]
                    if k == n:
                        
                        for a in [0, 2]:
                            # evaluate the posterior
                            deltlpos[a, 1] = retr_deltlpos(gdat, gdatmodi, array([k]), array([diffparaodim[a]]))
                    
                        gdatmodi.hess[indxstdpfrst, indxstdpseco] = 1. / 4. / deltparastep**2 * fabs(deltlpos[0, 1] + deltlpos[2, 1] - 2. * deltlpos[1, 1])
                        
                        if gdat.fittnumbtrap > 0:
                            if k in gdatmodi.thisindxsampcompconc:
                                indxtrapmoditemp = k - gdat.fittindxsampcompinit
                                indxpoplmoditemp = array([amin(where(indxtrapmoditemp // gdat.fittnumbtrappoplcumr == 0))])
                                indxregimoditemp = array([amin(where(indxtrapmoditemp // gdat.fittnumbtrapregipoplcumr[indxpoplmoditemp[0]] == 0))])
                                numbparapoplinittemp = indxtrapmoditemp - gdat.fittnumbtrapregipoplcuml[indxpoplmoditemp[0]][indxregimoditemp[0]]
                                indxelemmoditemp = [numbparapoplinittemp // gdat.fittnumbcomp[indxpoplmoditemp[0]]]
                                indxcompmoditemp = numbparapoplinittemp % gdat.fittnumbcomp[indxpoplmoditemp[0]]
                                
                                indxsampampltemp = k - indxcompmoditemp + gdat.fittindxcompampl[indxpoplmoditemp[0]]
                                amplfact = gdatmodi.thissampvarb[indxsampampltemp] / getattr(gdat, 'minm' + gdat.fittnamefeatampl[indxpoplmoditemp[0]])
                                
                                stdv = 1. / sqrt(gdatmodi.hess[indxstdpfrst, indxstdpseco])
                                if not isfinite(stdv):
                                    #raise Exception('Hessian is infinite.')
                                    print 'Hessian is infinite. Replacing with unity.'
                                    print 'deltlpos'
                                    print deltlpos
                                    print
                                    stdv = 1.
                                if strgcomp == gdat.fittnamefeatampl[indxpoplmoditemp[0]]:
                                    gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * amplfact**2. / \
                                                                                gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[indxpoplmoditemp[0]][indxregimoditemp[0]]]
                                else:
                                    gdatmodi.stdvstdpmatr[indxstdpfrst, indxstdpseco] += stdv * amplfact**0.5 / \
                                                                                gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[indxpoplmoditemp[0]][indxregimoditemp[0]]]
                                gdatmodi.dictmodi[l][d]['stdv' + strgcomp + 'indv'][indxelemmoditemp] = stdv
                                gdatmodi.dictmodi[l][d][gdat.fittnamefeatampl[l] + 'indv'][indxelemmoditemp] = gdatmodi.thissampvarb[indxsampampltemp]
        
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
    
    gdatmodi.stdvstdpmatr *= fudgstdv * 2.38 
    
    indx = where(logical_not(isfinite(gdatmodi.stdvstdpmatr)))
    for k in range(indx[0].size):
        if indx[0][k] == indx[1][k]:
            print 'Bad estimation of the proposal scale'
            print gdat.namestdp[indx[0][k]]
            print

        #gdatmodi.stdvstdpmatr[indx] = maxmstdv
        gdatmodi.stdvstdpmatr[indx] = gdat.stdvstdp[getattr(gdat, 'indxstdp' + gdat.namestdp[indx[0][k]])]
            
    proc_samp(gdat, gdatmodi, 'this', 'fitt')
    
    gdatmodi.stdvstdp = gdatmodi.stdvstdpmatr[gdat.indxstdp, gdat.indxstdp]
    
    if gdat.makeplot:
        
        xdat = gdat.indxstdp
        ydat = gdatmodi.stdvstdp
        
        pathopti = getattr(gdat, 'path' + gdat.namesampdist + 'opti')
        path = pathopti + 'stdv%d.pdf' % gdatmodi.indxprocwork
        tdpy.util.plot_gene(path, xdat, ydat, scalydat='logt', lablxdat='$i_{stdp}$', lablydat=r'$\sigma$', plottype='hist', limtydat=[amin(ydat) / 2., 2. * amax(ydat)])
        
        if gdat.fittnumbtrap > 0:
            for l in gdat.fittindxpopl:
                meanplot = getattr(gdat, 'mean' + gdat.fittnamefeatampl[l])
                minm = getattr(gdat, 'minm' + gdat.fittnamefeatampl[l])
                for d in gdat.fittindxregipopl[l]:
                    for strgcomp in gdat.fittliststrgcomp[l]:
                        path = pathopti + 'stdv' + strgcomp + 'pop%dreg%d.pdf' % (l, d)
                        factplot = getattr(gdat, 'fact' + strgcomp + 'plot')
                        factamplplot = getattr(gdat, 'fact' + gdat.fittnamefeatampl[l] + 'plot')
                        xdat = [gdatmodi.dictmodi[l][d][gdat.fittnamefeatampl[l] + 'indv'] * factamplplot, meanplot * factamplplot]
                        
                        if strgcomp == gdat.fittnamefeatampl[l]:
                            ydat = [gdatmodi.dictmodi[l][d]['stdv' + strgcomp + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)] / (meanplot / minm)**2.]
                        else:
                            ydat = [gdatmodi.dictmodi[l][d]['stdv' + strgcomp + 'indv'], gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)] / (meanplot / minm)**0.5]
                        lablxdat = getattr(gdat, 'labl' + gdat.fittnamefeatampl[l] + 'totl')
                        scalxdat = getattr(gdat, 'scal' + gdat.fittnamefeatampl[l] + 'plot')
                        limtxdat = array(getattr(gdat, 'limt' + gdat.fittnamefeatampl[l] + 'plot')) * factamplplot
                        tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                         lablydat=r'$\sigma_{%s}$' % getattr(gdat, 'labl' + strgcomp), plottype=['scat', 'lghtline'])
                        #tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                        #                                 lablydat=r'$\sigma_{%s}$%s' % (getattr(gdat, 'labl' + strgcomp), getattr(gdat, 'labl' + strgcomp + 'unit')), plottype=['scat', 'lghtline'])
                        
                        tdpy.util.plot_gene(path, xdat, ydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, limtxdat=limtxdat, \
                                                         lablydat=r'$\sigma_{%s}$' % getattr(gdat, 'labl' + strgcomp), plottype=['scat', 'lghtline'])

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


def retr_elemlist(gdat, gdatmodi):

    ## element parameters
    if gdat.fittnumbtrap > 0:
        ## lists of occupied and empty transdimensional parameters
        gdatmodi.thisindxelemfull = [[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl]
        cntr = 0
        for l in gdat.fittindxpopl:
            for d in gdat.fittindxregipopl[l]:
                cntr += gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]]
                gdatmodi.thisindxelemfull[l][d] = range(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]].astype(int))
        if cntr > 0:
            gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxelemfull, 'fitt')
        else: 
            gdatmodi.thisindxsampcomp = None
    else:
        gdatmodi.thisindxsampcomp = None
    

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
        print 'inittype'
        print gdat.inittype
  
    gdatmodi.thisamplpert = zeros(1)
    
    gdatmodi.thisstdpscalfact = 1.

    ## unit sample vector
    ## initialize randomly
    gdatmodi.thissamp = rand(gdat.fittnumbpara)
    gdatmodi.thissampvarb = zeros(gdat.fittnumbpara)
    
    if gdat.fittnumbtrap > 0:
        for l in gdat.fittindxpopl:
            meanelemtemp = icdf_fixp(gdat, 'fitt', gdatmodi.thissamp[gdat.fittindxfixpmeanelem[l]], gdat.fittindxfixpmeanelem[l])
            for d in gdat.fittindxregipopl[l]:
                gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] = poisson(meanelemtemp)
                print 'gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]]'
                print gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]]
                gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] = min(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]], gdat.fittmaxmnumbelem[l][d])
                gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] = max(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]], gdat.fittminmnumbelem[l][d])
                #gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] = random_integers(gdat.fittminmnumbelem[l][d], gdat.fittmaxmnumbelem[l][d])
    
    if gdat.fittnumbtrap > 0:
        print 'gdatmodi.thissamp[gdat.fittindxfixpnumbelemtotl]'
        print gdatmodi.thissamp[gdat.fittindxfixpnumbelem]
        print

    ## impose user-specified initial state
    if gdat.inittype == 'reco':
        if gdat.namerecostat != None:
            strgcnfg = gdat.namerecostat
        else:
            strgcnfg = gdat.strgcnfg
        path = gdat.pathoutp + 'stat_' + strgcnfg + '.h5'
        if os.path.exists(path):
            boolinitreco = True
            thisfile = h5py.File(path, 'r')
            if gdat.verbtype > 0:
                print 'Initializing from the state %s...' % path
                print 'Likelihood:'
                print thisfile['lliktotl'][...]
                
                # find the number of populations provided
                maxmindxpopl = 0
                for l in range(10):
                    for attr in thisfile:
                        if attr.startswith('lgalpop'):
                            indxpopl = int(attr[7])
                            if indxpopl > maxmindxpopl:
                                maxmindxpopl = indxpopl
                numbpoplinpt = maxmindxpopl + 1
                
                if numbpoplinpt != gdat.fittnumbpopl:
                    raise Exception('Cannot initialize with different number of populations.')
                # find the number of elements provided
                cntr = zeros(numbpoplinpt, dtype=int)
                for attr in thisfile:
                    if attr.startswith('lgalpop'):
                        indxpopl = int(attr[7])
                        cntr[indxpopl] += 1
                if gdat.verbtype > 0:
                    print 'Number of elements found:'
                    print cntr

            print 'attributes in the file'
            for attr in thisfile:
                print attr
            for attr in thisfile:
                for k, namefixp in enumerate(gdat.fittnamefixp):
                    if namefixp == attr:
                        gdatmodi.thissamp[k] = cdfn_fixp(gdat, 'fitt', thisfile[attr][()], k)
                        if gdatmodi.thissamp[k] == 0.:
                            print 'Warning CDF is zero.'
                        if not isfinite(thisfile[attr][()]):
                            raise Exception('Retreived state parameter is not finite.')
                        if not isfinite(gdatmodi.thissamp[k]) or gdatmodi.thissamp[k] < 0. or gdatmodi.thissamp[k] > 1.:
                            print 'namefixp'
                            print namefixp
                            print 'thisfile[attr][()]'
                            print thisfile[attr][()]
                            print 'gdatmodi.thissamp[k]'
                            print gdatmodi.thissamp[k]
                            raise Exception('CDF of the retreived state parameter is bad.')
            if gdat.fittnumbtrap > 0:
                for l in gdat.fittindxpopl:
                    for d in gdat.fittindxregipopl[l]:
                        maxmnumbelem = getattr(gdat, 'fittmaxmnumbelempop%dreg%d' % (l, d))
                        if gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] > maxmnumbelem:
                            gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] = maxmnumbelem
                            if gdat.verbtype > 0:
                                print 'Tapering off the element list...'

            if gdat.verbtype > 0 and gdat.fittnumbtrap > 0:
                print 'gdatmodi.thissamp[gdat.fittindxfixpnumbelem]'
                print gdatmodi.thissamp[gdat.fittindxfixpnumbelem]
            
            retr_elemlist(gdat, gdatmodi)
            gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)
            
            if (gdatmodi.thissamp == 0).all():
                raise Exception('Bad initialization.')
    
            if gdat.fittnumbtrap > 0 and gdatmodi.thisindxsampcomp != None:
                for strgcomp in gdat.fittliststrgcomptotl:
                    initcomp = [[[] for d in gdat.indxregi] for l in gdat.fittindxpopl]
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            initcomp[l][d] = empty(len(gdatmodi.thisindxelemfull[l][d]))
                            for k in range(len(gdatmodi.thisindxelemfull[l][d])):
                                namefiel = '%spop%dreg%d%04d' % (strgcomp, l, d, k)
                                for attr in thisfile:
                                    if namefiel == attr:
                                        initcomp[l][d][k] = thisfile[namefiel][()]
                    setattr(gdat, 'init' + strgcomp, initcomp)
                initcompfromstat(gdat, gdatmodi, 'init')
            thisfile.close()
        else:
            boolinitreco = False
            if gdat.verbtype > 0:
                print 'Initialization from previous state, %s, failed.' % path
        
    if gdat.inittype == 'refr' or gdat.inittype == 'pert':
        for k, namefixp in enumerate(gdat.fittnamefixp):
            if not (gdat.inittype == 'pert' and namefixp.startswith('numbelem')) and namefixp in gdat.truenamefixp:
                indxfixptrue = where(gdat.truenamefixp == namefixp)[0]
                gdatmodi.thissamp[k] = cdfn_fixp(gdat, 'fitt', gdat.truesampvarb[indxfixptrue], k)
                gdatmodi.thissampvarb[k] = icdf_fixp(gdat, 'fitt', gdatmodi.thissamp[k], k)

        retr_elemlist(gdat, gdatmodi)
        if gdat.fittnumbtrap > 0:
            if gdat.inittype == 'refr':
                initcompfromstat(gdat, gdatmodi, 'true')
    
    ## impose user-specified individual initial values
    for k, namefixp in enumerate(gdat.fittnamefixp):
        if gdat.inittype == 'reco' or  gdat.inittype == 'refr' or gdat.inittype == 'pert':
            try:
                getattr(gdat, 'init' + namefixp)
                print 'Conflicting initial state arguments detected, init keyword takes precedence.'
            except:
                pass
        try:
            initvalu = getattr(gdat, 'init' + namefixp)
            gdatmodi.thissamp[k] = cdfn_fixp(gdat, 'fitt', initvalu, k)
            if gdat.verbtype > 0:
                print 'Received initial condition for %s: %.3g' % (namefixp, initvalu)
            print
        except:
            pass

    if gdat.inittype == 'rand' or gdat.inittype == 'reco' and not boolinitreco:
        if gdat.verbtype > 0:
            print 'Initializing from a random state...'
        retr_elemlist(gdat, gdatmodi)
        gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)

    # check the initial unit sample vector for bad entries
    if gdat.fittnumbtrap > 0:
        indxsampdiff = setdiff1d(gdat.fittindxpara, gdat.fittindxfixpnumbelemtotl)
        indxsampbaddlowr = where(gdatmodi.thissamp[indxsampdiff] <= 0.)[0]
        indxsampbadduppr = where(gdatmodi.thissamp[indxsampdiff] >= 1.)[0]
        indxsampbaddlowr = indxsampdiff[indxsampbaddlowr]
        indxsampbadduppr = indxsampdiff[indxsampbadduppr]
    else:
        indxsampbaddlowr = where(gdatmodi.thissamp <= 0.)[0]
        indxsampbadduppr = where(gdatmodi.thissamp >= 1.)[0]
    indxsampbadd = concatenate((indxsampbaddlowr, indxsampbadduppr))
    if indxsampbadd.size > 0:
        print 'Initial value caused unit sample vector to go outside the unit interval...'
        print 'gdat.fittnamepara[indxsampbadd]'
        print gdat.fittnamepara[indxsampbadd]
        print 'gdatmodi.thissamp[indxsampbadd, None]'
        print gdatmodi.thissamp[indxsampbadd, None]
        gdatmodi.thissamp[indxsampbadd] = rand(indxsampbadd.size)
        print 'Initializing these parameters from the prior...'
        #raise Exception('')
    
    gdatmodi.thissampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.thissamp, gdatmodi.thisindxsampcomp)
    
    if gdat.verbtype > 1:
        show_samp(gdat, gdatmodi, thisonly=True)
    
    ## sample index
    gdatmodi.cntrswep = 0
   
    if gdat.diagmode:
        if gdat.indxswepsave.size != gdat.numbsamp:
            raise Exception('Inappropriate number of samples.')

    # initialize the worker sampler
    ## prepare gdatmodi
    #gdatmodi.thislliktotl = 0.
    #gdatmodi.thissbrtdfnc = zeros_like(gdat.expo)
    #gdatmodi.thissbrthost = zeros_like(gdat.expo)
    ##gdatmodi.thisdeflhost = zeros_like(gdat.expo)
    ##gdatmodi.thispsfnconv = zeros_like(gdat.expo)
    #gdatmodi.thisllik = zeros_like(gdat.expo)
    ##prep_gdatmodi(gdat, gdatmodi, gdatmodi, 'this')
    #gdatmodi.thisstdvsamp = zeros(gdat.fittnumbpara)
    
    gdatmodi.thisaccp = zeros(1, dtype=bool)
    gdatmodi.thisaccppsfn = zeros(1, dtype=bool)
    gdatmodi.thisaccpprio = zeros(1, dtype=bool)
    gdatmodi.thisaccpprop = zeros(1, dtype=bool)
    gdatmodi.thisindxproptype = zeros(1, dtype=int)
    gdatmodi.thisauxipara = zeros(gdat.fittmaxmnumbcomp)
    gdatmodi.thislcomfact = zeros(1)
    gdatmodi.thisljcbfact = zeros(1)
    gdatmodi.thisaccpprob = zeros(1)
    gdatmodi.thischro = zeros(gdat.numbchro)
    gdatmodi.thisdeltlliktotl = zeros(1)
    gdatmodi.thislpau = zeros(gdat.fittnumblpau)
    gdatmodi.thismemoresi = zeros(1)
    gdatmodi.thislfctasym = zeros(1)
    gdatmodi.thislpriprop = zeros(gdat.fittnumblpri)
    
    # dummy definitions required for logs
    gdatmodi.nextdeltlliktotl = 0.

    # log the initial state
    if gdat.verbtype > 1:
        tdpy.util.show_memo(gdatmodi, 'gdatmodi')

    # process the initial sample, define the variables to be processed in each sample
    proc_samp(gdat, gdatmodi, 'this', 'fitt')
    
    # dummy definitions
    #if ('lght' in gdat.fittelemtype or 'clus' in gdat.fittelemtype or 'lghtline' in gdat.fittelemtype) and gdat.fittnumbtrap > 0:
    #    gdatmodi.nextsbrtdfnc = zeros_like(gdatmodi.thissbrtdfnc)
    
    # enter interactive mode
    if gdat.intrevalcntpmodl:
        plot_genemaps(gdat, gdatmodi, 'this', 'cntpmodl', strgcbar='cntpdata', thisindxener=0, thisindxevtt=0, intreval=True)
    if gdat.intrevalcntpresi:
        plot_genemaps(gdat, gdatmodi, 'this', 'cntpresi', thisindxener=0, thisindxevtt=0, intreval=True)

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

    if gdat.verbtype > 0:
        print 'Sampling...'

    if gdat.verbtype > 1:
        print 'gdat.stdvstdp'
        print gdat.stdvstdp[:, None]
        print
    
    gdatmodi.thisstdvstdp = copy(gdat.stdvstdp)

    gdatmodi.optidone = False 

    while gdatmodi.cntrswep < gdat.numbswep:
        
        gdatmodi.thischro[:] = 0.
        
        if gdat.diagmode:
            gdatmodi.thissampdiag = copy(gdatmodi.thissamp)
            gdatmodi.thissampvarbdiag = copy(gdatmodi.thissampvarb)

        if gdat.optitypetemp == 'hess':
            optihess(gdat, gdatmodi)
            gdatmodi.optidone = True
            
        if gdatmodi.optidone:
            thisfile = h5py.File(gdat.pathoutpthis + 'opti.h5', 'w')
            thisfile.create_dataset('stdvstdp', data=gdatmodi.stdvstdp)
            thisfile.close()
            break

        if gdat.emptsamp:
            print 'Empty sampling. Sample number %d' % gdatmodi.cntrswep
            break

        initchro(gdat, gdatmodi, 'next', 'totl')
        
        if gdat.verbtype > 1:
            print
            print '-' * 10
            print 'Sweep %d' % gdatmodi.cntrswep

        thismakefram = (gdatmodi.cntrswep % gdat.numbswepplot == 0) and gdatmodi.indxprocwork == int(float(gdatmodi.cntrswep) / gdat.numbswep * gdat.numbproc) \
                                                                                   and gdat.makeplotfram and gdat.makeplot
        # propose the next sample
        if gdat.verbtype > 1:        
            print
            print '-----'
            print 'Proposing...'
            print 'thislliktotl'
            print gdatmodi.thislliktotl
            print 'thislpritotl'
            print gdatmodi.thislpritotl
            print 'thislpostotl'
            print gdatmodi.thislliktotl + gdatmodi.thislpritotl
            print
        
        if gdat.burntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
            gdatmodi.thisfacttmpr = ((gdatmodi.cntrswep + 1.) / gdat.numbburntmpr)**4
            #gdatmodi.thistmprfactdeltllik = gdatmodi.thisfacttmpr
            gdatmodi.thistmprfactdeltllik = 1. # gdatmodi.thisfacttmpr
            gdatmodi.thistmprfactstdv = 1. / gdatmodi.thisfacttmpr
            #gdatmodi.thistmprlposelem = -1000. * (1. - gdatmodi.thisfacttmpr) * concatenate(gdatmodi.thisindxsampcomp['comp']).size
            gdatmodi.thistmprlposelem = 0.
        else:
            gdatmodi.thistmprfactdeltllik = 1.
            gdatmodi.thistmprfactstdv = 1.
            gdatmodi.thistmprlposelem = 0. 
        
        initchro(gdat, gdatmodi, 'next', 'prop')
        prop_stat(gdat, gdatmodi, 'fitt')
        stopchro(gdat, gdatmodi, 'next', 'prop')
       
        if gdat.optitypetemp == 'auto' and gdatmodi.cntrswep == 0 or gdat.evoltype == 'maxmllik':
            gdatmodi.thisstdpscalfact *= 1.5**gdatmodi.nextdeltlliktotl
        else:
            gdatmodi.thisstdpscalfact = 1.
            
        if gdat.verbtype > 1:
            show_samp(gdat, gdatmodi)
    
        # diagnostics
        if gdat.diagmode:
            
            initchro(gdat, gdatmodi, 'next', 'diag')
            
            indxsampbadd = where((gdatmodi.thissamp[gdat.fittnumbpopl*gdat.numbregi:] > 1.) | (gdatmodi.thissamp[gdat.fittnumbpopl*gdat.numbregi:] < 0.))[0] + 1
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
            
            if gdatmodi.propwith and gdatmodi.nextlpautotl != 0.:
                raise Exception('Auxiliary variable PDF should is not zero during a within-model proposal.')
            
            if not isfinite(gdatmodi.thislliktotl):
                raise Exception('Log-likelihood is infinite!')
    
            if gdatmodi.cntrswep == 0:
                gdatmodi.thislliktotlprev = gdatmodi.thislliktotl
            if gdatmodi.thislliktotl - gdatmodi.thislliktotlprev < -10.:
                print 'gdatmodi.thislliktotlprev'
                print gdatmodi.thislliktotlprev
                print 'gdatmodi.thislliktotl'
                print gdatmodi.thislliktotl
                #print 'loglikelihood drop is very unlikely!'
                raise Exception('loglikelihood drop is very unlikely!')
            gdatmodi.thislliktotlprev = gdatmodi.thislliktotl
       
            for strgstat in ['this', 'next']:
                for strgvarb in ['samp', 'sampvarb']:
                    varb = getattr(gdatmodi, strgstat + strgvarb)
                    if not isfinite(varb).all():
                        print 'gdatmodi' + strgstat + strgvarb
                        for k in gdat.fittindxpara:
                            print varb[k]
                        raise Exception('Sample vector is not finite.')
            
            if gdat.fittnumbtrap > 0:
                if gdat.fittboolelemsbrtdfncanyy:
                    frac = amin(gdatmodi.thissbrtdfnc) / mean(gdatmodi.thissbrtdfnc)
                    cntppntschec = retr_cntp(gdat, gdatmodi.thissbrtdfnc, gdat.indxregi, gdat.indxcubeeval)
                    for d in gdat.indxregi:
                        if amin(cntppntschec[d]) < -0.1 and frac < -1e-3:
                            raise Exception('thissbrtdfnc went negative by %.3g percent.' % (100. * frac))
                    
            # check the population index
            try:
                if gdatmodi.indxpoplmodi < 0:
                    raise Exception('Bad population index')
            except:
                pass
            
            for d in gdat.indxregi:
                if (gdatmodi.thiscntpmodl[d] <= 0.).any() or not (isfinite(gdatmodi.thiscntpmodl[d])).all():
                    raise Exception('Current flux model is not positive')
            if gdatmodi.cntrswep > 0:
                if isnan(gdatmodi.thislpri).any():
                    raise Exception('Delta log-prior is not finite.')
                if not isfinite(gdatmodi.nextdeltlliktotl):
                    raise Exception('Delta log-likelihood is not finite.')

            if gdat.fittnumbtrap > 0:
                for l in gdat.fittindxpopl:
                    for d in gdat.fittindxregipopl[l]:
                        if gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l][d]] != len(gdatmodi.thisindxelemfull[l][d]):
                            print 'gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem]'
                            print gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem]
                            print 'gdatmodi.thisindxelemfull'
                            print gdatmodi.thisindxelemfull
                            raise Exception('Number of elements is inconsistent with the element index list.')

                        if gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[l][d]] != len(gdatmodi.thisindxelemfull[l][d]):
                            raise Exception('Number of elements is inconsistent across data structures.')
                        
                        for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                            if gdat.fittlistscalcomp[l][k] == 'gaus' or gdat.fittlistscalcomp[l][k] == 'igam':
                                continue
                            comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l][d]]
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
        
            stopchro(gdat, gdatmodi, 'next', 'diag')
    
        # save the sample
        if gdat.boolsave[gdatmodi.cntrswep]:
           
            initchro(gdat, gdatmodi, 'next', 'save')
        
            if gdat.savestat:
                
                if gdat.namesavestat != None:
                    strgcnfg = gdat.namesavestat
                else:
                    strgcnfg = gdat.strgcnfg
                path = gdat.pathoutp + 'stat_' + strgcnfg + '.h5'
                
                booltemp = False
                if os.path.isfile(path):
                    thisfilechec = h5py.File(path, 'r')
                    if thisfilechec['lliktotl'][...] > gdatmodi.thislliktotl:
                        if gdat.verbtype > 0:
                            print 'Not saving the state to %s because loglikelihood is lower...' % path
                            print 'Likelihood in the file:'
                            print thisfilechec['lliktotl'][...]
                    else:
                        booltemp = True
                    thisfilechec.close()
                else:
                    booltemp = True
                if gdat.forcsavestat:
                    booltemp = True
                if booltemp:
                    if gdatmodi.indxprocwork > 0:
                        continue
                    if gdat.verbtype > 0:
                        print 'Saving the state to %s...' % path
        
                    thisfile = h5py.File(path, 'w')
                    thisfile.create_dataset('lliktotl', data=gdatmodi.thislliktotl)
                    for namefixp in gdat.fittnamefixp:
                        indxfixp = getattr(gdat, 'fittindxfixp' + namefixp)
                        valu = gdatmodi.thissampvarb[indxfixp]
                        thisfile.create_dataset(namefixp, data=valu)
                    if gdat.fittnumbtrap > 0:
                        for l in gdat.fittindxpopl:
                            for strgcomp in gdat.fittliststrgcomp[l]:
                                for d in gdat.fittindxregipopl[l]:
                                    comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l][d]]
                                    for k in arange(comp.size):
                                        name = strgcomp + 'pop%dreg%d%04d' % (l, d, k)
                                        thisfile.create_dataset(name, data=comp[k])
                    thisfile.close()
            
            # preprocess the current sample to calculate variables that are not updated
            proc_samp(gdat, gdatmodi, 'this', 'fitt')
            
            indxsampsave = gdat.indxsampsave[gdatmodi.cntrswep]
            
            # fill the sample lists
            for strgvarb in gdat.liststrgvarbarrysamp:
                valu = getattr(gdatmodi, 'this' + strgvarb)
                workdict['list' + strgvarb][indxsampsave, ...] = valu
            for strgvarb in gdat.liststrgvarblistsamp:
                workdict['list' + strgvarb].append(deepcopy(getattr(gdatmodi, 'this' + strgvarb)))
            stopchro(gdat, gdatmodi, 'next', 'save')

        # plot the current sample
        if thismakefram:
            
            initchro(gdat, gdatmodi, 'next', 'plot')
            
            if gdat.verbtype > 0:
                print 'Process %d is in queue for making a frame.' % gdatmodi.indxprocwork
            
            proc_samp(gdat, gdatmodi, 'this', 'fitt')
            
            if gdat.numbproc > 1:
                gdatmodi.lock.acquire()
            
            if gdat.verbtype > 0:
                print 'Process %d started making a frame.' % gdatmodi.indxprocwork
            
            plot_samp(gdat, gdatmodi, 'this', 'fitt')
            
            if gdat.verbtype > 0:
                print 'Process %d finished making a frame.' % gdatmodi.indxprocwork
        
            if gdat.numbproc > 1:
                gdatmodi.lock.release()
        
            stopchro(gdat, gdatmodi, 'next', 'plot')
    
        # temp
        if gdat.fittpsfnevaltype != 'none':
            if gdat.fittpsfntype == 'doubking':
                if gdatmodi.nextsampvarb[gdat.fittindxfixppsfp[1]] >= gdatmodi.nextsampvarb[gdat.fittindxfixppsfp[3]]:
                    for k in range(20):
                        print 'Proposal rejected due to PSF'
                    gdatmodi.thisaccppsfn = False
                    print 'gdatmodi.nextsampvarb'
                    print gdatmodi.nextsampvarb
                    print 'gdatmodi.nextsampvarb[gdat.fittindxfixppsfp]'
                    print gdatmodi.nextsampvarb[gdat.fittindxfixppsfp]
                    print 'gdatmodi.propbrth'
                    print gdatmodi.propbrth
       
        # assertions
        if gdat.diagmode:
            if not gdat.calcllik and gdatmodi.propllik:
                raise Exception('')
                
            if gdat.fittnumbtrap == 0 and gdatmodi.propelem:
                raise Exception('')
    
        # determine the acceptance probability
        gdatmodi.thisaccpprop = gdatmodi.thisaccpprio and gdatmodi.thisaccppsfn
        if gdatmodi.thisaccpprop:
            
            proc_samp(gdat, gdatmodi, 'next', 'fitt')
            if not gdat.calcllik:
                gdatmodi.nextdeltlliktotl = 0.
        
            if gdat.verbtype > 1:
                print 'gdatmodi.nextdeltlliktotl'
                print gdatmodi.nextdeltlliktotl
                print 'gdatmodi.nextdeltlpritotl'
                print gdatmodi.nextdeltlpritotl
                print 'gdatmodi.nextlpautotl'
                print gdatmodi.nextlpautotl
                print 'gdatmodi.thislfctasym'
                print gdatmodi.thislfctasym
                print 'gdatmodi.thislcomfact'
                print gdatmodi.thislcomfact
                print 'gdatmodi.thisljcbfact'
                print gdatmodi.thisljcbfact
                print 'gdatmodi.thistmprlposelem'
                print gdatmodi.thistmprlposelem
                print
            
            if gdat.diagmode:
                if gdatmodi.nextdeltlliktotl == 0 and gdatmodi.nextdeltlpritotl == 0.:
                    #raise Exception('Both likelihood and prior did not change.')
                    print 'Both likelihood and prior did not change.'

            # evaluate the acceptance probability
            gdatmodi.thisaccpprob[0] = exp(gdatmodi.thistmprfactdeltllik * gdatmodi.nextdeltlliktotl + gdatmodi.thistmprlposelem + gdatmodi.nextdeltlpritotl)
            
        else:
            gdatmodi.thisaccpprob[0] = 0.
    
        # accept or reject the proposal
        if gdat.evoltype == 'maxmllik':
            booltemp = gdatmodi.thisaccpprop and gdatmodi.nextdeltlliktotl > 0.
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

            if False:
                print 'Rejected.'

            gdatmodi.thisaccp = False
        
        if False:
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

        ## variables to be saved for each sweep
        for strg in gdat.liststrgvarbarryswep:
            workdict['list' + strg][gdatmodi.cntrswep, ...] = getattr(gdatmodi, 'this' + strg)
        
        # save the execution time for the sweep
        stopchro(gdat, gdatmodi, 'next', 'totl')
        
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
                    
                    proc_samp(gdat, gdatmodi, 'this', 'fitt')
                    
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
                            binsampl = getattr(gdat, 'bins' + gdat.namefeatampl)
                            numbaccp = empty(gdat.numbbinsplot, dtype=int)
                            for a in gdat.indxbinsplot: 
                                boolbins = (binsampl[a] < workdict['listamplpert'][indxswepintv, 0]) & (workdict['listamplpert'][indxswepintv, 0]< binsampl[a+1])
                                numbaccp[a] = where(boolaccp & boolproptype & boolbins)[0].size
                            numbtotl = empty(gdat.numbbinsplot, dtype=int)
                            for a in gdat.indxbinsplot: 
                                boolbins = (binsampl[a] < workdict['listamplpert'][indxswepintv, 0]) & (workdict['listamplpert'][indxswepintv, 0]< binsampl[a+1])
                                numbtotl[a] = where(boolproptype & boolbins)[0].size
                            percaccp = zeros(gdat.numbbinsplot)
                            indx = where(numbtotl > 0)[0]
                            if indx.size > 0:
                                percaccp[indx] = 100. / numbtotl[indx].astype(float) * numbaccp[indx]
                        else:
                            numbaccp = where(boolaccp & boolproptype)[0].size
                            numbtotl = where(boolproptype)[0].size
                            if numbtotl > 0:
                                percaccp = 100. * numbaccp / float(numbtotl)
                            else:
                                percaccp = 0.
                        
                        if k in gdat.indxstdp:
                            strgstdvstdp = '%.3g' % gdat.stdvstdp[k]
                        else:
                            strgstdvstdp = ''
                        
                        if gdat.showmoreaccp and gdat.indxproptype[k] in gdat.indxproptypecomp:
                            for a in gdat.indxbinsplot:
                                print '%30s %40s %10s' % ('%s-%02d' % (gdat.legdproptype[k], a), 'acceptance rate: %3d%% (%5d out of %5d)' % \
                                                                                                            (percaccp[a], numbaccp[a], numbtotl[a]), strgstdvstdp)
                        else:
                            print '%30s %40s %10s' % (gdat.legdproptype[k], 'acceptance rate: %3d%% (%5d out of %5d)' % (percaccp, numbaccp, numbtotl), strgstdvstdp)
                        
                    if gdat.burntmpr and gdatmodi.cntrswep < gdat.numbburntmpr:
                        print 'Tempered burn-in'
                        print 'gdatmodi.thisfacttmpr'
                        print gdatmodi.thisfacttmpr
                    print 
                    numbpara = gdat.fittnumbfixp
                    print 'gdat.fittnumbfixp'
                    print gdat.fittnumbfixp
                    if gdat.fittnumbtrap > 0:
                        for l in gdat.fittindxpopl:
                            for d in gdat.fittindxregipopl[l]:
                                numbpara += gdatmodi.thisindxsampcomp['comp'][l][d].size
                    print 'Current number of parameters:'
                    print numbpara
                    if numbpara * 4 > gdat.factthin:
                        print 'Warning! Thinning factor is not enough!'
                        print 'gdat.factthin'
                        print gdat.factthin
                    print 'gdatmodi.thislliktotl'
                    print gdatmodi.thislliktotl
                    print 'gdatmodi.thislpritotl'
                    print gdatmodi.thislpritotl
                    for attr, valu in gdatmodi.__dict__.iteritems():
                        if isinstance(valu, ndarray):
                            #print 'attr'
                            #print attr
                            #print 'valu'
                            #print valu.shape
                            #print
                            if 8 * valu.size * gdat.numbsamptotl > 1e9:
                                print 'Warning! %s has total length %d and size %s' % (attr, valu.size * gdat.numbsamptotl, \
                                                                                                tdpy.util.retr_strgmemo(8 * valu.size * gdat.numbsamptotl))
                    print 'Backgrounds'
                    print gdatmodi.thissamp[gdat.fittindxfixpbacp]
                    print gdatmodi.thissampvarb[gdat.fittindxfixpbacp]
                    if gdat.fittnumbtrap > 0:
                        print 'Number of elements:'
                        for l in gdat.fittindxpopl:
                            print gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[l]].astype(int)
                        print 'Mean number of elements:'
                        print gdatmodi.thissampvarb[gdat.fittindxfixpmeanelem]
                        for l in gdat.fittindxpopl:
                            print 'Log-slope of the amplitude parameter distribution, population %d:' % l
                            indxfixp = getattr(gdat, 'fittindxfixp' + gdat.fittnamefeatampl[l] + 'distsloppop%d' % l)
                            print gdatmodi.thissampvarb[indxfixp]
                        print 'Log-prior penalization term: '
                        print gdatmodi.thislpri[0]
                        if gdat.allwrefr and gdat.asscrefr:
                            print 'Completeness'
                            for q in gdat.indxrefr:
                                for l in gdat.fittindxpopl:
                                    for d in gdat.indxregi:
                                        if gdat.refrnumbelem[q][d] == 0:
                                            continue
                                        namevarb = 'ref%dpop%dreg%d' % (q, l, d)
                                        print 'Region %d, Reference %d, Population %d' % (d, q, l)
                                        print 'Total'
                                        print getattr(gdatmodi, 'thiscmpl' + namevarb)
                                        refrfeat = getattr(gdat, 'refr' + gdat.namefeatsignrefr)
                                        if len(refrfeat[q][d]) > 0:
                                            print 'Binned in significance feature'
                                            print getattr(gdatmodi, 'thiscmpl' + gdat.namefeatsignrefr + namevarb)
                                            print 
                            print 'False discovery rate'
                            for q in gdat.indxrefr:
                                for l in gdat.fittindxpopl:
                                    for d in gdat.indxregi:
                                        if gdat.refrnumbelem[q][d] == 0:
                                            continue
                                        namevarb = 'ref%dpop%dreg%d' % (q, l, d)
                                        print 'Region %d, Reference %d, Population %d' % (d, q, l)
                                        print 'Total'
                                        print getattr(gdatmodi, 'thisfdis' + namevarb)
                                        refrfeat = getattr(gdat, 'refr' + gdat.namefeatsignrefr)
                                        if len(refrfeat[q][d]) > 0:
                                            print 'Binned in significance feature'
                                            print getattr(gdatmodi, 'thisfdis' + gdat.namefeatsignrefr + namevarb)
    
                    print 'Residual RMS'
                    for d in gdat.indxregi:
                        thiscntpresi = getattr(gdatmodi, 'thiscntpresireg%d' % d)
                        print sqrt(mean(thiscntpresi**2))

                    print 'Chronometers: '
                    print 'chro'
                    for k in range(gdat.numbchro):
                        for name, valu in gdat.indxchro.iteritems():
                            if valu == k:
                                print '%s: %.3g msec' % (name, gdatmodi.thischro[k] * 1e3)
                                if name == 'llik' and gdat.numbpixl > 1 and gdat.fittnumbtrap > 0 and 'loclhash' in gdat.commelemspatevaltype:
                                    print '%.3g per pixel' % (gdatmodi.thischro[k] * 1e3 / amin(gdat.numbpixlprox))
                   
                    print 

        if gdat.diagmode:
            indx = where(abs(gdatmodi.thissampdiag - gdatmodi.thissamp) > 0.)[0]
            if gdat.verbtype > 1:
                print 'Unit parameters that changed:'
                print indx
                print gdat.fittnamepara[indx]
            indx = where(abs(gdatmodi.thissampvarbdiag - gdatmodi.thissampvarb) > 0.)[0]
            if gdat.verbtype > 1:
                print 'Parameters that changed:'
                print indx
                print gdat.fittnamepara[indx]
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
        if gdat.optitypetemp != 'auto':
            gdatmodi.cntrswep += 1
    
    for strgvarb in gdat.liststrgvarbarry + gdat.liststrgvarblistsamp:
        valu = workdict['list' + strgvarb]
        setattr(gdatmodi, 'list' + strgvarb, valu)

    gdatmodi.timereal = time.time() - timereal
    gdatmodi.timeproc = time.clock() - timeproc
    
    path = gdat.pathoutpthis + 'gdatmodi%04d' % gdatmodi.indxprocwork
    writfile(gdatmodi, path) 
    
