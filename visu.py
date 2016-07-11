# common imports
from __init__ import *

# internal functions
from util import *

def plot_post(pathpcat):
    
    hdun = pf.open(pathpcat)

    print 'Loading PCAT output file %s...' % pathpcat
    print pf.info(pathpcat)

    gdat = gdatstrt()
    
    gdat.verbtype = 1

    gdat.numbener = hdun[0].header['numbener']
    gdat.numbevtt = hdun[0].header['numbevtt']
    gdat.numbpopl = hdun[0].header['numbpopl']
    gdat.numbpsfipara = hdun[0].header['numbpsfipara']
    gdat.numbformpara = hdun[0].header['numbformpara']
    
    gdat.indxener = arange(gdat.numbener)
    gdat.indxevtt = arange(gdat.numbevtt)
    gdat.indxpopl = arange(gdat.numbpopl)
    
    gdat.numbsamp = hdun[0].header['numbsamp']
    gdat.numbburn = hdun[0].header['numbburn']
    gdat.numbswep = hdun[0].header['numbswep']
    gdat.factthin = hdun[0].header['factthin']
    gdat.numbpopl = hdun[0].header['numbpopl']
    gdat.numbproc = hdun[0].header['numbproc']

    timeatcr = hdun[0].header['timeatcr']
    
    gdat.maxmgang = hdun[0].header['maxmgang']
    gdat.lgalcntr = hdun[0].header['lgalcntr']
    gdat.bgalcntr = hdun[0].header['bgalcntr']
    gdat.numbflux = hdun[0].header['numbflux']

    gdat.minmlgal = hdun[0].header['minmlgal']
    gdat.maxmlgal = hdun[0].header['maxmlgal']
    gdat.minmbgal = hdun[0].header['minmbgal']
    gdat.maxmbgal = hdun[0].header['maxmbgal']
    
    gdat.minmflux = hdun[0].header['minmflux']
    gdat.maxmflux = hdun[0].header['maxmflux']
    gdat.minmsind = hdun[0].header['minmsind']
    gdat.maxmsind = hdun[0].header['maxmsind']

    gdat.datatype = hdun[0].header['datatype']
    gdat.regitype = hdun[0].header['regitype']
    gdat.psfntype = hdun[0].header['psfntype']
    gdat.exprtype = hdun[0].header['exprtype']
    gdat.pixltype = hdun[0].header['pixltype']
    gdat.fdfntype = hdun[0].header['fdfntype']

    if gdat.pixltype == 'heal':
        gdat.numbsideheal = hdun[0].header['numbsideheal']
    else:
        gdat.numbsidecart = hdun[0].header['numbsidecart']

    gdat.trueinfo = hdun[0].header['trueinfo']
    
    gdat.margsize = hdun[0].header['margsize']
    gdat.strgtime = hdun[0].header['strgtime']
    gdat.strgcnfg = hdun[0].header['strgcnfg']

    gdat.stdvfdfnnorm = hdun[0].header['stdvfdfnnorm']
    gdat.stdvfdfnslop = hdun[0].header['stdvfdfnslop']
    gdat.stdvpsfipara = hdun[0].header['stdvpsfipara']
    gdat.stdvback = hdun[0].header['stdvback']
    gdat.stdvlbhl = hdun[0].header['stdvlbhl']
    gdat.stdvflux = hdun[0].header['stdvflux']
    gdat.radispmrlbhl = hdun[0].header['radispmrlbhl']
    gdat.fracrand = hdun[0].header['fracrand']
    
    gdat.specfraceval = hdun[0].header['specfraceval']

    if gdat.trueinfo and gdat.datatype == 'mock':
        gdat.mockpsfntype = hdun[0].header['mockpsfntype']

    gdat.strgexpr = hdun[0].header['strgexpr']
    gdat.strgback = []
    gdat.lablback = []
    gdat.nameback = []
    k = 0
    while True:
        try:
            gdat.strgback.append(hdun[0].header['strgback%04d' % k])
            gdat.lablback.append(hdun[0].header['lablback%04d' % k])
            gdat.nameback.append(hdun[0].header['nameback%04d' % k])
            k += 1
        except:
            break
    gdat.strgexpo = hdun[0].header['strgexpo']

    gdat.levi = hdun[0].header['levi']
    gdat.info = hdun[0].header['info']
    
    gdat.pathdata = hdun[0].header['pathdata']

    gdat.indxenerincl = hdun['indxenerincl'].data
    gdat.indxevttincl = hdun['indxevttincl'].data

    gdat.maxmnumbpnts = hdun['maxmnumbpnts'].data
        
    gdat.listlgalhist = hdun['lgalhist'].data
    gdat.listbgalhist = hdun['bgalhist'].data
    gdat.listganghist = hdun['ganghist'].data
    gdat.listaanghist = hdun['aanghist'].data
    gdat.listspechist = hdun['spechist'].data
    gdat.listsindhist = hdun['sindhist'].data
    
    pntsprob = hdun['pntsprob'].data
    
    listllik = hdun['llik'].data
    listlpri = hdun['lpri'].data
    
    atcr = hdun['atcr'].data
    gmrbstat = hdun['gmrbstat'].data
    
    listmodlcnts = hdun['modlcnts'].data
    
    if gdat.trueinfo and gdat.datatype == 'mock':
        gdat.mockspatdist = []
        for l in gdat.indxpopl:
            gdat.mockspatdist.append(hdun[0].header['mockspatdistpop%d' % l])
        gdat.mockfdfntype = hdun[0].header['mockfdfntype']
        if gdat.mockfdfntype == 'powr':
            gdat.mockfdfnslop = hdun['mockfdfnslop'].data
        if gdat.mockfdfntype == 'brok':
            gdat.mockfdfnsloplowr = hdun['mockfdfnsloplowr'].data
            gdat.mockfdfnslopuppr = hdun['mockfdfnslopuppr'].data
            gdat.mockfdfnbrek = hdun['mockfdfnbrek'].data
        gdat.mocknormback = hdun['mocknormback'].data
        gdat.mocknumbpnts = hdun['mocknumbpnts'].data

    # prior boundaries
    gdat.minmnormback = hdun['minmnormback'].data
    gdat.maxmnormback = hdun['maxmnormback'].data

    gdat.meansdfn = hdun['meansdfn'].data
    gdat.stdvsdfn = hdun['stdvsdfn'].data
        
    # bins
    gdat.binsener = hdun['binsener'].data
    gdat.indxenerfdfn = array([hdun[0].header['indxenerfdfn']])
   
    # hyperprior limits
    gdat.minmfdfnnorm = hdun['minmfdfnnorm'].data 
    gdat.maxmfdfnnorm = hdun['maxmfdfnnorm'].data
    if gdat.fdfntype == 'powr':
        gdat.minmfdfnslop = hdun['minmfdfnslop'].data
        gdat.maxmfdfnslop = hdun['maxmfdfnslop'].data
    if gdat.fdfntype == 'brok':
        gdat.minmfdfnbrek = hdun['minmfdfnbrek'].data
        gdat.maxmfdfnbrek = hdun['maxmfdfnbrek'].data
        gdat.minmfdfnsloplowr = hdun['minmfdfnsloplowr'].data
        gdat.maxmfdfnsloplowr = hdun['maxmfdfnsloplowr'].data
        gdat.minmfdfnslopuppr = hdun['minmfdfnslopuppr'].data
        gdat.maxmfdfnslopuppr = hdun['maxmfdfnslopuppr'].data
        
    # utilities
    gdat.listsamp = hdun['listsamp'].data
    gdat.probprop = hdun['probprop'].data
    gdat.listindxprop = hdun['indxprop'].data
    gdat.listchrototl = hdun['listchrototl'].data
    gdat.listchrollik = hdun['listchrollik'].data
    gdat.listaccp = hdun['accp'].data
    gdat.listindxparamodi = hdun['indxparamodi'].data
    gdat.listauxipara = hdun['auxipara'].data
    gdat.listlaccfrac = hdun['laccfrac'].data
    gdat.listnumbpair = hdun['numbpair'].data
    gdat.listjcbnfact = hdun['jcbnfact'].data
    gdat.listcombfact = hdun['combfact'].data
    gdat.listpntsfluxmean = hdun['listpntsfluxmean'].data

    # posterior distributions
    gdat.listnumbpnts = hdun['numbpnts'].data
    gdat.listfdfnnorm = hdun['fdfnnorm'].data
    if gdat.fdfntype == 'powr':
        gdat.listfdfnslop = hdun['fdfnslop'].data
    if gdat.fdfntype == 'brok':
        gdat.listfdfnsloplowr = hdun['fdfnsloplowr'].data
        gdat.listfdfnslopuppr = hdun['fdfnslopuppr'].data
        gdat.listfdfnbrek = hdun['fdfnbrek'].data
    gdat.listpsfipara = hdun['psfipara'].data
    gdat.listnormback = hdun['normback'].data

    gdat.makeplot = True

    # setup the sampler
    setp(gdat) 

    # truth gdat.information
    if gdat.datatype == 'mock':
        gdat.truenumbpnts = hdun['mocknumbpnts'].data
        gdat.truelgal = []
        gdat.truebgal = []
        gdat.truespec = []
        gdat.truesind = []
        for l in gdat.indxpopl:
            gdat.truelgal.append(hdun['mocklgalpop%d' % l].data)
            gdat.truebgal.append(hdun['mockbgalpop%d' % l].data)
            gdat.truespec.append(hdun['mockspecpop%d' % l].data)
            gdat.truesind.append(hdun['mocksindpop%d' % l].data)
        if gdat.mockfdfntype == 'powr':
            gdat.mockfdfnslop = hdun['mockfdfnslop'].data
        if gdat.mockfdfntype == 'brok':
            gdat.mockfdfnbrek = hdun['mockfdfnbrek'].data
            gdat.mockfdfnsloplowr = hdun['mockfdfnsloplowr'].data
            gdat.mockfdfnslopuppr = hdun['mockfdfnslopuppr'].data
        gdat.truenormback = hdun['mocknormback'].data
        gdat.truepsfipara = hdun['mockpsfipara'].data
        
    gdat.listlgal = [[] for l in gdat.indxpopl]
    gdat.listbgal = [[] for l in gdat.indxpopl]
    gdat.listspec = [[] for l in gdat.indxpopl]
    gdat.listsind = [[] for l in gdat.indxpopl]
    gdat.listgang = [[] for l in gdat.indxpopl]
    gdat.listaang = [[] for l in gdat.indxpopl]
    for l in gdat.indxpopl:
        lgal = hdun['lgalpop%d' % l].data
        bgal = hdun['bgalpop%d' % l].data
        spec = hdun['specpop%d' % l].data
        sind = hdun['sindpop%d' % l].data
        gang = hdun['gangpop%d' % l].data
        aang = hdun['aangpop%d' % l].data
        for j in gdat.indxsamptotl:
            indxpnts = where(spec[j, gdat.indxenerfdfn[0], :] > 0.)[0]
            numbpnts = indxpnts.size
            gdat.listlgal[l].append(lgal[j, indxpnts])
            gdat.listbgal[l].append(bgal[j, indxpnts])
            # temp
            gdat.listspec[l].append(spec[j, :, indxpnts].T)
            gdat.listsind[l].append(sind[j, indxpnts])
            gdat.listgang[l].append(gang[j, indxpnts])
            gdat.listaang[l].append(aang[j, indxpnts])

    # indices of the parameters to be plotted
    numbparaplot = min(gdat.numbpara, 50)
    size = numbparaplot - gdat.indxsampcompinit
    indxparaplot = concatenate([arange(gdat.indxsampcompinit), sort(choice(arange(gdat.indxsampcompinit, gdat.numbpara), size=size, replace=False))])

    # Gelman-Rubin test
    if gdat.numbproc > 1 and isfinite(gmrbstat).all():
        print 'Making the Gelman-Rubin TS plot...'
        tim0 = time.time()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        bins = linspace(1., amax(gmrbstat), 40)
        axis.hist(gmrbstat, bins=bins)
        axis.set_xlabel('PSRF')
        axis.set_ylabel('$N_{pix}$')
        plt.tight_layout()
        figr.savefig(gdat.pathplot + 'gmrbdist_' + gdat.rtag + '.pdf')
        plt.close(figr)
        
        tim1 = time.time()
        print 'Done in %.3g seconds.' % (tim1 - tim0)

    print 'Calculating the autocorrelation of the chain...'
    tim0 = time.time()
   
    # temp
    if False:
        print listllik[:, 0:8]
        numbproc = listllik.shape[0]
        listllik = listllik.flatten()
        levi = retr_levi(listllik)
        print levi
        print retr_info(listllik, levi)
        hist, bins = histogram(listllik)

        print 'hist'
        print hist
        print 'bins'
        print bins
        binsmaxm = bins[argmax(hist)]
        indxsamptemp = where((listllik > binsmaxm - 100) & (listllik < binsmaxm + 100))[0]
        listllik = listllik[indxsamptemp]
        levi = retr_levi(listllik)
        print levi
        print retr_info(listllik, levi)
        path = gdat.pathplot + 'llik_' + gdat.rtag
        print
        print
        tdpy.mcmc.plot_trac(path, listllik.flatten(), '$P(D|x)$')
        return

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    numbsampatcr = atcr.size
    
    axis.plot(arange(numbsampatcr), atcr)
    axis.set_xlabel('Sample index')
    axis.set_ylabel(r'$\tilde{\eta}$')
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'atcr_' + gdat.rtag + '.pdf')
    plt.close(figr)

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)

    print 'Calculating proposal and acceptance rates...'
    tim0 = time.time()

    # plot proposal efficiency
    numbtimemcmc = 20
    binstimemcmc = linspace(0., gdat.numbswep, numbtimemcmc)

    numbtick = 2

    figr, axgr = plt.subplots(gdat.numbprop, 1, figsize=(gdat.plotsize, gdat.numbprop * gdat.plotsize / 4.), sharex='all')
    for n, axis in enumerate(axgr):
        hist = axis.hist(where(gdat.listindxprop == n)[0], bins=binstimemcmc)[0]
        axis.hist(where((gdat.listindxprop == n) & (gdat.listaccp == True))[0], bins=binstimemcmc)
        axis.set_ylabel('%s' % gdat.strgprop[n])
        if n == gdat.numbprop - 1:
            axis.set_xlabel('$i_{samp}$')
        
        # define the y-axis
        maxm = amax(hist)
        axis.set_ylim([0., maxm])
        listtick = linspace(maxm / 2., maxm, numbtick)
        listlabltick = ['%.3g' % tick for tick in listtick]
        axis.set_yticks(listtick)
        axis.set_yticklabels(listlabltick)
    
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'propeffiprop_' + gdat.rtag + '.pdf')
    plt.close(figr)
    
    figr, axgr = plt.subplots(numbparaplot, 1, figsize=(gdat.plotsize, numbparaplot * gdat.plotsize / 4.), sharex='all')
    for n, axis in enumerate(axgr):
        hist = axis.hist(where(gdat.listindxparamodi == indxparaplot[n])[0], bins=binstimemcmc)[0]
        axis.hist(where((gdat.listindxparamodi == n) & (gdat.listaccp == True))[0], bins=binstimemcmc)
        axis.set_ylabel('$p_{%d}$' % indxparaplot[n])
        if n == gdat.numbprop - 1:
            axis.set_xlabel('$i_{samp}$')
        
        # define the y-axis
        maxm = amax(hist)
        axis.set_ylim([0., maxm])
        listtick = linspace(maxm / 2., maxm, numbtick)
        listlabltick = ['%.3g' % tick for tick in listtick]
        axis.set_yticks(listtick)
        axis.set_yticklabels(listlabltick)
    
    plt.subplots_adjust(hspace=0)
    figr.savefig(gdat.pathplot + 'propeffipara_' + gdat.rtag + '.pdf')
    plt.close(figr)
    
    # temp
    if False:
        print 'hey'
        print 'gdat.listindxprop'
        print gdat.listindxprop
    
    indxsampsplt = where(gdat.listindxprop == gdat.indxpropsplt)[0]
    indxsampmerg = where(gdat.listindxprop == gdat.indxpropmerg)[0]
            
    listname = ['laccfrac', 'numbpair', 'combfact', 'jcbnfact']
    listvarb = [gdat.listlaccfrac, gdat.listnumbpair, gdat.listcombfact, gdat.listjcbnfact]
    for k in range(4):
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
   
        # temp
        if False:
            print 'hey'
            print 'indxsampsplt'
            print indxsampsplt
            print 'indxsampmerg'
            print indxsampmerg
            print

        axis.hist(listvarb[k][indxsampsplt])
        axis.hist(listvarb[k][indxsampmerg])
        axis.set_ylabel('$N_{samp}$')
        axis.set_xlabel(listname[k])
        plt.tight_layout()
        figr.savefig(gdat.pathplot + listname[k] + gdat.rtag + '.pdf')
        plt.close(figr)
    
    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Calculating proposal execution times...'
    tim0 = time.time()

    try:
        plot_chro(gdat)
    except:
        print 'Time performance plots crashed'

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Making the grid posterior plot...'
    tim0 = time.time()
   
    # temp
    if False:
        path = gdat.pathplot + 'listsamp_' + gdat.rtag + '_'
        tdpy.mcmc.plot_grid(path, gdat.listsamp[:, indxparaplot], ['%d' % k for k in indxparaplot], numbplotside=10)

    numbpntspost = 3
    if gdat.trueinfo and gdat.datatype == 'mock' and gdat.mocknumbpnts[0] == numbpntspost and gdat.numbpopl == 1 and gdat.numbback == 1:
        numbpara = numbpntspost * gdat.numbcompcolr + gdat.numbener
        listpost = zeros((gdat.numbsamp, numbpara))
        for j in gdat.indxsamptotl:
            for k in range(numbpntspost):
                listpost[j, 0*numbpntspost+k] = gdat.listlgal[0][j][k]
                listpost[j, 1*numbpntspost+k] = gdat.listbgal[0][j][k]
                listpost[j, 2*numbpntspost+k] = gdat.listspec[0][j][gdat.indxenerfdfn[0], k]
                listpost[j, 3*numbpntspost+k] = gdat.listsind[0][j][k]
        for i in gdat.indxener:
            listpost[:, 4*numbpntspost+i] = gdat.listnormback[:, 0, i]
        truepost = zeros(numbpara)
        truepost[0*numbpntspost:1*numbpntspost] = gdat.truelgal[0][k]
        truepost[1*numbpntspost:2*numbpntspost] = gdat.truebgal[0][k]
        truepost[2*numbpntspost:3*numbpntspost] = gdat.truespec[0][0, gdat.indxenerfdfn[0], k]
        truepost[3*numbpntspost:4*numbpntspost] = gdat.truesind[0][k]
        truepost[4*numbpntspost:] = gdat.truenormback[0, :]
        path = gdat.pathplot + 'postdist_' + gdat.rtag
        strgpost = ['$%s_%d$' % (strg, indxpnts + 1) for strg in ['l', 'b', 'f', 's'] for indxpnts in arange(numbpnts)]
        strgpost += ['$A_{%d}$' % i for i in gdat.indxener]
        tdpy.mcmc.plot_grid(path, listpost, strgpost, truepara=truepost, numbtickbins=3)
           
    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Parsing the sample bundle and making frames...'
    tim0 = time.time()
    
    # flux match with the true catalog
    if gdat.trueinfo:
        for l in gdat.indxpopl:
            listindxmodl = []
            for p in gdat.indxsamptotl:
                indxmodl, indxtruepntsassc = corr_catl(gdat, l, gdat.listlgal[l][p], gdat.listbgal[l][p], gdat.listspec[l][p])
                listindxmodl.append(indxmodl)
            postspecmtch = zeros((3, gdat.numbener, gdat.truenumbpnts[l]))
            for i in gdat.indxener:
                gdat.listspecmtch = zeros((gdat.numbsamptotl, gdat.truenumbpnts[l]))
                for p in gdat.indxsamptotl:
                    indxpntstrue = where(listindxmodl[p] >= 0)[0]
                    gdat.listspecmtch[p, indxpntstrue] = gdat.listspec[l][p][i, listindxmodl[p][indxpntstrue]]
                postspecmtch[0, i, :] = percentile(gdat.listspecmtch, 16., axis=0)
                postspecmtch[1, i, :] = percentile(gdat.listspecmtch, 50., axis=0)
                postspecmtch[2, i, :] = percentile(gdat.listspecmtch, 84., axis=0)
            plot_scatspec(gdat, l, postspecmtch=postspecmtch)

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Making posterior distribution plots...'

    if gdat.pixltype == 'heal':
        reso = 60. / gdat.numbsideheal
        numbbinslgcr = int((gdat.maxmlgal - gdat.minmlgal) / reso)
        numbbinsbgcr = int((gdat.maxmbgal - gdat.minmbgal) / reso)
        pntsprobcart = zeros((numbbinslgcr, numbbinsbgcr, gdat.numbpopl, gdat.numbener, gdat.numbflux))
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for h in range(gdat.numbflux):
                    pntsprobcart[:, :, l, i, h] = tdpy.util.retr_cart(pntsprob[l, i, :, h], 
                                                                      indxpixlrofi=gdat.indxpixlrofi, \
                                                                      numbsideinpt=gdat.numbsideheal, \
                                                                      minmlgal=gdat.minmlgal, \
                                                                      maxmlgal=gdat.maxmlgal, \
                                                                      minmbgal=gdat.minmbgal, \
                                                                      maxmbgal=gdat.maxmbgal, \
                                                                      reso=reso)
    else:
        pntsprobcart = pntsprob.reshape((gdat.numbpopl, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbflux))
        pntsprobcart = swapaxes(swapaxes(pntsprobcart, 0, 2), 1, 3)
        
    # stacked posteiors binned in position and flux
    plot_pntsprob(gdat, pntsprobcart, ptag='quad')
    plot_pntsprob(gdat, pntsprobcart, ptag='full', full=True)
    plot_pntsprob(gdat, pntsprobcart, ptag='cumu', cumu=True)

    # compute median and 68% credible intervals of the inferred parameters
    gdat.postnormback = tdpy.util.retr_postvarb(gdat.listnormback)
    gdat.postpsfipara = tdpy.util.retr_postvarb(gdat.listpsfipara)
    
    # compute the medians of secondary variables
    gdat.medinormback = gdat.postnormback[0, :, :]
    gdat.medipsfipara = gdat.postpsfipara[0, :]
    gdat.medipsfn = retr_psfn(gdat, gdat.medipsfipara, gdat.indxener, gdat.binsangl, gdat.psfntype)
    gdat.medifwhm = 2. * retr_psfnwdth(gdat, gdat.medipsfn, 0.5)
    gdat.medibackfwhmcnts = retr_backfwhmcnts(gdat, gdat.medinormback, gdat.medifwhm)
    
    ## standard deviation axis
    gdat.medibinssigm = retr_sigm(gdat, gdat.binscnts, gdat.medibackfwhmcnts)

    # flux distribution
    for l in gdat.indxpopl:
        plot_histspec(gdat, l, listspechist=gdat.listspechist[:, l, :, :])

    # color distribution
    for l in gdat.indxpopl:
        plot_histsind(gdat, l, listsindhist=gdat.listsindhist[:, l, :])

    # fraction of emission components
    if gdat.numbback == 2:
        postpntsfluxmean = tdpy.util.retr_postvarb(gdat.listpntsfluxmean)
        plot_compfrac(gdat, postpntsfluxmean=postpntsfluxmean)

    # PSF parameters
    path = gdat.pathplot + 'psfipara_' + gdat.rtag
    if gdat.psfntype == 'singgaus' or gdat.psfntype == 'singking':
        gdat.listpsfipara[:, gdat.indxpsfiparainit] = rad2deg(gdat.listpsfipara[:, gdat.indxpsfiparainit])
        if gdat.trueinfo:
            gdat.truepsfipara[gdat.indxpsfiparainit] = rad2deg(gdat.truepsfipara[gdat.indxpsfiparainit])
    elif gdat.psfntype == 'doubgaus' or gdat.psfntype == 'gausking':
        gdat.listpsfipara[:, gdat.indxpsfiparainit+1] = rad2deg(gdat.listpsfipara[:, gdat.indxpsfiparainit+1])
        gdat.listpsfipara[:, gdat.indxpsfiparainit+2] = rad2deg(gdat.listpsfipara[:, gdat.indxpsfiparainit+2])
        if gdat.trueinfo:
            gdat.truepsfipara[gdat.indxpsfiparainit+1] = rad2deg(gdat.truepsfipara[gdat.indxpsfiparainit+1])
            gdat.truepsfipara[gdat.indxpsfiparainit+2] = rad2deg(gdat.truepsfipara[gdat.indxpsfiparainit+2])
    elif gdat.psfntype == 'doubking':
        gdat.listpsfipara[:, gdat.indxpsfiparainit+1] = rad2deg(gdat.listpsfipara[:, gdat.indxpsfiparainit+1])
        gdat.listpsfipara[:, gdat.indxpsfiparainit+3] = rad2deg(gdat.listpsfipara[:, gdat.indxpsfiparainit+3])
        if gdat.trueinfo:
            gdat.truepsfipara[gdat.indxpsfiparainit+1] = rad2deg(gdat.truepsfipara[gdat.indxpsfiparainit+1])
            gdat.truepsfipara[gdat.indxpsfiparainit+3] = rad2deg(gdat.truepsfipara[gdat.indxpsfiparainit+3])
    if gdat.trueinfo and gdat.psfntype == 'doubking':
        truepara = gdat.truepsfipara
    else:
        truepara = array([None] * gdat.numbpsfipara)
    tdpy.mcmc.plot_grid(path, gdat.listpsfipara, gdat.strgpsfipara, truepara=truepara, numbplotside=gdat.numbformpara, \
        numbbins=gdat.numbbins, numbtickbins=3)
    
    for k in range(gdat.numbpsfipara):
        path = gdat.pathplot + 'psfipara%d_' % k + gdat.rtag + '_'
        tdpy.mcmc.plot_trac(path, gdat.listpsfipara[:, k], gdat.strgpsfipara[k])
    
    # log-likelihood
    path = gdat.pathplot + 'llik_' + gdat.rtag
    tdpy.mcmc.plot_trac(path, listllik.flatten(), '$P(D|x)$')

    # log-prior
    path = gdat.pathplot + 'lpri_' + gdat.rtag
    tdpy.mcmc.plot_trac(path, listlpri.flatten(), '$P(x)$')

    # number, expected number of PS and flux conditional prior power law index 
    for l in range(gdat.numbpopl):
        
        # number of point sources
        path = gdat.pathplot + 'numbpntsdist_pop%d_' % l + gdat.rtag
        if gdat.trueinfo and gdat.truenumbpnts != None:
            truepara = gdat.truenumbpnts[l]
        else:
            truepara = None
        tdpy.mcmc.plot_trac(path, gdat.listnumbpnts[:, l], '$N$', truepara=truepara)

        # mean number of point sources
        path = gdat.pathplot + 'fdfnnorm_pop%d_' % l + gdat.rtag
        if gdat.trueinfo and gdat.datatype == 'mock':
            if gdat.mockfdfntype == 'powr':
                truepara = gdat.mocknumbpnts[l]
            else:
                truepara = None
        else:
            truepara = None
        tdpy.mcmc.plot_trac(path, gdat.listfdfnnorm[:, l], '$\mu$', truepara=truepara)

        # flux distribution
        titl = gdat.binsenerstrg[gdat.indxenerfdfn[0]]
        if gdat.fdfntype == 'powr':
            # power law index
            path = gdat.pathplot + 'fdfnslop_pop%d_' % l + gdat.rtag
            # temp
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfdfntype == 'powr':
                    truepara = gdat.mockfdfnslop[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha$'
            tdpy.mcmc.plot_trac(path, gdat.listfdfnslop[:, l], labl, truepara=truepara)
        
        if gdat.fdfntype == 'brok':
            # break flux
            path = gdat.pathplot + 'fdfnbrek_pop%d_' % l + gdat.rtag
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfdfntype == 'brok':
                    truepara = gdat.mockfdfnbrek[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$f_b$'
            tdpy.mcmc.plot_trac(path, gdat.listfdfnbrek[:, l], labl, truepara=truepara, scalpara='logt')
        
            # lower power law index
            path = gdat.pathplot + 'fdfnsloplowr_pop%d_' % l + gdat.rtag
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfdfntype == 'brok':
                    truepara = gdat.mockfdfnsloplowr[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha_1$'
            tdpy.mcmc.plot_trac(path, gdat.listfdfnsloplowr[:, l], labl, truepara=truepara)
        
            # uppr power law index
            path = gdat.pathplot + 'fdfnslopuppr_pop%d_' % l + gdat.rtag
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfdfntype == 'brok':
                    truepara = gdat.mockfdfnslopuppr[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha_2$'
            tdpy.mcmc.plot_trac(path, gdat.listfdfnslopuppr[:, l], labl, truepara=truepara)
        
    # background normalization
    for i in gdat.indxener:
        for c in gdat.indxback:
            path = gdat.pathplot + gdat.nameback[c] + '%d_' % i + gdat.rtag
            if gdat.trueinfo and gdat.datatype == 'mock':
                truepara = gdat.truenormback[c, i]
            else:
                truepara = None
            titl = gdat.binsenerstrg[i]
            labl = gdat.lablback[c] + '$_{%d}$' % i
            tdpy.mcmc.plot_trac(path, gdat.listnormback[:, c, i], labl, truepara=truepara, titl=titl)

    # plot log-likelihood
    figr, axrw = plt.subplots(2, 1, figsize=(gdat.plotsize, 2 * gdat.plotsize))
    figr.suptitle(r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi))
    for k, axis in enumerate(axrw):
        if k == 0:
            if amin(listllik) != amax(listllik):
                axis.hist(listllik.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(D|x)$')
        else:
            if amin(listlpri) != amax(listlpri):
                axis.hist(listlpri.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(x)$')
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'leviinfo_' + gdat.rtag + '.pdf')
    plt.close(figr)

    #make_anim()

    tim1 = time.time()
    print 'Plots are produced in %.3g seconds.' % (tim1 - tim0)


def plot_chro(gdat):

    binstime = logspace(log10(amin(gdat.listchrototl[where(gdat.listchrototl > 0)] * 1e3)), log10(amax(gdat.listchrototl * 1e3)), 50)
    figr, axcl = plt.subplots(gdat.numbprop, 1, figsize=(2 * gdat.plotsize, gdat.numbprop * gdat.plotsize / 2.))
    for k in range(gdat.numbprop):
        indxswepchro = where((gdat.listindxprop == k) & (gdat.listchrototl[:, 0] > 0))[0]
        if indxswepchro.size > 0:
            axcl[k].hist(gdat.listchrototl[indxswepchro, 0] * 1e3, binstime, log=True, label=gdat.strgprop[k])
        axcl[k].set_xlim([amin(binstime), amax(binstime)])
        axcl[k].set_ylim([0.5, None])
        axcl[k].set_ylabel(gdat.strgprop[k])
        axcl[k].set_xscale('log')
    axcl[-1].set_xlabel('$t$ [ms]')
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'chroprop_' + gdat.rtag + '.pdf')
    plt.close(figr)

    labl = ['Total', 'Proposal', 'Prior', 'Likelihood']
    numblabl = len(labl)
    figr, axcl = plt.subplots(2, 1, figsize=(2 * gdat.plotsize, numblabl * gdat.plotsize / 2.))
    for k in range(1, 4):
        axcl[0].hist(gdat.listchrototl[where(gdat.listchrototl[:, k] > 0)[0], k] * 1e3, binstime, log=True, label=labl[k])

    indxswepchro = where(gdat.listchrototl[:, 0])[0]
    if indxswepchro.size > 0:
        axcl[1].hist(gdat.listchrototl[indxswepchro, 0] * 1e3, binstime, log=True, label=labl[0], color='black')
    axcl[1].set_title(r'$\langle t \rangle$ = %.3g ms' % mean(gdat.listchrototl[where(gdat.listchrototl[:, 0] > 0)[0], 0] * 1e3))
    axcl[0].set_xlim([amin(binstime), amax(binstime)])
    axcl[1].set_xlabel('$t$ [ms]')
    axcl[0].set_xscale('log')
    axcl[1].set_xscale('log')
    axcl[0].set_ylim([0.5, None])
    axcl[1].set_ylim([0.5, None])
    axcl[0].legend(loc=1)
    axcl[1].legend(loc=2)
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'chrototl_' + gdat.rtag + '.pdf')
    plt.close(figr)

    listlabl = ['Setup', 'Pixel', 'Mesh', 'PS Flux', 'Total Flux', 'Counts', 'Likelihood']
    numblabl = len(listlabl)
    figr, axcl = plt.subplots(gdat.numbchrollik, 1, figsize=(2 * gdat.plotsize, gdat.plotsize * numblabl / 2.))
    
    maxmchrollik = amax(gdat.listchrollik * 1e3)
    if maxmchrollik > 0.:
        minmchrollik = amin(gdat.listchrollik[where(gdat.listchrollik > 0)] * 1e3)
        if maxmchrollik != minmchrollik:
            binstime = logspace(log10(minmchrollik), log10(maxmchrollik), 50)
            for k in range(gdat.numbchrollik):
                axcl[k].hist(gdat.listchrollik[where(gdat.listchrollik[:, k] > 0)[0], k] * 1e3, binstime, log=True, label=listlabl[k])
                axcl[k].set_xlim([amin(binstime), amax(binstime)])
                axcl[k].set_ylim([0.5, None])
                axcl[k].set_ylabel(listlabl[k])
                axcl[k].set_xscale('log')
            axcl[-1].set_xlabel('$t$ [ms]')
            plt.tight_layout()
            figr.savefig(gdat.pathplot + 'chrollik_' + gdat.rtag + '.pdf')
            plt.close(figr)


def plot_compfrac(gdat, gdatmodi=None, postpntsfluxmean=None):
    
    if postpntsfluxmean != None:
        post = True
    else:
        post = False
        
    listlinestyl = ['-', '--', '-.', ':']
    listcolr = ['black', 'b', 'b', 'b']
    listlabl = ['Data', 'PS', 'Iso', 'FDM']

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    
    listydat = empty((gdat.numbback + 2, gdat.numbener))
    listyerr = zeros((2, gdat.numbback + 2, gdat.numbener))
    
    listydat[0, :] = gdat.datafluxmean
    if post:
        listydat[1, :] = postpntsfluxmean[0, :]
        listyerr[:, 1, :] = tdpy.util.retr_errrvarb(postpntsfluxmean)
        for c in gdat.indxback:
            listydat[c+2, :] = gdat.postnormback[0, c, :] * gdat.backfluxmean[c]
            listyerr[:, c+2, :] = tdpy.util.retr_errrvarb(gdat.postnormback[:, c, :]) * gdat.backfluxmean[c]
    else:
        listydat[1, :] = mean(sum(gdatmodi.thispntsflux * gdat.expo, 2) / sum(gdat.expo, 2), 1)
        for c in gdat.indxback:
            listydat[c+2, :] = gdatmodi.thissampvarb[gdat.indxsampnormback[c, :]] * gdat.backfluxmean[c]
    
    xdat = gdat.meanener
    for k in range(gdat.numbback + 2):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [GeV]')
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [GeV/cm$^2$/s/sr]')
    axis.legend()

    if post:
        path = gdat.pathplot + 'compfracspec_' + gdat.rtag + '.pdf'
    else:
        path = gdat.pathplot + 'compfracspec_' + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
   
    listlablbackflux = ['PS', 'Iso', 'FDM']
    listexpl = [0.1, 0, 0]

    listsize = zeros(gdat.numbback + 1)
    for k in range(gdat.numbback + 1):
        if gdat.numbener == 1:
            listsize[k] = gdat.diffener * listydat[k+1, :]
        else:
            listsize[k] = trapz(listydat[k+1, :], gdat.meanener)
    listsize *= 100. / sum(listsize)
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))

    axis.pie(listsize, explode=listexpl, labels=listlablbackflux, autopct='%1.1f%%')
    axis.axis('equal')

    if post:
        path = gdat.pathplot + 'compfrac_' + gdat.rtag + '.pdf'
    else:
        path = gdat.pathplot + 'compfrac_' + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    plt.subplots_adjust(top=0.8, bottom=0.2, left=0.2, right=0.8)
    plt.savefig(path)
    plt.close(figr)
     

def plot_histsind(gdat, l, gdatmodi=None, listsindhist=None):
    
    if listsindhist == None:
        post = False
    else:
        post = True
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if post:
        xdat = gdat.meansind
        postsindhist = tdpy.util.retr_postvarb(listsindhist)
        ydat = postsindhist[0, :]
        yerr = tdpy.util.retr_errrvarb(postsindhist)
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    else:
        axis.hist(gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[l]], gdat.binssind, alpha=0.5, color='b', log=True, label='Sample')
    if gdat.trueinfo:
        axis.hist(gdat.truesind[l], gdat.binssind, alpha=0.5, color='g', log=True, label=gdat.truelabl)
        if gdat.datatype == 'mock':
            axis.hist(gdat.exprsind, gdat.binssind, alpha=0.1, color='red', log=True, label='3FGL')
    axis.set_yscale('log')
    axis.set_xlabel('$s$')
    axis.set_xlim([gdat.minmsind, gdat.maxmsind])
    axis.set_ylabel('$N$')
    axis.set_ylim([0.5, None])
    axis.legend(loc=2)
    if post:
        path = gdat.pathplot + 'histsind_pop%d_' % l + gdat.rtag + '.pdf'
    else:
        path = gdat.pathplot + 'histsind_pop%d_' % l + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    

def plot_histspec(gdat, l, gdatmodi=None, plotspec=False, listspechist=None):
  
    if listspechist == None:
        post = False
    else:
        post = True
    
    if post:
        binssigm = gdat.medibinssigm
    else:
        binssigm = gdatmodi.binssigm

    if plotspec:
        numbcols = gdat.numbener
    else:
        numbcols = 1

    figr, axcl = plt.subplots(1, numbcols, figsize=(gdat.plotsize * numbcols, gdat.plotsize))
    if plotspec:
        indxenertemp = gdat.indxener
    else:
        indxenertemp = gdat.indxenerfdfn
        axcltemp = []
        for i in gdat.indxener:
            if i == gdat.indxenerfdfn:
                axcltemp.append(axcl)
            else:
                axcltemp.append(None)
        axcl = axcltemp
    
    for i in indxenertemp:
        axis = axcl[i]
        if post:
            xdat = gdat.meanspec[i, :]
            postspechist = tdpy.util.retr_postvarb(listspechist[:, :, i])
            ydat = postspechist[0, :]
            yerr = tdpy.util.retr_errrvarb(postspechist)
            axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        else:
            spec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]][i, :]
            axis.hist(spec, gdat.binsspec[i, :], alpha=0.5, color='b', log=True, label='Sample')
            # superimpose the current prior flux distribution
            if i == gdat.indxenerfdfn[0]:
                fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[l]]
                if gdat.fdfntype == 'powr':
                    fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[l]]  
                    fluxhistmodl = fdfnnorm * pdfn_flux_powr(gdat, gdat.meanflux, fdfnslop) * gdat.diffflux
                if gdat.fdfntype == 'brok':
                    fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[l]]  
                    fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[l]]  
                    fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[l]]  
                    fluxhistmodl = fdfnnorm * pdfn_flux_brok(gdat, gdat.meanflux, fdfnbrek, fdfnsloplowr, fdfnslopuppr) * gdat.diffflux
                axis.plot(gdat.meanspec[i, :], fluxhistmodl, ls='--', alpha=0.5, color='b')

        # add horizontal axes for counts and fluctuation significance
        axiscnts = axis.twiny()
        axissigm = axis.twiny()
        axiscnts.xaxis.set_ticks_position('bottom')
        axiscnts.xaxis.set_label_position('bottom')
        axiscnts.set_xlim([gdat.binscnts[i, 0], gdat.binscnts[i, -1]])
        axissigm.set_xlim([binssigm[i, 0], binssigm[i, -1]])
        axiscnts.set_xscale('log')
        axissigm.set_xscale('log')
        axissigm.set_xlabel(r'$\sigma$')
        axiscnts.set_xlabel('$C$')
        axiscnts.spines['bottom'].set_position(('axes', 1.))
        axissigm.spines['top'].set_position(('axes', 1.05))
                
        # superimpose the true catalog
        if gdat.trueinfo:
            truehist = axis.hist(gdat.truespec[l][0, i, :], gdat.binsspec[i, :], alpha=0.5, color='g', log=True, label=gdat.truelabl)
            if gdat.datatype == 'mock':
                if gdat.exprtype == 'ferm':
                    axis.hist(gdat.exprspec[0, i, :], gdat.binsspec[i, :], color='red', alpha=0.1, log=True, label='3FGL')
        axis.set_yscale('log')
        axis.set_xlabel('$f$ ' + gdat.strgfluxunit)
        axis.set_xscale('log')
        axis.text(0.75, 0.65, gdat.binsenerstrg[i], ha='center', va='center', transform=axis.transAxes)
        if gdat.trueinfo:
            axis.set_ylim([0.5, None])
        axis.set_xlim([gdat.minmspec[i], gdat.maxmspec[i]])
        if plotspec:
            if i == 0:
                axis.set_ylabel('$N$')
            if i == numbcols / 2:
                axis.legend()
        else:
            axis.set_ylabel('$N$')
            axis.legend(loc=7)
    if plotspec:
        strg = 'spec'
    else:
        strg = 'flux'
    if post:
        path = gdat.pathplot + 'hist%s_pop%d_' % (strg, l) + gdat.rtag + '.pdf'
    else:
        path = gdat.pathplot + 'hist%s_pop%d_' % (strg, l) + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    #plt.subplots_adjust(top=0.9)
    plt.savefig(path)
    plt.close(figr)
    

def plot_scatspec(gdat, l, gdatmodi=None, postspecmtch=None):
    
    if postspecmtch != None:
        post = True
    else:
        post = False

    figr, axrw = plt.subplots(1, gdat.numbener, figsize=(gdat.plotsize * gdat.numbener, gdat.plotsize))
    if gdat.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        xdat = gdat.truespec[l][0, i, :]
        xerr = tdpy.util.retr_errrvarb(gdat.truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
 
        labl = '$f_{samp}$ ' + gdat.strgfluxunit
        if post:
            yerr[0, :] = postspecmtch[0, i, :]
            ydat = postspecmtch[1, i, :]
            yerr[1, :] = postspecmtch[2, i, :]
            yerr = fabs(yerr - ydat)
        else:
            ydat = gdatmodi.thisspecmtch[i, :]

        indx = where(ydat > 0.)[0]
        if indx.size > 0:
            axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='black')
        
        # tag multiple associations
        if not post:
            indx = gdatmodi.indxtruepntsassc[l].mult
            if len(indx) > 0:
                axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='red')
    
        # superimpose the bias line
        fluxbias = retr_fluxbias(gdat, gdat.binsspec, i)
        axis.plot(gdat.binsspec[i, :], gdat.binsspec[i, :], ls='--', alpha=0.2, color='black')
        axis.plot(gdat.binsspec[i, :], fluxbias[0, :], ls='--', alpha=0.2, color='black')
        axis.plot(gdat.binsspec[i, :], fluxbias[1, :], ls='--', alpha=0.2, color='black')
        
        # temp
        #if gdat.indxtruepntstimevari[l].size > 0:
        #    axis.errorbar(xdat[gdat.indxtruepntstimevari[l]], ydat[gdat.indxtruepntstimevari[l]], ls='', yerr=yerr[:, gdat.indxtruepntstimevari[l]], \
        #        lw=1, marker='o', markersize=5, color='red')
   
        if gdat.datatype == 'mock':
            strg = 'true'
        else:
            if gdat.exprtype == 'ferm':
                strg = '3FGL'
        axis.set_xlabel('$f_{%s}$ ' % strg + gdat.strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_title(gdat.enerstrg[i])

        ylim = [gdat.minmspec[i], gdat.maxmspec[i]]
        
        axis.set_ylim(ylim)
        axis.set_xlim([gdat.minmspec[i], gdat.maxmspec[i]])
        axis.set_title(gdat.binsenerstrg[i])

    if postspecmtch != None:
        path = gdat.pathplot + 'scatspec%d_' % l + gdat.rtag + '.pdf'
    elif gdatmodi.thisspecmtch != None:
        path = gdat.pathplot + 'scatspec%d_' % l + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep

    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)


def plot_scatpixl(gdat, gdatmodi, l):
    
    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            slop, intc, coef, pval, stdv = sp.stats.linregress(gdat.datacnts[i, :, m], gdatmodi.thismodlcnts[i, :, m])
            axis.scatter(gdat.datacnts[i, :, m], gdatmodi.thismodlcnts[i, :, m], alpha=0.5)

            axislimt = [0., amax(gdat.datacnts[i, :, m]) * 1.5]
            axis.set_xlim(axislimt)
            axis.set_ylim(axislimt)
            
            axis.text(0.5, 0.9, '$m = %.3g, b = %.3g, r = %.3g$' % (slop, intc, coef), va='center', ha='center', transform=axis.transAxes, fontsize=16)
            if m == gdat.numbevtt - 1:
                axis.set_xlabel('Data Counts')
            if i == 0:
                labl = 'Model Counts'
                if gdat.exprtype == 'ferm':
                    labl += ', ' + gdat.evttstrg[m]
                axis.set_ylabel(labl)
            if m == 0:
                axis.set_title(gdat.enerstrg[i])
            
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'scatpixl%d_' % l + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    
    
def plot_look(gdat):

    indxpixlproxsize = zeros(gdat.numbpixl)
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbfluxprox):
        for j in gdat.indxpixl:
            indxpixlproxsize[j] = indxpixlprox[h][j].size
        binspixlsize = logspace(log10(amin(indxpixlproxsize)), log10(amax(indxpixlproxsize)), 100)
        axis.hist(indxpixlproxsize, binspixlsize, log=True, label='Flux bin %d' % h)
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'look.pdf')
    plt.close()
    
    
def plot_psfn_type():
    
    devi = linspace(0., 5., 100)
    y = zeros((x.size, 5))

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    singgaus = retr_singgaus(devi, 0.25)
    axis.plot(devi, singgaus, label='Single Gaussian')

    singking = retr_singking(devi, 0.25, 10.)
    axis.plot(devi, singking, label='Single King')

    doubgaus = retr_doubgaus(devi, 0.1, 0.25, 1.)
    axis.plot(devi, doubgaus, label='Double Gaussian')

    gausking = retr_gausking(devi, 0.1, 0.25, 1., 10.)
    axis.plot(devi, gausking, label='Gaussian + King')

    doubking = retr_doubking(devi, 0.1, 0.25, 10., 1., 5.)
    axis.plot(devi, doubking, label='Double King')

    axis.legend()
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylim([1e-3, None])
    plt.show()
    
def plot_minmfluxinfo(minmfluxarry, listinfo, listlevi):
    
    figr, axis = plt.subplots()
    ax_ = axis.twinx()
    axis.plot(minmfluxarry, listinfo, label='Relative entropy')
    axis.legend(bbox_to_anchor=[0.2, 1.08], loc=2)
    
    ax_.plot(minmfluxarry, listlevi, label='Log-evidence', color='g')
    ax_.legend(bbox_to_anchor=[0.8, 1.08])

    axis.set_ylabel('$D_{KL}$ [nats]')
    ax_.set_ylabel(r'$\log P(D)$ [nats]')
    axis.set_xlabel('$f_{min}$ [1/cm$^2$/s/GeV]')
    axis.set_xscale('log')
    plt.tight_layout()
    figr.savefig(os.environ["PCAT_DATA_PATH"] + '/imag/minmfluxinfo.pdf')
    plt.close(figr)
    
    
def plot_evidtest():
    
    minmgain = -1.
    maxmgain = 5.
    minmdevi = 0.
    maxmdevi = 5.
    gain = linspace(minmgain, maxmgain, 100)
    devi = linspace(minmdevi, maxmdevi, 100)

    evid = log(sqrt(1. + exp(2. * gain[None, :])) * exp(-devi[:, None]**2 / 2. / (1. + 1. / exp(2. * gain[None, :]))))
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    figr.suptitle('Log-Bayesian Evidence For Lower-Dimension Model', fontsize=18)
    imag = axis.imshow(evid, extent=[minmgain, maxmgain, minmdevi, maxmdevi], cmap='winter', origin='lower')
    cset1 = plt.contourf(gain, devi, evid, cmap='winter')
    axis.set_xlabel('Information gain')
    axis.set_ylabel('Goodness of fit')
    plt.colorbar(imag, ax=axis, fraction=0.03)

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'evidtest_' + gdat.rtag + '.pdf')
    plt.close(figr)
    
    
def plot_pntsprob(gdat, pntsprobcart, ptag, full=False, cumu=False):
    
    if cumu:
        numbrows = 1
        numbcols = 1
    else:
        numbcols = 2
        if full:
            numbrows = gdat.numbflux / 2
        else:
            numbrows = 2
        
    if gdat.exprtype == 'ferm':
        strgvarb = '$f$'
        strgunit = ' [cm$^2$ s GeV]'
    if gdat.exprtype == 'sdss':
        strgvarb = '$C$'
        strgunit = ' [counts]'
    titl = strgvarb + strgunit

    for l in gdat.indxpopl:
        figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
        if numbrows == 1:
            axgr = [axgr]            
        figr.suptitle(titl, fontsize=18)
        for a, axrw in enumerate(axgr):
            if numbcols == 1:
                axrw = [axrw]
            for b, axis in enumerate(axrw):
                h = a * 2 + b
                if full:
                    indxlowr = h
                    indxuppr = h + 1
                elif cumu:
                    indxlowr = 0
                    indxuppr = gdat.numbflux
                else:
                    if h < 3:
                        indxlowr = 2 * h
                        indxuppr = 2 * (h + 1)
                    else:
                        indxlowr = 2 * h
                        indxuppr = gdat.numbflux
                temp = sum(pntsprobcart[:, :, l, gdat.indxenerfdfn[0], indxlowr:indxuppr], 2)
                imag = axis.imshow(temp, origin='lower', cmap='BuPu', extent=gdat.exttrofi)#, norm=mpl.colors.LogNorm(vmin=0.01, vmax=1))
                plt.colorbar(imag, fraction=0.05, ax=axis)

                # superimpose true PS
                # temp
                if gdat.trueinfo:
                    indxpnts = where((gdat.binsspec[gdat.indxenerfdfn[0], indxlowr] < gdat.truespec[l][0, gdat.indxenerfdfn[0], :]) & \
                        (gdat.truespec[l][0, gdat.indxenerfdfn[0], :] < gdat.binsspec[gdat.indxenerfdfn[0], indxuppr]))[0]
                    mar1 = axis.scatter(gdat.truelgal[l][indxpnts], gdat.truebgal[l][indxpnts], s=100, alpha=0.5, marker='*', lw=2, color='g')
                axis.set_xlabel(gdat.longlabl)
                axis.set_ylabel(gdat.latilabl)
                axis.set_xlim([gdat.frambndrmarg, -gdat.frambndrmarg])
                axis.set_ylim([-gdat.frambndrmarg, gdat.frambndrmarg])
                axis.axvline(gdat.frambndr, ls='--', alpha=0.3, color='black')
                axis.axvline(-gdat.frambndr, ls='--', alpha=0.3, color='black')
                axis.axhline(gdat.frambndr, ls='--', alpha=0.3, color='black')
                axis.axhline(-gdat.frambndr, ls='--', alpha=0.3, color='black')
                axis.set_title(tdpy.util.mexp(gdat.binsspec[gdat.indxenerfdfn, indxlowr]) + ' $<$ ' + \
                    strgvarb + ' $<$ ' + tdpy.util.mexp(gdat.binsspec[gdat.indxenerfdfn, indxuppr]))
        plt.tight_layout()
        figr.savefig(gdat.pathplot + 'pntsbind' + ptag + '%d%d' % (l, gdat.indxenerincl[gdat.indxenerfdfn]) + '_' + gdat.rtag + '.pdf')
        plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.binsangl)

    figr, axgr = plt.subplots(1, 2, figsize=(2 * gdat.plotsize, gdat.plotsize))
    figr.suptitle('King Function', fontsize=20)
    for k, axis in enumerate(axgr):
        if k == 0:
            sigmlist = [0.25]
            gammlist = [1.01, 2.5, 10.]
            lloc = 3
        else:
            sigmlist = [0.1, 0.25, 1.]
            gammlist = [2.]
            lloc = 1
        for sigm in sigmlist:
            for gamm in gammlist:
                axis.plot(angl, retr_singking(angl, sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        axis.legend(loc=lloc)
        axis.set_yscale('log')
        axis.set_xlabel(r'$\theta$ ' + gdat.strganglunit)
        axis.set_xlabel(r'$\mathcal{K}(\theta)$')
        
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'king.pdf')
    plt.close(figr)
    
    
def plot_psfn(gdat, gdatmodi):
    
    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            axis.plot(gdat.binsanglplot, gdatmodi.thispsfn[i, :, m], label='Sample')
            if gdat.trueinfo:
                if gdat.exprtype == 'ferm':
                    labl = 'Fermi-LAT'
                if gdat.exprtype == 'sdss':
                    labl = 'SDSS'
                axis.plot(gdat.binsanglplot, gdat.truepsfn[i, :, m], label=labl, color='g', ls='--')
            axis.set_yscale('log')
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$\theta$ ' + gdat.strganglunit)
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
            if m == 0:
                axis.set_title(gdat.enerstrg[i])  
            if i == gdat.numbener - 1 and m == gdat.numbevtt - 1:
                axis.legend(loc=2)
            indxsamp = gdat.indxsamppsfipara[i*gdat.numbformpara+m*gdat.numbener*gdat.numbformpara]
    
            if gdat.psfntype == 'singgaus':
                strg = r'$\sigma = %.3g$ ' % gdatmodi.thissampvarb[indxsamp]
            elif gdat.psfntype == 'singking':
                strg = r'$\sigma = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp]) + gdat.strganglunit + '\n'
                strg += r'$\gamma = %.3g$' % gdatmodi.thissampvarb[indxsamp+1]
            elif gdat.psfntype == 'doubgaus':
                strg = r'$f = %.3g$' % gdatmodi.thissampvarb[indxsamp] + '\n'
                strg += r'$\sigma = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+1]) + gdat.strganglunit + '\n'
                strg += r'$\sigma = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+2]) + gdat.strganglunit
            elif gdat.psfntype == 'gausking':
                strg = r'$f_G = %.3g$' % gdatmodi.thissampvarb[indxsamp] + '\n'
                strg += r'$\sigma_G = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+1]) + gdat.strganglunit + '\n'
                strg += r'$\sigma_K = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+2]) + gdat.strganglunit + '\n'
                strg += r'$\gamma = %.3g$' % gdatmodi.thissampvarb[indxsamp+3]
            elif gdat.psfntype == 'doubking':
                strg = r'$f_c = %.3g$' % gdatmodi.thissampvarb[indxsamp] + '\n'
                strg += r'$\sigma_c = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+1]) + gdat.strganglunit + '\n'
                strg += r'$\gamma_c = %.3g$' % gdatmodi.thissampvarb[indxsamp+2] + '\n'
                strg += r'$\sigma_t = %.3g$ ' % rad2deg(gdatmodi.thissampvarb[indxsamp+3]) + gdat.strganglunit + '\n'
                strg += r'$\gamma_t = %.3g$' % gdatmodi.thissampvarb[indxsamp+4]
            axis.text(0.75, 0.75, strg, va='center', ha='center', transform=axis.transAxes, fontsize=18)
            
            if gdat.exprtype == 'ferm':
                axis.set_ylim([1e-3, 1e6])
            if gdat.exprtype == 'sdss':
                axis.set_ylim([1e4, 1e11])

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'psfnprof_' + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    
    
def plot_fwhm(gdat, gdatmodi):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))

    tranfwhm = transpose(gdatmodi.thisfwhm)
    imag = axis.imshow(rad2deg(tranfwhm), origin='lower', extent=[gdat.binsener[0], gdat.binsener[-1], 0, 4], cmap='BuPu', interpolation='none', alpha=0.4)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.set_xscale('log')
    axis.set_ylabel('PSF Class')
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_yticks([0.5, 1.5, 2.5, 3.5])
    axis.set_yticklabels(['0', '1', '2', '3'])
    axis.set_xticks(gdat.binsener)
    axis.set_xticklabels(['%.2g' % gdat.binsener[i] for i in gdat.indxener])
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            axis.text(gdat.meanener[i], gdat.indxevttincl[m] + 0.5, r'$%.3g^\circ$' % rad2deg(tranfwhm[m, i]), ha='center', va='center', fontsize=14)

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'fwhmcnts_' + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    
    
def plot_datacntshist(gdat):

    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            datacntstemp = gdat.datacnts[i, :, m]
            
            maxmdatacntstemp = amax(datacntstemp)
            minmdatacntstemp = amin(datacntstemp[where(datacntstemp > 0.)])
            gdat.binscntstemp = linspace(minmdatacntstemp, maxmdatacntstemp, 20)
            meancntstemp = (gdat.binscntstemp[1:] + gdat.binscntstemp[0:-1]) / 2.
            diffcntstemp = gdat.binscntstemp[1:] - gdat.binscntstemp[0:-1]
            
            datacntshist = axis.hist(datacntstemp, gdat.binscntstemp, color='b')[0]

            init = [meancntstemp[argmax(datacntshist)], 1.]
            
            axis.set_xlim([amin(gdat.binscntstemp), amax(gdat.binscntstemp)])
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            #axis.set_xscale('log')
            if m == 0:
                axis.set_title(gdat.binsenerstrg[i])
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
        
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'datacntshist' + gdat.rtag + '.pdf')
    plt.close(figr)
    
    
def plot_intr():
    
    with plt.xkcd():

        from matplotlib import patheffects
        mpl.rcParams['path.effects'] = [patheffects.withStroke(linewidth=0)]

        figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))

        catl = arange(80)
        probcatl = pss.pmf(catl, 30.) + 0.5 * pss.pmf(catl, 60.)
        axis.plot(catl, probcatl)
        axis.set_xticks([10, 30, 60])
        axis.set_xticklabels(["Crackpot's Catalog", "Best-fit catalog", "Not-so-best-fit catalog"])
        axis.set_yticks([])
        text = axis.set_title("Exploring the catalog space with Probabilistic cataloging")
        text.set_position([.5, 1.05])
        axis.set_xlabel('Catalog index')
        axis.set_ylabel("Probability")
        
        axis.tick_params(axis='x', colors='#B6E954')
        axis.tick_params(axis='y', colors='#B6E954')
        axis.spines['bottom'].set_color('#B6E954')
        axis.spines['top'].set_color('#B6E954') 
        axis.spines['right'].set_color('#B6E954')
        axis.spines['left'].set_color('#B6E954')
        axis.yaxis.label.set_color('#B6E954')
        axis.xaxis.label.set_color('#B6E954')
        axis.title.set_color('#B6E954')

        axis.set_axis_bgcolor('black')
        figr.set_facecolor('black')
        plt.tight_layout()
        plt.savefig(os.environ["PCAT_DATA_PATH"] + '/imag/talkintr.pdf', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_eval(gdat):

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for k in range(gdat.numbfluxprox + 1):
        if k == 0 or k == gdat.numbfluxprox:
            alph = 1.
            if k == 0:
                labl = 'Dimmest PS'
                colr = 'b'
            else:
                labl = 'Brightest PS'
                colr = 'g'
        else:
            alph = 0.2
            labl = None
            colr = 'black'
        axis.plot(gdat.binsanglplot, gdat.binsfluxprox[k] * gdat.truepsfn[0, :, 0], label=labl, color=colr, alpha=alph)
        axis.set_xlim([amin(gdat.binsanglplot), amax(gdat.binsanglplot)])
        if k > 0:
            axis.axvline(gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
    axis.set_yscale('log')
    axis.set_xlabel(r'$\theta$ ' + gdat.strganglunit)
    axis.set_ylabel('$f$ [1/cm$^2$/s/sr/GeV]')
    axis.axhline(gdat.specfraceval * amax(gdat.binsfluxprox[0] * gdat.truepsfn[0, :, 0]), color='red', ls=':', label='Flux floor')
    axis.legend(handlelength=.5, loc=3)
    
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'eval_' + gdat.rtag + '.pdf')
    plt.close(figr)


def plot_3fgl_thrs(gdat):

    path = os.environ["PCAT_DATA_PATH"] + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = pf.getdata(path, 0)

    bgalfgl3 = linspace(-90., 90., 481)
    lgalfgl3 = linspace(-180., 180., 960)

    bgalexpo = linspace(-90., 90., 400)
    lgalexpo = linspace(-180., 180., 800)

    #fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)
    fluxthrs = griddata([lgalfgl3, bgalfgl3], fluxthrs, [gdat.lgalheal, gdat.bgalheal])

    cntsthrs = fluxthrs * gdat.expo

    jbgal = where(abs(bgalexpo) < 10.)[0]
    jlgal = where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel(gdat.longlabl)
    axis.set_ylabel(gdat.latilabl)

    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'thrs_' + gdat.rtag + '.pdf')
    plt.close(figr)
    

def make_anim():

    listname = ['errrpnts0A', 'datacnts0A', 'resicnts0A', 'modlcnts0A', 'histspec', 'scatspec', 'psfnprof', 'compfrac0', 'compfracspec', 'scatpixl']
    
    for name in listname:
    
        strg = '%s*0.pdf' % name
        listfile = fnmatch.filter(os.listdir(gdat.pathplot), strg)[int(numbburn/numbswepplot):]
        
        print fnmatch.filter(os.listdir(gdat.pathplot), strg)
        print listfile

        nfile = len(listfile)
        jfile = choice(arange(nfile), replace=False, size=nfile)

        cmnd = 'convert -delay 20 '
        for k in range(nfile):
            cmnd += '%s ' % listfile[jfile[k]]
        cmnd += ' %s.gif' % name
        os.system(cmnd)

        
def plot_histcnts(gdat, l, gdatmodi=None):

    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            if gdat.trueinfo:
                truehist = axis.hist(gdat.truecnts[l][i, :, m], gdat.binscnts[i, :], alpha=0.5, color='g', log=True, label=gdat.truelabl)
                if gdat.datatype == 'mock':
                    if gdat.exprtype == 'ferm':
                        axis.hist(gdat.exprcnts[i, :, m], gdat.binscnts[i, :], alpha=0.1, color='red', log=True, label='3FGL')
            axis.hist(gdatmodi.thiscnts[l][i, :, m], gdat.binscnts[i, :], alpha=0.5, color='b', log=True, label='Sample')
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            axis.set_xscale('log')
            if gdat.trueinfo:
                axis.set_ylim([0.5, 1e3])
            if m == 0:
                axis.set_title(gdat.binsenerstrg[i])
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
            if m == gdat.numbevtt / 2 and i == gdat.numbener / 2:
                axis.legend()
        
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'histcnts_pop%d_' % l + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    

def plot_diagfram(gdat, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'diagthis')
    axis, cbar = retr_imag(gdat, axis, gdatmodi.thispntsfluxmodi, indxenerplot, indxevttplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'diagnext')
    axis, cbar = retr_imag(gdat, axis, gdat.nextpntsfluxmodi, indxenerplot, indxevttplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'diagdiff')
    axis, cbar = retr_imag(gdat, axis, gdat.diffpntsfluxmodi, indxenerplot, indxevttplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)


def plot_nextstat(gdat, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'nextstat')
    axis, cbar = retr_imag(gdat, axis, gdat.datacnts, indxenerplot, indxevttplot)
   
    mrkrsize = retr_mrkrsize(gdat, fabs(gdatmodi.modispec), indxenerplot) * 10
    for k in range(gdat.numbmodipnts):
        if gdatmodi.modispec[indxenerplot, k] > 0:
            colr = 'yellow'
        else:
            colr = 'red'
        if gdat.exprtype == 'ferm':
            xaxi = gdatmodi.modilgal[k]
            yaxi = gdatmodi.modibgal[k]
        else:
            xaxi = gdatmodi.modilgal[k] * 3600.
            yaxi = gdatmodi.modibgal[k] * 3600.
        axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=gdat.mrkralph, marker='o', linewidth=2, color=colr)
        text = r'%s, lnL = %.3g' % (gdat.strgprop[gdatmodi.thisindxprop], gdat.deltllik)
        if gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg:
            text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (gdatmodi.laccfrac, gdatmodi.thisjcbnfact, gdatmodi.thiscombfact, gdatmodi.listaccp[j])
        axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)


def plot_datacnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'datacnts')
    axis, cbar = retr_imag(gdat, axis, gdat.datacnts, indxenerplot, indxevttplot, satuuppr=gdat.datacntssatu)
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_modlcnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'modlcnts')
    axis, cbar = retr_imag(gdat, axis, gdatmodi.thismodlcnts, indxenerplot, indxevttplot, satuuppr=gdat.datacntssatu)
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_resicnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'resicnts')
    axis, cbar = retr_imag(gdat, axis, gdatmodi.thisresicnts, indxenerplot, indxevttplot, \
        satulowr=-gdat.resicntssatu, satuuppr=gdat.resicntssatu, cmap='RdBu')
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_errrpnts(gdat, gdatmodi, indxenerplot, indxevttplot):

    # temp
    #satulowr = -1.
    #satuuppr = 1.
    satulowr = None
    satuuppr = None

    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'errrpnts')
    axis, cbar = retr_imag(gdat, axis, gdatmodi.thiserrrpnts, indxenerplot, indxevttplot, satulowr=satulowr, satuuppr=satuuppr, cmap='RdBu', mean=True, logt=True)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_catlfram(gdat, gdatmodi, indxpoplplot, indxenerplot, indxevttplot):
    
    figr, axis, path = init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, 'catlfram')

    if gdat.trueinfo:
        for l in gdat.indxpopl:
            numbpnts = int(gdat.truenumbpnts[l])
            for a in range(numbpnts):
                if indxevttplot == None:
                    cnts = sum(gdat.truecnts[l][indxenerplot, a, :])
                    sigm = sqrt(sum(gdat.truesigm[l][indxenerplot, a, :]**2))
                else:
                    cnts = gdat.truecnts[l][indxenerplot, a, indxevttplot]
                    sigm = gdat.truesigm[l][indxenerplot, a, indxevttplot]
                axis.text(gdat.truelgal[l][a] + 0.7, gdat.truebgal[l][a] - 0.7, '%d/%.2f' % (cnts, sigm), color='g', fontsize=13)

    for l in gdat.indxpopl:
        numbpnts = int(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[l]])
        for a in range(numbpnts):
            if indxevttplot == None:
                cnts = sum(gdatmodi.thiscnts[l][indxenerplot, a, :])
            else:
                cnts = gdatmodi.thiscnts[l][indxenerplot, a, indxevttplot]
            axis.text(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l][a]] - 0.5, gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l][a]] + 0.3, \
                                '%d' % cnts, color='b', fontsize=13)
    
    if indxevttplot == None:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % mean(gdatmodi.thisbackcntsmean[indxenerplot, :]), fontsize=18)
    else:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % gdatmodi.thisbackcntsmean[indxenerplot, indxevttplot], fontsize=18)
        
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)

    
def plot_pntsdiff():

    tempindxpntsfull = []
    tempnumbpnts = 10
    tempindxpntsfull.append(range(tempnumbpnts))
    tempsamp = rand(gdat.numbpara)
    tempsamp[gdat.indxsampnumbpnts] = array([tempnumbpnts])
    tempsamp[gdat.indxsampfdfnnorm] = cdfn_logt(array([tempnumbpnts]), minmfdfnnorm, factfdfnnorm)
    tempsamp[gdat.indxsampfdfnslop] = cdfn_atan(array([1.5]), minmfdfnslop, factfdfnslop)
    tempsamp[gdat.indxsamppsfipara] = 0.5
    for c in gdat.indxback:
        tempsamp[gdat.indxsampnormback[c, 0]] = cdfn_logt(array([1.]), minmnormback[c], factnormback[c])
    
    tempsampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, drmcsamp[:, 0])
    temppntsflux, temppntscnts, tempmodlflux, tempmodlcnts = retr_maps(gdat, tempindxpntsfull, tempsampvarb)

    pntscnts = tempcnts[0, :, 0] * gdat.expo[0, :, 0] * gdat.apix * gdat.diffener[0]
    isotcnts = isotflux[0, :, 0] * gdat.expo[0, :, 0] * gdat.apix * gdat.diffener[0]
    
    totlcnts0 = isotcnts
    totlcnts1 = pntscnts / 1e6 + isotcnts
    
    nrept = 100
    totlcntsheal0 = zeros((nrept, gdat.numbpixlheal))
    totlcntsheal1 = zeros((nrept, gdat.numbpixlheal))
    for k in range(nrept):
        for n in range(gdat.indxpixlrofi.size):
            totlcntsheal0[k, gdat.indxpixlrofi[n]] = poisson(totlcnts0[n])
            totlcntsheal1[k, gdat.indxpixlrofi[n]] = poisson(totlcnts1[n])
          
    maxmcnts = max(amax(totlcntsheal0), amax(totlcntsheal1))
    gdat.binscnts = linspace(0., maxmcnts, maxmcnts + 2)
    diffcnts = gdat.binscnts[1:] - gdat.binscnts[:-1]
    meancnts = (gdat.binscnts[1:] + gdat.binscnts[:-1]) / 2.
    
    hist0 = empty((nrept, maxmcnts + 1))
    hist1 = empty((nrept, maxmcnts + 1))
    for k in range(nrept):
        hist0[k, :] = histogram(totlcntsheal0[k, gdat.indxpixlrofi], gdat.binscnts)[0].astype(float)
        hist0[k, :] *= 1. / sum(hist0[k, :]) / diffcnts
        hist1[k, :] = histogram(totlcntsheal1[k, gdat.indxpixlrofi], gdat.binscnts)[0].astype(float)
        hist1[k, :] *= 1. / sum(hist1[k, :]) / diffcnts
        
    
    totlcntscart0 = tdpy.util.retr_cart(totlcntsheal0[0, :], minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
    totlcntscart1 = tdpy.util.retr_cart(totlcntsheal1[0, :], minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
    
    fig = plt.figure(figsize=(2 * gdat.plotsize, 2 * gdat.plotsize))
    axis = figr.add_subplot(221)
    imag = axis.imshow(totlcntscart0, origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    
    axis.set_xlabel(gdat.longlabl)
    axis.set_ylabel(gdat.latilabl)
    if gdat.exprtype == 'ferm':
        axis.set_xlim([gdat.maxmlgal, gdat.minmlgal])
        axis.set_ylim([gdat.minmbgal, gdat.maxmbgal])
    else:
        axis.set_xlim(array([gdat.maxmlgal, gdat.minmlgal]) * 3600.)
        axis.set_ylim(array([gdat.minmbgal, gdat.maxmbgal]) * 3600.)
    
    axis.set_title('Isotropic')
    
    axis = figr.add_subplot(222)
    imag = axis.imshow(totlcntscart1, origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.set_xlabel(r'$x$ [$^\circ$]')
    axis.set_xlim([gdat.maxmlgal, gdat.minmlgal])
    axis.set_ylim([gdat.minmbgal, gdat.maxmbgal])
    axis.set_title('Isotropic + Unresolved PS')
    
    axis.scatter(tempsampvarb[trueindxsamplgal], tempsampvarb[trueindxsampbgal], s=50, alpha=0.8, marker='x', color='g', linewidth=2)
    
    axis = figr.add_subplot(212)

    tdpy.mcmc.plot_braz(ax, meancnts, hist0,  lcol='lightgreen', alpha=0.3, dcol='green', mcol='darkgreen', labl='Isotropic')
    tdpy.mcmc.plot_braz(ax, meancnts, hist1, lcol='lightblue', alpha=0.3, dcol='blue', mcol='darkblue', labl='Isotropic + Unresolved PS')

    axis.set_xlabel('$k$')
    axis.set_ylabel('d$P$/d$k$')
    axis.set_yscale('log')
    axis.legend()

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'pntsdiff.pdf')
    plt.close(figr)
    
