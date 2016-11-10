# common imports
from __init__ import *

# internal functions
from util import *

def plot_post(pathpcat, verbtype=1, makeanim=False):
    
    gdat = tdpy.util.gdatstrt()
    gdat.verbtype = verbtype

    if gdat.verbtype > 0:
        print 'Producing plots for the PCAT output...'

    # read PCAT output file
    if gdat.verbtype > 0:
        print 'Reading %s...' % pathpcat
    hdun = pf.open(pathpcat)

    gdat.numbener = hdun[0].header['numbener']
    gdat.numbevtt = hdun[0].header['numbevtt']
    gdat.numbpopl = hdun[0].header['numbpopl']
    
    gdat.numbproc = hdun[0].header['numbproc']
    gdat.pathbase = hdun[0].header['pathbase']
    
    gdat.pntstype = hdun[0].header['pntstype']
    gdat.truevarioaxi = hdun[0].header['truevarioaxi']
    gdat.truepsfntype = hdun[0].header['truepsfntype']
    
    gdat.anglfact = hdun[0].header['anglfact']
    gdat.strganglunit = hdun[0].header['strganglunit']
    gdat.modlvarioaxi = hdun[0].header['modlvarioaxi']
    
    #gdat.numbpsfp = hdun[0].header['numbpsfp']
    #gdat.numbpsfpform = hdun[0].header['numbformpara']
    
    gdat.optiprop = hdun[0].header['optiprop']
    
    gdat.indxener = arange(gdat.numbener)
    gdat.indxevtt = arange(gdat.numbevtt)
    gdat.indxpopl = arange(gdat.numbpopl)
    
    gdat.numbsamp = hdun[0].header['numbsamp']
    gdat.numbburn = hdun[0].header['numbburn']
    gdat.numbswep = hdun[0].header['numbswep']
    gdat.numbswepplot = hdun[0].header['numbswepplot']
    gdat.factthin = hdun[0].header['factthin']
    gdat.numbpopl = hdun[0].header['numbpopl']
    gdat.numbproc = hdun[0].header['numbproc']
    
    gdat.strgcatl = hdun[0].header['strgcatl']
    gdat.strgfunctime = hdun[0].header['strgfunctime']

    gdat.numbfluxprox = hdun[0].header['numbfluxprox']
    timeatcr = hdun[0].header['timeatcr']
    
    gdat.anglassc = hdun[0].header['anglassc']
    gdat.maxmgang = hdun[0].header['maxmgang']
    gdat.numbfluxplot = hdun[0].header['numbfluxplot']

    gdat.minmlgal = hdun[0].header['minmlgal']
    gdat.maxmlgal = hdun[0].header['maxmlgal']
    gdat.minmbgal = hdun[0].header['minmbgal']
    gdat.maxmbgal = hdun[0].header['maxmbgal']
    
    # plotting
    gdat.anotcatl = hdun[0].header['anotcatl']
   
    gdat.numbsidepntsprob = hdun[0].header['numbsidepntsprob']

    gdat.minmflux = hdun[0].header['minmflux']
    gdat.maxmflux = hdun[0].header['maxmflux']

    gdat.datatype = hdun[0].header['datatype']
    gdat.modlpsfntype = hdun[0].header['modlpsfntype']
    gdat.exprtype = hdun[0].header['exprtype']
    gdat.pixltype = hdun[0].header['pixltype']
    
    gdat.lgalcntr = hdun[0].header['lgalcntr']
    gdat.bgalcntr = hdun[0].header['bgalcntr']
    gdat.scalmaps = hdun[0].header['scalmaps']
    gdat.satumaps = hdun[0].header['satumaps']
    
    gdat.probpsfp = hdun[0].header['probpsfp']
    
    gdat.exprinfo = hdun[0].header['exprinfo']
    gdat.priotype = hdun[0].header['priotype']
    gdat.priofactdoff = hdun[0].header['priofactdoff']

    gdat.spatdisttype = []
    gdat.fluxdisttype = []
    gdat.spectype = []
    for l in gdat.indxpopl:
        gdat.spatdisttype.append(hdun[0].header['spatdisttypepop%d' % l])
        gdat.fluxdisttype.append(hdun[0].header['fluxdisttypepop%d' % l])
        gdat.spectype.append(hdun[0].header['spectypepop%d' % l])

    if gdat.pixltype == 'heal':
        gdat.numbsideheal = hdun[0].header['numbsideheal']
    if gdat.pixltype == 'cart':
        gdat.numbsidecart = hdun[0].header['numbsidecart']
    if gdat.pixltype == 'unbd':
        gdat.numbdatasamp = hdun[0].header['numbdatasamp']

    # best-fit sample
    ## maximum likelihood
    gdat.indxswepmaxmllik = hdun[0].header['indxswepmaxmllik']
    gdat.maxmllikswep = hdun[0].header['maxmllikswep']
    gdat.indxswepmaxmlpos = hdun[0].header['indxswepmaxmlpos']
    gdat.maxmlposswep = hdun[0].header['maxmlposswep']
        
    gdat.strgenerunit = hdun[0].header['strgenerunit']
    gdat.strgflux = hdun[0].header['strgflux']
    gdat.strgfluxunit = hdun[0].header['strgfluxunit']
    gdat.maxmangl = hdun[0].header['maxmangl']
    gdat.margfactmodl = hdun[0].header['margfactmodl']
    gdat.margfactcomp = hdun[0].header['margfactcomp']
    gdat.strgtimestmp = hdun[0].header['strgtimestmp']
    gdat.strgxaxi = hdun[0].header['strgxaxi']
    gdat.strgyaxi = hdun[0].header['strgyaxi']
    gdat.strgcnfg = hdun[0].header['strgcnfg']
    gdat.numbback = hdun[0].header['numbback']

    gdat.stdvmeanpnts = hdun[0].header['stdvmeanpnts']
    gdat.stdvfluxdistslop = hdun[0].header['stdvfluxdistslop']
    gdat.stdvproppsfp = hdun[0].header['stdvproppsfp']
    gdat.stdvback = hdun[0].header['stdvback']
    gdat.stdvlbhl = hdun[0].header['stdvlbhl']
    gdat.stdvflux = hdun[0].header['stdvflux']
    gdat.radispmr = hdun[0].header['radispmr']
    gdat.fracrand = hdun[0].header['fracrand']
    
    gdat.psfninfoprio = hdun[0].header['psfninfoprio']
    
    gdat.specfraceval = hdun[0].header['specfraceval']

    if gdat.datatype == 'mock':
        gdat.mocknumbpopl = hdun[0].header['mocknumbpopl']
        gdat.mockminmflux = hdun[0].header['mockminmflux']
        gdat.mockmaxmflux = hdun[0].header['mockmaxmflux']

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
    gdat.listspephist = hdun['spephist'].data
    
    gdat.pntsprob = hdun['pntsprob'].data
    
    listllik = hdun['llik'].data
    listlpri = hdun['lpri'].data
    
    atcr = hdun['atcr'].data
    gmrbstat = hdun['gmrbstat'].data
    indxpixlsave = hdun['indxpixlsave'].data
    listmodlcnts = hdun['modlcnts'].data
    
    if gdat.datatype == 'mock':
        gdat.mockspatdisttype = []
        gdat.mockfluxdisttype = []
        gdat.mockspectype = []
        for l in gdat.indxpopl:
            gdat.mockspatdisttype.append(hdun[0].header['mockspatdisttypepop%d' % l])
            gdat.mockfluxdisttype.append(hdun[0].header['mockfluxdisttypepop%d' % l])
            gdat.mockspectype.append(hdun[0].header['mockspectypepop%d' % l])
        
        gdat.mockvarioaxi = hdun[0].header['mockvarioaxi']
        gdat.mockfluxdistslop = hdun['mockfluxdistslop'].data
        gdat.mockfluxdistsloplowr = hdun['mockfluxdistsloplowr'].data
        gdat.mockfluxdistslopuppr = hdun['mockfluxdistslopuppr'].data
        gdat.mockfluxdistbrek = hdun['mockfluxdistbrek'].data
        gdat.mocksinddistmean = hdun['mocksinddistmean'].data
        gdat.mocksinddiststdv = hdun['mocksinddiststdv'].data
        gdat.mockcurvdistmean = hdun['mockcurvdistmean'].data
        gdat.mockcurvdiststdv = hdun['mockcurvdiststdv'].data
        gdat.mocknormback = hdun['mocknormback'].data
        gdat.mocknumbpnts = hdun['mocknumbpnts'].data

    # prior boundaries
    gdat.minmnormback = hdun['minmnormback'].data
    gdat.maxmnormback = hdun['maxmnormback'].data

    gdat.sinddistmean = hdun['sinddistmean'].data
    gdat.sinddiststdv = hdun['sinddiststdv'].data
        
    gdat.curvdistmean = hdun['curvdistmean'].data
    gdat.curvdiststdv = hdun['curvdiststdv'].data
        
    # bins
    gdat.indxevttfull = hdun['indxevttfull'].data
    gdat.binsenerfull = hdun['binsenerfull'].data
    gdat.indxenerfluxdist = array([hdun[0].header['indxenerfluxdist']])
   
    # hyperprior limits
    gdat.minmmeanpnts = hdun['minmmeanpnts'].data 
    gdat.maxmmeanpnts = hdun['maxmmeanpnts'].data
    gdat.minmfluxdistslop = hdun['minmfluxdistslop'].data
    gdat.maxmfluxdistslop = hdun['maxmfluxdistslop'].data
    gdat.minmfluxdistbrek = hdun['minmfluxdistbrek'].data
    gdat.maxmfluxdistbrek = hdun['maxmfluxdistbrek'].data
    gdat.minmfluxdistsloplowr = hdun['minmfluxdistsloplowr'].data
    gdat.maxmfluxdistsloplowr = hdun['maxmfluxdistsloplowr'].data
    gdat.minmfluxdistslopuppr = hdun['minmfluxdistslopuppr'].data
    gdat.maxmfluxdistslopuppr = hdun['maxmfluxdistslopuppr'].data
        
    # utilities
    gdat.listsamp = hdun['listsamp'].data
    gdat.listsampvarb = hdun['listsampvarb'].data
    gdat.listindxpntsfull = hdun['listindxpntsfull'].data
    gdat.probprop = hdun['probprop'].data
    gdat.listindxprop = hdun['indxprop'].data
    gdat.listchrototl = hdun['listchrototl'].data
    gdat.listchrollik = hdun['listchrollik'].data
    gdat.listaccp = hdun['accp'].data
    gdat.listindxparamodi = hdun['indxparamodi'].data
    gdat.listauxipara = hdun['auxipara'].data
    gdat.listlaccfact = hdun['laccfact'].data
    gdat.listnumbpair = hdun['numbpair'].data
    gdat.listjcbnfact = hdun['jcbnfact'].data
    gdat.listcombfact = hdun['combfact'].data
    gdat.listboolreje = hdun['boolreje'].data
    gdat.listpntsfluxmean = hdun['listpntsfluxmean'].data
    gdat.listdeltllik = hdun['listdeltllik'].data
    gdat.listdeltlpri = hdun['listdeltlpri'].data
    gdat.listmemoresi = hdun['listmemoresi'].data

    # posterior distributions
    gdat.listnumbpnts = hdun['numbpnts'].data
    gdat.listmeanpnts = hdun['meanpnts'].data
    gdat.listfluxdistslop = hdun['fluxdistslop'].data
    gdat.listfluxdistsloplowr = hdun['fluxdistsloplowr'].data
    gdat.listfluxdistslopuppr = hdun['fluxdistslopuppr'].data
    gdat.listfluxdistbrek = hdun['fluxdistbrek'].data
    gdat.listpsfp = hdun['psfp'].data
    gdat.listnormback = hdun['normback'].data
    
    # best-fit sample
    ## maximum likelihood
    gdat.sampvarbmaxmllik = hdun['sampvarbmaxmllik'].data
    gdat.sampvarbmaxmlpos = hdun['sampvarbmaxmlpos'].data

    strgbest = 'ML'
    
    gdat.makeplot = True
    gdat.diagmode = False
    
    gdat.strgbinsener = []
    for i in range(gdat.numbener):
        gdat.strgbinsener.append(hdun[0].header['strgbinsener%04d' % i])
    
    # initial setup
    setpinit(gdat) 
    
    # mock catalog
    if gdat.datatype == 'mock':
        gdat.mockdatacnts = hdun['mockdatacnts'].data
        gdat.mocknumbpnts = hdun['mocknumbpnts'].data
        gdat.mocklgal = []
        gdat.mockbgal = []
        gdat.mockspec = []
        if gdat.numbener > 1:
            gdat.mockspep = []
        gdat.mockcnts = []
        for l in gdat.indxpopl:
            gdat.mocklgal.append(hdun['mocklgalpop%d' % l].data)
            gdat.mockbgal.append(hdun['mockbgalpop%d' % l].data)
            gdat.mockspec.append(hdun['mockspecpop%d' % l].data)
            if gdat.numbener > 1:
                gdat.mockspep.append(hdun['mockspeppop%d' % l].data)
            gdat.mockcnts.append(hdun['mockcntspop%d' % l].data)
        gdat.mockfluxdistslop = hdun['mockfluxdistslop'].data
        gdat.mockfluxdistbrek = hdun['mockfluxdistbrek'].data
        gdat.mockfluxdistsloplowr = hdun['mockfluxdistsloplowr'].data
        gdat.mockfluxdistslopuppr = hdun['mockfluxdistslopuppr'].data
        gdat.mocknormback = hdun['mocknormback'].data
        
    # final setup
    setpfinl(gdat) 
    
    gdat.listlgal = [[] for l in gdat.indxpopl]
    gdat.listbgal = [[] for l in gdat.indxpopl]
    gdat.listspec = [[] for l in gdat.indxpopl]
    gdat.listspep = [[] for l in gdat.indxpopl]
    gdat.listgang = [[] for l in gdat.indxpopl]
    gdat.listaang = [[] for l in gdat.indxpopl]
    for l in gdat.indxpopl:
        lgal = hdun['lgalpop%d' % l].data
        bgal = hdun['bgalpop%d' % l].data
        spec = hdun['specpop%d' % l].data
        spep = hdun['speppop%d' % l].data
        gang = hdun['gangpop%d' % l].data
        aang = hdun['aangpop%d' % l].data
        for j in gdat.indxsamptotl:
            indxpnts = where(spec[j, gdat.indxenerfluxdist[0], :] > 0.)[0]
            numbpnts = indxpnts.size
            gdat.listlgal[l].append(lgal[j, indxpnts])
            gdat.listbgal[l].append(bgal[j, indxpnts])
            # temp
            gdat.listspec[l].append(spec[j, :, indxpnts].T)
            gdat.listspep[l].append(spep[j, indxpnts, :])
            gdat.listgang[l].append(gang[j, indxpnts])
            gdat.listaang[l].append(aang[j, indxpnts])

    timetotlinit = gdat.functime()
    
    ## convert the table into a list of lists
    listindxpntsfulltemp = [[[] for l in gdat.indxpopl] for n in gdat.indxsamptotl]
    for n in gdat.indxsamptotl:
        for l in gdat.indxpopl:
            numbpnts = gdat.listnumbpnts[n, l]
            listindxpntsfulltemp[n][l] = gdat.listindxpntsfull[n, l, :numbpnts]
    gdat.listindxpntsfull = listindxpntsfulltemp
    
    # indices of the parameters to be plotted
    numbparaplot = min(gdat.numbpara, 50)
    size = numbparaplot - gdat.indxsampcompinit
    indxparaplot = concatenate([arange(gdat.indxsampcompinit), sort(choice(arange(gdat.indxsampcompinit, gdat.numbpara), size=size, replace=False))])

    # Gelman-Rubin test
    if gdat.numbproc > 1:
        if isfinite(gmrbstat).all():
            figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
            bins = linspace(1., amax(gmrbstat), 40)
            axis.hist(gmrbstat, bins=bins)
            axis.set_xlabel('PSRF')
            axis.set_ylabel('$N_{pix}$')
            plt.tight_layout()
            figr.savefig(gdat.pathplot + 'diag/gmrbhist.pdf')
            plt.close(figr)
            path = gdat.pathplot + 'diag/gmrbheal.pdf'
            maps = zeros(gdat.numbpixl)
            maps[indxpixlsave] = gmrbstat
            tdpy.util.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
        else:
            print 'Inappropriate Gelman-Rubin test statistics encountered.'

    # plot autocorrelation
    tdpy.mcmc.plot_atcr(gdat.pathdiag, atcr, timeatcr)
    
    # plot proposal efficiency
    if gdat.verbtype > 0:
        print 'Making proposal efficiency plots...'
    # temp
    if False:
        numbtimemcmc = 20
        binstimemcmc = linspace(0., gdat.numbswep, numbtimemcmc)
        numbtick = 2
        figr, axgr = plt.subplots(gdat.numbprop, 1, figsize=(gdat.plotsize, gdat.numbprop * gdat.plotsize / 4.), sharex='all')
        for n, axis in enumerate(axgr):
            histtotl = axis.hist(where(gdat.listindxprop == n)[0], bins=binstimemcmc)[0]
            histaccp = axis.hist(intersect1d(where(gdat.listindxprop == n)[0], where(gdat.listaccp == True)[0]), bins=binstimemcmc)[0]
            axis.set_ylabel('%s' % gdat.strgprop[n])
            if n == gdat.numbprop - 1:
                axis.set_xlabel('$i_{samp}$')
            # define the y-axis
            maxm = amax(histtotl)
            axis.set_ylim([0., maxm])
            try:
                axis.set_title('%03d' % int(sum(histaccp) / sum(histtotl) * 100.))
            except:
                pass
            listtick = linspace(maxm / 2., maxm, numbtick)
            listlabltick = ['%.3g' % tick for tick in listtick]
            axis.set_yticks(listtick)
            axis.set_yticklabels(listlabltick)
        plt.tight_layout()
        figr.savefig(gdat.pathdiag + 'accpratiprop.pdf')
        plt.close(figr)
        
        figr, axgr = plt.subplots(numbparaplot, 1, figsize=(gdat.plotsize, numbparaplot * gdat.plotsize / 4.), sharex='all')
        for n, axis in enumerate(axgr):
            hist = axis.hist(where(gdat.listindxparamodi == indxparaplot[n])[0], bins=binstimemcmc)[0]
            axis.hist(intersect1d(where(gdat.listindxparamodi == n)[0], where(gdat.listaccp)[0]), bins=binstimemcmc)
            axis.set_ylabel('$p_{%d}$' % indxparaplot[n])
            if n == gdat.numbprop - 1:
                axis.set_xlabel('$i_{samp}$')
            
            # define the y-axis
            maxm = amax(hist[0])
            axis.set_ylim([0., maxm])
            listtick = linspace(maxm / 2., maxm, numbtick)
            listlabltick = ['%.3g' % tick for tick in listtick]
            axis.set_yticks(listtick)
            axis.set_yticklabels(listlabltick)
        
        plt.subplots_adjust(hspace=0)
        figr.savefig(gdat.pathdiag + 'accpratipara.pdf')
        plt.close(figr)
   
        # plot split and merge diagnostics
        indxsampsplttotl = where(gdat.listindxprop == gdat.indxpropsplt)[0]
        indxsampsplt = intersect1d(where(gdat.listindxprop == gdat.indxpropsplt)[0], where(gdat.listboolreje == False)[0])
        indxsampmergtotl = where(gdat.listindxprop == gdat.indxpropmerg)[0]
        indxsampmerg = where((gdat.listindxprop == gdat.indxpropmerg) & (gdat.listboolreje == False))[0]
        indxsampspmr = union1d(indxsampsplt, indxsampmerg)
        indxsampspmrtotl = union1d(indxsampsplttotl, indxsampmergtotl)
        indxsampreje = where(gdat.listboolreje == False)
        if indxsampspmrtotl.size > 0:

            ## create plot subfolder
            os.system('mkdir -p %s' % gdat.pathplot + 'spmr')

            ## labels and names
            listlabl = ['$u_f$', '$u_r$', r'$u_\phi$', '$u_s$', '$N_{pair}$', \
                                            r'$\alpha_c\alpha_j$', r'$\alpha_P\alpha_c\alpha_j$', r'$\alpha_c$', r'$\alpha_j$', r'$\alpha_L$', r'$\alpha_P$']
            listname = ['fracauxi', 'radiauxi', 'anglauxi', 'sindauxi', 'numbpair', 'laccfact', 'laccfacttotl', 'combfact', 'jcbnfct', 'deltllik', 'deltlpri']
            
            ## variables
            listvarb = [gdat.listauxipara[:, 0], gdat.anglfact * gdat.listauxipara[:, 1], gdat.listauxipara[:, 2], gdat.listauxipara[:, 3], gdat.listnumbpair, \
                                             exp(gdat.listlaccfact), exp(gdat.listlaccfact + gdat.listdeltlpri), gdat.listcombfact, \
                                             gdat.listjcbnfact, exp(gdat.listdeltllik), exp(gdat.listdeltlpri)]
           
            print 'gdat.listauxipara'
            print gdat.listauxipara
            print 'gdat.listnumbpair'
            print gdat.listnumbpair
            print 'gdat.listlaccfact'
            print gdat.listlaccfact
            print 'gdat.listdeltlpri'
            print gdat.listdeltlpri
            print 'gdat.listcombfact'
            print gdat.listcombfact
            print 'gdat.listdeltllik'
            print gdat.listdeltllik[indxsampsplttotl]
            print gdat.listdeltllik[indxsampmergtotl]
            print 'gdat.listdeltlpri'
            print gdat.listdeltlpri[indxsampsplttotl]
            print gdat.listdeltlpri[indxsampmergtotl]
            print

            # rescale the Jacobian and the prior fraction to make them dimensionless
            listvarb[5][indxsampsplttotl] /= gdat.minmflux
            listvarb[5][indxsampmergtotl] *= gdat.minmflux
            listvarb[8][indxsampsplttotl] /= gdat.minmflux
            listvarb[8][indxsampmergtotl] *= gdat.minmflux
            listvarb[10][indxsampsplttotl] *= gdat.minmflux
            listvarb[10][indxsampmergtotl] /= gdat.minmflux

            for k in range(len(listlabl)):
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                maxm = amax(listvarb[k][indxsampspmrtotl])
                
                if k <= 4:
                    indxsampspmrtemp = indxsampspmrtotl
                    indxsampsplttemp = indxsampsplttotl
                else:
                    indxsampspmrtemp = indxsampspmr
                    indxsampsplttemp = indxsampsplt
                
                if k >= 5:
                    axis.set_xscale('log')
                    minm = amin(listvarb[k][indxsampspmrtemp][where(listvarb[k][indxsampspmrtemp] > 0.)])
                    bins = logspace(log10(minm), log10(maxm), 40)
                else:
                    minm = amin(listvarb[k][indxsampspmrtemp])
                    bins = linspace(minm, maxm, 40)
              
                try:
                    axis.hist(listvarb[k][indxsampsplttemp], bins=bins, label='Split', alpha=gdat.alphmrkr)
                    axis.hist(listvarb[k][indxsampmerg], bins=bins, label='Merge', alpha=gdat.alphmrkr)
                except:
                    if gdat.verbtype > 0:
                        print 'Skipping split merge diagnostic plot...'
                axis.set_ylabel('$N_{samp}$')
                axis.set_xlabel(listlabl[k])
                axis.legend()
                plt.tight_layout()
                figr.savefig(gdat.pathplot + 'spmr/' + listname[k] + '.pdf')
                plt.close(figr)
    
    if gdat.verbtype > 0:
        print 'Calculating proposal execution times...'
    timeinit = gdat.functime()
    plot_chro(gdat)
    timefinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # temp
    if False:
        path = gdat.pathplot + 'listsamp_'
        tdpy.mcmc.plot_grid(path, gdat.listsamp[:, indxparaplot], ['%d' % k for k in indxparaplot], numbplotside=10)

    numbpntspost = 3
    # temp -- posterior plots only work for numbcompcolr == 4
    if gdat.strgcnfg == 'test_post' and gdat.datatype == 'mock' and gdat.mocknumbpnts[0] == numbpntspost and gdat.numbpopl == 1 and gdat.numbback == 1:
        if gdat.verbtype > 0:
            print 'Making the grid posterior plot...'
        timeinit = gdat.functime()
   
        numbpara = numbpntspost * gdat.numbcompcolr + gdat.numbener
        listpost = zeros((gdat.numbsamp, numbpara))
        for j in gdat.indxsamptotl:
            for k in range(numbpntspost):
                listpost[j, 0*numbpntspost+k] = gdat.listlgal[0][j][k]
                listpost[j, 1*numbpntspost+k] = gdat.listbgal[0][j][k]
                listpost[j, 2*numbpntspost+k] = gdat.listspec[0][j][gdat.indxenerfluxdist[0], k]
                listpost[j, 3*numbpntspost+k] = gdat.listspep[0][j][k, :]
        for i in gdat.indxener:
            listpost[:, 4*numbpntspost+i] = gdat.listnormback[:, 0, i]
        truepost = zeros(numbpara)
        truepost[0*numbpntspost:1*numbpntspost] = gdat.truelgal[0][k]
        truepost[1*numbpntspost:2*numbpntspost] = gdat.truebgal[0][k]
        truepost[2*numbpntspost:3*numbpntspost] = gdat.truespec[0][0, gdat.indxenerfluxdist[0], k]
        truepost[3*numbpntspost:4*numbpntspost] = gdat.truespep[0][k, 0]
        truepost[4*numbpntspost:] = gdat.truenormback[0, :]
        path = gdat.pathplot + 'postdist'
        strgpara = ['$%s_%d$' % (strg, indxpnts + 1) for strg in ['l', 'b', 'f', 's'] for indxpnts in arange(numbpnts)]
        strgpara += ['$A_{%d}$' % i for i in gdat.indxener]
        
        tdpy.mcmc.plot_grid(path, listpost, strgpara, truepara=truepost, numbtickbins=3)

        # find the matrix of partial derivatives
        ## evaluate the likelihood at the sample
        # temp
        #gdatmodi = tdpy.util.gdatstrt()
        #retr_llik(gdat, gdatmodi, init=True)

        ## evaluate the likelihood at the vicinity of the sample
        #for a in range():
        #    for b in range():
        #        retr_llik(gdat, gdatmodi)

        timefinl = gdat.functime()
        if gdat.verbtype > 0:
            print 'Done in %.3g seconds.' % (timefinl - timeinit)

    # cross correlation with the reference catalog
    if gdat.verbtype > 0:
        print 'Cross-correlating the sample catalogs with the reference catalog...'
    timeinit = gdat.functime()
    if gdat.trueinfo:
        for l in gdat.indxpopl:
            listindxmodl = []
            for n in gdat.indxsamptotl:
                indxmodl, indxtruepntsassc = corr_catl(gdat, l, gdat.listlgal[l][n], gdat.listbgal[l][n], gdat.listspec[l][n])
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
    timefinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Done in %.3g seconds.' % (timefinl - timeinit)
    
    # stacked posteiors binned in position and flux
    plot_pntsprob(gdat, ptag='quad')
    plot_pntsprob(gdat, ptag='full', full=True)
    plot_pntsprob(gdat, ptag='cumu', cumu=True)

    # compute median and 68% credible intervals of the inferred parameters
    gdat.postnormback = tdpy.util.retr_postvarb(gdat.listnormback)
    gdat.postpsfp = tdpy.util.retr_postvarb(gdat.listpsfp)
    
    # compute the medians of secondary variables
    gdat.medinormback = gdat.postnormback[0, :, :]
    if gdat.pntstype == 'lght':
        gdat.medipsfp = gdat.postpsfp[0, :]
        gdat.medipsfn = retr_psfn(gdat, gdat.medipsfp, gdat.indxener, gdat.binsangl, gdat.modlpsfntype)
        gdat.medifwhm = 2. * retr_psfnwdth(gdat, gdat.medipsfn, 0.5)
        if gdat.correxpo:
            gdat.medicntsbackfwhm = retr_cntsbackfwhm(gdat, gdat.medinormback, gdat.medifwhm)
            gdat.medibinssigm = retr_sigm(gdat, gdat.binscnts, gdat.medicntsbackfwhm)

    # flux distribution
    for l in gdat.indxpopl:
        plot_histspec(gdat, l, listspechist=gdat.listspechist[:, l, :, :])

    # color distribution
    if gdat.pntstype == 'lght' and gdat.numbener > 1:
        for l in gdat.indxpopl:
            plot_histsind(gdat, l, listsindhist=gdat.listspephist[l][:, 0, :])

    # fraction of emission components
    if gdat.numbback == 2:
        postpntsfluxmean = tdpy.util.retr_postvarb(gdat.listpntsfluxmean)
        plot_compfrac(gdat, postpntsfluxmean=postpntsfluxmean)

    # mosaic plot of images of posterior catalogs 
    plot_mosa(gdat)
    
    # PSF parameters
    path = gdat.pathpost + 'psfp'
    if gdat.modlpsfntype == 'singgaus' or gdat.modlpsfntype == 'singking':
        gdat.listpsfp[:, gdat.indxpsfpinit] = gdat.listpsfp[:, gdat.indxpsfpinit] * gdat.anglfact
        if gdat.trueinfo:
            gdat.truepsfp[gdat.indxpsfpinit] = gdat.truepsfp[gdat.indxpsfpinit] * gdat.anglfact
    elif gdat.modlpsfntype == 'doubgaus' or gdat.modlpsfntype == 'gausking':
        gdat.listpsfp[:, gdat.indxpsfpinit+1] = gdat.listpsfp[:, gdat.indxpsfpinit+1] * gdat.anglfact
        gdat.listpsfp[:, gdat.indxpsfpinit+2] = gdat.listpsfp[:, gdat.indxpsfpinit+2] * gdat.anglfact
        if gdat.trueinfo:
            gdat.truepsfp[gdat.indxpsfpinit+1] = gdat.truepsfp[gdat.indxpsfpinit+1] * gdat.anglfact
            gdat.truepsfp[gdat.indxpsfpinit+2] = gdat.truepsfp[gdat.indxpsfpinit+2] * gdat.anglfact
    elif gdat.modlpsfntype == 'doubking':
        gdat.listpsfp[:, gdat.indxpsfpinit+1] = gdat.listpsfp[:, gdat.indxpsfpinit+1] * gdat.anglfact
        gdat.listpsfp[:, gdat.indxpsfpinit+3] = gdat.listpsfp[:, gdat.indxpsfpinit+3] * gdat.anglfact
        if gdat.trueinfo:
            gdat.truepsfp[gdat.indxpsfpinit+1] = gdat.truepsfp[gdat.indxpsfpinit+1] * gdat.anglfact
            gdat.truepsfp[gdat.indxpsfpinit+3] = gdat.truepsfp[gdat.indxpsfpinit+3] * gdat.anglfact

    if gdat.probpsfp != 0.:
        if gdat.trueinfo and gdat.modlpsfntype == 'doubking':
            truepara = gdat.truepsfp
        else:
            truepara = array([None] * gdat.numbpsfp)
        
        tdpy.mcmc.plot_grid(path, gdat.listpsfp, gdat.strgpsfp, truepara=truepara, numbplotside=gdat.numbpsfptotl, numbbinsplot=gdat.numbbins, numbtickbins=3)
        
        for k in range(gdat.numbpsfp):
            if std(gdat.listpsfp[:, k]) != 0:
                path = gdat.pathpost + 'psfp_%d' % k
                tdpy.mcmc.plot_trac(path, gdat.listpsfp[:, k], gdat.strgpsfp[k])
    
    # log-likelihood
    path = gdat.pathpost + 'llik'
    tdpy.mcmc.plot_trac(path, listllik.flatten(), '$P(D|x)$', varbdraw=[gdat.maxmllikswep], labldraw=['Maximum likelihood Sample'])

    # log-prior
    path = gdat.pathpost + 'lpri'
    tdpy.mcmc.plot_trac(path, listlpri.flatten(), '$P(x)$')

    # number, expected number of PS and flux conditional prior power law index 
    for l in gdat.indxpopl:
        
        # number of point sources
        path = gdat.pathpost + 'numbpntsdist_pop%d' % l
        if gdat.trueinfo and gdat.truenumbpnts != None:
            truepara = gdat.truenumbpnts[l]
        else:
            truepara = None
        tdpy.mcmc.plot_trac(path, gdat.listnumbpnts[:, l], '$N$', truepara=truepara)

        # temp
        if False:
            gdat.numbsamp = 1000
            gdat.mocknumbpnts = array([400])
            gdat.listmeanpnts = randn(gdat.numbsamp).reshape((gdat.numbsamp, 1))
            gdat.listmeanpnts[where(gdat.listmeanpnts[:, 0] > 0.)[0], 0] *= 30
            gdat.listmeanpnts[where(gdat.listmeanpnts[:, 0] < 0.)[0], 0] *= 45
            gdat.listmeanpnts += gdat.mocknumbpnts

        # mean number of point sources
        path = gdat.pathpost + 'meanpnts_pop%d' % l
        if gdat.trueinfo and gdat.datatype == 'mock':
            truepara = gdat.mocknumbpnts[l]
        else:
            truepara = None
   
        tdpy.mcmc.plot_trac(path, gdat.listmeanpnts[:, l], '$\mu$', truepara=truepara, )

        # temp
        if False:
            gdat.listfluxdistsloplowr = randn(gdat.numbsamp).reshape((gdat.numbsamp, 1))
            gdat.listfluxdistsloplowr[where(gdat.listfluxdistsloplowr[:, 0] > 0.)[0], 0] *= 0.15
            gdat.listfluxdistsloplowr[where(gdat.listfluxdistsloplowr[:, 0] < 0.)[0], 0] *= 0.1
            gdat.listfluxdistsloplowr += gdat.mockfluxdistsloplowr
            gdat.listfluxdistslopuppr = randn(gdat.numbsamp).reshape((gdat.numbsamp, 1))
            gdat.listfluxdistslopuppr[where(gdat.listfluxdistslopuppr[:, 0] > 0.)[0], 0] *= 0.03
            gdat.listfluxdistslopuppr[where(gdat.listfluxdistslopuppr[:, 0] < 0.)[0], 0] *= 0.05
            gdat.listfluxdistslopuppr += gdat.mockfluxdistslopuppr
            gdat.mockfluxdistsloplowr = 0.1 * rand() + gdat.mockfluxdistsloplowr
            gdat.mockfluxdistslopuppr = 0.05 * randn() + gdat.mockfluxdistslopuppr

        # flux distribution
        titl = gdat.strgbinsener[gdat.indxenerfluxdist[0]]
        if gdat.fluxdisttype[l] == 'powr':
            # power law index
            path = gdat.pathpost + 'fluxdistslop_pop%d' % l
            # temp
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfluxdisttype[l] == 'powr':
                    truepara = gdat.mockfluxdistslop[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha$'
            tdpy.mcmc.plot_trac(path, gdat.listfluxdistslop[:, l], labl, truepara=truepara)
        
        if gdat.fluxdisttype[l] == 'brok':
            # break flux
            path = gdat.pathpost + 'fluxdistbrek_pop%d' % l
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfluxdisttype[l] == 'brok':
                    truepara = gdat.mockfluxdistbrek[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$%s_b$' % gdat.strgflux
            tdpy.mcmc.plot_trac(path, gdat.listfluxdistbrek[:, l], labl, truepara=truepara, scalpara='logt')
        
            # lower power law index
            path = gdat.pathplot + 'fluxdistsloplowr_pop%d' % l
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfluxdisttype[l] == 'brok':
                    truepara = gdat.mockfluxdistsloplowr[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha_1$'
            tdpy.mcmc.plot_trac(path, gdat.listfluxdistsloplowr[:, l], labl, truepara=truepara)
        
            # uppr power law index
            path = gdat.pathpost + 'fluxdistslopuppr_pop%d' % l
            if gdat.trueinfo and gdat.datatype == 'mock':
                if gdat.mockfluxdisttype[l] == 'brok':
                    truepara = gdat.mockfluxdistslopuppr[l]
                else:
                    truepara = None
            else:
                truepara = None
            labl =  r'$\alpha_2$'
            tdpy.mcmc.plot_trac(path, gdat.listfluxdistslopuppr[:, l], labl, truepara=truepara)
        
    # background normalization
    for i in gdat.indxener:
        for c in gdat.indxback:
            path = gdat.pathpost + gdat.nameback[c] + '%d' % i
            if gdat.trueinfo and gdat.datatype == 'mock':
                truepara = gdat.truenormback[c, i]
            else:
                truepara = None
            titl = gdat.strgbinsener[i]
            labl = gdat.lablback[c] + '$_{%d}$' % i
            varbdraw=[gdat.sampvarbmaxmllik[gdat.indxsampnormback[c, i]]]
            tdpy.mcmc.plot_trac(path, gdat.listnormback[:, c, i], labl, truepara=truepara, titl=titl, scalpara='logt', varbdraw=varbdraw, labldraw=[strgbest])

    ## joint distribution of background component pairs
    for i in gdat.indxener:
        for c in gdat.indxback:
            for k in arange(c):
                path = gdat.pathplot + gdat.nameback[c] + gdat.nameback[k] + '%d' % i
                if gdat.trueinfo and gdat.datatype == 'mock':
                    truepara = array([gdat.truenormback[c, i], gdat.truenormback[k, i]])
                else:
                    truepara = None
                indxbackpair = array([c, k])
                print 'hey'
                print 'gdat.lablback'
                print gdat.lablback
                strgpara = ['%s$_{%d}$' % (gdat.lablback[a], i) for a in indxbackpair]
                print 'indxbackpair'
                print indxbackpair
                print 'strgpara'
                print strgpara
                print 
                tdpy.mcmc.plot_grid(path, gdat.listnormback[:, indxbackpair, i], strgpara, truepara=truepara, join=True)
    
    # plot log-likelihood
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    figr.suptitle(r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi))
    axis.hist(listllik.flatten())
    axis.set_ylabel(r'$N_{samp}$')
    axis.set_xlabel(r'$\ln P(D|x)$')
    axis.axvline(amax(listllik), label='Maximum saved log-likelihood')
    axis.axvline(amax(listllik), label='Overall maximum log-likelihood')
    plt.subplots_adjust(top=0.87, bottom=0.1, left=0.1)
    figr.savefig(gdat.pathpost + 'histllik.pdf')
    plt.close(figr)

    # plot log-prior
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    figr.suptitle(r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi))
    if amin(listlpri) != amax(listlpri):
        axis.hist(listlpri.flatten())
        axis.set_ylabel(r'$N_{samp}$')
        axis.set_xlabel(r'$\ln P(x)$')
    plt.subplots_adjust(top=0.87, bottom=0.1, left=0.1)
    figr.savefig(gdat.pathpost + 'histlpri.pdf')
    plt.close(figr)

    # animate the frame plots
    if makeanim:
        make_anim(gdat)

    # plot resident memory
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxsamp, mean(gdat.listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(gdat.pathdiag + 'memoresi.pdf')
    plt.close(figr)

    timetotlfinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Plots are produced in %.3g seconds.' % (timetotlfinl - timetotlinit)


def plot_chro(gdat):

    if amin(gdat.listchrototl <= 0.) or amin(gdat.listchrollik <= 0.):
        print 'Invalid chronometer...'

    gdat.listchrototl *= 1e3
    binstime = logspace(log10(amin(gdat.listchrototl[where(gdat.listchrototl > 0)])), log10(amax(gdat.listchrototl)), 50)
    figr, axcl = plt.subplots(gdat.numbprop, 1, figsize=(2 * gdat.plotsize, gdat.numbprop * gdat.plotsize / 3.))
    for k in range(gdat.numbprop):
        indxswepchro = intersect1d(where(gdat.listindxprop == k)[0], where(gdat.listchrototl[:, 0] > 0)[0])
        try:
            axcl[k].hist(gdat.listchrototl[indxswepchro, 0], binstime, log=True, label=gdat.strgprop[k])
        except:
            if gdat.verbtype > 0:
                print 'Skipping proposal timing plot...'
        axcl[k].set_xlim([amin(binstime), amax(binstime)])
        axcl[k].set_ylim([0.5, None])
        axcl[k].set_ylabel(gdat.strgprop[k])
        axcl[k].set_xscale('log')
        if k != gdat.numbprop - 1:
            axcl[k].set_xticklabels([])
    axcl[-1].set_xlabel('$t$ [ms]')
    plt.subplots_adjust(hspace=0.05)
    figr.savefig(gdat.pathdiag + 'chroprop.pdf')
    plt.close(figr)

    labl = ['Total', 'Proposal', 'Prior', 'Likelihood', 'Save', 'Rest']
    listcolr = ['black', 'b', 'g', 'r', 'y', 'm']
    numblabl = len(labl)
    figr, axcl = plt.subplots(2, 1, figsize=(2 * gdat.plotsize, gdat.plotsize))
    for k in range(numblabl - 1):
        if k == numblabl - 2:
            varb = gdat.listchrototl[:, 0] - sum(gdat.listchrototl[:, 1:], 1)
        else:
            varb = gdat.listchrototl[:, k+1]
        try:
            axcl[0].hist(varb, binstime, log=True, label=labl[k+1], alpha=0.5, color=listcolr[k+1])
        except:
            pass
    axcl[1].hist(gdat.listchrototl[:, 0], binstime, log=True, label=labl[0], color=listcolr[0], alpha=0.5)
    axcl[1].set_title(r'$\langle t \rangle$ = %.3g ms' % mean(gdat.listchrototl[where(gdat.listchrototl[:, 0] > 0)[0], 0]))
    axcl[0].set_xlim([amin(binstime), amax(binstime)])
    axcl[1].set_xlim([amin(binstime), amax(binstime)])
    axcl[0].set_xscale('log')
    axcl[1].set_xscale('log')
    axcl[0].set_ylim([0.5, None])
    axcl[1].set_ylim([0.5, None])
    axcl[0].legend(loc=9, ncol=4)
    axcl[1].legend(loc=2)
    axcl[0].set_xticklabels([])
    axcl[1].set_xlabel('$t$ [ms]')
    plt.subplots_adjust(bottom=0.15, hspace=0.18)
    figr.savefig(gdat.pathdiag + 'chrototl.pdf')
    plt.close(figr)

    gdat.listchrollik *= 1e3
    if gdat.pntstype == 'lght':
        listlabl = ['Reading sample', 'Gathering pixels', 'Meshing pixels', 'PS flux map', 'Total Flux map', 'Counts', 'Likelihood']
    if gdat.pntstype == 'lens':
        listlabl = ['Solving LE', 'Counts', 'Likelihood']
    numblabl = len(listlabl)
    figr, axcl = plt.subplots(gdat.numbchrollik, 1, figsize=(2 * gdat.plotsize, gdat.plotsize * numblabl / 3.))
    maxmchrollik = amax(gdat.listchrollik)
    minmchrollik = amin(gdat.listchrollik[where(gdat.listchrollik > 0)])
    binstime = logspace(log10(minmchrollik), log10(maxmchrollik), 50)
    for k in range(gdat.numbchrollik):
        chro = gdat.listchrollik[where(gdat.listchrollik[:, k] > 0)[0], k]
        try:
            axcl[k].hist(chro, binstime, log=True, label=listlabl[k])
        except:
            print 'Skipping likelihood timing plot...'
        axcl[k].set_xlim([amin(binstime), amax(binstime)])
        axcl[k].set_ylim([0.5, None])
        axcl[k].set_ylabel(listlabl[k])
        axcl[k].set_xscale('log')
        if k != gdat.numbchrollik - 1:
            axcl[k].set_xticklabels([])
        axcl[k].axvline(mean(chro), ls='--', alpha=0.2, color='black')
    axcl[-1].set_xlabel('$t$ [ms]')
    plt.subplots_adjust(hspace=0.05)
    figr.savefig(gdat.pathdiag + 'chrollik.pdf')
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
            listydat[c+2, :] = gdat.postnormback[0, c, :] * gdat.backfluxmean[c, :]
            listyerr[:, c+2, :] = tdpy.util.retr_errrvarb(gdat.postnormback[:, c, :]) * gdat.backfluxmean[c, :]
    else:
        if gdat.correxpo:
            listydat[1, :] = sum(sum(gdatmodi.thispntsflux * gdat.expo, 1), 1) / sum(sum(gdat.expo, 1), 1)
        else:
            listydat[1, :] = sum(gdatmodi.thispntsflux)
            
        for c in gdat.indxback:
            listydat[c+2, :] = gdatmodi.thissampvarb[gdat.indxsampnormback[c, :]] * gdat.backfluxmean[c, :]
   
    xdat = gdat.meanener
    for k in range(gdat.numbback + 2):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [%s]' % gdat.strgenerunit)
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [%s/cm$^2$/s/sr]' % gdat.strgenerunit)
    axis.legend()

    if post:
        path = gdat.pathplot + 'compfracspec.pdf'
    else:
        path = gdat.pathplot + 'fram/compfracspec_swep%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
   
    listexpl = zeros(gdat.numbback + 1)
    listexpl[0] = 0.1

    listsize = zeros(gdat.numbback + 1)
    for k in range(gdat.numbback + 1):
        if gdat.numbener == 1:
            listsize[k] = gdat.diffener * listydat[k+1, :]
        else:
            listsize[k] = trapz(listydat[k+1, :], gdat.meanener)
    listsize *= 100. / sum(listsize)
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))

    labl = ['PS'] + gdat.lablback
    
    axis.pie(listsize, explode=listexpl, labels=labl, autopct='%1.1f%%')
    axis.axis('equal')

    if post:
        path = gdat.pathplot + 'compfrac.pdf'
    else:
        path = gdat.pathplot + 'fram/compfrac_swep%09d.pdf' % gdatmodi.cntrswep
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
        try:
            axis.hist(gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[l][gdatmodi.indxmodlpntscomp[l]]], gdat.binssind, alpha=gdat.alphmrkr, color='b', \
                                                                                                                                                log=True, label='Sample')
        except:
            print 'Skipping the plot of color histogram...'
            print 'gdat.binssind'
            print gdat.binssind
            print 'thissind'
            print gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[l][gdatmodi.indxmodlpntscomp[l]]]
    if gdat.trueinfo and gdat.indxtruepntscomp[l].size > 0:
        axis.hist(gdat.truespep[l][gdat.indxtruepntscomp[l], 0], gdat.binssind, alpha=gdat.alphmrkr, color='g', log=True, label=gdat.truelabl)
        if gdat.datatype == 'mock' and gdat.exprinfo:
            axis.hist(gdat.exprsind, gdat.binssind, alpha=gdat.alphmrkr, color='red', log=True, label='3FGL')
    axis.set_yscale('log')
    axis.set_xlabel('$s$')
    axis.set_xlim([gdat.minmsind, gdat.maxmsind])
    axis.set_ylabel('$N$')
    axis.set_ylim(gdat.limshist)
    axis.legend(loc=2)
    if post:
        path = gdat.pathpost + 'histsind_pop%d.pdf' % l
    else:
        path = gdat.pathfram + 'histsind_pop%d_swep%09d.pdf' % (l, gdatmodi.cntrswep)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    

def plot_brgt(gdat, gdatmodi=None, listspecbrgthist=None):
    
    if listspecbrgthist == None:
        post = False
    else:
        post = True
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if post:
        xdat = gdat.meansind
        postsindhist = tdpy.util.retr_postvarb(listspecbrgthist)
        ydat = postsindhist[0, :]
        yerr = tdpy.util.retr_errrvarb(postsindhist)
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    else:   
        fluxbrgt, fluxbrgtassc = retr_fluxbrgt(gdat, gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsamplgal)], \
                                                                         gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampbgal)], \
                                                                         gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec)[gdat.indxenerfluxdist[0], :]])
        print 'hey'
        print 'fluxbrgt'
        print fluxbrgt
        print 'fluxbrgtassc'
        print fluxbrgtassc
        axis.scatter(fluxbrgt, fluxbrgtassc, alpha=gdat.alphmrkr, color='b', label='Sample')
        axis.scatter(fluxbrgt[0], sum(fluxbrgtassc), alpha=gdat.alphmrkr, color='b', label='Sample - Total')
    if gdat.trueinfo:
        axis.scatter(gdat.truefluxbrgt, gdat.truefluxbrgtassc, alpha=gdat.alphmrkr, color='g', label=gdat.truelabl)
        if gdat.datatype == 'mock' and gdat.exprinfo:
            axis.scatter(gdat.exprfluxbrgt, gdat.exprfluxbrgtassc, alpha=gdat.alphmrkr, color='red', label=gdat.strgcatl)
    axis.set_xscale('log')
    axis.set_xlabel(r'$%s_{assc}$' % gdat.strgflux)
    axis.set_xlim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_ylim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_ylabel(r'$%s_{max}$' % gdat.strgflux)
    axis.legend(loc=2)
    if post:
        path = gdat.pathpost + 'scatbrgt.pdf'
    else:
        path = gdat.pathfram + 'scatbrgt_swep%09d.pdf' % (gdatmodi.cntrswep)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    

def plot_fluxsind(gdat, l, strgtype='scat', gdatmodi=None, listspechist=None, listsindhist=None):
    
    if listsindhist == None:
        post = False
    else:
        post = True
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if post:
        postsindhist = tdpy.util.retr_postvarb(listsindhist)
        ydat = postsindhist[0, :]
        yerr = tdpy.util.retr_errrvarb(postsindhist)
        #axis.hist2d(flux, sind, gdat.binsflux, gdat.binssind, color='b', alpha=gdat.alphmrkr)
        sns.kdeplot(flux, sind, ax=axis, cmap="Blues")

    else:
        flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]]
        sind = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[l]]
        if strgtype == 'hist':
            #hist = histogram2d(flux, sind, bins=[gdat.binsflux, gdat.binssind])[0]
            #axis.pcolor(gdat.binsflux, gdat.binssind, hist, cmap='Blues', label='Sample',  alpha=gdat.alphmrkr)
            try:
                sns.kdeplot(flux, sind, ax=axis, cmap='Blues', label='Sample', legend=True)
            except:
                print 'Skipping flux-color scatter plot...'
        else:
            axis.scatter(flux, sind, alpha=gdat.alphmrkr, color='b', label='Sample')
    
    if gdat.trueinfo:
        flux = gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]
        if strgtype == 'hist':
            #hist = histogram2d(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truespep[l][:, 0], bins=[gdat.binsflux, gdat.binssind], )[0]
            #axis.pcolor(gdat.binsflux, gdat.binssind, hist, cmap='Greens', label=gdat.truelabl,  alpha=gdat.alphmrkr)
            sns.kdeplot(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truespep[l][:, 0], ax=axis, cmap='Greens', label=gdat.truelabl)
        else:
            axis.scatter(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :], gdat.truespep[l][:, 0], alpha=gdat.alphmrkr, color='g', label=gdat.truelabl)
        if gdat.datatype == 'mock' and gdat.exprinfo:
            if strgtype == 'hist':
                #hist = histogram2d(gdat.exprspec[0, gdat.indxenerfluxdist[0], :], gdat.exprsind, bins=[gdat.binsflux, gdat.binssind])[0]
                #axis.pcolor(gdat.binsflux, gdat.binssind, hist, color='Reds', label=gdat.nameexpr, alpha=gdat.alphmrkr)
                sns.kdeplot(gdat.exprspec[0, gdat.indxenerfluxdist[0], :], gdat.exprsind, ax=axis, cmap='Reds', label=gdat.nameexpr)
            else:
                axis.scatter(gdat.exprspec[0, gdat.indxenerfluxdist[0], :], gdat.exprsind, alpha=gdat.alphmrkr, color='r', label=gdat.nameexpr)
    axis.set_xscale('log')
    axis.set_xlabel('$%s$ [%s]' % (gdat.strgflux, gdat.strgfluxunit))
    axis.set_ylabel('$s$')
    axis.set_xlim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_ylim([gdat.minmsind, gdat.maxmsind])
    axis.legend(loc=2)
    if post:
        path = gdat.pathplot + '%sfluxsind_pop%d' % (strgtype, l) + '.pdf'
    else:
        path = gdat.pathplot + 'fram/%sfluxsind_pop%d' % (strgtype, l) + '_swep%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    

def plot_histspec(gdat, l, gdatmodi=None, plotspec=False, listspechist=None):
  
    if listspechist == None:
        post = False
    else:
        post = True
    
    if gdat.pixltype != 'unbd' and gdat.pntstype == 'lght':
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
        indxenertemp = gdat.indxenerfluxdist
        axcltemp = []
        for i in gdat.indxener:
            if i == gdat.indxenerfluxdist:
                axcltemp.append(axcl)
            else:
                axcltemp.append(None)
        axcl = axcltemp
    
    for i in indxenertemp:
        axis = axcl[i]
        if post:
            xdat = gdat.factfluxconv * gdat.meanspecplot[i, :]
            postspechist = tdpy.util.retr_postvarb(listspechist[:, :, i])
            ydat = postspechist[0, :]
            yerr = tdpy.util.retr_errrvarb(postspechist)
            axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        else:
            spec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]][i, gdatmodi.indxmodlpntscomp[l]]
            # temp -- gdat.binsspecplot may be too narrow if there are many energy bins or PS colors are extreme
            try:
                axis.hist(gdat.factfluxconv * spec, gdat.factfluxconv * gdat.binsspecplot[i, :], alpha=gdat.alphmrkr, color='b', log=True, label='Sample')
            except:
                print 'Spectral bins are inappropriate. Skipping the spectral histogram...'
            
            # superimpose the current prior flux distribution
            if i == gdat.indxenerfluxdist[0]:
                meanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]]
                if gdat.fluxdisttype[l] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]]  
                    fluxhistmodl = meanpnts * pdfn_flux_powr(gdat, gdat.meanfluxplot, fluxdistslop) * gdat.difffluxplot
                if gdat.fluxdisttype[l] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]]  
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]]  
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]]  
                    fluxhistmodl = meanpnts * pdfn_flux_brok(gdat, gdat.meanfluxplot, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr) * gdat.difffluxplot
                
                axis.plot(gdat.factfluxconv * gdat.meanspecplot[i, :], fluxhistmodl, ls='--', alpha=gdat.alphmrkr, color='b')

        if gdat.pixltype != 'unbd' and gdat.pntstype == 'lght':
            # add another horizontal axis for counts
            axiscnts = axis.twiny()
            axiscnts.set_xscale('log')
            axiscnts.set_xlabel('$C$')
            axiscnts.spines['bottom'].set_position(('axes', 1.))
            axiscnts.set_xlim([gdat.binscnts[i, 0], gdat.binscnts[i, -1]])
            axiscnts.xaxis.set_ticks_position('bottom')
            axiscnts.xaxis.set_label_position('bottom')
           
            # temp
            if False:
                # add yet another horizontal axis for fluctuation significance
                axissigm = axis.twiny()
                axissigm.set_xscale('log')
                axissigm.set_xlabel(r'$\sigma$')
                axissigm.spines['top'].set_position(('axes', 1.05))
                axissigm.set_xlim([binssigm[i, 0], binssigm[i, -1]])
                axissigm.axvline(1., ls='--', alpha=0.1)
                axissigm.axvline(5., ls='--', alpha=0.1)
    
        # superimpose the true catalog
        if gdat.trueinfo and gdat.indxtruepntscomp[l].size > 0:
            truehist = axis.hist(gdat.factfluxconv * gdat.truespec[l][0, i, gdat.indxtruepntscomp[l]], gdat.factfluxconv * gdat.binsspecplot[i, :], alpha=gdat.alphmrkr, \
                                                                                                                            color='g', log=True, label=gdat.truelabl)
            if gdat.datatype == 'mock' and gdat.exprinfo:
                if gdat.exprtype == 'ferm':
                    axis.hist(gdat.factfluxconv * gdat.exprspec[0, i, :], gdat.factfluxconv * gdat.binsspecplot[i, :], color='red', alpha=gdat.alphmrkr, log=True, label='3FGL')
        
        axis.set_yscale('log')
        axis.set_xlabel('$%s$ [%s]' % (gdat.strgflux, gdat.strgfluxunit))
        axis.set_xscale('log')
        axis.text(0.75, 0.65, gdat.strgbinsener[i], ha='center', va='center', transform=axis.transAxes)
        axis.set_ylim(gdat.limshist)
        axis.set_xlim([gdat.factfluxconv * gdat.minmspecplot[i], gdat.factfluxconv * gdat.maxmspecplot[i]])
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
        path = gdat.pathpost + 'hist%s_pop%d' % (strg, l) + '.pdf'
    else:
        path = gdat.pathfram + 'hist%s_pop%d' % (strg, l) + '_swep%09d.pdf' % gdatmodi.cntrswep
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    

def plot_scatspec(gdat, l, gdatmodi=None, postspecmtch=None, plotdiff=False):
    
    if postspecmtch != None:
        post = True
    else:
        post = False

    figr, axrw = plt.subplots(1, gdat.numbener, figsize=(gdat.plotsize * gdat.numbener, gdat.plotsize))
    if gdat.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        # prepare data to be plotted
        xdat = gdat.truespec[l][0, i, :]
        xerr = tdpy.util.retr_errrvarb(gdat.truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
        if post:
            yerr[0, :] = postspecmtch[0, i, :]
            ydat = postspecmtch[1, i, :]
            yerr[1, :] = postspecmtch[2, i, :]
            yerr = fabs(yerr - ydat)
        else:
            ydat = gdatmodi.thisspecmtch[i, :]

        xdat *= gdat.factfluxconv
        xerr *= gdat.factfluxconv
        ydat *= gdat.factfluxconv
        yerr *= gdat.factfluxconv
        
        if plotdiff:
            ydat = 100. * (ydat - xdat) / xdat
       
        # plot all associations
        indx = where(ydat > 0.)[0]
        if indx.size > 0:
            axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], \
                                                                                                    lw=1, marker='o', markersize=5, color='black', alpha=0.1)
       
        # plot associations inside the comparison area
        indx = intersect1d(where(ydat > 0.)[0], gdat.indxtruepntscomp[l])
        if indx.size > 0:
            axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], \
                                                                                                    lw=1, marker='o', markersize=5, color='black')
        
        # plot associations to multiple model point sources
        if not post:
            indx = intersect1d(gdatmodi.indxtruepntsassc[l].mult, gdat.indxtruepntscomp[l])
            if len(indx) > 0:
                axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], \
                                                                                                                        lw=1, marker='o', markersize=5, color='r')
    
        if plotdiff:
            axis.axhline(0., ls='--', alpha=gdat.alphmrkr, color='black')
        else:
            # superimpose the bias line
            fluxbias = retr_fluxbias(gdat, gdat.binsspecplot[i, :], i)
            axis.plot(gdat.factfluxconv * gdat.binsspecplot[i, :], gdat.factfluxconv * gdat.binsspecplot[i, :], ls='--', alpha=gdat.alphmrkr, color='black')
            axis.plot(gdat.factfluxconv * gdat.binsspecplot[i, :], gdat.factfluxconv * fluxbias[0, :], ls='--', alpha=gdat.alphmrkr, color='black')
            axis.plot(gdat.factfluxconv * gdat.binsspecplot[i, :], gdat.factfluxconv * fluxbias[1, :], ls='--', alpha=gdat.alphmrkr, color='black')
        
        if gdat.datatype == 'mock':
            strg = 'true'
        else:
            if gdat.exprtype == 'ferm':
                strg = '3FGL'
        axis.set_xlabel('$%s_{%s}$ [%s]' % (gdat.strgflux, strg, gdat.strgfluxunit))
        if i == 0:
            axis.set_ylabel('$%s_{samp}$ [%s]' % (gdat.strgflux, gdat.strgfluxunit))
        if i == gdat.numbener / 2:
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='black', alpha=0.1)
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='black')
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='r')
            axis.legend(loc=2)
        if not plotdiff:
            axis.set_yscale('log')
        axis.set_xscale('log')
        axis.set_title(gdat.strgbinsener[i])
        if plotdiff:
            limsyaxi = array([-100., 100.])
        else:
            limsyaxi = gdat.factfluxconv * array([gdat.minmspecplot[i], gdat.maxmspecplot[i]])
        axis.set_ylim(limsyaxi)
        axis.set_xlim([gdat.factfluxconv * gdat.minmspecplot[i], gdat.factfluxconv * gdat.maxmspecplot[i]])
        axis.set_title(gdat.strgbinsener[i])
   
    if plotdiff:
        strg = 'diff'
    else:
        strg = ''
    if postspecmtch != None:
        path = gdat.pathpost + 'scatspec%s_pop%d' % (strg, l) + '.pdf'
    elif gdatmodi.thisspecmtch != None:
        path = gdat.pathfram + 'scatspec%s_pop%d' % (strg, l) + '_swep%09d.pdf' % gdatmodi.cntrswep

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
            axis.scatter(gdat.datacnts[i, :, m], gdatmodi.thismodlcnts[i, :, m], alpha=gdat.alphmrkr)

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
                axis.set_title(gdat.strgbinsener[i])
            
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'fram/scatpixl_pop%d' % l + '_swep%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 41
    numbfluxprox = len(gdat.indxpixlprox)
    bins = empty((numbfluxprox, numbbins))
    indxpixlproxsize = empty((numbfluxprox, gdat.numbpixl))
    for h in range(gdat.numbfluxprox):
        for j in gdat.indxpixl:
            indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
        bins[h, :] = logspace(log10(amin(indxpixlproxsize[h, :])), log10(amax(indxpixlproxsize[h, :])), numbbins)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbfluxprox):
        axis.hist(indxpixlproxsize[h, :], bins=bins[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphmrkr)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixl, label='ROI', ls='--')
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    axis.legend(loc=2)
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'init/indxprox.pdf')
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
    plt.savefig(gdat.pathplot + 'evidtest.pdf')
    plt.close(figr)
    
    
def plot_pntsprob(gdat, ptag, full=False, cumu=False):
    
    if cumu:
        numbrows = 1
        numbcols = 1
    else:
        numbcols = 2
        if full:
            numbrows = gdat.numbfluxplot / 2
        else:
            numbrows = 2
        
    for l in gdat.indxpopl:
        figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
        if numbrows == 1:
            axgr = [axgr]            
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
                    indxuppr = gdat.numbfluxplot
                else:
                    if h < 3:
                        indxlowr = 2 * h
                        indxuppr = 2 * (h + 1)
                    else:
                        indxlowr = 2 * h
                        indxuppr = gdat.numbfluxplot
                temp = sum(gdat.pntsprob[l, :, :, indxlowr:indxuppr], 2)
                if where(temp > 0.)[0].size > 0:
                    imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
                else:
                    imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                    
                # superimpose true PS
                # temp
                if gdat.trueinfo:
                    indxpnts = where((gdat.binsspecplot[gdat.indxenerfluxdist[0], indxlowr] < gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]) & \
                        (gdat.truespec[l][0, gdat.indxenerfluxdist[0], :] < gdat.binsspecplot[gdat.indxenerfluxdist[0], indxuppr]))[0]
                    mrkrsize = retr_mrkrsize(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist[0], indxpnts])
                    axis.scatter(gdat.anglfact * gdat.truelgal[l][indxpnts], gdat.anglfact * gdat.truebgal[l][indxpnts], \
                                                                                            s=mrkrsize, alpha=gdat.alphmrkr, marker='*', lw=2, color='g')
                
                if a == numbrows - 1:
                    axis.set_xlabel(gdat.strgxaxi)
                else:
                    axis.set_xticklabels([])
                if b == 0:
                    axis.set_ylabel(gdat.strgyaxi)
                else:
                    axis.set_yticklabels([])

                draw_frambndr(gdat, axis)
                
                titl = tdpy.util.mexp(gdat.binsspecplot[gdat.indxenerfluxdist, indxlowr]) + ' $< %s <$' % gdat.strgflux + \
                                                                                            tdpy.util.mexp(gdat.binsspecplot[gdat.indxenerfluxdist, indxuppr])
                axis.set_title(titl)
        
        plt.figtext(0.5, 0.97, '$%s$ [%s]' % (gdat.strgflux, gdat.strgfluxunit), ha='center', va='center')
        axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
        cbar = figr.colorbar(imag, cax=axiscomm)
        #cbar.set_ticks(gdat.tickdatacnts[i, :])
        #cbar.set_ticklabels(gdat.labldatacnts[i, :])
        plt.subplots_adjust(left=0.1, top=.92, hspace=0.05, wspace=0.03, bottom=0.08)
        figr.savefig(gdat.pathpost + 'pntsbind' + ptag + '%d%d' % (l, gdat.indxenerincl[gdat.indxenerfluxdist]) + '.pdf')
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
    if gdat.modlvarioaxi:
        psfntemp = gdatmodi.thispsfn[:, :, :, 0]
    else:
        psfntemp = gdatmodi.thispsfn

    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            # angular range to be plotted
            indxangltemp = where(psfntemp[i, :, m] > 1e-6 * amax(psfntemp[i, :, m]))[0]
            
            axis.plot(gdat.binsanglplot[indxangltemp], psfntemp[i, indxangltemp, m], label='Sample')
            
            if gdat.truepsfn != None:
                if gdat.truevarioaxi:
                    truepsfntemp = gdat.truepsfn[i, indxangltemp, m, 0]
                    truepsfntotltemp = gdat.truepsfn[i, :, m, 0]
                else:
                    truepsfntemp = gdat.truepsfn[i, indxangltemp, m]
                    truepsfntotltemp = gdat.truepsfn[i, :, m]
                axis.plot(gdat.binsanglplot[indxangltemp], truepsfntemp, label=gdat.nameexpr, color='g', ls='--')
                axis.set_xlim([0., gdat.binsanglplot[amax(where(truepsfntotltemp > 1e-6 * amax(truepsfntotltemp))[0])]])
            axis.set_yscale('log')
                
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$\theta$ [%s]' % gdat.strganglunit)
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
            if m == 0:
                axis.set_title(gdat.strgbinsener[i])  
            if i == gdat.numbener - 1 and m == gdat.numbevtt - 1:
                axis.legend()
  
            indxsamptemp = gdat.indxsamppsfp[0] + i * gdat.numbpsfptotl + m * gdat.numbener * gdat.numbpsfptotl
            if gdat.modlpsfntype == 'singgaus':
                strg = r'$\sigma = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp]) + gdat.strganglunit
                # temp
            elif gdat.modlpsfntype == 'singking':
                strg = r'$\sigma = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp]) + gdat.strganglunit + '\n'
                strg += r'$\gamma = %.3g$' % gdatmodi.thissampvarb[indxsamptemp+1]
            elif gdat.modlpsfntype == 'doubgaus':
                strg = r'$f = %.3g$' % gdatmodi.thissampvarb[indxsamptemp] + '\n'
                strg += r'$\sigma = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+1]) + gdat.strganglunit + '\n'
                strg += r'$\sigma = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+2]) + gdat.strganglunit
            elif gdat.modlpsfntype == 'gausking':
                strg = r'$f_G = %.3g$' % gdatmodi.thissampvarb[indxsamptemp] + '\n'
                strg += r'$\sigma_G = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+1]) + gdat.strganglunit + '\n'
                strg += r'$\sigma_K = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+2]) + gdat.strganglunit + '\n'
                strg += r'$\gamma = %.3g$' % gdatmodi.thissampvarb[indxsamptemp+3]
            elif gdat.modlpsfntype == 'doubking':
                strg = r'$f_c = %.3g$' % gdatmodi.thissampvarb[indxsamptemp] + '\n'
                strg += r'$\sigma_c = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+1]) + gdat.strganglunit + '\n'
                strg += r'$\gamma_c = %.3g$' % gdatmodi.thissampvarb[indxsamptemp+2] + '\n'
                strg += r'$\sigma_t = %.3g$ ' % (gdat.anglfact * gdatmodi.thissampvarb[indxsamptemp+3]) + gdat.strganglunit + '\n'
                strg += r'$\gamma_t = %.3g$' % gdatmodi.thissampvarb[indxsamptemp+4]
            axis.text(0.2, 0.2, strg, va='center', ha='center', transform=axis.transAxes, fontsize=18)
            
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'fram/psfnprof_swep%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    
    
def plot_fwhm(gdat, gdatmodi):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))

    tranfwhm = transpose(gdatmodi.thisfwhm)
    imag = plt.pcolor(gdat.enerfact * gdat.meanener, gdat.indxevtt, rad2deg(tranfwhm), cmap='BuPu')
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.set_xscale('log')
    axis.set_ylabel('PSF Class')
    axis.set_xlabel(r'$E$ [%s]' % gdat.strgenerunit)
    axis.set_yticks(gdat.indxevtt + 0.5)
    axis.set_yticklabels(['%d' % m for m in gdat.evttfull[gdat.indxevttincl]])
    axis.set_xticks(gdat.binsener * gdat.enerfact)
    axis.set_xticklabels(['%.2g' % (gdat.binsener[i] * gdat.enerfact) for i in arange(gdat.numbener + 1)])
    axis.set_xlim([gdat.binsener[0] * gdat.enerfact, gdat.binsener[-1] * gdat.enerfact])
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            axis.text(gdat.meanener[i], gdat.indxevttincl[m] + 0.5, r'$%.3g^\circ$' % rad2deg(tranfwhm[m, i]), ha='center', va='center', fontsize=14)

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'fram/fwhmcnts_swep%09d.pdf' % gdatmodi.cntrswep)
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
                axis.set_title(gdat.strgbinsener[i])
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
        
    plt.tight_layout()
    plt.savefig(gdat.pathinit + 'datacntshist.pdf')
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

    if gdat.truevarioaxi:
        psfntemp = copy(gdat.truepsfn[0, :, 0, 0])
    else:
        psfntemp = copy(gdat.truepsfn[0, :, 0])
    psfntemp /= amax(psfntemp)
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
        axis.plot(gdat.binsanglplot, gdat.binsfluxprox[k] * psfntemp, label=labl, color=colr, alpha=alph)
        axis.set_xlim([amin(gdat.binsanglplot), amax(gdat.binsanglplot)])
        if k > 0:
            axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
    axis.set_yscale('log')
    axis.set_xlabel(r'$\theta$ [%s]' % gdat.strganglunit)
    axis.set_ylabel('$%s$ [%s]' % (gdat.strgflux, gdat.strgfluxunit))

    limt = gdat.specfraceval * amax(gdat.binsfluxprox[0] * psfntemp)
    maxmangltemp = interp(1e-1 * limt, gdat.binsfluxprox[k] * psfntemp[::-1], gdat.binsanglplot[::-1])
    
    if limt > 0:
        axis.axhline(limt, color='red', ls=':', label='Flux floor')
    axis.set_xlim([None, maxmangltemp])
    axis.set_ylim([1e-3 * limt, None])
    legd = axis.legend(frameon=True, shadow=True, fancybox=True, framealpha=1.)
    legd.get_frame().set_facecolor('white')
    
    plt.tight_layout()
    plt.savefig(gdat.pathinit + 'eval.pdf')
    plt.close(figr)


def plot_mosa(gdat):

    # empty global object
    gdatmodi = tdpy.util.gdatstrt()
    
    # data structure to hold the indices of model PS to be compared to the reference catalog 
    gdatmodi.indxmodlpntscomp = [[] for l in gdat.indxpopl]
    
    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamp:
        indxsampmosa = choice(gdat.indxsamp, size=numbsampmosa, replace=False)
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            gdatmodi.thissampvarb = gdat.listsampvarb[indxsampmosa[numbcols*a+b], :].flatten()
                            gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampspec, \
                                                            gdatmodi.thisindxsampspep, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdat.listindxpntsfull[l])
                            gdatmodi.indxmodlpntscomp[l] = retr_indxpntscomp(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]], \
                                                                                                            gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]])
                            # cross-correlate with the reference catalog
                            if gdat.trueinfo:
                                gdatmodi.indxtruepntsassc = []
                                indxmodl, indxtruepntsassc = corr_catl(gdat, l, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]], \
                                                                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]])
                                gdatmodi.indxtruepntsassc.append(indxtruepntsassc)
         
                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.strgxaxi)
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.strgyaxi)
                            else:
                                axis.set_yticklabels([])
                            
                            imag = retr_imag(gdat, axis, gdat.datacnts, i, m, vmin=gdat.minmdatacnts[i], vmax=gdat.maxmdatacnts[i], scal=gdat.scalmaps)
                            supr_fram(gdat, gdatmodi, axis, i, m)
                   
                    plt.figtext(0.5, 0.93, gdat.strgbinsener[i], ha='center', va='center')
                    axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
                    cbar = figr.colorbar(imag, cax=axiscomm)
                    cbar.set_ticks(gdat.tickdatacnts[i, :])
                    cbar.set_ticklabels(gdat.labldatacnts[i, :])
                    plt.subplots_adjust(left=0.1, top=.91, hspace=0.01, wspace=0.03, bottom=0.09)
                    if l == 1:
                        strg = ''
                    else:
                        strg = '_pop%d' % l
                    if m == None:
                        path = gdat.pathpost + 'mosa' + strg + '_%dA.pdf' % (gdat.indxenerincl[i])
                    else:
                        path = gdat.pathpost + 'mosa' + strg + '_%d%d.pdf' % (gdat.indxenerincl[i], gdat.indxevttincl[m])
                    plt.savefig(path)
                    plt.close(figr)
    else:
        if gdat.verbtype > 0:
            print 'Skipping the mosaic image plot...'


def plot_grap(igal=False):
        
    figr, axis = plt.subplots(figsize=(6, 6))

    grap = nx.DiGraph()
    if igal:
        listcolr = ['black', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'black', 'black']
    else:
        listcolr = ['black', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'olive', 'magenta']
    grap.add_edges_from([ \
                         ('fluxdistslop', 'flux'), \
                         ('meanpnts', 'numbpnts'), \
                         ('modl','data'), \
                         ('psfp', 'modl'), \
                         ('normback', 'modl'), \
                         ('lgal','modl'), \
                         ('bgal','modl'), \
                         ('flux','modl'), \
                         ('numbpnts','lgal'), \
                         ('numbpnts','bgal'), \
                         ('numbpnts','flux') \
                        ])
    if igal:
        grap.add_edges_from([ \
                             ('spatdistslop', 'lgal'), \
                             ('spatdistslop', 'bgal') \
                            ])
    else:
        grap.add_edges_from([ \
                             ('sind', 'modl'), \
                             ('numbpnts', 'sind') \
                            ])

        
    labl = {}
    labl['fluxdistslop'] = r'$\alpha$'
    labl['meanpnts'] = r'$\mu$'
    labl['data'] = r'$\mathcal{D}$'
    labl['numbpnts'] = '$N$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['normback'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{l}$'
    labl['bgal'] = r'$\vec{b}$'
    if igal:
        labl['spatdistslop'] = r'$\beta$'
        labl['flux'] = r'$\vec{f_i}$'
    else:
        labl['flux'] = r'$\vec{f}$'
        labl['sind'] = r'$\vec{s}$'
    labl['modl'] = r'$\mathcal{M}$'
    
        
    
    posi = nx.circular_layout(grap)
    posi['meanpnts'] = array([0., 0.15])
    posi['numbpnts'] = array([0., 0.075])
    posi['normback'] = array([0.5, -0.0])
    posi['psfp'] = array([0.7, -0.0])
    posi['fluxdistslop'] = array([0.3, 0.15])
    posi['lgal'] = array([-0.3, -0.0])
    posi['bgal'] = array([-0.1, -0.0])
    posi['flux'] = array([0.1, -0.0])
    if igal:
        posi['spatdistslop'] = array([-0.2, 0.0])
    else:
        posi['sind'] = array([0.3, -0.0])
    posi['modl'] = array([0., -0.075])
    posi['data'] = array([0., -0.15])
    
    size = 500
    nx.draw(grap, posi, labels=labl, ax=axis, edgelist=[], nodelist=[])
    nx.draw_networkx_edges(grap, posi, ax=axis, labels=labl, edge_color=listcolr)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['modl', 'data'], node_color='grey', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbpnts'], node_color='b', node_size=size)
    if igal:
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'fluxdistslop', 'spatdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'flux'], node_color='g', node_size=size)
    else:
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'fluxdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'flux', 'sind'], node_color='g', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['psfp', 'normback'], node_color='y', node_size=size)
    
    pathplot = os.environ["PCAT_DATA_PATH"] + '/imag/'
    #plt.subplots_adjust(top)
    plt.tight_layout()
    if igal:
        strg = 'igal'
    else:
        strg = ''
    plt.savefig(pathplot + 'grap%s.pdf' % strg)
    plt.close(figr)


#plot_grap(True)
plot_grap()
def plot_3fgl_thrs(gdat):

    path = os.environ["PCAT_DATA_PATH"] + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = pf.getdata(path, 0)

    bgalfgl3 = linspace(-90., 90., 481)
    lgalfgl3 = linspace(-180., 180., 960)

    bgalexpo = linspace(-90., 90., 400)
    lgalexpo = linspace(-180., 180., 800)

    #fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)
    fluxthrs = griddata([lgalfgl3, bgalfgl3], fluxthrs, [gdat.lgalheal])

    cntsthrs = fluxthrs * gdat.expo

    jbgal = where(abs(bgalexpo) < 10.)[0]
    jlgal = where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel(gdat.strgxaxi)
    axis.set_ylabel(gdat.strgyaxi)

    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'thrs.pdf')
    plt.close(figr)
    

def make_anim(gdat):

    listname = ['psfnprof', 'compfrac', 'compfracspec', 'scatpixl']
    listnamepopl = ['histspec', 'histflux', 'scatfluxsind', 'histfluxsind', 'histcnts', 'histsind', 'scatspec']
    listnamepoplenerevtt = ['datacnts', 'resicnts', 'modlcnts']
    listnameenerevtt = ['factoaxi', 'errr', 'errrcnts']
    for l in gdat.indxpopl:
        for k in range(len(listnamepopl)):
            listname.append('%s_pop%d' % (listnamepopl[k], l))
        for i in gdat.indxener:
            for m in gdat.indxevttplot:
                if m == None:
                    strg = 'A'
                else:
                    strg = '%d' % m
                for k in range(len(listnamepoplenerevtt)):
                    listname.append('%s_pop%d_%d%s' % (listnamepoplenerevtt[k], l, i, strg))
                for k in range(len(listnameenerevtt)):
                    listname.append('%s_%d%s' % (listnameenerevtt[k], i, strg))

    pathanim = gdat.pathplot + 'anim/'

    os.system('mkdir -p %s' % pathanim)
    if gdat.verbtype > 0:
        print 'Making animations...'
    for name in listname:
    
        strg = '%s_swep*.pdf' % name
        listfile = fnmatch.filter(os.listdir(gdat.pathfram), strg)[int(gdat.numbburn / gdat.numbswepplot):]
        
        numbfile = len(listfile)
        if numbfile == 0:
            if gdat.verbtype > 0:
                print 'Skipping animation for %s...' % name
            continue
        else:
            if gdat.verbtype > 0:
                print 'Producing animation for %s...' % name
            
        indxfileanim = choice(arange(numbfile), replace=False, size=numbfile)

        cmnd = 'convert -delay 20 -density 300 -quality 100 '
        for k in range(numbfile):
            cmnd += '%s%s ' % (gdat.pathfram, listfile[indxfileanim[k]])
        cmnd += ' %s%s.gif' % (pathanim, name)
        os.system(cmnd)

       
def plot_factoaxi(gdat, i, m, gdatmodi=None):

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if gdat.modlvarioaxi:
        if gdatmodi == None:
            psfpoaxinorm = gdat.postpsfpoaxinorm
            psfpoaxiindx = gdat.postpsfpoaxiindx
        else:
            psfpoaxinorm = gdatmodi.thissampvarb[gdat.indxsamppsfpoaxinorm]
            psfpoaxiindx = gdatmodi.thissampvarb[gdat.indxsamppsfpoaxiindx]
        factoaxi = retr_factoaxi(gdat, gdat.binsoaxi, psfpoaxinorm, psfpoaxiindx)
        axis.plot(gdat.binsoaxi * gdat.anglfact, factoaxi[i, :], label='Sample', color='b')
    
    if gdat.truefactoaxi != None:
        axis.plot(gdat.binsoaxi * gdat.anglfact, gdat.truefactoaxi[i, :], label=gdat.truelabl, color='g', ls='--')
        axis.set_ylim([0., 2. * amax(gdat.truefactoaxi[i, :])])
    axis.set_xlabel(r'$\phi$ [%s]' % gdat.strganglunit)
    axis.set_ylabel(r'f($\phi$)')
    axis.legend(loc=2)
    plt.tight_layout()
    if gdatmodi == None:
        strg = 'post/factoaxi_%d%d.pdf' % (i, m)
    else:
        strg = 'fram/factoaxi_%d%d_swep%09d.pdf' % (i, m, gdatmodi.cntrswep)
    plt.savefig(gdat.pathplot + strg)
    plt.close(figr)
    

def plot_histcnts(gdat, l, gdatmodi=None):

    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            try:
                axis.hist(gdatmodi.thiscnts[l][i, :, m], gdat.binscnts[i, :], alpha=gdat.alphmrkr, color='b', log=True, label='Sample')
            except:
                if gdat.verbtype > 0:
                    print 'Skipping PS count histogram plot...'
            if gdat.trueinfo:
                try:
                    truehist = axis.hist(gdat.truecnts[l][i, :, m], gdat.binscnts[i, :], alpha=gdat.alphmrkr, color='g', log=True, label=gdat.truelabl)
                except:
                    if gdat.verbtype > 0:
                        print 'Skipping PS count histogram plot...'
                if gdat.datatype == 'mock' and gdat.exprinfo:
                    if gdat.exprtype == 'ferm':
                        axis.hist(gdat.exprcnts[i, :, m], gdat.binscnts[i, :], alpha=gdat.alphmrkr, color='red', log=True, label='3FGL')
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            axis.set_xscale('log')
            axis.set_ylim(gdat.limshist)
            if m == 0:
                axis.set_title(gdat.strgbinsener[i])
            if i == 0 and gdat.exprtype == 'ferm':
                axis.set_ylabel(gdat.evttstrg[m])
            if m == gdat.numbevtt / 2 and i == gdat.numbener / 2:
                axis.legend()
        
    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'fram/histcnts_pop%d' % l + '_swep%09d.pdf' % gdatmodi.cntrswep)
    plt.close(figr)
    

def plot_datacnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_figr(gdat, indxenerplot, indxevttplot, 'datacnts', indxpoplplot=indxpoplplot, gdatmodi=gdatmodi)
    make_catllabl(gdat, axis)
    if gdat.correxpo:
        imag = retr_imag(gdat, axis, gdat.datacnts, indxenerplot, indxevttplot, vmin=gdat.minmdatacnts[indxenerplot], vmax=gdat.maxmdatacnts[indxenerplot], scal=gdat.scalmaps)
        make_cbar(gdat, axis, imag, indxenerplot, tick=gdat.tickdatacnts[indxenerplot, :], labl=gdat.labldatacnts[indxenerplot, :])
    else:
        imag = retr_scat(gdat, axis, gdat.datacnts, indxenerplot, indxevttplot)
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_modlcnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_figr(gdat, indxenerplot, indxevttplot, 'modlcnts', indxpoplplot=indxpoplplot, gdatmodi=gdatmodi)
    make_catllabl(gdat, axis)
    if gdat.pixltype == 'unbd':
        modltemp = gdatmodi.thismodlflux
    else:
        modltemp = gdatmodi.thismodlcnts
    imag = retr_imag(gdat, axis, modltemp, indxenerplot, indxevttplot, vmin=gdat.minmdatacnts[indxenerplot], \
                                                                                                                vmax=gdat.maxmdatacnts[indxenerplot], scal=gdat.scalmaps)
    make_cbar(gdat, axis, imag, indxenerplot, tick=gdat.tickdatacnts[indxenerplot, :], labl=gdat.labldatacnts[indxenerplot, :])
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_resicnts(gdat, indxpoplplot, gdatmodi, indxenerplot, indxevttplot):

    figr, axis, path = init_figr(gdat, indxenerplot, indxevttplot, 'resicnts', indxpoplplot=indxpoplplot, gdatmodi=gdatmodi)
    make_catllabl(gdat, axis)
    imag = retr_imag(gdat, axis, gdatmodi.thisresicnts, indxenerplot, indxevttplot, vmax=gdat.maxmresicnts[indxenerplot], cmap='RdBu', scal=gdat.scalmaps)
    supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot)
    make_cbar(gdat, axis, imag, indxenerplot, tick=gdat.tickresicnts[indxenerplot, :], labl=gdat.lablresicnts[indxenerplot, :])
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
def plot_errrcnts(gdat, gdatmodi, indxenerplot, indxevttplot):
   
    listtupl = [ \
                #['errrcnts', False, gdatmodi.thiserrrcnts, gdat.maxmerrrcnts[indxenerplot], gdat.tickerrrcnts[indxenerplot, :], gdat.lablerrrcnts[indxenerplot, :]], \
                #[    'errr',  True,     gdatmodi.thiserrr,     gdat.maxmerrr[indxenerplot],     gdat.tickerrr[indxenerplot, :],     gdat.lablerrr[indxenerplot, :]], \
                ['errrcnts', False, gdatmodi.thiserrrcnts, None, None, None], \
                [    'errr',  True,     gdatmodi.thiserrr, None, None, None], \
               ]
    numbiter = len(listtupl)
    for k in range(numbiter):
        strg, boolmean, varb, vmax, tick, labl = listtupl[k] 
        figr, axis, path = init_figr(gdat, indxenerplot, indxevttplot, strg, gdatmodi=gdatmodi)
        imag = retr_imag(gdat, axis, varb, indxenerplot, indxevttplot, vmax=vmax, cmap='RdBu', mean=boolmean)
        make_cbar(gdat, axis, imag, indxenerplot, tick=tick, labl=labl)
        plt.tight_layout()
        plt.savefig(path)
        plt.close(figr)
    
    
def plot_catlfram(gdat, gdatmodi, indxpoplplot, indxenerplot, indxevttplot):
    
    figr, axis, path = init_figr(gdat, indxenerplot, indxevttplot, 'catlfram', indxpoplplot=indxpoplplot, gdatmodi=gdatmodi)
    make_framlabl(gdat, axis)

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
    tempsamp[gdat.indxsampmeanpnts] = cdfn_logt(array([tempnumbpnts]), minmmeanpnts, factmeanpnts)
    tempsamp[gdat.indxsampfluxdistslop] = cdfn_atan(array([1.5]), minmfluxdistslop, factfluxdistslop)
    tempsamp[gdat.indxsamppsfp] = 0.5
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
        
    
    totlcntscart0 = tdpy.util.retr_cart(totlcntsheal0[0, :], minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                            minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
    totlcntscart1 = tdpy.util.retr_cart(totlcntsheal1[0, :], minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                            minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
    
    fig = plt.figure(figsize=(2 * gdat.plotsize, 2 * gdat.plotsize))
    axis = figr.add_subplot(221)
    imag = axis.imshow(totlcntscart0, origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    
    axis.set_xlabel(gdat.strgxaxi)
    axis.set_ylabel(gdat.strgyaxi)
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
    
    axis.scatter(tempsampvarb[trueindxsamplgal], tempsampvarb[trueindxsampbgal], s=50, alpha=gdat.alphmrkr, marker='x', color='g', linewidth=2)
    
    axis = figr.add_subplot(212)

    tdpy.mcmc.plot_braz(ax, meancnts, hist0,  lcol='lightgreen', alpha=gdat.alphmrkr, dcol='green', mcol='darkgreen', labl='Isotropic')
    tdpy.mcmc.plot_braz(ax, meancnts, hist1, lcol='lightblue', alpha=gdat.alphmrkr, dcol='blue', mcol='darkblue', labl='Isotropic + Unresolved PS')

    axis.set_xlabel('$k$')
    axis.set_ylabel('d$P$/d$k$')
    axis.set_yscale('log')
    axis.legend()

    plt.tight_layout()
    plt.savefig(gdat.pathplot + 'pntsdiff.pdf')
    plt.close(figr)
    
