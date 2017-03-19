# common imports
from __init__ import *

# internal functions
from util import *

def plot_samp(gdat, gdatmodi, strg):
    
    if strg == 'true':
        strgtype = strg
    else:
        strgtype = ''
    
    if gdatmodi != None:
        gdatobjt = gdatmodi
    else:
        gdatobjt = gdat
    
    spatdisttype = getattr(gdat, strgtype + 'spatdisttype')
    
    if strg != 'post':
    
        spectype = getattr(gdat, strgtype + 'spectype')
        indxpopl = getattr(gdat, strgtype + 'indxpopl')
        sampvarb = getattr(gdatobjt, strg + 'sampvarb')
        numbpnts = sampvarb[getattr(gdat, strgtype + 'indxfixpnumbpnts')].astype(int)
        if gdat.elemtype == 'lght' and gdat.numbener > 1:
            specplot = getattr(gdatobjt, strg + 'specplot') 
        indxpntsfull = list(getattr(gdatobjt, strg + 'indxpntsfull'))
        indxsampcomp = retr_indxsampcomp(gdat, indxpntsfull, spectype)
        
    # plots
    ## frame-only
    if gdatmodi != None:
        ## brightest PS
        # temp
        if False and gdat.elemtype == 'lght':
            if gdatmodi == None or sum(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]) != 0:
                plot_brgt(gdat, gdatmodi, strg)
    
    if gdatmodi != None:
        strgswep = '_%09d' % gdatmodi.cntrswep
    else:
        strgswep = ''
        
    # element features projected onto the data axes
    if strg == 'this' or strg  == 'true' and gdat.datatype == 'mock':
        alph = 0.1
        if strg == 'this':
            pathtemp = gdat.pathfram
            colr = 'b'
        else:
            pathtemp = gdat.pathinit
            colr = 'g'

        # PS spectra
        if gdat.numbener > 1 and gdat.elemtype == 'lght':
            for l in gdat.indxpopl:
                listxdat = []
                listplottype = []
                listydat = []
                for k in range(numbpnts[l]):
                    spec = copy(specplot[l][:, k])
                    if gdat.enerdiff:
                        spec *= gdat.meanenerplot**2
                    listplottype.append('line')
                    listxdat.append(gdat.meanenerplot)
                    listydat.append(spec)
                
                path = pathtemp + 'specpop%d%s.pdf' % (l, strgswep)
                tdpy.util.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, colr=colr, alph=alph, plottype=listplottype, \
                      limtxdat=[gdat.minmener, gdat.maxmener], lablydat='$E^2dN/dE$ [%s]' % gdat.lablfeatunit['flux'], limtydat=[amin(gdat.minmspecplot), amax(gdat.maxmspecplot)])
        
        # deflection profiles
        if gdat.elemtype == 'lens':
            xdat = gdat.binsanglplot * gdat.anglfact
            lablxdat = gdat.lablfeattotl['gang']
            for l in indxpopl:
                listydat = []
                listvlinfrst = []
                listvlinseco = []
                for k in arange(numbpnts[l]):
                    listydat.append(retr_deflcutf(gdat.binsanglplot, sampvarb[indxsampcomp['flux'][l][k]], sampvarb[indxsampcomp['asca'][l][k]], \
                                                                                                                    sampvarb[indxsampcomp['acut'][l][k]]) * gdat.anglfact)
                    listvlinfrst.append(sampvarb[indxsampcomp['asca'][l][k]] * gdat.anglfact) 
                    listvlinseco.append(sampvarb[indxsampcomp['acut'][l][k]] * gdat.anglfact)
                    
                    if sampvarb[indxsampcomp['acut'][l][k]] * gdat.anglfact < 0.:
                        print 'sampvarb[indxsampcomp[acut][l][k]]'
                        print sampvarb[indxsampcomp['acut'][l][k]] * gdat.anglfact
                        raise Exception('')

                listydat.append(xdat * 0. + gdat.anglfact * sampvarb[getattr(gdat, strgtype + 'indxfixpbeinhost')])
                path = pathtemp + 'deflsubhpop%d%s.pdf' % (l, strgswep)
                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, \
                                                        colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
    
    ## PSF radial profile
    if gdat.elemtype == 'lght':
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.oaxitype:
                    indxydat = [i, slice(None), m, 0]
                else:
                    indxydat = [i, slice(None), m]
                strgindxydat = '%d%d' % (i, m)
                lablxaxi = gdat.lablfeattotl['gang']
                plot_gene(gdat, gdatmodi, strg, 'psfn', 'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalyaxi='logt', \
                                                                            factxdat=gdat.anglfact, lablxaxi=lablxaxi, lablyaxi=r'$\mathcal{P}$')
                if gdat.oaxitype or gdat.trueoaxitype:
                    plot_gene(gdat, gdatmodi, strg, 'factoaxi', 'binsoaxi', indxydat=[i, m, slice(None)], strgindxydat=strgindxydat, \
                                                                                    factxdat=gdat.anglfact, lablxaxi=lablxaxi, lablyaxi=r'$f(\phi)$')

    # number of background counts per PSF
    if False:
        for i in gdat.indxener:
            path = gdat.pathplot + 'cntsbackfwhmflux%d_%09d.pdf' % (i, gdatmodi.cntrswep)
            tdpy.util.plot_maps(path, sum(gdatmodi.thiscntsbackfwhm, 2)[i, :], pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                                                minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                                minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)

    # element feature correlations
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            if strg != 'true' and gdat.trueinfo:
                plot_scatspec(gdat, l, gdatmodi)
                plot_scatspec(gdat, l, gdatmodi, plotdiff=True)
            for a, strgfrst in enumerate(gdat.liststrgfeatodim[l]):
                for b, strgseco in enumerate(gdat.liststrgfeatodim[l]):
                    if a < b:
                        for strgplottype in ['hist', 'scat']:
                            plot_elemtdim(gdat, gdatmodi, strg, l, strgplottype, strgfrst, strgseco)

    # temp -- restrict compfrac and other plots to indxmodlpntscomp
    plot_compfrac(gdat, gdatmodi, strg)
    
    for i in gdat.indxener:
        for m in gdat.indxevttplot:
            plot_scatpixl(gdat, gdatmodi, strg)
   
    # temp -- RGB data counts
    if gdat.numbener > 1:
        pass
   
    # element features
    if not (strg == 'true' and gdat.datatype == 'inpt'):
        for l in gdat.indxpopl:
            indxydat = [l, slice(None)]
            strgindxydat = '%d' % l

            for strgfeat in gdat.liststrgfeatodim[l]:
                # temp
                if not (strgfeat == 'gang' and gdat.spatdisttype[l] != 'gang' or strgfeat == 'curv' and gdat.spectype[l] != 'curv' or \
                                                                                                                  strgfeat == 'expo' and gdat.spectype[l] != 'expo'):
                
                    if strgfeat == 'flux':
                        if gdat.elemtype == 'lens':
                            strgxaxitwin = 'mass'
                            lablxaxitwin = r'$M$ [$M_{\odot}$]'
                        if gdat.elemtype == 'lght':
                            strgxaxitwin = 'cntsplot'
                            lablxaxitwin = '$C$'
                        offs = [0.9, 0.8]
                    else:
                        strgxaxitwin = None
                        lablxaxitwin = None
                        offs = None
                    scalxaxi = gdat.dictglob['scal' + strgfeat + 'plot']
                    factxaxi = gdat.dictglob['fact' + strgfeat + 'plot']
                    lablxaxi = gdat.lablfeattotl[strgfeat]
                    limtxdat = [getattr(gdat, 'minm' + strgfeat + 'plot') * factxaxi, getattr(gdat, 'maxm' + strgfeat + 'plot') * factxaxi]
                    
                    limtydat = gdat.limtpntshist
                    plot_gene(gdat, gdatmodi, strg, strgfeat + 'hist', 'mean' + strgfeat, scalyaxi='logt', lablxaxi=lablxaxi, lablyaxi=r'$N$', factxdat=factxaxi, \
                                          scalxaxi=scalxaxi, limtxdat=limtxdat, limtydat=limtydat, offslegd=offs, indxydat=indxydat, strgindxydat=strgindxydat, histodim=True, \
                                          strgxaxitwin=strgxaxitwin, lablxaxitwin=lablxaxitwin)
           
    if gdatmodi != None:
        gdatobjt = gdatmodi
    else:
        gdatobjt = gdat

    if gdat.elemtype == 'lght' and gdat.calcerrr:
        if strg != 'true':
            for i in gdat.indxener:
                plot_genemaps(gdat, gdatmodi, strg, 'errrcnts', strgcbar='resicnts', thisindxener=i)
    
    if strg != 'true':
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                plot_genemaps(gdat, gdatmodi, strg, 'datacnts', thisindxener=i, thisindxevtt=m)
                if gdat.numbpopl > 1:
                    for l in gdat.indxpopl:
                        plot_genemaps(gdat, gdatmodi, strg, 'datacnts', thisindxener=i, thisindxevtt=m, thisindxpopl=l)
    
    # temp
    if 'gaus' in spatdisttype:
        plot_genemaps(gdat, gdatmodi, strg, 'lpdfspatpriointp', tdim=True)

    if strg == 'post':
        listbool = [False, True]
    else:
        listbool = [False]
    for stdv in listbool:
        
        if gdat.elemtype == 'lens':
            plot_genemaps(gdat, gdatmodi, strg, 'conv', tdim=True, stdv=stdv)
            plot_genemaps(gdat, gdatmodi, strg, 'convelem', strgcbar='conv', tdim=True, stdv=stdv)
            for i in gdat.indxener:
                plot_genemaps(gdat, gdatmodi, strg, 'hostcntsmaps', strgcbar='datacnts', thisindxener=i, tdim=True, stdv=stdv)
                plot_genemaps(gdat, gdatmodi, strg, 'lenscnts', strgcbar='datacnts', thisindxener=i, tdim=True, stdv=stdv)
                plot_genemaps(gdat, gdatmodi, strg, 'modlfluxuncv', strgcbar='datacnts', thisindxener=i, tdim=True, stdv=stdv)
            if strg != 'true':
                plot_genemaps(gdat, gdatmodi, strg, 'deflcomp', tdim=True, stdv=stdv)
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if strg != 'true':
                    plot_genemaps(gdat, gdatmodi, strg, 'llik', strgcbar='llik', thisindxener=i, thisindxevtt=m, stdv=stdv)
                plot_genemaps(gdat, gdatmodi, strg, 'modlcnts', strgcbar='datacnts', thisindxener=i, thisindxevtt=m, stdv=stdv)
                plot_genemaps(gdat, gdatmodi, strg, 'resicnts', thisindxener=i, thisindxevtt=m, stdv=stdv)
        
    if gdat.elemtype == 'lens':
        plot_gene(gdat, gdatmodi, strg, 'convpsecelemodim', 'meanmpolodim', lablxaxi='$k$ [1/kpc]', lablyaxi='$P_{subs}(k)$', ylim=[1., 1e3], scalxaxi='logt', scalyaxi='logt') #\
        plot_gene(gdat, gdatmodi, strg, 'convpsecodim', 'meanmpolodim', lablxaxi='$k$ [1/kpc]', lablyaxi='$P(k)$', ylim=[1., 1e3], scalxaxi='logt', scalyaxi='logt') #\
#                                                                                                    strgxaxitwin='meananglodim', lablxaxitwin=gdat.lablgang, scalyaxi='logt')
        plot_gene(gdat, gdatmodi, strg, 'convpsecodim', 'meanwvecodim', lablxaxi='$k$ [1/kpc]', lablyaxi='$P(k)$', ylim=[0.1, 1e4], scalxaxi='logt', scalyaxi='logt') #\
#                                                                                                   strgxaxitwin='meanwlenodim', lablxaxitwin='$l$', scalyaxi='logt')
        plot_gene(gdat, gdatmodi, strg, 'histdefl', 'meandefl', scal='self', lablxaxi=r'$\alpha$ [arcsec]', lablyaxi=r'$N_{pix}$', factxdat=gdat.anglfact, histodim=True)

    for l in gdat.indxpopl:
        if gdat.elemtype == 'lens':

            # overall deflection field
            plot_defl(gdat, gdatmodi, strg)
            
            # deflection field due to individual lenses
            for k in range(gdat.numbdeflsingplot):
                plot_defl(gdat, gdatmodi, strg, indxdefl=k)
            
            # residual deflection field
            if strg != 'true':
                plot_defl(gdat, gdatmodi, strg, strgcomp='resi')
                for k in range(gdat.numbdeflsingplot):
                    plot_defl(gdat, gdatmodi, strg, strgcomp='resi', indxdefl=k)

        if strg != 'true':
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    #if strg == 'this':
                    #    plot_datacnts(gdat, gdatmodi, i, m, l)
                    pass
                    

def plot_post(gdat=None, pathpcat=None, verbtype=1, makeanim=False, writ=True, prio=False):
    
    if writ:
        gdat = tdpy.util.gdatstrt()
        gdat.verbtype = verbtype
        gdat.makeanim = makeanim

        if gdat.verbtype > 0:
            print 'Reading %s...' % pathpcat

        # read PCAT output file
        shel = shelve.open(pathpcat)
        for keyy in shel:
            if keyy == 'gdat':
                gdat = shel[keyy]
        shel.close()
    
    if gdat.verbtype > 0:
        print 'Producing postprocessing plots...'

    timetotlinit = gdat.functime()
    
    gdat.strgbest = 'ML'
  
    # prior components
    gdat.indxlpri = arange(gdat.numblpri)
    gdat.indxlpau = arange(gdat.numblpau)
    indxbadd = where(isfinite(gdat.listlpri) == False)
    if indxbadd[0].size > 0:
        gdat.listlpri[indxbadd] = 5.
    indxbadd = where(isfinite(gdat.listlpriprop) == False)
    if indxbadd[0].size > 0:
        gdat.listlpriprop[indxbadd] = 5.
    
    for n in gdat.indxproptype:
        for k in gdat.indxlpri:
            if (gdat.listlpri[:, k] != 0.).any():
                path = gdat.pathpostlpri + 'lpri%04d' % k + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlpri[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
            if (gdat.listlpriprop[:, k] != 0.).any():
                path = gdat.pathpostlpri + 'lpriprop%04d' % k + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlpriprop[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
            if (gdat.listlpriprop[:, k] - gdat.listlpri[:, k] != 0.).any():
                path = gdat.pathpostlpri + 'lpridelt%04d' % k + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlpriprop[gdat.listindxsamptotl[n], k] - gdat.listlpri[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
        for k in gdat.indxlpau:
            if (gdat.listlpau[:, k] != 0.).any():
                path = gdat.pathpostlpri + 'lpau%04d' % k + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlpau[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
        if (gdat.listlfctprop[:, 0] != 0.).any():
            path = gdat.pathpostlpri + 'lfctprop' + gdat.nameproptype[n]
            tdpy.mcmc.plot_trac(path, gdat.listlfctprop[gdat.listindxsamptotl[n], 0], '$q_{lfct}$', logthist=True)
    
    # Gelman-Rubin test
    # temp
    if False and gdat.numbproc > 1:
        if isfinite(gdat.gmrbstat).all():
            if gdat.verbtype > 0:
                print 'Gelman-Rubin TS...'
    
            figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
            minm = min(amin(gdat.gmrbstat), amin(gdat.gmrbfixp))
            maxm = max(amax(gdat.gmrbstat), amax(gdat.gmrbfixp))
            bins = linspace(minm, maxm, 40)
            axis.hist(gdat.gmrbstat.flatten(), bins=bins, label='Data proj.')
            axis.hist(gdat.gmrbfixp, bins=bins, label='Fixed dim.')
            axis.set_xlabel('PSRF')
            axis.set_ylabel('$N_{stat}$')
            plt.tight_layout()
            figr.savefig(gdat.pathplot + 'diag/gmrbhist.pdf')
            plt.close(figr)
            
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    maps = gdat.gmrbstat[i, :, m]
                    path = gdat.pathdiag + 'gmrbmaps_%d%d.pdf' % (i, m)
                    tdpy.util.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                            minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                            minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
        else:
            print 'Inappropriate Gelman-Rubin test statistics encountered.'

    # plot autocorrelation
    if gdat.verbtype > 0:
        print 'Autocorrelation...'
    tdpy.mcmc.plot_atcr(gdat.pathdiag, gdat.atcr, gdat.timeatcr)
    
    # plot proposal efficiency
    if gdat.verbtype > 0:
        print 'Acceptance ratio...'
    numbtimemcmc = 20
    binstimemcmc = linspace(0., gdat.numbswep, numbtimemcmc)
    numbtick = 2
    sizefigryaxi = max(gdat.numbproptype * gdat.plotsize / 4., gdat.plotsize / 2.)
    figr, axgr = plt.subplots(gdat.numbproptype, 1, figsize=(gdat.plotsize, sizefigryaxi), sharex='all')
    if gdat.numbproptype == 1:
        axgr = [axgr]
    for n, axis in enumerate(axgr):
        histtotl = axis.hist(gdat.listindxsamptotl[n], bins=binstimemcmc)[0]
        histaccp = axis.hist(gdat.listindxsamptotlaccp[n], bins=binstimemcmc)[0]
        axis.set_ylabel('%s' % gdat.lablproptype[n])
        if n == gdat.numbproptype - 1:
            axis.set_xlabel('$i_{samp}$')
        # define the y-axis
        maxm = amax(histtotl)
        axis.set_ylim([0., maxm])
        try:
            axis.set_title('%.3g%%' % int(sum(histaccp) / sum(histtotl) * 100.))
        except:
            pass
        listtick = linspace(maxm / 2., maxm, numbtick)
        listlabltick = ['%.3g' % tick for tick in listtick]
        axis.set_yticks(listtick)
        axis.set_yticklabels(listlabltick)
    plt.tight_layout()
    figr.savefig(gdat.pathdiag + 'accpratiproptype.pdf')
    plt.close(figr)
   
    # post-processing frame plots
    if gdat.numbframpost != None:
        gdatmodi = tdpy.util.gdatstrt()
        gdat.indxsamptotlfram = arange(gdat.numbframpost) * (gdat.indxsamptotl - gdat.indxsamptotl % gdat.numbframpost) / gdat.numbframpost
        for n in gdat.indxsamptotlfram:
            gdatmodi.cntrswep = n
            gdatmodi.thisindxpntsfull = deepcopy(gdat.listindxpntsfull[n])
            gdatmodi.thissampvarb = copy(gdat.listsampvarb[n, :])
            proc_samp(gdat, gdatmodi, 'this')
            plot_samp(gdat, gdatmodi, 'this')

    # plot split and merge diagnostics
    if gdat.numbtrap > 0 and gdat.probtran > 0. and gdat.probbrde < 1.:
        if gdat.verbtype > 0:
            print 'Split and merge related plots...'
    
        indxsampsplttotl = where(gdat.listindxprop == gdat.indxpropsplt)[0]
        indxsampsplt = intersect1d(where(gdat.listindxprop == gdat.indxpropsplt)[0], where(gdat.listaccpprop)[0])
        indxsampmergtotl = where(gdat.listindxprop == gdat.indxpropmerg)[0]
        indxsampmerg = intersect1d(where(gdat.listindxprop == gdat.indxpropmerg)[0], where(gdat.listaccpprop)[0])
        indxsampspmr = union1d(indxsampsplt, indxsampmerg)
        indxsampspmrtotl = union1d(indxsampsplttotl, indxsampmergtotl)
        if indxsampspmrtotl.size > 0:
    
            ## labels and names
            listlabl = ['$u_f$', '$u_r$', r'$u_\phi$', '$u_s$', '$N_{pair}$', \
                                            r'$\alpha_c\alpha_j$', r'$\alpha_P\alpha_c\alpha_j$', r'$\alpha_c$', r'$\alpha_j$', r'$\alpha_L$', r'$\alpha_P$']
            listname = ['fracauxi', 'radiauxi', 'anglauxi', 'sindauxi', 'numbpair', 'laccfact', 'laccfacttotl', 'combfact', 'jcbnfct']
            
            ## variables
            listvarb = [gdat.listauxipara[:, 0], gdat.anglfact * gdat.listauxipara[:, 1], gdat.listauxipara[:, 2], gdat.listauxipara[:, 3], gdat.listnumbpair, \
                                             exp(gdat.listlaccfact), exp(gdat.listlaccfact + gdat.listdeltlpri), gdat.listcombfact, \
                                             gdat.listjcbnfact]
           
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
    
                make_legd(axis)
                plt.tight_layout()
                figr.savefig(gdat.pathpostspmr + listname[k] + '.pdf')
                plt.close(figr)
    
    if gdat.verbtype > 0:
        print 'Proposal execution times...'
    plot_chro(gdat)

    if gdat.verbtype > 0:
        print 'Plotting the posterior...'
    # posterior versions of the frame plots
    plot_samp(gdat, None, 'post')
    
    if gdat.verbtype > 0:
        print 'Plotting a mosaic of samples...'
    
    ## mosaic of images of posterior catalogs 
    plot_mosa(gdat)
    
    if gdat.verbtype > 0:
        print 'Fixed dimensional parameter traces...'
    
    if gdat.elemtype == 'lens':
        path = gdat.pathpost + 'fracsubh'
        tdpy.mcmc.plot_trac(path, gdat.listfracsubh, gdat.lablfracsubh, truepara=gdat.truefracsubh)
    
    ## fixed-dimensional parameters
    ### trace and marginal distribution of each parameter
    for k in gdat.indxfixp:
        path = gdat.pathpostfixp + gdat.namefixp[k]
        tdpy.mcmc.plot_trac(path, gdat.listfixp[:, k], gdat.lablfixp[k], truepara=gdat.corrfixp[k], scalpara=gdat.scalfixp[k])
    
    if gdat.checprio and not prio:
        # this works only for scalar variables -- needs to be generalized to all variables
        for namevarbscal in gdat.listnamevarbscal:
            titl = '$D_{KL} = %.3g$' % getattr(gdat, 'infototl' + namevarbscal)
            xdat = getattr(gdat, 'mean' + namevarbscal)
            lablxdat = getattr(gdat, 'labl' + namevarbscal)
            
            path = gdat.pathpost + 'info' + namevarbscal
            ydat = getattr(gdat, 'info' + namevarbscal)
            lablydat = r'$\Delta D_{KL}$'
            tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, titl=titl)

            path = gdat.pathpost + 'pdfn' + namevarbscal
            ydat = [getattr(gdat, 'pdfnpost' + namevarbscal), getattr(gdat, 'pdfnprio' + namevarbscal)]
            lablydat = r'$P$'
            legd = ['$P$(%s|$D$)' % lablxdat, '$P$(%s)' % lablxdat]
            tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, titl=titl, colr=['k', 'k'], linestyl=['-', '--'], legd=legd)
    
    if gdat.verbtype > 0:
        print 'Fixed dimensional parameter covariance...'
    
    ### covariance
    #### overall
    path = gdat.pathpostfixp + 'fixp'
    tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.indxfixpprop] * gdat.factfixpplot[None, gdat.indxfixpprop], gdat.lablfixp[gdat.indxfixpprop], \
                                                                                          truepara=gdat.corrfixp[gdat.indxfixpprop] * gdat.factfixpplot[gdat.indxfixpprop])
    
    #### individual processes
    if gdat.numbproc > 1:
        for k in gdat.indxproc:
            path = gdat.pathpostfixpproc + 'proc%d' % k
            tdpy.mcmc.plot_grid(path, gdat.listsampvarbproc[:, k, gdat.indxfixpprop] * gdat.factfixpplot[None, gdat.indxfixpprop], gdat.lablfixp[gdat.indxfixpprop], \
                                                                                      truepara=gdat.corrfixp[gdat.indxfixpprop] * gdat.factfixpplot[gdat.indxfixpprop])
    
    ### grouped covariance plots
    if gdat.verbtype > 0:
        print 'Hyperparameters...'
    
    #### hyperparameters
    path = gdat.pathpostfixp + 'hypr'
    tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.indxfixphypr], gdat.lablfixp[gdat.indxfixphypr], truepara=[gdat.corrfixp[k] for k in gdat.indxfixphypr])
    
    if gdat.verbtype > 0:
        print 'PSF parameters...'
    
    #### PSF
    if gdat.proppsfp:
        path = gdat.pathpostfixp + 'psfp'
        tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.indxfixppsfp], gdat.lablfixp[gdat.indxfixppsfp], truepara=[gdat.corrfixp[k] for k in gdat.indxfixppsfp], \
                                                                                                                                        numbplotside=gdat.numbpsfptotl)
    if gdat.verbtype > 0:
        print 'Background parameters...'
    
    #### backgrounds
    if gdat.propbacp:
        path = gdat.pathpostfixp + 'bacp'
        tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.indxfixpbacp], gdat.lablfixp[gdat.indxfixpbacp], \
                                                                                                        truepara=[gdat.corrfixp[k] for k in gdat.indxfixpbacp])
        if gdat.numbback == 2 and gdat.specback == [None, None]:
            for i in gdat.indxener:
                indx = gdat.indxfixpbacp[gdat.indxback*gdat.numbener+i]
                path = gdat.pathpostfixp + 'bacpene%d' % i
                tdpy.mcmc.plot_grid(path, gdat.listfixp[:, indx], gdat.lablfixp[indx], truepara=[gdat.corrfixp[k] for k in indx], join=True)
    
    if gdat.verbtype > 0:
        print 'Transdimensional parameters...'
    
    ## randomly selected trandimensional parameters
    if gdat.numbtrap > 0:
        # choose the parameters based on persistence
        stdvlistsamptran = std(gdat.listsamp[:, gdat.indxsamptrap], axis=0)
        indxtrapgood = where(stdvlistsamptran > 0.)[0]
        numbtrapgood = indxtrapgood.size
        numbtrapplot = min(10, numbtrapgood)
        if numbtrapplot > 0:
            indxtrapplot = sort(choice(gdat.indxsamptrap[indxtrapgood], size=numbtrapplot, replace=False))
            path = gdat.pathpost + 'listsamp'
            tdpy.mcmc.plot_grid(path, gdat.listsamp[:, indxtrapplot], ['%d' % k for k in indxtrapplot])
            path = gdat.pathpost + 'listsampvarb'
            tdpy.mcmc.plot_grid(path, gdat.listsampvarb[:, indxtrapplot], ['%d' % k for k in indxtrapplot])

    if gdat.verbtype > 0:
        print 'Binned transdimensional parameters...'
   
    # stacked posteiors binned in position and flux
    if gdat.numbtrap > 0:
        liststrgbins = ['quad', 'full']
        for l in gdat.indxpopl:
            plot_postbindmaps(gdat, l, 'cumu')
            for strgbins in liststrgbins:
                for strgfeatsign in gdat.liststrgfeatsign[l]:
                    plot_postbindmaps(gdat, l, strgbins, strgfeatsign)

    if gdat.verbtype > 0:
        print 'Prior and likelihood...'
    
    ## log-prior and log-likelihood
    gdat.listlliktotl = sum(gdat.listllik, axis=(1, 2, 3))
    gdat.listlpritotl = sum(gdat.listlpri, axis=1)

    for strg in ['llik']:

        if strg == 'lpri':
            labltemp = '\ln P(x)'
        if strg == 'llik':
            labltemp = '\ln P(D|x)'
        
        setattr(gdat, 'list' + strg + 'flat', getattr(gdat, 'list' + strg + 'totl').flatten())
        setattr(gdat, 'listdelt' + strg + 'totlflat', getattr(gdat, 'listdelt' + strg + 'totl').flatten())
        
        labl = r'$%s$' % labltemp
        labldelt = r'$\Delta%s$' % labltemp

        path = gdat.pathpost + strg
        titl = r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi)
        tdpy.mcmc.plot_hist(path, getattr(gdat, 'list' + strg + 'flat'), labl, titl)
        if strg == 'llik':
            varbdraw = [gdat.maxmllikswep]
            labldraw = ['Maximum likelihood Sample']
        else:
            varbdraw = None
            labldraw = None
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'list' + strg + 'flat'), labl, varbdraw=varbdraw, labldraw=labldraw)

        path = getattr(gdat, 'pathpostdelt%s' % strg) + 'delt%s' % strg
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat'), labldelt)
        if gdat.numbproc > 1:
            for k in gdat.indxproc:
                path = getattr(gdat, 'pathpostdelt%s' % strg) + 'delt%s_proc%04d' % (strg, k)
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totl')[:, k], labldelt, titl='Process %d' % k)
        for n in gdat.indxproptype:
            path = getattr(gdat, 'pathpostdelt%s' % strg) + 'delt%s_%s' % (strg, gdat.nameproptype[n])
            tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotl[n]], labldelt, titl=gdat.nameproptype[n])
            path = getattr(gdat, 'pathpostdelt%saccp' % strg) + 'delt%s_%s_accp' % (strg, gdat.nameproptype[n])
            tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotlaccp[n]], labldelt, titl=gdat.nameproptype[n] + ', Accepted')
            path = getattr(gdat, 'pathpostdelt%sreje' % strg) + 'delt%s_%s_reje' % (strg, gdat.nameproptype[n])
            tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotlreje[n]], labldelt, titl=gdat.nameproptype[n] + ', Rejected')
        
    # plot resident memory
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxswep, mean(gdat.listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(gdat.pathdiag + 'memoresi.pdf')
    plt.close(figr)

    # animate the frame plots
    if gdat.makeanim:
        make_anim(gdat)

    timetotlfinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Plots and animations are produced in %.3g seconds.' % (timetotlfinl - timetotlinit)


def plot_chro(gdat):

    gdat.listchrototl *= 1e3
    indxchro = array([0, 1, 2, 4])
    binstime = logspace(log10(amin(gdat.listchrototl[where(gdat.listchrototl > 0)])), log10(amax(gdat.listchrototl[:, indxchro])), 50)

    labl = ['Choice', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Posterior', 'Advance', 'Total']
    listcolr = ['b', 'g', 'r', 'm', 'orange', 'cyan', 'yellow', 'black']
    numblabl = len(labl)
    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    for k in range(numblabl):
        varb = gdat.listchrototl[:, k]
        axis.hist(varb, binstime, log=True, label=labl[k], color=listcolr[k], alpha=0.3)
    
    for k in range(numblabl):
        varb = gdat.listchrototl[:, k]
        axis.hist(varb, binstime, log=True, edgecolor=listcolr[k], linewidth=5, facecolor='none')

    axis.set_title(r'$\langle t \rangle$ = %.3g ms' % mean(gdat.listchrototl[where(gdat.listchrototl[:, 0] > 0)[0], 0]))
    axis.set_xlim([amin(binstime), amax(binstime)])
    axis.set_xscale('log')
    axis.set_ylim([0.5, None])
    make_legd(axis)
    axis.set_xlabel('$t$ [ms]')
    
    plt.tight_layout()
    figr.savefig(gdat.pathdiag + 'chrototl.pdf')
    plt.close(figr)

    #for k in range(numblabl):
    #    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    #    if k == numblabl - 1:
    #        varb = gdat.listchrototl[:, 0] - sum(gdat.listchrototl[:, 1:], 1)
    #    else:
    #        varb = gdat.listchrototl[:, k]
    #    axis.hist(varb, binstime, log=True)
    #    axis.set_title(labl[k])
    #    axis.set_xlim([amin(binstime), amax(binstime)])
    #    axis.set_xscale('log')
    #    axis.set_xlabel('$t$ [ms]')
    #    axis.set_ylim([0.5, None])
    #    figr.savefig(gdat.pathdiag + 'chro%04d.pdf' % k)
    #    plt.close(figr)

    gdat.listchroproc *= 1e3
   
    if (gdat.listchroproc != 0).any():
        listlabl = ['PSF Intp.', 'Variables', 'Pixel mesh', 'Energy mesh', 'PS flux', 'Total flux', 'Lens host', 'Lens source', 'PSF conv.', 'Counts', 'Unbinned', 'Likelihood']
        numblabl = len(listlabl)
        figr, axcl = plt.subplots(gdat.numbchroproc, 1, figsize=(2 * gdat.plotsize, gdat.plotsize * numblabl / 3.))
        maxmchroproc = amax(gdat.listchroproc)
        minmchroproc = amin(gdat.listchroproc[where(gdat.listchroproc > 0)])
        binstime = logspace(log10(minmchroproc), log10(maxmchroproc), 50)
    
        for k in range(gdat.numbchroproc):
            chro = gdat.listchroproc[where(gdat.listchroproc[:, k] > 0)[0], k]
            try:
                axcl[k].hist(chro, binstime, log=True, label=listlabl[k])
            except:
                pass
            axcl[k].set_xlim([amin(binstime), amax(binstime)])
            axcl[k].set_ylim([0.5, None])
            axcl[k].set_ylabel(listlabl[k])
            axcl[k].set_xscale('log')
            if k != gdat.numbchroproc - 1:
                axcl[k].set_xticklabels([])
            axcl[k].axvline(mean(chro), ls='--', alpha=0.2, color='black')
        axcl[-1].set_xlabel('$t$ [ms]')
        plt.subplots_adjust(hspace=0.05)
        figr.savefig(gdat.pathdiag + 'chroproc.pdf')
        plt.close(figr)


def plot_compfrac(gdat, gdatmodi, strg):
    
    listydat = zeros((gdat.numblablcompfracspec, gdat.numbener))
    listyerr = zeros((2, gdat.numblablcompfracspec, gdat.numbener))
   
    cntr = 0
        
    if strg == 'true':
        strgtype = strg
    else:
        strgtype = ''
    
    specback = getattr(gdat, strgtype + 'specback')

    ## background templates
    for c in gdat.indxback:
        temp = retr_varb(gdat, gdatmodi, strg, 'sampvarb', indx=[gdat.indxfixpbacp[gdat.indxbacpback[c]]])
        if specback[c] != None:
            norm = temp * specback[c]
        else:
            norm = temp
        listydat[cntr, :] = norm * gdat.backfluxmean[c, :]

        if strg == 'post':
            temp = retr_varb(gdat, gdatmodi, strg, 'sampvarb', indx=[gdat.indxfixpbacp[gdat.indxbacpback[c]]], perc='errr')
            if specback[c] != None:
                norm = temp * specback[c]
            else:
                norm = temp
            listyerr[:, cntr, :] = norm * gdat.backfluxmean[None, c, :]

        if gdat.elemtype == 'lens':
            listydat[cntr, :] *= 4. * gdat.maxmgang**2
        
        cntr += 1
    
    ## PS
    if gdat.elemtype == 'lght':
        listydat[cntr, :] = retr_varb(gdat, gdatmodi, strg, 'pntsfluxmean')
        if strg == 'post':
            listyerr[:, cntr, :] = retr_varb(gdat, gdatmodi, strg, 'pntsfluxmean', perc='errr')
        cntr += 1

    if gdat.elemtype == 'lens':
        for strgtemp in ['sour', 'host']:
            indxvarb = getattr(gdat, 'indxfixpspec' + strgtemp)
            if gdatmodi == None:
                if strg == 'post':
                    listydat[cntr, :] = getattr(gdat, 'medifixp')[indxvarb]
                    # temp -- indxfixpspec**** should be multidimensional
                    listyerr[:, cntr, :] = getattr(gdat, 'errrfixp')[:, indxvarb]
                else:
                    listydat[cntr, :] = gdat.truesampvarb[indxvarb]
            else:
                listydat[cntr, :] = gdatmodi.thissampvarb[indxvarb]
            cntr += 1
            
    cntrdata = cntr

    ## data
    listydat[cntr, :] = gdat.datafluxmean
    cntr += 1
    
    ## total model
    if gdat.numblablcompfrac > 1:
        if specback[c] != None:
            norm = temp * specback[c]
        else:
            norm = temp
        listydat[cntr, :] = sum(listydat[:cntrdata, :], 0)
        if gdatmodi == None:
            listyerr[:, cntr, :] = mean(listydat[:cntrdata, :], 0)
        cntr += 1

    # plot energy spectra of the data, background model components and total background
    if gdat.numbener > 1:
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        
        listmrkr = [(2 + k/2, 1 + k % 2, 0) for k in range(16)]

        # plot reference spectra
        if gdat.listspecrefrplot != None:
            for k in range(len(gdat.listspecrefrplot)):
                axis.plot(gdat.listenerrefrplot[k], gdat.listspecrefrplot[k], label=gdat.listlablrefrplot[k], color='m')

        xdat = gdat.meanener
        for k in range(gdat.numblablcompfracspec):
            if k == cntrdata:
                colr = 'black'
                linestyl = '-'
                mrkr = None
            else:
                if strg == 'true':
                    colr = 'g'
                    linestyl = '-.'
                    mrkr = listmrkr[k]
                elif strg == 'this':
                    colr = 'b'
                    linestyl = '-.'
                    mrkr = listmrkr[k]
                elif strg == 'post':
                    colr = 'black'
                    linestyl = ''
                    mrkr = 'o'
            ydat = listydat[k, :]
            yerr = listyerr[:, k, :]
            if gdat.enerdiff:
                ydat *= gdat.meanener**2
                yerr *= gdat.meanener**2
            axis.errorbar(xdat, ydat, yerr=yerr, label=gdat.listlablcompfracspec[k], color=colr, marker=mrkr, ls=linestyl, markersize=15)
    
        axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])
        try:
            gdat.ylimcompfracspec
        except:
            gdat.ylimcompfracspec = axis.get_ylim()
        axis.set_ylim([1e-4 * amax(gdat.datafluxmean), 1e1 * amax(gdat.datafluxmean)])
        axis.set_yscale('log')
        axis.set_xlabel('$E$ [%s]' % gdat.strgenerunit)
        axis.set_xscale('log')
        axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [%s/cm$^2$/s/sr]' % gdat.strgenerunit)
        make_legd(axis, numbcols=3)
        
        if gdat.exprtype == 'chan':
            axis.set_ylim([None, 1e4])

        plt.tight_layout()
        path = retr_plotpath(gdat, gdatmodi, strg, 'compfracspec')
        figr.savefig(path)
        plt.close(figr)
    

    # pie plot illustrating contribution of each background template (and PS) to the total model
    listexpl = zeros(gdat.numblablcompfrac)
    listsize = zeros(gdat.numblablcompfrac)
    
    for k in range(gdat.numblablcompfrac):
        if gdat.enerdiff:
            listsize[k] = sum(listydat[k, :] * gdat.deltener)
        else:
            listsize[k] = listydat[k, :]
    

    listsize *= 100. / sum(listsize)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))

    axis.pie(listsize, explode=listexpl, labels=gdat.listlablcompfrac, autopct='%1.1f%%')
    axis.axis('equal')

    plt.subplots_adjust(top=0.8, bottom=0.2, left=0.2, right=0.8)
    path = retr_plotpath(gdat, gdatmodi, strg, 'compfrac')
    figr.savefig(path)
    plt.close(figr)
     

def plot_pdfntotlflux():

    minm = 1e-9
    maxm = 10e-9
    numbvarb = 90
    numbsamp = 100000
    numbbins = 40
    alph = 0.5
    
    binssing = linspace(minm, maxm, numbvarb + 1)
    meansing = (binssing[:-1] + binssing[1:]) / 2.
    deltsing = binssing[1:] - binssing[:-1]
    
    binsdoub = linspace(2. * minm, 2. * maxm, 2 * numbvarb)
    meandoub = (binsdoub[:-1] + binsdoub[1:]) / 2.
    deltdoub = binsdoub[1:] - binsdoub[:-1]
    
    bins = linspace(minm, 2. * maxm, 2 * numbvarb + 1)
    
    arry = empty((2, numbsamp))
    
    minmslop = 1.5
    maxmslop = 3.
    numbslop = 4
    sloparry = linspace(minmslop, maxmslop, numbslop)
    for n in range(numbslop):
        slop = sloparry[n]
        for k in range(2):
            arry[k, :] = (rand(numbsamp) * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
        
        totl = sum(arry, 0)
        
        powrprob = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop)) * meansing**(-slop)
        
        convprob = convolve(powrprob, powrprob) * deltdoub[0]
        
        indxdoub = where(meandoub <= maxm)[0]
        convprobpoly = polyval(polyfit(meandoub[indxdoub], convprob[indxdoub], 8), meandoub[indxdoub])
        
        figr, axis = plt.subplots()
        axis.hist(arry[k, :], bins=bins, alpha=alph, label='$f_1$ (Sampled)', color='b')
        axis.hist(totl, bins=bins, alpha=alph, label='$f_0$ (Sampled)', color='g')
        axis.plot(meansing, powrprob * numbsamp * deltsing, label='$f_1$ (Analytic)', color='b')
        axis.plot(meandoub, convprob * numbsamp * deltdoub[0], label='$f_0$ (Numerically convolved)', color='g')
        
        axis.plot(meandoub[indxdoub], convprobpoly * numbsamp * deltdoub[indxdoub], label='$f_0$ (Fit)', color='r')
    
        axis.set_ylim([0.5, numbsamp])
        axis.set_xlabel('$f$')
        axis.set_xlim([amin(bins), amax(bins)])
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_ylabel('$N_{samp}$')
        make_legd(axis)
        plt.tight_layout()
        pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/powrpdfn/'
        os.system('mkdir -p ' + pathfold)
        figr.savefig(pathfold + 'powrpdfn%04d.pdf' % n)
        plt.close(figr)
        

def plot_brgt(gdat, gdatmodi, strg):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if gdatmodi == None:
        # temp
        pass
    else:   
        if fluxbrgt.size > 0:
            axis.scatter(fluxbrgt, fluxbrgtassc, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
            axis.scatter(fluxbrgt[0], sum(fluxbrgtassc), alpha=gdat.alphmrkr, color='b', label='Sample - Total')
    if gdat.truefluxbrgt.size > 0:
        axis.scatter(gdat.truefluxbrgt, gdat.truefluxbrgtassc, alpha=gdat.alphmrkr, color='g', label=gdat.truelabl)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_ylim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_xlabel(r'$%s_{max}$%s' % (gdat.lablfeat['flux'], gdat.lablfeatunit['flux']))
    axis.set_ylabel(r'$%s_{asc}$%s' % (gdat.lablfeat['flux'], gdat.lablfeatunit['flux']))
    make_legd(axis, loca=2)
    
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatbrgt')
    figr.savefig(path)
    plt.close(figr)
    

def plot_elemtdim(gdat, gdatmodi, strg, l, strgplottype, strgfrst, strgseco, strgmome='medi'):
    
    sizelarg = 10
    sizesmll = 1
    
    legdmome = getattr(gdat, 'legd' + strgmome)

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if strg == 'post':
        varb = getattr(gdat, strgmome + strgfrst + strgseco + 'hist')[l, :, :]
        varbfrst = getattr(gdat, 'bins' + strgfrst + 'plot') * gdat.dictglob['fact' + strgfrst + 'plot']
        varbseco = getattr(gdat, 'bins' + strgseco + 'plot') * gdat.dictglob['fact' + strgseco + 'plot']
        labl = gdat.legdsampdist + ' ' + legdmome
        imag = axis.pcolor(varbfrst, varbseco, varb.T, cmap='Greys', label=labl)
        #imag = axis.imshow(gdat.medifluxsindhist[l], cmap='Purples', extent=[])
        plt.colorbar(imag) 

    if strg == 'this':
        if strgplottype == 'hist':
            meanfrst = getattr(gdat, 'bins' + strgfrst + 'plot') * gdat.dictglob['fact' + strgfrst + 'plot']
            meanseco = getattr(gdat, 'bins' + strgseco + 'plot') * gdat.dictglob['fact' + strgseco + 'plot']
            hist = getattr(gdatmodi, strg + strgfrst + strgseco + 'hist')[l, :, :]
            imag = axis.pcolor(meanfrst, meanseco, hist.T, cmap='Blues', label=gdat.legdsamp, alpha=gdat.alphmrkr)
        else:
            varbfrst = getattr(gdatmodi, 'this' + strgfrst)[l] * gdat.dictglob['fact' + strgfrst + 'plot']
            varbseco = getattr(gdatmodi, 'this' + strgseco)[l] * gdat.dictglob['fact' + strgseco + 'plot']
            axis.scatter(varbfrst, varbseco, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
    
    # true
    if gdat.trueinfo:
        try:
            truevarbfrst = getattr(gdat, 'true' + strgfrst)[l] * gdat.dictglob['fact' + strgfrst + 'plot']
            truevarbseco = getattr(gdat, 'true' + strgseco)[l] * gdat.dictglob['fact' + strgseco + 'plot']
            axis.scatter(truevarbfrst, truevarbseco, alpha=gdat.alphmrkr, color='g', label=gdat.truelabl, s=sizelarg)
        except:
            pass

    # experimental
    if gdat.datatype == 'mock':
        try:
            exprvarbfrsttotl = getattr(gdat, 'expr' + strgfrst + 'totl')[l] * gdat.dictglob['fact' + strgfrst + 'plot']
            exprvarbsecototl = getattr(gdat, 'expr' + strgseco + 'totl')[l] * gdat.dictglob['fact' + strgseco + 'plot']
            axis.scatter(exprvarbfrsttotl, exprvarbsecototl, alpha=0.05, color='r', label=gdat.nameexpr + ' All', s=sizesmll)
        except:
            pass
        try:
            exprvarbfrst = getattr(gdat, 'expr' + strgfrst)[l] * gdat.dictglob['fact' + strgfrst + 'plot']
            exprvarbseco = getattr(gdat, 'expr' + strgseco)[l] * gdat.dictglob['fact' + strgseco + 'plot']
            axis.scatter(exprvarbfrst, exprvarbseco, alpha=0.3, color='r', label=gdat.nameexpr, s=sizelarg)
        except:
            pass

    scalfrst = gdat.dictglob['scal' + strgfrst + 'plot']
    scalseco = gdat.dictglob['scal' + strgseco + 'plot']
    if scalfrst == 'logt':
        axis.set_xscale('log')
    if scalseco == 'logt':
        axis.set_yscale('log')

    axis.set_xlabel(gdat.lablfeattotl[strgfrst])
    axis.set_ylabel(gdat.lablfeattotl[strgseco])
    limtfrst = gdat.dictglob['limt' + strgfrst + 'plot']
    limtseco = gdat.dictglob['limt' + strgseco + 'plot']
    axis.set_xlim(limtfrst * gdat.dictglob['fact' + strgfrst + 'plot'])
    axis.set_ylim(limtseco * gdat.dictglob['fact' + strgseco + 'plot'])
    make_legd(axis, loca=2)

    plt.tight_layout()
    if strg == 'post':
        strgmometemp = strgmome
    else:
        strgmometemp = ''
    path = retr_plotpath(gdat, gdatmodi, strg, '%s%s%s%spop%d' % (strgmometemp, strgplottype, strgfrst, strgseco, l))
    figr.savefig(path)
    plt.close(figr)
    

def plot_gene(gdat, gdatmodi, strg, strgydat, strgxdat, indxydat=None, strgindxydat=None, indxxdat=None, strgindxxdat=None, \
                     scal=None, scalxaxi=None, scalyaxi=None, limtxdat=None, limtydat=None, \
                     lablxaxi='', lablyaxi='', factxdat=1., factydat=1., histodim=False, ylim=None, strgxaxitwin=None, lablxaxitwin=None, offslegd=None, tdim=False):
   
    if scal == None:
        if scalxaxi == None:
            scalxaxi = 'linr'
        if scalyaxi == None:
            scalyaxi = 'linr'
    else:
        scalxaxi = scal
        scalyaxi = scal

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    
    if tdim:
        xdat = retr_fromgdat(gdat, gdatmodi, strg, strgxdat) * factxdat
        ydat = retr_fromgdat(gdat, gdatmodi, strg, strgydat) * factydat
    else:
        xdat = getattr(gdat, strgxdat + 'plot') * factxdat
        ydat = retr_fromgdat(gdat, gdatmodi, strg, strgydat) * factydat
    
    if indxxdat != None:
        xdat = xdat[indxxdat]
    if indxydat != None:
        ydat = ydat[indxydat]
    
    if tdim:
        axis.scatter(xdat, ydat, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
    else:
        if histodim:
            deltxdat = getattr(gdat, 'delt' + strgxdat[4:] + 'plot') * factxdat

            if scalxaxi == 'logt':
                xdattemp = 10**(log10(xdat) - diff(log10(xdat))[0] / 2.)
            else:
                xdattemp = xdat - deltxdat / 2.
    
    if strg == 'post':
        yerr = retr_fromgdat(gdat, gdatmodi, strg, strgydat, errr=True) * factydat
        
        if indxydat != None:
            yerr = yerr[[slice(None)] + indxydat]
        
        # label
        if strgydat.endswith('hist'):
            ##  element distribution
            labl = gdat.legdsampdist + ' distribution'
        else:
            ##  other
            labl = gdat.legdsampdist + ' distribution'
            
        temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color='black', label=labl, lw=1, capsize=10)
        for caps in listcaps:
            caps.set_markeredgewidth(1)

    else:
        if histodim:
            axis.bar(xdattemp, ydat, deltxdat, label=gdat.legdsamp, alpha=0.5)
        else:
            axis.plot(xdat, ydat, label=gdat.legdsamp, alpha=0.5)

    try:
        if strgxdat[4:] == strgydat[:4]:
            strgtemp = strgydat[:4]
            bins = copy(getattr(gdat, 'bins' + strgtemp + 'plot'))
            varb = copy(getattr(gdat, 'expr' + strgtemp))
            if strgtemp == 'lgal' or strgtemp == 'bgal':
            	bins *= gdat.anglfact
            	varb *= gdat.anglfact
            axis.hist(varb, bins, label=gdat.strgcatl, alpha=0.5, color='r')
    except:
        pass
    
    if gdat.datatype == 'mock':
        ydat = getattr(gdat, 'true' + strgydat)
        if indxydat != None:
            ydat = ydat[indxydat]
        if histodim:
            axis.bar(xdattemp, ydat, deltxdat, color='g', label='True', alpha=0.5)
        else:
            axis.plot(xdat, ydat, color='g', label='True', alpha=0.5)

    if scalxaxi == 'logt':
        axis.set_xscale('log')
    if scalyaxi == 'logt':
        if where(ydat > 0.)[0].size > 0:
            axis.set_yscale('log')
    
    if ylim != None:
        axis.set_ylim(ylim)

    axis.set_xlabel(lablxaxi)
    axis.set_ylabel(lablyaxi)

    # add another horizontal axis for counts
    if strgxaxitwin != None:
        axistwin = axis.twiny()
        if scalxaxi == 'logt':
            axistwin.set_xscale('log')
        axistwin.set_xlabel(lablxaxitwin)
        axistwin.spines['bottom'].set_position(('axes', 1.))
        binstwin = getattr(gdat, 'bins' + strgxaxitwin)
        minm = amin(binstwin) 
        maxm = amax(binstwin) 
        axistwin.set_xlim([minm, maxm])
        axistwin.xaxis.set_ticks_position('bottom')
        axistwin.xaxis.set_label_position('bottom')
       
    # superimpose prior on the feature
    if strgydat.endswith('hist') and strgydat[:4] in gdat.liststrgfeatprio[indxydat[0]]:
        xdatprio = getattr(gdat, strgxdat) * factxdat
        if strg == 'true' or strg == 'post':
            gdattemp = gdat
        else:
            gdattemp = gdatmodi
    
        if gdat.datatype == 'mock':
            truexdatprio = getattr(gdat, 'true' + strgxdat) * factxdat
            trueydatsupr = getattr(gdat, 'true' + strgydat + 'prio')[indxydat]
            axis.plot(truexdatprio, trueydatsupr, ls='-', alpha=gdat.alphmrkr, color='g')

        if strg == 'post':
            ydatsupr = getattr(gdattemp, 'medi' + strgydat + 'prio')[indxydat]
            yerrsupr = getattr(gdattemp, 'errr' + strgydat + 'prio')[[slice(None)] + indxydat]
            labl = gdat.legdsampdist + ' hyper-distribution'
            tdpy.util.plot_braz(axis, xdatprio, ydatsupr, yerr=yerrsupr, lcol='lightgrey', dcol='grey', labl=labl)
        else:
            ydatsupr = getattr(gdattemp, strg + strgydat + 'prio')[indxydat]
            axis.plot(xdatprio, ydatsupr, ls='--', alpha=gdat.alphmrkr, color='b')

    if indxydat != None:
        strgydat += strgindxydat
    
    if indxxdat != None:
        strgxdat += strgindxxdat
    
    if limtxdat != None:
        axis.set_xlim(limtxdat)
    else:
        axis.set_xlim([amin(xdat), amax(xdat)])
    if limtydat != None:
        axis.set_ylim(limtydat)
    else:
        if histodim:
            axis.set_ylim(gdat.histylim)
        else:
            axis.set_ylim([amin(ydat), amax(ydat)])
        
    make_legd(axis, offs=offslegd)
    
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, strgydat)
    figr.savefig(path)
    plt.close(figr)


def plot_scatspec(gdat, l, gdatmodi, plotdiff=False):
    
    figr, axrw = plt.subplots(1, gdat.numbener, figsize=(gdat.plotsize * gdat.numbener, gdat.plotsize))
    if gdat.numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        # prepare data to be plotted
        xdat = copy(gdat.truespec[l][0, i, :])
        
        xerr = tdpy.util.retr_errrvarb(gdat.truespec[l][:, i, :])
        
        yerr = zeros((2, xdat.size))
        if gdatmodi == None:
            ydat = copy(gdat.medispecassc[l][i, :])
            yerr = copy(gdat.errrspecassc[l][:, i, :])
        else:
            ydat = gdatmodi.thisspecassc[l][i, :]

        if (ydat == 0.).all():
            continue
            
        xdat *= gdat.fluxfactplot
        xerr *= gdat.fluxfactplot
        if plotdiff:
            ydat = 100. * (ydat - xdat) / xdat
        else: 
            ydat *= gdat.fluxfactplot
            yerr *= gdat.fluxfactplot
        
        # plot all associations
        indx = where(ydat > 0.)[0]
        if indx.size > 0:
            axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='black', alpha=0.1)
       
        # plot associations inside the comparison area
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, xerr=xerr, lw=1, marker='o', markersize=5, color='black')
        
        # plot associations to multiple model point sources
        if gdatmodi != None:
            indxpntsasscmult = gdatmodi.trueindxpntsasscmult[l]
            if len(indxpntsasscmult) > 0:
                axis.errorbar(xdat[indxpntsasscmult], ydat[indxpntsasscmult], ls='', yerr=yerr[:, indxpntsasscmult], xerr=xerr[:, indxpntsasscmult], \
                                                                                                                        lw=1, marker='o', markersize=5, color='r')
        else:
            indxpntsasscmult = array([])

        if plotdiff:
            axis.axhline(0., ls='--', alpha=gdat.alphmrkr, color='black')
        else:
            if gdat.elemtype == 'lght':
                # superimpose the bias line
                fluxbias = retr_fluxbias(gdat, gdat.binsspecplot[i, :], i)
                axis.plot(gdat.fluxfactplot * gdat.binsspecplot[i, :], gdat.fluxfactplot * gdat.binsspecplot[i, :], ls='--', alpha=gdat.alphmrkr, color='black')
                axis.plot(gdat.fluxfactplot * gdat.binsspecplot[i, :], gdat.fluxfactplot * fluxbias[0, :], ls='--', alpha=gdat.alphmrkr, color='black')
                axis.plot(gdat.fluxfactplot * gdat.binsspecplot[i, :], gdat.fluxfactplot * fluxbias[1, :], ls='--', alpha=gdat.alphmrkr, color='black')
        
        axis.set_xlabel(r'$%s^{%s}$%s' % (gdat.lablfeat['flux'], gdat.strgcatl, gdat.lablfeatunit['flux']))
        if i == 0:
            axis.set_ylabel('$%s^{samp}$%s' % (gdat.lablfeat['flux'], gdat.lablfeatunit['flux']))
        if i == gdat.numbener / 2:
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='black', alpha=0.1)
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='black')
            axis.errorbar(0., 0., ls='', yerr=0., xerr=0., lw=1, marker='o', markersize=5, color='r')
            # temp
            #make_legd(axis, loca=2)
        if not plotdiff:
            axis.set_yscale('log')
        axis.set_xscale('log')
        if gdat.enerbins:
            axis.set_title(gdat.strgener[i])
        if plotdiff:
            limsyaxi = array([-100., 100.])
        else:
            limsyaxi = gdat.fluxfactplot * array([gdat.minmspecplot[i], gdat.maxmspecplot[i]])
        axis.set_ylim(limsyaxi)
        axis.set_xlim([gdat.fluxfactplot * gdat.minmspecplot[i], gdat.fluxfactplot * gdat.maxmspecplot[i]])
   
    plt.tight_layout()
    if plotdiff:
        strg = 'diff'
    else:
        strg = ''
    path = retr_plotpath(gdat, gdatmodi, 'post', 'scatspec%spop%d' % (strg, l))
    figr.savefig(path)
    plt.close(figr)


def make_legd(axis, offs=None, loca=1, numbcols=1):

    legd = axis.legend(fancybox=True, frameon=True, bbox_to_anchor=offs, bbox_transform=axis.transAxes, ncol=numbcols, loc=loca)
    
    legd.get_frame().set_fill(True)
    legd.get_frame().set_facecolor('white')


def plot_scatpixl(gdat, gdatmodi, strg):
    
    figr, axgr = plt.subplots(gdat.numbevtt, gdat.numbener, figsize=(gdat.numbener * gdat.plotsize, gdat.numbevtt * gdat.plotsize))
    if gdat.numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if gdat.numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            ydat = retr_fromgdat(gdat, gdatmodi, strg, 'modlcnts')[i, gdat.indxpixlsave, m]
            
            axis.scatter(gdat.datacnts[i, gdat.indxpixlsave, m], ydat, alpha=gdat.alphmrkr)

            axislimt = [0., amax(gdat.datacnts[i, gdat.indxpixlsave, m]) * 1.1]
            axis.set_xlim(axislimt)
            axis.set_ylim(axislimt)
            
            slop, intc, coef, pval, stdv = sp.stats.linregress(gdat.datacnts[i, gdat.indxpixlsave, m], ydat)
            axis.text(0.5, 0.9, '$m = %.3g, b = %.3g, r = %.3g$' % (slop, intc, coef), va='center', ha='center', transform=axis.transAxes, fontsize=16)
            
            if m == gdat.numbevtt - 1:
                axis.set_xlabel('Data Counts')
            if i == 0:
                labl = 'Model Counts'
                if gdat.evttbins:
                    labl += ', ' + gdat.strgevtt[m]
                axis.set_ylabel(labl)
            if gdat.enerbins and m == 0:
                axis.set_title(gdat.strgener[i])
            
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatpixl')
    figr.savefig(path)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 41
    numbfluxprox = len(gdat.indxpixlprox)
    bins = empty((numbfluxprox, numbbins))
    indxpixlproxsize = empty((numbfluxprox, gdat.numbpixl))
    for h in range(gdat.numbfluxprox):
        for j in gdat.indxpixl:
            try:
                indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
            except:
                indxpixlproxsize[h, j] = gdat.numbpixl
        bins[h, :] = logspace(log10(amin(indxpixlproxsize[h, :])), log10(amax(indxpixlproxsize[h, :])), numbbins)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbfluxprox):
        axis.hist(indxpixlproxsize[h, :], bins=bins[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphmrkr)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixl, label='ROI', ls='--')
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    make_legd(axis)
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'init/indxprox.pdf')
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

    make_legd(axis)
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
    figr.savefig(gdat.pathplot + 'evidtest.pdf')
    plt.close(figr)
    
    
def plot_postbindmaps(gdat, indxpopltemp, strgbins, strgfeatsign=None):
    
    if strgfeatsign != None:
        numbparaplot = getattr(gdat, 'numb' + strgfeatsign + 'plot')
    else:
        numbparaplot = 1

    if strgbins == 'cumu':
        numbrows = 1
        numbcols = 1
    else:
        numbcols = 2
        if strgbins == 'full':
            numbrows = numbparaplot / 2
        else:
            numbrows = 2
    
    if strgfeatsign != None:
        indxfeat = gdat.liststrgfeat.index(strgfeatsign)
    else:
        indxfeat = arange(len(gdat.liststrgfeat))

    print 'gdat.pntsprob'
    print gdat.pntsprob.shape
    print 'strgbins'
    print strgbins
    print 'strgfeatsign'
    print strgfeatsign

    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
    if numbrows == 1:
        axgr = [axgr]            
    for a, axrw in enumerate(axgr):
        if numbcols == 1:
            axrw = [axrw]
        for b, axis in enumerate(axrw):
            h = a * 2 + b
            if strgbins == 'full':
                indxlowr = h
                indxuppr = h + 1
            elif strgbins == 'cumu':
                indxlowr = 0
                indxuppr = numbparaplot
            else:
                if h < 3:
                    indxlowr = 2 * h
                    indxuppr = 2 * (h + 1)
                else:
                    indxlowr = 2 * h
                    indxuppr = numbparaplot
            
            print 'indxlowr'
            print indxlowr
            print 'indxuppr'
            print indxuppr
            if strgfeatsign != None:
                print 'indxfeat'
                print indxfeat
                temp = sum(gdat.pntsprob[indxpopltemp, :, :, indxlowr:indxuppr, indxfeat], 2).T
            else:
                temp = sum(sum(gdat.pntsprob[indxpopltemp, :, :, indxlowr:indxuppr, :], 2), 2).T
                
            if where(temp > 0.)[0].size > 0:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
            else:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                
            if strgfeatsign != None:
                bins = getattr(gdat, 'bins' + strgfeatsign + 'plot')
            
            # superimpose true PS
            # temp
            try:
                if strgfeatsign != None:
                    truefeat = getattr(gdat, 'true' + strgfeatsign)[indxpopltemp] 
                    indxpnts = where((bins[indxlowr] < truefeat) & (truefeat < bins[indxuppr]))[0]
                else:
                    indxpnts = arange(gdat.trueflux[indxpopltemp].size)

                print 'true'
                print 'indxpnts'
                print indxpnts
                print 'bins[indxlowr]'
                print bins[indxlowr]
                print 'bins[indxuppr]'
                print bins[indxuppr]
                print

                mrkrsize = retr_mrkrsize(gdat, gdat.trueflux[indxpopltemp][indxpnts])
                axis.scatter(gdat.anglfact * gdat.truelgal[indxpopltemp][indxpnts], gdat.anglfact * gdat.truebgal[indxpopltemp][indxpnts], \
                                                                                        s=mrkrsize, alpha=gdat.alphmrkr, marker='*', lw=2, color='g')
            except:
                pass

            if a == numbrows - 1:
                axis.set_xlabel(gdat.lablfeattotl['lgal'])
            else:
                axis.set_xticklabels([])
            if b == 0:
                axis.set_ylabel(gdat.lablfeattotl['bgal'])
            else:
                axis.set_yticklabels([])

            draw_frambndr(gdat, axis)
            
            if strgbins != 'cumu':
                lablfeat = gdat.lablfeat[strgfeatsign]
                factfeat = gdat.dictglob['fact' + strgfeatsign + 'plot']
                titl = tdpy.util.mexp(factfeat * bins[indxlowr]) + ' < $%s$ < ' % lablfeat + tdpy.util.mexp(factfeat * bins[indxuppr])
                axis.set_title(titl)
    
    if strgfeatsign != None:
        lablfeattotl = gdat.lablfeattotl[strgfeatsign]
        plt.figtext(0.5, 0.95, '%s' % lablfeattotl, ha='center', va='center')
    axiscomm = figr.add_axes([0.9, 0.2, 0.02, 0.6])
    cbar = figr.colorbar(imag, cax=axiscomm)

    plt.subplots_adjust(left=0.18, top=.9, right=0.82, bottom=0.15, hspace=0.08, wspace=0.08)
    if strgbins == 'cumu':
        strgtemp = ''
    else:
        strgtemp = strgfeatsign
    path = gdat.pathpost + 'postbindmaps%s%s%d' % (strgbins, strgtemp, indxpopltemp) + '.pdf'
    figr.savefig(path)
    plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.binsanglplot)

    figr, axgr = plt.subplots(1, 2, figsize=(2 * gdat.plotsize, gdat.plotsize))
    figr.suptitle('King Function', fontsize=20)
    for k, axis in enumerate(axgr):
        if k == 0:
            sigmlist = [0.25]
            gammlist = [1.01, 2.5, 10.]
        else:
            sigmlist = [0.1, 0.25, 1.]
            gammlist = [2.]
        for sigm in sigmlist:
            for gamm in gammlist:
                axis.plot(angl, retr_singking(angl, sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        make_legd(axis)
        axis.set_yscale('log')
        axis.set_xlabel(gdat.lablfeattotl['gang'])
        axis.set_xlabel(r'$\mathcal{K}$')
        
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'king.pdf')
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
            minmdatacntstemp = 0.
            gdat.binscntstemp = linspace(minmdatacntstemp, maxmdatacntstemp, 20)
            meancntstemp = (gdat.binscntstemp[1:] + gdat.binscntstemp[0:-1]) / 2.
            diffcntstemp = gdat.binscntstemp[1:] - gdat.binscntstemp[0:-1]
            
            datacntshist = axis.hist(datacntstemp, gdat.binscntstemp, color='b')[0]

            init = [meancntstemp[argmax(datacntshist)], 1.]
            
            axis.set_xlim([amin(gdat.binscntstemp), amax(gdat.binscntstemp)])
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            if m == 0 and gdat.enerbins:
                axis.set_title(gdat.strgener[i])
            if i == 0 and gdat.evttbins:
                axis.set_ylabel(gdat.strgevtt[m])
            axis.set_yscale('log')

    plt.tight_layout()
    figr.savefig(gdat.pathinit + 'datacntshist.pdf')
    plt.close(figr)
    
    
def plot_intr(gdat):
    
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
        figr.savefig(gdat.pathimag + 'talkintr.pdf', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_eval(gdat):

    if gdat.exproaxitype:
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
        axis.plot(gdat.binsanglplot * gdat.anglfact, gdat.binsfluxprox[k] * psfntemp, label=labl, color=colr, alpha=alph)
        axis.set_xlim([amin(gdat.binsanglplot) * gdat.anglfact, amax(gdat.binsanglplot) * gdat.anglfact])
        if k > 0:
            axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
    axis.set_yscale('log')
    axis.set_xlabel(gdat.lablfeattotl['gang'])
    axis.set_ylabel(gdat.lablfeattotl['flux'])

    limt = gdat.specfraceval * amax(gdat.binsfluxprox[0] * psfntemp)
    maxmangltemp = interp(1e-1 * limt, gdat.binsfluxprox[k] * psfntemp[::-1], gdat.binsanglplot[::-1] * gdat.anglfact)
    
    axis.set_ylim([1e-3 * limt, None])
    if limt > 0:
        axis.axhline(limt, color='red', ls=':', label='Flux floor')
    axis.set_xlim([None, maxmangltemp])
    
    make_legd(axis)

    plt.tight_layout()
    figr.savefig(gdat.pathinit + 'eval.pdf')
    plt.close(figr)


def plot_mosa(gdat):

    # empty global object
    gdatmodi = tdpy.util.gdatstrt()
    gdatmodi.thischroproc = zeros(gdat.numbchroproc)

    # data structure to hold the indices of model PS to be compared to the reference catalog 
    gdatmodi.indxmodlpntscomp = [[] for l in gdat.indxpopl]
    
    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamptotl:
        indxsampmosa = choice(gdat.indxsamptotl, size=numbsampmosa, replace=False)
        for l in gdat.indxpopl:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            
                            n = indxsampmosa[numbcols*a+b]
                            gdatmodi.thissampvarb = gdat.listsampvarb[n, :].flatten()
                            gdatmodi.thissamp = gdat.listsamp[n, :].flatten()
                            gdatmodi.thisindxpntsfull = gdat.listindxpntsfull[n]
                            
                            if gdat.numbtrap > 0:
                                # temp
                                #gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
                                proc_samp(gdat, gdatmodi, 'this')

                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.lablfeattotl['lgal'])
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.lablfeattotl['bgal'])
                            else:
                                axis.set_yticklabels([])
                            
                            imag = retr_imag(gdat, axis, gdat.datacnts, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                            supr_fram(gdat, gdatmodi, 'this', axis, l)
                    
                    if gdat.enerbins:
                        plt.figtext(0.5, 0.93, gdat.strgener[i], ha='center', va='center')
                    axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
                    cbar = figr.colorbar(imag, cax=axiscomm)
                    cbar.set_ticks(gdat.tickdatacnts)
                    cbar.set_ticklabels(gdat.labldatacnts)
                    plt.subplots_adjust(left=0.1, top=.91, hspace=0.01, wspace=0.05, bottom=0.09)
                    if l == 1:
                        strg = ''
                    else:
                        strg = 'pop%d' % l
                    if m == None:
                        path = gdat.pathpost + 'mosa' + strg + 'ene%dA.pdf' % (gdat.indxenerincl[i])
                    else:
                        path = gdat.pathpost + 'mosa' + strg + 'ene%devtt%d.pdf' % (gdat.indxenerincl[i], gdat.indxevttincl[m])
                    figr.savefig(path)
                    plt.close(figr)
    else:
        if gdat.verbtype > 0:
            print 'Skipping the mosaic plot...'


def plot_grap(plottype='igal', verbtype=0):
        
    figr, axis = plt.subplots(figsize=(6, 6))

    grap = nx.DiGraph()
    if plottype == 'ngal':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']

    if plottype == 'igal':
        listcolr = ['black', 'black', 'black', 'olive', 'black', 'black', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'olive', \
                                                                                                    'magenta', 'magenta', 'magenta', 'magenta', 'magenta']
    if plottype == 'chan':
        listcolr = ['black', 'black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']

    if plottype == 'lens':
        listcolr = ['olive', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', 'magenta']

    if plottype == 'lensprim':
        listcolr = ['olive', 'black', 'black', 'olive', 'black', 'olive', 'olive', 'olive', 'olive', 'olive', 'black', 'black', 'olive', 'magenta', 'magenta', \
                                                                                                                                        'magenta', 'magenta', 'magenta']

    grap.add_edges_from([ \
                         ('meanpnts', 'numbpnts'), \
                         ('modl','data'), \
                         ('psfp', 'modl'), \
                         ('bacp', 'modl'), \
                         ('lgal','modl'), \
                         ('bgal','modl'), \
                         ('numbpnts','lgal'), \
                         ('numbpnts','bgal'), \
                        ])
    
    if plottype == 'igal' or plottype == 'ngal' or plottype == 'chan':
        grap.add_edges_from([ \
                             ('ampldistslop', 'ampl'), \
                             ('ampl', 'modl'), \
                             ('numbpnts','ampl'), \
                             ('numbpnts', 'sind'), \
                             ('sind','modl'), \
                            ])
    if plottype == 'lens' or plottype == 'lensprim':
        grap.add_edges_from([ \
                             ('lenp', 'modl'), \
                             ('defsdistslop', 'defs'), \
                             ('defs', 'modl'), \
                             ('numbpnts','defs'), \
                            ])
    if plottype == 'lensprim':
        grap.add_edges_from([ \
                             ('ascadistslop', 'asca'), \
                             ('asca', 'modl'), \
                             ('numbpnts','asca'), \
                             ('acutdistslop', 'acut'), \
                             ('acut', 'modl'), \
                             ('numbpnts','acut'), \
                            ])
    
    if plottype == 'igal' or plottype == 'chan':
        grap.add_edges_from([ \
                             ('sinddistslop', 'sind'), \
                             ('expodistslop', 'expo'), \
                            ])
    
    if plottype == 'igal':
        grap.add_edges_from([ \
                             ('numbpnts', 'expo'), \
                             ('expo', 'modl'), \
                             ('spatdistslop', 'lgal'), \
                             ('spatdistslop', 'bgal'), \
                            ])
        
    labl = {}
    labl['numbpnts'] = '$N$'
    
    if plottype == 'igal':
        labl['ampldistslop'] = r'$\vec{\alpha}$'
        labl['meanpnts'] = r'$\vec{\mu}$'
    else:
        labl['meanpnts'] = r'$\mu$'
    if plottype == 'chan' or plottype == 'ngal':
        labl['ampldistslop'] = r'$\alpha$'
    if plottype == 'lens' or plottype == 'lensprim':
        labl['defsdistslop'] = r'$\alpha$'
    if plottype == 'lensprim':
        labl['ascadistslop'] = r'$\lambda_s$'
        labl['acutdistslop'] = r'$\lambda_\tau$'
    
    if plottype == 'igal':
        labl['expodistslop'] = r'$\vec{\tau_{E_c}}$'
        labl['sinddistslop'] = r'$\vec{\tau_s}$'
    if plottype == 'chan':
        labl['sinddistslop'] = r'$\beta$'
    
    if plottype == 'igal':
        labl['spatdistslop'] = r'$\vec{\gamma}$'
    if plottype == 'lens' or plottype == 'lensprim':
        labl['lenp'] = r'$\vec{\chi}$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['bacp'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{\theta_1}$'
    labl['bgal'] = r'$\vec{\theta_2}$'
    if plottype == 'igal' or plottype == 'ngal' or plottype == 'chan':
        labl['sind'] = r'$\vec{s}$'
        labl['ampl'] = r'$\vec{f}$'
    else:
        labl['defs'] = r'$\vec{\alpha_s}$'
    if plottype == 'lensprim':
        labl['asca'] = r'$\vec{\theta_s}$'
        labl['acut'] = r'$\vec{\tau}$'
        
    if plottype == 'igal':
        labl['expo'] = r'$\vec{E_c}$'
    labl['modl'] = r'$\mathcal{M}$'
    labl['data'] = r'$\mathcal{D}$'
    
    posi = nx.circular_layout(grap)
    posi['sinddistslop'] = array([0.4, 0.15])
    if plottype == 'igal':
        posi['expodistslop'] = array([0.6, 0.15])
        posi['spatdistslop'] = array([-0.2, 0.15])
    if plottype == 'igal' or plottype == 'ngal' or plottype == 'chan':
        posi['numbpnts'] = array([0., 0.075])
        posi['meanpnts'] = array([0., 0.15])
        posi['ampldistslop'] = array([0.2, 0.15])
    if plottype == 'lens' or plottype == 'lensprim':
        posi['numbpnts'] = array([-0.1, 0.075])
        posi['meanpnts'] = array([-0.1, 0.15])
        posi['defsdistslop'] = array([0.1, 0.15])
    if plottype == 'lensprim':
        posi['ascadistslop'] = array([0.3, 0.15])
        posi['acutdistslop'] = array([0.5, 0.15])
    if plottype == 'igal':
        posi['psfp'] = array([0.9, -0.0])
        posi['bacp'] = array([0.7, -0.0])
    if plottype == 'ngal':
        posi['psfp'] = array([0.7, -0.0])
        posi['bacp'] = array([0.5, -0.0])
    if plottype == 'chan':
        posi['psfp'] = array([0.7, -0.0])
        posi['bacp'] = array([0.5, -0.0])
    if plottype == 'lens':
        posi['psfp'] = array([0.3, -0.0])
        posi['bacp'] = array([0.5, -0.0])
        posi['lenp'] = array([0.7, -0.0])
    if plottype == 'lensprim':
        posi['psfp'] = array([0.7, -0.0])
        posi['bacp'] = array([0.9, -0.0])
        posi['lenp'] = array([1.1, -0.0])
    posi['lgal'] = array([-0.3, -0.0])
    posi['bgal'] = array([-0.1, -0.0])
    if plottype == 'igal' or plottype == 'ngal' or plottype == 'chan':
        posi['sind'] = array([0.3, -0.0])
        posi['ampl'] = array([0.1, -0.0])
    if plottype == 'igal':
        posi['expo'] = array([0.5, -0.0])
    if plottype == 'lensprim' or plottype == 'lens':
        posi['defs'] = array([0.1, -0.0])
    if plottype == 'lensprim':
        posi['asca'] = array([0.3, -0.0])
        posi['acut'] = array([0.5, -0.0])
    posi['modl'] = array([0., -0.075])
    posi['data'] = array([0., -0.15])
   
    if verbtype > 0:
        numb = max(len(grap.edges()), len(listcolr))
        for k in range(numb):
            try:
                print '%15s %15s %15s' % (grap.edges()[k][0], grap.edges()[k][1], listcolr[k])
            except:
                print 'unequal'
        print

    size = 500
    nx.draw(grap, posi, labels=labl, ax=axis, edgelist=[], nodelist=[])
    nx.draw_networkx_edges(grap, posi, ax=axis, labels=labl, edge_color=listcolr)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['modl', 'data'], node_color='grey', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbpnts'], node_color='b', node_size=size)
    if plottype == 'igal':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'spatdistslop', 'ampldistslop', 'sinddistslop', 'expodistslop'], \
                                                                                                                                                node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind', 'expo'], node_color='g', node_size=size)
    if plottype == 'chan':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'ampldistslop', 'sinddistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind'], node_color='g', node_size=size)
    if plottype == 'ngal':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'ampldistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind'], node_color='g', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['psfp', 'bacp'], node_color='y', node_size=size)
    if plottype == 'lens' or plottype == 'lensprim':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'defsdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lenp'], node_color='y', node_size=size)
    if plottype == 'lens':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs'], node_color='g', node_size=size)
    if plottype == 'lensprim':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['ascadistslop', 'acutdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs', 'asca', 'acut'], node_color='g', node_size=size)
    
    pathplot = os.environ["PCAT_DATA_PATH"] + '/imag/'
    plt.tight_layout()
    figr.savefig(pathplot + 'grap%s.pdf' % plottype)
    plt.close(figr)


#plot_grap(verbtype=1)
#plot_grap(plottype='ngal', verbtype=1)
#plot_grap(plottype='lens', verbtype=1)
#plot_grap(plottype='lensprim', verbtype=1)
#plot_grap(plottype='chan', verbtype=1)


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
    axis.set_xlabel(gdat.lablfeattotl['lgal'])
    axis.set_ylabel(gdat.lablfeattotl['bgal'])

    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'thrs.pdf')
    plt.close(figr)
    

def make_anim(gdat):
    
    listfile = fnmatch.filter(os.listdir(gdat.pathfram), '*_swep*.pdf')
    listfiletemp = []
    for thisfile in listfile:
        listfiletemp.extend((thisfile.split('_')[0]).rsplit('/', 1))
    
    listname = list(set(listfiletemp))

    if gdat.verbtype > 0:
        print 'Making animations of frame plots...'
    for name in listname:
        
        strgtemp = '%s*_swep*.pdf' % name
        listfile = fnmatch.filter(os.listdir(gdat.pathfram), strgtemp)
        numbfile = len(listfile)
        liststrgextn = []
        for k in range(numbfile):
            liststrgextn.append((listfile[k].split(name)[1]).split('_')[0])
        
        liststrgextn = list(set(liststrgextn))
        
        for k in range(len(liststrgextn)):
    
            listfile = fnmatch.filter(os.listdir(gdat.pathfram), name + liststrgextn[k] + '_swep*.pdf')
            numbfile = len(listfile)
            
            indxfilelowr = int(ceil(numbfile * float(gdat.numbburn) / gdat.numbswep))
            if indxfilelowr < numbfile:
                indxfileanim = arange(indxfilelowr, numbfile)
            else:
                continue
            
            if gdat.verbtype > 0:
                print 'Making %s animation...' % name
                
            indxfileanim = choice(indxfileanim, replace=False, size=indxfileanim.size)
            
            cmnd = 'convert -delay 20 -density 300 -quality 100 '
            for n in range(indxfileanim.size):
                cmnd += '%s%s ' % (gdat.pathfram, listfile[indxfileanim[n]])
            cmnd += ' %s%s.gif' % (gdat.pathanim, name + liststrgextn[k])
            os.system(cmnd)

       
def plot_init(gdat):
        
    # make initial plots
    if gdat.makeplot:
        #plot_3fgl_thrs(gdat)
        plot_datacntshist(gdat)
        if gdat.evalcirc != 'full':
            plot_indxprox(gdat)
        #if gdat.exprtype == 'ferm':
        #    plot_fgl3(gdat)
        
        # temp
        if gdat.makeplotintr:
            plot_intr(gdat)
            #plot_pert()
            #plot_king(gdat)
    
            if gdat.elemtype == 'lens':
                xdat = gdat.binsanglplot * gdat.anglfact
                lablxdat = gdat.lablfeattotl['gang']
                
                # temp
                #angl = linspace(0.1, 2., 20)
                #defl = retr_deflcutf(angl, 0.1, 0.05, 0.3)
                #coef = polyfit(angl, defl, 4)

                listdeflscal = array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1]) / gdat.anglfact
                listanglscal = array([0.05, 0.05, 0.1, 0.05, 0.05, 0.05]) / gdat.anglfact
                listanglcutf = array([0.3, 0.3, 0.3, 1., 2., 0.]) / gdat.anglfact
                listasym = [False, False, False, False, False, True]
                listydat = []
                for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
                    listydat.append(retr_deflcutf(gdat.binsanglplot, deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
                
                listydat.append(xdat * 0. + 0.05)
                
                path = gdat.pathinitintr + 'deflcutf.pdf'
                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=r'$\alpha$ [$^{\prime\prime}$]')

                spec = 1e-19 # [erg/cm^2/s]
                listsize = array([0.1, 0.1, 0.1]) / gdat.anglfact
                listindx = array([3., 4., 5.])
                listydat = []
                listlegd = []
                for size, indx in zip(listsize, listindx):
                    listydat.append(retr_sersprof(spec, gdat.binsanglplot, size, indx=indx))
                    listlegd.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
                path = gdat.pathinitintr + 'sersprof.pdf'
                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=r'$f$', listlegd=listlegd)
            
                minmredshost = 0.01
                maxmredshost = 0.4
                minmredssour = 0.01
                maxmredssour = 2.
                numbreds = 200
                retr_axis(gdat, 'redshost', minmredshost, maxmredshost, numbreds)
                retr_axis(gdat, 'redssour', minmredssour, maxmredssour, numbreds)
                
                gdat.meanadishost = empty(numbreds)
                for k in range(numbreds):
                    gdat.meanadishost[k] = gdat.adisobjt(gdat.meanredshost[k]) * 1e3
                
                minmmass = zeros((numbreds + 1, numbreds + 1))
                maxmmass = zeros((numbreds + 1, numbreds + 1))
                for k, redshost in enumerate(gdat.binsredshost):
                    for n, redssour in enumerate(gdat.binsredssour):
                        if redssour > redshost:
                            adishost = gdat.adisobjt(redshost) * 1e3
                            adissour = gdat.adisobjt(redssour) * 1e3
                            adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
                            massfromdeflscal = retr_massfromdeflscal(gdat, adissour, adishost, adishostsour, gdat.anglscal, gdat.anglcutf)
                            minmmass[n, k] = log10(massfromdeflscal * gdat.minmflux)
                            maxmmass[n, k] = log10(massfromdeflscal * gdat.maxmflux)
                
                valulevl = linspace(6., 12., 20)
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                imag = axis.imshow(minmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=6, vmax=9)
                cont = axis.contour(gdat.binsredshost, gdat.binsredssour, minmmass, 10, colors='g', levels=valulevl)
                axis.clabel(cont, inline=1, fontsize=10)
                axis.set_xlabel('$z_h$')
                axis.set_ylabel('$z_s$')
                axis.set_title(r'$M$ [$M_{\odot}$]')
                path = gdat.pathinitintr + 'massredsminm.pdf'
                plt.colorbar(imag) 
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=6, vmax=10)
                cont = axis.contour(gdat.binsredshost, gdat.binsredssour, maxmmass, 10, colors='g', levels=valulevl)
                axis.clabel(cont, inline=1, fontsize=10)
                axis.set_xlabel('$z_h$')
                axis.set_ylabel('$z_s$')
                axis.set_title(r'$M$ [$M_{\odot}$]')
                path = gdat.pathinitintr + 'massredsmaxm.pdf'
                plt.colorbar(imag) 
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.plot(gdat.meanredshost, gdat.meanadishost * gdat.sizepixl)
                axis.plot(gdat.meanredshost, gdat.meanadishost * 2. * gdat.maxmgangdata)
                axis.set_xlabel('$z_h$')
                axis.set_yscale('log')
                axis.set_ylabel(r'$\lambda$ [kpc]')
                path = gdat.pathinitintr + 'wlenreds.pdf'
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
                
            if gdat.evalcirc and gdat.elemtype == 'lght':
                plot_eval(gdat)
            return
        
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            
            # temp
            if False and gdat.pixltype == 'cart' and gdat.elemtype == 'lght':
                figr, axis, path = init_figr(gdat, None, 'datacntspeak', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.datacnts, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                axis.scatter(gdat.anglfact * gdat.lgalcart[gdat.indxxaximaxm], gdat.anglfact * gdat.bgalcart[gdat.indxyaximaxm], alpha=0.6, s=20, edgecolor='none')
                
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
    
            if gdat.correxpo:
                figr, axis, path = init_figr(gdat, None, 'expo', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.expo, '', 'expo', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
        
            for c in gdat.indxback:
                figr, axis, path = init_figr(gdat, None, 'backcnts', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.backcnts[c], '', 'datacnts', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
    
            if gdat.numbback > 1:
                figr, axis, path = init_figr(gdat, None, 'backcntstotl', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.backcntstotl, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
        
            figr, axis, path = init_figr(gdat, None, 'diffcntstotl', '', indxenerplot=i, indxevttplot=m)
            imag = retr_imag(gdat, axis, gdat.datacnts - gdat.backcntstotl, '', 'datacnts', thisindxener=i, thisindxevtt=m)
            make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
            plt.tight_layout()
            figr.savefig(path)
            plt.close(figr)
    
            if gdat.elemtype == 'lens':
                
                figr, axis, path = init_figr(gdat, None, 'modlcntsraww', 'true', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.truemodlcntsraww, '', 'datacnts', thisindxener=i, thisindxevtt=m, tdim=True)
                make_cbar(gdat, axis, imag, 0, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)


def plot_defl(gdat, gdatmodi, strg, strgcomp='', indxdefl=None, thisindxpopl=-1):

    strgvarb = 'defl'
    if indxdefl != None:
        strgvarb += 'sing'
    strgvarb += strgcomp
    
    defl = retr_fromgdat(gdat, gdatmodi, strg, strgvarb)
   
    # temp
    #if strgcomp == 'resi':
    #    defl *= 1e2

    strgplot = strg + strgvarb

    if indxdefl != None:
        defl = defl[:, :, :, indxdefl]
        strgplot += '%04d' % indxdefl

    figr, axis, path = init_figr(gdat, gdatmodi, strgplot, strg, indxpoplplot=thisindxpopl)
    make_catllabl(gdat, strg, axis)
    draw_frambndr(gdat, axis)
  
    defllgal = defl[:, :, 0]
    deflbgal = defl[:, :, 1]
    fact = 4
    deflmagn = sqrt(defllgal[::fact, ::fact]**2 + deflbgal[::fact, ::fact]**2)
    axis.imshow(zeros((10, 10)))
    ptch = axis.quiver(gdat.anglfact * gdat.lgalgridcart[::fact, ::fact], gdat.anglfact * gdat.bgalgridcart[::fact, ::fact], \
                  gdat.fluxfactplot * defllgal[::fact, ::fact], gdat.fluxfactplot * deflbgal[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strg, axis)
    #plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    figr.savefig(path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strg, strgvarb, strgcbar=None, thisindxener=None, thisindxevtt=-1, tdim=False, thisindxpopl=-1, stdv=False):
    
    if strgcbar == None:
        strgcbar = strgvarb
   
    if strgvarb == 'datacnts':
        strgplot = strgvarb
    else:
        if strg == 'post':
            if stdv:
                strgtemp = 'stdv'
            else:
                strgtemp = 'medi'
        else:
            strgtemp = ''
        strgplot = strgtemp + strgvarb
    
    maps = retr_fromgdat(gdat, gdatmodi, strg, strgvarb, stdv=stdv)
    
    figr, axis, path = init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=thisindxener, indxevttplot=thisindxevtt, indxpoplplot=thisindxpopl)
    
    imag = retr_imag(gdat, axis, maps, strg, strgcbar, thisindxener=thisindxener, thisindxevtt=thisindxevtt, tdim=tdim)
    
    tick = getattr(gdat, 'tick' + strgcbar) 
    labl = getattr(gdat, 'labl' + strgcbar) 
    make_cbar(gdat, axis, imag, tick=tick, labl=labl)
    
    make_catllabl(gdat, strg, axis)
    supr_fram(gdat, gdatmodi, strg, axis, thisindxpopl)

    plt.tight_layout()
    figr.savefig(path)
    plt.close(figr)
    
