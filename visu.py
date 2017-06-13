# common imports
from __init__ import *

# internal functions
from util import *

def plot_samp(gdat, gdatmodi, strg):
   
    if gdat.shrtfram and strg == 'this':
        plot_genemaps(gdat, gdatmodi, strg, 'datacnts', thisindxener=0, thisindxevtt=0)
        plot_genemaps(gdat, gdatmodi, strg, 'resicnts', thisindxener=0, thisindxevtt=0)
    else:    
        if strg == 'true':
            strgmodl = strg
        else:
            strgmodl = 'fitt'
        
        if gdatmodi != None:
            gdatobjt = gdatmodi
        else:
            gdatobjt = gdat
        
        numbtrap = getattr(gdat, strgmodl + 'numbtrap')
        sampvarb = getattr(gdatobjt, strg + 'sampvarb')
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        numbdeflsubhplot = getattr(gdat, strgmodl + 'numbdeflsubhplot')
        numbdeflsingplot = getattr(gdat, strgmodl + 'numbdeflsingplot')
        if numbtrap > 0:
            liststrgfeatcorr = getattr(gdat, strgmodl + 'liststrgfeatcorr')
            spectype = getattr(gdat, strgmodl + 'spectype')
            indxpopl = getattr(gdat, strgmodl + 'indxpopl')
            numbpnts = sampvarb[getattr(gdat, strgmodl + 'indxfixpnumbpnts')].astype(int)
            liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
            if strg != 'post':
                indxelemfull = list(getattr(gdatobjt, strg + 'indxelemfull'))
                indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)
            
        if gdatmodi != None:
            strgswep = '_%09d' % gdatmodi.cntrswep
        else:
            strgswep = ''
        
        if gdatmodi != None:
            gdatobjt = gdatmodi
        else:
            gdatobjt = gdat
    
        # plots
        ## highest amplitude element
        # temp
        if numbtrap > 0:
            if False:
                if gdat.elemtype == 'lght':
                    plot_brgt(gdat, gdatmodi, strg)
        
            # completeness and false discovery rate
            if strg != 'true':
                for strgclas in ['cmpl', 'fdis']:   
                    nameinte = strgclas + '/'
                    for l in gdat.trueindxpopl:
                        lablydat = getattr(gdat, 'labl' + strgclas + 'pop%d' % l)
                        indxydat = [l, slice(None)]
                        strgindxydat = 'pop%d' % l
                        for strgfeat in gdat.trueliststrgfeatodim[l]:
                            if strgfeat in gdat.fittliststrgfeatodim[l]:
                                if not (gdat.datatype == 'inpt' and getattr(gdat, 'expr' + strgfeat) == None):
                                    factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                                    lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                                    scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                                    limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                                    plot_gene(gdat, gdatmodi, strg, strgclas + strgfeat, 'mean' + strgfeat, lablxaxi=lablxdat, \
                                              lablyaxi=lablydat, factxdat=factxdat, plottype='errr', \
                                              scalxaxi=scalxdat, limtydat=[0., 1.], limtxdat=limtxdat, indxydat=indxydat, \
                                              strgindxydat=strgindxydat, omittrue=True, nameinte=nameinte)
                
            # element features projected onto the data axes
            if gdat.datatype == 'mock':
                alph = 0.1
                if strg == 'this':
                    pathtemp = gdat.pathplot + gdat.namesampdist + '/fram/'
                    colr = 'b'
                elif strg == 'mlik':
                    pathtemp = gdat.pathplot + gdat.namesampdist + '/finl/'
                    colr = 'r'
                elif strg == 'true':
                    pathtemp = gdat.pathinit
                    colr = 'g'
                elif strg == 'post':
                    pathtemp = gdat.pathplot + gdat.namesampdist + '/finl/'
                    colr = 'black'
    
                # PS spectra
                if gdat.numbener > 1 and gdat.elemtype == 'lght':
                    
                    # temp
                    if not (strg == 'post' and not gdat.condcatl) and not strg == 'post':
                        if strg == 'post':
                            specplot = [empty((gdat.numbener, gdat.numbstkscond)) for l in indxpopl]
                            for l in indxpopl:
                                for r in gdat.indxstkscond:
                                    specplot[l][:, r] = gdat.dictglob['poststkscond'][r]['specplot'][l]
                        else:
                            specplot = getattr(gdatobjt, strg + 'specplot')
                        for l in indxpopl:
                            listxdat = []
                            listplottype = []
                            listydat = []
                            for k in range(numbpnts[l]):
                                specplottemp = specplot[l]
                                if strgmodl == 'true':
                                    specplottemp = specplottemp[0, :, k]
                                else:
                                    specplottemp = specplottemp[:, k]
                                if gdat.enerdiff:
                                    specplottemp *= gdat.meanenerplot**2
                                listplottype.append('line')
                                listxdat.append(gdat.meanenerplot)
                                listydat.append(specplottemp)
                            
                            path = pathtemp + strg + 'specpop%d%s.pdf' % (l, strgswep)
                            tdpy.util.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, colr=colr, alph=alph, \
                                  plottype=listplottype, \
                                  limtxdat=[gdat.minmener, gdat.maxmener], lablydat='$E^2dN/dE$ [%s]' % gdat.lablfluxunit, limtydat=[amin(gdat.minmspec), amax(gdat.maxmspec)])
                
                if gdat.elemtype == 'lens':
                    if not (strg == 'post' and not gdat.condcatl):
                        ## deflection profiles
                        if gdat.variasca and gdat.variacut:
                            xdat = gdat.meanangl * gdat.anglfact
                            lablxdat = gdat.lablgangtotl
                            deflprof = retr_fromgdat(gdat, gdatmodi, strg, 'deflprof')
                            for l in indxpopl:
                                listydat = []
                                listvlinfrst = []
                                listvlinseco = []
    
                                for k in arange(numbpnts[l]):
                                    if strgmodl == 'true':
                                        deflproftemp = deflprof[l][0, :, :]
                                    else:
                                        deflproftemp = deflprof[l]
                                    listydat.append(deflproftemp[:, k] * gdat.anglfact)
                                    listvlinfrst.append(sampvarb[indxsampcomp['asca'][l][k]] * gdat.anglfact) 
                                    listvlinseco.append(sampvarb[indxsampcomp['acut'][l][k]] * gdat.anglfact)
                                    
                                listydat.append(xdat * 0. + gdat.anglfact * sampvarb[getattr(gdat, strgmodl + 'indxfixpbeinhost')])
                                path = pathtemp + strg + 'deflsubhpop%d%s.pdf' % (l, strgswep)
                                limtydat = [1e-3, 1.]
                                limtxdat = [1e-3, 1.]
                                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, drawdiag=True, limtydat=limtydat, \
                                        limtxdat=limtxdat, colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
                    
                    ## radial mass budget
                    xdat = gdat.meananglhalf * gdat.anglfact
                    lablxdat = gdat.lablanglfromhosttotl
                    for namecalc in ['delt', 'intg']:
                        
                        # host and subhalo masses
                        limtydat = [gdat.minmmcut, getattr(gdat, 'maxmmasshost' + namecalc + 'bein')]
                        lablydat = getattr(gdat, 'lablmass' + namecalc + 'totl')

                        path = pathtemp + strg + 'mass%s%s.pdf' % (namecalc, strgswep)
                        listydat = [retr_fromgdat(gdat, gdatmodi, strg, 'masshost' + namecalc), \
                                                        retr_fromgdat(gdat, gdatmodi, strg, 'masssubh' + namecalc)]
                        if strg == 'post':
                            listyerr = [retr_fromgdat(gdat, gdatmodi, strg, 'masshost' + namecalc, mometype='errr'), \
                                                        retr_fromgdat(gdat, gdatmodi, strg, 'masssubh' + namecalc, mometype='errr')]
                        else:
                            listyerr = None
                        tdpy.util.plot_gene(path, xdat, listydat, yerr=listyerr, scalydat='logt', lablxdat=lablxdat, limtydat=limtydat, \
                                                                                                colr=colr, lablydat=lablydat, listlinestyl=['-', '--'])

                        # subhalo mass fraction
                        limtydat = [1e-3, 0.1]
                        lablydat = getattr(gdat, 'lablfracsubh' + namecalc + 'totl')

                        path = pathtemp + strg + 'fracsubh%s%s.pdf' % (namecalc, strgswep)
                        listydat = retr_fromgdat(gdat, gdatmodi, strg, 'fracsubh' + namecalc)
                        if strg == 'post':
                            listyerr = retr_fromgdat(gdat, gdatmodi, strg, 'fracsubh' + namecalc, mometype='errr')
                        else:
                            listyerr = None
                        tdpy.util.plot_gene(path, xdat, listydat, yerr=listyerr, scalydat='logt', lablxdat=lablxdat, limtydat=limtydat, \
                                                                                                                                  colr=colr, lablydat=lablydat)

            if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                ## PSF radial profile
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        if gdat.fittoaxitype:
                            indxydat = [i, slice(None), m, 0]
                        else:
                            indxydat = [i, slice(None), m]
                        strgindxydat = '%d%d' % (i, m)
                        lablxaxi = gdat.lablgangtotl
                        limtydat= array([1e-3, 1e3]) * gdat.anglfact**2
                        plot_gene(gdat, gdatmodi, strg, 'psfn', 'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalyaxi='logt', \
                                                                                       factxdat=gdat.anglfact, lablxaxi=lablxaxi, lablyaxi=r'$\mathcal{P}$', limtydat=limtydat)
                        if gdat.fittoaxitype or gdat.trueoaxitype:
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
            for l in indxpopl:
                if strg != 'true' and gdat.trueinfo:
                    for strgfeat in gdat.fittliststrgfeatodim[l]:
                        if gdat.datatype == 'mock':
                            plot_scatassc(gdat, gdatmodi, strg, l, strgfeat)
                            plot_scatassc(gdat, gdatmodi, strg, l, strgfeat, plotdiff=True)
                        else:
                            # temp
                            try:
                                plot_scatassc(gdat, gdatmodi, strg, l, strgfeat)
                                plot_scatassc(gdat, gdatmodi, strg, l, strgfeat, plotdiff=True)
                            except:
                                pass
                for a, strgfrst in enumerate(liststrgfeatcorr[l]):
                    for b, strgseco in enumerate(liststrgfeatcorr[l]):
                        if a < b:
                            for strgplottype in ['hist', 'scat']:
                                if strgmodl == 'true' and strgplottype == 'hist':
                                    continue
                                plot_elemtdim(gdat, gdatmodi, strg, l, strgplottype, strgfrst, strgseco)
    
            # element feature histograms
            if not (strg == 'true' and gdat.datatype == 'inpt'):
                for l in indxpopl:
                    indxydat = [l, slice(None)]
                    strgindxydat = 'pop%d' % l
    
                    for strgfeat in liststrgfeatodim[l] + gdat.liststrgfeatplot:
                        
                        scalxaxi = getattr(gdat, 'scal' + strgfeat + 'plot')
                        factxaxi = getattr(gdat, 'fact' + strgfeat + 'plot')
                        lablxaxi = getattr(gdat, 'labl' + strgfeat + 'totl')
                        limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxaxi, getattr(gdat, 'maxm' + strgfeat) * factxaxi]
                        
                        limtydat = gdat.limtpntshist
                    
                        # for true model, also plot the significant elements only
                        # temp
                        #if strgmodl == 'true' and strgfeat in gdat.listnamefeatpars:
                        #    listname = ['hist' + strgfeat, 'hist' + strgfeat + 'pars']
                        #else:
                        #    listname = ['hist' + strgfeat]
                        listname = ['hist' + strgfeat]
                        for name in listname:
                            plot_gene(gdat, gdatmodi, strg, name, 'mean' + strgfeat, scalyaxi='logt', lablxaxi=lablxaxi, lablyaxi=r'$N$', factxdat=factxaxi, histodim=True, \
                                              scalxaxi=scalxaxi, limtydat=limtydat, limtxdat=limtxdat, indxydat=indxydat, strgindxydat=strgindxydat, nameinte='histodim/')
               
        # temp -- restrict compfrac and other plots to indxmodlpntscomp
        plot_compfrac(gdat, gdatmodi, strg)
        
        for i in gdat.indxener:
            for m in gdat.indxevttplot:
                plot_scatpixl(gdat, gdatmodi, strg)
       
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and gdat.calcerrr:
            if strg != 'true':
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strg, 'errrcnts', strgcbar='resicnts', thisindxener=i)
        
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                plot_genemaps(gdat, gdatmodi, strg, 'datacnts', thisindxener=i, thisindxevtt=m)
                if numbpopl > 1:
                    for l in indxpopl:
                        plot_genemaps(gdat, gdatmodi, strg, 'datacnts', thisindxener=i, thisindxevtt=m, thisindxpopl=l)
        
        # temp
        if numbtrap > 0:
            for l in indxpopl:
                liststrgfeatmodu = getattr(gdat, strgmodl + 'liststrgfeatmodu')
                liststrgpdfnmodu = getattr(gdat, strgmodl + 'liststrgpdfnmodu')
                for strgfeat, strgpdfn in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]):
                    if strgpdfn == 'tmplreln':
                        plot_genemaps(gdat, gdatmodi, 'fitt', 'lpdfspatpriointp', tdim=True)
                    if strgpdfn == 'tmplgaum':
                        plot_genemaps(gdat, gdatmodi, strg, 'lpdfspatpriointp', tdim=True)
    
        if gdat.elemtype == 'lens':
            if strgmodl == 'true':
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strg, 's2nr', thisindxener=i)
                
            plot_genemaps(gdat, gdatmodi, strg, 'magn', tdim=True)
            plot_genemaps(gdat, gdatmodi, strg, 'conv', tdim=True)
            if numbtrap > 0:
                plot_genemaps(gdat, gdatmodi, strg, 'convelem', tdim=True)
            for i in gdat.indxener:
                plot_genemaps(gdat, gdatmodi, strg, 'hostcntsmaps', strgcbar='datacnts', thisindxener=i, tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'lenscnts', strgcbar='datacnts', thisindxener=i, tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'gradlenscntsmgtd', strgcbar='datacnts', thisindxener=i, tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'modlcntsuncv', strgcbar='datacnts', thisindxener=i, tdim=True)
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if strg != 'true':
                    plot_genemaps(gdat, gdatmodi, strg, 'llik', strgcbar='llikmaps', thisindxener=i, thisindxevtt=m)
                plot_genemaps(gdat, gdatmodi, strg, 'modlcnts', strgcbar='datacnts', thisindxener=i, thisindxevtt=m)
                plot_genemaps(gdat, gdatmodi, strg, 'resicnts', thisindxener=i, thisindxevtt=m)
            
        if gdat.elemtype == 'lens':
            plot_gene(gdat, gdatmodi, strg, 'convpsecodim', 'meanwvecodim', lablxaxi='$k$ [1/kpc]', lablyaxi='$P(k)$', limtydat=[1e-1, 1e2], scalxaxi='logt', scalyaxi='logt')
            plot_gene(gdat, gdatmodi, strg, 'histdefl', 'meandefl', scal='self', lablxaxi=r'$\alpha$ [arcsec]', lablyaxi=r'$N_{pix}$', factxdat=gdat.anglfact, \
                                                                                                                                            histodim=True)
            if numbtrap > 0:
                plot_gene(gdat, gdatmodi, strg, 'convpsecelemodim', 'meanwvecodim', lablxaxi='$k$ [1/kpc]', lablyaxi='$P_{sub}(k)$', \
                                                                                                                    limtydat=[1e-5, 1e-1], scalxaxi='logt', scalyaxi='logt')
                plot_gene(gdat, gdatmodi, strg, 'histdeflelem', 'meandeflelem', scal='self', lablxaxi=r'$\alpha$ [arcsec]', \
                                                                                                         lablyaxi=r'$N_{pix}$', factxdat=gdat.anglfact, histodim=True)
    
        if gdat.elemtype == 'lens':
    
            # gradient of the lens emission
            plot_defl(gdat, gdatmodi, strg, 'gradlenscnts')
            
            # overall deflection field
            plot_defl(gdat, gdatmodi, strg, multfact=0.1)
            
            # deflection field due to individual lenses
            for k in range(numbdeflsingplot):  
                if k == 0:
                    multfact = 0.1
                elif k == 1:
                    multfact = 1.
                elif k >= 2:
                    multfact = 10.
                plot_defl(gdat, gdatmodi, strg, indxdefl=k, multfact=multfact)
            
            # residual deflection field
            if strg != 'true' and gdat.datatype == 'mock':
                plot_defl(gdat, gdatmodi, strg, strgcomp='resi', multfact=100.)
                for k in range(numbdeflsingplot):
                    plot_defl(gdat, gdatmodi, strg, strgcomp='resi', indxdefl=k, multfact=100.)
                
                if numbtrap > 0:
                    plot_genemaps(gdat, gdatmodi, strg, 'resiconvelem', tdim=True)
                    plot_genemaps(gdat, gdatmodi, strg, 'resiconvelemperc', tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'resimagn', tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'resimagnperc', tdim=True)
    

def plot_post(gdat=None, pathpcat=None, verbtype=1, prio=False):
    
    if gdat.verbtype > 0:
        print 'Producing postprocessing plots...'

    timetotlinit = gdat.functime()
    
    gdat.strgbest = 'ML'
  
    # prior components
    if gdat.diagmode:
        if gdat.verbtype > 0:
            print 'Plotting the prior distribution...'

        gdat.indxlpri = arange(gdat.numblpri)
        gdat.indxlpau = arange(gdat.numblpau)
        indxbadd = where(isfinite(gdat.listlpri) == False)
        if indxbadd[0].size > 0:
            gdat.listlpri[indxbadd] = 5.
        indxbadd = where(isfinite(gdat.listlpriprop) == False)
        if indxbadd[0].size > 0:
            gdat.listlpriprop[indxbadd] = 5.
        
        for n in gdat.indxproptype:
            pathpostlpri = getattr(gdat, 'path' + gdat.namesampdist + 'finllpri')
            for k in gdat.indxlpri:
                if (gdat.listlpri[:, k] != 0.).any():
                    path = pathpostlpri + 'lpri%04d' % k + gdat.nameproptype[n]
                    tdpy.mcmc.plot_trac(path, gdat.listlpri[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
                if (gdat.listlpriprop[:, k] != 0.).any():
                    path = pathpostlpri + 'lpriprop%04d' % k + gdat.nameproptype[n]
                    tdpy.mcmc.plot_trac(path, gdat.listlpriprop[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
                if (gdat.listlpriprop[:, k] - gdat.listlpri[:, k] != 0.).any():
                    path = pathpostlpri + 'lpridelt%04d' % k + gdat.nameproptype[n]
                    tdpy.mcmc.plot_trac(path, gdat.listlpriprop[gdat.listindxsamptotl[n], k] - gdat.listlpri[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
            for k in gdat.indxlpau:
                if (gdat.listlpau[:, k] != 0.).any():
                    path = pathpostlpri + 'lpau%04d' % k + gdat.nameproptype[n]
                    tdpy.mcmc.plot_trac(path, gdat.listlpau[gdat.listindxsamptotl[n], k], '%04d' % k, logthist=True)
            if (gdat.listlfctasym != 0.).any():
                path = pathpostlpri + 'lfctprop' + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlfctasym[gdat.listindxsamptotl[n]], '$q_{lfct}$', logthist=True)
    
    # Gelman-Rubin test
    pathdiag = getattr(gdat, 'path' + gdat.namesampdist + 'finldiag')
    if gdat.numbproc > 1:
        if isfinite(gdat.gmrbstat).all():
            if gdat.verbtype > 0:
                print 'Gelman-Rubin TS...'
    
            figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
            minm = min(amin(gdat.gmrbstat), amin(gdat.gmrbfixp))
            maxm = max(amax(gdat.gmrbstat), amax(gdat.gmrbfixp))
            bins = linspace(minm, maxm, 40)
            # temp
            try:
                axis.hist(gdat.gmrbstat.flatten(), bins=bins, label='Data proj.')
            except:
                print 'gdat.gmrbstat'
                print gdat.gmrbstat
                print 'bins'
                print bins
            try:
                axis.hist(gdat.gmrbfixp, bins=bins, label='Fixed dim.')
            except:
                print 'gdat.gmrbfixp'
                print gdat.gmrbfixp
                print 'bins'
                print bins
            axis.set_xlabel('PSRF')
            axis.set_ylabel('$N_{stat}$')
            plt.tight_layout()
            figr.savefig(pathdiag + 'gmrbhist.pdf')
            plt.close(figr)
            
            try:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        maps = gdat.gmrbstat[i, :, m]
                        path = pathdiag + 'gmrbmaps_%d%d.pdf' % (i, m)
                        tdpy.util.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                                minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                                minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
            except:
                print 'gdat.gmrbstat'
                print gdat.gmrbstat
        else:
            print 'Inappropriate Gelman-Rubin test statistics encountered.'

    # plot autocorrelation
    if gdat.verbtype > 0:
        print 'Autocorrelation...'
    tdpy.mcmc.plot_atcr(pathdiag, gdat.atcr, gdat.timeatcr)
    
    # plot proposal efficiency
    if gdat.verbtype > 0:
        print 'Acceptance ratio...'
    numbtimemcmc = 20
    binstimemcmc = linspace(0., gdat.numbswep, numbtimemcmc)
    numbtick = 2
    numbplotfram = 8
    cntr = 0
    while cntr < gdat.numbproptype:
        thisnumbplot = min(numbplotfram, gdat.numbproptype - cntr)
        sizefigryaxi = max(thisnumbplot * gdat.plotsize / 4., gdat.plotsize / 2.)
        figr, axgr = plt.subplots(thisnumbplot, 1, figsize=(gdat.plotsize, sizefigryaxi), sharex='all')
        if thisnumbplot == 1:
            axgr = [axgr]
        for n, axis in enumerate(axgr):
            histtotl = axis.hist(gdat.listindxsamptotl[n+cntr], bins=binstimemcmc)[0]
            histaccp = axis.hist(gdat.listindxsamptotlaccp[n+cntr], bins=binstimemcmc)[0]
            # temp
            axis.set_ylabel('%s' % gdat.lablproptype[n+cntr])
            if n + cntr == gdat.numbproptype - 1:
                axis.set_xlabel('$i_{samp}$')
            maxm = amax(histtotl)
            axis.set_ylim([0., maxm])
            if sum(histtotl) != 0:
                axis.set_title('%.3g%%' % int(sum(histaccp) / sum(histtotl) * 100.))
            listtick = linspace(maxm / 2., maxm, numbtick)
            listlabltick = ['%.3g' % tick for tick in listtick]
            axis.set_yticks(listtick)
            axis.set_yticklabels(listlabltick)
           
            if False:
                print 'gdat.listindxsamptotl[n]'
                print gdat.listindxsamptotl[n]
                print 'gdat.listindxsamptotlaccp[n]'
                print gdat.listindxsamptotlaccp[n]
                print 'binstimemcmc'
                print binstimemcmc
                print 'gdat.lablproptype[n]'
                print gdat.lablproptype[n]
                print 'listtick'
                print listtick
                print 'listlabltick'
                print listlabltick
                print
    
        plt.tight_layout()
        figr.savefig(pathdiag + 'accpratiproptype%04d.pdf' % (cntr / numbplotfram))
        plt.close(figr)
        cntr += numbplotfram
   
    # post-processing frame plots
    if gdat.numbframpost != None:
        gdatmodi = tdpy.util.gdatstrt()
        gdat.indxsamptotlfram = arange(gdat.numbframpost) * (gdat.indxsamptotl - gdat.indxsamptotl % gdat.numbframpost) / gdat.numbframpost
        for n in gdat.indxsamptotlfram:
            gdatmodi.cntrswep = n
            gdatmodi.thisindxelemfull = deepcopy(gdat.listindxelemfull[n])
            gdatmodi.thissampvarb = copy(gdat.listsampvarb[n, :])
            proc_samp(gdat, gdatmodi, 'this')
            plot_samp(gdat, gdatmodi, 'this')

    # plot split and merge diagnostics
    if gdat.fittnumbtrap > 0 and gdat.probtran > 0. and gdat.probbrde < 1.:
        if gdat.verbtype > 0:
            print 'Split and merge related plots...'
    
        indxsampsplttotl = where(gdat.listindxproptype == gdat.indxproptypesplt)[0]
        indxsampsplt = intersect1d(where(gdat.listindxproptype == gdat.indxproptypesplt)[0], where(gdat.listaccpprop)[0])
        indxsampmergtotl = where(gdat.listindxproptype == gdat.indxproptypemerg)[0]
        indxsampmerg = intersect1d(where(gdat.listindxproptype == gdat.indxproptypemerg)[0], where(gdat.listaccpprop)[0])
        indxsampspmrtotl = concatenate((indxsampsplttotl, indxsampmergtotl))
        if indxsampspmrtotl.size > 0:
    
            ## labels and names
            listlabl = ['$u_{frac}$', '$u_l$', r'$u_b$', '$N_{pair}$', r'$\log\alpha_j$', r'$\log\alpha_c$']
            listname = ['fracauxi', 'lgalauxi', 'bgalauxi', 'numbpair', 'ljcbfact', 'lcomfact']
            
            ## variables
            listvarb = [gdat.listauxipara[:, 0], gdat.anglfact * gdat.listauxipara[:, 1], gdat.anglfact * gdat.listauxipara[:, 2], gdat.listnumbpair, \
                                                                                                                                  gdat.listljcbfact, gdat.listlcomfact]
           
            for k in range(len(listlabl)):
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                maxm = amax(listvarb[k][indxsampspmrtotl])
                minm = amin(listvarb[k][indxsampspmrtotl])
                bins = linspace(minm, maxm, 40)
              
                axis.hist(listvarb[k][indxsampsplttotl], bins=bins, label='Split', alpha=gdat.alphmrkr, color='c')
                axis.hist(listvarb[k][indxsampmergtotl], bins=bins, label='Merge', alpha=gdat.alphmrkr, color='y')
                axis.hist(listvarb[k][indxsampsplt], bins=bins, label='Split/Accepted', alpha=gdat.alphmrkr, edgecolor='c', lw=20, facecolor='none')
                axis.hist(listvarb[k][indxsampmerg], bins=bins, label='Merge/Accepted', alpha=gdat.alphmrkr, edgecolor='y', lw=20, facecolor='none')
                
                axis.set_ylabel('$N_{samp}$')
                axis.set_xlabel(listlabl[k])
    
                make_legd(axis)
                plt.tight_layout()
                path = getattr(gdat, 'path' + gdat.namesampdist + 'finlspmr') + listname[k] + '.pdf'
                figr.savefig(path)
                plt.close(figr)
    
    if gdat.verbtype > 0:
        print 'Proposal execution times...'
    plot_chro(gdat)

    if gdat.verbtype > 0:
        print 'Derived quantities...'

    # posterior versions of the frame plots
    plot_samp(gdat, None, 'post')
    proc_samp(gdat, None, 'mlik')

    if gdat.verbtype > 0:
        print 'A mosaic of samples...'
    
    ## mosaic of images of posterior catalogs 
    plot_mosa(gdat)
    
    if gdat.verbtype > 0:
        print 'Fixed dimensional parameter traces...'
    
    ## randomly selected trandimensional parameters
    if gdat.fittnumbtrap > 0:
        if gdat.verbtype > 0:
            print 'Transdimensional parameters...'
    
        # choose the parameters based on persistence
        stdvlistsamptran = std(gdat.listsamp[:, gdat.fittindxsamptrap], axis=0)
        indxtrapgood = where(stdvlistsamptran > 0.)[0]
        numbtrapgood = indxtrapgood.size
        numbtrapplot = min(10, numbtrapgood)
        indxtrapplot = sort(choice(gdat.fittindxsamptrap[indxtrapgood], size=numbtrapplot, replace=False))

        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'listelemfrst'
        tdpy.mcmc.plot_grid(path, gdat.listsampvarb[:, gdat.fittindxsamptrap[:3]] * gdat.fittfactplotpara[gdat.fittindxsamptrap[:3]], \
                                                                                            [gdat.fittlablpara[k] for k in gdat.fittindxsamptrap[:3]])

        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'listsamp'
        tdpy.mcmc.plot_grid(path, gdat.listsamp[:, indxtrapplot], ['%d' % k for k in indxtrapplot])
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'listsampvarb'
        tdpy.mcmc.plot_grid(path, gdat.listsampvarb[:, indxtrapplot] * gdat.fittfactplotpara[indxtrapplot], [gdat.fittlablpara[k] for k in indxtrapplot])

    ## scalar variables
    ### trace and marginal distribution of each parameter
    for name in gdat.listnamevarbscal:
        scal = getattr(gdat, 'scal' + name) 
        factplot = getattr(gdat, 'fact' + name + 'plot')
        corr = getattr(gdat, 'corr' + name)
        if corr == None:
            truepara = None
        else:
            truepara = getattr(gdat, 'corr' + name) * factplot
        labltotl = getattr(gdat, 'labl' + name + 'totl')
        listvarb = getattr(gdat, 'list' + name) * factplot
        mlik = getattr(gdat, 'mlik' + name) * factplot
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscaltrac') + name
        tdpy.mcmc.plot_trac(path, listvarb, labltotl, truepara=truepara, scalpara=scal, varbdraw=[mlik], labldraw=[''], colrdraw=['r'])
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalhist') + name
        tdpy.mcmc.plot_hist(path, listvarb, labltotl, truepara=truepara, scalpara=scal, varbdraw=[mlik], labldraw=[''], colrdraw=['r'])
       
        for nameseco in gdat.listnamevarbscal:
            if name == nameseco:
                continue
            pathjoin = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscaljoin') + name + nameseco
            scalseco = getattr(gdat, 'scal' + nameseco) 
            factplotseco = getattr(gdat, 'fact' + nameseco + 'plot')
            corrseco = getattr(gdat, 'corr' + nameseco)
            if corrseco == None:
                trueparaseco = None
            else:
                trueparaseco = getattr(gdat, 'corr' + nameseco) * factplotseco
            labltotlseco = getattr(gdat, 'labl' + nameseco + 'totl')
            listvarbseco = getattr(gdat, 'list' + nameseco) * factplotseco
            mlikseco = getattr(gdat, 'mlik' + nameseco) * factplotseco
            
            listjoin = vstack((listvarb, listvarbseco)).T
            tdpy.mcmc.plot_grid(pathjoin, listjoin, [labltotl, labltotlseco], scalpara=[scal, scalseco], truepara=[truepara, trueparaseco], \
                                                                                                                                  join=True, varbdraw=[mlik, mlikseco])
        
    if gdat.checprio and not prio:
        # this works only for scalar variables -- needs to be generalized to all variables
        for namevarbscal in gdat.listnamevarbscal:
            titl = '$D_{KL} = %.3g$' % getattr(gdat, 'infototl' + namevarbscal)
            xdat = getattr(gdat, 'mean' + namevarbscal) * getattr(gdat, 'fact' + namevarbscal + 'plot')
            lablxdat = getattr(gdat, 'labl' + namevarbscal + 'totl')
            
            path = gdat.pathinfo + 'info' + namevarbscal
            ydat = getattr(gdat, 'info' + namevarbscal)
            lablydat = r'$\Delta D_{KL}$'
            tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, titl=titl)

            path = gdat.pathinfo + 'pdfn' + namevarbscal
            ydat = [getattr(gdat, 'pdfnpost' + namevarbscal), getattr(gdat, 'pdfnprio' + namevarbscal)]
            lablydat = r'$P$'
            legd = ['$P$(%s|$D$)' % lablxdat, '$P$(%s)' % lablxdat]
            tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, titl=titl, colr=['k', 'k'], linestyl=['-', '--'], legd=legd)
    
    if gdat.verbtype > 0:
        print 'Fixed dimensional parameter covariance...'
    
    ### covariance
    ## overall
    path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'fixp'
    truepara = gdat.fittcorrfixp * gdat.fittfactfixpplot
    mlikpara = gdat.mlikfixp * gdat.fittfactfixpplot
    tdpy.mcmc.plot_grid(path, gdat.listfixp * gdat.fittfactfixpplot[None, :], gdat.fittlablfixptotl, truepara=truepara, varbdraw=mlikpara)
    
    ## individual processes
    if gdat.numbproc > 1:
        for k in gdat.indxproc:
            path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalproc') + 'proc%04d' % k
            tdpy.mcmc.plot_grid(path, gdat.listsampvarbproc[:, k, gdat.fittindxfixp] * gdat.fittfactfixpplot[None, gdat.fittindxfixp], \
                                gdat.fittlablfixptotl[gdat.fittindxfixp], truepara=gdat.fittcorrfixp[gdat.fittindxfixp] * gdat.fittfactfixpplot[gdat.fittindxfixp])
    
    ## grouped covariance plots
    if gdat.verbtype > 0:
        print 'Hyperparameters...'
    
    if gdat.fittnumbtrap > 0:
        #### hyperparameters
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'hypr'
        tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.fittindxfixphypr] * gdat.fittfactfixpplot[None, gdat.fittindxfixphypr], gdat.fittlablfixptotl[gdat.fittindxfixphypr], \
                                                                                           truepara=[gdat.fittcorrfixp[k] * gdat.fittfactfixpplot[k] for k in gdat.fittindxfixphypr])
    
    if gdat.verbtype > 0:
        print 'PSF parameters...'
    
    #### PSF
    if gdat.proppsfp:
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'psfp'
        tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.fittindxfixppsfp] * gdat.fittfactfixpplot[None, gdat.fittindxfixppsfp], gdat.fittlablfixptotl[gdat.fittindxfixppsfp], \
                                          truepara=[gdat.fittcorrfixp[k] * gdat.fittfactfixpplot[k] for k in gdat.fittindxfixppsfp], numbplotside=gdat.fittnumbpsfptotl)
    if gdat.verbtype > 0:
        print 'Background parameters...'
    
    #### backgrounds
    if gdat.propbacp:
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'bacp'
        tdpy.mcmc.plot_grid(path, gdat.listfixp[:, gdat.fittindxfixpbacp], gdat.fittlablfixptotl[gdat.fittindxfixpbacp], \
                                                                                                        truepara=[gdat.fittcorrfixp[k] for k in gdat.fittindxfixpbacp])
        if gdat.fittnumbback == 2 and gdat.fittspecback == [None, None]:
            for i in gdat.indxener:
                indx = gdat.fittindxfixpbacp[gdat.fittindxback*gdat.numbener+i]
                path = getattr(gdat, 'path' + gdat.namesampdist + 'finlvarbscalcova') + 'bacpene%d' % i
                tdpy.mcmc.plot_grid(path, gdat.listfixp[:, indx], gdat.fittlablfixptotl[indx], truepara=[gdat.fittcorrfixp[k] for k in indx], join=True)
        
    if gdat.verbtype > 0:
        print 'Binned transdimensional parameters...'
   
    # stacked posteiors binned in position and flux
    if gdat.fittnumbtrap > 0:
        liststrgbins = ['quad', 'full']
        for l in gdat.fittindxpopl:
            plot_postbindmaps(gdat, l, 'cumu')
            for strgbins in liststrgbins:
                for strgfeatsign in gdat.fittliststrgfeatsign[l]:
                    plot_postbindmaps(gdat, l, strgbins, strgfeatsign)

    if gdat.verbtype > 0:
        print 'Prior and likelihood...'
    
    ## log-prior and log-likelihood
    gdat.listlliktotl = sum(gdat.listllik, axis=(1, 2, 3))
    gdat.listlpritotl = sum(gdat.listlpri, axis=1)

    for strg in ['lpri', 'llik']:

        if strg == 'lpri':
            labltemp = '\ln P(M)'
        if strg == 'llik':
            labltemp = '\ln P(D|M)'
        labl = r'$%s$' % labltemp

        setattr(gdat, 'list' + strg + 'flat', getattr(gdat, 'list' + strg + 'totl').flatten())
        
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finl') + strg
        titl = r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi)
        print 'strg'
        print strg
        print getattr(gdat, 'list' + strg + 'flat')
        print
        tdpy.mcmc.plot_hist(path, getattr(gdat, 'list' + strg + 'flat'), labl, titl)
        varbdraw = []
        labldraw = []
        colrdraw = []
        if gdat.datatype == 'mock':
            varbdraw += [getattr(gdat, 'true' + strg + 'totl')]
            labldraw += ['True model']
            colrdraw += ['g']
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'list' + strg + 'flat'), labl, varbdraw=varbdraw, labldraw=labldraw, colrdraw=colrdraw)
    
        if strg == 'llik':
            labldelt = r'$\Delta%s$' % labltemp
            setattr(gdat, 'listdelt' + strg + 'totlflat', getattr(gdat, 'listdelt' + strg + 'totl').flatten())
            pathbase = getattr(gdat, 'path' + gdat.namesampdist + 'finldelt%s' % strg)
            path = pathbase + 'delt%s' % strg
    
            tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat'), labldelt)
            
            if gdat.numbproc > 1:
                for k in gdat.indxproc:
                    path = pathbase + 'delt%s_proc%04d' % (strg, k)
                    tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totl')[:, k], labldelt, titl='Process %d' % k)
            for n in gdat.indxproptype:
                path = pathbase + 'delt%s_%s' % (strg, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotl[n]], labldelt, titl=gdat.nameproptype[n])
                path = pathbase + 'delt%s_%s_accp' % (strg, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotlaccp[n]], labldelt, titl=gdat.nameproptype[n] + ', Accepted')
                path = pathbase + 'delt%s_%s_reje' % (strg, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strg + 'totlflat')[gdat.listindxsamptotlreje[n]], labldelt, titl=gdat.nameproptype[n] + ', Rejected')
        
    # plot resident memory
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxswep, mean(gdat.listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(pathdiag + 'memoresi.pdf')
    plt.close(figr)

    # animate the frame plots
    if gdat.makeanim:
        make_anim(gdat.strgcnfg, full=True)

    timetotlfinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Plots and animations are produced in %.3g seconds.' % (timetotlfinl - timetotlinit)


def plot_chro(gdat):
    
    pathdiag = getattr(gdat, 'path' + gdat.namesampdist + 'finldiag')

    gdat.listchro *= 1e3
    indxchro = array([0, 1, 2, 4])
    binstime = logspace(log10(amin(gdat.listchro[where(gdat.listchro > 0)])), log10(amax(gdat.listchro[:, indxchro])), 50)

    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    for k in range(gdat.numbchro):
        varb = gdat.listchro[:, k]
        axis.hist(varb, binstime, log=True, label=gdat.listlegdchro[k], linewidth=2, alpha=0.3)

    axis.set_title(r'$\langle t \rangle$ = %.3g ms' % mean(gdat.listchro[where(gdat.listchro[:, 0] > 0)[0], 0]))
    axis.set_xlim([amin(binstime), amax(binstime)])
    axis.set_xscale('log')
    axis.set_ylim([0.5, None])
    make_legd(axis)
    axis.set_xlabel('$t$ [ms]')
    
    plt.tight_layout()
    figr.savefig(pathdiag + 'chro.pdf')
    plt.close(figr)

    gdat.listchro *= 1e3
   
    if (gdat.listchro != 0).any() and False:
        figr, axcl = plt.subplots(gdat.numbchro, 1, figsize=(2 * gdat.plotsize, gdat.plotsize * gdat.numbchro / 3.))
        maxmchro = amax(gdat.listchro)
        minmchro = amin(gdat.listchro[where(gdat.listchro > 0)])
        binstime = logspace(log10(minmchro), log10(maxmchro), 50)
    
        for k in range(gdat.numbchro):
            chro = gdat.listchro[where(gdat.listchro[:, k] > 0)[0], k]
            axcl[k].hist(chro, binstime, log=True, label=listlabl[k])
            axcl[k].set_xlim([amin(binstime), amax(binstime)])
            axcl[k].set_ylim([0.5, None])
            axcl[k].set_ylabel(listlabl[k])
            axcl[k].set_xscale('log')
            if k != gdat.numbchro - 1:
                axcl[k].set_xticklabels([])
            axcl[k].axvline(mean(chro), ls='--', alpha=0.2, color='black')
        axcl[-1].set_xlabel('$t$ [ms]')
        plt.subplots_adjust(hspace=0.05)
        figr.savefig(pathdiag + 'chro.pdf')
        plt.close(figr)


def plot_compfrac(gdat, gdatmodi, strg):
    
    if strg == 'true':
        strgmodl = strg
    else:
        strgmodl = 'fitt'
    
    specback = getattr(gdat, strgmodl + 'specback')
    indxback = getattr(gdat, strgmodl + 'indxback')
    maxmgang = getattr(gdat, strgmodl + 'maxmgang')
    indxfixpbacp = getattr(gdat, strgmodl + 'indxfixpbacp')
    backfluxmean = getattr(gdat, strgmodl + 'backfluxmean')
    indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
    numblablcompfrac = getattr(gdat, strgmodl + 'numblablcompfrac')
    numblablcompfracspec = getattr(gdat, strgmodl + 'numblablcompfracspec')

    listydat = zeros((numblablcompfracspec, gdat.numbener))
    listyerr = zeros((2, numblablcompfracspec, gdat.numbener))
   
    cntr = 0
        
    ## background templates
    for c in indxback:
        temp = retr_varb(gdat, gdatmodi, strg, 'sampvarb', indx=[indxfixpbacp[indxbacpback[c]]])
        if specback[c] != None:
            norm = temp * specback[c]
        else:
            norm = temp
        listydat[cntr, :] = norm * backfluxmean[c, :]

        if strg == 'post':
            temp = retr_varb(gdat, gdatmodi, strg, 'sampvarb', indx=[indxfixpbacp[indxbacpback[c]]], perc='errr')
            if specback[c] != None:
                norm = temp * specback[c]
            else:
                norm = temp
            listyerr[:, cntr, :] = norm * backfluxmean[None, c, :]

        if gdat.elemtype == 'lens':
            listydat[cntr, :] *= 4. * maxmgang**2
        
        cntr += 1
    
    ## PS
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        listydat[cntr, :] = retr_varb(gdat, gdatmodi, strg, 'pntsfluxmean')
        if strg == 'post':
            listyerr[:, cntr, :] = retr_varb(gdat, gdatmodi, strg, 'pntsfluxmean', perc='errr')
        cntr += 1

    if gdat.elemtype == 'lens':
        for strgtemp in ['sour', 'host']:
            indxvarb = getattr(gdat, strgmodl + 'indxfixpspec' + strgtemp)
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
    if numblablcompfrac > 1:
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
        
        listlablcompfracspec = getattr(gdat, strgmodl + 'listlablcompfracspec')

        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        
        listmrkr = [(2 + k/2, 1 + k % 2, 0) for k in range(16)]

        # plot reference spectra
        if gdat.listspecrefrplot != None:
            for k in range(len(gdat.listspecrefrplot)):
                axis.plot(gdat.listenerrefrplot[k], gdat.listspecrefrplot[k], label=gdat.listlablrefrplot[k], color='m')

        xdat = gdat.meanener
        for k in range(numblablcompfracspec):
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
            axis.errorbar(xdat, ydat, yerr=yerr, label=listlablcompfracspec[k], color=colr, marker=mrkr, ls=linestyl, markersize=15)
    
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
        axis.scatter(gdat.truefluxbrgt, gdat.truefluxbrgtassc, alpha=gdat.alphmrkr, color='g', label=gdat.legdtrue)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_ylim([gdat.minmfluxplot, gdat.maxmfluxplot])
    axis.set_xlabel(r'$%s_{max}$%s' % (gdat.lablflux, gdat.lablfluxunit))
    axis.set_ylabel(r'$%s_{asc}$%s' % (gdat.lablflux, gdat.lablfluxunit))
    make_legd(axis, loca=2)
    
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatbrgt')
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def savefigr(gdat, gdatmodi, figr, path):
    
    #if gdatmodi != None and gdat.numbproc > 1:
    #    gdatmodi.lock.acquire()
    #    print 'Process %d acquiring the lock...' % gdatmodi.indxprocwork 
    
    plt.savefig(path)
    
    #if gdatmodi != None and gdat.numbproc > 1:
    #    gdatmodi.lock.release()
    #    print 'Process %d releasing the lock...' % gdatmodi.indxprocwork 
        

def plot_elemtdim(gdat, gdatmodi, strg, l, strgplottype, strgfrst, strgseco, strgmome='medi'):
    
    sizelarg = 10
    sizesmll = 1
    
    legdmome = getattr(gdat, 'legd' + strgmome)

    if strg == 'true':  
        strgmodl = strg
    else:
        strgmodl = 'fitt'

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if strg == 'post':
        labl = gdat.legdsampdist + ' ' + legdmome
        if strgplottype == 'hist':
            varb = getattr(gdat, strgmome + 'hist' + strgfrst + strgseco)[l, :, :]
            varbfrst = getattr(gdat, 'bins' + strgfrst) * getattr(gdat, 'fact' + strgfrst + 'plot')
            varbseco = getattr(gdat, 'bins' + strgseco) * getattr(gdat, 'fact' + strgseco + 'plot')
            imag = axis.pcolor(varbfrst, varbseco, varb.T, cmap='Greys', label=labl)
            plt.colorbar(imag)
        else:
            if gdat.condcatl:
                varbfrst = zeros(gdat.numbprvlhigh)
                varbseco = zeros(gdat.numbprvlhigh)
                cntr = 0
                for r in gdat.indxstkscond:
                    if r in gdat.indxprvlhigh:
                        varbfrst[cntr] = gdat.dictglob['poststkscond'][r][strgfrst][l] * getattr(gdat, 'fact' + strgfrst + 'plot')
                        varbseco[cntr] = gdat.dictglob['poststkscond'][r][strgseco][l] * getattr(gdat, 'fact' + strgseco + 'plot')
                        cntr += 1
                axis.scatter(varbfrst, varbseco, alpha=gdat.alphmrkr, color='black', label=gdat.legdsamp)

    if strg == 'this' or strg == 'mlik':
        if strgplottype == 'hist':
            meanfrst = getattr(gdat, 'bins' + strgfrst) * getattr(gdat, 'fact' + strgfrst + 'plot')
            meanseco = getattr(gdat, 'bins' + strgseco) * getattr(gdat, 'fact' + strgseco + 'plot')
            hist = getattr(gdatmodi, strg + 'hist' + strgfrst + strgseco)[l, :, :]
            imag = axis.pcolor(meanfrst, meanseco, hist.T, cmap='Blues', label=gdat.legdsamp, alpha=gdat.alphmrkr)
        else:
            varbfrst = getattr(gdatmodi, 'this' + strgfrst)[l] * getattr(gdat, 'fact' + strgfrst + 'plot')
            varbseco = getattr(gdatmodi, 'this' + strgseco)[l] * getattr(gdat, 'fact' + strgseco + 'plot')
            axis.scatter(varbfrst, varbseco, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
    
    # true
    if gdat.trueinfo:
        try:
            truevarbfrst = getattr(gdat, 'true' + strgfrst)[l] * getattr(gdat, 'fact' + strgfrst + 'plot')
            truevarbseco = getattr(gdat, 'true' + strgseco)[l] * getattr(gdat, 'fact' + strgseco + 'plot')
            axis.scatter(truevarbfrst, truevarbseco, alpha=gdat.alphmrkr, color='g', label=gdat.legdtrue, s=sizelarg)
        except:
            pass

    # experimental
    if gdat.datatype == 'mock':
        try:
            exprvarbfrsttotl = getattr(gdat, 'expr' + strgfrst + 'totl')[l] * getattr(gdat, 'fact' + strgfrst + 'plot')
            exprvarbsecototl = getattr(gdat, 'expr' + strgseco + 'totl')[l] * getattr(gdat, 'fact' + strgseco + 'plot')
            axis.scatter(exprvarbfrsttotl, exprvarbsecototl, alpha=0.05, color='r', label=gdat.nameexpr + ' All', s=sizesmll)
        except:
            pass
        try:
            exprvarbfrst = getattr(gdat, 'expr' + strgfrst)[l] * getattr(gdat, 'fact' + strgfrst + 'plot')
            exprvarbseco = getattr(gdat, 'expr' + strgseco)[l] * getattr(gdat, 'fact' + strgseco + 'plot')
            axis.scatter(exprvarbfrst, exprvarbseco, alpha=0.3, color='r', label=gdat.nameexpr, s=sizelarg)
        except:
            pass

    plot_sigmcont(gdat, axis, l, strgfrst=strgfrst, strgseco=strgseco)
    
    scalfrst = getattr(gdat, 'scal' + strgfrst + 'plot')
    scalseco = getattr(gdat, 'scal' + strgseco + 'plot')
    if scalfrst == 'logt':
        axis.set_xscale('log')
    if scalseco == 'logt':
        axis.set_yscale('log')

    axis.set_xlabel(getattr(gdat, 'labl' + strgfrst + 'totl'))
    axis.set_ylabel(getattr(gdat, 'labl' + strgseco + 'totl'))
    limtfrst = getattr(gdat, 'limt' + strgfrst + 'plot')
    limtseco = getattr(gdat, 'limt' + strgseco + 'plot')
    factplotfrst = getattr(gdat, 'fact' + strgfrst + 'plot')
    factplotseco = getattr(gdat, 'fact' + strgseco + 'plot')
    axis.set_xlim(limtfrst * factplotfrst)
    axis.set_ylim(limtseco * factplotseco)

    plt.tight_layout()
    if strg == 'post':
        strgmometemp = strgmome
    else:
        strgmometemp = ''
    path = retr_plotpath(gdat, gdatmodi, strg, '%s%s%s%spop%d' % (strgmometemp, strgplottype, strgfrst, strgseco, l), nameinte=strgplottype + 'tdim/')
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_sigmcont(gdat, axis, l, strgfrst=None, strgseco=None):
    
    if strgfrst == 'deltllik' or strgseco == 'deltllik':
        for pval in gdat.pvalcont:
            deltlliksigm = scipy.stats.chi2.ppf(1. - pval, gdat.fittnumbcomp[l])
            if strgfrst == 'deltllik':
                axis.axvline(deltlliksigm, ls='--', color='black', alpha=0.2) 
            if strgseco == 'deltllik':
                axis.axhline(deltlliksigm, ls='--', color='black', alpha=0.2) 
    

def plot_gene(gdat, gdatmodi, strg, strgydat, strgxdat, indxydat=None, strgindxydat=None, indxxdat=None, strgindxxdat=None, plottype='none', \
                     scal=None, scalxaxi=None, scalyaxi=None, limtxdat=None, limtydat=None, omittrue=False, nameinte='', \
                     lablxaxi='', lablyaxi='', factxdat=1., factydat=1., histodim=False, offslegd=None, tdim=False):
   
    if strg == 'true':
        strgmodl = strg
    else:
        strgmodl = 'fitt'

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
        # temp
        if strgxdat[4:] in gdat.fittliststrgfeattotl:
            xdat = getattr(gdat, strgxdat) * factxdat
        else:
            xdat = getattr(gdat, strgxdat) * factxdat
        ydat = retr_fromgdat(gdat, gdatmodi, strg, strgydat) * factydat
    
    if indxxdat != None:
        xdat = xdat[indxxdat]
    if indxydat != None:
        ydat = ydat[indxydat]
    
    if tdim:
        axis.scatter(xdat, ydat, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
    else:
        if histodim:
            # temp
            if strgxdat[4:] in gdat.fittliststrgfeattotl:
                deltxdat = getattr(gdat, 'delt' + strgxdat[4:]) * factxdat
                binsxdat = getattr(gdat, 'bins' + strgxdat[4:]) * factxdat
            else:
                deltxdat = getattr(gdat, 'delt' + strgxdat[4:]) * factxdat
                binsxdat = getattr(gdat, 'bins' + strgxdat[4:]) * factxdat

            xdattemp = binsxdat[:-1] + deltxdat / 2.
    
    if strg == 'post':
        yerr = retr_fromgdat(gdat, gdatmodi, strg, strgydat, mometype='errr') * factydat
        
        if indxydat != None:
            yerr = yerr[[slice(None)] + indxydat]
        
        # label
        if strgydat.startswith('hist'):
            ##  element distribution
            labl = gdat.legdsampdist + ' distribution'
        else:
            ##  other
            labl = gdat.legdsampdist + ' distribution'
        
        temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color='black', label=labl, lw=1, capsize=10)
        for caps in listcaps:
            caps.set_markeredgewidth(1)

    elif strg == 'this' or strg == 'mlik':
        
        if strg == 'this':
            legd = gdat.legdsamp
        else:
            legd = gdat.legdmlik

        if histodim:
            axis.bar(xdattemp, ydat, deltxdat, label=gdat.legdsamp, alpha=0.5)
        else:
            if plottype == 'errr':
                yerr = retr_fromgdat(gdat, gdatmodi, strg, strgydat, mometype='errr') * factydat

                if indxydat != None:
                    yerr = yerr[[slice(None)] + indxydat]
                temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color='b', label=legd, lw=1, capsize=10)
                for caps in listcaps:
                    caps.set_markeredgewidth(1)
            else:
                axis.plot(xdat, ydat, label=gdat.legdsamp, alpha=0.5)
    
    # reference
    try:
        if strgydat.startswith('hist'):
            strgtemp = strgydat[:4]

            bins = copy(getattr(gdat, 'bins' + strgtemp))
            varb = copy(getattr(gdat, 'expr' + strgtemp))
            if strgtemp == 'lgal' or strgtemp == 'bgal':
            	bins *= gdat.anglfact
            	varb *= gdat.anglfact
            
            if gdat.strgcnfg == 'pcat_chan_inpt':
                print 'plot_gene'
                print 'strgxdat'
                print strgxdat
                print 'strgydat'
                print strgydat
                print 'strgtemp'
                print strgtemp
                print 'bins'
                print bins
                print 'varb'
                print varb
                print
            axis.hist(varb, bins, label=gdat.strgcatl, alpha=0.5, color='r')
    except:
        pass
    
    # true
    if gdat.datatype == 'mock' and not omittrue:
        try:
            ydat = getattr(gdat, 'true' + strgydat)
            if indxydat != None:
                ydat = ydat[indxydat]
            if histodim:
                axis.bar(xdattemp, ydat, deltxdat, color='g', label=gdat.legdtrue, alpha=0.5)
            else:
                axis.plot(xdat, ydat, color='g', label=gdat.legdtrue, alpha=0.5)
        except:
            if gdat.verbtype > 0:
                print 'Skipping truth plot for %s...' % strgydat
    
    # axis scales
    if scalxaxi == 'logt':
        axis.set_xscale('log')
    if scalyaxi == 'logt':
        if where(ydat > 0.)[0].size > 0:
            axis.set_yscale('log')
    
    # axis labels
    axis.set_xlabel(lablxaxi)
    axis.set_ylabel(lablyaxi)

    # superimpose prior on the feature
    liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
    liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
    
    if strgydat.startswith('hist') and indxydat != None and strgydat[4:] in liststrgfeatprio[indxydat[0]]:
        xdatprio = getattr(gdat, strgmodl + strgxdat + 'prio') * factxdat
        if strg == 'true' or strg == 'post':
            gdattemp = gdat
        else:
            gdattemp = gdatmodi
    
        if gdat.datatype == 'mock' and not omittrue:
            truexdatprio = getattr(gdat, 'true' + strgxdat + 'prio') * factxdat
            trueydatsupr = getattr(gdat, 'true' + strgydat + 'prio')[indxydat]
            axis.plot(truexdatprio, trueydatsupr, ls='-', alpha=gdat.alphmrkr, color='g')
        
        if strg != 'true':
            if strg == 'post':
                ydatsupr = getattr(gdattemp, 'medi' + strgydat + 'prio')[indxydat]
                yerrsupr = getattr(gdattemp, 'errr' + strgydat + 'prio')[[slice(None)] + indxydat]
                labl = gdat.legdsampdist + ' hyper-distribution'
                tdpy.util.plot_braz(axis, xdatprio, ydatsupr, yerr=yerrsupr, lcol='lightgrey', dcol='grey', labl=labl)
            else:
                ydatsupr = getattr(gdattemp, strg + strgydat + 'prio')[indxydat]
                axis.plot(xdatprio, ydatsupr, ls='--', alpha=gdat.alphmrkr, color='b')

    if strgydat.startswith('hist') and indxydat != None and strgydat[4:] == 'deltllik':
        plot_sigmcont(gdat, axis, indxydat[0], strgfrst=strgxdat[4:])
    
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
        axis.set_ylim([amin(ydat), amax(ydat)])
    
    make_legd(axis, offs=offslegd)
    
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, strgydat, nameinte=nameinte)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatassc(gdat, gdatmodi, strg, l, strgfeat, plotdiff=False):
    
    figr, axis = plt.subplots(1, 1, figsize=(gdat.plotsize * gdat.numbener, gdat.plotsize))

    # prepare data to be plotted
    xdat = copy(getattr(gdat, 'true' + strgfeat)[l][0, :])
    xerr = tdpy.util.retr_errrvarb(getattr(gdat, 'true' + strgfeat)[l])
    
    minmplot = getattr(gdat, 'minm' + strgfeat)
    maxmplot = getattr(gdat, 'maxm' + strgfeat)
    binsplot = getattr(gdat, 'bins' + strgfeat)
    
    yerr = zeros((2, xdat.size))
    if gdatmodi == None:
        ydat = copy(getattr(gdat, 'medi' + strgfeat + 'assc')[l])
        yerr = copy(getattr(gdat, 'errr' + strgfeat + 'assc')[l])
    else:
        ydat = getattr(gdatmodi, 'this' + strgfeat + 'assc')[l]

    factplot = getattr(gdat, 'fact' + strgfeat + 'plot')
    
    xdat *= factplot
    xerr *= factplot
    if plotdiff:
        ydat = 100. * (ydat - xdat) / xdat
    else: 
        ydat *= factplot
        yerr *= factplot
    
    # plot all associations
    if plotdiff:
        indx = where(ydat > -100.)[0]
    else:
        indx = where(ydat > 0.)[0]
    if indx.size > 0:
        axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='black')
        
    # temp -- plot associations inside the comparison area
    if plotdiff:
        axis.axhline(0., ls='--', alpha=gdat.alphmrkr, color='black')
    else:
        axis.plot(factplot * binsplot, factplot * binsplot, ls='--', alpha=gdat.alphmrkr, color='black')
    
    lablfeat = getattr(gdat, 'labl' + strgfeat)
    lablfeatunit = getattr(gdat, 'labl' + strgfeat + 'unit')
    lablxdat = r'$%s^{%s}$%s' % (lablfeat, gdat.strgcatl, lablfeatunit)
    lablydat = '$%s^{samp}$%s' % (lablfeat, lablfeatunit)
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)
    if indx.size > 0 and getattr(gdat, 'scal' + strgfeat + 'plot') == 'logt':
        if not plotdiff:
            axis.set_yscale('log')
        axis.set_xscale('log')
    if plotdiff:
        limsyaxi = array([-100., 100.])
    else:
        limsyaxi = factplot * array([minmplot, maxmplot])
    axis.set_ylim(limsyaxi)
    axis.set_xlim([factplot * minmplot, factplot * maxmplot])
   
    plt.tight_layout()
    if plotdiff:
        strgtype = 'diff'
    else:
        strgtype = ''
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatassc' + strgfeat + '%spop%d' % (strgtype, l), nameinte='assc/')
    savefigr(gdat, gdatmodi, figr, path)
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
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 40
    numbfluxprox = len(gdat.indxpixlprox)
    bins = empty((gdat.numbprox, numbbins + 1))
    indxpixlproxsize = empty((numbfluxprox, gdat.numbpixl))
    for h in range(gdat.numbprox):
        for j in gdat.indxpixl:
            try:
                indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
            except:
                indxpixlproxsize[h, j] = gdat.numbpixl
        bins[h, :] = logspace(log10(amin(indxpixlproxsize[h, :])), log10(amax(indxpixlproxsize[h, :])), numbbins + 1)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbprox):
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
    
    
def plot_postbindmaps(gdat, indxpopltemp, strgbins, strgfeat=None):
    
    if strgfeat != None:
        numbparaplot = gdat.numbbinsplot
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
    
    if strgfeat != None:
        indxfeatsign = gdat.fittliststrgfeatsign[indxpopltemp].index(strgfeat)
    else:
        indxfeatsign = arange(len(gdat.fittliststrgfeatsign))

    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
    if numbrows == 1:
        axgr = [axgr]            
    for a, axrw in enumerate(axgr):
        if numbcols == 1:
            axrw = [axrw]
        for b, axis in enumerate(axrw):
            if strgfeat != None:
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
                temp = sum(gdat.pntsprob[indxpopltemp][:, :, indxlowr:indxuppr, indxfeatsign], 2).T
            else:
                temp = sum(sum(gdat.pntsprob[indxpopltemp], 2), 2).T
                
            if where(temp > 0.)[0].size > 0:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
            else:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                
            if strgfeat != None:
                bins = getattr(gdat, 'bins' + strgfeat)
            
            # superimpose true PS
            if gdat.trueinfo:
                
                if strgfeat != None:
                    truefeat = getattr(gdat, 'true' + strgfeat)[indxpopltemp][0, :] 
                    indxpnts = where((bins[indxlowr] < truefeat) & (truefeat < bins[indxuppr]))[0]
                    
                else:
                    indxpnts = arange(gdat.truenumbpnts[indxpopltemp])
                
                truecompsign = getattr(gdat, 'true' + gdat.namefeatsign)[indxpopltemp]
                
                mrkrsize = retr_mrkrsize(gdat, truecompsign[0, indxpnts])
                axis.scatter(gdat.anglfact * gdat.truelgal[indxpopltemp][0, indxpnts], gdat.anglfact * gdat.truebgal[indxpopltemp][0, indxpnts], \
                                                                                        s=mrkrsize, alpha=gdat.alphmrkr, marker='*', lw=2, color='g')

            if a == numbrows - 1:
                axis.set_xlabel(gdat.labllgaltotl)
            else:
                axis.set_xticklabels([])
            if b == 0:
                axis.set_ylabel(gdat.lablbgaltotl)
            else:
                axis.set_yticklabels([])

            draw_frambndr(gdat, axis)
            
            if strgbins != 'cumu':
                lablfeat = getattr(gdat, 'labl' + strgfeat)
                factfeatplot = getattr(gdat, 'fact' + strgfeat + 'plot')
                titl = tdpy.util.mexp(factfeatplot * bins[indxlowr]) + ' < $%s$ < ' % lablfeat + tdpy.util.mexp(factfeatplot * bins[indxuppr])
                axis.set_title(titl)
    
    if strgfeat != None:
        lablfeattotl = getattr(gdat, 'labl' + strgfeat + 'totl')
        plt.figtext(0.5, 0.95, '%s' % lablfeattotl, ha='center', va='center')
    axiscomm = figr.add_axes([0.87, 0.2, 0.02, 0.6])
    cbar = figr.colorbar(imag, cax=axiscomm)

    plt.subplots_adjust(left=0.18, top=.9, right=0.82, bottom=0.15, hspace=0.08, wspace=0.08)
    if strgbins == 'cumu':
        strgtemp = ''
    else:
        strgtemp = strgfeat
    path = getattr(gdat, 'path' + gdat.namesampdist + 'finl') + 'postbindmaps%s%s%d' % (strgbins, strgtemp, indxpopltemp) + '.pdf'
    figr.savefig(path)
    plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.binsangl)

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
        axis.set_xlabel(gdat.lablgangtotl)
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
            axis.hist(datacntstemp, gdat.binsdatacnts, alpha=0.5, color='black')
            axis.axvline(mean(datacntstemp), ls='--', color='black')
            axis.set_xscale('log')
            if m == gdat.numbevtt - 1:
                axis.set_xlabel(r'$k$')
            if m == 0 and gdat.numbener > 1:
                axis.set_title(gdat.strgener[i])
            if i == 0 and gdat.numbevtt > 1:
                axis.set_ylabel(gdat.strgevtt[m])
            axis.set_yscale('log')

    plt.tight_layout()
    figr.savefig(gdat.pathinit + 'datacntshist.pdf')
    plt.close(figr)
    
    
def plot_intr(gdat):
    
    if gdat.verbtype > 0:
        print 'Making PCAT introductory plots...'

    plot_grap(verbtype=1)
    plot_grap(plottype='ngal', verbtype=1)
    #plot_grap(plottype='lens', verbtype=1)
    #plot_grap(plottype='lensprim', verbtype=1)
    plot_grap(plottype='lensseco', verbtype=1)
    #plot_grap(plottype='chan', verbtype=1)
    
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
        psfntemp = copy(gdat.exprpsfn[0, :, 0, 0])
    else:
        psfntemp = copy(gdat.exprpsfn[0, :, 0])
    psfntemp /= amax(psfntemp)
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for k in range(gdat.numbprox + 1):
        if k == 0 or k == gdat.numbprox:
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
        axis.plot(gdat.binsangl * gdat.anglfact, gdat.binsprox[k] * psfntemp, label=labl, color=colr, alpha=alph)
        axis.set_xlim([amin(gdat.binsangl) * gdat.anglfact, amax(gdat.binsangl) * gdat.anglfact])
        if k > 0:
            axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
    axis.set_yscale('log')
    axis.set_xlabel(gdat.lablgangtotl)
    axis.set_ylabel(gdat.lablfluxtotl)

    limt = gdat.specfraceval * amax(gdat.binsprox[0] * psfntemp)
    maxmangltemp = interp(1e-1 * limt, gdat.binsprox[k] * psfntemp[::-1], gdat.binsangl[::-1] * gdat.anglfact)
    
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
    gdatmodi.thischro = zeros(gdat.numbchro)

    # data structure to hold the indices of model PS to be compared to the reference catalog 
    gdatmodi.indxmodlpntscomp = [[] for l in gdat.fittindxpopl]
    
    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamptotl:
        indxsampmosa = choice(gdat.indxsamptotl, size=numbsampmosa, replace=False)
        for l in gdat.fittindxpopl:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            
                            n = indxsampmosa[numbcols*a+b]
                            gdatmodi.thissampvarb = gdat.listsampvarb[n, :].flatten()
                            gdatmodi.thissamp = gdat.listsamp[n, :].flatten()
                            
                            if gdat.fittnumbtrap > 0:
                                gdatmodi.thisindxelemfull = gdat.listindxelemfull[n]
                                proc_samp(gdat, gdatmodi, 'this')

                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.labllgaltotl)
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.lablbgaltotl)
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
                    plt.subplots_adjust(left=0.1, top=.91, hspace=0.03, wspace=0.1, bottom=0.09)
                    if l == 1:
                        strg = ''
                    else:
                        strg = 'pop%d' % l
                    pathfinl = getattr(gdat, 'path' + gdat.namesampdist + 'finl')
                    if m == None:
                        path = pathfinl + 'mosa' + strg + 'ene%dA.pdf' % (gdat.indxenerincl[i])
                    else:
                        path = pathfinl + 'mosa' + strg + 'ene%devtt%d.pdf' % (gdat.indxenerincl[i], gdat.indxevttincl[m])
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

    if plottype == 'lensseco':
        listcolr = ['olive', 'black', 'olive', 'olive', 'olive', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', \
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
    if plottype == 'lens' or plottype == 'lensprim' or plottype == 'lensseco':
        grap.add_edges_from([ \
                             ('lenp', 'modl'), \
                             ('defsdistslop', 'defs'), \
                             ('defs', 'modl'), \
                             ('numbpnts','defs'), \
                            ])
    if plottype == 'lensprim' or plottype == 'lensseco':
        grap.add_edges_from([ \
                             ('asca', 'modl'), \
                             ('numbpnts','asca'), \
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
    if plottype.startswith('lens'):
        nameelem = r'\rm{sub}'
    else:
        nameelem = r'\rm{pts}'
    labl['numbpnts'] = '$N_{%s}$' % nameelem
    
    if plottype == 'igal':
        labl['ampldistslop'] = r'$\vec{\beta}$'
        labl['meanpnts'] = r'$\vec{\mu}_{%s}$' % nameelem
    else:
        labl['meanpnts'] = r'$\mu_{%s}$' % nameelem
    if plottype == 'chan' or plottype == 'ngal':
        labl['ampldistslop'] = r'$\beta$'
    if plottype == 'lens' or plottype == 'lensprim' or plottype == 'lensseco':
        labl['defsdistslop'] = r'$\beta$'
    
    if plottype == 'igal':
        labl['expodistslop'] = r'$\vec{\tau_\rm{E}}$'
        labl['sinddistslop'] = r'$\vec{\tau_\rm{s}}$'
    if plottype == 'chan':
        labl['sinddistslop'] = r'$\beta$'
    
    if plottype == 'igal':
        labl['spatdistslop'] = r'$\vec{\gamma}$'
    if plottype == 'lens' or plottype == 'lensprim' or plottype == 'lensseco':
        labl['lenp'] = r'$\vec{\chi}$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['bacp'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{\theta_1}$'
    labl['bgal'] = r'$\vec{\theta_2}$'
    if plottype == 'igal' or plottype == 'ngal' or plottype == 'chan':
        labl['sind'] = r'$\vec{s}$'
        labl['ampl'] = r'$\vec{f}$'
    else:
        labl['defs'] = r'$\vec{\alpha_\rm{s}}$'
    if plottype == 'lensprim' or plottype == 'lensseco':
        labl['asca'] = r'$\vec{\theta_\rm{s}}$'
        labl['acut'] = r'$\vec{\theta_\rm{c}}$'
        
    if plottype == 'igal':
        labl['expo'] = r'$\vec{E_{\rm{c}}}$'
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
    if plottype == 'lens' or plottype == 'lensprim' or plottype == 'lensseco':
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
    if plottype == 'lensprim' or plottype == 'lensseco':
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
    if plottype == 'lensprim' or plottype == 'lens' or plottype == 'lensseco':
        posi['defs'] = array([0.1, -0.0])
    if plottype == 'lensprim' or plottype == 'lensseco':
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

    size = 1000
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
    if plottype == 'lens' or plottype == 'lensprim' or plottype == 'lensseco':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanpnts', 'defsdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lenp'], node_color='y', node_size=size)
    if plottype == 'lens':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs'], node_color='g', node_size=size)
    if plottype == 'lensprim':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['ascadistslop', 'acutdistslop'], node_color='r', node_size=size)
    if plottype == 'lensprim' or plottype == 'lensseco':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs', 'asca', 'acut'], node_color='g', node_size=size)
    
    pathplot = os.environ["PCAT_DATA_PATH"] + '/imag/'
    plt.tight_layout()
    figr.savefig(pathplot + 'grap%s.pdf' % plottype)
    plt.close(figr)


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
    axis.set_xlabel(gdat.labllgaltotl)
    axis.set_ylabel(gdat.lablbgaltotl)

    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    figr.savefig(gdat.pathplot + 'thrs.pdf')
    plt.close(figr)
    

def plot_init(gdat):
        
    # make initial plots
    if gdat.makeplot:
        
        # ferm
        #plot_3fgl_thrs(gdat)
        #if gdat.exprtype == 'ferm':
        #    plot_fgl3(gdat)
        
        plot_datacntshist(gdat)
        if gdat.evalcirc != 'full':
            plot_indxprox(gdat)
        if gdat.evalcirc != 'full' and (gdat.elemtype == 'lght' or gdat.elemtype == 'clus'):
            plot_eval(gdat)
        
        # temp
        if gdat.makeplotintr:
            plot_intr(gdat)
            #plot_pert()
            #plot_king(gdat)
    
            if gdat.elemtype == 'lens':
                xdat = gdat.binsangl[1:] * gdat.anglfact
                lablxdat = gdat.lablgangtotl
                
                listdeflscal = array([4e-2, 4e-2, 4e-2]) / gdat.anglfact
                listanglscal = array([0.05, 0.1, 0.05]) / gdat.anglfact
                listanglcutf = array([1.,    1.,  10.]) / gdat.anglfact
                listasym = [False, False, False]
                listydat = []
                for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
                    listydat.append(retr_deflcutf(gdat.binsangl[1:], deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
                
                for scalxdat in ['self', 'logt']:
                    path = gdat.pathinitintr + 'deflcutf' + scalxdat + '.pdf'
                    tdpy.util.plot_gene(path, xdat, listydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, \
                                                                                lablydat=r'$\alpha_n$ [$^{\prime\prime}$]', limtydat=[1e-3, 1.5e-2], limtxdat=[None, 2.])
               
                # pixel-convoltuion of the Sersic profile
                xdat = gdat.binslgalsers * gdat.anglfact
                for n in range(gdat.numbindxsers + 1):
                    for k in range(gdat.numbhalfsers + 1):
                        path = gdat.pathinitintr + 'sersprofconv%04d%04d.pdf' % (gdat.binshalfsers[n], gdat.binsindxsers[k])
                        tdpy.util.plot_gene(path, xdat, gdat.sersprof[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxsoldtotl, \
                                                                                                                                                        limtydat=[1e-8, 1e010])
                        path = gdat.pathinitintr + 'sersprofcntr%04d%04d.pdf' % (gdat.binshalfsers[n], gdat.binsindxsers[k])
                        tdpy.util.plot_gene(path, xdat, gdat.sersprofcntr[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxsoldtotl)
               
                xdat = gdat.binsangl * gdat.anglfact
                listspec = array([1e-19, 1e-18, 1e-18, 1e-18]) / gdat.anglfact
                listsize = array([0.3, 1., 1., 1.]) / gdat.anglfact
                listindx = array([4., 2., 4., 10.])
                listydat = []
                listlegd = []
                for spec, size, indx in zip(listspec, listsize, listindx):
                    listydat.append(retr_sersprof(gdat, spec, gdat.binsangl, size, indx=indx) * gdat.anglfact**2 / (4. * pi)**2)
                    listlegd.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
                path = gdat.pathinitintr + 'sersprof.pdf'
                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxsoldtotl, \
                                                                 listlegd=listlegd, listhlin=1e-7, limtydat=[1e-8, 1e0])
            
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
                
                asca = 0.1 / gdat.anglfact
                acut = 1. / gdat.anglfact
    
                minmmass = zeros((numbreds + 1, numbreds + 1))
                maxmmass = zeros((numbreds + 1, numbreds + 1))
                for k, redshost in enumerate(gdat.binsredshost):
                    for n, redssour in enumerate(gdat.binsredssour):
                        if redssour > redshost:
                            adishost = gdat.adisobjt(redshost) * 1e3
                            adissour = gdat.adisobjt(redssour) * 1e3
                            adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
                            factmcutfromdefs = retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut)
                            minmmass[n, k] = log10(factmcutfromdefs * gdat.minmdefs)
                            maxmmass[n, k] = log10(factmcutfromdefs * gdat.maxmdefs)
               
                valulevl = linspace(7.5, 9., 10)
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                cont = axis.contour(gdat.binsredshost, gdat.binsredssour, minmmass, 10, colors='g', levels=valulevl)
                axis.clabel(cont, inline=1, fontsize=10, fmt='%.3g')
                axis.set_xlabel(r'$z_{\rm{hst}}$')
                axis.set_ylabel(r'$z_{\rm{src}}$')
                axis.set_title(r'$M_{c,min}$ [$M_{\odot}$]')
                path = gdat.pathinitintr + 'massredsminm.pdf'
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
                
                valulevl = linspace(9., 11., 20)
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=9., vmax=11.)
                cont = axis.contour(gdat.binsredshost, gdat.binsredssour, maxmmass, 10, colors='g', levels=valulevl)
                axis.clabel(cont, inline=1, fontsize=10, fmt='%.3g')
                axis.set_xlabel('$z_{hst}$')
                axis.set_ylabel('$z_{src}$')
                axis.set_title(r'$M_{c,max}$ [$M_{\odot}$]')
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
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                fracacutasca = logspace(-1., 2., 20)
                mcut = 1e8 * retr_mcutfrommscl(fracacutasca)
                axis.loglog(fracacutasca, mcut)
                axis.set_xlabel(r'$\tau_n$')
                axis.set_ylabel(gdat.lablmcuttotl)
                axis.axhline(1e8, ls='--')
                path = gdat.pathinitintr + 'mcut.pdf'
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)

    for i in gdat.indxener:
        for m in gdat.indxevtt:
            
            # temp
            if False and gdat.pixltype == 'cart' and (gdat.elemtype == 'lght' or gdat.elemtype == 'clus'):
                figr, axis, path = init_figr(gdat, None, 'datacntspeak', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.datacnts, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                axis.scatter(gdat.anglfact * gdat.meanlgalcart[gdat.indxxaximaxm], gdat.anglfact * gdat.meanbgalcart[gdat.indxyaximaxm], alpha=0.6, s=20, edgecolor='none')
                
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
    
            if gdat.correxpo:
                figr, axis, path = init_figr(gdat, None, 'expo', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.expo, '', 'expomaps', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
            
            for strgmodl in gdat.liststrgmodl:
                numbback = getattr(gdat, strgmodl + 'numbback')
                indxback = getattr(gdat, strgmodl + 'indxback')
                backcnts = getattr(gdat, strgmodl + 'backcnts')
                backcntstotl = getattr(gdat, strgmodl + 'backcntstotl')
                for c in indxback:
                    figr, axis, path = init_figr(gdat, None, 'backcnts', '', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, backcnts[c], '', 'datacnts', thisindxener=i, thisindxevtt=m)
                    make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
        
                if numbback > 1:
                    figr, axis, path = init_figr(gdat, None, 'backcntstotl', '', indxenerplot=i, indxevttplot=m)
                    imag = retr_imag(gdat, axis, backcntstotl, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                    make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
            
                figr, axis, path = init_figr(gdat, None, 'diffcntstotl', '', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.datacnts - backcntstotl, '', 'datacnts', thisindxener=i, thisindxevtt=m)
                make_cbar(gdat, axis, imag, i, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)
    
            if gdat.elemtype == 'lens' and gdat.datatype == 'mock':
                figr, axis, path = init_figr(gdat, None, 'modlcntsraww', 'true', indxenerplot=i, indxevttplot=m)
                imag = retr_imag(gdat, axis, gdat.truemodlcntsraww, '', 'datacnts', thisindxener=i, thisindxevtt=m, tdim=True)
                make_cbar(gdat, axis, imag, 0, tick=gdat.tickdatacnts, labl=gdat.labldatacnts)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)


def plot_defl(gdat, gdatmodi, strg, strgvarb='defl', strgcomp='', indxdefl=None, thisindxpopl=-1, multfact=1.):

    if indxdefl != None:
        strgvarb += 'sing'
    strgvarb = strgcomp + strgvarb
    
    defl = retr_fromgdat(gdat, gdatmodi, strg, strgvarb)
    
    defl *= multfact

    if indxdefl != None:
        defl = defl[:, :, :, indxdefl]
        strgvarb += '%04d' % indxdefl

    figr, axis, path = init_figr(gdat, gdatmodi, strgvarb, strg, indxpoplplot=thisindxpopl)
    make_catllabl(gdat, strg, axis)
    draw_frambndr(gdat, axis)
  
    defllgal = defl[:, :, 0]
    deflbgal = defl[:, :, 1]
    fact = 4
    axis.imshow(zeros((10, 10)))
    ptch = axis.quiver(gdat.anglfact * gdat.lgalgridcart[::fact, ::fact], gdat.anglfact * gdat.bgalgridcart[::fact, ::fact], \
                  gdat.fluxfactplot * defllgal[::fact, ::fact], gdat.fluxfactplot * deflbgal[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strg, axis)
    #plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strg, strgvarb, strgcbar=None, thisindxener=None, thisindxevtt=-1, tdim=False, thisindxpopl=-1, mometype='medi', intreval=False):
    
    if strgcbar == None:
        strgcbar = strgvarb
   
    if strgvarb == 'datacnts':
        strgplot = strgvarb
    else:
        if strg == 'post':
            strgtemp = mometype
        else:
            strgtemp = ''
        strgplot = strgtemp + strgvarb
    
    figr, axis, path = init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=thisindxener, indxevttplot=thisindxevtt, indxpoplplot=thisindxpopl, intreval=intreval)
    
    maps = retr_fromgdat(gdat, gdatmodi, strg, strgvarb, mometype=mometype)
    imag = retr_imag(gdat, axis, maps, strg, strgcbar, thisindxener=thisindxener, thisindxevtt=thisindxevtt, tdim=tdim)
    tick = getattr(gdat, 'tick' + strgcbar) 
    labl = getattr(gdat, 'labl' + strgcbar) 
    make_cbar(gdat, axis, imag, tick=tick, labl=labl)
    make_catllabl(gdat, strg, axis)
    supr_fram(gdat, gdatmodi, strg, axis, thisindxpopl)

    # Add two sliders for tweaking the parameters
    if intreval:
        print 'Interactive session began...'
        
        freq_slider = []
        for k, namefixp in enumerate(gdat.fittnamefixp):
            indxfixp = getattr(gdat, 'fittindxfixp' + namefixp)

            initvalu = gdatmodi.thissampvarb[indxfixp] * gdat.fittfactfixpplot[k]
            labl = getattr(gdat, 'labl' + namefixp)
            minm = getattr(gdat, 'minm' + namefixp) * gdat.fittfactfixpplot[k]
            maxm = getattr(gdat, 'maxm' + namefixp) * gdat.fittfactfixpplot[k]
            freq_slider_ax = figr.add_axes([0.08, 0.04 + 0.05 * k, 0.2, 0.04])
            freq_slider.append(Slider(freq_slider_ax, labl, minm, maxm, valinit=initvalu))
        def sliders_on_changed(val):
            print 'Slider changed.'
            for k, namefixp in enumerate(gdat.fittnamefixp):
                print namefixp
                print freq_slider[k].val
                gdatmodi.thissampvarb[k] = freq_slider[k].val / gdat.fittfactfixpplot[k]
            print
            proc_samp(gdat, gdatmodi, 'this')
            maps = retr_fromgdat(gdat, gdatmodi, strg, strgvarb, mometype=mometype)
            retr_imag(gdat, axis, maps, strg, strgcbar, thisindxener=thisindxener, thisindxevtt=thisindxevtt, tdim=tdim, imag=imag)
            for ptch in axis.get_children():
                if isinstance(ptch, mpl.patches.Circle) or isinstance(ptch, mpl.collections.PathCollection):#isinstance(ptch, mpl.lines.Line2D):
                    ptch.remove()
            supr_fram(gdat, gdatmodi, strg, axis, thisindxpopl)
            plt.show(block=False)

        for k, namefixp in enumerate(gdat.fittnamefixp):
            freq_slider[k].on_changed(sliders_on_changed)

        plt.show()
       
        inpt = raw_input("Press enter to continue...")
        plt.close()
        raise Exception('Interactive session ended...')
    else:
        plt.tight_layout()
        savefigr(gdat, gdatmodi, figr, path)
        plt.close(figr)
    
