# common imports
from __init__ import *

# internal functions
from util import *

def plot_samp(gdat, gdatmodi, strg):
   
    if gdat.shrtfram and strg == 'this':
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_genemaps(gdat, gdatmodi, strg, 'cntpdata', d, i, m)
                    plot_genemaps(gdat, gdatmodi, strg, 'cntpresi', d, i, m)
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
            
            if strg != 'post':
                numbelem = sampvarb[getattr(gdat, strgmodl + 'indxfixpnumbelem')].astype(int)
            liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
            if strg != 'post':
                indxelemfull = list(getattr(gdatobjt, strg + 'indxelemfull'))
                indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)
        numbback = getattr(gdat, strgmodl + 'numbback')
        backtype = getattr(gdat, strgmodl + 'backtype')
        indxback = getattr(gdat, strgmodl + 'indxback')
        convdiff = getattr(gdat, strgmodl + 'convdiff')
        listnamediff = getattr(gdat, strgmodl + 'listnamediff')
        listnameecomtotl = getattr(gdat, strgmodl + 'listnameecomtotl')
        unifback = getattr(gdat, strgmodl + 'unifback')
        listnameback = getattr(gdat, strgmodl + 'listnameback')
        lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
        
        #cntpback = getattr(gdatobjt, strg + 'cntpback')
        #cntpbacktotl = getattr(gdatobjt, strg + 'cntpbacktotl')
            
        if gdatmodi != None:
            strgswep = '_%09d' % gdatmodi.cntrswep
        else:
            strgswep = ''
        
        if gdatmodi != None:
            gdatobjt = gdatmodi
        else:
            gdatobjt = gdat
   
        # plots
        ## histograms of the number of counts per pixel
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    for nameecom in listnameecomtotl:
                        name = 'histcntp' + nameecom + 'reg%dene%devt%d' % (d, i, m)
                        plot_gene(gdat, gdatmodi, strg, name, 'meancntpdata', scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, \
                                                                                                    lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbpixl])

        ## highest amplitude element
        # temp
        if numbtrap > 0:
            if False:
                if gdat.elemtype == 'lght':
                    plot_brgt(gdat, gdatmodi, strg)
        
            # completeness and false discovery rate
            if strg != 'true' and gdat.allwrefr:
                for strgclas in ['cmpl', 'fdis']:   
                    nameinte = strgclas + '/'
                    for d in gdat.indxregi:
                        for q in gdat.indxrefr:
                            for l in gdat.fittindxpopl:
                                lablydat = getattr(gdat, 'labl' + strgclas + 'reg%dref%dpop%d' % (d, q, l))
                                strgindxydat = 'reg%dref%dpop%d' % (d, q, l)
                                for strgfeat in gdat.refrliststrgfeat[q]:
                                    if strgclas == 'fdis' and not strgfeat in liststrgfeatodim[l]:
                                        continue
                                    if not strgfeat.startswith('spec') and not strgfeat.startswith('defl') and not strgfeat in gdat.refrliststrgfeatonly[q][l] and \
                                                                        not (gdat.datatype == 'mock' and (strgfeat.endswith('pars') or strgfeat.endswith('nrel'))):
                                        factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                                        lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                                        scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                                        limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                                        plot_gene(gdat, gdatmodi, strg, strgclas + strgfeat + strgindxydat, 'mean' + strgfeat, lablxdat=lablxdat, \
                                                  lablydat=lablydat, factxdat=factxdat, plottype='errr', \
                                                  scalxdat=scalxdat, limtydat=[0., 1.], limtxdat=limtxdat, \
                                                  omittrue=True, nameinte=nameinte)
                
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
    
            # transdimensional element features projected onto the data axes
            if not (strg == 'post' and not gdat.condcatl):
                
                if gdat.numbener > 1 and gdat.elemtype == 'lght':
                    
                    # PS spectra
                    if strg == 'post':
                        specplot = [empty((gdat.numbenerplot, gdat.numbstkscond))]
                        for r in gdat.indxstkscond:
                            specplot[0][:, r] = gdat.dictglob['poststkscond'][r]['specplot'][0, :]
                    else:
                        specplot = getattr(gdatobjt, strg + 'specplot')
                    for d in gdat.indxregi:
                        for l in range(len(specplot)):
                            listxdat = []
                            listplottype = []
                            
                            for k in range(specplot[d][l].shape[-1]):
                                listxdat.append(gdat.meanenerplot)
                                listplottype.append('line')
                            
                            for specconvunit in gdat.listspecconvunit:
                                listydat = []
                                
                                for k in range(specplot[d][l].shape[-1]):
                                    specplottemp = specplot[d][l]
                                    if strgmodl == 'true':
                                        specplottemp = copy(specplottemp[0, :, k])
                                    else:
                                        specplottemp = copy(specplottemp[:, k])
                                    if specconvunit[0] == 'ene1':
                                        specplottemp *= gdat.meanenerplot
                                    if specconvunit[0] == 'ene2':
                                        specplottemp *= gdat.meanenerplot**2
                                    if specconvunit[0] == 'ene3':
                                        # temp
                                        pass
                                    listydat.append(specplottemp)
                                
                                lablydat = getattr(gdat, 'lablflux' + specconvunit[0] + specconvunit[1] + 'totl')
                                if specconvunit[1] == gdat.nameenerunit:
                                    factydat = 1.
                                else:
                                    factydat = getattr(gdat, 'fact' + specconvunit[1] + gdat.nameenerunit)
                                strgtemp = specconvunit[0] + specconvunit[1]
                                if specconvunit[0] == 'ene3':
                                    strgtemp += specconvunit[2]
                                path = pathtemp + strg + 'specreg%dpop%d%s%s.pdf' % (d, l, strgtemp, strgswep)
                                for ydat in listydat:
                                    ydat *= factydat
                                limtydat = [amin(gdat.minmspec) * factydat, amax(gdat.maxmspec) * factydat]
                                tdpy.util.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, colr=colr, alph=alph, \
                                                           plottype=listplottype, limtxdat=[gdat.minmener, gdat.maxmener], lablydat=lablydat, \
                                                           limtydat=limtydat)
                
                if gdat.elemtype == 'lens':

                    ## deflection profiles
                    if gdat.variasca and gdat.variacut:
                        xdat = gdat.meanangl * gdat.anglfact
                        lablxdat = gdat.lablgangtotl
                        if strg == 'post':
                            deflprof = [empty((gdat.numbangl, gdat.numbstkscond))]
                            asca = [empty(gdat.numbstkscond)]
                            acut = [empty(gdat.numbstkscond)]
                            for r in gdat.indxstkscond:
                                deflprof[0][:, r] = gdat.dictglob['poststkscond'][r]['deflprof'][0, :]
                                asca[0][r] = gdat.dictglob['poststkscond'][r]['asca'][0]
                                acut[0][r] = gdat.dictglob['poststkscond'][r]['acut'][0]
                        else:
                            deflprof = getattr(gdatobjt, strg + 'deflprof')
                            asca = getattr(gdatobjt, strg + 'asca')
                            acut = getattr(gdatobjt, strg + 'acut')

                        for d in gdat.indxregi:
                            for l in range(len(deflprof[d])):
                                listydat = []
                                listvlinfrst = []
                                listvlinseco = []
    
                                if strgmodl == 'true':
                                    print 'dl'
                                    print d, l
                                    print 'deflprof'
                                    print len(deflprof)
                                    print len(deflprof[0])
                                    print 
                                    deflproftemp = deflprof[d][l][0, :, :]
                                else:
                                    deflproftemp = deflprof[d][l]
                                for k in range(deflprof[d][l].shape[-1]):
                                    listydat.append(deflproftemp[:, k] * gdat.anglfact)
                                    if strg == 'true':
                                        ascatemp = asca[d][l][0, k]
                                        acuttemp = acut[d][l][0, k]
                                    else:
                                        ascatemp = asca[d][l][k]
                                        acuttemp = acut[d][l][k]
                                    listvlinfrst.append(ascatemp * gdat.anglfact) 
                                    listvlinseco.append(acuttemp * gdat.anglfact)
                                
                                indxfixpbeinhost = getattr(gdat, strgmodl + 'indxfixpbeinhost')
                                beinhost = retr_fromgdat(gdat, gdatmodi, strg, 'sampvarb', indxvarb=indxfixpbeinhost)
                                listydat.append(xdat * 0. + gdat.anglfact * beinhost[d])
                                path = pathtemp + strg + 'deflsubhreg%dpop%d%s.pdf' % (d, l, strgswep)
                                limtydat = [1e-3, 1.]
                                limtxdat = [1e-3, 1.]
                                tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, drawdiag=True, limtydat=limtydat, \
                                        limtxdat=limtxdat, colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
                    
            if gdat.datatype == 'mock':
                if gdat.elemtype == 'lens':
                    ## radial mass budget
                    xdat = gdat.meananglhalf * gdat.anglfact
                    lablxdat = gdat.lablanglfromhosttotl
                    for namecalc in ['delt', 'intg']:
                        
                        # host and subhalo masses
                        limtydat = [gdat.minmmcut, getattr(gdat, 'maxmmasshost' + namecalc + 'bein')]
                        lablydat = getattr(gdat, 'lablmass' + namecalc + 'totl')

                        path = retr_plotpath(gdat, gdatmodi, strg, 'mass%s%s.pdf' % (namecalc, strgswep))
                        listydat = [retr_fromgdat(gdat, gdatmodi, strg, 'masshost' + namecalc), \
                                                        retr_fromgdat(gdat, gdatmodi, strg, 'masssubh' + namecalc)]
                        if strg == 'post':
                            listyerr = [retr_fromgdat(gdat, gdatmodi, strg, 'masshost' + namecalc, mometype='errr'), \
                                                        retr_fromgdat(gdat, gdatmodi, strg, 'masssubh' + namecalc, mometype='errr')]
                        else:
                            listyerr = None
                        tdpy.util.plot_gene(path, xdat, listydat, yerr=listyerr, scalydat='logt', lablxdat=lablxdat, limtydat=limtydat, \
                                                                                                lablydat=lablydat, listlinestyl=['-', '--'], colr=colr)

                        # subhalo mass fraction
                        limtydat = [1e-3, 0.1]
                        lablydat = getattr(gdat, 'lablfracsubh' + namecalc + 'totl')

                        path = retr_plotpath(gdat, gdatmodi, strg, 'fracsubh%s%s.pdf' % (namecalc, strgswep))
                        listydat = retr_fromgdat(gdat, gdatmodi, strg, 'fracsubh' + namecalc)
                        if strg == 'post':
                            listyerr = retr_fromgdat(gdat, gdatmodi, strg, 'fracsubh' + namecalc, mometype='errr')
                        else:
                            listyerr = None
                        tdpy.util.plot_gene(path, xdat, listydat, yerr=listyerr, scalydat='logt', lablxdat=lablxdat, limtydat=limtydat, lablydat=lablydat, colr=colr)

            if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                ## PSF radial profile
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        if gdat.fittoaxitype:
                            indxydat = [i, slice(None), m, 0]
                        else:
                            indxydat = [i, slice(None), m]
                        strgindxydat = 'ene%d%evtd' % (i, m)
                        lablxdat = gdat.lablgangtotl
                        limtydat= array([1e-3, 1e3]) * gdat.anglfact**2
                        plot_gene(gdat, gdatmodi, strg, 'psfn', 'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalydat='logt', \
                                                                                       factxdat=gdat.anglfact, lablxdat=lablxdat, lablydat=r'$\mathcal{P}$', limtydat=limtydat)
                        if gdat.fittoaxitype or gdat.trueoaxitype:
                            plot_gene(gdat, gdatmodi, strg, 'factoaxi', 'binsoaxi', indxydat=[i, m, slice(None)], strgindxydat=strgindxydat, \
                                                                                            factxdat=gdat.anglfact, lablxdat=lablxdat, lablydat=r'$f(\phi)$')
    
                # number of background counts inside PSF FWHM
                # temp
                if False:
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            strgindxydat = '%d%d' % (i, m)
                            plot_gene(gdat, gdatmodi, strg, 'cntsbackfwhm', 'meancnts', indxydat=[i, m], strgindxydat=strgindxydat, \
                                                                                                                    lablxdat=gdat.lablcnts, lablydat=gdat.lablcntsbackfwhm)
    
            # element feature correlations
            for d in gdat.indxregi:
                for l in indxpopl:
                    if strg != 'true' and gdat.refrinfo and gdat.allwrefr:
                        for strgfeat in gdat.fittliststrgfeatodim[l]:
                            for q in gdat.indxrefr:
                                if not strgfeat in gdat.refrliststrgfeat[q][l]:
                                    continue
                                plot_scatassc(gdat, gdatmodi, strg, q, l, strgfeat, d)
                                plot_scatassc(gdat, gdatmodi, strg, q, l, strgfeat, d, plotdiff=True)
                    for a, strgfrst in enumerate(liststrgfeatcorr[l]):
                        for b, strgseco in enumerate(liststrgfeatcorr[l]):
                            if a < b:
                                for strgplottype in ['hist', 'scat']:
                                    if strgmodl == 'true' and strgplottype == 'hist':
                                        continue
                                    if not gdat.condcatl and strgplottype == 'scat':
                                        continue
                                    plot_elemtdim(gdat, gdatmodi, strg, d, l, strgplottype, strgfrst, strgseco)
    
            # element feature histograms
            if not (strg == 'true' and gdat.datatype == 'inpt'):
                limtydat = gdat.limtydathistfeat
                for l in indxpopl:
                    strgindxydat = 'pop%d' % l
    
                    for strgfeat in liststrgfeatodim[l] + gdat.liststrgfeatplot:
                        
                        if strgfeat in gdat.liststrgfeatplot:
                            indxydat = [l, slice(None), gdat.indxenerpivt[0]]
                        else:
                            indxydat = [l, slice(None)]
                        
                        scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                        factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                        lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                        limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                        
                        # for true model, also plot the significant elements only
                        # temp
                        #if strgmodl == 'true' and strgfeat in gdat.listnamefeatpars:
                        #    listname = ['hist' + strgfeat, 'hist' + strgfeat + 'pars']
                        #else:
                        #    listname = ['hist' + strgfeat]

                        listname = ['hist' + strgfeat + 'reg%dpop%d' % (d, l)]
                        for name in listname:
                            listydattype = ['totl', 'sden']
                            for ydattype in listydattype:
                                
                                # plot the surface density of elements only for the amplitude feature
                                if strgfeat != gdat.namefeatampl and ydattype == 'sden':
                                    continue
        
                                ## plot the surface density of elements
                                if ydattype == 'sden':
                                    if gdat.sdenunit == 'degr':
                                        factydat = (pi / 180.)**2 / (2. * gdat.maxmgang)**2
                                        lablydat = r'$\Sigma_{%s}$ [deg$^{-2}$]' % gdat.lablelemextn
                                    if gdat.sdenunit == 'ster':
                                        factydat = 1. / (2. * gdat.maxmgang)**2
                                        lablydat = r'$\Sigma_{%s}$ [sr$^{-2}$]' % gdat.lablelemextn
                                ## plot the total number of elements
                                if ydattype == 'totl':
                                    factydat = 1.
                                    lablydat = r'$N_{%s}$' % gdat.lablelemextn
                                
                                print 'plot_samp()'
                                print 'element feature'
                                print 'strg'
                                print strg
                                print 'factydat'
                                print factydat
                                print 'name'
                                print name
                                print 'ydattype'
                                print ydattype
                                print 'strgfeat'
                                print strgfeat
                                print
                                plot_gene(gdat, gdatmodi, strg, name, 'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                  lablydat=lablydat, factxdat=factxdat, histodim=True, factydat=factydat, ydattype=ydattype, \
                                                  scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                  #indxydat=indxydat, strgindxydat=strgindxydat, \
                                                  nameinte='histodim/')
        
        # data and model count scatter
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    plot_scatcntp(gdat, gdatmodi, strg, d, i, m)
       
        # temp -- restrict other plots to indxmodlelemcomp
        if gdat.enerbins:
            for specconvunit in gdat.listspecconvunit:
                if not (isinstance(backtype[0], str) and backtype[0].startswith('mpol')):
                    for d in gdat.indxregi:
                        plot_sbrt(gdat, gdatmodi, strg, d, specconvunit)
       
        # data count maps
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_genemaps(gdat, gdatmodi, strg, 'cntpdata', d, i, m)
                    if numbpopl > 1:
                        for l in indxpopl:
                            plot_genemaps(gdat, gdatmodi, strg, 'cntpdata', d, i, m, indxpoplplot=l)
        
        ## spatial priors
        # temp
        if numbtrap > 0:
            liststrgfeatmodu = getattr(gdat, strgmodl + 'liststrgfeatmodu')
            liststrgpdfnmodu = getattr(gdat, strgmodl + 'liststrgpdfnmodu')
            for l in indxpopl:
                for strgfeat, strgpdfn in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]):
                    if strgpdfn == 'tmplreln':
                        plot_genemaps(gdat, gdatmodi, 'fitt', 'lpdfspatpriointp', tdim=True)
                    if strgpdfn == 'tmplgaum':
                        plot_genemaps(gdat, gdatmodi, strg, 'lpdfspatpriointp', tdim=True)
    
        # model count maps
        ## backgrounds
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    for c in indxback:
                        if isinstance(backtype[c], str) and backtype[c].startswith('mpol'):
                            continue
                        if not unifback[c]:
                            plot_genemaps(gdat, gdatmodi, strg, 'cntpback%04d' % c, d, i, m, strgcbar='cntpdata')
        
        ## count error
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and gdat.calcerrr and numbtrap > 0:
            if strg != 'true':
                for d in gdat.indxregi:
                    for i in gdat.indxener:
                        plot_genemaps(gdat, gdatmodi, strg, 'cntperrr', d, i, -1, strgcbar='cntpresi')
        
        ## diffuse components 
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for k, name in enumerate(listnamediff):
                    plot_genemaps(gdat, gdatmodi, strg, 'cntp' + name, d, i, strgcbar='cntpdata')
       
        ## model and residual
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_genemaps(gdat, gdatmodi, strg, 'cntpmodl', d, i, m, strgcbar='cntpdata')
                    plot_genemaps(gdat, gdatmodi, strg, 'cntpresi', d, i, m)
       
        # likelihood
        if strg != 'true':
            for d in gdat.indxregi:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_genemaps(gdat, gdatmodi, strg, 'llik', d, i, m, strgcbar='llikmaps')
        
        if lensmodltype != 'none':
            ## lensing signal to noise
            if strgmodl == 'true':
                for d in gdat.indxregi:
                    for i in gdat.indxener:
                        plot_genemaps(gdat, gdatmodi, strg, 's2nr', d, i, -1)
            for d in gdat.indxregi:
                plot_genemaps(gdat, gdatmodi, strg, 'magn', d, tdim=True)
                plot_genemaps(gdat, gdatmodi, strg, 'conv', d, tdim=True)
                if numbtrap > 0:
                    plot_genemaps(gdat, gdatmodi, strg, 'convelem', d, tdim=True)
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strg, 'cntplens', d, i, strgcbar='cntpdata', tdim=True)
                    plot_genemaps(gdat, gdatmodi, strg, 'cntplensgradmgtd', d, i, strgcbar='cntpdata', tdim=True)
        
        if gdat.penalpridiff:
            for d in gdat.indxregi:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_gene(gdat, gdatmodi, strg, 'psecodimdatapntsreg%dene%devt%d' % (d, i, m), 'meanmpolodim', lablxdat='$l$', lablydat='$P_{resi}(l)$', \
                                                                                                                 limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                        plot_gene(gdat, gdatmodi, strg, 'psecodimdatapntsprioreg%dene%devt%d' % (d, i, m), 'meanmpolodim', lablxdat='$l$', \
                                                                                           lablydat='$P_{prio}(l)$', limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
            
        if lensmodltype != 'none':
            for d in gdat.indxregi:
                indxydat = [d, slice(None)]
                strgindxydat = 'reg%d' % d
                plot_gene(gdat, gdatmodi, strg, 'convpsecodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P(k)$', limtydat=[1e-1, 1e2], \
                                                                                               scalxdat='logt', scalydat='logt', indxydat=indxydat, strgindxydat=strgindxydat)
                plot_gene(gdat, gdatmodi, strg, 'histdefl', 'meandefl', scal='self', lablxdat=r'$\alpha$ [arcsec]', lablydat=r'$N_{pix}$', factxdat=gdat.anglfact, \
                                                                                                             strgindxydat=strgindxydat, indxydat=indxydat, histodim=True)
        if gdat.elemtype == 'lens':
            if numbtrap > 0:
                for d in gdat.indxregi:
                    indxydat = [d, slice(None)]
                    strgindxydat = 'reg%d' % d
                    plot_gene(gdat, gdatmodi, strg, 'convpsecelemodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P_{sub}(k)$', \
                                                                      strgindxydat=strgindxydat, indxydat=indxydat, limtydat=[1e-5, 1e-1], scalxdat='logt', scalydat='logt')
                    plot_gene(gdat, gdatmodi, strg, 'histdeflelem', 'meandeflelem', scal='self', lablxdat=r'$\alpha$ [arcsec]', \
                                                                      strgindxydat=strgindxydat, indxydat=indxydat, lablydat=r'$N_{pix}$', factxdat=gdat.anglfact, histodim=True)
    
        if lensmodltype != 'none':
            for d in gdat.indxregi:
                # gradient of the lens emission
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_defl(gdat, gdatmodi, strg, d, 'cntplensgrad', indxenerplot=i, indxevttplot=m)
            
                # overall deflection field
                plot_defl(gdat, gdatmodi, strg, d, multfact=0.1)
                
                # deflection field due to individual lenses
                for k in range(numbdeflsingplot):  
                    if k == 0:
                        multfact = 0.1
                    elif k == 1:
                        multfact = 1.
                    elif k >= 2:
                        multfact = 10.
                    plot_defl(gdat, gdatmodi, strg, d, indxdefl=k, multfact=multfact)
                
                # residual deflection field
                if strg != 'true' and gdat.datatype == 'mock':
                    plot_defl(gdat, gdatmodi, strg, d, strgcomp='resi', multfact=100.)
                    for k in range(numbdeflsingplot):
                        plot_defl(gdat, gdatmodi, strg, d, strgcomp='resi', indxdefl=k, multfact=100.)
                    
                    if numbtrap > 0:
                        plot_genemaps(gdat, gdatmodi, strg, 'convelemresi', d, tdim=True)
                        plot_genemaps(gdat, gdatmodi, strg, 'convelemresiperc', d, tdim=True)
                    plot_genemaps(gdat, gdatmodi, strg, 'magnresi', d, tdim=True)
                    plot_genemaps(gdat, gdatmodi, strg, 'magnresiperc', d, tdim=True)
    

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
        sizefigrydat = max(thisnumbplot * gdat.plotsize / 4., gdat.plotsize / 2.)
        figr, axgr = plt.subplots(thisnumbplot, 1, figsize=(gdat.plotsize, sizefigrydat), sharex='all')
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
    for d in gdat.indxregi:
        plot_mosa(gdat, d)
    
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
        if gdat.fittnumbback == 2 and gdat.fittspecback == [False, False]:
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


def plot_sbrt(gdat, gdatmodi, strg, indxregiplot, specconvunit):
    
    if strg == 'true':
        strgmodl = strg
    else:
        strgmodl = 'fitt'
    
    backtype = getattr(gdat, strgmodl + 'backtype')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    specback = getattr(gdat, strgmodl + 'specback')
    indxback = getattr(gdat, strgmodl + 'indxback')
    maxmgang = getattr(gdat, strgmodl + 'maxmgang')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
    numblablsbrt = getattr(gdat, strgmodl + 'numblablsbrt')
    numblablsbrtspec = getattr(gdat, strgmodl + 'numblablsbrtspec')
    
    for b, namespatmean in enumerate(gdat.listnamespatmean):
    
        listydat = zeros((numblablsbrtspec, gdat.numbener))
        listyerr = zeros((2, numblablsbrtspec, gdat.numbener))
        
        cntr = 0
        indxvarb = [b, indxregiplot, slice(None)]

        for c in indxback:
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtback%04dmean' % c, indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtback%04dmean' % c, indxvarb=indxvarb, mometype='errr')
            cntr += 1
        
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and numbtrap > 0:
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtpntsmean', indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtpntsmean', indxvarb=indxvarb, mometype='errr')
            cntr += 1
        
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtpntssubtmean', indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtpntssubtmean', indxvarb=indxvarb, mometype='errr')
            cntr += 1
        
        if hostemistype != 'none':
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrthostmean', indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrthostmean', indxvarb=indxvarb, mometype='errr')
            cntr += 1
        
        if lensmodltype != 'none':
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtlensmean', indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtlensmean', indxvarb=indxvarb, mometype='errr')
            cntr += 1
        
        cntrdata = cntr

        ## data
        listydat[cntr, :] = gdat.sbrtdatamean[b, indxregiplot, :]
        cntr += 1
        
        ## total model
        if numblablsbrt > 1:
            listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtmodlmean', indxvarb=indxvarb)
            if strg == 'post':
                listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strg, 'sbrtmodlmean', indxvarb=indxvarb, mometype='errr')
            cntr += 1

        # plot energy spectra of the data, background model components and total background
        if gdat.numbener > 1:
            
            listlegdsbrtspec = getattr(gdat, strgmodl + 'listlegdsbrtspec')

            figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
            
            listmrkr = ['o', '>', 's', 'h', '*', 'p', 'x']

            # determine the energy scaling factor
            if specconvunit[0] == 'ene0':
                factener = 1.
            if specconvunit[0] == 'ene1':
                factener = gdat.meanener
            if specconvunit[0] == 'ene2':
                factener = gdat.meanener**2
            if specconvunit[0] == 'ene3':
                # temp
                pass
                factener = 1.
                #indxenerintv = where((gdat.meanener < specconvunit[4]) & (gdat.meanener > specconvunit[3]))[0]
                #ener = concatenate((array([specconvunit[3]]), gdat.meanener[indxenerintv], array([specconvunit[4]])))
                #
                #for k in range(3):
                #    if k == 0:
                #        ydattemp = 
                #    ydatminmener = interp(specconvunit[3], gdat.meanener, ydat)
                #    ydatmaxmener = interp(specconvunit[4], gdat.meanener, ydat)
                #    ydat = concatenate((array([ydatminmener]), ydat[indxenerintv], array([ydatmaxmener])))
                #    ydat = trapz(ydat, gdat.meanener)
                #
                #yerrminmener = interp(specconvunit[3], gdat.meanener, yerr, axis=1)
                #yerrmaxmener = interp(specconvunit[4], gdat.meanener, yerr, axis=1)
                #ydat = stack((array([yerrminmener]), ydat[indxenerintv], array([yerrmaxmener])))
                #
                #
                #yerr = trapz(yerr, gdat.meanener)

            # plot reference spectra
            if gdat.listspecrefrplot != None:
                for k in range(len(gdat.listspecrefrplot)):
                    axis.plot(gdat.listenerrefrplot[k], gdat.listspecrefrplot[k] * factener, label=gdat.listlablrefrplot[k], color='m')

            if specconvunit[1] == gdat.nameenerunit:
                factydat = 1.
            else:
                factydat = getattr(gdat, 'fact' + specconvunit[1] + gdat.nameenerunit)
            
            xdat = gdat.meanener
            for k in range(numblablsbrtspec):
                mrkr = listmrkr[k]
                if k == cntrdata:
                    colr = 'black'
                    linestyl = '-'
                else:
                    linestyl = '--'
                    if strg == 'true':
                        colr = 'g'
                    elif strg == 'this':
                        colr = 'b'
                    elif strg == 'post':
                        colr = 'black'
                
                ydat = copy(listydat[k, :])
                yerr = copy(listyerr[:, k, :])
                
                ydat *= factener
                yerr *= factener
                
                ydat *= factydat
                yerr *= factydat

                axis.errorbar(xdat, ydat, yerr=yerr, label=listlegdsbrtspec[k], color=colr, marker=mrkr, ls=linestyl, markersize=15)
        
            axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])

            limtydat = array([1e-4 * amin(listydat[cntrdata, :] * factener), 1e3 * amax(listydat[cntrdata, :] * factener)]) * factydat
            axis.set_ylim(limtydat)
            axis.set_yscale('log')
            axis.set_xlabel(gdat.lablenertotl)
            axis.set_xscale('log')
            labl = getattr(gdat, 'lablsbrt' + specconvunit[0] + specconvunit[1] + 'stertotl')
            axis.set_ylabel(labl)
            make_legd(axis, numbcols=2)
            
            plt.tight_layout()
            path = retr_plotpath(gdat, gdatmodi, strg, 'sdenmean%s%s%s' % (namespatmean, specconvunit[0], specconvunit[1]))
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
        axis.scatter(gdat.truefluxbrgt, gdat.truefluxbrgtassc, alpha=gdat.alphmrkr, color=gdat.listcolrrefr[0], label=gdat.legdrefr[0])
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
        

def plot_elemtdim(gdat, gdatmodi, strg, indxregiplot, indxpoplplot, strgplottype, strgfrst, strgseco, strgmome='medi'):
    
    sizelarg = 10
    sizesmll = 1
    
    legdmome = getattr(gdat, 'legd' + strgmome)
    
    limtfrst = getattr(gdat, 'limt' + strgfrst + 'plot')
    limtseco = getattr(gdat, 'limt' + strgseco + 'plot')
    factplotfrst = getattr(gdat, 'fact' + strgfrst + 'plot')
    factplotseco = getattr(gdat, 'fact' + strgseco + 'plot')

    if strg == 'true':  
        strgmodl = strg
    else:
        strgmodl = 'fitt'

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if strg == 'post':
        labl = gdat.legdsampdist + ' ' + legdmome
        if strgplottype == 'hist':
            varb = getattr(gdat, strgmome + 'hist' + strgfrst + strgseco + 'reg%dpop%d' % (indxregiplot, indxpoplplot))
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
                        varbfrst[cntr] = gdat.dictglob['poststkscond'][r][strgfrst][indxpoplplot] * getattr(gdat, 'fact' + strgfrst + 'plot')
                        varbseco[cntr] = gdat.dictglob['poststkscond'][r][strgseco][indxpoplplot] * getattr(gdat, 'fact' + strgseco + 'plot')
                        cntr += 1
                axis.scatter(varbfrst, varbseco, alpha=gdat.alphmrkr, color='black', label=gdat.legdsamp)
    
    if strg == 'this' or strg == 'mlik':
        if strgplottype == 'hist':
            meanfrst = getattr(gdat, 'bins' + strgfrst) * getattr(gdat, 'fact' + strgfrst + 'plot')
            meanseco = getattr(gdat, 'bins' + strgseco) * getattr(gdat, 'fact' + strgseco + 'plot')
            hist = getattr(gdatmodi, strg + 'hist' + strgfrst + strgseco + 'reg%dpop%d' % (indxregiplot, indxpoplplot))
            imag = axis.pcolor(meanfrst, meanseco, hist.T, cmap='Blues', label=gdat.legdsamp, alpha=gdat.alphmrkr)
        else:
            try:
                varbfrst = getattr(gdatmodi, 'this' + strgfrst)[indxpoplplot] * getattr(gdat, 'fact' + strgfrst + 'plot')
                varbseco = getattr(gdatmodi, 'this' + strgseco)[indxpoplplot] * getattr(gdat, 'fact' + strgseco + 'plot')
            except:
                varbfrst = array([limtfrst[0] * factplotfrst * 0.1])
                varbseco = array([limtseco[0] * factplotseco * 0.1])
            axis.scatter(varbfrst, varbseco, alpha=gdat.alphmrkr, color='b', label=gdat.legdsamp)
    
    # reference elements
    if gdat.allwrefr:
        for q in gdat.indxrefr:
            try:
                refrvarbfrst = getattr(gdat, 'refrfeat')[indxregiplot][q][strgfrst] * getattr(gdat, 'fact' + strgfrst + 'plot')
                refrvarbseco = getattr(gdat, 'refrfeat')[indxregiplot][q][strgseco] * getattr(gdat, 'fact' + strgseco + 'plot')
            except:
                refrvarbfrst = array([limtfrst[0] * factplotfrst * 0.1])
                refrvarbseco = array([limtseco[0] * factplotseco * 0.1])
            axis.scatter(refrvarbfrst, refrvarbseco, alpha=gdat.alphmrkr, color=gdat.listcolrrefr[q], label=gdat.legdrefr[q], s=sizelarg)

    plot_sigmcont(gdat, axis, indxpoplplot, strgfrst=strgfrst, strgseco=strgseco)
    
    scalfrst = getattr(gdat, 'scal' + strgfrst + 'plot')
    scalseco = getattr(gdat, 'scal' + strgseco + 'plot')
    if scalfrst == 'logt':
        axis.set_xscale('log')
    if scalseco == 'logt':
        axis.set_yscale('log')
    
    axis.set_xlabel(getattr(gdat, 'labl' + strgfrst + 'totl'))
    axis.set_ylabel(getattr(gdat, 'labl' + strgseco + 'totl'))
    axis.set_xlim(limtfrst * factplotfrst)
    axis.set_ylim(limtseco * factplotseco)

    plt.tight_layout()
    if strg == 'post':
        strgmometemp = strgmome
    else:
        strgmometemp = ''
    path = retr_plotpath(gdat, gdatmodi, strg, '%s%s%s%sreg%dpop%d' % (strgmometemp, strgplottype, strgfrst, strgseco, indxregiplot, indxpoplplot), nameinte=strgplottype + 'tdim/')
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
                     scal=None, scalxdat=None, scalydat=None, limtxdat=None, limtydat=None, omittrue=False, nameinte='', \
                     lablxdat='', lablydat='', factxdat=1., factydat=1., histodim=False, offslegd=None, tdim=False, ydattype='totl'):
   
    if strg == 'true':
        strgmodl = strg
    else:
        strgmodl = 'fitt'

    if scal == None:
        if scalxdat == None:
            scalxdat = 'linr'
        if scalydat == None:
            scalydat = 'linr'
    else:
        scalxdat = scal
        scalydat = scal

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
    
    # temp
    xerr = zeros((2, xdat.size))
    
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
        
        temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, marker='o', ls='', markersize=5, color='black', label=labl, lw=1, capsize=10)
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
                temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, marker='o', ls='', markersize=5, color='b', label=legd, lw=1, capsize=10)
                for caps in listcaps:
                    caps.set_markeredgewidth(1)
            else:
                axis.plot(xdat, ydat, label=gdat.legdsamp, alpha=0.5)
    
    # reference histogram
    try:
        if not omittrue and gdat.allwrefr:
            for q in gdat.indxrefr:
                if strgydat.startswith('psfn'):
                    ydattemp = getattr(gdat, 'psfnexpr')
                    print 'hey'
                    print 'strgydat'
                    print strgydat
                    print 'ydattemp'
                    summgene(ydattemp)
                    print 'indxydat'
                    print indxydat
                    print ''
                    #ydattemp = ydattemp[indxydat]
                    print 'ydattemp'
                    print ydattemp
                else:
                    ydattemp = getattr(gdat, 'refr' + strgydat[:-4] + 'ref%d' % q)
                ydat = ydattemp * factydat
                
                if indxydat != None:
                    ydat = ydat[indxydat]
                
                if strgydat.startswith('psfn'):
                    print 'ydat'
                    print ydat
                    print

                if histodim:
                    axis.bar(xdattemp, ydat, deltxdat, color=gdat.listcolrrefr[q], label=gdat.legdrefr[q], alpha=gdat.alphmrkr)
                else:
                    
                    if strgydat.startswith('psfn'):
                        print 'xdat'
                        print xdat
                        print 'ydat'
                        print ydat
                    axis.plot(xdat, ydat, color=gdat.listcolrrefr[q], label=gdat.legdrefr[q], alpha=gdat.alphmrkr)
    except:
        print 'Reference overplotting failed for %s' % strgydat
    # axis scales
    if scalxdat == 'logt':
        axis.set_xscale('log')
    if scalydat == 'logt':
        if where(ydat > 0.)[0].size > 0:
            axis.set_yscale('log')
    
    # axis labels
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)

    # superimpose prior on the feature
    liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
    liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
    
    if strgydat.startswith('hist') and strgydat != 'histdefl' and strgydat != 'histdeflelem':
        print 'strgydat'
        print strgydat
        print
        if strgydat[-4:-1] == 'pop':
            if strgydat[4:-4] in liststrgfeatprio[int(strgydat[-1])]:
                xdatprio = getattr(gdat, strgmodl + strgxdat + 'prio') * factxdat
                if strg == 'true' or strg == 'post':
                    gdattemp = gdat
                else:
                    gdattemp = gdatmodi
    
                if gdat.datatype == 'mock' and not omittrue:
                    truexdatprio = getattr(gdat, 'true' + strgxdat + 'prio') * factxdat
                    trueydatsupr = getattr(gdat, 'true' + strgydat + 'prio') * factydat
                    axis.plot(truexdatprio, trueydatsupr, ls='-', alpha=gdat.alphmrkr, color='g')
                
                if strg != 'true':
                    if strg == 'post':
                        ydatsupr = getattr(gdattemp, 'medi' + strgydat + 'prio') * factydat
                        yerrsupr = getattr(gdattemp, 'errr' + strgydat + 'prio') * factydat
                        labl = gdat.legdsampdist + ' hyper-distribution'
                        tdpy.util.plot_braz(axis, xdatprio, ydatsupr, yerr=yerrsupr, lcol='lightgrey', dcol='grey', labl=labl)
                    else:
                        ydatsupr = getattr(gdattemp, strg + strgydat + 'prio') * factydat
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
        axis.set_ylim([limtydat[0] * factydat, limtydat[1] * factydat])
    else:
        axis.set_ylim([amin(ydat), amax(ydat)])
    
    if ydattype != 'totl':
        strgydat += ydattype

    # temp
    try:
        make_legd(axis, offs=offslegd)
    except:
        pass

    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, strgydat, nameinte=nameinte)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatassc(gdat, gdatmodi, strg, q, l, strgfeat, indxregiplot, plotdiff=False):
    
    figr, axis = plt.subplots(1, 1, figsize=(gdat.plotsize, gdat.plotsize))

    # prepare data to be plotted
    xdat = copy(getattr(gdat, 'refr' + strgfeat)[indxregiplot][q][0, :])
    xerr = tdpy.util.retr_errrvarb(getattr(gdat, 'refr' + strgfeat)[indxregiplot][q])
    
    minmplot = getattr(gdat, 'minm' + strgfeat)
    maxmplot = getattr(gdat, 'maxm' + strgfeat)
    binsplot = getattr(gdat, 'bins' + strgfeat)
    factplot = getattr(gdat, 'fact' + strgfeat + 'plot')
    
    ydat = retr_fromgdat(gdat, gdatmodi, strg, strgfeat + 'asscreg%dref%dpop%d' % (indxregiplot, q, l))
    yerr = zeros((2, ydat.size))
    
    if strg == 'post':
        yerr = retr_fromgdat(gdat, gdatmodi, strg, strgfeat + 'asscreg%dref%dpop%d' % (indxregiplot, q, l), mometype='errr')
    
    xdat *= factplot
    xerr *= factplot
    if plotdiff:
        ydat = 100. * (ydat - xdat) / xdat
    else: 
        ydat *= factplot
        if strg == 'post':
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
    
    lablxdat = getattr(gdat, 'labl' + strgfeat + 'refr')
    lablydat = getattr(gdat, 'labl' + strgfeat + 'samp')
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)
    if indx.size > 0 and getattr(gdat, 'scal' + strgfeat + 'plot') == 'logt':
        if not plotdiff:
            axis.set_yscale('log')
        axis.set_xscale('log')
    if plotdiff:
        limsydat = array([-100., 100.])
    else:
        limsydat = factplot * array([minmplot, maxmplot])
    axis.set_ylim(limsydat)
    axis.set_xlim([factplot * minmplot, factplot * maxmplot])
   
    plt.tight_layout()
    if plotdiff:
        strgtype = 'diff'
    else:
        strgtype = ''
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatassc' + strgfeat + '%sref%dpop%d' % (strgtype, q, l), nameinte='assc/')
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def make_legd(axis, offs=None, loca=1, numbcols=1):
    
    legd = axis.legend(fancybox=True, frameon=True, bbox_to_anchor=offs, bbox_transform=axis.transAxes, ncol=numbcols, loc=loca, labelspacing=1, handlelength=2)
    
    legd.get_frame().set_fill(True)
    legd.get_frame().set_facecolor('white')


def plot_scatcntp(gdat, gdatmodi, strg, indxregiplot, indxenerplot, indxevttplot):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    ydat = retr_fromgdat(gdat, gdatmodi, strg, 'cntpmodl')[indxregiplot, indxenerplot, :, indxevttplot]
    axis.scatter(gdat.cntpdata[indxregiplot, indxenerplot, :, indxevttplot], ydat, alpha=gdat.alphmrkr)
    gdat.limtcntpdata = [gdat.binscntpdata[0], gdat.binscntpdata[-1]]
    axis.set_xlim(gdat.limtcntpdata)
    axis.set_ylim(gdat.limtcntpdata)
    axis.set_ylabel('$k^{modl}$')
    axis.set_xlabel('$k^{data}$')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.tight_layout()
    path = retr_plotpath(gdat, gdatmodi, strg, 'scatcntpreg%dene%devt%d' % (indxregiplot, indxenerplot, indxevttplot))
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 40
    numbfluxprox = len(gdat.indxpixlprox)
    bins = empty((gdat.numbprox, numbbins + 1))
    indxpixlproxsize = empty((numbfluxprox, gdat.numbpixlfull))
    for h in range(gdat.numbprox):
        for j in gdat.indxpixlfull:
            try:
                indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
            except:
                indxpixlproxsize[h, j] = gdat.numbpixlfull
        bins[h, :] = logspace(log10(amin(indxpixlproxsize[h, :])), log10(amax(indxpixlproxsize[h, :])), numbbins + 1)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbprox):
        axis.hist(indxpixlproxsize[h, :], bins=bins[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphmrkr)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixlfull, label='ROI', ls='--')
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
            
            # superimpose reference elements

            if gdat.allwrefr:
                refrsign = getattr(gdat, 'refr' + gdat.namefeatsign)
                for q in gdat.indxrefr:
                    if len(refrsign[q]) > 0:
                        if strgfeat != None:
                            refrfeat = getattr(gdat, 'refr' + strgfeat)[q][0, :] 
                            indxelem = where((bins[indxlowr] < refrfeat) & (refrfeat < bins[indxuppr]))[0]
                        else:
                            indxelem = arange(gdat.refrnumbelem[q])
                        
                        mrkrsize = retr_mrkrsize(gdat, refrsign[q][0, indxelem], gdat.namefeatsign)
                        axis.scatter(gdat.anglfact * gdat.refrlgal[q][0, indxelem], gdat.anglfact * gdat.refrbgal[q][0, indxelem], \
                                                                                                    s=mrkrsize, alpha=gdat.alphmrkr, marker='*', lw=2, color=gdat.listcolrrefr[q])

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
    
    
def plot_intr(gdat):
    
    if gdat.verbtype > 0:
        print 'Making PCAT introductory plots...'

    plot_grap(plottype='lght0000', verbtype=1)
    plot_grap(plottype='lght0001', verbtype=1)
    plot_grap(plottype='lght0002', verbtype=1)
    plot_grap(plottype='lght0003', verbtype=1)
    plot_grap(plottype='lens0000', verbtype=1)
    plot_grap(plottype='lens0001', verbtype=1)
    
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
        psfntemp = copy(gdat.psfnexpr[0, :, 0, 0])
    else:
        psfntemp = copy(gdat.psfnexpr[0, :, 0])
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
    axis.set_ylabel(gdat.lablsbrttotl)
    
    limt = gdat.specfraceval * amax(gdat.binsprox[0] * psfntemp)

    if limt != 0.:
        axis.axhline(limt, color='red', ls=':', label='Flux floor')
    
    #axis.set_ylim([1e-3 * limt, None])
    #maxmangltemp = interp(1e-1 * limt, gdat.binsprox[k] * psfntemp[::-1], gdat.binsangl[::-1] * gdat.anglfact)
    #axis.set_xlim([None, maxmangltemp])
    
    make_legd(axis)

    plt.tight_layout()
    figr.savefig(gdat.pathinit + 'eval.pdf')
    plt.close(figr)


def plot_mosa(gdat, indxregiplot):

    # empty global object
    gdatmodi = tdpy.util.gdatstrt()
    gdatmodi.thischro = zeros(gdat.numbchro)

    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamptotl:
        indxsampmosa = choice(gdat.indxsamptotl, size=numbsampmosa, replace=False)
        for l in gdat.fittindxpopl:
            for d in gdat.indxregi:
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
                                
                                imag = retr_imag(gdat, axis, gdat.cntpdata, '', 'cntpdata', indxregiplot, i, m)
                                supr_fram(gdat, gdatmodi, 'this', axis, indxregiplot, l)
                        
                        if gdat.enerbins:
                            plt.figtext(0.5, 0.93, gdat.strgener[i], ha='center', va='center')
                        axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
                        cbar = figr.colorbar(imag, cax=axiscomm)
                        cbar.set_ticks(gdat.tickcntpdata)
                        cbar.set_ticklabels(gdat.lablcntpdata)
                        plt.subplots_adjust(left=0.1, top=.91, hspace=0.03, wspace=0.1, bottom=0.09)
                        if l == 1:
                            strg = ''
                        else:
                            strg = 'pop%d' % l
                        pathfinl = getattr(gdat, 'path' + gdat.namesampdist + 'finl')
                        if m == None:
                            path = pathfinl + 'mosa' + strg + 'reg%dene%dA.pdf' % (indxregiplot, gdat.indxenerincl[i])
                        else:
                            path = pathfinl + 'mosa' + strg + 'reg%dene%devtt%d.pdf' % (indxregiplot, gdat.indxenerincl[i], gdat.indxevttincl[m])
                        figr.savefig(path)
                        plt.close(figr)
    else:
        if gdat.verbtype > 0:
            print 'Skipping the mosaic plot...'


def plot_grap(plottype, verbtype=0):
        
    figr, axis = plt.subplots(figsize=(6, 6))

    grap = nx.DiGraph()
    if plottype == 'lght0000':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']

    if plottype == 'lght0001':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta', 'black']

    if plottype == 'lght0002':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'olive', 'magenta', \
                                                                                                    'magenta', 'magenta', 'magenta', 'magenta', 'black']
    if plottype == 'lght0003':
        listcolr = ['black', 'black', 'black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', \
                                                                                                    'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']
    
    if plottype == 'lens0000':
        listcolr = ['olive', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', 'magenta']

    if plottype == 'lens0001':
        listcolr = ['olive', 'black', 'olive', 'olive', 'olive', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', \
                                                                                                                                        'magenta', 'magenta', 'magenta']

    grap.add_edges_from([ \
                         ('meanelem', 'numbelem'), \
                         ('modl','data'), \
                         ('psfp', 'modl'), \
                         ('bacp', 'modl'), \
                         ('lgal','modl'), \
                         ('bgal','modl'), \
                         ('numbelem','lgal'), \
                         ('numbelem','bgal'), \
                        ])
    
    if plottype.startswith('lght'):
        grap.add_edges_from([ \
                             ('ampldistslop', 'ampl'), \
                             ('ampl', 'modl'), \
                             ('numbelem','ampl'), \
                             ('numbelem', 'sind'), \
                             ('sind','modl'), \
                            ])
    
    if plottype.startswith('lens'):
        grap.add_edges_from([ \
                             ('lenp', 'modl'), \
                             ('defsdistslop', 'defs'), \
                             ('defs', 'modl'), \
                             ('numbelem','defs'), \
                            ])
    
    if plottype == 'lens0001':
        grap.add_edges_from([ \
                             ('asca', 'modl'), \
                             ('numbelem','asca'), \
                             ('acut', 'modl'), \
                             ('numbelem','acut'), \
                            ])
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('sinddistmean', 'sind'), \
                            ])
    
    if plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('numbelem', 'expc'), \
                             ('expc', 'modl'), \
                            ])
    
    if plottype == 'lght0003':
        grap.add_edges_from([ \
                             ('spatdistcons', 'lgal'), \
                             ('spatdistcons', 'bgal'), \
                            ])
        
    labl = {}
    if plottype.startswith('lens'):
        nameelem = r'\rm{sub}'
    else:
        nameelem = r'\rm{pts}'
    if plottype.startswith('lght') and (plottype == 'lght0001' or plottype == 'lght0002'):
        labl['numbelem'] = r'$\vec{N}_{%s}$' % nameelem
        labl['meanelem'] = r'$\vec{\mu}_{%s}$' % nameelem
    else:
        labl['numbelem'] = '$N_{%s}$' % nameelem
        labl['meanelem'] = r'$\mu_{%s}$' % nameelem
    
    if plottype.startswith('lght'):
        if plottype == 'lght0000' or plottype == 'lght0003':
            labl['ampldistslop'] = r'$\alpha$'
        else:
            labl['ampldistslop'] = r'$\vec{\alpha}$'
    if plottype.startswith('lens'):
        labl['defsdistslop'] = r'$\beta$'
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        labl['sinddistmean'] = r'$\vec{\beta}$'
    
    if plottype == 'lght0003':
        labl['spatdistcons'] = r'$\gamma$'
    if plottype.startswith('lens'):
        labl['lenp'] = r'$\vec{\chi}$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['bacp'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{\theta_1}$'
    labl['bgal'] = r'$\vec{\theta_2}$'
    if plottype.startswith('lght'):
        labl['sind'] = r'$\vec{s}$'
        labl['ampl'] = r'$\vec{f}$'
    else:
        labl['defs'] = r'$\vec{\alpha_{\rm{s}}}$'
    if plottype == 'lens0001':
        labl['asca'] = r'$\vec{\theta_{\rm{s}}}$'
        labl['acut'] = r'$\vec{\theta_{\rm{c}}}$'
        
    if plottype == 'lght0002':
        labl['expc'] = r'$\vec{E_{\rm{c}}}$'
    labl['modl'] = r'$\mathcal{M}$'
    labl['data'] = r'$\mathcal{D}$'
    
    posi = nx.circular_layout(grap)
    posi['sinddistmean'] = array([0.4, 0.15])
    if plottype == 'lght0003':
        posi['spatdistcons'] = array([-0.2, 0.15])
    if plottype.startswith('lght'):
        posi['numbelem'] = array([0., 0.075])
        posi['meanelem'] = array([0., 0.15])
        posi['ampldistslop'] = array([0.2, 0.15])
    if plottype.startswith('lens'):
        posi['numbelem'] = array([-0.1, 0.075])
        posi['meanelem'] = array([-0.1, 0.15])
        posi['defsdistslop'] = array([0.1, 0.15])
    
    if plottype.startswith('lght'):
        if plottype == 'lght0002':
            posi['psfp'] = array([0.9, -0.0])
            posi['bacp'] = array([0.7, -0.0])
        else:
            posi['psfp'] = array([0.7, -0.0])
            posi['bacp'] = array([0.5, -0.0])
    if plottype == 'lens0000':
        posi['psfp'] = array([0.3, -0.0])
        posi['bacp'] = array([0.5, -0.0])
        posi['lenp'] = array([0.7, -0.0])
    if plottype == 'lens0001':
        posi['psfp'] = array([0.7, -0.0])
        posi['bacp'] = array([0.9, -0.0])
        posi['lenp'] = array([1.1, -0.0])
    posi['lgal'] = array([-0.3, -0.0])
    posi['bgal'] = array([-0.1, -0.0])
    if plottype.startswith('lght'):
        posi['ampl'] = array([0.1, -0.0])
        posi['sind'] = array([0.3, -0.0])
    if plottype == 'lght0002':
        posi['expc'] = array([0.5, -0.0])

    if plottype.startswith('lens'):
        posi['defs'] = array([0.1, -0.0])
    if plottype == 'lens0001':
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
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbelem'], node_color='b', node_size=size)
    if plottype.startswith('lght'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'ampldistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind'], node_color='g', node_size=size)
    if plottype == 'lght0001' or plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['sinddistmean'], node_color='r', node_size=size)
    if plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['expc'], node_color='g', node_size=size)
    if plottype == 'lght0003':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['spatdistcons'], node_color='r', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['psfp', 'bacp'], node_color='y', node_size=size)
    if plottype.startswith('lens'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'defsdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lenp'], node_color='y', node_size=size)
    if plottype == 'lens0000':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs'], node_color='g', node_size=size)
    if plottype == 'lens0001':
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
        
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    name = 'histcntpdatareg%dene%devt%d' % (d, i, m)
                    plot_gene(gdat, None, '', name, 'meancntpdata', scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbpixl])
        if gdat.numbpixl > 1 and gdat.evalcirc != 'full':
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

    for d in gdat.indxregi:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                
                # temp
                if False and gdat.pixltype == 'cart' and (gdat.elemtype == 'lght' or gdat.elemtype == 'clus'):
                    figr, axis, path = init_figr(gdat, None, 'cntpdatapeak', '', d, i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.cntpdata, '', 'cntpdata', d, i, m)
                    make_cbar(gdat, axis, imag, i, tick=gdat.tickcntpdata, labl=gdat.lablcntpdata)
                    axis.scatter(gdat.anglfact * gdat.meanlgalcart[gdat.indxxdatmaxm], gdat.anglfact * gdat.meanbgalcart[gdat.indxydatmaxm], alpha=0.6, s=20, edgecolor='none')
                    
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
        
                if gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
                    figr, axis, path = init_figr(gdat, None, 'cntpmodlraww', 'true', d, i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.truecntpmodlraww, '', 'cntpdata', d, i, m, tdim=True)
                    make_cbar(gdat, axis, imag, 0, tick=gdat.tickcntpdata, labl=gdat.lablcntpdata)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)

        if gdat.correxpo:
            gdat.lablnumbpixl = r'$N_{\rm{pix}}$'
            gdat.limtexpo = [gdat.minmexpo, gdat.maxmexpo]
            if gdat.enerbins:
                path = gdat.pathinit + 'expototlmeanreg%d.pdf' % d
                tdpy.util.plot_gene(path, gdat.meanener, gdat.expototlmean[d, :], scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, \
                                                                                                lablydat=gdat.lablexpototl, limtydat=gdat.limtexpo)
            
            gdat.histexpototl = empty((gdat.numbener, gdat.numbtickcbar - 1))
            for i in gdat.indxener:
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.hist(gdat.expototl[d, i, :], gdat.binsexpo)
                axis.set_xlabel(gdat.lablexpototl)
                axis.set_ylabel(gdat.lablnumbpixl)
                axis.set_xscale('log')
                axis.set_yscale('log')
                plt.tight_layout()
                path = gdat.pathinit + 'histexporeg%dene%d.pdf' % (d, i)
                figr.savefig(path)
                plt.close(figr)

                for m in gdat.indxevtt:
                    figr, axis, path = init_figr(gdat, None, 'expo', '', d, i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.expo, '', 'expo', d, i, m)
                    make_cbar(gdat, axis, imag, i)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
                

def plot_defl(gdat, gdatmodi, strg, indxregiplot, strgvarb='defl', strgcomp='', indxdefl=None, indxpoplplot=-1, multfact=1., indxenerplot=None, indxevttplot=None):

    if indxdefl != None:
        strgvarb += 'sing'
    strgvarb = strgvarb + strgcomp
    
    defl = retr_fromgdat(gdat, gdatmodi, strg, strgvarb)
    
    defl *= multfact
   
    print 'defl'
    summgene(defl)
    print 

    if indxenerplot != None:
        defl = defl[indxregiplot, indxenerplot, :, indxevttplot, ...]
    else:
        defl = defl[indxregiplot, ...]
    print 'defl'
    summgene(defl)
    print 

    if indxdefl != None:
        defl = defl[..., indxdefl]
        strgvarb += '%04d' % indxdefl
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))

    figr, axis, path = init_figr(gdat, gdatmodi, strgvarb, strg, indxregiplot, indxenerplot, indxevttplot, indxpoplplot)
    make_catllabl(gdat, strg, axis)
    draw_frambndr(gdat, axis)
  
    defllgal = defl[:, :, 0]
    deflbgal = defl[:, :, 1]
    fact = 4
    axis.imshow(zeros((10, 10)))
    ptch = axis.quiver(gdat.anglfact * gdat.lgalgridcart[::fact, ::fact], gdat.anglfact * gdat.bgalgridcart[::fact, ::fact], \
                  gdat.fluxfactplot * defllgal[::fact, ::fact], gdat.fluxfactplot * deflbgal[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strg, axis, indxregiplot)
    #plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strg, strgvarb, indxregiplot, indxenerplot=None, indxevttplot=-1, strgcbar=None, \
                                                                tdim=False, indxpoplplot=-1, mometype='medi', intreval=False, indxmaps=None):
    
    if strgcbar == None:
        strgcbar = strgvarb
  
    # construct the string for the map
    if strgvarb == 'cntpdata':
        strgplot = strgvarb
    else:
        if strg == 'post':
            strgtemp = mometype
        else:
            strgtemp = ''
        strgplot = strgtemp + strgvarb
    if indxmaps != None:
        strgplot += '%04d' % indxmaps

    figr, axis, path = init_figr(gdat, gdatmodi, strgplot, strg, indxregiplot, indxenerplot, indxevttplot, indxpoplplot, intreval=intreval)
   
    maps = retr_fromgdat(gdat, gdatmodi, strg, strgvarb, mometype=mometype, indxlist=indxmaps)
    print 'maps'
    summgene(maps)
    imag = retr_imag(gdat, axis, maps, strg, strgcbar, indxregiplot, indxenerplot, indxevttplot, tdim=tdim)
    tick = getattr(gdat, 'tick' + strgcbar) 
    labl = getattr(gdat, 'labl' + strgcbar) 

    make_cbar(gdat, axis, imag, tick=tick, labl=labl)
    make_catllabl(gdat, strg, axis)
    if gdat.suprelem:
        supr_fram(gdat, gdatmodi, strg, axis, indxregiplot, indxpoplplot)

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
            retr_imag(gdat, axis, maps, strg, strgcbar, indxregiplot, indxenerplot, indxevttplot, tdim=tdim, imag=imag)
            for ptch in axis.get_children():
                if isinstance(ptch, mpl.patches.Circle) or isinstance(ptch, mpl.collections.PathCollection):#isinstance(ptch, mpl.lines.Line2D):
                    ptch.remove()
            supr_fram(gdat, gdatmodi, strg, axis, indxregiplot, indxpoplplot)
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
    
