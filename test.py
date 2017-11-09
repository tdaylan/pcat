def plot_post(gdat=None, pathpcat=None, verbtype=1, prio=False):
    
    if gdat.verbtype > 0:
        print 'Producing postprocessing plots...'

    timetotlinit = gdat.functime()
    
    gdat.strgbest = 'ML'
    
    for name, valu in gdat.__dict__.iteritems():
        if name.startswith('list'):
            print name
            print

    # prior components
    if gdat.diagmode:
        if gdat.verbtype > 0:
            print 'Plotting the prior distribution...'

        for n in gdat.indxproptype:
            pathpostlpri = getattr(gdat, 'path' + gdat.namesampdist + 'finllpri')
            #path = pathpostlpri + 'deltlliktotl' + gdat.nameproptype[n]
            #tdpy.mcmc.plot_trac(path, gdat.listdeltlliktotl[gdat.listindxsamptotl[n]], r'$\log \Delta P(D|M)$', logthist=True)
            #for k in gdat.fittindxlpri:
            #    path = pathpostlpri + 'deltlpritotl%04d' % k + gdat.nameproptype[n]
            #    tdpy.mcmc.plot_trac(path, gdat.listdeltlpritotl[gdat.listindxsamptotl[n], k], r'$\log \Delta P_{%d}(M)$' % k, logthist=True)
            if gdat.nameproptype[n] == 'brth' or gdat.nameproptype[n] == 'deth' or gdat.nameproptype[n] == 'splt' or gdat.nameproptype[n] == 'merg':
                path = pathpostlpri + 'lpautotl%04d' + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlpautotl[gdat.listindxsamptotl[n]], '', logthist=True)
                for k in gdat.fittindxlpau:
                    path = pathpostlpri + 'lpau%04d' % k + gdat.nameproptype[n]
                    tdpy.mcmc.plot_trac(path, gdat.listlpau[gdat.listindxsamptotl[n], k], r'$\log u_{%d}$' % k, logthist=True)
                
            if gdat.nameproptype[n] == 'splt' or gdat.nameproptype[n] == 'merg':
                path = pathpostlpri + 'lrpp' + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listlrpp[gdat.listindxsamptotl[n]], r'$\log \alpha_p$', logthist=True)

                path = pathpostlpri + 'ljcb' + gdat.nameproptype[n]
                tdpy.mcmc.plot_trac(path, gdat.listljcb[gdat.listindxsamptotl[n]], r'$\log \alpha_c$', logthist=True)

    # Gelman-Rubin test
    pathdiag = getattr(gdat, 'path' + gdat.namesampdist + 'finldiag')
    if gdat.numbproc > 1:
        for d in gdat.indxregi:
            if isfinite(gdat.gmrbstat[d]).all():
                if gdat.verbtype > 0:
                    print 'Gelman-Rubin TS...'
    
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                minm = min(amin(gdat.gmrbstat[d]), amin(gdat.gmrbfixp))
                maxm = max(amax(gdat.gmrbstat[d]), amax(gdat.gmrbfixp))
                bins = linspace(minm, maxm, 40)
                axis.hist(gdat.gmrbstat[d].flatten(), bins=bins, label='Data proj.')
                axis.hist(gdat.gmrbfixp, bins=bins, label='Fixed dim.')
                axis.set_xlabel('PSRF')
                axis.set_ylabel('$N_{stat}$')
                plt.tight_layout()
                figr.savefig(pathdiag + 'gmrbhist.pdf')
                plt.close(figr)
                
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.plot(gdat.fittindxfixp, gdat.gmrbfixp)
                axis.set_xticklabels(gdat.fittlablfixp)
                axis.set_ylabel('PSRF')
                plt.tight_layout()
                figr.savefig(pathdiag + 'gmrbfixp.pdf')
                plt.close(figr)
                
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        maps = gdat.gmrbstat[d][i, :, m]
                        path = pathdiag + 'gmrbdataene%devt%d.pdf' % (i, m)
                        tdpy.util.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                                minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                                minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
            else:
                print 'Inappropriate Gelman-Rubin test statistics encountered for region %d.' % d

    # plot autocorrelation
    if gdat.verbtype > 0:
        print 'Autocorrelation...'
    for d in gdat.indxregi:
        tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrcntp[d][0, 0, 0, 0, :], gdat.timeatcrcntp[d][0, 0, 0, 0], strgextn='cntp')
    tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrpara[0, 0, :], gdat.timeatcrpara[0, 0], strgextn='para')
    
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
            proc_samp(gdat, gdatmodi, 'this', 'fitt')
            plot_samp(gdat, gdatmodi, 'this', 'fitt')

    # plot split and merge diagnostics
    if gdat.fittnumbtrap > 0 and gdat.probspmr > 0.:
        if gdat.verbtype > 0:
            print 'Split and merge related plots...'
    
        indxsampsplttotl = where(gdat.listindxproptype == gdat.indxproptypesplt)[0]
        indxsampsplt = intersect1d(where(gdat.listindxproptype == gdat.indxproptypesplt)[0], where(gdat.listaccpprop)[0])
        indxsampmergtotl = where(gdat.listindxproptype == gdat.indxproptypemerg)[0]
        indxsampmerg = intersect1d(where(gdat.listindxproptype == gdat.indxproptypemerg)[0], where(gdat.listaccpprop)[0])
        indxsampspmrtotl = concatenate((indxsampsplttotl, indxsampmergtotl))
        if indxsampspmrtotl.size > 0:
            for l in gdat.fittindxpopl:
                ## labels and names
                if gdat.fittelemtype[l] == 'lghtline':
                    listlabl = ['$u_e$', '$u_f$', r'$\log\alpha_j$', r'$\log\alpha_p$']
                    listname = ['elinauxi', 'fracauxi', 'ljcb', 'lrpp']
                    listvarb = [gdat.listauxipara[:, 0], gdat.listauxipara[:, 1], gdat.listljcb, gdat.listlrpp]
                else:
                    listlabl = ['$u_x$', r'$u_y$', r'$u_f$', r'$\log\alpha_j$', r'$\log\alpha_p$']
                    listname = ['lgalauxi', 'bgalauxi', 'fracauxi', 'ljcb', 'lrpp']
                    listvarb = [gdat.anglfact * gdat.listauxipara[:, 0], gdat.anglfact * gdat.listauxipara[:, 1], gdat.listauxipara[:, 2], gdat.listljcb, gdat.listlrpp]
                
                ## variables
                print 'gdat.listauxipara'
                summgene(gdat.listauxipara)
                print 'gdat.listauxipara[:, 0]'
                summgene(gdat.listauxipara[:, 0])
                print 'gdat.listauxipara[:, 1]'
                summgene(gdat.listauxipara[:, 1])
                if gdat.fittelemtype[l] != 'lghtline':
                    print 'gdat.listauxipara[:, 2]'
                    summgene(gdat.listauxipara[:, 2])
                
                for k in range(len(listlabl)):
                    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                    maxm = amax(listvarb[k][indxsampspmrtotl])
                    minm = amin(listvarb[k][indxsampspmrtotl])
                    bins = linspace(minm, maxm, 40)
                  
                    axis.hist(listvarb[k][indxsampsplttotl], bins=bins, label='Split', alpha=gdat.alphmrkr, color='c')
                    axis.hist(listvarb[k][indxsampmergtotl], bins=bins, label='Merge', alpha=gdat.alphmrkr, color='y')
                    axis.hist(listvarb[k][indxsampsplt], bins=bins, label='Split, Accepted', edgecolor='c', facecolor='none')
                    axis.hist(listvarb[k][indxsampmerg], bins=bins, label='Merge, Accepted', edgecolor='y', facecolor='none')
                    
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
    plot_samp(gdat, None, 'post', 'fitt')
    #proc_samp(gdat, None, 'mlik', 'fitt')
    
    if gdat.verbtype > 0:
        print 'A mosaic of samples...'
    
    ## mosaic of images of posterior catalogs
    if gdat.numbpixl > 1:
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
        if numbtrapplot > 0:
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
    
            # temp -- this comes up in lens_syst
            try:
                tdpy.mcmc.plot_grid(pathjoin, listjoin, [labltotl, labltotlseco], scalpara=[scal, scalseco], truepara=[truepara, trueparaseco], \
                                                                                                                                  join=True, varbdraw=[mlik, mlikseco])
            except:
                pass

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
    if gdat.proppsfp and gdat.numbpixl > 1 and gdat.fittpsfnevaltype != 'none':
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
    if gdat.fittnumbtrap > 0 and gdat.numbpixl > 1:
        liststrgbins = ['quad', 'full']
        for d in gdat.indxregi:
            for l in gdat.fittindxpopl:
                plot_posthistlgalbgalelemstkd(gdat, d, l, 'cumu')
                for strgbins in liststrgbins:
                    for strgfeatsign in gdat.fittliststrgfeatsign[l]:
                        plot_posthistlgalbgalelemstkd(gdat, d, l, strgbins, strgfeatsign)

    if gdat.verbtype > 0:
        print 'Prior and likelihood...'
    
    for strgpdfn in ['lpri', 'llik']:

        if strgpdfn == 'lpri':
            labltemp = '\ln P(M)'
        if strgpdfn == 'llik':
            labltemp = '\ln P(D|M)'
        labl = r'$%s$' % labltemp

        setattr(gdat, 'list' + strgpdfn + 'flat', getattr(gdat, 'list' + strgpdfn + 'totl').flatten())
        
        path = getattr(gdat, 'path' + gdat.namesampdist + 'finl') + strgpdfn
        titl = r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (gdat.info, gdat.levi)
        tdpy.mcmc.plot_hist(path, getattr(gdat, 'list' + strgpdfn + 'flat'), labl, titl)
        varbdraw = []
        labldraw = []
        colrdraw = []
        if gdat.datatype == 'mock':
            varbdraw += [getattr(gdat, 'true' + strgpdfn + 'totl')]
            labldraw += ['True model']
            colrdraw += [gdat.refrcolr]
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'list' + strgpdfn + 'flat'), labl, varbdraw=varbdraw, labldraw=labldraw, colrdraw=colrdraw)
    
        if strgpdfn == 'llik':
            labldelt = r'$\Delta%s$' % labltemp
            setattr(gdat, 'listdelt' + strgpdfn + 'totlflat', getattr(gdat, 'listdelt' + strgpdfn + 'totl').flatten())
            pathbase = getattr(gdat, 'path' + gdat.namesampdist + 'finldelt%s' % strgpdfn)
            path = pathbase + 'delt%s' % strgpdfn
    
            tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strgpdfn + 'totlflat'), labldelt)
            
            if gdat.numbproc > 1:
                for k in gdat.indxproc:
                    path = pathbase + 'delt%s_proc%04d' % (strgpdfn, k)
                    tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strgpdfn + 'totl')[:, k], labldelt, titl='Process %d' % k)
            for n in gdat.indxproptype:
                path = pathbase + 'delt%s_%s' % (strgpdfn, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strgpdfn + 'totlflat')[gdat.listindxsamptotl[n]], labldelt, titl=gdat.nameproptype[n])
                path = pathbase + 'delt%s_%s_accp' % (strgpdfn, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strgpdfn + 'totlflat')[gdat.listindxsamptotlaccp[n]], labldelt, titl=gdat.nameproptype[n] + ', Accepted')
                path = pathbase + 'delt%s_%s_reje' % (strgpdfn, gdat.nameproptype[n])
                tdpy.mcmc.plot_trac(path, getattr(gdat, 'listdelt' + strgpdfn + 'totlflat')[gdat.listindxsamptotlreje[n]], labldelt, titl=gdat.nameproptype[n] + ', Rejected')
        
    # plot resident memory
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxswep, mean(gdat.listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(pathdiag + 'memoresi.pdf')
    plt.close(figr)

    timetotlfinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Plots and animations are produced in %.3g seconds.' % (timetotlfinl - timetotlinit)



