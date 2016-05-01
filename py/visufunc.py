
# coding: utf-8

# In[ ]:

def plot_look():

    ypixlsize = zeros(npixl)
    figr, axis = plt.subplots(figsize=(10, 6))
    for h in range(numbspecprox):
        for j in ipixl:
            ypixlsize[j] = ypixl[h][j].size
        binspixlsize = logspace(log10(amin(ypixlsize)), log10(amax(ypixlsize)), 100)
        hist, bins = histogram(ypixlsize, binspixlsize)
        mean = sqrt(bins[:-1] * bins[1:])
        axis.loglog(mean, hist, label='Flux bin %d' % h)
    axis.set_title("Number of pixels in the pixel lookup tables")
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    axis.legend()
    
    plt.savefig(plotpath + 'look.png')
    plt.close()


# In[ ]:

def plot_datacntshist():

    
    figr, axgr = plt.subplots(numbevtt, numbener, figsize=(7 * numbener, 7 * numbevtt))
    if numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            datacntstemp = datacnts[i, :, m]
            
            maxmdatacntstemp = amax(datacntstemp)
            minmdatacntstemp = amin(datacntstemp[where(datacntstemp > 0.)])
            binscntstemp = linspace(minmdatacntstemp, maxmdatacntstemp, 20)
            meancntstemp = (binscntstemp[1:] + binscntstemp[0:-1]) / 2.
            diffcntstemp = binscntstemp[1:] - binscntstemp[0:-1]
            
            datacntshist = axis.hist(datacntstemp, binscntstemp, color='b')[0]

            init = [meancntstemp[argmax(datacntshist)], 1.]
            
            global datacntshistglob
            datacntshistglob = diffcntstemp * sum(datacntshist)
    
            axis.set_xlim([amin(binscntstemp), amax(binscntstemp)])
            if m == numbevtt - 1:
                axis.set_xlabel(r'$k$')
            #axis.set_xscale('log')
            if m == 0:
                axis.set_title(binsenerstrg[i])
            if i == 0 and exprtype == 'ferm':
                axis.set_ylabel(evttstrg[m])
        
    figr.subplots_adjust(wspace=0.3, hspace=0.2)
    plt.savefig(plotpath + 'datacntshist' + rtag + '.png')
    plt.close(figr)


# In[ ]:

def plot_intr():
    
    with plt.xkcd():

        figr, axis = plt.subplots(figsize=(14, 6))

        catl = arange(80)
        probcatl = pss.pmf(catl, 30.) + 0.5 * pss.pmf(catl, 60.)
        axis.plot(catl, probcatl)
        axis.set_xticks([10, 30, 60])
        axis.set_xticklabels(["Crackpot's Catalog", "Best-fit catalog", "Not-so-best-fit catalog"])
        axis.set_yticks([])
        axis.set_title("Exploring the catalog space with Probabilistic cataloging")
        axis.set_xlabel('Catalog index')
        axis.set_ylabel("Probability")
        plt.savefig(os.environ["PNTS_TRAN_DATA_PATH"] + '/png/talkintr.png')
        plt.close()  


# In[ ]:

def plot_heal(heal, rofi=True, titl=''):
    
    if rofi:
        healtemp = copy(heal)
        heal = zeros(npixlheal)
        heal[jpixlrofi] = healtemp

    cart = tdpy_util.retr_cart(heal, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
    imag = plt.imshow(cart, origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, fraction=0.05)
    plt.title(titl)
    plt.show()


# In[ ]:

def plot_3fgl_thrs():


    expoheal = sum(sum(expoheal, 2)[1:3, :], axis=0)


    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/catl/3fgl_thrs.fits'
    fluxthrs = pf.getdata(path, 0)

    bgalfgl3 = linspace(-90., 90., 481)
    lgalfgl3 = linspace(-180., 180., 960)

    bgalexpo = linspace(-90., 90., 400)
    lgalexpo = linspace(-180., 180., 800)

    fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)

    cntsthrs = fluxthrs * expo

    jbgal = where(abs(bgalexpo) < 10.)[0]
    jlgal = where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    

    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)

    axis.set_title('3FGL Detection Flux Threshold [1/cm$^2$/s], 1.0 GeV - 10. GeV')
    imag = plt.imshow(fluxthrs[amin(jbgal):amax(jbgal)+1, amin(jlghprofi):amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, fraction=0.05)
    plt.savefig(plotpath + 'thrs_' + rtag + '.png')
    plt.close(figr)


# In[ ]:

def plot_fgl3():
    
    figr, axis = plt.subplots()
    bins = logspace(log10(amin(fgl3timevari[where(fgl3timevari > 0.)[0]])), log10(amax(fgl3timevari)), 100)
    axis.hist(fgl3timevari, bins=bins, label='All', log=True)
    axis.hist(fgl3timevari[indxfgl3rofi], bins=bins, label='ROI', log=True)
    axis.axvline(72.44, ls='--', alpha=0.5, color='black')
    axis.set_xlabel('3FGL time variability index')
    axis.set_xscale('log')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(plotpath + 'fgl3timevari.png')
    plt.close(figr)
    

    figr, axis = plt.subplots()
    indxfgl3scut = where(isfinite(fgl3scut))[0]
    bins = linspace(amin(fgl3scut[indxfgl3scut]), amax(fgl3scut[indxfgl3scut]), 100)
    axis.hist(fgl3scut[indxfgl3scut], bins=bins, label='All', log=True)
    axis.hist(fgl3scut[intersect1d(indxfgl3scut, indxfgl3rofi)], bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral cutoff')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(plotpath + 'fgl3scut.png')
    plt.close(figr)
    
    figr, axis = plt.subplots()
    indxfgl3scur = where(isfinite(fgl3scur))[0]
    bins = linspace(amin(fgl3scur[indxfgl3scur]), amax(fgl3scur[indxfgl3scur]), 100)
    axis.hist(fgl3scur[indxfgl3scur], bins=bins, label='All', log=True)
    axis.hist(fgl3scur[intersect1d(indxfgl3scur, indxfgl3rofi)], bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral curvature')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(plotpath + 'fgl3scur.png')
    plt.close(figr)
    
    figr, axis = plt.subplots()
    indxfgl3sind = where(isfinite(fgl3sind))[0]
    bins = linspace(amin(fgl3sind[indxfgl3sind]), amax(fgl3sind[indxfgl3sind]), 100)
    axis.hist(fgl3sind[indxfgl3sind], bins=bins, label='All', log=True)
    axis.hist(fgl3sind[intersect1d(indxfgl3sind, indxfgl3rofi)], bins=bins, label='ROI', log=True)
    axis.set_xlabel('3FGL spectral index')
    axis.legend()
    axis.set_ylim([0.1, None])
    plt.savefig(plotpath + 'fgl3sind.png')
    plt.close(figr)
    
    strgfgl3spectype = ['LogParabola', 'PLExpCutoff', 'PLSuperExpCutoff', 'PowerLaw']



# In[ ]:

def plot_samp():

    global thisresicnts, errrmodlcnts
    
    thisresicnts = datacnts - thismodlcnts

    thispsfn = retr_psfn(thissampvarb[indxsamppsfipara], indxener, angldisp, psfntype=modlpsfntype)
    if pixltype == 'cart':
        thispsfn = thispsfn.reshape((numbener, -1, numbevtt))
            

    thisfwhm = retr_fwhm(thispsfn)
    
    plot_psfn(thispsfn)
    
    

    global thisbackcntsmean
    thisbackcntsmean = empty((numbener, numbevtt))
    for c in indxback:
        thisbackcntsmean += mean(thissampvarb[indxsampnormback[c, :, None, None]] * backflux[c] * expo *                                     diffener[:, None, None] * pi * thisfwhm[:, None, :]**2 / 4., 1)

    
    
    thiscnts = []
    for l in indxpopl:
        ppixl = retr_pixl(thissampvarb[thisindxsampbgal[l]], thissampvarb[thisindxsamplgal[l]])
        cntstemp = thissampvarb[thisindxsampspec[l]][:, :, None] * expo[:, ppixl, :] * diffener[:, None, None]
        thiscnts.append(cntstemp)

    global jtruepntsbias, jtruepntsmiss, specmtch
    indxmodl, jtruepntsbias, jtruepntsmiss = pair_catl(truelgal[l],                          truebgal[l],                          truespec[l][0, :, :],                          thissampvarb[thisindxsamplgal[l]],                          thissampvarb[thisindxsampbgal[l]],                          thissampvarb[thisindxsampspec[l]])

    thisspecmtch = thissampvarb[thisindxsampspec[l]][:, indxmodl]
    thisspecmtch[:, jtruepntsmiss] = 0.


    if thissampvarb[indxsampnumbpnts[l]] > 1:
        for l in indxpopl:
            if colrprio:
                plot_histsind(l)
            plot_scatpixl(l)
            if trueinfo:
                plot_scatspec(l, thisspecmtch=thisspecmtch)
            plot_histspec(l)
            plot_histcnts(l, thiscnts)
            plot_compfrac()

            

    
    for i in indxener:
        
        plot_datacnts(i, None)
        #plot_catl(i, None, thiscnts)
        plot_modlcnts(i, None)
        plot_resicnts(i, None, thisresicnts)

        #for m in indxevtt:
        #    plot_datacnts(i, m)
        #    plot_catl(i, m, thiscnts)
        #    plot_modlcnts(i, m)
        #    plot_resicnts(i, m, thisresicnts)
        
    #if numbener == 3:
    #    plot_datacnts(None, None)
        
    #plot_fwhm(thisfwhm)
    #plot_backcntsmean(thisbackcntsmean)
    
        
    tempsampvarb, tempppixl, tempcnts, temppntsflux,         tempmodlflux, tempmodlcnts = pars_samp(thisindxpntsfull, drmcsamp[:, 0])
    errrmodlcnts = thismodlcnts - tempmodlcnts
    
    
    for i in indxener:
        plot_errrcnts(i, None, errrmodlcnts)
    


# In[ ]:

def plot_histcnts(l, thiscnts):

    figr, axgr = plt.subplots(numbevtt, numbener, figsize=(7 * numbener, 7 * numbevtt))
    if numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            if trueinfo:
                truehist = axis.hist(truecnts[l][i, :, m], binscnts[i, :], alpha=0.5, color='g', log=True, label=truelabl)
                if datatype == 'mock':
                    axis.hist(fgl3cnts[i, :, m], binscnts[i, :], alpha=0.5, color='red', log=True, label='3FGL')
            axis.hist(thiscnts[l][i, :, m], binscnts[i, :], alpha=0.5, color='b', log=True, label='Sample')
            if m == numbevtt - 1:
                axis.set_xlabel(r'$k$')
            axis.set_xscale('log')
            if trueinfo:
                axis.set_ylim([0.1, 1e3])
            if m == 0:
                axis.set_title(binsenerstrg[i])
            if i == 0 and exprtype == 'ferm':
                axis.set_ylabel(evttstrg[m])
            if m == numbevtt / 2 and i == numbener / 2:
                axis.legend()
        
    figr.subplots_adjust(wspace=0.3)
    plt.savefig(plotpath + 'histcnts%d_' % l + rtag + '_%09d.png' % j)
    plt.close(figr)
     


# In[ ]:

def plot_datacnts(pener, pevtt, nextstat=False):

    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        if pener == None:
            axis.set_title('')
        else:
            axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)


        
    # plot the model count map
    if pevtt == None:
        if pener == None:
            imag = axis.imshow(sum(datacntscart, axis=3), origin='lower', extent=extt, interpolation='none')
        else:
            imag = axis.imshow(sum(datacntscart[:, :, pener, :], axis=2), origin='lower',                              interpolation='none', cmap='Reds', extent=extt)
    else:
        imag = axis.imshow(datacntscart[:, :, pener, pevtt], origin='lower', interpolation='none',                          cmap='Reds', extent=extt)
    
    if pevtt != None or pener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
    
    # superimpose catalogs
    for l in indxpopl:


        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal[jtruepntsmiss], bgal[jtruepntsmiss], s=mrkrsize[jtruepntsmiss],                        alpha=mrkralph, label=truelabl + ', missed', marker='x', linewidth=2, color='g')
            axis.scatter(lgal[jtruepntsbias], bgal[jtruepntsbias], s=mrkrsize[jtruepntsbias],                        alpha=mrkralph, label=truelabl + ', biased', marker='o', linewidth=2, color='g')
            indxpnts = setdiff1d(arange(truenumbpnts, dtype=int), concatenate((jtruepntsbias, jtruepntsmiss)))
            axis.scatter(lgal[indxpnts], bgal[indxpnts], s=mrkrsize[indxpnts],                        alpha=mrkralph, label=truelabl + ', hit', marker='D', linewidth=2, color='g')
            for l in indxpopl:
                if jtruepntstimevari[l].size > 0:
                    axis.scatter(lgal[jtruepntstimevari[l]], bgal[jtruepntstimevari[l]], s=100,                                label=truelabl + ', variable', marker='*', linewidth=2, color='y')


        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
            
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')


    if nextstat:
        
        if False:
            print 'indxprop'
            print indxprop
            print 'reje'
            print reje
        
        if indxprop == indxpropsplt:
            print 'thislgal'
            print thislgal
            print 'thisbgal'
            print thisbgal
            print 'thisspec'
            print thisspec
            print 'nextlgal0'
            print nextlgal0
            print 'nextlgal1'
            print nextlgal1
            print 'nextbgal0'
            print nextbgal0
            print 'nextbgal1'
            print nextbgal1
            print 'nextspec0'
            print nextspec0
            print 'nextspec1'
            print nextspec1
            
        if indxprop == indxpropmerg:
            print 'thislgal0'
            print thislgal0
            print 'thisbgal0'
            print thisbgal0
            print 'thisspec0'
            print thisspec0
            print 'thislgal1'
            print thislgal1
            print 'thisbgal1'
            print thisbgal1
            print 'thisspec1'
            print thisspec1
            print 'nextlgal'
            print nextlgal
            print 'nextbgal'
            print nextbgal
            print 'nextspec'
            print nextspec  
        
        
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(abs(modispec[pener, k]), pener)
            
            if exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=mrkralph,                        marker='o', linewidth=2, color=colr)
            
            if False:
                print 'modilgal[k]'
                print modilgal[k]
                print 'modibgal[k]'
                print modibgal[k]
                print
                
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == indxpropsplt or indxprop == indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        if pener == None:
            path = plotpath + 'datacntsAA_' + rtag + '_%09d.png' % j
        else:
            path = plotpath + 'datacnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = plotpath + 'datacnts%d%d_' % (pener, jevtt[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)


# In[ ]:

def plot_modlcnts(pener, pevtt):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        modlcntstemp = sum(thismodlcnts[pener, :, :], axis=1)
    else:
        modlcntstemp = thismodlcnts[pener, :, pevtt]
    if pixltype == 'heal':
        modlcntstemp = tdpy_util.retr_cart(modlcntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=minmlgal, maxmlgal=maxmlgal,                                            minmbgal=minmbgal, maxmbgal=maxmbgal)
    else:
        modlcntstemp = modlcntstemp.reshape((nsidecart, nsidecart)).T
    modlcntstemp[where(modlcntstemp > datacntssatu[pener])] = datacntssatu[pener]
    
    imag = plt.imshow(modlcntstemp, origin='lower', cmap='Reds', extent=extt)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    
    # superimpose catalogs
    for l in indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = truelgal[l]
            bgal = truebgal[l]
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = plotpath + 'modlcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = plotpath + 'modlcnts%d%d_' % (pener, jevtt[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    


# In[ ]:

def plot_resicnts(pener, pevtt, resicnts, nextstat=False):

    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        
    # plot the model count map
    if pevtt == None:
        resicntstemp = sum(resicnts[pener, :, :], axis=1)
    else:
        resicntstemp = resicnts[pener, :, pevtt]
    if pixltype == 'heal':
        resicntstemp = tdpy_util.retr_cart(resicntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=minmlgal, maxmlgal=maxmlgal,                                            minmbgal=minmbgal, maxmbgal=maxmbgal)
    else:
        resicntstemp = resicntstemp.reshape((nsidecart, nsidecart))
    resicntstemp[where(resicntstemp > resicntssatu[pener])] = resicntssatu[pener]
    resicntstemp[where(resicntstemp < -resicntssatu[pener])] = -resicntssatu[pener]
    
    imag = axis.imshow(resicntstemp, origin='lower', cmap='RdGy', extent=extt)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    # superimpose catalogs
    for l in indxpopl:

        # model catalog
        mrkrsize = retr_mrkrsize(thissampvarb[thisindxsampspec[l]][pener, :], pener)
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            mrkrsize = retr_mrkrsize(truespec[l][0, pener, :], pener)
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=mrkrsize, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

        
    if nextstat:
        for k in range(modilgal.size):
            if modispec[pener, k] > 0:
                colr = 'yellow'
            else:
                colr = 'red'
            mrkrsize = retr_mrkrsize(abs(modispec[pener, k]), pener)
            
            
            if exprtype == 'ferm':
                xaxi = modilgal[k]
                yaxi = modibgal[k]
            else:
                xaxi = modilgal[k] * 3600.
                yaxi = modibgal[k] * 3600.
            axis.scatter(xaxi, yaxi, s=mrkrsize, alpha=mrkralph,                        marker='o', linewidth=2, color=colr)

            
            text = r'indxprop = %d, lnL = %.3g' % (indxprop, deltllik)
            if indxprop == indxpropsplt or indxprop == indxpropmerg:
                text += ', lnF = %.3g, lnJ = %.3g, lnC = %.3g, %d' % (laccfrac, thisjcbnfact, thiscombfact, listaccp[j])
            axis.text(0.6, 1.1, text, va='center', ha='center', transform=axis.transAxes, fontsize=18)

            
    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')

    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        path = plotpath + 'resicnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = plotpath + 'resicnts%d%d_' % (pener, jevtt[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    
    plt.close(figr)
    


# In[ ]:

def plot_errrcnts(pener, pevtt, errrcntsrofi):

    if pevtt == None:
        errrcntstemp = sum(errrcntsrofi[pener, :, :], axis=1)
    else:
        errrcntstemp = errrcntsrofi[pener, :, pevtt]
    
    if pixltype == 'heal':
        errrcntstemp = tdpy_util.retr_cart(errrcntstemp, jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                            minmlgal=minmlgal, maxmlgal=maxmlgal,                                            minmbgal=minmbgal, maxmbgal=maxmbgal)
    else:
        errrcntstemp = errrcntstemp.reshape((nsidecart, nsidecart))
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)

    # plot the error count map
    imag = axis.imshow(errrcntstemp, origin='lower', cmap='Reds', extent=extt)
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    if pevtt == None:
        path = plotpath + 'errrcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = plotpath + 'errrcnts%d%d_' % (pener, jevtt[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    


# In[ ]:

def plot_catl(pener, pevtt, thiscnts):
    
    # begin figure
    figr, axis = plt.subplots(figsize=(12, 12))
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    axis.set_xlim([frambndrmarg, -frambndrmarg])
    axis.set_ylim([-frambndrmarg, frambndrmarg])
    if pevtt == None:
        axis.set_title(binsenerstrg[pener])
    else:
        titl = binsenerstrg[pener]
        if exprtype == 'ferm':
            titl += ', ' + evttstrg[pevtt]
        axis.set_title(titl)
        

    # superimpose catalogs
    for l in indxpopl:
        
        # model catalog
        lgal = thissampvarb[thisindxsamplgal[l]]
        bgal = thissampvarb[thisindxsampbgal[l]]
        if exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        axis.scatter(lgal, bgal, s=300, alpha=mrkralph, label='Sample', marker='+', linewidth=2, color='b')

        # true catalog
        if trueinfo:
            lgal = copy(truelgal[l])
            bgal = copy(truebgal[l])
            if exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal, bgal, s=300, alpha=mrkralph, label=truelabl, marker='x', linewidth=2, color='g')

    
    
    if trueinfo:
        for l in indxpopl:
            numbpnts = int(truenumbpnts[l])
            for a in range(numbpnts):
                if pevtt == None:
                    cnts = sum(truecnts[l][pener, a, :])
                    sigm = sqrt(sum(truesigm[l][pener, a, :]**2))
                else:
                    cnts = truecnts[l][pener, a, pevtt]
                    sigm = truesigm[l][pener, a, pevtt]
                axis.text(truelgal[l][a] + 0.7, truelgal[l][a] - 0.7,                         '%d/%.2f' % (cnts, sigm), color='g', fontsize=13)

    for l in indxpopl:
        numbpnts = int(thissampvarb[indxsampnumbpnts[l]])
        for a in range(numbpnts):
            if pevtt == None:
                cnts = sum(thiscnts[l][pener, a, :])
            else:
                cnts = thiscnts[l][pener, a, pevtt]
            axis.text(thissampvarb[thisindxsamplgal[l][a]] - 0.5, thissampvarb[thisindxsampbgal[l][a]] + 0.3,                     '%d' % cnts, color='b', fontsize=13)
    
    axis.axvline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-frambndr, ls='--', alpha=0.3, color='black')
    
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)
    
    if pevtt == None:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % mean(thisbackcntsmean[pener, :]), fontsize=18)
    else:
        plt.figtext(0.76, 0.92, '$C_{back} = %d$' % thisbackcntsmean[pener, pevtt], fontsize=18)
        
    if pevtt == None:
        path = plotpath + 'catlcnts%dA_' % pener + rtag + '_%09d.png' % j
    else:
        path = plotpath + 'catlcnts%d%d_' % (pener, jevtt[pevtt]) + rtag + '_%09d.png' % j
    plt.savefig(path)
    
    plt.close(figr)



    
    
    



# In[ ]:

def plot_topo():
    
        
    figr, axis = plt.subplots()
    datacnts = zeros(npixlheal)
    datacnts[jpixlrofi] = datacnts[0,:,3]
    testcart = tdpy_util.retr_cart(datacnts, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
    imag = axis.imshow(testcart, origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    axis.scatter(mocklgal, mockbgal, c='g')
    axis.scatter(lgal, bgal, c='b', s=50)
    axis.set_xlim(extent[0:2])
    axis.set_ylim(extent[2:4])
    axis.set_title('PSF3 Mock count map')
    plt.show()
    
    data = linspace(0., 1., numbbins)

    # parameters of interest
    indxpnts = choice(arange(numbpnts), size=3, replace=False)
    jpara = sort(concatenate((indxsampnormback[1, :], indxsamplgal[indxpnts], indxsampbgal[indxpnts], indxsampspec[0, indxpnts])))
    
    
    labl = ['A', '$\psi_1$', '$\phi_1$', '$f_1$', '$\psi_2$', '$\phi_2$', '$f_2$', '$\psi_3$', '$\phi_3$', '$f_3$']
    labllist = list(itertools.combinations(labl, 2))
    indxlist = list(itertools.combinations(jpara, 2))
    ncomb = len(labllist)
    

    llik = zeros((ncomb, numbbins, numbbins))        
    indxsamp = arange(5, npara, dtype=int)
    xlog = ones(ncomb, dtype=bool)
    ylog = ones(ncomb, dtype=bool)
    extt = zeros((ncomb, 4))
    thiscntr = -1
    for k in range(ncomb):
        
        for i in range(2):
            if indxlist[k][i] == indxsampnormback[1, :]:
                extt[k, 2 * i] = minmnfdm
                extt[k, 2 * i + 1] = maxmnfdm
            elif where(indxlist[k][i] == indxsamplgal[indxpnts])[0].size > 0:
                extt[k, 2 * i] = minmgang
                extt[k, 2 * i + 1] = maxmgang
            elif where(indxlist[k][i] == indxsampbgal[indxpnts])[0].size > 0:
                extt[k, 2 * i] = 0.
                extt[k, 2 * i + 1] = 2. * pi
                if i == 0:
                    xlog[k] = False
                else:
                    ylog[k] = False
            else:
                extt[k, 2 * i] = hyprmflx[0]
                extt[k, 2 * i + 1] = maxmflux

            

    for k in range(ncomb):

        # set all the other parameters to random values
        indxsamprand = setdiff1d(indxsamp, array(indxlist[k]))
        samp[indxsamprand] = rand()

        # scan over the two parameters of interest
        for n in range(numbbins):
            for m in range(numbbins):
                samp[indxlist[k][0]] = data[n]
                samp[indxlist[k][1]] = data[m]
                llik[k, n, m] = retr_llik(samp)

        thiscntr = tdpy_util.show_prog(k, ncomb, thiscntr)
            
    
    llik *= 1e-7
    nrows = 15
    ncoln = 3
    figr, axgr = plt.subplots(nrows, ncoln, figsize=(24, 90))
    figr.suptitle('Likelihood Topology', fontsize=40)
    for x, axrw in enumerate(axgr):
        for y, axis in enumerate(axrw):  
            k = x * ncoln + y
            im = axis.imshow(llik[k, :, :], origin='lower', cmap='Reds', extent=extt[k,:], aspect='auto')
            figr.colorbar(im, ax=axis, fraction=0.04)
            axis.set_xlabel(labllist[k][0])
            axis.set_ylabel(labllist[k][1])
            if xlog[k]:
                axis.set_xscale('log')
            if ylog[k]:
                axis.set_yscale('log')
                  
    figr.subplots_adjust(top=0.97, hspace=0.4, wspace=0.4)
    plt.show()


# In[ ]:

def plot_pntsdiff():

    tempindxpntsfull = []
    tempnumbpnts = 10
    tempindxpntsfull.append(range(tempnumbpnts))
    tempsamp = rand(maxmsampsize)
    tempsamp[indxsampnumbpnts] = array([tempnumbpnts])
    tempsamp[indxsampfdfnnorm] = cdfn_logt(array([tempnumbpnts]), minmfdfnnorm, factfdfnnorm)
    tempsamp[indxsampfdfnslop] = cdfn_atan(array([1.5]), minmfdfnslop, factfdfnslop)
    tempsamp[indxsamppsfipara] = 0.5
    for c in indxback:
        tempsamp[indxsampnormback[c, 0]] = cdfn_logt(array([1.]), minmnormback[c], factnormback[c])
    tempsampvarb, tempppixl, tempcnts,         temppntsflux, tempflux, tempcntsrofi = pars_samp(tempindxpntsfull, tempsamp)


    pntscnts = tempcnts[0, :, 0] * expo[0, :, 0] * apix * diffener[0]
    isotcnts = isotflux[0, :, 0] * expo[0, :, 0] * apix * diffener[0]
    
    totlcnts0 = isotcnts
    totlcnts1 = pntscnts / 1e6 + isotcnts
    
    nrept = 100
    totlcntsheal0 = zeros((nrept, npixlheal))
    totlcntsheal1 = zeros((nrept, npixlheal))
    for k in range(nrept):
        for n in range(jpixlrofi.size):
            totlcntsheal0[k, jpixlrofi[n]] = poisson(totlcnts0[n])
            totlcntsheal1[k, jpixlrofi[n]] = poisson(totlcnts1[n])
          
    maxmcnts = max(amax(totlcntsheal0), amax(totlcntsheal1))
    binscnts = linspace(0., maxmcnts, maxmcnts + 2)
    diffcnts = binscnts[1:] - binscnts[:-1]
    meancnts = (binscnts[1:] + binscnts[:-1]) / 2.
    
    hist0 = empty((nrept, maxmcnts + 1))
    hist1 = empty((nrept, maxmcnts + 1))
    for k in range(nrept):
        hist0[k, :] = histogram(totlcntsheal0[k, jpixlrofi], binscnts)[0].astype(float)
        hist0[k, :] *= 1. / sum(hist0[k, :]) / diffcnts
        hist1[k, :] = histogram(totlcntsheal1[k, jpixlrofi], binscnts)[0].astype(float)
        hist1[k, :] *= 1. / sum(hist1[k, :]) / diffcnts
        
    
    totlcntscart0 = tdpy_util.retr_cart(totlcntsheal0[0, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
    totlcntscart1 = tdpy_util.retr_cart(totlcntsheal1[0, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
    
    fig = plt.figure(figsize=(12, 12))
    axis = figr.add_subplot(221)
    imag = axis.imshow(totlcntscart0, origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, ax=axis, fraction=0.05)
    
    axis.set_xlabel(longlabl)
    axis.set_ylabel(latilabl)
    if exprtype == 'ferm':
        axis.set_xlim([maxmlgal, minmlgal])
        axis.set_ylim([minmbgal, maxmbgal])
    else:
        axis.set_xlim(array([maxmlgal, minmlgal]) * 3600.)
        axis.set_ylim(array([minmbgal, maxmbgal]) * 3600.)
    
    axis.set_title('Isotropic')

    
    axis = figr.add_subplot(222)
    imag = axis.imshow(totlcntscart1, origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, ax=axis, fraction=0.05)

    axis.set_xlabel(r'$x$ [$^\circ$]')
    axis.set_xlim([maxmlgal, minmlgal])
    axis.set_ylim([minmbgal, maxmbgal])
    axis.set_title('Isotropic + Unresolved PS')
    
    axis.scatter(tempsampvarb[trueindxsamplgal], tempsampvarb[trueindxsampbgal],                s=50, alpha=0.8, marker='x', color='g', linewidth=2)
    
    axis = figr.add_subplot(212)


    tdpy_util.plot_braz(ax, meancnts, hist0,  lcol='lightgreen',                         alpha=0.3, dcol='green', mcol='darkgreen', labl='Isotropic')
    tdpy_util.plot_braz(ax, meancnts, hist1, lcol='lightblue',                         alpha=0.3, dcol='blue', mcol='darkblue', labl='Isotropic + Unresolved PS')

    axis.set_title('Count PDF')
    axis.set_xlabel('$k$')
    axis.set_ylabel('d$P$/d$k$')
    axis.set_yscale('log')
    axis.legend()

    plt.savefig(plotpath + 'pntsdiff.png')
    plt.close(figr)
    



# In[ ]:

def plot_psfn_type():
    
    devi = linspace(0., 5., 100)
    y = zeros((x.size, 5))

    figr, axis = plt.subplots(figsize=(10, 6))
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


# In[ ]:

def plot_post(pathprobcatl):
          
    hdun = pf.open(pathprobcatl)

    hdun.info()

    global numbener, numbevtt, numbpopl
    numbener = hdun[0].header['numbener']
    numbevtt = hdun[0].header['numbevtt']
    numbpopl = hdun[0].header['numbpopl']
    numbpsfipara = hdun[0].header['numbpsfipara']
    nformpara = hdun[0].header['nformpara']
    
    global indxener, indxpopl
    indxener = arange(numbener)
    indxevtt = arange(numbevtt)
    indxpopl = arange(numbpopl)
    
    numbsamp = hdun[0].header['numbsamp']
    numbburn = hdun[0].header['numbburn']
    numbswep = hdun[0].header['numbswep']
    factthin = hdun[0].header['factthin']
    numbpopl = hdun[0].header['numbpopl']
    numbproc = hdun[0].header['numbproc']

    global numbspec
    maxmgang = hdun[0].header['maxmgang']
    lgalcntr = hdun[0].header['lgalcntr']
    bgalcntr = hdun[0].header['bgalcntr']
    numbspec = hdun[0].header['numbspec']
    nsideheal = hdun[0].header['nsideheal']
    nsidecart = hdun[0].header['nsidecart']

    minmlgal = hdun[0].header['minmlgal']
    maxmlgal = hdun[0].header['maxmlgal']
    minmbgal = hdun[0].header['minmbgal']
    maxmbgal = hdun[0].header['maxmbgal']

    global datatype, regitype, modlpsfntype, exprtype, pixltype
    datatype = hdun[0].header['datatype']
    regitype = hdun[0].header['regitype']
    modlpsfntype = hdun[0].header['modlpsfntype']
    if datatype == 'mock':
        mockpsfntype = hdun[0].header['mockpsfntype']
    exprtype = hdun[0].header['exprtype']
    pixltype = hdun[0].header['pixltype']
    
    global trueinfo
    trueinfo = hdun[0].header['trueinfo']
    colrprio = hdun[0].header['colrprio']
    
    margsize = hdun[0].header['margsize']
    datetimestrg = hdun[0].header['datetimestrg']

    
    levi = hdun[0].header['levi']
    info = hdun[0].header['info']
    
    global rtag
    rtag = retr_rtag(None)
    
    global plotpath
    if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
        plotfold = '/n/pan/www/tansu/png/pnts_tran/'
    else:
        plotfold = os.environ["PNTS_TRAN_DATA_PATH"] + '/png/'
    plotpath = plotfold + datetimestrg + '_' + rtag + '/'
    cmnd = 'mkdir -p ' + plotpath
    os.system(cmnd)


    lgalheal, bgalheal, nsideheal, npixlheal, apix = tdpy_util.retr_heal(nsideheal)
    jpixlrofi = where((abs(lgalheal) < maxmgang) & (abs(bgalheal) < maxmgang))[0]


    global jener, jevtt
    jener = hdun['jener'].data
    jevtt = hdun['jevtt'].data
    
    
    minmpsfipara, maxmpsfipara, factpsfipara, strgpsfipara,         scalpsfipara, jpsfiparainit = retr_psfimodl(nformpara, exprtype, psfntype, numbener, numbevtt, jener, jevtt)

        
    listlgal = []
    listbgal = []
    listspec = []
    if colrprio:
        listsind = []
    for l in indxpopl:
        listlgal.append(hdun['lgalpopl%d' % l].data)
        listbgal.append(hdun['bgalpopl%d' % l].data)
        listspec.append(hdun['specpopl%d' % l].data)
        if colrprio:
            listsind.append(hdun['sindpopl%d' % l].data)
        
    listspechist = hdun['spechist'].data
    pntsprob = hdun['pntsprob'].data
    
    listllik = hdun['llik'].data
    listlpri = hdun['lpri'].data
    
    gmrbstat = hdun['gmrbstat'].data
    listmodlcnts = hdun['modlcnts'].data
    

    # truth information
    global truelgal, truebgal, truespec
    if trueinfo:
        truenumbpnts = hdun['truenumbpnts'].data
        
        truelgal = []
        truebgal = []
        truespec = []
        truesind = []
        truespec = []
        for l in indxpopl:
            truelgal.append(hdun['truelgalpopl%d' % l].data)
            truebgal.append(hdun['truebgalpopl%d' % l].data)
            truespec.append(hdun['truespecpopl%d' % l].data)
        truefdfnnorm = hdun['truefdfnnorm'].data
        truefdfnslop = hdun['truefdfnslop'].data
        truenormback = hdun['truenormback'].data
        truepsfipara = hdun['truepsfipara'].data
        
    
    # prior boundaries
    if colrprio:
        global minmsind, maxmsind, binssind, meansind
        minmsind = hdun['minmsind'].data
        maxmsind = hdun['maxmsind'].data
        binssind = hdun['binssind'].data
        meansind = hdun['meansind'].data
        
    global minmspec, maxmspec, binsspec, meanspec
    minmspec = hdun['minmspec'].data
    maxmspec = hdun['maxmspec'].data
    binsspec = hdun['binsspec'].data
    meanspec = hdun['meanspec'].data
    
    global binsener, meanener
    binsener = hdun['binsener'].data
    meanener = hdun['meanener'].data
    
    # utilities
    probprop = hdun['probprop'].data
    listindxprop = hdun['indxprop'].data
    listchro = hdun['chro'].data
    listaccp = hdun['accp'].data
    listindxsampmodi = hdun['sampmodi'].data
    listauxipara = hdun['auxipara'].data
    listlaccfrac = hdun['laccfrac'].data
    listnumbpair = hdun['numbpair'].data
    listjcbnfact = hdun['jcbnfact'].data
    listcombfact = hdun['combfact'].data
    listpntsfluxmean = hdun['listpntsfluxmean'].data

    listnumbpnts = hdun['numbpnts'].data
    listfdfnnorm = hdun['fdfnnorm'].data
    listfdfnslop = hdun['fdfnslop'].data
    listpsfipara = hdun['psfipara'].data
    listnormback = hdun['normback'].data
    
    numbprop = len(probprop)
    
    global maxmgangmarg
    maxmgangmarg = maxmgang + margsize
     
        
    global extt, frambndr, frambndrmarg
    extt = [minmlgal, maxmlgal, minmbgal, maxmbgal]
    if exprtype == 'sdss':
        extt *= 3600.
        frambndr = maxmgang * 3600.
        frambndrmarg = maxmgangmarg * 3600.
    else:
        frambndr = maxmgang
        frambndrmarg = maxmgangmarg
        
    global strgfluxunit
    strgfluxunit = retr_strgfluxunit(exprtype)

    # energy bin string
    global binsenerstrg, enerstrg
    enerstrg, binsenerstrg = retr_enerstrg(exprtype)
    
    global longlabl, latilabl
    if regitype == 'igal':
        longlabl = '$l$'
        latilabl = '$b$'
    else:
        longlabl = r'$\nu$'
        latilabl = r'$\mu$'
        
    if exprtype == 'ferm':
        longlabl += r' [$^\circ$]'
        latilabl += r' [$^\circ$]'
    if exprtype == 'sdss':
        longlabl += ' [arcsec]'
        latilabl += ' [arcsec]'
        
    if trueinfo:
        global truelabl
        if datatype == 'mock':
            truelabl = 'Mock data'
        else:
            truelabl = '3FGL'

        
    # Gelman-Rubin test
    if numbproc > 1:
        print 'Making the Gelman-Rubin TS plot...'
        tim0 = time.time()
        
        figr, axis = plt.subplots()
        axis.hist(gmrbstat, bins=linspace(1., amax(gmrbstat), 40))
        axis.set_title('Gelman-Rubin Convergence Test')
        axis.set_xlabel('PSRF')
        axis.set_ylabel('$N_{pix}$')
        figr.savefig(plotpath + 'gmrbdist_' + rtag + '.png')
        plt.close(figr)
        
        tim1 = time.time()
        print 'Done in %.3g seconds.' % (tim1 - tim0)

            
    if False:
        print 'Regressing the data map against model maps from the posterior...'
        tim0 = time.time()

        # Linear regression of data and model maps
        slop = zeros(numbsamp * numbproc)
        intc = zeros(numbsamp * numbproc)
        coef = zeros(numbsamp * numbproc)
        for p in range(numbsamp * numbproc):
            slop[p], intc[p], coef[p], pval, stdv =                 sp.stats.linregress(datacnts[0, gpixl, 0], listmodlcnts.reshape((numbsamp * numbproc), ngpixl)[p, :])
        figr, axcl = plt.subplots(1, 3, figsize=(12, 4))
        for a, axis in enumerate(axcl):
            if a == 0:
                axis.hist(slop, 40)
                axis.set_ylabel('$N_{samp}$')
                axis.set_xlabel('$m$')
            if a == 1:
                axis.hist(intc, 40)
                axis.set_xlabel(r'$\beta$')
            if a == 2:
                axis.hist(coef, 40)
                axis.set_xlabel('$R$')
        figr.savefig(plotpath + 'lregdist_' + rtag + '.png')
        plt.close(figr)

        tim1 = time.time()
        print 'Done in %.3g seconds.' % (tim1 - tim0)

        
    print 'Calculating proposal and acceptance rates...'
    tim0 = time.time()

    mcmctimebins = linspace(0., numbswep, 100)
    figr, axgr = plt.subplots(numbprop, 1, figsize=(10, 15), sharex='all')
    for g, axis in enumerate(axgr):
        axis.hist(where(listindxprop == g)[0], bins=mcmctimebins)
        axis.hist(where((listindxprop == g) & (listaccp == True))[0], bins=mcmctimebins)
        axis.set_ylabel('%d' % g)
        if g == numbprop - 1:
            axis.set_xlabel('$i_{samp}$')
    figr.savefig(plotpath + 'propeffi_' + rtag + '.png')
    figr.subplots_adjust(hspace=0.)
    plt.close(figr)
     

    indxsampsplt = where(listindxprop == indxpropsplt)[0]
    indxsampmerg = where(listindxprop == indxpropmerg)[0]
            
    listname = ['laccfrac', 'numbpair', 'combfact', 'jcbnfact']
    listvarb = [listlaccfrac, listnumbpair, listcombfact, listjcbnfact]
    for k in range(4):
        figr, axis = plt.subplots()
        axis.hist(listvarb[k][indxsampsplt])
        axis.hist(listvarb[k][indxsampmerg])
        axis.set_ylabel('$N_{samp}$')
        axis.set_xlabel(listname[k])
        figr.subplots_adjust(bottom=0.2)
        figr.savefig(plotpath + listname[k] + rtag + '.png')
        plt.close(figr)
    
    
    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Calculating proposal execution times...'
    tim0 = time.time()

    binstime = logspace(log10(amin(listchro[where(listchro > 0)] * 1e3)), log10(amax(listchro * 1e3)), 50)
    with sns.color_palette("nipy_spectral", numbprop):
        figr, axis = plt.subplots(figsize=(14, 12))
        
        axis.hist(listchro[where(listchro[:, 0] > 0)[0], 0] * 1e3, binstime, facecolor='none', log=True, lw=1, ls='--', edgecolor='black')
        for g in range(numbprop):
            indxlistchro = where((listindxprop == g) & (listchro[:, 0] > 0))[0]
            # temp
            axis.hist(listchro[indxlistchro, 0] * 1e3, binstime, edgecolor='none', log=True, alpha=0.3)#, label=propstrg[g])
        axis.set_xlabel('$t$ [ms]')
        axis.set_xscale('log')
        axis.set_xlim([amin(binstime), amax(binstime)])
        axis.set_ylim([0.5, None])
        axis.legend(loc=2)
        figr.savefig(plotpath + 'chroprop_' + rtag + '.png')
        plt.close(figr)

    labl = ['Total', 'Proposal', 'Prior', 'Likelihood']
    figr, axcl = plt.subplots(2, 1, figsize=(14, 10))
    for k in range(1, 4):
        axcl[0].hist(listchro[where(listchro[:, k] > 0)[0], k] * 1e3, binstime, log=True, alpha=0.5, label=labl[k])
    axcl[1].hist(listchro[where(listchro[:, 0] > 0)[0], 0] * 1e3, binstime, log=True, label=labl[0], color='black')
    axcl[1].set_title(r'$\langle t \rangle$ = %.3g ms' % mean(listchro[where(listchro[:, 0] > 0)[0], 0] * 1e3))
    axcl[0].set_xlim([amin(binstime), amax(binstime)])
    axcl[1].set_xlabel('$t$ [ms]')
    axcl[0].set_xscale('log')
    axcl[1].set_xscale('log')
    axcl[0].set_ylim([0.5, None])
    axcl[1].set_ylim([0.5, None])
    axcl[0].legend(loc=1)
    axcl[1].legend(loc=2)
    figr.savefig(plotpath + 'chrototl_' + rtag + '.png')
    plt.close(figr)

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Parsing the sample bundle and making frames...'
    tim0 = time.time()

    
    if False:
        if maxmnumbpnts[0] == 4 and numbener == 1 and numbpopl == 1:
            for k in range(maxmnumbpnts[0]):
                listpostdist = zeros((numbsamp, 3 * k)) 
                listpostdist[:, 0*k:1*k] = listlgalpnts[0][k]
                listpostdist[:, 1*k:2*k] = listbgalpnts[0][k]
                listpostdist[:, 2*k:3*k] = listspecpnts[0][k][0, :]
                path = plotpath + 'postdist_%d_' % k + rtag
                if trueinfo:
                    truepara = truesampvarb[indxsamppsfipara[ipsfipara]]
                else:
                    truepara = None
                tdpy_util.plot_mcmc(listpostdist, parastrgpsfipara[ipsfipara],                                     truepara=truepara, path=path, numbbins=numbbins, quan=True)
                
                
            
    # flux match with the true catalog
    if trueinfo:
        for l in indxpopl:
            
            discspecmtch = zeros(truenumbpnts) + numbsamp
            listindxmodl = []
            for k in range(numbsamp):
                indxmodl, jtruepntsbias, jtruepntsmiss = pair_catl(truelgal[l],                                      truebgal[l],                                      truespec[l][0, :, :],                                      listlgal[l][k, :],                                      listbgal[l][k, :],                                      listspec[l][k, :, :])
                listindxmodl.append(indxmodl)
                discspecmtch[jtruepntsmiss] -= 1.
            discspecmtch /= numbsamp

            
            postspecmtch = zeros((3, numbener, truenumbpnts))
            for i in indxener:
                listspecmtch = zeros((numbsamp, truenumbpnts))
                for k in range(numbsamp):
                    indxpntstrue = where(listindxmodl[k] >= 0)[0]
                    listspecmtch[k, indxpntstrue] = listspec[l][k][i, listindxmodl[k][indxpntstrue]]
                postspecmtch[0, i, :] = percentile(listspecmtch, 16., axis=0)
                postspecmtch[1, i, :] = percentile(listspecmtch, 50., axis=0)
                postspecmtch[2, i, :] = percentile(listspecmtch, 84., axis=0)
            
            plot_scatspec(l, postspecmtch=postspecmtch)

            # store the comparison
            path = os.environ["PNTS_TRAN_DATA_PATH"] + '/pcatcomp_popl%d_' % l + rtag + '.fits'
            compbund = stack((truespec[l], postspecmtch))

    tim1 = time.time()
    print 'Done in %.3g seconds.' % (tim1 - tim0)
    print 'Making posterior distribution plots...'

    if pixltype == 'heal':
        
        reso = 60. / nsideheal
        numbbinslgcr = (maxmlgal - minmlgal) / reso
        numbbinsbgcr = (maxmbgal - minmbgal) / reso

        pntsprobcart = zeros((numbbinslgcr, numbbinsbgcr, numbpopl, numbener, numbspec))
        for l in indxpopl:
            for i in indxener:
                for h in range(numbspec):
                    pntsprobcart[:, :, l, i, h] = tdpy_util.retr_cart(pntsprob[l, i, :, h], jpixlrofi=jpixlrofi, nsideinpt=nsideheal,                                                                       minmlgal=minmlgal, maxmlgal=maxmlgal,                                                                       minmbgal=minmbgal, maxmbgal=maxmbgal, reso=reso)
    else:
        pntsprobcart = pntsprob.reshape((numbpopl, numbener, nsidecart, nsidecart, numbspec))
        pntsprobcart = swapaxes(swapaxes(pntsprobcart, 0, 2), 1, 3)
        
    # stacked posteiors binned in position and flux
    plot_pntsprob(pntsprobcart, ptag='quad')
    plot_pntsprob(pntsprobcart, ptag='full', full=True)
    plot_pntsprob(pntsprobcart, ptag='cumu', cumu=True)
   
    # flux distribution
    for l in indxpopl:
        plot_histspec(l, listspechist=listspechist[:, l, :, :])

    # fraction of emission components
    postpntsfluxmean = retr_postvarb(listpntsfluxmean)
    postnormback = retr_postvarb(listnormback)
    plot_compfrac(postpntsfluxmean=postpntsfluxmean, postnormback=postnormback)


    # PSF parameters
    path = plotpath + 'psfipara_' + rtag
    if modlpsfntype == 'singgaus' or modlpsfntype == 'singking':
        listpsfipara[:, jpsfiparainit] = rad2deg(listpsfipara[:, jpsfiparainit])
        if trueinfo:
            truepsfipara[jpsfiparainit] = rad2deg(truepsfipara[jpsfiparainit])
    elif modlpsfntype == 'doubgaus' or modlpsfntype == 'gausking':
        listpsfipara[:, jpsfiparainit+1] = rad2deg(listpsfipara[:, jpsfiparainit+1])
        listpsfipara[:, jpsfiparainit+2] = rad2deg(listpsfipara[:, jpsfiparainit+2])
        if trueinfo:
            truepsfipara[jpsfiparainit+1] = rad2deg(truepsfipara[jpsfiparainit+1])
            truepsfipara[jpsfiparainit+2] = rad2deg(truepsfipara[jpsfiparainit+2])
    elif modlpsfntype == 'doubking':
        listpsfipara[:, jpsfiparainit+1] = rad2deg(listpsfipara[:, jpsfiparainit+1])
        listpsfipara[:, jpsfiparainit+3] = rad2deg(listpsfipara[:, jpsfiparainit+3])
        if trueinfo:
            truepsfipara[jpsfiparainit+1] = rad2deg(truepsfipara[jpsfiparainit+1])
            truepsfipara[jpsfiparainit+3] = rad2deg(truepsfipara[jpsfiparainit+3])
    if trueinfo and psfntype == 'doubking':
        truepara = truepsfipara
    else:
        truepara = array([None] * numbpsfipara)
    tdpy_util.plot_mcmc(listpsfipara, strgpsfipara, truepara=truepara,                         nplot=nformpara, path=path, numbbins=numbbins, quan=True, ntickbins=3)
    
    for k in range(numbpsfipara):
        path = plotpath + 'psfipara%d_' % k + rtag + '.png'
        tdpy_util.plot_trac(listpsfipara[:, k], strgpsfipara[k], path=path, quan=True)
    
    
    # log-likelihood
    path = plotpath + 'llik_' + rtag + '.png'
    tdpy_util.plot_trac(listllik.flatten(), '$P(D|y)$', path=path)

    # log-prior
    path = plotpath + 'lpri_' + rtag + '.png'
    tdpy_util.plot_trac(listlpri.flatten(), '$P(y)$', path=path)
    

    # number, expected number of PS and flux conditional prior power law index 
    for l in range(numbpopl):
        
        # number of point sources
        path = plotpath + 'numbpntsdist_popl%d_' % l + rtag + '.png'
        if trueinfo and truenumbpnts != None:
            truepara = truenumbpnts[l]
        else:
            truepara = None
        tdpy_util.plot_trac(listnumbpnts[:, l], '$N$', truepara=truepara, path=path)

        # mean number of point sources
        path = plotpath + 'fdfnnorm_popl%d_' % l + rtag + '.png'
        if trueinfo and truefdfnnorm != None:
            truepara = truefdfnnorm[l]
        else:
            truepara = None
        tdpy_util.plot_trac(listfdfnnorm[:, l], '$\mu$', truepara=truepara, path=path)

        # flux distribution power law index
        for i in indxenerfdfn:
            path = plotpath + 'fdfnslopdist_popl%d%d_' % (l, i) + rtag + '.png'
            if trueinfo and truefdfnslop != None:
                truepara = truefdfnslop[l, i]
            else:
                truepara = None
            titl = binsenerstrg[i]
            labl =  r'$\alpha_{%d}$' % i
            tdpy_util.plot_trac(listfdfnslop[:, l, i], labl, truepara=truepara, path=path, titl=titl)
        
        
    # isotropic background normalization
    for i in indxener:
        path = plotpath + 'nisodist%d_' % i + rtag + '.png'
        if trueinfo:
            if datatype == 'mock':
                truepara = truenormback[0, i]
            else:
                truepara = None
        titl = binsenerstrg[i]
        labl = r'$\mathcal{I}_{%d}$' % i
        tdpy_util.plot_trac(listnormback[:, 0, i], labl, truepara=truepara, path=path, titl=titl)
       
    if exprtype == 'ferm':
        # diffuse model normalization
        for i in indxener:
            path = plotpath + 'nfdmdist%d_' % i + rtag + '.png'
            if trueinfo:
                if datatype == 'mock':
                    truepara = truenormback[1, i]
                else:
                    truepara = None
            else:
                truepara = None
            titl = binsenerstrg[i]
            labl = r'$\mathcal{D}_{%d}$' % i
            tdpy_util.plot_trac(listnormback[:, 1, i], labl, truepara=truepara, path=path, titl=titl)

    # plot log-likelihood
    figr, axrw = plt.subplots(2, 1, figsize=(7, 12))
    figr.suptitle(r'$D_{KL} = %.5g, \ln P(D) = %.5g$' % (info, levi))
    for k, axis in enumerate(axrw):
        if k == 0:
            if amin(listllik) != amax(listllik):
                axis.hist(listllik.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(D|x)$')
        else:
            if amin(listlpri) != amax(listllik):
                axis.hist(listlpri.flatten())
                axis.set_ylabel(r'$N_{samp}$')
                axis.set_xlabel(r'$\ln P(x)$')
    figr.savefig(plotpath + 'leviinfo_' + rtag + '.png')
    plt.close(figr)

    #make_anim()

    tim1 = time.time()
    print 'Plots are produced in %.3g seconds.' % (tim1 - tim0)



# In[ ]:

def plot_compfrac(postpntsfluxmean=None, postnormback=None):
    
    if postpntsfluxmean != None:
        post = True
    else:
        post = False
        
        
    listlinestyl = ['-', '--', '-.', ':']
    listcolr = ['black', 'b', 'b', 'b']
    listlabl = ['Data', 'PS', 'Iso', 'FDM']

    figr, axis = plt.subplots(figsize=(7 * numbener, 7))
    
    listydat = empty((numbback + 2, numbener))
    listyerr = zeros((2, numbback + 2, numbener))
    
    listydat[0, :] = datafluxmean
    if post:
        listydat[1, :] = postpntsfluxmean[0, :]
        listyerr[:, 1, :] = retr_errrvarb(postpntsfluxmean)
        for c in indxback:
            listydat[c+2, :] = postnormback[0, c, :] * backfluxmean[c]
            listyerr[:, c+2, :] = retr_errrvarb(postnormback[:, c, :]) * backfluxmean[c]
    else:
        listydat[1, :] = mean(sum(thispntsflux * expo, 2) / sum(expo, 2), 1)
        for c in indxback:
            listydat[c+2, :] = thissampvarb[indxsampnormback[c, :]] * backfluxmean[c]

    
    xdat = meanener
    for k in range(numbback + 2):
        ydat = meanener**2 * listydat[k, :]
        yerr = meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5,                     ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    # Fermi-LAT results
    if trueinfo:
        if datatype == 'mock':
            
            pass
            
        else:
            
            if exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PNTS_TRAN_DATA_PATH"] + '/fermspec' + name + '.csv'
                    data = loadtxt(path)
                    enertemp = data[:, 0] # [GeV]
                    fluxtemp = data[:, 1] * 1e-3 # [GeV/cm^2/s/sr]
                    fluxtemp = interp(meanener, enertemp, fluxtemp)
                    #fluxtemp = interpolate.interp1d(enertemp, fluxtemp)(meanener)
                    axis.plot(meanener, fluxtemp, marker=listmrkr[k], color=listcolr[k], label=listlabl[k])


    axis.set_xlim([amin(binsener), amax(binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [GeV]')
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [GeV/cm$^2$/s/sr]')
    axis.legend()

    if post:
        path = plotpath + 'compfracspec' + rtag + '.png'
    else:
        path = plotpath + 'compfracspec' + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    
    listlabl = ['PS', 'Iso', 'FDM']
    listexpl = [0.1, 0, 0]

    listsize = zeros(numbback + 1)
    for k in range(numbback + 1):
        if numbener == 1:
            listsize[k] = diffener * listydat[k+1, :]
        else:
            listsize[k] = trapz(listydat[k+1, :], meanener)
    listsize *= 100. / sum(listsize)
        
    figr, axis = plt.subplots()

    axis.pie(listsize, explode=listexpl, labels=listlabl, autopct='%1.1f%%')
    axis.axis('equal')

    if post:
        path = plotpath + 'compfrac' + rtag + '.png'
    else:
        path = plotpath + 'compfrac' + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
     

    


# In[ ]:

def plot_histsind(l, postsindhist=None):
    
    if postsindhist == None:
        post = False
    else:
        post = True
        
    figr, axis = plt.subplots()
    if post:
        xdat = meansind[i, :]
        ydat = postsindhist[0, :]
        yerr = retr_errrvarb(postsindhist)
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    else:
        axis.hist(thissampvarb[thisindxsampsind[l]], binssind, alpha=0.5, color='b', log=True, label='Sample')
    if trueinfo:
        axis.hist(truesind[l], binssind, alpha=0.5, color='g', log=True, label=truelabl)
        if datatype == 'mock':
            axis.hist(fgl3sind, binssind, alpha=0.1, color='red', log=True, label='3FGL')
    axis.set_yscale('log')
    axis.set_xlabel('$s$')
    axis.set_xlim([minmsind, maxmsind])
    axis.set_ylabel('$N$')
    if post:
        path = plotpath + 'histsind_popl%d' % l + rtag + '.png'
    else:
        path = plotpath + 'histsind_popl%d' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    



# In[ ]:

def plot_histspec(l, listspechist=None):
    
    if listspechist == None:
        post = False
    else:
        post = True
        
    figr, axcl = plt.subplots(1, numbener, figsize=(7 * numbener, 7))
    if numbener == 1:
        axcl = [axcl]
    for i, axis in enumerate(axcl):
        if post:
            xdat = meanspec[i, :]
            yerr = empty((2, numbspec))
            ydat = percentile(listspechist[:, i, :], 50., axis=0)
            yerr[0, :] = percentile(listspechist[:, i, :], 16., axis=0)
            yerr[1, :] = percentile(listspechist[:, i, :], 84., axis=0)
            yerr = abs(yerr - ydat)
            axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        else:
            axis.hist(thissampvarb[thisindxsampspec[l]][i, :], binsspec[i, :], alpha=0.5, color='b',                     log=True, label='Sample')
            
            if not colrprio or i == indxenerfdfn:
                fdfnslop = thissampvarb[indxsampfdfnslop[l, i]]  
                fluxhistmodl = retr_fdfn(thissampvarb[indxsampfdfnnorm[l]], fdfnslop, i)
                axis.plot(meanspec[i, :], fluxhistmodl, ls='--', alpha=0.5, color='b')
            
        if trueinfo:
            truehist = axis.hist(truespec[l][0, i, :], binsspec[i, :],                                alpha=0.5, color='g', log=True, label=truelabl)
            if datatype == 'mock':
                axis.hist(fgl3spec[0, i, :], binsspec[i, :], color='red', alpha=0.1, log=True, label='3FGL')

        axis.set_yscale('log')
        axis.set_xlabel('$f$ ' + strgfluxunit)
        axis.set_xscale('log')
        axis.set_title(enerstrg[i])
        if trueinfo:
            axis.set_ylim([0.1, 1e3])
        axis.set_xlim([minmspec[i], maxmspec[i]])
        if i == 0:
            axis.set_ylabel('$N$')
        if i == numbener / 2:
            axis.legend()
        
    figr.subplots_adjust(wspace=0.3, bottom=0.2)
    
    if post:
        path = plotpath + 'histspec%d' % l + rtag + '.png'
    else:
        path = plotpath + 'histspec%d' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)
    


# In[ ]:


def plot_scatspec(l, postspecmtch=None, thisspecmtch=None):
    
    figr, axrw = plt.subplots(1, numbener, figsize=(7 * numbener, 6))
    if numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):

        xdat = truespec[l][0, i, :]
        xerr = retr_errrvarb(truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
 
        if postspecmtch != None or thisspecmtch != None:
            labl = '$f_{samp}$ ' + strgfluxunit
            axis.plot(meanspec[i, :], meanspec[i, :], ls='--', color='black', alpha=0.2)
            if postspecmtch != None:
                yerr[0, :] = postspecmtch[0, i, :]
                ydat = postspecmtch[1, i, :]
                yerr[1, :] = postspecmtch[2, i, :]
                yerr = abs(yerr - ydat)
            else:
                ydat = thisspecmtch[i, :]

        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
    
        if jtruepntstimevari[l].size > 0:
            axis.errorbar(xdat[jtruepntstimevari[l]], ydat[jtruepntstimevari[l]],                           ls='', yerr=yerr[:, jtruepntstimevari[l]],                           lw=1, marker='o', markersize=5, color='red')
            
    
        axis.set_xlabel('$f_{true}$ ' + strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_title(enerstrg[i])

        ylim = [minmspec[i], maxmspec[i]]
        
        axis.set_ylim(ylim)
        axis.set_xlim([minmspec[i], maxmspec[i]])
        axis.set_title(binsenerstrg[i])

    figr.subplots_adjust(wspace=0.4, bottom=0.2)

    if postspecmtch != None:
        path = plotpath + 'scatspec%d_' % l + rtag + '.png'
    elif thisspecmtch != None:
        path = plotpath + 'scatspec%d_' % l + rtag + '_%09d.png' % j

    plt.savefig(path)
    plt.close(figr)



# In[ ]:

def plot_scatpixl(l):
    
    figr, axgr = plt.subplots(numbevtt, numbener, figsize=(numbener * 7, numbevtt * 7))
    if numbevtt == 1:
        axgr = [axgr]
    for m, axrw in enumerate(axgr):
        if numbener == 1:
            axrw = [axrw]
        for i, axis in enumerate(axrw):
            
            slop, intc, coef, pval, stdv = sp.stats.linregress(datacnts[i, :, m], thismodlcnts[i, :, m])
            axis.scatter(datacnts[i, :, m], thismodlcnts[i, :, m], alpha=0.5)

            axislimt = [0., amax(datacnts[i, :, m]) * 1.5]
            axis.set_xlim(axislimt)
            axis.set_ylim(axislimt)
            
            axis.text(0.5, 0.9, '$m = %.3g, b = %.3g, r = %.3g$' % (slop, intc, coef),                     va='center', ha='center', transform=axis.transAxes, fontsize=16)
            if m == numbevtt - 1:
                axis.set_xlabel('Data Counts')
            if i == 0:
                labl = 'Model Counts'
                if exprtype == 'ferm':
                    labl += ', ' + evttstrg[m]
                axis.set_ylabel(labl)
            if m == 0:
                axis.set_title(enerstrg[i])


            
    figr.subplots_adjust(hspace=0.4, wspace=0.4, top=0.8)
    plt.savefig(plotpath + 'scatpixl%d_' % l + rtag + '_%09d.png' % j)
    plt.close(figr)
    
    


# In[ ]:


def plot_discspec(l, discspecmtch=None):
    
    figr, axrw = plt.subplots(1, numbener, figsize=(7 * numbener, 6))
    if numbener == 1:
        axrw = [axrw]
    for i, axis in enumerate(axrw):
        xdat = truespec[l][0, i, :]
        xerr = retr_errrvarb(truespec[l][:, i, :])
        yerr = zeros((2, xdat.size))
        if discspecmtch != None:
            ydat = discspecmtch[i, :]
            labl = '$f_{hit}$'
            ylim = [0., 1.]
        if postspecmtch != None or thisspecmtch != None:
            labl = '$f_{samp}$ ' + strgfluxunit
            axis.set_yscale('log')
            axis.plot(meanspec[i, :], meanspec[i, :], ls='--', color='black', alpha=0.2)
            if postspecmtch != None:
                yerr[0, :] = postspecmtch[0, i, :]
                ydat = postspecmtch[1, i, :]
                yerr[1, :] = postspecmtch[2, i, :]
                yerr = abs(yerr - ydat)
            else:
                ydat = thisspecmtch[i, :]
        axis.errorbar(xdat, ydat, ls='', yerr=yerr, lw=1, marker='o', markersize=5, color='black')
        axis.set_xlabel('$f_{true}$ ' + strgfluxunit)
        if i == 0:
            axis.set_ylabel(labl)
        axis.set_xscale('log')
        axis.set_title(enerstrg[i])
        if colrprio:
            if i == indxenerfdfn[0]:
                ylim = [minmspec[i], maxmspec[i]]
        else:
            ylim = [minmspec[i], maxmspec[i]]      
        axis.set_ylim(ylim)
        axis.set_xlim([minmspec[i], maxmspec[i]])
        axis.set_title(binsenerstrg[i])
    figr.subplots_adjust(wspace=0.4, bottom=0.2)
    path = plotpath + 'discspec%d_' % l + rtag + '_%09d.png' % j
    plt.savefig(path)
    plt.close(figr)


