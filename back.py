def proc_samp_back(gdat, gdatmodi, strg, raww=False):

    if gdatmodi != None:
        gdatobjt = gdatmodi
    else:
        gdatobjt = gdat

    if strg == 'true':
        strgtype = 'true'
        strgindx = 'true'
    else:
        strgtype = ''
        strgindx = 'this'
    
    
    proptrannext = (strg == 'next') and gdatmodi.proptran

    # grab the sample vector
    if proptrannext:
        sampvarb = getattr(gdatobjt, strgindx + 'sampvarb')
    else:
        sampvarb = getattr(gdatobjt, strg + 'sampvarb')


    if gdat.pntstype == 'lens':
        lgalsour = sampvarb[getattr(gdat, strgtype + 'indxfixplgalsour')]
        bgalsour = sampvarb[getattr(gdat, strgtype + 'indxfixpbgalsour')]
        specsour = sampvarb[getattr(gdat, strgtype + 'indxfixpspecsour')]
        sizesour = sampvarb[getattr(gdat, strgtype + 'indxfixpsizesour')]
        ellpsour = sampvarb[getattr(gdat, strgtype + 'indxfixpellpsour')]
        anglsour = sampvarb[getattr(gdat, strgtype + 'indxfixpanglsour')]
        lgalhost = sampvarb[getattr(gdat, strgtype + 'indxfixplgalhost')]
        bgalhost = sampvarb[getattr(gdat, strgtype + 'indxfixpbgalhost')]
        spechost = sampvarb[getattr(gdat, strgtype + 'indxfixpspechost')]
        if raww:
            beinhost = 0.
        else:
            beinhost = sampvarb[getattr(gdat, strgtype + 'indxfixpbeinhost')]
        ellphost = sampvarb[getattr(gdat, strgtype + 'indxfixpellphost')]
        anglhost = sampvarb[getattr(gdat, strgtype + 'indxfixpanglhost')]
        if raww:
            sherhost = 0.
        else:
            sherhost = sampvarb[getattr(gdat, strgtype + 'indxfixpsherhost')]
        sanghost = sampvarb[getattr(gdat, strgtype + 'indxfixpsanghost')]
       
        # host halo and external shear
        defl = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, sherhost, sanghost, 0.)
    
    if gdat.pntstype == 'lght':
        varioaxi = getattr(gdat, strgtype + 'varioaxi')
           



    if strg == 'next' and gdatmodi.proptran:
        if gdat.pntstype == 'lght':
            psfnintp = gdatmodi.thispsfnintp
        numbpntstemp = array([gdatmodi.numbpntsmodi])
        if gdat.pntstype == 'lens':
            defl = gdatmodi.thisdefl
        else:
            pntsflux = gdatmodi.thispntsflux
        numbpopl = 1
        lgalconc = gdatmodi.modilgal[:gdatmodi.numbpntsmodi]
        bgalconc = gdatmodi.modibgal[:gdatmodi.numbpntsmodi]
        specconc = gdatmodi.modispec[:, :gdatmodi.numbpntsmodi]
    
    else:
        indxpntsfull = list(getattr(gdatobjt, strgindx + 'indxpntsfull'))
        indxsamplgal, indxsampbgal, indxsampspec, indxsampspep, indxsampcompcolr = retr_indx(gdat, indxpntsfull)
        numbpntstemp = getattr(gdatobjt, strgindx + 'sampvarb')[gdat.indxfixpnumbpnts].astype(int)
        indxsamplgal = getattr(gdatobjt, strgindx + 'indxsamplgal')
        indxsampbgal = getattr(gdatobjt, strgindx + 'indxsampbgal')
        lgal = []
        bgal = []
        flux = []
        spec = []
        numbpopl = numbpntstemp.size
        for l in range(numbpopl):
            lgal.append(sampvarb[indxsamplgal[l]])
            bgal.append(sampvarb[indxsampbgal[l]])
            flux.append(sampvarb[indxsampspec[l]][gdat.indxenerfluxdist[0], :])
            spec.append(sampvarb[indxsampspec[l]])
        lgalconc = concatenate(lgal)
        bgalconc = concatenate(bgal)
        specconc = concatenate(spec, axis=1)
        if gdat.numbener > 1:
            indxsampspec = getattr(gdatobjt, strgindx + 'indxsampspec')
            spep = []
            for l in range(numbpopl):
                spep.append(sampvarb[indxsampspep[l]])
        
        # process a sample vector and the occupancy list to calculate secondary variables
	    ## secondary variables needed for likelihood evaluation    
        psfp = sampvarb[getattr(gdat, 'indxfixppsfp')]
        if gdat.pntstype == 'lens':
            psfnkern = AiryDisk2DKernel(psfp[0] / gdat.sizepixl)
            
        if gdat.pntstype == 'lght':
            ### PSF off-axis factor
            if varioaxi:
                oaxinorm = sampvarb[getattr(gdat, 'indxfixppsfpoaxinorm')]
                oaxiindx = sampvarb[getattr(gdat, 'indxfixppsfpoaxiindx')]
                factoaxi = retr_factoaxi(gdat, gdat.binsoaxi, oaxinorm, oaxiindx)
        
            psfntype = getattr(gdat, strgtype + 'psfntype')
            
            psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, varioaxi)
            if varioaxi:
                psfnintp = []
                for p in gdat.indxoaxi:
                    psfnintp.append(interp1d(gdat.binsangl, psfn[:, :, :, p], axis=1))
            else:
                psfnintp = interp1d(gdat.binsangl, psfn, axis=1)
            setattr(gdatobjt, strg + 'psfnintp', psfnintp)
    
    if strg == 'this' or strg == 'next' and not gdatmodi.proptran:
        setattr(gdatobjt, strg + 'lgal', lgal)
        setattr(gdatobjt, strg + 'bgal', bgal)
        setattr(gdatobjt, strg + 'flux', flux)
        setattr(gdatobjt, strg + 'spec', spec)
        if gdat.numbener > 1:
            setattr(gdatobjt, strg + 'spep', spep)
    
    bacp = sampvarb[getattr(gdat, 'indxfixpbacp')]

    if gdat.pntstype == 'lens':
        
        ## components
        if numbpntstemp > 0 and not raww:
            for k in range(numbpntstemp):
                # temp -- fix truncation
                defl += retr_defl(gdat, lgalconc[k], bgalconc[k], specconc[0, k], 0., 0., 0., 0., 0.)
                        
        # lensed image
        lensflux = retr_mapsraww(gdat, gdat.lgalgridcart - defl[:, :, 0], gdat.bgalgridcart - defl[:, :, 1], bacp, lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)

        # host emission
        mapshost = retr_mapssers(gdat, gdat.lgalgridcart, gdat.bgalgridcart, lgalhost, bgalhost, spechost, beinhost, ellphost, anglhost)
        
        # total emission
        modlfluxuncv = lensflux + bacp * gdat.backfluxcart[0][:, :, :, 0] + mapshost
        
        # convolve the lensed image with the PSF
        # temp
        if False:
            modlflux = empty_like(gdat.datacnts)
            for i in gdat.indxener:
                modlflux[i, :, :, 0] = convolve(modlflux[i, :, :, 0], psfnkern[i]).flatten()
        else:
            modlflux = modlfluxuncv.reshape((gdat.numbener, gdat.numbpixl, 1))
   
        if False:
            # temp
            mapssour = retr_mapsraww(gdat, gdat.lgalgridcart, gdat.bgalgridcart, bacp, lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
        
        setattr(gdatobjt, strg + 'defl', defl)
        
        if gdat.verbtype > 1:
            print 'modlflux'
            summgene(modlflux)
            print 'lensflux'
            summgene(lensflux)
            print 'defl'
            summgene(defl)

    if gdat.pntstype == 'lght':
        ### PS flux map
        pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, varioaxi, evalcirc=gdat.evalcirc)
        setattr(gdatobjt, strg + 'pntsflux', pntsflux)
        
        ### model flux map
        modlflux = retr_mapslght(gdat, bacp, pntsflux, gdat.indxcube)
   
    ### count map
    if gdat.pixltype != 'unbd':
        modlcnts = modlflux * gdat.expo * gdat.apix
        if gdat.enerbins:
            modlcnts *= gdat.diffener[:, None, None]
        setattr(gdatobjt, strg + 'modlcnts', modlcnts)
    
    if strg != 'true':
        
        ### log-prior
        if gdat.numbtrap > 0:
                
            if strg == 'next' and gdatmodi.proptran:
                if gdatmodi.propbrth or gdatmodi.propsplt:
                    fact = -1.
                else:
                    fact = 1.
                gdatmodi.templpri[0] = gdatmodi.thislpri[0] + fact * gdat.priofactdoff * gdat.numbcompcolr[gdatmodi.indxpoplmodi]
                gdatmodi.templpri[1:] = gdatmodi.thislpri[1:]
            else:
                for l in gdat.indxpopl:
                    gdatmodi.templpri[0] = gdat.priofactdoff * gdat.numbcompcolr[l] * numbpntstemp[l]
                    
                    gdatmodi.templpri[1+l] = retr_probpois(sampvarb[gdat.indxfixpnumbpnts[l]], sampvarb[gdat.indxfixpmeanpnts[l]])
                    
                    if gdat.fluxdisttype[l] == 'powr':
                        gdatmodi.tempfluxdistslop[l] = sampvarb[gdat.indxfixpfluxdistslop[l]]
                        lpriflux = sum(log(pdfn_flux_powr(gdat, spec[l][gdat.indxenerfluxdist[0], :], gdatmodi.tempfluxdistslop[l])))
                    if gdat.fluxdisttype[l] == 'brok':
                        gdatmodi.tempfluxdistbrek[l] = sampvarb[gdat.indxfixpfluxdistbrek[l]]
                        gdatmodi.tempfluxdistsloplowr[l] = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                        gdatmodi.tempfluxdistslopuppr[l] = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                        lpriflux = sum(log(pdfn_flux_brok(gdat, spec[l][gdat.indxenerfluxdist[0], :], gdatmodi.tempfluxdistbrek[l], \
                                                                                                                gdatmodi.tempfluxdistsloplowr[l], gdatmodi.tempfluxdistslopuppr[l])))
                    gdatmodi.templpri[1+gdat.numbpopl+l] = lpriflux
                    
                    if gdat.numbener > 1:
                        gdatmodi.tempsinddistmean[l] = sampvarb[gdat.indxfixpsinddistmean[l]]
                        gdatmodi.tempsinddiststdv[l] = sampvarb[gdat.indxfixpsinddiststdv[l]]
                        gdatmodi.templpri[1+2*gdat.numbpopl+l] = sum(lpdf_gaus(spep[l][:, 0], gdatmodi.tempsinddistmean[l], gdatmodi.tempsinddiststdv[l])) 
            
            setattr(gdatobjt, strg  + 'fluxdistslop', gdatmodi.tempfluxdistslop)
            setattr(gdatobjt, strg  + 'fluxdistslop', gdatmodi.tempfluxdistslop)
            setattr(gdatobjt, strg  + 'fluxdistbrek', gdatmodi.tempfluxdistbrek)
            setattr(gdatobjt, strg  + 'fluxdistsloplowr', gdatmodi.tempfluxdistsloplowr)
            setattr(gdatobjt, strg  + 'fluxdistslopuppr', gdatmodi.tempfluxdistslopuppr)
            if gdat.numbener > 1:
                setattr(gdatobjt, strg  + 'sinddistmean', gdatmodi.tempsinddistmean)
                setattr(gdatobjt, strg  + 'sinddiststdv', gdatmodi.tempsinddiststdv)
        
        ### log-likelihood
        if gdat.pixltype == 'unbd':
            gdatmodi.templlik = gdat.numbdatasamp * log(modlfluxtotl) - modlfluxtotl + log(modlflux)
        else:
            if gdat.liketype == 'pois':
                gdatmodi.templlik = gdat.datacnts * log(modlcnts) - modlcnts
            if gdat.liketype == 'gaus':
                gdatmodi.templlik = -0.5 * (gdat.datacnts - modlcnts)**2 / gdat.datacnts
            
        if strg == 'next':
            print 'hey'
            print 'gdatmodi.templlik'
            summgene(gdatmodi.templlik)
            print 'gdatmodi.thisllik'
            summgene(gdatmodi.thisllik)
            print
            deltlpri = sum(gdatmodi.templpri - gdatmodi.thislpri)
            deltllik = sum(gdatmodi.templlik - gdatmodi.thisllik)
            setattr(gdatmodi, 'deltlpri', deltlpri) 
            setattr(gdatmodi, 'deltllik', deltllik)
            
            # temp
            laccfact = 0.
            setattr(gdatmodi, 'laccfact', laccfact)
        setattr(gdatmodi, strg + 'lpri', gdatmodi.templpri) 
        setattr(gdatmodi, strg + 'gdatmodi.templlik', gdatmodi.templlik) 
    
    # mock data specific
    if strg == 'true':
        if raww:
            strgvarb = 'truemodlcntsraww'
        else:
            strgvarb = 'datacnts'
        
        # generate count data
        cntstemp = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxevtt:
                    cntstemp[i, j, m] = poisson(modlcnts[i, j, m])
        setattr(gdat, strgvarb, cntstemp)
        
        if raww:
            return

    ## tertiary variables that are not needed for likelihood evaluation
    if strg != 'next':
       
        setattr(gdatobjt, strg + 'psfp', psfp)
        if gdat.pixltype != 'unbd':
            resicnts = gdat.datacnts - modlcnts
            setattr(gdatobjt, strg + 'resicnts', resicnts)
        if gdat.pntstype == 'lght':
            if varioaxi:
                setattr(gdatobjt, strg + 'factoaxi', factoaxi)
            
            setattr(gdatobjt, strg + 'psfn', psfn)

            ### PSF FWHM
            fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strg + 'fwhm', fwhm)
            
	    	### mean PS flux map 
            pntsfluxmean = sum(sum(pntsflux * gdat.expo, 2), 1) / sum(sum(gdat.expo, 2), 1)
            setattr(gdatobjt, strg + 'pntsfluxmean', pntsfluxmean)

            if gdat.pixltype != 'unbd':
                ### number of background counts per PSF
                cntsbackfwhm = retr_cntsbackfwhm(gdat, bacp, fwhm)
            
                ### number of counts and standard deviation of each PS
                cnts = []
                sigm = []
                for l in gdat.indxpopl:
                    # temp -- zero exposure pixels will give zero counts
                    indxpixltemp = retr_indxpixl(gdat, bgal[l], lgal[l])
                    cntstemp = sampvarb[indxsampspec[l]][:, :, None] * gdat.expofull[:, indxpixltemp, :]
                    if gdat.enerbins:
                        cntstemp *= gdat.diffener[:, None, None]
                    cnts.append(cntstemp)
                    if gdat.varioaxi:
                        sigmtemp = retr_sigm(gdat, cntstemp, cntsbackfwhm, lgal=lgal[l], bgal=bgal[l])
                    else:
                        sigmtemp = retr_sigm(gdat, cntstemp, cntsbackfwhm)
                    sigm.append(sigmtemp)
                
            if gdat.calcerrr and gdat.numbtrap > 0:
                pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, gdat.varioaxi, evalcirc=False)
                pntscnts = pntsflux * gdat.expo * gdat.apix
                if gdat.enerbins:
                    pntscnts *= gdat.diffener[:, None, None]
                errrcnts = pntscnts - temppntscnts
                errr = zeros_like(errrcnts)
                indxcubegood = where(temppntscnts > 1e-10)
                errr[indxcubegood] = 100. * errrcnts[indxcubegood] / temppntscnts[indxcubegood]
                setattr(gdatobjt, strg + 'errrcnts', errrcnts)
                setattr(gdatobjt, strg + 'errr', errr)
                if False and amax(fabs(errr)) > 0.1:
                    raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')

        if gdat.pntstype == 'lens':

            setattr(gdatobjt, strg + 'psfnkern', psfnkern)
            setattr(gdatobjt, strg + 'modlfluxuncv', modlfluxuncv)
            setattr(gdatobjt, strg + 'lensflux', lensflux)
            setattr(gdatobjt, strg + 'mapshost', mapshost)

            ### deflection
            #### number of deflection components
            numbdeflpnts = min(3, numbpntstemp)
            if numbpntstemp > 0:
                indxpntssortbrgt = argsort(spec[0][0, :])[::-1]
                lgalsort = lgal[0][indxpntssortbrgt][:numbdeflpnts]
                bgalsort = bgal[0][indxpntssortbrgt][:numbdeflpnts]
                beinsort = spec[0][0, indxpntssortbrgt][:numbdeflpnts]
            
            numbdeflsing = numbdeflpnts + 2
            deflsing = empty((gdat.numbsidecart, gdat.numbsidecart, 2, numbdeflsing))
            for k in range(numbdeflsing):
                if k == 0:
                    deflsing[:, :, :, k] = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, 0., 0., 0.)
                elif k == 1:
                    deflsing[:, :, :, k] = retr_defl(gdat, 0., 0., 0., 0., 0., sherhost, sanghost, 0.)
                else:
                    deflsing[:, :, :, k] = retr_defl(gdat, lgalsort[k-2], bgalsort[k-2], beinsort[k-2], 0., 0., 0., 0., 0.)

            ### convergence
            conv = retr_conv(gdat, defl) 
            convpsec = retr_psec(gdat, conv)
            convpsecodim = retr_psecodim(gdat, convpsec) 
            histdefl = histogram(defl, bins=gdat.binsdefl)[0]
            setattr(gdatobjt, strg + 'conv', conv)
            setattr(gdatobjt, strg + 'convpsec', convpsec)
            setattr(gdatobjt, strg + 'convpsecodim', convpsecodim)
            setattr(gdatobjt, strg + 'numbdeflsing', numbdeflsing)
            setattr(gdatobjt, strg + 'histdefl', histdefl)
            setattr(gdatobjt, strg + 'deflsing', deflsing)
     
	    ### PS indices to compare with the reference catalog
        if strg == 'this' and gdat.numbtrap > 0 and gdat.trueinfo:
            for l in gdat.indxpopl:
                gdatmodi.indxmodlpntscomp[l] = retr_indxpntscomp(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]], gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]])
            gdatmodi.trueindxpntsassc = []
            gdatmodi.thisspecassc = []
            for l in gdat.indxpopl:
                indxmodl, trueindxpntsassc = corr_catl(gdat, gdatmodi, l, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]], gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]], \
                                                                                                                            gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]])
                gdatmodi.trueindxpntsassc.append(trueindxpntsassc)
                gdatmodi.thisspecassc.append(zeros((gdat.numbener, gdat.truefixp[gdat.indxfixpnumbpnts[l]])))
                temp = where(indxmodl >= 0)[0]
                gdatmodi.thisspecassc[l][:, temp] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l]][:, indxmodl[temp]]
                gdatmodi.thisspecassc[l][:, gdatmodi.trueindxpntsassc[l].miss] = 0.

            if gdat.pntstype == 'lens' and gdat.trueinfo and gdat.datatype == 'mock':
                gdatmodi.thisdeflresi = gdatmodi.thisdefl - gdat.truedefl
                gdatmodi.thisdeflcomp = sum(gdatmodi.thisdefl * gdat.truedefl, 2) / sqrt(gdatmodi.thisdefl[:, :, 0]**2 + gdatmodi.thisdefl[:, :, 1]**2) / \
                                                                                                           sqrt(gdat.truedefl[:, :, 0]**2 + gdat.truedefl[:, :, 1]**2)
            
        if gdatmodi != None:
            ### corrected prior
            gdatmodi.thislprinorm = 0.
            for l in gdat.indxpopl:
                # temp -- brok terms are not complete
                break
                numbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]]
                meanpnts = gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[l]]
                gdatmodi.thislprinorm += numbpnts * gdat.priofactlgalbgal + gdat.priofactfluxdistslop + gdat.priofactmeanpnts - log(meanpnts)
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist[0], :]]
                if gdat.fluxdisttype[l] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]]
                    gdatmodi.thislprinorm -= log(1. + fluxdistslop**2)
                    gdatmodi.thislprinorm += sum(log(pdfn_flux_powr(gdat, flux, fluxdistslop)))
                if gdat.fluxdisttype[l] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[l]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                    gdatmodi.thislprinorm += sum(log(pdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))
           
        	# temp 
            if False:
                numbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[0]]
                diffllikdiffpara = empty(numbpnts)
                for k in range(numbpnts):
                    diffllikdiffpara[k]
                gdatmodi.listdiffllikdiffpara.append(diffllikdiffpara)
        
                tranmatr = diffllikdiffpara[:, None] * gdatmodi.listdiffllikdiffpara[gdatmodi.cntrswep-1][None, :]
                gdatmodi.listtranmatr.append(tranmatr)



