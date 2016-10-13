# common imports
from __init__ import *

def retr_psfnwdth(gdat, psfn, frac):
    
    wdth = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            indxanglgood = argsort(psfn[i, :, m])
            intpwdth = max(frac * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpwdth > amin(psfn[i, indxanglgood, m]) and intpwdth < amax(psfn[i, indxanglgood, m]):
                wdth[i, m] = interp1d(psfn[i, indxanglgood, m], gdat.binsangl[indxanglgood])(intpwdth)
    return wdth


def retr_spec(gdat, flux, sind):
    
    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])

    spec = flux[None, :] * (gdat.meanener[:, None] / gdat.enerfluxdist)**(-sind[None, :])
    
    return spec


def retr_indx(gdat, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampsind = []
    indxsampcomp = []
    for l in gdat.indxpopl:
        indxsamplgaltemp = gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:l]) + array(indxpntsfull[l], dtype=int) * gdat.numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], gdat.numbener, 0) +repeat(arange(gdat.numbener), len(indxpntsfull[l])).reshape(gdat.numbener, -1))
        indxsampsind.append(indxsamplgaltemp + 2 + gdat.numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, gdat.numbcomp) + tile(arange(gdat.numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(gdat, lgal, bgal, spec, psfipara, psfntype):
    
    numbpnts = lgal.size
    pntsfluxsing = empty((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for k in range(numbpnts):
        
        # calculate the distance to all pixels from each point source
        dist = retr_angldistunit(gdat, lgal[k], bgal[k], gdat.indxpixl)

        indx = argsort(dist)
        dist = dist[indx]
        indxpixltemp = gdat.indxpixl[indx]
        
        # evaluate the PSF
        psfn = retr_psfn(gdat, psfipara, gdat.indxener, dist, psfntype)
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                pntsfluxsing[k, i, indxpixltemp, m] = spec[i, k] * psfn[i, :, m]

    # sum contributions from all PS
    pntsflux = sum(pntsfluxsing, 0) 

    return pntsflux


def retr_rofi_flux(gdat, normback, pntsflux, tempindx):
    
    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += normback[c, :, None, None] * gdat.backflux[c][tempindx]        
    
    return modlflux


def cdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):

    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (gdat.maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunit = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (flux**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
    indxflux = where(flux >= fluxdistbrek)[0]
    
    if indxflux.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
        fluxunit[indxflux] = temp + norm / (1. - fluxdistslopuppr) * fluxdistbrek**fluxdistslopuppr * \
                                                                        (flux[indxflux]**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr))

    return fluxunit


def pdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):

    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (gdat.maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    pdfn = norm * (flux / fluxdistbrek)**(-fluxdistsloplowr)
    indxflux = where(flux >= fluxdistbrek)[0]
    
    if indxflux.size > 0:
        pdfn[indxflux] = norm * (flux[indxflux] / fluxdistbrek)**(-fluxdistslopuppr)
        
    return pdfn


def icdf_flux_brok(gdat, fluxunit, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):
   
    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (gdat.maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunitbrek = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
    flux = (fluxunit * (1. - fluxdistsloplowr) / norm / fluxdistbrek**fluxdistsloplowr + gdat.minmflux**(1. - fluxdistsloplowr))**(1. / (1. - fluxdistsloplowr))
    indxfluxunit = where(fluxunit >= fluxunitbrek)[0]
    
    if indxfluxunit.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
        flux[indxfluxunit] = ((fluxunit[indxfluxunit] - temp) * (1. - fluxdistslopuppr) / norm / fluxdistbrek**fluxdistslopuppr + \
                                                                                            fluxdistbrek**(1. - fluxdistslopuppr))**(1. / (1. - fluxdistslopuppr))

    return flux


def cdfn_flux_powr(gdat, flux, fluxdistslop):
        
    fluxunit = (flux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) / (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop))
        
    return fluxunit


def icdf_flux_powr(gdat, fluxunit, fluxdistslop):

    # temp
    try:
        flux = (fluxunit * (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) + gdat.minmflux**(1. - fluxdistslop))**(1. / (1. - fluxdistslop))
    except:
        print '(fluxunit * (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) + gdat.minmflux**(1. - fluxdistslop))'
        print (fluxunit * (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) + gdat.minmflux**(1. - fluxdistslop))
        print '(1. / (1. - fluxdistslop))'
        print (1. / (1. - fluxdistslop))
        raise

    return flux


def pdfn_flux_powr(gdat, flux, fluxdistslop):
  
    norm = (1. - fluxdistslop) / (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop))
    
    pdfn = norm * flux**(-fluxdistslop)
    
    return pdfn


def icdf_self(paraunit, minmpara, factpara):
    para = factpara * paraunit + minmpara
    return para


def cdfn_self(para, minmpara, factpara):
    paraunit = (para - minmpara) / factpara
    return paraunit


def icdf_eerr(paraunit, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    cdfnnormpara = paraunit * cdfnnormdiff + cdfnnormminm
    tranpara = sp.special.erfinv(2. * cdfnnormpara - 1.) * sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def cdfn_eerr(para, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    tranpara = (para - meanpara) / stdvpara
    cdfnnormpara = 0.5 * (sp.special.erf(tranpara / sqrt(2.)) + 1.)
    paraunit = (cdfnnormpara - cdfnnormminm) / cdfnnormdiff
    return paraunit


def icdf_logt(paraunit, minmpara, factpara):
    para = exp(paraunit * factpara) * minmpara
    return para


def cdfn_logt(para, minmpara, factpara):
    paraunit = log(para / minmpara) / factpara
    return paraunit


def icdf_atan(paraunit, minmpara, factpara):
    para = tan(factpara * paraunit + arctan(minmpara))
    return para


def cdfn_atan(para, minmpara, factpara):
    paraunit = (arctan(para) - arctan(minmpara)) / factpara
    return paraunit


def icdf_psfipara(gdat, psfiparaunit, thisindxpsfipara):

    minmpsfipara = gdat.minmpsfipara[thisindxpsfipara]
    factpsfipara = gdat.factpsfipara[thisindxpsfipara]
    scalpsfipara = gdat.scalpsfipara[thisindxpsfipara]
        
    if scalpsfipara == 'self':
        psfipara = icdf_self(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfipara = icdf_logt(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfipara = icdf_atan(psfiparaunit, minmpsfipara, factpsfipara)

    return psfipara


def cdfn_psfipara(gdat, psfipara, thisindxpsfipara):
    
    minmpsfipara = gdat.minmpsfipara[thisindxpsfipara]
    factpsfipara = gdat.factpsfipara[thisindxpsfipara]
    scalpsfipara = gdat.scalpsfipara[thisindxpsfipara]

    if scalpsfipara == 'self':
        psfiparaunit = cdfn_self(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfiparaunit = cdfn_logt(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfiparaunit = cdfn_atan(psfipara, minmpsfipara, factpsfipara)
    
    return psfiparaunit


def retr_thisindxprop(gdat, gdatmodi):
    
    # choose the population to be modified
    gdatmodi.indxpoplmodi = choice(gdat.indxpopl)
    
    numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
    if numbpnts == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probpropmaxm)
    elif numbpnts == gdat.minmnumbpnts:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probpropminm)
    else:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probprop)

    if gdat.verbtype > 1:
        print inspect.stack()[0][3]
        print 'indxpoplmodi'
        print gdatmodi.indxpoplmodi
        print 'thisnumbpnts'
        print numbpnts
        print 'maxmnumbpnts'
        print gdat.maxmnumbpnts
        print


def retr_indxpixl(gdat, bgal, lgal):

    if gdat.pixltype == 'heal':
        indxpixl = gdat.pixlcnvt[ang2pix(gdat.numbsideheal, pi / 2. - bgal, lgal)]
        if gdat.diagmode:
            if (indxpixl == -1).any():  
                print 'pixlcnvt went negative!'
                raise Exception

    if gdat.pixltype == 'cart':
        indxlgcr = floor(gdat.numbsidecart * (lgal - gdat.minmlgal) / 2. / gdat.maxmgang).astype(int)
        indxbgcr = floor(gdat.numbsidecart * (bgal - gdat.minmbgal) / 2. / gdat.maxmgang).astype(int)

        if isscalar(indxlgcr):
            if indxlgcr < 0:
                indxlgcr = 0
            if indxlgcr >= gdat.numbsidecart:
                indxlgcr = gdat.numbsidecart - 1
        else:
            indxlgcr[where(indxlgcr < 0)] = 0
            indxlgcr[where(indxlgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        if isscalar(indxbgcr):
            if indxbgcr < 0:
                indxbgcr = 0
            if indxbgcr >= gdat.numbsidecart:
                indxbgcr = gdat.numbsidecart - 1
        else:
            indxbgcr[where(indxbgcr < 0)] = 0
            indxbgcr[where(indxbgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        indxpixl = indxlgcr * gdat.numbsidecart + indxbgcr

    return indxpixl


def retr_elpsfrac(elpsaxis):
    
    distnorm = sum(((listsamp - gdat.elpscntr[None, :]) / elpsaxis[None, :])**2, axis=1)
    indxsampregu = where(distnorm < 1.)[0]
    thissampfrac = indxsampregu.size / gdat.numbsamp
    vari = (thissampfrac / 0.05 - 1.)**2
    
    return vari


def retr_llik(gdat, gdatmodi, init=False):

    if init:
        gdatmodi.thisllik = gdat.datacnts * log(gdatmodi.thismodlcnts) - gdatmodi.thismodlcnts
        gdatmodi.thislliktotl = sum(gdatmodi.thisllik)
        
    elif gdatmodi.thisindxprop >= gdat.indxproppsfipara:

        # load convenience variables
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            lgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsamplgal)]
            bgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampbgal)]
            spec = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdatmodi.indxenermodi, :]]
        if gdatmodi.thisindxprop >= gdat.indxpropbrth:
            lgal = gdatmodi.modilgal[:gdatmodi.numbmodipnts]
            bgal = gdatmodi.modibgal[:gdatmodi.numbmodipnts]
            spec = gdatmodi.modispec[meshgrid(gdatmodi.indxenermodi, arange(gdatmodi.numbmodipnts), indexing='ij')]
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timeinit
        
        # determine pixels over which to evaluate the log-likelihood
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            gdatmodi.indxpixlmodi = gdat.indxpixl
        if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
            thisindxpixlprox = []
            for k in range(gdatmodi.numbmodipnts):
                # temp -- this may not work for extreme color PS!
                # take the flux at the pivot energy
                if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                    fluxtemp = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenerfluxdist, k]]
                else:
                    fluxtemp = gdatmodi.modispec[gdat.indxenerfluxdist, k]
                # find the flux index
                indxfluxproxtemp = amin(where(gdat.binsfluxprox - fabs(fluxtemp) > 0.)[0]) - 1
                indxpixltemp = retr_indxpixl(gdat, bgal[k], lgal[k])
                thisindxpixlprox.append(gdat.indxpixlprox[indxfluxproxtemp][indxpixltemp])
            gdatmodi.indxpixlmodi = unique(concatenate(thisindxpixlprox))
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timeinit

        # construct the mesh grid for likelihood evaluation
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            gdatmodi.indxcubemodi = meshgrid(gdatmodi.indxenermodi, gdatmodi.indxpixlmodi, gdat.indxevtt, indexing='ij')
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 2] = timefinl - timeinit

        # update the model PS flux map, if needed
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            # temp
            #if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            #    gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] = 0.
            #else:
            gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] = copy(gdatmodi.thispntsflux[gdatmodi.indxcubemodi])
                
            # evaluate the PSF for each PS over the set of data pixels to be updated
            # temp
            if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                numbrept = 2
            else:
                numbrept = 1
            for n in range(numbrept):
                # grab the PSF
                if n == 0:
                    psfnintp = gdatmodi.thispsfnintp
                else:
                    psfnintp = gdatmodi.nextpsfnintp
                for k in range(gdatmodi.numbmodipnts):
                    # calculate the distance to the pixels to be updated
                    dist = retr_angldistunit(gdat, lgal[k], bgal[k], thisindxpixlprox[k])
                    # interpolate the PSF
                    psfn = psfnintp(dist)
                    # add the contribution of the PS to the the proposed flux map
                    for i in range(gdatmodi.indxenermodi.size):
                        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                            if n == 0:
                                spectemp = -spec[i, k]
                            else:
                                spectemp = spec[i, k]
                        else:
                            spectemp = spec[i, k]
                        gdatmodi.nextpntsflux[gdatmodi.indxenermodi[i], thisindxpixlprox[k], :] += spectemp * psfn[gdatmodi.indxenermodi[i], :, :]
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 3] = timefinl - timeinit
        
        # update the total model flux map
        timeinit = gdat.functime()
        indxtemp = meshgrid(gdat.indxback, gdatmodi.indxenermodi, indexing='ij')
        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            normback = gdatmodi.nextsampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.thispntsflux
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            normback = gdatmodi.thissampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.nextpntsflux
        gdatmodi.nextmodlflux[gdatmodi.indxcubemodi] = retr_rofi_flux(gdat, normback, pntsflux, gdatmodi.indxcubemodi)
        
        
        # temp
        if True and gdat.strgcnfg == 'test_uppr':
            
            temppntsflux, temppntscnts, tempmodlflux, tempmodlcnts = retr_maps(gdat, list(gdatmodi.thisindxpntsfull), copy(gdatmodi.thissampvarb))
            gdatmodi.thispntscnts = gdatmodi.thispntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
            gdatmodi.thiserrrpnts = gdatmodi.thispntscnts - temppntscnts
            
            facttemp = gdat.expo[gdatmodi.indxcubemodi].flatten() * gdat.diffener[gdatmodi.indxenermodi] * gdat.apix

            print 'facttemp'
            print facttemp.shape
            print amin(facttemp)
            print amax(facttemp)

            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = temppntsflux[gdatmodi.indxcubemodi] * facttemp
            path = gdat.pathdiag + 'temppntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
            print 'temppntsflux[gdatmodi.indxcubemodi]'
            print temppntsflux[gdatmodi.indxcubemodi].shape
            print amin(temppntsflux[gdatmodi.indxcubemodi])
            print amax(temppntsflux[gdatmodi.indxcubemodi])

            print 'temp[gdatmodi.indxpixlmodi]'
            print temp[gdatmodi.indxpixlmodi].shape
            print amin(temp[gdatmodi.indxpixlmodi])
            print amax(temp[gdatmodi.indxpixlmodi])

            
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = (gdatmodi.thispntscnts[gdatmodi.indxcubemodi] - temppntscnts[gdatmodi.indxcubemodi]) * facttemp
            path = gdat.pathdiag + 'errrpntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, resi=True, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)

            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.thispntsflux[gdatmodi.indxcubemodi] * facttemp
            path = gdat.pathdiag + 'thispntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
        
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] * facttemp
            path = gdat.pathdiag + 'nextpntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
        
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = (gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] - gdatmodi.thispntsflux[gdatmodi.indxcubemodi]) * facttemp
            path = gdat.pathdiag + 'diffpntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
            
            if amax(fabs(gdatmodi.thiserrrpnts)) > 0.:
                print 
                raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')

            
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 4] = timefinl - timeinit

        # calculate the total model count map
        timeinit = gdat.functime()
        gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi] = gdatmodi.nextmodlflux[gdatmodi.indxcubemodi] * gdat.expo[gdatmodi.indxcubemodi] * \
                                                                                                                gdat.apix * gdat.diffener[gdatmodi.indxenermodi, None, None] # [1]
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 5] = timefinl - timeinit

        # calculate the likelihood
        timeinit = gdat.functime()
        gdatmodi.nextllik[gdatmodi.indxcubemodi] = gdat.datacnts[gdatmodi.indxcubemodi] * log(gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi]) \
                                                                                                                - gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi]
        gdatmodi.deltllik = sum(gdatmodi.nextllik[gdatmodi.indxcubemodi] - gdatmodi.thisllik[gdatmodi.indxcubemodi])
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 6] = timefinl - timeinit

        if gdat.diagmode:
            if not isfinite(gdatmodi.nextllik[gdatmodi.indxcubemodi]).any():
                warnings.warn('Log-likelihood went NAN!')
            
    else:
        gdatmodi.deltllik = 0.
        
    
def retr_backfwhmcnts(gdat, normback, fwhm):

    backfwhmcnts = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:
        backfwhmcnts += normback[c, :, None, None] * gdat.backflux[c] * gdat.expo * gdat.diffener[:, None, None] * pi * fwhm[:, None, :]**2 / 4.

    return backfwhmcnts


def retr_sigm(gdat, cnts, backfwhmcnts):
    
    sigm = cnts / sum(mean(backfwhmcnts, 1), 1)[:, None]

    return sigm


def retr_fluxlpribind(gdatmodi, gdat, l):
    
    if gdat.fluxdisttype[l] == 'powr':
        fluxhistmodl = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]] * gdat.diffflux * pdfn_flux_powr(gdat, gdat.meanflux, \
                                                                               gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]])
    if gdat.fluxdisttype[l] == 'brok':
        fluxhistmodl = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]] * gdat.diffflux * pdfn_flux_brok(gdat, gdat.meanflux, \
                                                                               gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]], \
                                                                               gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]], \
                                                                               gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]])
    fluxhistdata = histogram(thisflux, gdat.binsflux)[0]
    fluxlpribind = retr_probpois(fluxhistdata, fluxhistmodl)
    
    return fluxlpribind


def retr_fluxlpri(gdatmodi, gdat, l):
    
    if gdat.fluxdisttype[l] == 'powr':
        fluxlpri = sum(log(pdfn_flux_powr(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]], \
                                                gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]])))
    if gdat.fluxdisttype[l] == 'brok':
        fluxlpri = sum(log(pdfn_flux_brok(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]], \
                                                gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]], \
                                                gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]], \
                                                gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]])))
    
    return fluxlpri


def retr_probpois(data, modl):
    
    prob = data * log(modl) - modl - sp.special.gammaln(data + 1)

    return prob


def retr_lpri(gdat, gdatmodi, init=False):
        
    if init:
        for l in gdat.indxpopl:
            thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]]
            if gdat.bindprio:
                gdatmodi.thislpri[l, :] = retr_fluxlpribind(gdatmodi, gdat, l)
            else:
                gdatmodi.thislpri[l, 0] = retr_probpois(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[l]], gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]])
                gdatmodi.thislpri[l, 1] = retr_fluxlpri(gdatmodi, gdat, l)
        gdatmodi.thislpritotl = sum(gdatmodi.thislpri)
        gdatmodi.nextlpri = copy(gdatmodi.thislpri)
    else:
        
        # initialize the delta log-prior
        gdatmodi.deltlpri = 0.

        # determine which type of delta log-prior is to be calculated
        boolupdtmeanpnts = gdatmodi.thisindxprop == gdat.indxpropmeanpnts
        boolupdtfluxdist = gdatmodi.thisindxprop == gdat.indxpropfluxdistslop or gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or \
                                                                                 gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or \
                                                                                 gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr
        boolupdtnumbpnts = gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg

        # calculate contributions to the delta log-prior
        if boolupdtmeanpnts or boolupdtfluxdist or boolupdtnumbpnts:

            # penalty term due to the number of degrees of freedom
            if boolupdtnumbpnts:
                if gdatmodi.thisindxprop == gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxpropsplt:
                    deltdoff = -gdat.numbcompcolr
                else:
                    deltdoff = gdat.numbcompcolr
                gdatmodi.deltlpri += gdat.priofactdoff * deltdoff

            # binned flux prior
            # temp -- binned prior currently does not work with splits and merges!
            if gdat.bindprio:
               
                # mean number of PS
                if boolupdtmeanpnts:
                    meanpnts = gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
                else:
                    meanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
                
                # hyperparameters on the flux distribution
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        fluxdistslop = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                            fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                            fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                            fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                
                # number of PS
                if boolupdtnumbpnts:
                    fluxhist = histogram(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :]], gdat.binsflux)[0] 
                    if gdatmodi.thisindxprop == gdat.indxpropbrth:
                        fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                    elif gdatmodi.thisindxprop == gdat.indxpropdeth:
                        fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                    elif gdatmodi.thisindxprop == gdat.indxpropsplt:
                        fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                        fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 1:3], gdat.binsflux)[0]
                    elif gdatmodi.thisindxprop == gdat.indxpropmerg:
                        fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, :2], gdat.binsflux)[0]
                        fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 2], gdat.binsflux)[0]
                
                # construct the flux prior
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxhistmodl = meanpnts * gdat.diffflux * pdfn_flux_powr(gdat, gdat.meanflux, fluxdistslop)
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxhistmodl = meanpnts * gdat.diffflux * pdfn_flux_brok(gdat, gdat.meanflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
            
                gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] = retr_probpois(fluxhist, fluxhistmodl)
            
                gdatmodi.deltlpri = sum(gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, :])

            # unbinned flux prior
            else:
    
                if boolupdtnumbpnts or boolupdtmeanpnts:
                    if boolupdtnumbpnts:
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = retr_probpois(gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]])
                    else:
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = retr_probpois(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]])
                    gdatmodi.deltlpri += gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 0]
                    
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_powr(gdat, \
                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]], \
                                                                    gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])))
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                        fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                            fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                            fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                            fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_brok(gdat, \
                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]], \
                                                    fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))
                    gdatmodi.deltlpri += gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 1]
       
            if gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg:
                
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    lprbfrst = log(pdfn_flux_powr(gdat, gdatmodi.fluxfrst, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]))
                    lprbseco = log(pdfn_flux_powr(gdat, gdatmodi.fluxseco, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]))
                    lprbpare = log(pdfn_flux_powr(gdat, gdatmodi.fluxpare, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]))
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    lprbfrst += log(pdfn_flux_brok(gdat, gdatmodi.fluxfrst, \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
                    
                    lprbseco += log(pdfn_flux_brok(gdat, gdatmodi.fluxseco, \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]))

                    lprbpare -= log(pdfn_flux_brok(gdat, gdatmodi.fluxpare, \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
                if gdatmodi.thisindxprop == gdat.indxpropsplt:
                    gdatmodi.deltlpri += lprbfrst
                    gdatmodi.deltlpri += lprbseco
                    gdatmodi.deltlpri -= lprbpare
                else:
                    gdatmodi.deltlpri += lprbpare
                    gdatmodi.deltlpri -= lprbfrst
                    gdatmodi.deltlpri -= lprbseco
                    
                    # temp
                    ## split
                    # P(f1)P(l1)P(b1)P(s1)P(f2)P(l2)P(b2)P(s2) / P(f0)P(l0)P(b0)P(s0)P(uf)P(ur)P(up)P(us)
                    # P(f1)P(f2)P(l2)P(b2) / P(f0)P(uf)P(ur)P(up)
                    # P(f1)P(f2) / P(f0)

        
def retr_sampvarb(gdat, indxpntsfull, samp):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxsampnumbpnts] = samp[gdat.indxsampnumbpnts]
    for l in gdat.indxpopl:
        sampvarb[gdat.indxsampmeanpnts[l]] = icdf_logt(samp[gdat.indxsampmeanpnts[l]], gdat.minmmeanpnts[l], gdat.factmeanpnts[l])
        if gdat.fluxdisttype[l] == 'powr':
            sampvarb[gdat.indxsampfluxdistslop[l]] = icdf_atan(samp[gdat.indxsampfluxdistslop[l]], gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
        if gdat.fluxdisttype[l] == 'brok':
            sampvarb[gdat.indxsampfluxdistbrek[l]] = icdf_logt(samp[gdat.indxsampfluxdistbrek[l]], gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
            sampvarb[gdat.indxsampfluxdistsloplowr[l]] = icdf_atan(samp[gdat.indxsampfluxdistsloplowr[l]], gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
            sampvarb[gdat.indxsampfluxdistslopuppr[l]] = icdf_atan(samp[gdat.indxsampfluxdistslopuppr[l]], gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])

    for k in gdat.indxpsfipara:
        sampvarb[gdat.indxsamppsfipara[k]] = icdf_psfipara(gdat, samp[gdat.indxsamppsfipara[k]], k)
    for c in gdat.indxback:
        sampvarb[gdat.indxsampnormback[c, :]] = icdf_logt(samp[gdat.indxsampnormback[c, :]], gdat.minmnormback[c], gdat.factnormback[c])
    
    for l in gdat.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
        if gdat.fluxdisttype[l] == 'powr':
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_powr(gdat, samp[indxsampspec[l][gdat.indxenerfluxdist, :]], sampvarb[gdat.indxsampfluxdistslop[l]])
        if gdat.fluxdisttype[l] == 'brok':
            fluxunit = samp[indxsampspec[l][gdat.indxenerfluxdist[0], :]]
            fluxdistbrek = sampvarb[gdat.indxsampfluxdistbrek[l]]
            fluxdistsloplowr = sampvarb[gdat.indxsampfluxdistsloplowr[l]]
            fluxdistslopuppr = sampvarb[gdat.indxsampfluxdistslopuppr[l]]
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_brok(gdat, fluxunit, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        sampvarb[indxsampsind[l]] = icdf_eerr(samp[indxsampsind[l]], gdat.sinddistmean[l], gdat.sinddiststdv[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
        sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampspec[l][gdat.indxenerfluxdist[0], :]], sampvarb[indxsampsind[l]])
    
    return sampvarb
    

def retr_maps(gdat, indxpntsfull, sampvarb):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
     
    listspectemp = []
    for l in gdat.indxpopl:
        listspectemp.append(sampvarb[indxsampspec[l]])

    pntsflux = retr_pntsflux(gdat, sampvarb[concatenate(indxsamplgal)], sampvarb[concatenate(indxsampbgal)], \
        concatenate(listspectemp, axis=1), sampvarb[gdat.indxsamppsfipara], gdat.modlpsfntype)
    
    totlflux = retr_rofi_flux(gdat, sampvarb[gdat.indxsampnormback], pntsflux, gdat.indxcube)
    
    pntscnts = pntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    totlcnts = totlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    
    return pntsflux, pntscnts, totlflux, totlcnts


def retr_mrkrsize(gdat, flux):

    mrkrsize = (log(flux) - log(gdat.minmflux)) / (log(gdat.maxmflux) - log(gdat.minmflux)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    if False:
        print 'flux'
        print flux
        print 'gdat.maxmflux'
        print gdat.maxmflux
        print 'gdat.minmflux'
        print gdat.minmflux
        print 'mrkrsize'
        print mrkrsize

    return mrkrsize


def test_intp():

    pass 

    return 


def retr_chanpsfn(gdat):

    numbformpara = 1
    gdat.chanpsfipara = array([1., 1.]) * pi / 3600. / 180.
    psfn = retr_psfn(gdat, gdat.chanpsfipara, gdat.indxener, gdat.binsangl, psfntype='singgaus')
    return psfn


def retr_fermpsfn(gdat):
   
    if False:
        reco = 8
    else:
        reco = 7

    if reco == 8:
        path = gdat.pathdata + 'expr/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
    irfn = pf.getdata(path, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    numbfermscalpara = 3
    numbfermformpara = 5
    
    fermscal = zeros((gdat.numbevtt, numbfermscalpara))
    fermform = zeros((gdat.numbener, gdat.numbevtt, numbfermformpara))
    
    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']
    for m in gdat.indxevtt:
        if reco == 8:
            irfn = pf.getdata(path, 1 + 3 * gdat.indxevttincl[m])
            fermscal[m, :] = pf.getdata(path, 2 + 3 * gdat.indxevttincl[m])['PSFSCALE']
        else:
            if m == 1:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 0:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
            else:
                continue
            irfn = pf.getdata(path, 1)
            fermscal[m, :] = pf.getdata(path, 2)['PSFSCALE']
        for k in range(numbfermformpara):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
        
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)

    # store the fermi PSF parameters
    gdat.fermpsfipara = zeros((gdat.numbener * numbfermformpara * gdat.numbevtt))
    for m in gdat.indxevtt:
        for k in range(numbfermformpara):
            indxfermpsfiparatemp = m * numbfermformpara * gdat.numbener + gdat.indxener * numbfermformpara + k
            gdat.fermpsfipara[indxfermpsfiparatemp] = fermform[:, m, k]

    # calculate the scale factor
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # evaluate the PSF
    gdat.fermpsfn = retr_psfn(gdat, gdat.fermpsfipara, gdat.indxener, gdat.binsangl, 'doubking')


def retr_sdsspsfn(gdat):
   
    numbpsfiparaevtt = gdat.numbener * 3
    frac = array([[0.5, 0.5, 0.5]]).T
    sigc = array([[0.5, 1., 3.]]).T
    sigt = array([[0.5, 1., 3.]]).T
    gdat.sdsspsfipara = empty(numbpsfiparaevtt)
    gdat.sdsspsfn = retr_doubgaus(gdat.binsangl[None, :, None], frac[:, None, :], sigc[:, None, :], sigt[:, None, :])


def updt_samp(gdat, gdatmodi):
    
    if gdatmodi.thisindxprop == gdat.indxpropmeanpnts:
        gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
    
    # determine if a hyperparameter is to be updated
    thisbool = False
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslop or gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or \
                                                gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        thisbool = True
    
    # update the hyperparameters
    if thisbool:
        
        ## update the sample vector
        flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]]
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_powr(gdat, flux, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
            elif gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
            else:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)

        ## update the unit sample vector -- this is unique for hyperparameter updates
        gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :], -1] = fluxunit
        
        ## update the prior register
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
        gdatmodi.thislpritotl = sum(gdatmodi.thislpri)

    # proposals that change the likelihood
    if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
        gdatmodi.thisllik[gdatmodi.indxcubemodi] = gdatmodi.nextllik[gdatmodi.indxcubemodi]
        gdatmodi.thislliktotl = sum(gdatmodi.thisllik)
        gdatmodi.thismodlcnts[gdatmodi.indxcubemodi] = gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi]
        
    # PSF
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        # temp
        gdatmodi.thispsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
     
    # background normalization
    if gdatmodi.thisindxprop == gdat.indxpropnormback:
        gdatmodi.thissampvarb[gdat.indxsampnormback[gdatmodi.indxbackmodi, gdatmodi.indxenermodi]] = \
            gdatmodi.nextsampvarb[gdat.indxsampnormback[gdatmodi.indxbackmodi, gdatmodi.indxenermodi]]
        
    # proposals that change the PS flux map
    if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
        gdatmodi.thispntsflux[gdatmodi.indxcubemodi] = copy(gdatmodi.nextpntsflux[gdatmodi.indxcubemodi])

    # transdimensinal updates
    if gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
        
    ## birth
    if gdatmodi.thisindxprop == gdat.indxpropbrth:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].append(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0])
        del gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]

        ### update the components
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcomplgal]] = gdatmodi.modilgal[0]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompbgal]] = gdatmodi.modibgal[0]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspec]] = gdatmodi.modispec[:, 0]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompsind]] = gdatmodi.modisind[0]
        
    ## death
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdatmodi.killindxpnts)
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdatmodi.killindxpnts)

    ## split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:

        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].append(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0])
        del gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]
        
        ### update the components
        #### first component
        gdatmodi.thissampvarb[gdatmodi.indxsampfrst] = gdatmodi.modilgal[1]
        gdatmodi.thissampvarb[gdatmodi.indxsampfrst+1] = gdatmodi.modibgal[1]
        gdatmodi.thissampvarb[gdatmodi.indxsampfrst+2:gdatmodi.indxsampfrst+2+gdat.numbener] = gdatmodi.modispec[:, 1]
        gdatmodi.thissampvarb[gdatmodi.indxsampfrst+2+gdat.numbener] = gdatmodi.modisind[1]
        #### second component
        gdatmodi.thissampvarb[gdatmodi.indxsampseco] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+1] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+2:gdatmodi.indxsampseco+2+gdat.numbener] = gdatmodi.modispec[:, 2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+2+gdat.numbener] = gdatmodi.modisind[2]
        
    ## merge
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdatmodi.mergindxseco)
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdatmodi.mergindxseco)

        ### update the component
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcomplgal]] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompbgal]] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspec]] = gdatmodi.modispec[:, 2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompsind]] = gdatmodi.modisind[2]
        
    ## PS parameter proposals
    if gdatmodi.thisindxprop >= gdat.indxproplgal:  
        if gdatmodi.thisindxprop == gdat.indxproplgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modilgal[1]
        elif gdatmodi.thisindxprop == gdat.indxpropbgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modibgal[1]
        else:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodispec] = gdatmodi.modispec[:, 1]
            if gdatmodi.thisindxprop == gdat.indxpropsind:
                gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modisind[1]

    # update the unit sample vector
    if gdatmodi.indxsampmodi.size > 0:
        gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 0] = gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 1]


def retr_listpair(gdat, lgal, bgal):
    
    if gdat.verbtype > 1:
        print 'Finding PS pairs inside the linking length...'
    
    listpair = []
    for k in range(lgal.size):
        # temp -- linking uses the Cartesian approximation, which is accurate enough for splits and merges inside a small circle
        indxpnts = k + 1 + where(sqrt((bgal[k+1:] - bgal[k])**2 + (lgal[k+1:] - lgal[k])**2) < gdat.radispmr)[0]
        for n in range(indxpnts.size):
            listpair.append([k, indxpnts[n]])
    
    if gdat.diagmode:
        boolgood = True
        for n in range(len(listpair)):
            if sqrt((lgal[listpair[n][0]] - lgal[listpair[n][1]])**2 + (bgal[listpair[n][0]] - bgal[listpair[n][1]])**2) >= gdat.radispmr:
                boolgood = False
        if not boolgood:
            Exception('Inappropriate list of pairs')

    return listpair


def retr_fermdata(gdat):
    
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'

    fgl3 = pf.getdata(path)
    
    fgl3spectemp = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], fgl3['Flux10000_100000']))
    fgl3spectemp = fgl3spectemp[gdat.indxenerincl, :] / gdat.diffener[:, None]
    
    # sort the catalog in decreasing flux
    indxfgl3sort = argsort(fgl3spectemp[gdat.indxenerfluxdist[0], :])[::-1]

    fgl3spectemp = fgl3spectemp[:, indxfgl3sort]

    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], fgl3['Unc_Flux10000_100000']))
    
    # temp
    fgl3specstdvtemp[where(isfinite(fgl3specstdvtemp) == False)] = 0.

    fgl3specstdvtemp = fgl3specstdvtemp[gdat.indxenerincl, :, :] / gdat.diffener[:, None, None]
    fgl3specstdvtemp = fgl3specstdvtemp[:, indxfgl3sort, :]
   
    fgl3numbpntsfull = fgl3['glon'].size
    
    fgl3lgal = deg2rad(fgl3['glon'][indxfgl3sort])
    fgl3lgal = ((fgl3lgal - pi) % (2. * pi)) - pi

    fgl3bgal = deg2rad(fgl3['glat'][indxfgl3sort])
            
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'][indxfgl3sort] + fgl3['Conf_68_SemiMajor'][indxfgl3sort]) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng'][indxfgl3sort]) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3sind = fgl3['Spectral_Index'][indxfgl3sort]
    
    fgl3strg = fgl3['Source_Name'][indxfgl3sort]
    fgl3strgclss = fgl3['CLASS1'][indxfgl3sort]
    fgl3strgassc = fgl3['ASSOC1'][indxfgl3sort]
    
    fgl3spectype = fgl3['SpectrumType'][indxfgl3sort]
    fgl3scur = fgl3['beta'][indxfgl3sort]
    fgl3scut = fgl3['Cutoff'][indxfgl3sort] * 1e-3
    
    fgl3timevari = fgl3['Variability_Index'][indxfgl3sort]
    
    fgl3spec = zeros((3, gdat.numbener, fgl3numbpntsfull))
    fgl3spec[0, :, :] = fgl3spectemp
    fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdvtemp[:, :, 0]
    fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdvtemp[:, :, 1]
    
    if not isfinite(fgl3specstdvtemp).all():
        raise Exception('fgl3specstdvtemp')

    if not isfinite(fgl3spec).all():
        raise Exception('fgl3spec')

    # adjust 3FGL positions according to the ROI center
    if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
        rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='Y')
        # temp
        #rttr = hp.rotator.Rotator(rot=[0., rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='Y')
        fgl3bgal, fgl3lgal = rttr(pi / 2. - fgl3bgal, fgl3lgal)
        fgl3bgal = pi / 2. - fgl3bgal

    # select the 3FGL point sources in the ROI
    indxfgl3rofi = arange(fgl3lgal.size, dtype=int)
    for i in gdat.indxener:
        indxfgl3rofi = intersect1d(where((fgl3spec[0, i, :] > gdat.minmspec[i]) & (fgl3spec[0, i, :] < gdat.maxmspec[i]))[0], indxfgl3rofi)
    indxfgl3rofi = intersect1d(where((fabs(fgl3lgal) < gdat.maxmgangmarg) & (fabs(fgl3bgal) < gdat.maxmgangmarg))[0], indxfgl3rofi)

    # time variability
    indxfgl3vari = where(fgl3timevari[indxfgl3rofi] > 72.44)[0]
    
    # number of 3FGL PS in the ROI
    fgl3numbpnts = indxfgl3rofi.size

    # compute the 3FGL counts
    fgl3cnts = empty((gdat.numbener, fgl3numbpnts, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            # temp -- zero exposure pixels will give zero counts
            indxpixltemp = retr_indxpixl(gdat, fgl3bgal[indxfgl3rofi], fgl3lgal[indxfgl3rofi])
            fgl3cnts[i, :, m] = fgl3spec[0, i, indxfgl3rofi] * gdat.expo[i, indxpixltemp, m] * gdat.diffener[i]

    gdat.exprnumbpnts = fgl3numbpnts
    gdat.exprlgal = fgl3lgal[indxfgl3rofi]
    gdat.exprbgal = fgl3bgal[indxfgl3rofi]
    gdat.exprspec = fgl3spec[:, :, indxfgl3rofi]
    gdat.exprcnts = fgl3cnts
    gdat.exprsind = fgl3sind[indxfgl3rofi]
    gdat.exprstrg = fgl3strg[indxfgl3rofi]
    gdat.exprstrgclss = fgl3strgclss[indxfgl3rofi]
    gdat.exprstrgassc = fgl3strgassc[indxfgl3rofi]
    gdat.indxexprvari = indxfgl3vari
    
    path = gdat.pathimag + '3fgl/'
    os.system('mkdir -p %s' % path)
    path += '3fglspecaxisstdv.pdf' 
    if not os.path.isfile(path):
        print 'fgl3axisstdv'
        print fgl3axisstdv
        print 'gl3spec[0, gdat.indxenerfluxdist[0], :]'
        print fgl3spec[0, gdat.indxenerfluxdist[0], :]
        print 

        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.scatter(fgl3spec[0, gdat.indxenerfluxdist[0], :], fgl3axisstdv)
        
        #axis.set_xlim([amin(fgl3spec[0, gdat.indxenerfluxdist[0], :]), amax(fgl3spec[0, gdat.indxenerfluxdist[0], :])])
        axis.set_ylim([0., .5])
        axis.set_xlim([1e-11, 1e-7])
        #axis.set_ylim([amin(fgl3axisstdv), amax(fgl3axisstdv)])
        axis.set_xscale('log')
        axis.set_xlabel('$f$ [1/cm^2/s/GeV]')
        axis.set_ylabel('$\sigma_{r}$ [deg]')
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)


def retr_rtag(gdat):
    
    rtag = '%d' % (gdat.numbswep)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxsamp, stdv):
    
    if gdat.fracrand > 0.:
        if rand() < gdat.fracrand:
            gdatmodi.drmcsamp[indxsamp, 1] = rand()
        else:
            gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_angldistunit(gdat, lgal1, bgal1, indxpixltemp, retranglcosi=False):
    
    xaxi, yaxi, zaxi = retr_unit(lgal1, bgal1)
    anglcosi = gdat.xaxigrid[indxpixltemp] * xaxi + gdat.yaxigrid[indxpixltemp] * yaxi + gdat.zaxigrid[indxpixltemp] * zaxi
    if retranglcosi:
        return anglcosi
    else:
        angldist = arccos(anglcosi)
        return angldist


def retr_singgaus(scaldevi, sigc):
    
    psfn = 1. / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_singking(scaldevi, sigc, gamc):
    
    psfn = 1. / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc)

    return psfn


def retr_doubgaus(scaldevi, frac, sigc, sigt):
    
    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_gausking(scaldevi, frac, sigc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) + \
        (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_lgalbgal(gang, aang):
    
    lgal = gang * cos(aang)
    bgal = gang * sin(aang)

    return lgal, bgal


def retr_gang(lgal, bgal):
    
    gang = rad2deg(arccos(cos(lgal) * cos(bgal)))

    return gang


def retr_aang(lgal, bgal):

    aang = arctan2(bgal, lgal)

    return aang


def retr_enerstrg(gdat):
    
    binsenerstrg = []
    if gdat.exprtype == 'ferm':
        enerstrg = []
        for i in gdat.indxener:
            binsenerstrg.append('%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]))
            enerstrg.append('%.3g GeV' % gdat.meanener[i])
            
    if gdat.exprtype == 'sdss':
        gdat.binsenerstrg = ['i-band', 'r-band', 'g-band']
        enerstrg = ['i-band', 'r-band', 'g-band']
    
    if gdat.exprtype == 'chan':
        enerstrg = []
        for i in gdat.indxener:
            binsenerstrg.append('%.3g KeV - %.3g KeV' % (gdat.binsener[i] * 1e3, gdat.binsener[i+1] * 1e3))
            enerstrg.append('%.3g KeV' % (gdat.meanener[i] * 1e3))
    
    return enerstrg, binsenerstrg


def retr_prop(gdat, gdatmodi):
 
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal,  gdatmodi.thisindxsampspec, \
            gdatmodi.thisindxsampsind, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    if gdat.verbtype > 1:
        print 'retr_prop(): '

        print 'drmcsamp, thissampvarb'
        for k in range(gdatmodi.thissampvarb.size):
            if k == gdat.indxsampcompinit:
                print
            if k > gdat.indxsampcompinit and (k - gdat.indxsampcompinit) % gdat.numbcomp == 0:
                print
            print '%14.4g %14.4g %14.4g' % (gdatmodi.drmcsamp[k, 0], gdatmodi.drmcsamp[k, 1], gdatmodi.thissampvarb[k])
        print
            
        print 'thisindxpntsfull: ', gdatmodi.thisindxpntsfull
        print 'thisindxpntsempt: ', gdatmodi.thisindxpntsempt  
        print 'thisindxsamplgal: ', gdatmodi.thisindxsamplgal
        print 'thisindxsampbgal: ', gdatmodi.thisindxsampbgal
        print 'thisindxsampspec: '
        print gdatmodi.thisindxsampspec
        print 'thisindxsampsind: ', gdatmodi.thisindxsampsind
        print 'thisindxsampcomp: ', gdatmodi.thisindxsampcomp
        print
        
    # hyperparameter changes
    # mean number of point sources
    if gdatmodi.thisindxprop == gdat.indxpropmeanpnts:
        gdatmodi.indxsampmodi = gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvmeanpnts)
        gdatmodi.nextsampvarb[gdat.indxsampmeanpnts] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
            gdat.minmmeanpnts[gdatmodi.indxpoplmodi], gdat.factmeanpnts[gdatmodi.indxpoplmodi])
        
    # flux distribution function shape
    ## flux distribution power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslop:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistslop[gdatmodi.indxpoplmodi], gdat.factfluxdistslop[gdatmodi.indxpoplmodi])
        gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]))
        if gdat.verbtype > 1:
            print 'indxsampvarbmodi'
            print gdatmodi.indxsampvarbmodi
            print 'thissampvarb[gdat.indxsampfluxdistslop]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistslop]
            print 'nextsampvarb[gdat.indxsampfluxdistslop]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop]
    
    ## FDF break flux
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvflux)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistbrek[gdatmodi.indxpoplmodi], gdat.factfluxdistbrek[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistbrek]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek]
            print 'nextsampvarb[gdat.indxsampfluxdistbrek]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek]
    
    ## FDF lower power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistsloplowr[gdatmodi.indxpoplmodi], gdat.factfluxdistsloplowr[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistsloplowr]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr]
            print 'nextsampvarb[gdat.indxsampfluxdistsloplowr]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr]
    
    ## FDF upper power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistslopuppr[gdatmodi.indxpoplmodi], gdat.factfluxdistslopuppr[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistslopuppr]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr]
            print 'nextsampvarb[gdat.indxsampfluxdistslopuppr]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr]
   
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]))

    # PSF parameter change 
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        
        # index of the PSF parameter to change
        gdat.indxpsfiparamodi = choice(gdat.indxpsfipara)

        # the energy bin of the PS flux map to be modified
        gdatmodi.indxenermodi = array([(gdat.indxpsfiparamodi % gdat.numbpsfiparaevtt) // gdat.numbformpara])
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdat.indxsamppsfipara[gdat.indxpsfiparamodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvpsfipara)
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara] = copy(gdatmodi.thissampvarb[gdat.indxsamppsfipara])
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara[gdat.indxpsfiparamodi]] = \
            icdf_psfipara(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.indxpsfiparamodi)
            
        gdatmodi.numbmodipnts = int(sum(gdatmodi.thissampvarb[gdat.indxsampnumbpnts]))
                    
        # construct the proposed PSF
        gdatmodi.nextpsfn = retr_psfn(gdat, gdatmodi.nextsampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.modlpsfntype)
        
        temp = retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.modlpsfntype)
        gdatmodi.nextpsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        
        if gdat.verbtype > 1:
           
            print 'indxpsfiparamodi: ', gdat.indxpsfiparamodi
            print 'indxenermodi: ', gdatmodi.indxenermodi
            print 'indxsampmodi: ', gdatmodi.indxsampmodi
            print 'thispsfipara'
            print gdatmodi.thissampvarb[gdat.indxsamppsfipara]
            print 'nextpsfipara'
            print gdatmodi.nextsampvarb[gdat.indxsamppsfipara]
            print 'thissampvarb[indxsampmodi]: ', gdatmodi.thissampvarb[gdatmodi.indxsampmodi]
            print 'nextpsfipara[indxsampmodi]: ', gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
            print

    # background changes
    if gdatmodi.thisindxprop == gdat.indxpropnormback:

        ## determine the sample index to be changed
        gdatmodi.indxenermodi = choice(gdat.indxener)
        gdatmodi.indxbackmodi = choice(gdat.indxback)
        gdatmodi.indxsampmodi = gdat.indxsampnormback[gdatmodi.indxbackmodi, gdatmodi.indxenermodi]
        
        ## propose
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvback)

        ## transform back from the unit space
        gdatmodi.nextsampvarb[gdat.indxsampnormback] = copy(gdatmodi.thissampvarb[gdat.indxsampnormback])
        gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
                                                                                    gdat.minmnormback[gdatmodi.indxbackmodi], gdat.factnormback[gdatmodi.indxbackmodi])

        if gdat.verbtype > 1:
            print 'indxenermodi: ', gdatmodi.indxenermodi
            print 'indxbackmodi: ', gdatmodi.indxbackmodi
            print 'indxsampmodi: ', gdatmodi.indxsampmodi
            print 'thissampvarb[gdat.indxsampnormback]: ', gdatmodi.thissampvarb[gdat.indxsampnormback]
            print 'nextsampvarb[gdat.indxsampnormback]: ', gdatmodi.nextsampvarb[gdat.indxsampnormback]
            print

    # birth
    if gdatmodi.thisindxprop == gdat.indxpropbrth:

        # temp -- modi
        gdatmodi.numbmodipnts = 1
        #thisnumbpntsmodi = gdat.maxmnumbpnts[gdatmodi.indxpoplmodi] - int(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]])
        #gdatmodi.numbmodipnts = choice(gdat.listnumbpntsmodi[thisnumbpntsmodi], p=gdat.probnumbpntsmodi[thisnumbpntsmodi])

        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + gdatmodi.numbmodipnts
    
        # initial sample index to add the new PS
        # temp -- modi
        indxsampbrth = gdat.indxsampcompinit + gdat.maxmnumbcompcumu[gdatmodi.indxpoplmodi] + \
                                                                    array(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][:gdatmodi.numbmodipnts]) * gdat.numbcomp
        if False:
            print 'hey'
            print 'thisnumbpntsmodi'
            print thisnumbpntsmodi
            print 'gdatmodi.numbmodipnts'
            print gdatmodi.numbmodipnts
            print 'gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][:gdatmodi.numbmodipnts]'
            print gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][:gdatmodi.numbmodipnts]
            print 'indxsampbrth'
            print indxsampbrth
            print 
        # sample auxiliary variables
        numbcompmodi = gdatmodi.numbmodipnts * gdat.numbcomp
        numbcompcolrmodi = gdatmodi.numbmodipnts * gdat.numbcompcolr
        gdatmodi.auxipara = rand(numbcompcolrmodi)
        gdatmodi.indxsampmodi = empty(numbcompmodi, dtype=int)
        for k in range(gdatmodi.numbmodipnts):
            gdatmodi.drmcsamp[indxsampbrth[k]+gdat.indxcomplbhl, -1] = gdatmodi.auxipara[k*gdatmodi.numbmodipnts:k*gdatmodi.numbmodipnts+2]
            gdatmodi.drmcsamp[indxsampbrth[k]+gdat.indxcompflux, -1] = gdatmodi.auxipara[k*gdatmodi.numbmodipnts+2]
            gdatmodi.drmcsamp[indxsampbrth[k]+gdat.indxcompsind, -1] = gdatmodi.auxipara[k*gdatmodi.numbmodipnts+3]
            # sample indices to be modified
            gdatmodi.indxsampmodi[k*gdat.numbcomp:(k+1)*gdat.numbcomp] = indxsampbrth[k] + gdat.indxcomp

        # modification catalog
        gdatmodi.modilgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcomplgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        gdatmodi.modibgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompbgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        fluxunit = gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1]
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_powr(gdat, fluxunit, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
            fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_brok(gdat, array([fluxunit]), fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        gdatmodi.modisind[0] = icdf_eerr(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                    gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        gdatmodi.modispec[:, 0] = retr_spec(gdat, modiflux, gdatmodi.modisind[0]).flatten()
    
        if gdat.verbtype > 1:
            print 'numbmodipnts'
            print gdatmodi.numbmodipnts
            print 'auxipara: ', gdatmodi.auxipara
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
    # kill
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        gdatmodi.killindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        gdatmodi.indxsampmodi = array([])
            
        # modification catalog
        gdatmodi.numbmodipnts = 1
        gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modispec[:, 0] = -gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, killindxindxpnts]]

        if gdat.verbtype > 1:
            print 'killindxpnts: ', gdatmodi.killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
  
    # split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:
        
        gdatmodi.numbmodipnts = 3
        
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + \
                                                            int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.numbcomp
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp
        gdatmodi.indxsampseco = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + int(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]) * gdat.numbcomp
        indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, indxfinlfrst, dtype=int), arange(gdatmodi.indxsampseco, indxfinlseco, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, gdatmodi.spltindxindxpnts]]
        gdatmodi.fluxpare = thisspec[gdat.indxenerfluxdist[0]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        
        if gdat.verbtype > 1:
            print 'spltindxindxpnts: ', gdatmodi.spltindxindxpnts
            print 'indxsampfrst: ', gdatmodi.indxsampfrst
            print 'indxfinlfrst: ', indxfinlfrst
            print 'indxsampseco: ', gdatmodi.indxsampseco
            print 'indxfinlseco: ', indxfinlseco
            print 'thislgal: ', gdat.anglfact * thislgal
            print 'thisbgal: ', gdat.anglfact * thisbgal
            print 'thisspec: ', thisspec
            print 'thisflux: ', gdatmodi.fluxpare
            print 'thissind: ', thissind
            
        # determine the new components
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.numbcompcolr)
        gdatmodi.auxipara[0] = rand()
        gdatmodi.auxipara[1] = rand() * gdat.radispmr
        gdatmodi.auxipara[2] = rand() * 2. * pi
        gdatmodi.auxipara[3] = icdf_eerr(rand(), gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                    gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        
        if gdat.verbtype > 1:
            print 'auxipara[0]: ', gdatmodi.auxipara[0]
            print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
            print 'auxipara[2]: ', gdatmodi.auxipara[2]
            print 'auxipara[3]: ', gdatmodi.auxipara[3]
            print
            
        gdatmodi.fluxfrst = gdatmodi.auxipara[0] * gdatmodi.fluxpare
        gdatmodi.spltlgalfrst = thislgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        gdatmodi.spltbgalfrst = thisbgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        gdatmodi.spltsindfrst = thissind
        
        gdatmodi.fluxseco = (1. - gdatmodi.auxipara[0]) * gdatmodi.fluxpare
        gdatmodi.spltlgalseco = thislgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        gdatmodi.spltbgalseco = thisbgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        gdatmodi.spltsindseco = gdatmodi.auxipara[3]
        
        if gdat.verbtype > 1:
            print 'spltlgalfrst: ', gdat.anglfact * gdatmodi.spltlgalfrst
            print 'spltlgalseco: ', gdat.anglfact * gdatmodi.spltlgalseco
            print 'spltbgalfrst: ', gdat.anglfact * gdatmodi.spltbgalfrst
            print 'spltbgalseco: ', gdat.anglfact * gdatmodi.spltbgalseco
            print 'spltfluxfrst: ', gdatmodi.fluxfrst
            print 'spltfluxseco: ', gdatmodi.fluxseco
            print 'spltsindfrst: ', gdatmodi.spltsindfrst
            print 'spltsindseco: ', gdatmodi.spltsindseco

        if fabs(gdatmodi.spltlgalfrst) > gdat.maxmgangmarg or fabs(gdatmodi.spltlgalseco) > gdat.maxmgangmarg or \
                                            fabs(gdatmodi.spltbgalfrst) > gdat.maxmgangmarg or fabs(gdatmodi.spltbgalseco) > gdat.maxmgangmarg or \
                                            gdatmodi.fluxfrst < gdat.minmflux or gdatmodi.fluxseco < gdat.minmflux:
            gdatmodi.boolreje = True

        if gdat.verbtype > 1:
            print 'boolreje'
            print gdatmodi.boolreje

        if not gdatmodi.boolreje:

            lgal = concatenate((array([gdatmodi.spltlgalfrst, gdatmodi.spltlgalseco]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([gdatmodi.spltbgalfrst, gdatmodi.spltbgalseco]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]], thisbgal)))
            
            listpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.numbpair = len(listpair)

            if gdatmodi.numbpair == 0:
                print 'Number of pairs should not be zero in the reverse proposal of a split'
                raise

            ## first new component
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcomplgal, -1] = cdfn_self(gdatmodi.spltlgalfrst, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompbgal, -1] = cdfn_self(gdatmodi.spltbgalfrst, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, gdatmodi.fluxfrst, \
                                                                                    gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompsind, -1] = cdfn_eerr(gdatmodi.spltsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspecfrst = retr_spec(gdat, gdatmodi.fluxfrst, gdatmodi.spltsindfrst)

            ## second new component
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcomplgal, -1] = cdfn_self(gdatmodi.spltlgalseco, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompbgal, -1] = cdfn_self(gdatmodi.spltbgalseco, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, gdatmodi.fluxseco, \
                                                                                    gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompsind, -1] = cdfn_eerr(gdatmodi.spltsindseco, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspecseco = retr_spec(gdat, gdatmodi.fluxseco, gdatmodi.spltsindseco)

            ## component to be removed
            gdatmodi.modilgal[0] = thislgal
            gdatmodi.modibgal[0] = thisbgal
            gdatmodi.modispec[:, 0] = -thisspec.flatten()
            gdatmodi.modisind[0] = thissind
            
            ## first component to be added
            gdatmodi.modilgal[1] = gdatmodi.spltlgalfrst
            gdatmodi.modibgal[1] = gdatmodi.spltbgalfrst
            gdatmodi.modispec[:, 1] = nextspecfrst.flatten()
            gdatmodi.modisind[1] = gdatmodi.spltsindfrst

            # second component to be added
            gdatmodi.modilgal[2] = gdatmodi.spltlgalseco
            gdatmodi.modibgal[2] = gdatmodi.spltbgalseco
            gdatmodi.modispec[:, 2] = nextspecseco.flatten()
            gdatmodi.modisind[2] = gdatmodi.spltsindseco

    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        # number of point sources to be modified
        gdatmodi.numbmodipnts = 3
        
        # proposed number of point sources
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # list of point source pairs available for merge proposal
        listpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]])
        gdatmodi.numbpair = len(listpair)
        
        if gdat.verbtype > 1:
            print 'listpair'
            print listpair
           
        # check if merge will be proposed
        if gdatmodi.numbpair == 0:
            gdatmodi.boolreje = True
        else:

            # sample a pair
            indxpairtemp = choice(arange(gdatmodi.numbpair))

            # determine PS indices to be merged
            mergindxindxpntsfrst = listpair[indxpairtemp][0]
            mergindxindxpntsseco = listpair[indxpairtemp][1]
  
            ## first PS index to be merged
            gdatmodi.mergindxfrst = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]

            ## second PS index to be merged
            gdatmodi.mergindxseco = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsseco]

            # determine indices of the modified elements in the sample vector
            ## first PS
            gdatmodi.indxsampfrst = gdat.indxsampcompinit + gdat.numbcomp * gdatmodi.mergindxfrst
            indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp
            
            ## second PS
            gdatmodi.indxsampseco = gdat.indxsampcompinit + gdat.numbcomp * gdatmodi.mergindxseco
            indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp

            # indices of the sample vector elements to be modified
            gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, indxfinlfrst)

            # indices of the PS to be merged
            mergindxpnts = sort(array([gdatmodi.mergindxfrst, gdatmodi.mergindxseco], dtype=int))

            # PS parameters to be merged
            ## first PS
            gdatmodi.lgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.bgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.specfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsfrst]]
            gdatmodi.sindfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.fluxfrst = gdatmodi.specfrst[gdat.indxenerfluxdist[0]]

            ## second PS
            gdatmodi.lgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.bgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.specseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsseco]]
            gdatmodi.sindseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.fluxseco = gdatmodi.specseco[gdat.indxenerfluxdist[0]]

            # auxiliary parameters
            auxifrac = gdatmodi.fluxfrst / (gdatmodi.fluxfrst + gdatmodi.fluxseco) 
            auxiradi = sqrt((gdatmodi.lgalseco - gdatmodi.lgalfrst)**2 + (gdatmodi.bgalseco - gdatmodi.bgalfrst)**2)
            auxiangl = pi + arctan2(gdatmodi.bgalseco - gdatmodi.bgalfrst, gdatmodi.lgalseco - gdatmodi.lgalfrst)
            auxisind = gdatmodi.sindseco

            # temp
            gdatmodi.auxipara = zeros(gdat.numbcompcolr)
            gdatmodi.auxipara[0] = auxifrac
            gdatmodi.auxipara[1] = auxiradi
            gdatmodi.auxipara[2] = auxiangl
            gdatmodi.auxipara[3] = gdatmodi.sindseco
            
            # merged PS
            gdatmodi.fluxpare = gdatmodi.fluxfrst + gdatmodi.fluxseco
            if gdatmodi.fluxpare > gdat.maxmflux:
                gdatmodi.boolreje = True
            gdatmodi.lgalpare = gdatmodi.lgalfrst + (1. - auxifrac) * (gdatmodi.lgalseco - gdatmodi.lgalfrst)
            gdatmodi.bgalpare = gdatmodi.bgalfrst + (1. - auxifrac) * (gdatmodi.bgalseco - gdatmodi.bgalfrst)
            gdatmodi.sindpare = gdatmodi.sindfrst
            gdatmodi.specpare = retr_spec(gdat, gdatmodi.fluxpare, gdatmodi.sindpare)

            # determine the unit variables for the merged PS
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst, -1] = cdfn_self(gdatmodi.lgalpare, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+1, -1] = cdfn_self(gdatmodi.bgalpare, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+2, -1] = cdfn_flux_powr(gdat, gdatmodi.fluxpare, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+3, -1] = gdatmodi.drmcsamp[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsfrst], -2]

            # PSs to be added to the PS flux map
            ## first PS
            gdatmodi.modilgal[0] = gdatmodi.lgalfrst
            gdatmodi.modibgal[0] = gdatmodi.bgalfrst
            gdatmodi.modispec[:, 0] = -gdatmodi.specfrst.flatten()
            gdatmodi.modisind[0] = gdatmodi.sindfrst

            ## second PS
            gdatmodi.modilgal[1] = gdatmodi.lgalseco
            gdatmodi.modibgal[1] = gdatmodi.bgalseco
            gdatmodi.modispec[:, 1] = -gdatmodi.specseco.flatten()
            gdatmodi.modisind[1] = gdatmodi.sindseco

            ## parent PS
            gdatmodi.modilgal[2] = gdatmodi.lgalpare
            gdatmodi.modibgal[2] = gdatmodi.bgalpare
            gdatmodi.modispec[:, 2] = gdatmodi.specpare.flatten()
            gdatmodi.modisind[2] = gdatmodi.sindpare

            if gdat.verbtype > 1:
                print 'mergindxfrst: ', gdatmodi.mergindxfrst
                print 'mergindxindxpntsfrst: ', mergindxindxpntsfrst
                print 'mergindxseco: ', gdatmodi.mergindxseco
                print 'mergindxindxpntsseco: ', mergindxindxpntsseco
                print 'indxsampfrst: ', gdatmodi.indxsampfrst
                print 'indxfinlfrst: ', indxfinlfrst
                print 'indxsampseco: ', gdatmodi.indxsampseco
                print 'indxfinlseco: ', indxfinlseco
                print 'merglgalfrst: ', gdat.anglfact * gdatmodi.lgalfrst
                print 'mergbgalfrst: ', gdat.anglfact * gdatmodi.bgalfrst
                print 'mergfluxfrst: ', gdatmodi.fluxfrst
                print 'mergsindfrst: ', gdatmodi.sindfrst
                print 'merglgalseco: ', gdat.anglfact * gdatmodi.lgalseco
                print 'mergbgalseco: ', gdat.anglfact * gdatmodi.bgalseco
                print 'mergfluxseco: ', gdatmodi.fluxseco
                print 'mergsindseco: ', gdatmodi.sindseco
                print 'merglgalpare: ', gdat.anglfact * gdatmodi.lgalpare
                print 'mergbgalpare: ', gdat.anglfact * gdatmodi.bgalpare
                print 'mergspecpare: ', gdatmodi.specpare
                print 'mergfluxpare: ', gdatmodi.fluxpare
                print 'mergsindpare: ', gdatmodi.sindpare
                print 'auxipara[0]: ', gdatmodi.auxipara[0]
                print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
                print 'auxipara[2]: ', gdatmodi.auxipara[2]
                print 'auxipara[3]: ', gdatmodi.auxipara[3]
                print

            if auxiradi > gdat.radispmr:
                print 'Auxiliary radius during a merge cannot be larger than the linking length of %.3g %s.' % (gdat.anglfact * gdat.radispmr, gdat.strganglunit)
                raise

    # PS parameter change
    if gdatmodi.thisindxprop >= gdat.indxproplgal:     
        
        gdatmodi.indxenermodi = gdat.indxener
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdatmodi.thisindxprop == gdat.indxproplgal:
                gdatmodi.indxcompmodi = gdat.indxcomplgal
            else:
                gdatmodi.indxcompmodi = gdat.indxcompbgal
        else:
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                gdatmodi.indxcompmodi = gdat.indxcompflux
            elif gdatmodi.thisindxprop == gdat.indxpropsind:
                gdatmodi.indxcompmodi = gdat.indxcompsind
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        gdatmodi.indxsampmodiinit = gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:gdatmodi.indxpoplmodi]) + modiindxpnts * gdat.numbcomp
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdatmodi.indxsampmodiinit + gdatmodi.indxcompmodi
        gdatmodi.indxsampmodispec = gdatmodi.indxsampmodiinit + 2 + gdat.indxener
        
        thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], modiindxindxpnts]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
        thisspec = retr_spec(gdat, thisflux, thissind)
            
        # propose
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdat.stdvlbhlvari:
                retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl * gdat.minmflux / thisflux)
            else:
                retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl) 
        if gdatmodi.thisindxprop == gdat.indxpropflux:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvflux)
        if gdatmodi.thisindxprop == gdat.indxpropsind:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvsind)

        gdatmodi.numbmodipnts = 2
        gdatmodi.modispec[:, 0] = -thisspec.flatten()
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdatmodi.indxcompmodi == 0:
                gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modilgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
                gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            else:
                gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modibgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.modispec[:, 1] = thisspec.flatten()
        else:
            gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_powr(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fluxdistslop)
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_brok(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
                gdatmodi.modisind[1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            else:
                modiflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, modiindxindxpnts]]
                gdatmodi.modisind[1] = icdf_eerr(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                        gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            gdatmodi.modispec[:, 1] = retr_spec(gdat, modiflux, gdatmodi.modisind[1]).flatten()

        if gdat.verbtype > 1:
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print 'indxcompmodi: ', gdatmodi.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts

    if gdat.verbtype > 1:
        print 'indxpoplmodi'
        print gdatmodi.indxpoplmodi
        if not (gdatmodi.boolreje or gdatmodi.thisindxprop == gdat.indxpropdeth):
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
            print 'drmcsamp[gdatmodi.indxsampmodi, :]'
            print gdatmodi.drmcsamp[gdatmodi.indxsampmodi, :]
        
    # energy bin in which to evaluate the log-likelihood
    if gdat.indxpropbrth <= gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdatmodi.indxenermodi = gdat.indxener

    if gdat.verbtype > 1:
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            print 'indxenermodi: ', gdatmodi.indxenermodi

    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg) and not gdatmodi.boolreje:
        
        ## Jacobian
        jcbnfacttemp = gdatmodi.fluxpare * fabs(gdatmodi.auxipara[1] * (sin(gdatmodi.auxipara[2]) * cos(gdatmodi.auxipara[2]) + cos(gdatmodi.auxipara[2])**2))
        if gdatmodi.thisindxprop == gdat.indxpropsplt:
            gdatmodi.jcbnfact = jcbnfacttemp
        else:
            gdatmodi.jcbnfact = 1. / jcbnfacttemp
        
        ## combinatorial factor
        combfacttemp = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.numbpair
        if gdatmodi.thisindxprop == gdat.indxpropsplt:
            gdatmodi.combfact = combfacttemp
        else:
            gdatmodi.combfact = 1. / combfacttemp
        
        gdatmodi.laccfact = log(gdatmodi.jcbnfact * gdatmodi.combfact)

        if gdat.verbtype > 1:
            print 'jcbnfact'
            print gdatmodi.jcbnfact
            print 'combfact'
            print gdatmodi.combfact
            print 'laccfact'
            print gdatmodi.laccfact
            print

    else:
        gdatmodi.jcbnfact = 0.
        gdatmodi.combfact = 0.
        gdatmodi.laccfact = 0.  
   
    # temp
    # define the index to be changed in the sample vector if a hyperparameter is being updated
    if not gdatmodi.boolreje and (gdatmodi.thisindxprop != gdat.indxpropfluxdistslop and gdatmodi.thisindxprop != gdat.indxpropfluxdistbrek and \
                                            gdatmodi.thisindxprop != gdat.indxpropfluxdistsloplowr and gdatmodi.thisindxprop != gdat.indxpropfluxdistslopuppr):
        gdatmodi.indxsampvarbmodi = gdatmodi.indxsampmodi

        
def retr_psfn(gdat, psfipara, indxenertemp, thisangl, psfntype):

    if psfntype == 'singgaus':
        numbformpara = 1
    elif psfntype == 'singking':
        numbformpara = 2
    elif psfntype == 'doubgaus':
        numbformpara = 3
    elif psfntype == 'gausking':
        numbformpara = 4
    elif psfntype == 'doubking':
        numbformpara = 5
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        scalangl = thisangl[None, :, None]

    indxpsfiparatemp = numbformpara * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    
    if psfntype == 'singgaus':
        sigc = psfipara[indxpsfiparatemp]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    elif psfntype == 'singking':
        sigc = psfipara[indxpsfiparatemp]
        gamc = psfipara[indxpsfiparatemp+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(scalangl, sigc, gamc)
        
    elif psfntype == 'doubgaus':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(scalangl, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        gamt = psfipara[indxpsfiparatemp+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(scalangl, frac, sigc, sigt, gamt)
        
    elif psfntype == 'doubking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        gamc = psfipara[indxpsfiparatemp+2]
        sigt = psfipara[indxpsfiparatemp+3]
        gamt = psfipara[indxpsfiparatemp+4]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
    
    # normalize the PSF
    if gdat.exprtype == 'ferm':
        psfn /= 2. * pi * trapz(psfn * sin(thisangl[None, :, None]), thisangl, axis=1)[:, None, :]
    
    # temp
    if True and (gdat.strgcnfg == 'pcat_ferm_expr_ngal' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp1' or \
                                                            gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp2' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp3'):
        print 'CORRECTING THE PSF.'
        tempcorr = array([1., 0.8, 0.8])
        psfn *= tempcorr[:, None, None]

    return psfn


def retr_psfimodl(gdat):
    
    # PSF parameters
    if gdat.modlpsfntype == 'singgaus':
        gdat.numbformpara = 1
    elif gdat.modlpsfntype == 'singking':
        gdat.numbformpara = 2 
    elif gdat.modlpsfntype == 'doubgaus':
        gdat.numbformpara = 3
    elif gdat.modlpsfntype == 'gausking':
        gdat.numbformpara = 4
    elif gdat.modlpsfntype == 'doubking':
        gdat.numbformpara = 5
       
    gdat.indxformpara = arange(gdat.numbformpara) 
    gdat.numbpsfiparaevtt = gdat.numbener * gdat.numbformpara
    gdat.numbpsfipara = gdat.numbpsfiparaevtt * gdat.numbevtt
    gdat.indxpsfipara = arange(gdat.numbpsfipara)

    minmformpara = zeros(gdat.numbformpara)
    maxmformpara = zeros(gdat.numbformpara)
    factformpara = zeros(gdat.numbformpara)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    if gdat.exprtype == 'ferm':
        minmanglpsfn = 0.1 # deg2rad(0.0001)
        maxmanglpsfn = 5. # deg2rad(5.)
        #minmanglpsfn = 0.01
        #maxmanglpsfn = 3.
        # temp
        minmgamm = 2.
        maxmgamm = 20.
    if gdat.exprtype == 'chan' or gdat.exprtype == 'sdss':
        minmanglpsfn = 0.5 / gdat.anglfact
        maxmanglpsfn = 2. / gdat.anglfact 
    minmpsfnfrac = 0.
    maxmpsfnfrac = 1.
    
    if gdat.modlpsfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif gdat.modlpsfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif gdat.modlpsfntype == 'doubgaus':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmanglpsfn
        maxmformpara[2] = maxmanglpsfn
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'logt'
        strgformpara = ['$f_c', r'$\sigma_c', r'$\sigma_t']
    elif gdat.modlpsfntype == 'gausking':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmanglpsfn
        maxmformpara[2] = maxmanglpsfn
        minmformpara[3] = minmgamm
        maxmformpara[3] = maxmgamm
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'logt'
        scalformpara[3] = 'atan'
        strgformpara = ['$f_g', r'$\sigma_g', r'$\sigma_k', r'$\gamma']
    elif gdat.modlpsfntype == 'doubking':
        minmformpara[0] = minmpsfnfrac
        maxmformpara[0] = maxmpsfnfrac
        minmformpara[1] = minmanglpsfn
        maxmformpara[1] = maxmanglpsfn
        minmformpara[2] = minmgamm
        maxmformpara[2] = maxmgamm
        minmformpara[3] = minmanglpsfn
        maxmformpara[3] = maxmanglpsfn
        minmformpara[4] = minmgamm
        maxmformpara[4] = maxmgamm
        scalformpara[0] = 'self'
        scalformpara[1] = 'logt'
        scalformpara[2] = 'atan'
        scalformpara[3] = 'logt'
        scalformpara[4] = 'atan'
        strgformpara = ['$f_c', r'$\sigma_c', r'$\gamma_c', r'$\sigma_t', r'$\gamma_t']

    for k in gdat.indxformpara:
        if scalformpara[k] == 'self':
            factformpara[k] = maxmformpara[k] - minmformpara[k]
        if scalformpara[k] == 'logt':
            factformpara[k] = log(maxmformpara[k] / minmformpara[k])
        if scalformpara[k] == 'atan':
            factformpara[k] = arctan(maxmformpara[k]) - arctan(minmformpara[k])
            
    gdat.minmpsfipara = tile(tile(minmformpara, gdat.numbener), gdat.numbevtt)
    gdat.maxmpsfipara = tile(tile(maxmformpara, gdat.numbener), gdat.numbevtt)
    gdat.scalpsfipara = tile(tile(scalformpara, gdat.numbener), gdat.numbevtt)
    gdat.factpsfipara = tile(tile(factformpara, gdat.numbener), gdat.numbevtt)
   
    # PSF parameter strings
    gdat.strgpsfipara = [strgformpara[k] + '^{%d%d}$' % (gdat.indxenerincl[i], gdat.indxevttincl[m]) \
        for m in gdat.indxevtt for i in gdat.indxener for k in gdat.indxformpara]
    gdat.indxpsfiparainit = (gdat.indxevtt[:, None] * gdat.numbpsfiparaevtt + gdat.indxener[None, :] * gdat.numbformpara).flatten()
    for k in arange(gdat.indxpsfiparainit.size):
        if gdat.modlpsfntype == 'singgaus' or gdat.modlpsfntype == 'singking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]] += ' ' + gdat.strganglunit
        elif gdat.modlpsfntype == 'doubgaus' or gdat.modlpsfntype == 'gausking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+2] += ' ' + gdat.strganglunit
        elif gdat.modlpsfntype == 'doubking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+3] += ' ' + gdat.strganglunit
    gdat.indxpsfipara = arange(gdat.numbpsfipara)


def retr_unit(lgal, bgal):

    xaxi = cos(bgal) * cos(lgal)
    yaxi = -cos(bgal) * sin(lgal)
    zaxi = sin(bgal)

    return xaxi, yaxi, zaxi


def retr_propmodl(gdat):

    gdat.strgprop = []
    cntr = tdpy.util.cntr()
    
    ## mean number of PS
    gdat.indxpropmeanpnts = cntr.incr()
    gdat.strgprop.append('meanpnts')

    ## flux distribution
    gdat.indxpropfluxdistslop = cntr.incr()
    gdat.strgprop.append('fluxdistslop')
    gdat.indxpropfluxdistbrek = cntr.incr()
    gdat.strgprop.append('fluxdistbrek')
    gdat.indxpropfluxdistsloplowr = cntr.incr()
    gdat.strgprop.append('fluxdistsloplowr')
    gdat.indxpropfluxdistslopuppr = cntr.incr()
    gdat.strgprop.append('fluxdistslopuppr')

    # PSF parameters
    gdat.indxproppsfipara = cntr.incr()
    gdat.strgprop.append('psfipara')
    
    # background normalization
    gdat.indxpropnormback = cntr.incr()
    gdat.strgprop.append('normback')
    
    # birth
    gdat.indxpropbrth = cntr.incr()
    gdat.strgprop.append('brth')
    
    # death
    gdat.indxpropdeth = cntr.incr()
    gdat.strgprop.append('deth')
    
    # split
    gdat.strgprop.append('splt')
    gdat.indxpropsplt = cntr.incr()
    
    # merge
    gdat.strgprop.append('merg')
    gdat.indxpropmerg = cntr.incr()
    
    # lgal
    gdat.strgprop.append('lgal')
    gdat.indxproplgal = cntr.incr()
    
    # bgal
    gdat.strgprop.append('bgal')
    gdat.indxpropbgal = cntr.incr()
    
    # spec
    gdat.strgprop.append('flux')
    gdat.indxpropflux = cntr.incr()
    
    # sind
    gdat.strgprop.append('sind')
    gdat.indxpropsind = cntr.incr()

    gdat.numbprop = len(gdat.strgprop)
    gdat.indxprop = arange(gdat.numbprop)

    if gdat.probprop == None:
            
        if gdat.boolpropfluxdistbrek:
            probfluxdistbrek = array([1.])
        else:
            probfluxdistbrek = array([0.])
            
        if gdat.boolpropfluxdist:
            probmeanpnts = array([1.])
            # temp
            if gdat.fluxdisttype[0] == 'powr':
                probfluxdistslop = array([1.])
                probfluxdistsloplowr = array([0.])
                probfluxdistslopuppr = array([0.])
            if gdat.fluxdisttype[0] == 'brok':
                probfluxdistslop = array([0.])
                probfluxdistsloplowr = array([1.])
                probfluxdistslopuppr = array([1.])
        else:
            probmeanpnts = array([0.])
            probfluxdistslop = array([0.])
            probfluxdistsloplowr = array([0.])
            probfluxdistslopuppr = array([0.])

        if gdat.boolproppsfn:
            gdat.probpsfipara = gdat.numbpsfipara
        else:
            gdat.probpsfipara = 0.
        probnormback = array([1.])
       
        probbrth = 0.2 * gdat.maxmnumbpntstotl
        probdeth = 0.2 * gdat.maxmnumbpntstotl
        probsplt = 0.  * gdat.maxmnumbpntstotl
        probmerg = 0.  * gdat.maxmnumbpntstotl
        
        problgal = gdat.maxmnumbpntstotl
        probbgal = gdat.maxmnumbpntstotl
        probspec = gdat.maxmnumbpntstotl
        if gdat.boolpropsind:
            probsind = gdat.maxmnumbpntstotl
        else:
            probsind = 0.
           
        gdat.probprop = empty(gdat.numbprop)
        gdat.probprop[0] = probmeanpnts
        gdat.probprop[1] = probfluxdistslop
        gdat.probprop[2] = probfluxdistbrek
        gdat.probprop[3] = probfluxdistsloplowr
        gdat.probprop[4] = probfluxdistslopuppr
        gdat.probprop[5] = gdat.probpsfipara
        gdat.probprop[6] = probnormback
        gdat.probprop[7] = probbrth
        gdat.probprop[8] = probdeth
        gdat.probprop[9] = probsplt
        gdat.probprop[10] = probmerg
        gdat.probprop[11] = problgal
        gdat.probprop[12] = probbgal
        gdat.probprop[13] = probspec
        gdat.probprop[14] = probsind
        
        probproptotl = sum(gdat.probprop)
        gdat.probprop /= probproptotl
        gdat.probpsfipara /= probproptotl
    else:
        gdat.probpsfipara = gdat.probprop[5]
       

def retr_randunitpsfipara(gdat):

    while True:
        randunitpsfipara = rand(gdat.numbpsfipara)
        if gdat.modlpsfntype == 'singgaus' or gdat.modlpsfntype == 'singking':
            break
        else:
            indxpar0 = 1
            if gdat.modlpsfntype == 'doubgaus' or gdat.modlpsfntype == 'gausking':
                indxpar1 = 2
            if gdat.modlpsfntype == 'doubking':
                indxpar1 = 3
            thisbool = True
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    indx = m * gdat.numbpsfiparaevtt + i * gdat.numbformpara
                    thisbool = thisbool and randunitpsfipara[indx+indxpar1] > randunitpsfipara[indx+indxpar0]
            if thisbool:
                break

    return randunitpsfipara


def setpinit(gdat):

    # paths
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    
    # process index
    gdat.indxproc = arange(gdat.numbproc)

    # axes
    ## energy
    gdat.numbener = gdat.indxenerincl.size
    gdat.numbenerfull = gdat.binsenerfull.size - 1
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.binsener = gdat.binsenerfull[gdat.indxenerinclbins]
    gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]
    gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
    gdat.indxener = arange(gdat.numbener, dtype=int)
    gdat.indxenerfluxdist = array([gdat.numbener / 2])
    gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
    factener = (gdat.meanener[gdat.indxenerfluxdist] / gdat.meanener)**2

    ## PSF class
    gdat.numbevttfull = gdat.indxevttfull.size
    gdat.numbevtt = gdat.indxevttincl.size
    gdat.indxevtt = arange(gdat.numbevtt)
    
    # angular deviation
    # temp -- check that gdat.numbangl does not degrade the performance
    gdat.numbangl = 1000
    gdat.binsangl = linspace(0., gdat.maxmangl, gdat.numbangl) # [rad]
    gdat.binsanglcosi = sort(cos(gdat.binsangl))
    
    # half size of the image where the sample catalog is compared against the reference
    gdat.maxmgangcomp = gdat.maxmgang * gdat.margfactcomp
    # half size of the spatial prior
    gdat.maxmgangmarg = gdat.maxmgang * gdat.margfactmodl

    # population index vector
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)
    if gdat.datatype == 'mock':
        gdat.mockindxpopl = arange(gdat.mocknumbpopl, dtype=int)

    if gdat.datatype == 'mock':
        gdat.mocksindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.mocksinddistmean) / gdat.mocksinddiststdv / sqrt(2.)) + 1.)
        gdat.mocksindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.mocksinddistmean) / gdat.mocksinddiststdv / sqrt(2.)) + 1.)
        gdat.mocksindcdfnnormdiff = gdat.mocksindcdfnnormmaxm - gdat.mocksindcdfnnormminm
    
    # construct the PSF
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
    if gdat.exprtype == 'chan':
        retr_chanpsfn(gdat)
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)

    if gdat.datatype == 'mock':

        # mock PSF parameters
        if gdat.mockpsfntype == 'singgaus':
            gdat.numbmockformpara = 1
        elif gdat.mockpsfntype == 'singking':
            gdat.numbmockformpara = 2 
        elif gdat.mockpsfntype == 'doubgaus':
            gdat.numbmockformpara = 3
        elif gdat.mockpsfntype == 'gausking':
            gdat.numbmockformpara = 4
        elif gdat.mockpsfntype == 'doubking':
            gdat.numbmockformpara = 5

        # temp
        #gdat.numbmockpsfiparaevtt = gdat.numbener * gdat.numbmockformpara
        #gdat.numbmockpsfipara = gdat.numbmockpsfiparaevtt * gdat.numbevtt
        #gdat.indxmockpsfipara = arange(gdat.numbmockpsfipara)   

        if gdat.exprtype == 'ferm':
            gdat.mockpsfipara = gdat.fermpsfipara
        if gdat.exprtype == 'sdss':
            gdat.mockpsfipara = gdat.sdsspsfipara
        if gdat.exprtype == 'chan':
            gdat.mockpsfipara = gdat.chanpsfipara
  
    # pixelization
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2

    if gdat.exprtype == 'ferm':
        gdat.nameexpr = 'Fermi-LAT'
    if gdat.exprtype == 'sdss':
        gdat.nameexpr = 'SDSS'
    if gdat.exprtype == 'chan':
        gdat.nameexpr = 'Chandra'
    
    gdat.numbchrototl = 4
    gdat.numbchrollik = 7

    # temp
    gdat.boolintpanglcosi = False

    # number of bins
    gdat.numbbins = 10

    gdat.minmnumbpnts = 1
    
    # the function to measure time
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

    gdat.indxback = arange(gdat.numbback)
    
    gdat.maxmnumbpntstotl = sum(gdat.maxmnumbpnts)

    # spatial priors
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
   
    # axes
    ## longitude
    gdat.numblgalpntsprob = 400
    gdat.numbbgalpntsprob = 400
    gdat.binslgalpntsprob = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbsidepntsprob + 1)
    gdat.binsbgalpntsprob = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbsidepntsprob + 1)
    gdat.indxlgalpntsprob = arange(gdat.numblgalpntsprob)
    gdat.indxbgalpntsprob = arange(gdat.numbbgalpntsprob)

    gdat.binslgal, gdat.meanlgal, gdat.difflgal, gdat.numblgal, gdat.indxlgal = tdpy.util.retr_axis(gdat.minmlgal, gdat.maxmlgal, 10)
    gdat.binslgal, gdat.meanlgal, gdat.difflgal, gdat.numblgal, gdat.indxlgal = tdpy.util.retr_axis(gdat.minmlgal, gdat.maxmlgal, 10)
    gdat.binsbgal, gdat.meanbgal, gdat.diffbgal, gdat.numbbgal, gdat.indxbgal = tdpy.util.retr_axis(gdat.minmbgal, gdat.maxmbgal, 10)

    ## radial
    gdat.numbgang = 10
    gdat.binsgang = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbgang + 1)

    ## azimuthal
    gdat.numbaang = 10
    gdat.binsaang = linspace(0., 2. * pi, gdat.numbaang + 1)

    ## flux
    gdat.numbflux = 20
    gdat.indxflux = arange(gdat.numbflux)
    gdat.binsflux = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbflux + 1)
    gdat.meanflux = sqrt(gdat.binsflux[1:] * gdat.binsflux[:-1])
    gdat.diffflux = gdat.binsflux[1:] - gdat.binsflux[:-1]
    ### pivot flux bin
    gdat.indxfluxpivt = gdat.numbflux / 2
    gdat.pivtflux = gdat.meanflux[gdat.indxfluxpivt]

    ## color
    gdat.numbsind = 20
    gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbsind + 1)
    gdat.meansind = (gdat.binssind[1:] + gdat.binssind[:-1]) / 2.
    gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    gdat.minmspec = gdat.minmflux * factener
    gdat.maxmspec = gdat.maxmflux * factener
    gdat.binsspec = gdat.binsflux[None, :] * factener[:, None]
    gdat.meanspec = empty((gdat.numbener, gdat.numbflux))
    for i in gdat.indxener:
        gdat.meanspec[i, :] = sqrt(gdat.binsspec[i, 1:] * gdat.binsspec[i, :-1])

    # input data
    if gdat.datatype == 'inpt':
        
        path = gdat.pathdata + 'inpt/' + gdat.strgexpr
        gdat.exprflux = pf.getdata(path)
        
        if gdat.pixltype == 'heal':
            if gdat.exprflux.ndim != 3:
                Exception('exprflux should be a 3D numpy array if pixelization is HealPix.')
        else:
            if gdat.exprflux.ndim != 4:
                Exception('exprflux should be a 4D numpy array if pixelization is Cartesian.')
        
        if gdat.pixltype == 'cart':
            gdat.numbsidecart = gdat.exprflux.shape[1]
            gdat.exprflux = gdat.exprflux.reshape((gdat.exprflux.shape[0], gdat.numbsidecart**2, gdat.exprflux.shape[3]))

        gdat.numbenerfull = gdat.exprflux.shape[0]
        gdat.numbpixlfull = gdat.exprflux.shape[1]
        gdat.numbevttfull = gdat.exprflux.shape[2]
        gdat.indxenerfull = arange(gdat.numbenerfull)
        gdat.indxevttfull = arange(gdat.numbevttfull)
        
        if gdat.pixltype == 'heal':
            gdat.numbsideheal = int(sqrt(gdat.numbpixlfull / 12))
    
    if gdat.pixltype == 'cart':
        gdat.binslgalcart = linspace(gdat.minmlgal, gdat.maxmlgal, gdat.numbsidecart + 1)
        gdat.binsbgalcart = linspace(gdat.minmbgal, gdat.maxmbgal, gdat.numbsidecart + 1)
        gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
        gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
        gdat.apix = (2. * gdat.maxmgang / gdat.numbsidecart)**2

    if gdat.pixltype == 'heal':
        
        lgalheal, bgalheal, gdat.numbpixlfull, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)
        lgalheal = deg2rad(lgalheal)
        bgalheal = deg2rad(bgalheal)
   
        gdat.indxpixlrofi = where((fabs(lgalheal) < gdat.maxmgang) & (fabs(bgalheal) < gdat.maxmgang))[0]
        
        gdat.indxpixlrofimargextd = where((fabs(lgalheal) < 1.2 * gdat.maxmgangmarg) & (fabs(bgalheal) < 1.2 * gdat.maxmgangmarg))[0]
        gdat.indxpixlrofimarg = where((fabs(lgalheal) < gdat.maxmgangmarg) & (fabs(bgalheal) < gdat.maxmgangmarg))[0]

        gdat.lgalgrid = lgalheal
        gdat.bgalgrid = bgalheal

    else:
        gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
        indxsidecart = arange(gdat.numbsidecart)
        temp = meshgrid(indxsidecart, indxsidecart, indexing='ij')
        gdat.bgalgrid = gdat.bgalcart[temp[1].flatten()]
        gdat.lgalgrid = gdat.lgalcart[temp[0].flatten()]
    
    gdat.indxpixlfull = arange(gdat.numbpixlfull)

    # plotting
    ## the normalized offset for text annotation of point sources in the frames
    gdat.offstext = gdat.maxmgang * 0.05
    ## figure size
    gdat.plotsize = 7
    ## text
    if gdat.datatype == 'mock':
        gdat.truelabl = 'Mock'
    if gdat.datatype == 'inpt':
        if gdat.exprtype == 'ferm':
            gdat.truelabl = '3FGL'
        if gdat.exprtype == 'sdss':
            gdat.truelabl = 'Hubble'
        if gdat.exprtype == 'chan':
            gdat.truelabl = 'Chandra'

    if gdat.exprtype == 'ferm':
        gdat.strganglunit = '$^o$'
        gdat.strganglunittext = 'degree'
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.strganglunit = '$^{\prime\prime}$'
        gdat.strganglunittext = 'arcsec'
            
    if gdat.exprtype == 'ferm':
        gdat.enerfact = 1.
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.enerfact = 1e3
    
    if gdat.strganglunit == '$^o$':
        gdat.anglfact = 180. / pi
    if gdat.strganglunit == '$^{\prime\prime}$':
        gdat.anglfact = 3600 * 180. / pi

    gdat.binsanglplot = gdat.anglfact * gdat.binsangl

    if gdat.lgalcntr == 0. and gdat.bgalcntr == 0.:
        gdat.longlabl = '$l$'
        gdat.latilabl = '$b$'
    else:
        gdat.longlabl = r'$\nu$'
        gdat.latilabl = r'$\mu$'
    
    gdat.longlabl += ' [%s]' % gdat.strganglunit
    gdat.latilabl += ' [%s]' % gdat.strganglunit
    
    gdat.truelablvari = gdat.truelabl + ' variable'
    gdat.truelablmiss = gdat.truelabl + ' miss'
    gdat.truelablbias = gdat.truelabl + ' off'
    gdat.truelablhits = gdat.truelabl + ' hit'
    gdat.truelablmult = gdat.truelabl + ' mult'

    # minimum angular distance from the center of the ROI
    gdat.minmgang = 1e-3

    # energy bin string
    gdat.enerstrg, gdat.binsenerstrg = retr_enerstrg(gdat)
    
    # PSF class string
    gdat.evttstrg = []
    for m in gdat.indxevtt:
        gdat.evttstrg.append('PSF%d' % gdat.indxevttincl[m])
        
    # number of components
    gdat.numbcomp = 3 + gdat.numbener
    gdat.indxcomp = arange(gdat.numbcomp)
    gdat.numbcompcolr = 4
    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
    
    # component indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcomplbhl = arange(2)
    gdat.indxcompspec = 2 + gdat.indxener
    gdat.indxcompflux = 2 + gdat.indxenerfluxdist
    gdat.indxcompsind = 2 + gdat.numbener
    
    # convenience factors for CDF and ICDF transforms
    ## mean number of PS
    gdat.factmeanpnts = log(gdat.maxmmeanpnts / gdat.minmmeanpnts)
    
    ## background normalization
    gdat.factnormback = log(gdat.maxmnormback / gdat.minmnormback)
    
    ## PS parameters
    gdat.factgang = log(gdat.maxmgangmarg / gdat.minmgang)

    gdat.factfluxdistslop = arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)
    gdat.factfluxdistbrek = log(gdat.maxmfluxdistbrek / gdat.minmfluxdistbrek)
    gdat.factfluxdistsloplowr = arctan(gdat.maxmfluxdistsloplowr) - arctan(gdat.minmfluxdistsloplowr)
    gdat.factfluxdistslopuppr = arctan(gdat.maxmfluxdistslopuppr) - arctan(gdat.minmfluxdistslopuppr)
    gdat.factsind = arctan(gdat.maxmsind) - arctan(gdat.minmsind)
    
    # temp -- these must be updated when updating hyperparameters on the color distribution
    gdat.sindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.sinddistmean) / gdat.sinddiststdv / sqrt(2.)) + 1.)
    gdat.sindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.sinddistmean) / gdat.sinddiststdv / sqrt(2.)) + 1.)
    gdat.sindcdfnnormdiff = gdat.sindcdfnnormmaxm - gdat.sindcdfnnormminm

    # exposure
    if gdat.strgexpo == 'unit':
        if gdat.datatype == 'mock':
            if gdat.pixltype == 'heal':
                gdat.expo= ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            if gdat.pixltype == 'cart':
                gdat.expo = ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
        if gdat.datatype == 'inpt':
            gdat.expo = ones_like(gdat.exprflux)
    else:
        path = gdat.pathdata + 'inpt/' + gdat.strgexpo
        gdat.expo = pf.getdata(path)
        if amin(gdat.expo) == amax(gdat.expo):
            print 'Bad input exposure map.'
            return
        if gdat.pixltype == 'cart':
            gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))

    # backgrounds
    gdat.backflux = []
    for c in gdat.indxback:
        if gdat.strgback[c] == 'unit':
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    backfluxtemp = ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.pixltype == 'cart':
                    backfluxtemp = ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
            if gdat.datatype == 'inpt':
                backfluxtemp = ones_like(gdat.exprflux)
        else:
            path = gdat.pathdata + 'inpt/' + gdat.strgback[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'cart':
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        gdat.backflux.append(backfluxtemp)
    
    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    if gdat.datatype == 'inpt':
        gdat.exprflux = gdat.exprflux[gdat.indxcubeincl]
    gdat.expo = gdat.expo[gdat.indxcubeincl]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcubeincl]
    
    # exclude voxels with vanishing exposure
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            gdat.indxpixlrofi = intersect1d(gdat.indxpixlrofi, where(gdat.expo[i, :, m] > 0.)[0])
    gdat.indxcuberofi = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.indxcube = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
       
    gdat.lgalgrid = gdat.lgalgrid[gdat.indxpixlrofi]
    gdat.bgalgrid = gdat.bgalgrid[gdat.indxpixlrofi]
    
    # store pixels as unit vectors
    gdat.xaxigrid, gdat.yaxigrid, gdat.zaxigrid = retr_unit(gdat.lgalgrid, gdat.bgalgrid)
   
    # construct a lookup table for converting HealPix pixels to ROI pixels
    if gdat.pixltype == 'heal':
        path = gdat.pathdata + 'pixlcnvt/'
        os.system('mkdir -p %s' % path)
        path += 'pixlcnvt_%09g.p' % gdat.maxmgang

        if os.path.isfile(path):
            if gdat.verbtype > 0:
                print 'Reading %s...' % path
            fobj = open(path, 'rb')
            gdat.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            gdat.pixlcnvt = zeros(gdat.numbpixlfull, dtype=int) - 1

            numbpixlmarg = gdat.indxpixlrofimargextd.size
            for k in range(numbpixlmarg):
                dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimargextd[k]], bgalheal[gdat.indxpixlrofimargextd[k]], gdat.indxpixl)
                gdat.pixlcnvt[gdat.indxpixlrofimargextd[k]] = argmin(dist)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
   
    if gdat.datatype == 'inpt':
        # temp
        #gdat.datacnts = gdat.datacnts[gdat.indxcuberofi]
        gdat.exprflux = gdat.exprflux[gdat.indxcuberofi]
    
    gdat.expofull = copy(gdat.expo)
    gdat.expo = gdat.expo[gdat.indxcuberofi]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcuberofi]

    # construct the PSF model
    retr_psfimodl(gdat)


def setpfinl(gdat):

    # set sample vector indices
    cntr = tdpy.util.cntr()
    gdat.indxsampnumbpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampmeanpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistslop = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistbrek = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistsloplowr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistslopuppr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsamppsfipara = arange(gdat.numbpsfipara) + cntr.incr(gdat.numbpsfipara)
    gdat.indxsampnormback = arange(gdat.numbback * gdat.numbener).reshape((gdat.numbback, gdat.numbener)) + cntr.incr(gdat.numbback * gdat.numbener)
    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.maxmnumbcompcumu = empty(gdat.numbpopl, dtype=int)
    for l in gdat.indxpopl:
        gdat.maxmnumbcompcumu[l] = sum(gdat.maxmnumbcomp[:l])
    gdat.indxsampcompinit = amax(gdat.indxsampnormback) + 1
    gdat.numbpara = int(gdat.indxsampcompinit + sum(gdat.maxmnumbcomp))
    gdat.indxpara = arange(gdat.numbpara)

    # number of burned sweeps
    if gdat.numbburn == None:
        if gdat.datatype == 'mock' and not gdat.randinit:
            minmnumbburn = 0
        else:
            minmnumbburn = 200000
        gdat.numbburn = min(minmnumbburn, gdat.numbswep - 1)

    # factor by which to thin the sweeps to get samples
    if gdat.factthin == None:
        # temp
        # gdat.factthin = min(2 * gdat.numbpara, gdat.numbswep - gdat.numbburn)
        gdat.factthin = int(min(0.1 * gdat.numbpara, gdat.numbswep - gdat.numbburn))

    # run tag
    gdat.rtag = retr_rtag(gdat)
    
    # plot paths
    gdat.pathplot = gdat.pathimag + 'pcat_' + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    gdat.pathdiag = gdat.pathplot + 'diag/'
    gdat.pathfram = gdat.pathplot + 'fram/'
    if gdat.makeplot:
        os.system('mkdir -p %s' % gdat.pathdiag)
        os.system('mkdir -p %s' % gdat.pathplot)
        os.system('mkdir -p %s' % gdat.pathfram)

    # get the experimental catalog
    if gdat.exprinfo:
        if gdat.exprtype == 'ferm':
            retr_fermdata(gdat)
        gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
        gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)

    # flag to indicate whether information from a deterministic catalog will be used or not
    gdat.trueinfo = gdat.datatype == 'mock' or gdat.exprinfo
    
    # get count data
    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = gdat.exprflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]
    
    # load mock catalog into the reference catalog data structure
    if gdat.datatype == 'mock':
        if gdat.trueinfo:
            gdat.truemeanpnts = gdat.mocknumbpnts
        
            gdat.truelgal = []
            gdat.truebgal = []
            gdat.truesind = []
            for l in gdat.indxpopl:
                gdat.truelgal.append(gdat.mocklgal[l])
                gdat.truebgal.append(gdat.mockbgal[l])
                gdat.truesind.append(gdat.mocksind[l])
                    
            gdat.truestrg = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgclss = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgassc = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truenumbpnts = gdat.mocknumbpnts
            
            if gdat.mockfluxdisttype[l] == 'powr':
                gdat.truefluxdistslop = gdat.mockfluxdistslop
            if gdat.mockfluxdisttype[l] == 'brok':
                gdat.truefluxdistbrek = gdat.mockfluxdistbrek
                gdat.truefluxdistsloplowr = gdat.mockfluxdistsloplowr
                gdat.truefluxdistslopuppr = gdat.mockfluxdistslopuppr
            gdat.truecnts = gdat.mockcnts
               
            gdat.truespec = []
            for l in gdat.indxpopl:
                gdat.truespectemp = empty((3, gdat.numbener, gdat.mocknumbpnts[l]))
                gdat.truespectemp[:] = gdat.mockspec[l][None, :, :]
                gdat.truespec.append(gdat.truespectemp)
            
            gdat.truepsfipara = gdat.mockpsfipara
            gdat.truepsfntype = gdat.mockpsfntype
            gdat.truenormback = gdat.mocknormback
            gdat.datacnts = gdat.mockdatacnts
        
    # spatially averaged background flux 
    gdat.backfluxmean = zeros((gdat.numbback, gdat.numbener))
    for c in gdat.indxback:
        for i in gdat.indxener:
            gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])

    # proposals
    retr_propmodl(gdat)
    
    # factors in the prior expression
    gdat.priofactmeanpnts = log(1. / (log(gdat.maxmmeanpnts) - log(gdat.minmmeanpnts)))
    gdat.priofactlgalbgal = 2. * log(1. / 2. / gdat.maxmgang)
    gdat.priofactfluxdistslop = gdat.numbener * log(1. / (arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)))
    gdat.priofactsplt = -2. * log(2. * gdat.maxmgangmarg) + log(gdat.radispmr) + log(2. * pi)
    # temp -- brok terms are missing

    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc

    # determine proposal probabilities
    gdat.probpropminm = copy(gdat.probprop)
    gdat.probpropmaxm = copy(gdat.probprop)
   
    ## do not allow death or merge when the number of PS is at its minimum or 
    ## births or splits when the number of PS is at its maximum
    gdat.probpropminm[[gdat.indxpropdeth, gdat.indxpropmerg]] = 0.
    gdat.probpropmaxm[[gdat.indxpropbrth, gdat.indxpropsplt]] = 0.
    gdat.probprop /= sum(gdat.probprop)
    gdat.probpropmaxm /= sum(gdat.probpropmaxm)
    gdat.probpropminm /= sum(gdat.probpropminm)
    
    gdat.minmlgalmarg = -gdat.maxmgangmarg
    gdat.maxmlgalmarg = gdat.maxmgangmarg
    gdat.minmbgalmarg = -gdat.maxmgangmarg
    gdat.maxmbgalmarg = gdat.maxmgangmarg

    # construct lists of possible changes to the number of PS for each PS model and the associated probabilities
    gdat.listnumbpntsmodi = []
    gdat.probnumbpntsmodi = []
    for k in range(sum(gdat.maxmnumbpnts)):
        gdat.listnumbpntsmodi.append(arange(1, k + 1))
        gdat.probnumbpntsmodi.append(1. / gdat.listnumbpntsmodi[k])
        gdat.probnumbpntsmodi[k] /= sum(gdat.probnumbpntsmodi[k])
   
    if gdat.verbtype > 1:
        print 'listnumbpntsmodi'
        print gdat.listnumbpntsmodi
        print 'probnumbpntsmodi'
        print gdat.probnumbpntsmodi
        print

    # temp
    gdat.tracsamp = False
    
    # plot settings
    ## marker opacity
    gdat.alphmrkr = 0.5
    gdat.alphpnts = 1.
    gdat.alphmaps = 1.

    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 500
    ## ROI
    gdat.exttrofi = array([gdat.minmlgal, gdat.maxmlgal, gdat.minmbgal, gdat.maxmbgal])
    gdat.exttrofi *= gdat.anglfact 
    gdat.frambndr = gdat.maxmgang * gdat.anglfact
    gdat.frambndrmarg = gdat.maxmgangmarg * gdat.anglfact
    
    # convenience variables
    gdat.fluxhistmodl = empty(gdat.numbflux)
    
    if gdat.exprtype == 'ferm':
        gdat.numbfluxprox = 3
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.numbfluxprox = 1
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
    gdat.maxmangleval = empty(gdat.numbfluxprox)
    for h in gdat.indxfluxprox:
        if gdat.specfraceval == 0:
            gdat.maxmangleval[h] = 3. * gdat.maxmgangmarg
        else:
            if gdat.exprtype == 'ferm':
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.fermpsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
            if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
                gdat.maxmangleval[h] = 10. / gdat.anglfact

    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(1000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)

    factener = (gdat.meanener / gdat.meanener[gdat.indxenerfluxdist[0]])**(-2.)
    
    # limits on counts, which are used to bin or overplot PS counts 
    gdat.minmcnts = gdat.minmflux * mean(mean(gdat.expo, 1), 1) * gdat.diffener * factener
    gdat.maxmcnts = gdat.maxmflux * mean(mean(gdat.expo, 1), 1) * gdat.diffener * factener
    gdat.binscnts = zeros((gdat.numbener, gdat.numbflux + 1))
    for i in gdat.indxener:
        gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbflux + 1) # [1]
       
    ## Real data
    # true data
    if gdat.trueinfo:
        if gdat.datatype == 'inpt':
            gdat.truenumbpnts = None
            gdat.truemeanpnts = None
            gdat.truenormback = None
    
            if gdat.exprtype == 'ferm':
                gdat.truenumbpnts = array([gdat.exprnumbpnts], dtype=int)
                gdat.truelgal = [gdat.exprlgal]
                gdat.truebgal = [gdat.exprbgal]
                gdat.truespec = [gdat.exprspec]
                gdat.truecnts = [gdat.exprcnts]
                gdat.truesind = [gdat.exprsind]
                gdat.truestrg = [gdat.exprstrg]
                gdat.truestrgclss = [gdat.exprstrgclss]
                gdat.truestrgassc = [gdat.exprstrgassc]
                gdat.indxtruevari = [gdat.indxexprvari]
                gdat.truepsfipara = gdat.fermpsfipara
                gdat.truepsfntype = 'doubking'
                gdat.truepsfn = gdat.fermpsfn
            elif gdat.exprtype == 'sdss':
                gdat.truepsfn = gdat.sdsspsfn

        if gdat.datatype == 'mock':
            gdat.truepsfn = retr_psfn(gdat, gdat.truepsfipara, gdat.indxener, gdat.binsangl, gdat.mockpsfntype)
            
        if gdat.truepsfn != None:
            gdat.truefwhm = 2. * retr_psfnwdth(gdat, gdat.truepsfn, 0.5)
        
        truebackcnts = []
        gdat.truesigm = []
        for l in gdat.indxpopl:
            indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
            truebackcntstemp = zeros((gdat.numbener, gdat.truenumbpnts[l], gdat.numbevtt))
            for c in gdat.indxback:
                truebackcntstemp += gdat.backflux[c][:, indxpixltemp, :] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None] * pi * \
                                                                                                                            gdat.truefwhm[:, None, :]**2 / 4.
            truebackcnts.append(truebackcntstemp)
            gdat.truesigm.append(gdat.truecnts[l] / sqrt(truebackcntstemp))
    
        for l in gdat.indxpopl:
            if not (isfinite(gdat.truespec[l]).all() and isfinite(gdat.truesind[l]).all()):
                print 'gdat.truespec'
                print gdat.truespec
                print 'gdat.truesind'
                print gdat.truesind
                raise Exception('True PS parameters are not finite.')

    else:
        gdat.truepsfn = None
    
    # determine the indices of true point sources, which will be compared againts the model sources
    if gdat.trueinfo:
        gdat.indxtruepntscomp = []
        for l in gdat.indxpopl:
            indxtruepntstemp = where((fabs(gdat.truelgal[l]) < gdat.maxmgangcomp) & (fabs(gdat.truebgal[l]) < gdat.maxmgangcomp))[0]
            gdat.indxtruepntscomp.append(indxtruepntstemp)

    # sanity checks
    # temp
    if (fabs(gdat.datacnts - rint(gdat.datacnts)) > 1e-3).any():
        print 'Fractional counts!'

    if amin(gdat.datacnts) < 0.:
        print 'Negative counts!'

    gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix / gdat.diffener
    gdat.datacntsmean = mean(sum(gdat.datacnts, 2), 1)
    gdat.datacntssatu = ceil((amax(sum(gdat.datacnts, 2), 1) - gdat.datacntsmean) * 0.05 + gdat.datacntsmean)
    gdat.resicntssatu = ceil(gdat.datacntssatu * 0.2)
    
    # auxiliary variables for plots
    # temp
    if False:
        if gdat.pixltype == 'heal':
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    gdat.datacntscarttemp = tdpy.util.retr_cart(gdat.datacnts[i, :, m], gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                        minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                    if i == 0 and m == 0:
                        gdat.datacntscart = zeros((gdat.datacntscarttemp.shape[0], gdat.datacntscarttemp.shape[1], gdat.numbener, gdat.numbevtt))
                    gdat.datacntscart[:, :, i, m] = gdat.datacntscarttemp
        if gdat.pixltype == 'cart':
            gdat.datacntscart = zeros((gdat.numbener, gdat.numbpixlfull, gdat.numbevtt)) 
            gdat.datacntscart[gdat.indxcuberofi] = gdat.datacnts
            gdat.datacntscart = gdat.datacntscart.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
            gdat.datacntscart = swapaxes(swapaxes(gdat.datacntscart, 0, 2), 0, 1)

        for i in gdat.indxener:
            indxdatacntscartsatu = where(gdat.datacntscart[:, :, i, :] > gdat.datacntssatu[i])
        gdat.datacntscart[indxdatacntscartsatu[0], indxdatacntscartsatu[1], i, indxdatacntscartsatu[2]] = gdat.datacntssatu[i]

    # make a look-up table of nearby pixels for each pixel
    path = gdat.pathdata + 'indxpixlprox/'
    os.system('mkdir -p %s' % path)
    path += 'indxpixlprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), 1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
    if gdat.verbtype > 0:
        print 'PSF evaluation will be performed up to %4.3g %s for the largest flux.' % (amax(gdat.maxmangleval) * gdat.anglfact, gdat.strganglunittext)
    if os.path.isfile(path):
        if gdat.verbtype > 0:
            print 'Previously computed nearby pixel look-up table will be used.'
            print 'Reading %s...' % path
        fobj = open(path, 'rb')
        gdat.indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        if gdat.verbtype > 0:
            print 'Computing the look-up table...'
        gdat.indxpixlprox = [[] for h in range(gdat.numbfluxprox)]
        for j in gdat.indxpixl:
            dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
            dist[j] = 0.
            for h in range(gdat.numbfluxprox):
                indxpixlproxtemp = where(dist < gdat.maxmangleval[h])[0]
                gdat.indxpixlprox[h].append(indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()

    if gdat.verbtype > 1:
        print 'Memory budget: indxpixlprox'
        totl = 0.
        for h in gdat.indxfluxprox:
            for n in gdat.indxpixl:
                totl += sys.getsizeof(gdat.indxpixlprox[h][n]) / 2.**20
        print '%.4g MB' % totl


def init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, strgplot, indxpoplplot=None):

    figr, axis = plt.subplots(figsize=(1.3 * gdat.plotsize, 1.3 * gdat.plotsize))
    axis.set_xlabel(gdat.longlabl)
    axis.set_ylabel(gdat.latilabl)
    axis.set_xlim([gdat.frambndrmarg, -gdat.frambndrmarg])
    axis.set_ylim([-gdat.frambndrmarg, gdat.frambndrmarg])
    if indxevttplot == None:
        if indxenerplot == None:
            axis.set_title('')
        else:
            axis.set_title(gdat.binsenerstrg[indxenerplot])
    else:
        titl = gdat.binsenerstrg[indxenerplot]
        if gdat.exprtype == 'ferm':
            titl += ', ' + gdat.evttstrg[indxevttplot]
        axis.set_title(titl)

    axis.axvline(gdat.frambndr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axvline(-gdat.frambndr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axhline(gdat.frambndr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axhline(-gdat.frambndr, ls='--', alpha=gdat.alphmrkr, color='black')

    if indxpoplplot == None:
        strg = ''
    else:
        strg = '_pop%d' % indxpoplplot

    if indxevttplot == None:
        path = gdat.pathplot + 'fram/' + strgplot + strg + '_%dA_swep%09d.pdf' % (gdat.indxenerincl[indxenerplot], gdatmodi.cntrswep)
    else:
        path = gdat.pathplot + 'fram/' + strgplot + strg + '_%d%d_swep%09d.pdf' % (gdat.indxenerincl[indxenerplot], gdat.indxevttincl[indxevttplot], gdatmodi.cntrswep)
    
    return figr, axis, path


def supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot):

    # temp
    if True:
        if gdat.trueinfo:
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                        label='Sample', marker='+', linewidth=2, color='b')
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                        label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                        label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                        label=gdat.truelablhits, marker='*', linewidth=2, color='g')
        axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=2)
        
    # true catalog
    if gdat.trueinfo:
        ## get the true catalog
        mrkrsize = retr_mrkrsize(gdat, gdat.truespec[indxpoplplot][0, gdat.indxenerfluxdist, :].flatten())
        lgal = copy(gdat.truelgal[indxpoplplot])
        bgal = copy(gdat.truebgal[indxpoplplot])
        numbpnts = int(gdat.truenumbpnts[indxpoplplot])
        
        ## associations
        ### missed
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].miss
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
        
        ### biased
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].bias[indxenerplot]
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                    label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
        
        ### hit
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].hits[indxenerplot]
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablhits, marker='*', linewidth=2, color='g')
        
        ## annotate
        if gdat.anotcatl:
            for a in range(numbpnts):
                strg = ''
                if gdat.truestrg[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrg[indxpoplplot][a]
                if gdat.truestrgassc[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrgassc[indxpoplplot][a]
                if gdat.truestrgclss[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrgclss[indxpoplplot][a]
                if strg != '':
                    axis.text(gdat.anglfact * gdat.truelgal[indxpoplplot][a], gdat.anglfact * gdat.truebgal[indxpoplplot][a] - gdat.offstext, strg, \
                                                                                                                        ha='center', va='center', color='g', fontsize=6)

    # model catalog
    mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[indxpoplplot][gdat.indxenerfluxdist, :]])
    lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[indxpoplplot]]
    bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[indxpoplplot]]
    axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphpnts, label='Sample', marker='+', linewidth=2, color='b')


def retr_levi(listllik):
    
    minmlistllik = amin(listllik)
    levi = log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    return levi


def retr_info(listllik, levi):
    
    info = mean(listllik) - levi

    return info


def retr_imag(gdat, axis, maps, thisindxener, thisindxevtt, cmap='Reds', mean=False, satuuppr=None, satulowr=None, titl='', minmcbar=None, maxmcbar=None):

    # filter the map
    if thisindxevtt == None:
        if thisindxener == None:
            if mean:
                maps = sum(maps * gdat.expo, axis=2) / sum(gdat.expo, axis=2)
            else:
                maps = sum(maps, axis=2)
        else:
            if mean:
                maps = sum(maps[thisindxener, :, :] * gdat.expo[thisindxener, :, :], axis=1) / sum(gdat.expo[thisindxener, :, :], axis=1)
            else:
                maps = sum(maps[thisindxener, :, :], axis=1)
    else:
        maps = maps[thisindxener, :, thisindxevtt]
    
    # temp
    maps = maps.squeeze()
    
    # project the map to 2D
    if gdat.pixltype == 'heal':
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
    
    if gdat.pixltype == 'cart':
        mapstemp = empty(gdat.numbsidecart**2)
        mapstemp[gdat.indxpixlrofi] = maps
        maps = mapstemp.reshape((gdat.numbsidecart, gdat.numbsidecart)).T
    
    # saturate the map
    if gdat.scalmaps == 'linrsatu':
        if satulowr != None:
            maps[where(maps < satulowr[thisindxener])] = satulowr[thisindxener]
        if satuuppr != None:
            maps[where(maps > satuuppr[thisindxener])] = satuuppr[thisindxener]
    if gdat.scalmaps == 'asnh':
        maps = arcsinh(maps)
   
    # temp
    #imag = axis.imshow(maps)
    if minmcbar != None and maxmcbar != None:
        clim = [minmcbar, maxmcbar]
    else:
        clim = None
    imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', clim=clim, alpha=gdat.alphmaps)
    axis.set_title(titl)

    # make a color bar
    if thisindxevtt != None or thisindxener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    return axis, cbar


def retr_jcbn():
   
    fluxpare, lgalpare, bgalpare, sindpare, fluxauxi, radiauxi, anglauxi, sindauxi \
                                                            = sympy.symbols('fluxpare lgalpare bgalpare sindpare fluxauxi radiauxi anglauxi sindauxi')
    
    matr = sympy.Matrix([[     fluxauxi, 0, 0 , 0,                         fluxpare,                                    0,                                               0, 0], \
                         [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (fluxauxi - 1) * radiauxi * sympy.sin(anglauxi), 0], \
                         [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (1 - fluxauxi) * radiauxi * sympy.cos(anglauxi), 0], \
                         [            0, 0, 0 , 1,                                0,                                    0,                                               0, 0], \
                         [ 1 - fluxauxi, 0, 0 , 0,                        -fluxpare,                                    0,                                               0, 0], \
                         [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), fluxauxi * radiauxi * sympy.sin(anglauxi), 0], \
                         [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), -fluxauxi * radiauxi * sympy.cos(anglauxi), 0], \
                         [            0, 0, 0 , 0,                               0 ,                                    0,                                               0, 1]])

    jcbn = matr.det()
    print jcbn

    return jcbn


def corr_catl(gdat, thisindxpopl, modllgal, modlbgal, modlspec):

    indxtruepntsassc = tdpy.util.gdatstrt()
    indxtruepntsassc.miss = []
    indxtruepntsassc.bias = [[] for i in gdat.indxener]
    indxtruepntsassc.hits = [[] for i in gdat.indxener]
    indxtruepntsassc.mult = []
        
    indxmodlpnts = zeros_like(gdat.truelgal[thisindxpopl], dtype=int) - 1
    specassc = zeros((gdat.numbener, gdat.truenumbpnts[thisindxpopl]), dtype=float)
    numbassc = zeros_like(gdat.truelgal[thisindxpopl], dtype=int)
    distassc = zeros_like(gdat.truelgal[thisindxpopl]) + 3 * gdat.maxmgang
    dir1 = array([gdat.truelgal[thisindxpopl], gdat.truebgal[thisindxpopl]])

    for k in range(modllgal.size):
        dir2 = array([modllgal[k], modlbgal[k]])
        dist = angdist(dir1, dir2, lonlat=True)
        thisindxtruepnts = where(dist < gdat.anglassc)[0]
        
        if thisindxtruepnts.size > 0:
            
            # if there are multiple associated true PS, sort them
            indx = argsort(dist[thisindxtruepnts])
            dist = dist[thisindxtruepnts][indx]
            thisindxtruepnts = thisindxtruepnts[indx]
                
            # store the index of the model PS
            numbassc[thisindxtruepnts[0]] += 1
            if dist[0] < distassc[thisindxtruepnts[0]]:
                specassc[:, thisindxtruepnts[0]] = modlspec[:, k]
                distassc[thisindxtruepnts[0]] = dist[0]
                indxmodlpnts[thisindxtruepnts[0]] = k

    # get the flux limit that delineates the biased associations and hits 
    fluxbias = empty((2, gdat.numbener, gdat.truenumbpnts[thisindxpopl]))
    for i in gdat.indxener:
        fluxbias[:, i, :] = retr_fluxbias(gdat, gdat.truespec[thisindxpopl][0, i, :], i)

    # divide associations into subgroups
    for k in range(gdat.truenumbpnts[thisindxpopl]):
        if numbassc[k] == 0:
            indxtruepntsassc.miss.append(k)
        else:
            if numbassc[k] > 1:
                indxtruepntsassc.mult.append(k)
    
            ## check whether the flux of the associated model point source matches well with the flux of the deterministic point source
            for i in gdat.indxener:
                boolbias = specassc[i, k] > fluxbias[1, i, k] or specassc[i, k] < fluxbias[0, i, k]
                if boolbias:
                    indxtruepntsassc.bias[i].append(k)
                else:
                    indxtruepntsassc.hits[i].append(k)
   
    if gdat.verbtype > 1:
        print 'Correlating catalogs...'
        print 'thisindxpopl'
        print thisindxpopl
        print 'indxtruepntsassc.hits'
        print indxtruepntsassc.hits
        print 'indxtruepntsassc.bias'
        print indxtruepntsassc.bias
        print 'indxtruepntsassc.mult'
        print indxtruepntsassc.mult
        print 'indxtruepntsassc.miss'
        print indxtruepntsassc.miss
        print 

    return indxmodlpnts, indxtruepntsassc


def retr_indxpntscomp(gdat, lgal, bgal):

    indxpntscomp = where((fabs(lgal) < gdat.maxmgangcomp) & (fabs(bgal) < gdat.maxmgangcomp))[0]

    return indxpntscomp


def retr_fluxbias(gdat, spec, indxenerthis):

    # convenience variables
    numbpnts = spec.size
    minmflux = gdat.minmspec[indxenerthis]
    maxmflux = gdat.maxmspec[indxenerthis]

    # tolerance factor at the minimum flux
    factlowr = 5. * ones(numbpnts)

    # tolerance factor at the maximum flux
    factuppr = 1.1 * ones(numbpnts)
    
    # if the flux of interest is above the maximum, i.e., for the non-pivot energy bins, extrapolate the bias lines parallel to the diagonal
    indxspeccons = where(spec > maxmflux)[0]
    factlowr[indxspeccons] = 1.1
    factuppr[indxspeccons] = 1.1

    # calculate the bias lines
    slop = (log(factuppr) - log(factlowr)) / (log(minmflux) - log(maxmflux))
    offs = log(factuppr) + slop * log(maxmflux)

    fluxbias = empty((2, numbpnts))
    fluxbias[0, :] = exp((1. + slop) * log(spec) - offs)
    fluxbias[1, :] = exp((1. - slop) * log(spec) + offs)

    return fluxbias


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
        axis.legend()
        plt.tight_layout()
        pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/powrpdfn/'
        os.system('mkdir -p ' + pathfold)
        figr.savefig(pathfold + 'powrpdfn%04d.pdf' % n)
        plt.close(figr)
        

