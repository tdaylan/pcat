# common imports
from __init__ import *

class gdatstrt(object):
    
    def __init__(self):
        pass


def retr_psfnwdth(gdat, psfn, frac):
    
    fwhm = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            indxanglgood = argsort(psfn[i, :, m])
            intpfwhm = max(frac * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, indxanglgood, m]) and intpfwhm < amax(psfn[i, indxanglgood, m]):
                fwhm[i, m] = interp1d(psfn[i, indxanglgood, m], gdat.binsangl[indxanglgood])(intpfwhm)
    return fwhm


def retr_spec(gdat, flux, sind):
    
    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])

    spec = flux[None, :] * (gdat.meanener[:, None] / gdat.enerfdfn)**(-sind[None, :])
    
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
    
    pntsflux = empty((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
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
                pntsflux[k, i, indxpixltemp, m] = spec[i, k] * psfn[i, :, m]

    # sum contributions from all PS
    pntsfluxtemp = sum(pntsflux, 0) 

    return pntsfluxtemp


def retr_rofi_flux(gdat, normback, pntsflux, tempindx):

    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += normback[c, :, None, None] * gdat.backflux[c][tempindx]        
    
    return modlflux


def cdfn_flux_brok(gdat, flux, fdfnbrek, fdfnsloplowr, fdfnslopuppr):

    norm = 1. / (fdfnbrek**fdfnsloplowr * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) + \
                 fdfnbrek**fdfnslopuppr * (gdat.maxmflux**(1. - fdfnslopuppr) - fdfnbrek**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    fluxunit = norm / (1. - fdfnsloplowr) * fdfnbrek**fdfnsloplowr * (flux**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr))
    indxflux = where(flux >= fdfnbrek)[0]
    
    if indxflux.size > 0:
        temp = norm * fdfnbrek**fdfnsloplowr / (1. - fdfnsloplowr) * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr))
        fluxunit[indxflux] = temp + norm / (1. - fdfnslopuppr) * fdfnbrek**fdfnslopuppr * (flux[indxflux]**(1. - fdfnslopuppr) - fdfnbrek**(1. - fdfnslopuppr))

    return fluxunit


def pdfn_flux_brok(gdat, flux, fdfnbrek, fdfnsloplowr, fdfnslopuppr):

    norm = 1. / (fdfnbrek**fdfnsloplowr * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) + \
                 fdfnbrek**fdfnslopuppr * (gdat.maxmflux**(1. - fdfnslopuppr) - fdfnbrek**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    pdfn = norm * (flux / fdfnbrek)**(-fdfnsloplowr)
    indxflux = where(flux >= fdfnbrek)[0]
    
    if indxflux.size > 0:
        pdfn[indxflux] = norm * (flux[indxflux] / fdfnbrek)**(-fdfnslopuppr)
        
    return pdfn


def icdf_flux_brok(gdat, fluxunit, fdfnbrek, fdfnsloplowr, fdfnslopuppr):
   
    norm = 1. / (fdfnbrek**fdfnsloplowr * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) + \
                 fdfnbrek**fdfnslopuppr * (gdat.maxmflux**(1. - fdfnslopuppr) - fdfnbrek**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    fluxunitbrek = norm / (1. - fdfnsloplowr) * fdfnbrek**fdfnsloplowr * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr))
    flux = (fluxunit * (1. - fdfnsloplowr) / norm / fdfnbrek**fdfnsloplowr + gdat.minmflux**(1. - fdfnsloplowr))**(1. / (1. - fdfnsloplowr))
    indxfluxunit = where(fluxunit >= fluxunitbrek)[0]
    
    if indxfluxunit.size > 0:
        temp = norm * fdfnbrek**fdfnsloplowr / (1. - fdfnsloplowr) * (fdfnbrek**(1. - fdfnsloplowr) - gdat.minmflux**(1. - fdfnsloplowr))
        flux[indxfluxunit] = ((fluxunit[indxfluxunit] - temp) * (1. - fdfnslopuppr) / norm / fdfnbrek**fdfnslopuppr + fdfnbrek**(1. - fdfnslopuppr))**(1. / (1. - fdfnslopuppr))

    return flux


def cdfn_flux_powr(gdat, flux, fdfnslop):
        
    fluxunit = (flux**(1. - fdfnslop) - gdat.minmflux**(1. - fdfnslop)) / (gdat.maxmflux**(1. - fdfnslop) - gdat.minmflux**(1. - fdfnslop))
        
    return fluxunit


def icdf_flux_powr(gdat, fluxunit, fdfnslop):

    flux = (fluxunit * (gdat.maxmflux**(1. - fdfnslop) - gdat.minmflux**(1. - fdfnslop)) + gdat.minmflux**(1. - fdfnslop))**(1. / (1. - fdfnslop))
    
    return flux


def pdfn_flux_powr(gdat, flux, fdfnslop):
  
    norm = (1. - fdfnslop) / (gdat.maxmflux**(1. - fdfnslop) - gdat.minmflux**(1. - fdfnslop))
    
    pdfn = norm * flux**(-fdfnslop)
    
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
        indxpixl = gdat.pixlcnvt[ang2pix(gdat.numbsideheal, deg2rad(90. - bgal), deg2rad(lgal))]
        
        if (indxpixl == -1).any():  
            print 'pixlcnvt went negative!'
    else:
        
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
            
        indxpixl = indxbgcr * gdat.numbsidecart + indxlgcr

    return indxpixl


def show_memo(objt, name):

    if isinstance(objt, list):
        for k in len(objt):
            size = sys.getsizeof(objt[k]) / 2.**20
    else:
        listsize = []
        listattr = []
        for attr, valu in objt.__dict__.iteritems():
            listsize.append(sys.getsizeof(valu) / 2.**20)
            listattr.append(attr)
        size = array(listsize)
        attr = array(listattr)
        sizetotl = sum(size) 
        print 'Memory budget: %s' % name
        print 'Total size: %.4g MB' % sizetotl
        
        # sort the sizes to get the largest tail
        indxsizetemp = argsort(size)[::-1]
        size = size[indxsizetemp]
        attr = attr[indxsizetemp]
        print 'Largest 5:'
        for k in range(5):
            print '%s: %.4g MB' % (attr[k], size[k])
        print 


def retr_llik(gdat, gdatmodi, init=False):

    if init:

        gdatmodi.thisllik = gdat.datacnts * log(gdatmodi.thismodlcnts) - gdatmodi.thismodlcnts
        
    elif gdatmodi.thisindxprop >= gdat.indxproppsfipara:

        # load convenience variables
        timebegn = time.time()
        
        ## save time by not reloading all PS into modi structure for PSF parameter updates
        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            lgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsamplgal)]
            bgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampbgal)]
            spec = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenermodi, :]]
        if gdatmodi.thisindxprop >= gdat.indxpropbrth:
            lgal = gdatmodi.modilgal[:gdat.numbmodipnts]
            bgal = gdatmodi.modibgal[:gdat.numbmodipnts]
            spec = gdatmodi.modispec[meshgrid(gdat.indxenermodi, arange(gdat.numbmodipnts), indexing='ij')]
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timebegn
        
        # determine pixels over which to evaluate the log-likelihood
        timebegn = time.time()
        
        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            gdat.indxpixlmodi = gdat.indxpixl
        if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
            thisindxpixlprox = []
            for k in range(gdat.numbmodipnts):
                
                # temp -- this may not work for extreme color PS!
                # take the flux at the pivot energy
                if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                    fluxtemp = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenerfdfn, k]]
                else:
                    fluxtemp = gdatmodi.modispec[gdat.indxenerfdfn, k]

                # find the flux index
                indxfluxproxtemp = amin(where(gdat.binsfluxprox - fabs(fluxtemp) > 0.)[0]) - 1
                indxpixltemp = retr_indxpixl(gdat, bgal[k], lgal[k])
                thisindxpixlprox.append(indxpixlprox[indxfluxproxtemp][indxpixltemp])
            gdat.indxpixlmodi = unique(concatenate(thisindxpixlprox))
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timebegn

        # construct the mesh grid for likelihood evaluation
        timebegn = time.time()
        
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            gdat.indxcubemodi = meshgrid(gdat.indxenermodi, gdat.indxpixlmodi, gdat.indxevtt, indexing='ij')

        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 2] = timefinl - timebegn

        # update the model point source flux map, if needed
        timebegn = time.time()
        
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            
            # temp
            #if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            #    gdatmodi.nextpntsflux[gdat.indxcubemodi] = 0.
            #else:
            gdatmodi.nextpntsflux[gdat.indxcubemodi] = copy(gdatmodi.thispntsflux[gdat.indxcubemodi])
                
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

                for k in range(gdat.numbmodipnts):
                    
                    # calculate the distance to the pixels to be updated
                    dist = retr_angldistunit(gdat, lgal[k], bgal[k], thisindxpixlprox[k])

                    # interpolate the PSF
                    psfn = psfnintp(dist)
    
                    # temp
                    if False and gdat.strgcnfg == 'cnfg_test':
                        print 'hey'
                        print 'n, k'
                        print n, k
                        print 'thispsfn'
                        print mean(gdatmodi.thispsfnintp(dist), 1)
                        print 'psfn'
                        print mean(psfn[:, :], 1)
                        print

                    # add the contribution of the PS to the the proposed flux map
                    for i in range(gdat.indxenermodi.size):
                        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                            if n == 0:
                                spectemp = -spec[i, k]
                            else:
                                spectemp = spec[i, k]
                        else:
                            spectemp = spec[i, k]
                        gdatmodi.nextpntsflux[gdat.indxenermodi[i], thisindxpixlprox[k], :] += spectemp * psfn[gdat.indxenermodi[i], :, :]

                        if False and gdat.strgcnfg == 'cnfg_test':
                            print 'i: ', i
                            print 'mean(spectemp)'
                            print spectemp
                            print 'mean(psfn[gdat.indxenermodi[i], :, :], 0)'
                            print mean(psfn[gdat.indxenermodi[i], :, :], 0)
                
                if False and gdat.strgcnfg == 'cnfg_test':
                    print 'mean(gdatmodi.nextpntsflux[0, thisindxpixlprox[k], :], 0)'
                    print mean(gdatmodi.nextpntsflux[0, thisindxpixlprox[k], :], 0)
                    print 'mean(gdatmodi.thispntsflux[0, thisindxpixlprox[k], :], 0)'
                    print mean(gdatmodi.thispntsflux[0, thisindxpixlprox[k], :], 0)

        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 3] = timefinl - timebegn

        # update the total model flux map
        timebegn = time.time()
       
        indxtemp = meshgrid(gdat.indxback, gdat.indxenermodi, indexing='ij')

        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            normback = gdatmodi.nextsampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.thispntsflux
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            normback = gdatmodi.thissampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.nextpntsflux
        gdatmodi.nextmodlflux[gdat.indxcubemodi] = retr_rofi_flux(gdat, normback, pntsflux, gdat.indxcubemodi)

        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 4] = timefinl - timebegn

        # calculate the total model count map
        timebegn = time.time()
        
        gdatmodi.nextmodlcnts[gdat.indxcubemodi] = gdatmodi.nextmodlflux[gdat.indxcubemodi] * gdat.expo[gdat.indxcubemodi] * \
            gdat.apix * gdat.diffener[gdat.indxenermodi, None, None] # [1]
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 5] = timefinl - timebegn

        # calculate the likelihood
        timebegn = time.time()
        
        gdatmodi.nextllik[gdat.indxcubemodi] = gdat.datacnts[gdat.indxcubemodi] * log(gdatmodi.nextmodlcnts[gdat.indxcubemodi]) - gdatmodi.nextmodlcnts[gdat.indxcubemodi]
            
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 6] = timefinl - timebegn

        if not isfinite(gdatmodi.nextllik[gdat.indxcubemodi]).any():
            warnings.warn('Log-likelihood went NAN!')
            
        gdatmodi.deltllik = sum(gdatmodi.nextllik[gdat.indxcubemodi] - gdatmodi.thisllik[gdat.indxcubemodi])
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


def retr_lpri(gdat, gdatmodi, init=False):
        
    if init:
        for l in gdat.indxpopl:
            fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[l]]
            if gdat.fdfntype == 'powr':
                fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[l]]
                fluxhistmodl = fdfnnorm * pdfn_flux_powr(gdat, gdat.meanflux, fdfnslop) * gdat.diffflux
            if gdat.fdfntype == 'brok':
                fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[l]]
                fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[l]]
                fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[l]]
                fluxhistmodl = fdfnnorm * pdfn_flux_brok(gdat, gdat.meanflux, fdfnbrek, fdfnsloplowr, fdfnslopuppr) * gdat.diffflux 
            spec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfdfn[0], :]]
            
            fluxhist = histogram(spec, gdat.binsflux)[0]
            lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)

            gdatmodi.thislpri[l, :] = lprbpois
           
            if gdat.strgcnfg == 'cnfg_test':
                print 'hey'
                print 'init'
                print 'fluxhist'
                print fluxhist
                print 'fluxhistmodl'
                print fluxhistmodl
                print 'gdatmodi.thislpri'
                print gdatmodi.thislpri
                print
    
        gdatmodi.nextlpri = copy(gdatmodi.thislpri)
                
    else:
        
        # determine if either the number of PS or any of the hyperpriors is being updated
        if gdatmodi.thisindxprop == gdat.indxpropfdfnnorm or gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg:
            thisbool = True
        else:
            thisbool = False
        if gdat.fdfntype == 'powr':
            if gdatmodi.thisindxprop == gdat.indxpropfdfnslop:
                thisbool = True
        if gdat.fdfntype == 'brok':
            if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek or gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr or gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
                thisbool = True
        
        # comptue the delta log-prior for the associated update
        if thisbool:

            # normalization of the FDF
            if gdatmodi.thisindxprop == gdat.indxpropfdfnnorm:
                fdfnnorm = gdatmodi.nextsampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]]
                if gdat.fdfntype == 'powr':
                    fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
                if gdat.fdfntype == 'brok':
                    fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                    fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                    fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            
            # shape of the FDF
            if gdat.fdfntype == 'powr':
                if  gdatmodi.thisindxprop == gdat.indxpropfdfnslop:
                    fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]]
                    fdfnslop = gdatmodi.nextsampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
            if gdat.fdfntype == 'brok':
                if  gdatmodi.thisindxprop == gdat.indxpropfdfnbrek or gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr or gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
                    fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]]
                    if  gdatmodi.thisindxprop == gdat.indxpropfdfnbrek:
                        fdfnbrek = gdatmodi.nextsampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                        fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                        fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                    if  gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr:
                        fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                        fdfnsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                        fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                    if  gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
                        fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                        fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                        fdfnslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            
            # number of PS
            if gdatmodi.thisindxprop == gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxpropdeth or \
                gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg:
                fdfnnorm = gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]]
                if gdat.fdfntype == 'powr':
                    fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
                if gdat.fdfntype == 'brok':
                    fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                    fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                    fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            
            # flux prior
            if gdat.fdfntype == 'powr':
                fluxhistmodl = fdfnnorm * pdfn_flux_powr(gdat, gdat.meanflux, fdfnslop) * gdat.diffflux
            if gdat.fdfntype == 'brok':
                fluxhistmodl = fdfnnorm * pdfn_flux_brok(gdat, gdat.meanflux, fdfnbrek, fdfnsloplowr, fdfnslopuppr) * gdat.diffflux
    
            # model flux distribution
            fluxhist = histogram(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn, :]], gdat.binsflux)[0] 
            if gdatmodi.thisindxprop == gdat.indxpropbrth:
                fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfdfn, 0], gdat.binsflux)[0]
            elif gdatmodi.thisindxprop == gdat.indxpropdeth:
                fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfdfn, 0], gdat.binsflux)[0]
            elif gdatmodi.thisindxprop == gdat.indxpropsplt:
                fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfdfn, 0], gdat.binsflux)[0]
                fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfdfn, 1:3], gdat.binsflux)[0]
            elif gdatmodi.thisindxprop == gdat.indxpropmerg:
                fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfdfn, 0:2], gdat.binsflux)[0]
                fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfdfn, 2], gdat.binsflux)[0]
            
            lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
            gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] = lprbpois

            # temp
            if gdat.strgcnfg == 'cnfg_test' and gdatmodi.cntrswep % 1000 == 0:
                print 'hey'
                print 'fluxhist'
                print fluxhist
                print 'fluxhistmodl'
                print fluxhistmodl
                print 'gdatmodi.thislpri'
                print gdatmodi.thislpri
                print 'gdatmodi.nextlpri'
                print gdatmodi.nextlpri
                print
                print 
                print

            gdatmodi.deltlpri = sum(gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, :])
        else:
            gdatmodi.deltlpri = 0.
        
        
def retr_sampvarb(gdat, indxpntsfull, samp):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxsampnumbpnts] = samp[gdat.indxsampnumbpnts]
    for l in gdat.indxpopl:
        sampvarb[gdat.indxsampfdfnnorm[l]] = icdf_logt(samp[gdat.indxsampfdfnnorm[l]], gdat.minmfdfnnorm[l], gdat.factfdfnnorm[l])
        if gdat.fdfntype == 'powr':
            sampvarb[gdat.indxsampfdfnslop[l]] = icdf_atan(samp[gdat.indxsampfdfnslop[l]], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
        if gdat.fdfntype == 'brok':
            sampvarb[gdat.indxsampfdfnbrek[l]] = icdf_logt(samp[gdat.indxsampfdfnbrek[l]], gdat.minmfdfnbrek[l], gdat.factfdfnbrek[l])
            sampvarb[gdat.indxsampfdfnsloplowr[l]] = icdf_atan(samp[gdat.indxsampfdfnsloplowr[l]], gdat.minmfdfnsloplowr[l], gdat.factfdfnsloplowr[l])
            sampvarb[gdat.indxsampfdfnslopuppr[l]] = icdf_atan(samp[gdat.indxsampfdfnslopuppr[l]], gdat.minmfdfnslopuppr[l], gdat.factfdfnslopuppr[l])

    for k in gdat.indxpsfipara:
        sampvarb[gdat.indxsamppsfipara[k]] = icdf_psfipara(gdat, samp[gdat.indxsamppsfipara[k]], k)
    for c in gdat.indxback:
        sampvarb[gdat.indxsampnormback[c, :]] = icdf_logt(samp[gdat.indxsampnormback[c, :]], gdat.minmnormback[c], gdat.factnormback[c])
    
    for l in gdat.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
        if gdat.fdfntype == 'powr':
            sampvarb[indxsampspec[l][gdat.indxenerfdfn, :]] = icdf_flux_powr(gdat, samp[indxsampspec[l][gdat.indxenerfdfn, :]], sampvarb[gdat.indxsampfdfnslop[l]])
        if gdat.fdfntype == 'brok':
            fluxunit = samp[indxsampspec[l][gdat.indxenerfdfn[0], :]]
            fdfnbrek = sampvarb[gdat.indxsampfdfnbrek[l]]
            fdfnsloplowr = sampvarb[gdat.indxsampfdfnsloplowr[l]]
            fdfnslopuppr = sampvarb[gdat.indxsampfdfnslopuppr[l]]
            sampvarb[indxsampspec[l][gdat.indxenerfdfn, :]] = icdf_flux_brok(gdat, fluxunit, fdfnbrek, fdfnsloplowr, fdfnslopuppr)
        sampvarb[indxsampsind[l]] = icdf_eerr(samp[indxsampsind[l]], gdat.meansdfn[l], gdat.stdvsdfn[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
        sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampspec[l][gdat.indxenerfdfn[0], :]], sampvarb[indxsampsind[l]])
    
    return sampvarb
    

def retr_maps(gdat, indxpntsfull, sampvarb):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
     
    listspectemp = []
    for l in gdat.indxpopl:
        listspectemp.append(sampvarb[indxsampspec[l]])

    pntsflux = retr_pntsflux(gdat, sampvarb[concatenate(indxsamplgal)], sampvarb[concatenate(indxsampbgal)], \
        concatenate(listspectemp, axis=1), sampvarb[gdat.indxsamppsfipara], gdat.psfntype)
    totlflux = retr_rofi_flux(gdat, sampvarb[gdat.indxsampnormback], pntsflux, gdat.indxcubefull)
    
    pntscnts = pntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    totlcnts = totlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    
    return pntsflux, pntscnts, totlflux, totlcnts


def retr_mrkrsize(gdat, flux):

    mrkrsize = (log(flux) - log(gdat.minmflux)) / (log(gdat.maxmflux) - log(gdat.minmflux)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
        
    return mrkrsize


def retr_fermpsfn(gdat):
   
    if False:
        reco = 8
    else:
        reco = 7

    if reco == 8:
        path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_back.fits'
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
                path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 0:
                path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_back.fits'
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
    factener = (10. * gdat.meanener[:, None])**fermscal[None, :, 2]
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * factener)**2 + fermscal[None, :, 1]**2)
    
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
    
    if gdatmodi.thisindxprop == gdat.indxpropfdfnnorm:
        gdatmodi.thissampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]]
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
    
    # determine if a hyperparameter is to be updated
    thisbool = False
    if gdat.fdfntype == 'powr':
        if gdatmodi.thisindxprop == gdat.indxpropfdfnslop:
            thisbool = True
    if gdat.fdfntype == 'brok':
        if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek or gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr or gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
            thisbool = True
    
    # update the hyperparameters
    if thisbool:
        
        ## update the sample vector
        flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn[0], :]]
        if gdat.fdfntype == 'powr':
            fdfnslop = gdatmodi.nextsampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
            gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_powr(gdat, flux, fdfnslop)
        if gdat.fdfntype == 'brok':
            if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek:
                fdfnbrek = gdatmodi.nextsampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
            elif gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr:
                fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                fdfnsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
            else:
                fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                fdfnslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_brok(gdat, flux, fdfnbrek, fdfnsloplowr, fdfnslopuppr)

        ## update the unit sample vector -- this is unique for hyperparameter updates
        gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn, :], -1] = fluxunit
        
        # update the prior register
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]

    # proposals that change the likelihood
    if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
        gdatmodi.thisllik[gdat.indxcubemodi] = gdatmodi.nextllik[gdat.indxcubemodi]
        gdatmodi.thismodlcnts[gdat.indxcubemodi] = gdatmodi.nextmodlcnts[gdat.indxcubemodi]
        
    # PSF
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        # temp
        gdatmodi.thispsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
     
    # background normalization
    if gdatmodi.thisindxprop == gdat.indxpropnormback:
        gdatmodi.thissampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]] = \
            gdatmodi.nextsampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]]
        
    # proposals that change the PS flux map
    if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
        gdatmodi.thispntsflux[gdat.indxcubemodi] = copy(gdatmodi.nextpntsflux[gdat.indxcubemodi])

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
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdat.killindxpnts)
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdat.killindxpnts)

    ## split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:

        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].append(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0])
        del gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]
        
        ### update the components
        #### first component
        gdatmodi.thissampvarb[gdat.indxsampchd0] = gdatmodi.modilgal[1]
        gdatmodi.thissampvarb[gdat.indxsampchd0+1] = gdatmodi.modibgal[1]
        gdatmodi.thissampvarb[gdat.indxsampchd0+2:gdat.indxsampchd0+2+gdat.numbener] = gdatmodi.modispec[:, 1]
        #### second component
        gdatmodi.thissampvarb[gdat.indxsampchd1] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdat.indxsampchd1+1] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdat.indxsampchd1+2:gdat.indxsampchd1+2+gdat.numbener] = gdatmodi.modispec[:, 2]
        
    ## merge
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdat.mergindxchd1)
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdat.mergindxchd1)

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
    
    pairlist = []
    for k in range(lgal.size):
        indxpnts = k + 1 + where((lgal[k+1:] < lgal[k] + gdat.radispmrlbhl) & \
            (lgal[k+1:] > lgal[k] - gdat.radispmrlbhl) & (bgal[k+1:] < bgal[k] + gdat.radispmrlbhl) & (bgal[k+1:] > bgal[k] - gdat.radispmrlbhl))[0]
        for l in range(indxpnts.size):
            pairlist.append([k, indxpnts[l]])
            
    return pairlist


def retr_expr(gdat):
    
    if gdat.exprtype == 'ferm':

        path = os.environ["PCAT_DATA_PATH"] + '/gll_psc_v16.fit'

        fgl3 = pf.getdata(path)
        
        fgl3spectemp = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], fgl3['Flux10000_100000']))
        fgl3spectemp = fgl3spectemp[gdat.indxenerincl, :] / gdat.diffener[:, None]
        
        # sort the catalog in decreasing flux
        indxfgl3sort = argsort(fgl3spectemp[gdat.indxenerfdfn[0], :])[::-1]

        fgl3spectemp = fgl3spectemp[:, indxfgl3sort]

        fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], fgl3['Unc_Flux10000_100000']))
        fgl3specstdvtemp = fgl3specstdvtemp[gdat.indxenerincl, :, :] / gdat.diffener[:, None, None]
        fgl3specstdvtemp = fgl3specstdvtemp[:, indxfgl3sort, :]

        fgl3numbpntsfull = fgl3['glon'].size
        
        fgl3lgal = fgl3['glon'][indxfgl3sort]
        fgl3lgal = ((fgl3lgal - 180.) % 360.) - 180.

        fgl3bgal = fgl3['glat'][indxfgl3sort]
                
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
        
        # adjust 3FGL positions according to the ROI center
        if gdat.regitype == 'ngal':
            rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
            fgl3bgal, fgl3lgal = rad2deg(rttr(deg2rad(90. - fgl3bgal), deg2rad(fgl3lgal)))
            fgl3bgal = 90. - fgl3bgal

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
    
    gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
    gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)
        


def retr_rtag(gdat, indxprocwork):
    
    
    if indxprocwork == None:
        rtag = 'AA_%d_%d_%d_%d_%s_%s_%s' % (gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
    else:
        rtag = '%02d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxsamp, stdv):
    
    if gdat.fracrand > 0.:
        if rand() < gdat.fracrand:
            gdatmodi.drmcsamp[indxsamp, 1] = rand()
        else:
            gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_angldistunit(gdat, lgal1, bgal1, indxpixltemp):
    
    xaxi, yaxi, zaxi = retr_unit(lgal1, bgal1)
    cositemp = gdat.xaxigrid[indxpixltemp] * xaxi + gdat.yaxigrid[indxpixltemp] * yaxi + gdat.zaxigrid[indxpixltemp] * zaxi
    angldist = arccos(cositemp)

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


def retr_strgfluxunit(gdat):
    
    if gdat.exprtype == 'ferm':
        strgfluxunit = r'[1/cm$^2$/s/GeV]'
    if gdat.exprtype == 'sdss':
        strgfluxunit = '[nMgy]'
        
    return strgfluxunit
     
   
def retr_gang(lgal, bgal):
    
    gang = rad2deg(arccos(cos(deg2rad(lgal)) * cos(deg2rad(bgal))))

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
    
    return enerstrg, binsenerstrg


def retr_prop(gdat, gdatmodi):
 
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal,  gdatmodi.thisindxsampspec, \
            gdatmodi.thisindxsampsind, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    if gdat.verbtype > 1:
        print 'retr_prop(): '

        print 'drmcsamp'
        print gdatmodi.drmcsamp
        
        print 'thissampvarb: '
        for k in range(gdatmodi.thissampvarb.size):
            if k == gdat.indxsampcompinit:
                print
            if k > gdat.indxsampcompinit and (k - gdat.indxsampcompinit) % gdat.numbcomp == 0:
                print
            print gdatmodi.thissampvarb[k]
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
    if gdatmodi.thisindxprop == gdat.indxpropfdfnnorm:
        gdatmodi.indxsampmodi = gdat.indxsampfdfnnorm[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvfdfnnorm)
        gdatmodi.nextsampvarb[gdat.indxsampfdfnnorm] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
            gdat.minmfdfnnorm[gdatmodi.indxpoplmodi], gdat.factfdfnnorm[gdatmodi.indxpoplmodi])
        
    # flux distribution function shape
    if gdat.fdfntype == 'powr':
        ## FDF power law slope
        if gdatmodi.thisindxprop == gdat.indxpropfdfnslop:
            gdatmodi.indxsampvarbmodi = gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfdfnslop)
            gdatmodi.nextsampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
                gdat.minmfdfnslop[gdatmodi.indxpoplmodi], gdat.factfdfnslop[gdatmodi.indxpoplmodi])
            gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn[0], :]))
            if gdat.verbtype > 1:
                print 'gdatmodi.indxsampvarbmodi'
                print gdatmodi.indxsampvarbmodi
                print 'thissampvarb[gdat.indxsampfdfnslop]'
                print gdatmodi.thissampvarb[gdat.indxsampfdfnslop]
                print 'nextsampvarb[gdat.indxsampfdfnslop]'
                print gdatmodi.nextsampvarb[gdat.indxsampfdfnslop]
    
    if gdat.fdfntype == 'brok':
        
        ## FDF break flux
        if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek:
            gdatmodi.indxsampvarbmodi = gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvflux)
            gdatmodi.nextsampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
                gdat.minmfdfnbrek[gdatmodi.indxpoplmodi], gdat.factfdfnbrek[gdatmodi.indxpoplmodi])
            if gdat.verbtype > 1:
                print 'thissampvarb[gdat.indxsampfdfnbrek]'
                print gdatmodi.thissampvarb[gdat.indxsampfdfnbrek]
                print 'nextsampvarb[gdat.indxsampfdfnbrek]'
                print gdatmodi.nextsampvarb[gdat.indxsampfdfnbrek]
        
        ## FDF lower power law slope
        if gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr:
            gdatmodi.indxsampvarbmodi = gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfdfnslop)
            gdatmodi.nextsampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
                gdat.minmfdfnsloplowr[gdatmodi.indxpoplmodi], gdat.factfdfnsloplowr[gdatmodi.indxpoplmodi])
            if gdat.verbtype > 1:
                print 'thissampvarb[gdat.indxsampfdfnsloplowr]'
                print gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr]
                print 'nextsampvarb[gdat.indxsampfdfnsloplowr]'
                print gdatmodi.nextsampvarb[gdat.indxsampfdfnsloplowr]
        
        ## FDF upper power law slope
        if gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
            gdatmodi.indxsampvarbmodi = gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfdfnslop)
            gdatmodi.nextsampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
                gdat.minmfdfnslopuppr[gdatmodi.indxpoplmodi], gdat.factfdfnslopuppr[gdatmodi.indxpoplmodi])
            if gdat.verbtype > 1:
                print 'thissampvarb[gdat.indxsampfdfnslopuppr]'
                print gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr]
                print 'nextsampvarb[gdat.indxsampfdfnslopuppr]'
                print gdatmodi.nextsampvarb[gdat.indxsampfdfnslopuppr]
   
        if gdatmodi.thisindxprop == gdat.indxpropfdfnbrek or gdatmodi.thisindxprop == gdat.indxpropfdfnsloplowr or gdatmodi.thisindxprop == gdat.indxpropfdfnslopuppr:
            gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn[0], :]))

    # PSF parameter change 
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        
        # index of the PSF parameter to change
        gdat.indxpsfiparamodi = choice(gdat.indxpsfipara)

        # the energy bin of the PS flux map to be modified
        gdat.indxenermodi = array([(gdat.indxpsfiparamodi % gdat.numbpsfiparaevtt) // gdat.numbformpara])
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdat.indxsamppsfipara[gdat.indxpsfiparamodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvpsfipara)
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara] = copy(gdatmodi.thissampvarb[gdat.indxsamppsfipara])
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara[gdat.indxpsfiparamodi]] = \
            icdf_psfipara(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.indxpsfiparamodi)
            
        gdat.numbmodipnts = int(sum(gdatmodi.thissampvarb[gdat.indxsampnumbpnts]))
                    
        # construct the proposed PSF
        gdatmodi.nextpsfn = retr_psfn(gdat, gdatmodi.nextsampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype)
        
        temp = retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype)
       
        if False:
            print 'hey'
            print 'gdatmodi.thispsfn'
            print mean(gdatmodi.thispsfn, 1)
            print 'mean'
            print mean(temp, 1)
            print 'gdatmodi.nextpsfn'
            print mean(gdatmodi.nextpsfn, 1)
            print
            print


        gdatmodi.nextpsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        
        if gdat.verbtype > 1:
           
            print 'indxpsfiparamodi: ', gdat.indxpsfiparamodi
            print 'indxenermodi: ', gdat.indxenermodi
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
        gdat.indxenermodi = choice(gdat.indxener)
        gdat.indxbackmodi = choice(gdat.indxback)
        gdatmodi.indxsampmodi = gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]
        
        ## propose
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvback)

        ## transform back from the unit space
        gdatmodi.nextsampvarb[gdat.indxsampnormback] = copy(gdatmodi.thissampvarb[gdat.indxsampnormback])
        gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
                                                                                    gdat.minmnormback[gdat.indxbackmodi], gdat.factnormback[gdat.indxbackmodi])

        if gdat.verbtype > 1:
            print 'indxenermodi: ', gdat.indxenermodi
            print 'indxbackmodi: ', gdat.indxbackmodi
            print 'indxsampmodi: ', gdatmodi.indxsampmodi
            print 'thissampvarb[gdat.indxsampnormback]: ', gdatmodi.thissampvarb[gdat.indxsampnormback]
            print 'nextsampvarb[gdat.indxsampnormback]: ', gdatmodi.nextsampvarb[gdat.indxsampnormback]
            print

    # birth
    if gdatmodi.thisindxprop == gdat.indxpropbrth:

        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxsampbrth = int(gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:gdatmodi.indxpoplmodi]) + gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0] * gdat.numbcomp)
        
        # sample auxiliary variables
        numbauxipara = gdat.numbcompcolr
        gdatmodi.auxipara = rand(numbauxipara)

        gdatmodi.drmcsamp[indxsampbrth:indxsampbrth+2, -1] = gdatmodi.auxipara[0:2]
        gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1] = gdatmodi.auxipara[-2]
        gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1] = gdatmodi.auxipara[-1]

        # sample indices to be modified
        gdatmodi.indxsampmodi = arange(indxsampbrth, indxsampbrth + gdat.numbcomp, dtype=int)

        # modification catalog
        gdat.numbmodipnts = 1
        gdatmodi.modilgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcomplgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        gdatmodi.modibgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompbgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        fluxunit = gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1]
        if gdat.fdfntype == 'powr':
            fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_powr(gdat, fluxunit, fdfnslop)
        if gdat.fdfntype == 'brok':
            fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
            fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_brok(gdat, array([fluxunit]), fdfnbrek, fdfnsloplowr, fdfnslopuppr)
        gdatmodi.modisind[0] = icdf_eerr(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1], gdat.meansdfn[gdatmodi.indxpoplmodi], \
                    gdat.stdvsdfn[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        gdatmodi.modispec[:, 0] = retr_spec(gdat, modiflux, gdatmodi.modisind[0]).flatten()
    
        if gdat.verbtype > 1:
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
        gdat.killindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        gdatmodi.indxsampmodi = array([])
            
        # modification catalog
        gdat.numbmodipnts = 1
        gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modispec[:, 0] = -gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, killindxindxpnts]]

        if gdat.verbtype > 1:
            print 'killindxpnts: ', gdat.killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
  
    # split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:
        
        gdat.numbmodipnts = 3
        
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        spltindxindxpnts = choice(thisindxindxpnts)
        spltindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        gdat.indxsampchd0 = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][spltindxindxpnts] * gdat.numbcomp
        indxfinl0 = gdat.indxsampchd0 + gdat.numbcomp
        gdat.indxsampchd1 = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0] * gdat.numbcomp
        indxfinl1 = gdat.indxsampchd1 + gdat.numbcomp
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdat.indxsampchd0, indxfinl0, dtype=int), arange(gdat.indxsampchd1, indxfinl1, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, spltindxindxpnts]]
        thisflux = thisspec[gdat.indxenerfdfn[0]]
        # temp
        #thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn[0], spltindxindxpnts]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        
        if gdat.verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxsampchd0: ', gdat.indxsampchd0
            print 'indxfinl0: ', indxfinl0
            print 'indxsampchd1: ', gdat.indxsampchd1
            print 'indxfinl1: ', indxfinl1
            if gdat.pixltype == 'heal':
                print 'thislgal: ', thislgal
                print 'thisbgal: ', thisbgal
            else:
                print 'thislgal: ', 3600. * thislgal
                print 'thisbgal: ', 3600. * thisbgal
            print 'thisflux: ', thisflux
            print 'thissind: ', thissind
            
        # determine the new components
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.numbcompcolr)
        gdatmodi.auxipara[0:2] = rand(2) * gdat.radispmrlbhl
        if gdat.fdfntype == 'powr':
            fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
            flux = icdf_flux_powr(gdat, rand(), fdfnslop)
        if gdat.fdfntype == 'brok':
            fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
            fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
            fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
            flux = icdf_flux_brok(gdat, rand(), fdfnbrek, fdfnsloplowr, fdfnslopuppr)
        gdatmodi.auxipara[2] = flux
        gdatmodi.auxipara[3] = icdf_eerr(rand(), gdat.meansdfn[gdatmodi.indxpoplmodi], \
                                    gdat.stdvsdfn[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        
        if gdat.verbtype > 1:
            if gdat.pixltype == 'heal':
                print 'auxipara[0]: ', gdatmodi.auxipara[0]
                print 'auxipara[1]: ', gdatmodi.auxipara[1]
            else:
                print 'auxipara[0]: ', 3600. * gdatmodi.auxipara[0]
                print 'auxipara[1]: ', 3600. * gdatmodi.auxipara[1]
            print 'auxipara[2:]'
            print gdatmodi.auxipara[2:]
            print
            
        nextflux0 = thisflux / (thisflux + gdatmodi.auxipara[2])
        nextflux1 = gdatmodi.auxipara[2] / (thisflux + gdatmodi.auxipara[2])
        nextsind0 = thissind
        nextsind1 = gdatmodi.auxipara[3]

        nextlgal0 = thislgal + gdatmodi.auxipara[0]
        nextlgal1 = thislgal - gdatmodi.auxipara[0]
        nextbgal0 = thisbgal + gdatmodi.auxipara[1]
        nextbgal1 = thisbgal - gdatmodi.auxipara[1]
        
        if gdat.verbtype > 1:
            if gdat.pixltype == 'heal':
                print 'nextlgal0: ', nextlgal0
                print 'nextlgal1: ', nextlgal1
                print 'nextbgal0: ', nextbgal0
                print 'nextbgal1: ', nextbgal1
            else:
                print 'nextlgal0: ', 3600. * nextlgal0
                print 'nextlgal1: ', 3600. * nextlgal1
                print 'nextbgal0: ', 3600. * nextbgal0
                print 'nextbgal1: ', 3600. * nextbgal1
            print 'nextflux0: ', nextflux0
            print 'nextflux1: ', nextflux1

        if fabs(nextlgal0) > gdat.maxmgangmarg or fabs(nextlgal1) > gdat.maxmgangmarg or fabs(nextbgal0) > gdat.maxmgangmarg or fabs(nextbgal1) > gdat.maxmgangmarg or \
                                    nextflux0 < gdat.minmflux or nextflux1 < gdat.minmflux:
            gdatmodi.boolreje = True
                
        if not gdatmodi.boolreje:

            lgal = concatenate((array([nextlgal0, nextlgal1]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgal0, nextbgal1]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]], thisbgal)))
            gdat.listpair = retr_listpair(gdat, lgal, bgal)

            ## first new component
            gdatmodi.drmcsamp[gdat.indxsampchd0+gdat.indxcomplgal, -1] = cdfn_self(nextlgal0, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampchd0+gdat.indxcompbgal, -1] = cdfn_self(nextbgal0, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampchd0+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, nextflux0, gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdat.indxsampchd0+gdat.indxcompsind, -1] = cdfn_eerr(nextsind0, gdat.meansdfn[gdatmodi.indxpoplmodi], gdat.stdvsdfn[gdatmodi.indxpoplmodi], \
                    gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspec0 = retr_spec(gdat, nextflux0, nextsind0)

            ## second new component
            gdatmodi.drmcsamp[gdat.indxsampchd1+gdat.indxcomplgal, -1] = cdfn_self(nextlgal1, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampchd1+gdat.indxcompbgal, -1] = cdfn_self(nextbgal1, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampchd1+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, nextflux1, gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdat.indxsampchd1+gdat.indxcompsind, -1] = cdfn_eerr(nextsind1, gdat.meansdfn[gdatmodi.indxpoplmodi], gdat.stdvsdfn[gdatmodi.indxpoplmodi], \
                    gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspec1 = retr_spec(gdat, nextflux1, nextsind1)

            ## component to be removed
            gdatmodi.modilgal[0] = thislgal
            gdatmodi.modibgal[0] = thisbgal
            gdatmodi.modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            gdatmodi.modilgal[1] = nextlgal0
            gdatmodi.modibgal[1] = nextbgal0
            gdatmodi.modispec[:, 1] = nextspec0.flatten()

            # second component to be added
            gdatmodi.modilgal[2] = nextlgal1
            gdatmodi.modibgal[2] = nextbgal1
            gdatmodi.modispec[:, 2] = nextspec1.flatten()

            tempindxindxpnts = setdiff1d(thisindxindxpnts, spltindxindxpnts)
            lgal = concatenate((gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][tempindxindxpnts]], array([nextlgal0]), array([nextlgal1]))) 
            bgal = concatenate((gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][tempindxindxpnts]], array([nextbgal0]), array([nextbgal1]))) 
            bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]
            listpair = retr_listpair(gdat, lgal, bgal)
            numbpair = len(listpair)
        
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        gdat.numbmodipnts = 3
        
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]])
            
        lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]]
        bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]
        listpair = retr_listpair(gdat, lgal, bgal)
        numbpair = len(listpair)
        
        if gdat.verbtype > 1:
            print 'lgal'
            print lgal
            print 'bgal'
            print bgal
            print 'pairlist'
            print gdat.listpair
           
        if numbpair == 0:
            gdatmodi.boolreje = True
        else:
            indxpairtemp = choice(arange(numbpair))
            mergindxindxpnts0 = gdat.listpair[indxpairtemp][0]
            mergindxindxpnts1 = gdat.listpair[indxpairtemp][1]
  
        if not gdatmodi.boolreje:

            # fisrt PS index to be merged
            mergindxchd0 = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpnts0]
            mergindxsampinit0 = gdat.indxsampcompinit + mergindxchd0 * gdat.numbcomp

            # second PS index to be merged
            gdat.mergindxchd1 = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpnts1]
            mergindxsampinit1 = gdat.indxsampcompinit + gdat.mergindxchd1 * gdat.numbcomp

            # determine the modified sample vector indices
            gdat.indxsampchd0 = gdat.indxsampcompinit + gdat.numbcomp * mergindxchd0
            indxfinl0 = gdat.indxsampchd0 + gdat.numbcomp
            gdat.indxsampchd1 = gdat.indxsampcompinit + gdat.numbcomp * gdat.mergindxchd1
            indxfinl1 = gdat.indxsampchd1 + gdat.numbcomp

            gdatmodi.indxsampmodi = arange(gdat.indxsampchd0, indxfinl0)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxchd0, gdat.mergindxchd1], dtype=int))

            thislgal0 = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpnts0]]
            thisbgal0 = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpnts0]]
            thisflux0 = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpnts0]]
            thissind0 = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpnts0]]

            thislgal1 = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpnts1]]
            thisbgal1 = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpnts1]]
            thisflux1 = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpnts1]]
            thissind1 = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpnts1]]

            # auxiliary component
            gdatmodi.auxipara = zeros(gdat.numbcomp)
            gdatmodi.auxipara[0] = (thislgal0 - thislgal1) / 2.
            gdatmodi.auxipara[1] = (thisbgal0 - thisbgal1) / 2.
            gdatmodi.auxipara[2:] = thisflux0 - thisflux1

            # merged PS
            nextlgal = (thislgal0 + thislgal1) / 2.
            nextbgal = (thisbgal0 + thisbgal1) / 2.
            nextflux = thisflux0 + thisflux1
            
            gdatmodi.drmcsamp[gdat.indxsampchd0, -1] = cdfn_self(nextlgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampchd0+1, -1] = cdfn_self(nextbgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            for i in gdat.indxener:
                gdatmodi.drmcsamp[gdat.indxsampchd0+2+i, -1] = cdfn_flux_powr(gdat, nextflux[i], gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]])

            gdat.numbmodipnts = 3
            ## first component to be merged
            gdatmodi.modilgal[0] = thislgal0
            gdatmodi.modibgal[0] = thisbgal0
            gdatmodi.modispec[:, 0] = -thisspec0.flatten()

            ## first component to be merged
            gdatmodi.modilgal[1] = thislgal1
            gdatmodi.modibgal[1] = thisbgal1
            gdatmodi.modispec[:, 1] = -thisspec1.flatten()

            ## component to be added
            gdatmodi.modilgal[2] = nextlgal
            gdatmodi.modibgal[2] = nextbgal
            gdatmodi.modispec[:, 2] = nextspec.flatten()

            if gdat.verbtype > 1:
                print 'mergindxchd0: ', mergindxchd0
                print 'mergindxindxpnts0: ', mergindxindxpnts0
                print 'mergindxchd1: ', gdat.mergindxchd1
                print 'mergindxindxpnts1: ', mergindxindxpnts1
                print 'indxsampchd0: ', gdat.indxsampchd0
                print 'indxfinl0: ', indxfinl0
                print 'indxsampchd1: ', gdat.indxsampchd1
                print 'indxfinl1: ', indxfinl1
                if gdat.pixltype == 'heal':
                    print 'thislgal0: ', thislgal0
                    print 'thisbgal0: ', thisbgal0
                    print 'thislgal1: ', thislgal1
                    print 'thisbgal1: ', thisbgal1
                else:
                    print 'thislgal0: ', 3600. * thislgal0
                    print 'thisbgal0: ', 3600. * thisbgal0
                    print 'thislgal1: ', 3600. * thislgal1
                    print 'thisbgal1: ', 3600. * thisbgal1 
                print 'thisflux0: ', thisflux0
                print 'thisflux1: ', thisflux1

                if gdat.pixltype == 'heal':
                    print 'nextlgal: ', nextlgal
                    print 'nextbgal: ', nextbgal
                    print 'auxipara[0]: ', gdatmodi.auxipara[0]
                    print 'auxipara[1]: ', gdatmodi.auxipara[1]
                else:
                    print 'nextlgal: ', 3600. * nextlgal
                    print 'nextbgal: ', 3600. * nextbgal
                    print 'auxipara[0]: ', 3600. * gdatmodi.auxipara[0]
                    print 'auxipara[1]: ', 3600. * gdatmodi.auxipara[1]
                print 'nextflux: ', nextflux
                print 'auxipara[2:]: ', gdatmodi.auxipara[2:]
                print

    # component change
    if gdatmodi.thisindxprop >= gdat.indxproplgal:     
        
        gdat.indxenermodi = gdat.indxener
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdatmodi.thisindxprop == gdat.indxproplgal:
                gdat.indxcompmodi = gdat.indxcomplgal
            else:
                gdat.indxcompmodi = gdat.indxcompbgal
        else:
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                gdat.indxcompmodi = gdat.indxcompflux
            elif gdatmodi.thisindxprop == gdat.indxpropsind:
                gdat.indxcompmodi = gdat.indxcompsind
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        gdatmodi.indxsampmodiinit = gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:gdatmodi.indxpoplmodi]) + modiindxpnts * gdat.numbcomp
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdatmodi.indxsampmodiinit + gdat.indxcompmodi
        gdatmodi.indxsampmodispec = gdatmodi.indxsampmodiinit + 2 + gdat.indxener
        
        # propose
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl) 
        if gdatmodi.thisindxprop == gdat.indxpropflux:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvflux)
        if gdatmodi.thisindxprop == gdat.indxpropsind:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvsind)

        thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn, modiindxindxpnts]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
        thisspec = retr_spec(gdat, thisflux, thissind)
            
        gdat.numbmodipnts = 2
        gdatmodi.modispec[:, 0] = -thisspec.flatten()
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdat.indxcompmodi == 0:
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
                if gdat.fdfntype == 'powr':
                    fdfnslop = gdatmodi.thissampvarb[gdat.indxsampfdfnslop[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_powr(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fdfnslop)
                if gdat.fdfntype == 'brok':
                    fdfnbrek = gdatmodi.thissampvarb[gdat.indxsampfdfnbrek[gdatmodi.indxpoplmodi]]
                    fdfnsloplowr = gdatmodi.thissampvarb[gdat.indxsampfdfnsloplowr[gdatmodi.indxpoplmodi]]
                    fdfnslopuppr = gdatmodi.thissampvarb[gdat.indxsampfdfnslopuppr[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_brok(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fdfnbrek, fdfnsloplowr, fdfnslopuppr)
                gdatmodi.modisind[1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            else:
                modiflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfdfn, modiindxindxpnts]]
                gdatmodi.modisind[1] = icdf_eerr(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.meansdfn[gdatmodi.indxpoplmodi], \
                        gdat.stdvsdfn[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            gdatmodi.modispec[:, 1] = retr_spec(gdat, modiflux, gdatmodi.modisind[1]).flatten()

        if gdat.verbtype > 1:
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print 'indxcompmodi: ', gdat.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts

    if gdat.verbtype > 1:
        print 'indxpoplmodi'
        print gdatmodi.indxpoplmodi
        if not (gdatmodi.boolreje or gdatmodi.thisindxprop == gdat.indxpropdeth):
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
            print 'gdatmodi.drmcsamp[gdatmodi.indxsampmodi, :]'
            print gdatmodi.drmcsamp[gdatmodi.indxsampmodi, :]
        
    # energy bin in which to evaluate the log-likelihood
    if gdat.indxpropbrth <= gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdat.indxenermodi = gdat.indxener

    if gdat.verbtype > 1:
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            print 'indxenermodi: ', gdat.indxenermodi

    # auxiliary variable density fraction and jacobian
    if (gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg) and not gdatmodi.boolreje:

        combfact = log(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]**2 / numbpair)
        if gdatmodi.thisindxprop == gdat.indxpropsplt:
            thisjcbnfact = gdat.spltjcbnfact
            thiscombfact = combfact 
        else:
            thisjcbnfact = -gdat.spltjcbnfact
            thiscombfact = -combfact 

        gdatmodi.laccfrac = thisjcbnfact + thiscombfact
        gdatmodi.listnumbpair[gdatmodi.cntrswep] = numbpair
        gdatmodi.listjcbnfact[gdatmodi.cntrswep] = thisjcbnfact
        gdatmodi.listcombfact[gdatmodi.cntrswep] = thiscombfact
        gdatmodi.listauxipara[gdatmodi.cntrswep, :] = gdatmodi.auxipara
        gdatmodi.listlaccfrac[gdatmodi.cntrswep] = gdatmodi.laccfrac

        if gdat.verbtype > 1:
            print 'thisjcbnfact'
            print thisjcbnfact
            print 'thiscombfact'
            print thiscombfact
            print 'laccfrac'
            print gdatmodi.laccfrac
            print 'listpair'
            print listpair
            print

    else:
        gdatmodi.laccfrac = 0.  
   
    # temp
    # define the index to be changed in the sample vector if a hyperparameter is being updated
    if gdat.fdfntype == 'powr':
        if gdatmodi.thisindxprop != gdat.indxpropfdfnslop:
            gdatmodi.indxsampvarbmodi = gdatmodi.indxsampmodi
    if gdat.fdfntype == 'brok':
        if gdatmodi.thisindxprop != gdat.indxpropfdfnbrek and gdatmodi.thisindxprop != gdat.indxpropfdfnsloplowr and gdatmodi.thisindxprop != gdat.indxpropfdfnslopuppr:
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
 
    scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]

    indxpsfiparatemp = numbformpara * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    
    if psfntype == 'singgaus':
        sigc = psfipara[indxpsfiparatemp]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(thisangltemp, sigc)

    elif psfntype == 'singking':
        sigc = psfipara[indxpsfiparatemp]
        gamc = psfipara[indxpsfiparatemp+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(thisangltemp, sigc, gamc)
        
    elif psfntype == 'doubgaus':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(thisangltemp, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        gamt = psfipara[indxpsfiparatemp+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(thisangltemp, frac, sigc, sigt, gamt)
        
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
    psfn /= 2. * pi * trapz(psfn * sin(thisangl[None, :, None]), thisangl, axis=1)[:, None, :]

    return psfn


def retr_psfimodl(gdat):
    
    # PSF parameters
    if gdat.psfntype == 'singgaus':
        gdat.numbformpara = 1
    elif gdat.psfntype == 'singking':
        gdat.numbformpara = 2 
    elif gdat.psfntype == 'doubgaus':
        gdat.numbformpara = 3
    elif gdat.psfntype == 'gausking':
        gdat.numbformpara = 4
    elif gdat.psfntype == 'doubking':
        gdat.numbformpara = 5
       
    gdat.indxformpara = arange(gdat.numbformpara) 
    gdat.numbpsfiparaevtt = gdat.numbener * gdat.numbformpara
    gdat.numbpsfipara = gdat.numbpsfiparaevtt * gdat.numbevtt
    gdat.indxpsfipara = arange(gdat.numbpsfipara)

    minmformpara = zeros(gdat.numbformpara)
    maxmformpara = zeros(gdat.numbformpara)
    factformpara = zeros(gdat.numbformpara)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    if gdat.exprtype == 'ferm':
        minmanglpsfn = 0.1 # deg2rad(0.0001)
        maxmanglpsfn = 5. # deg2rad(5.)
        #minmanglpsfn = 0.01
        #maxmanglpsfn = 3.
        # temp
        minmgamm = 2.
        maxmgamm = 20.
    if gdat.exprtype == 'sdss':
        minmanglpsfn = deg2rad(0.01 / 3600.)
        maxmanglpsfn = deg2rad(3. / 3600.)
    minmpsfnfrac = 0.
    maxmpsfnfrac = 1.
    if gdat.exprtype == 'sdss':
        minmanglpsfn = deg2rad(0.01 / 3600.)
        maxmanglpsfn = deg2rad(2. / 3600.)
    if gdat.psfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif gdat.psfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif gdat.psfntype == 'doubgaus':
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
    elif gdat.psfntype == 'gausking':
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
    elif gdat.psfntype == 'doubking':
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
        if gdat.psfntype == 'singgaus' or gdat.psfntype == 'singking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]] += ' ' + gdat.strganglunit
        elif gdat.psfntype == 'doubgaus' or gdat.psfntype == 'gausking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+2] += ' ' + gdat.strganglunit
        elif gdat.psfntype == 'doubking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+3] += ' ' + gdat.strganglunit
    gdat.indxpsfipara = arange(gdat.numbpsfipara)


def retr_unit(lgal, bgal):

    lgaltemp = deg2rad(lgal)
    bgaltemp = deg2rad(bgal)

    xaxi = cos(bgaltemp) * cos(lgaltemp)
    yaxi = -cos(bgaltemp) * sin(lgaltemp)
    zaxi = sin(bgaltemp)

    return xaxi, yaxi, zaxi


def retr_propmodl(gdat):

    gdat.strgprop = []
    cntr = tdpy.util.cntr()
    
    ## normalization of the FDF
    gdat.indxpropfdfnnorm = cntr.incr()
    gdat.strgprop.append('fdfnnorm')

    ## FDF shape
    if gdat.fdfntype == 'powr':
        gdat.indxpropfdfnslop = cntr.incr()
        gdat.strgprop.append('fdfnslop')
    if gdat.fdfntype == 'brok':
        gdat.indxpropfdfnbrek = cntr.incr()
        gdat.indxpropfdfnsloplowr = cntr.incr()
        gdat.indxpropfdfnslopuppr = cntr.incr()
        gdat.strgprop.append('fdfnbrek')
        gdat.strgprop.append('fdfnsloplowr')
        gdat.strgprop.append('fdfnslopuppr')

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
            
        if gdat.boolpropfdfn:
            probfdfnnorm = array([1.])
            if gdat.fdfntype == 'powr':
                probfdfnslop = array([1.])
            if gdat.fdfntype == 'brok':
                probfdfnbrek = array([1.])
                probfdfnsloplowr = array([1.])
                probfdfnslopuppr = array([1.])
        else:
            probfdfnnorm = array([0.])
            if gdat.fdfntype == 'powr':
                probfdfnslop = array([0.])
            if gdat.fdfntype == 'brok':
                probfdfnbrek = array([0.])
                probfdfnsloplowr = array([0.])
                probfdfnslopuppr = array([0.])

        if gdat.boolproppsfn:
            probpsfipara = array([1.]) * gdat.numbpsfipara
        else:
            probpsfipara = array([0.])
        probnormback = array([1.])
        
        probbrth = array([0.2 * sum(gdat.maxmnumbpnts) / 2.])
        probdeth = array([0.2 * sum(gdat.maxmnumbpnts) / 2.])
        probsplt = array([0. * sum(gdat.maxmnumbpnts) / 2.])
        probmerg = array([0. * sum(gdat.maxmnumbpnts) / 2.])
        
        problgal = array([sum(gdat.maxmnumbpnts) / 2.])
        probbgal = array([sum(gdat.maxmnumbpnts) / 2.])
        probspec = array([sum(gdat.maxmnumbpnts) / 2.])
        probsind = array([sum(gdat.maxmnumbpnts) / 2.])
           
        if gdat.fdfntype == 'powr':
            gdat.probprop = concatenate((probfdfnnorm, probfdfnslop, probpsfipara, probnormback, probbrth, probdeth, \
                probsplt, probmerg, problgal, probbgal, probspec, probsind))
        if gdat.fdfntype == 'brok':
            gdat.probprop = concatenate((probfdfnnorm, probfdfnbrek, probfdfnsloplowr, probfdfnslopuppr, probpsfipara, probnormback, probbrth, probdeth, \
                probsplt, probmerg, problgal, probbgal, probspec, probsind))
        gdat.probprop /= sum(gdat.probprop)
       

def retr_randunitpsfipara(gdat):

    while True:
        randunitpsfipara = rand(gdat.numbpsfipara)
        if gdat.psfntype == 'singgaus' or gdat.psfntype == 'singking':
            break
        else:
            indxpar0 = 1
            if gdat.psfntype == 'doubgaus' or gdat.psfntype == 'gausking':
                indxpar1 = 2
            if gdat.psfntype == 'doubking':
                indxpar1 = 3
            thisbool = True
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    indx = m * gdat.numbpsfiparaevtt + i * gdat.numbformpara
                    thisbool = thisbool and randunitpsfipara[indx+indxpar1] > randunitpsfipara[indx+indxpar0]
            if thisbool:
                break

    return randunitpsfipara


def retr_indxcube(gdat):

    gdat.indxcubefull = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
    gdat.indxcubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixl, gdat.indxevttincl, indexing='ij')
    if gdat.pixltype == 'heal':
        gdat.indxcubeheal = meshgrid(gdat.indxenerinclfull, gdat.indxpixlrofi, gdat.indxevttinclfull, indexing='ij')


def retr_expo(gdat):
 
    if gdat.strgexpo == 'unit':
        gdat.expo = ones((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    else:
        path = gdat.pathdata + '/' + gdat.strgexpo
        gdat.expo = pf.getdata(path)

        if gdat.pixltype == 'heal':
            gdat.expo = gdat.expo[gdat.indxcubeheal]
        else:
            gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))
            
        gdat.expo = gdat.expo[gdat.indxcubefilt]
    

def setp(gdat):
   
    # number of processes
    gdat.strgproc = os.uname()[1]
    if gdat.numbproc == None:
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu':
            gdat.numbproc = 20
        else:
            gdat.numbproc = 1
    gdat.indxproc = arange(gdat.numbproc) 

    if gdat.exprtype == 'ferm':
        if gdat.lablback == None:
            gdat.lablback = [r'$\mathcal{I}$', r'$\mathcal{D}$']
        if gdat.nameback == None:
            gdat.nameback = ['normisot', 'normfdfm']

    gdat.strgfluxunit = retr_strgfluxunit(gdat)
    
    gdat.numbchrototl = 4
    gdat.numbchrollik = 7

    # number of bins
    gdat.numbbins = 10

    gdat.minmnumbpnts = 1

    # the normalized offset for text annotation of point sources in the frames
    gdat.offstext = gdat.maxmgang * 0.05
    
    gdat.numbback = len(gdat.nameback)
    gdat.indxback = arange(gdat.numbback)
    
    gdat.numbevtt = gdat.indxevttincl.size
    gdat.indxevtt = arange(gdat.numbevtt)
    
    # axes
    ## longitude
    gdat.numblgal = 10
    gdat.binslgal = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numblgal + 1)

    ## latitude
    gdat.numbbgal = 10
    gdat.binsbgal = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbbgal + 1)

    ## radial
    gdat.numbgang = 10
    gdat.binsgang = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbgang + 1)

    ## azimuthal
    gdat.numbaang = 10
    gdat.binsaang = linspace(0., 2. * pi, gdat.numbaang + 1)

    ## flux
    gdat.numbflux = 10
    gdat.indxflux = arange(gdat.numbflux)
    gdat.binsflux = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbflux + 1)
    gdat.meanflux = sqrt(gdat.binsflux[1:] * gdat.binsflux[:-1])
    gdat.diffflux = gdat.binsflux[1:] - gdat.binsflux[:-1]
    ### pivot flux bin
    gdat.indxfluxpivt = gdat.numbflux / 2
    gdat.pivtflux = gdat.meanflux[gdat.indxfluxpivt]

    ## color
    gdat.numbsind = 10
    gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbsind + 1)
    gdat.meansind = (gdat.binssind[1:] + gdat.binssind[:-1]) / 2.
    gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    ## energy
    gdat.numbener = gdat.indxenerincl.size
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.binsener = array([0.1, 0.3, 1., 3., 10., 100.])[gdat.indxenerinclbins]
    gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]

    gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
    gdat.indxener = arange(gdat.numbener, dtype=int)
    
    gdat.indxenerfdfn = array([gdat.numbener / 2])
    gdat.enerfdfn = gdat.meanener[gdat.indxenerfdfn]
        
    factener = (gdat.meanener[gdat.indxenerfdfn] / gdat.meanener)**2

    gdat.minmspec = gdat.minmflux * factener
    gdat.maxmspec = gdat.maxmflux * factener
    gdat.binsspec = gdat.binsflux[None, :] * factener[:, None]
    gdat.meanspec = empty((gdat.numbener, gdat.numbflux))
    for i in gdat.indxener:
        gdat.meanspec[i, :] = sqrt(gdat.binsspec[i, 1:] * gdat.binsspec[i, :-1])

    # temp
    if gdat.exprtype == 'sdss':
        gdat.diffener = ones(gdat.numbener)
    
    # half-size of the ROI including the margins
    gdat.maxmgangmarg = gdat.maxmgang + gdat.margsize
        
    # angular gdat.deviation
    gdat.numbangl = 100
    if gdat.exprtype == 'sdss':
        gdat.maxmangl = deg2rad(15. / 3600.) # [rad]
    if gdat.exprtype == 'ferm':
        if gdat.specfraceval == 0.:
            gdat.maxmangl = deg2rad(3. * gdat.maxmgangmarg) # [rad]
        else:
            gdat.maxmangl = deg2rad(20.) # [rad]
    gdat.binsangl = linspace(0., gdat.maxmangl, gdat.numbangl) # [rad]
    
    # plotting
    ## figure size
    gdat.plotsize = 7
    ## text
    if gdat.exprtype == 'sdss':
        gdat.binsanglplot = rad2deg(gdat.binsangl) * 3600.
    if gdat.exprtype == 'ferm':
        gdat.binsanglplot = rad2deg(gdat.binsangl)

    if gdat.trueinfo:
        if gdat.datatype == 'mock':
            gdat.truelabl = 'Mock'
        if gdat.datatype == 'inpt':
            if gdat.exprtype == 'ferm':
                gdat.truelabl = '3FGL'
            if gdat.exprtype == 'sdss':
                gdat.truelabl = 'Hubble'

    if gdat.exprtype == 'ferm':
        gdat.strganglunit = '[deg]'
    if gdat.exprtype == 'sdss':
        gdat.strganglunit = '[arcsec]'

    if gdat.regitype == 'igal':
        gdat.longlabl = '$l$'
        gdat.latilabl = '$b$'
    if gdat.regitype == 'ngal':
        gdat.longlabl = r'$\nu$'
        gdat.latilabl = r'$\mu$'
    if gdat.exprtype == 'ferm':
        gdat.longlabl += r' [$^\circ$]'
        gdat.latilabl += r' [$^\circ$]'
    if gdat.exprtype == 'sdss':
        gdat.longlabl += ' [arcsec]'
        gdat.latilabl += ' [arcsec]'
    
    gdat.truelablvari = gdat.truelabl + ' variable'
    gdat.truelablmiss = gdat.truelabl + ' miss'
    gdat.truelablbias = gdat.truelabl + ' off'
    gdat.truelablhits = gdat.truelabl + ' hit'
    gdat.truelablmult = gdat.truelabl + ' mult'

    # energy bin string
    gdat.enerstrg, gdat.binsenerstrg = retr_enerstrg(gdat)
    
    # PSF class string
    gdat.evttstrg = []
    for m in gdat.indxevtt:
        gdat.evttstrg.append('PSF%d' % gdat.indxevttincl[m])
        
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)

    gdat.spltjcbnfact = log(2.**(2 - gdat.numbener))
    
    # number of components
    # temp -- 3->4 4->5
    gdat.numbcomp = 4 + gdat.numbener
    gdat.numbcompcolr = 5
    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
    
    # component indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompspec = 2 + gdat.indxener
    gdat.indxcompflux = 2 + gdat.indxenerfdfn
    gdat.indxcompsind = 2 + gdat.numbener
    
    # population index vector
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)

    # convenience factors for CDF and ICDF transforms
    gdat.factfdfnnorm = log(gdat.maxmfdfnnorm / gdat.minmfdfnnorm)
    if gdat.fdfntype == 'powr':
        gdat.factfdfnslop = arctan(gdat.maxmfdfnslop) - arctan(gdat.minmfdfnslop)
    if gdat.fdfntype == 'brok':
        gdat.factfdfnbrek = log(gdat.maxmfdfnbrek / gdat.minmfdfnbrek)
        gdat.factfdfnsloplowr = arctan(gdat.maxmfdfnsloplowr) - arctan(gdat.minmfdfnsloplowr)
        gdat.factfdfnslopuppr = arctan(gdat.maxmfdfnslopuppr) - arctan(gdat.minmfdfnslopuppr)
    gdat.factsind = arctan(gdat.maxmsind) - arctan(gdat.minmsind)
    gdat.sindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.meansdfn) / gdat.stdvsdfn / sqrt(2.)) + 1.)
    gdat.sindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.meansdfn) / gdat.stdvsdfn / sqrt(2.)) + 1.)
    gdat.sindcdfnnormdiff = gdat.sindcdfnnormmaxm - gdat.sindcdfnnormminm

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

        gdat.numbmockpsfiparaevtt = gdat.numbener * gdat.numbmockformpara
        gdat.numbmockpsfipara = gdat.numbmockpsfiparaevtt * gdat.numbevtt
        gdat.indxmockpsfipara = arange(gdat.numbmockpsfipara)   

    # maximum number of point sources that can be modified at once
    gdat.numbmodipnts = int(max(3, sum(gdat.maxmnumbpnts)))
    
    # construct the PSF model
    retr_psfimodl(gdat)

    # proposals
    retr_propmodl(gdat)
    
    # factors in the prior expression
    gdat.priofactlgalbgal = 2. * log(1. / 2. / gdat.maxmgang)
    if gdat.fdfntype == 'powr':
        gdat.priofactfdfnslop = gdat.numbener * log(1. / (arctan(gdat.maxmfdfnslop) - arctan(gdat.minmfdfnslop)))
    # temp
    if gdat.fdfntype == 'brok':
        pass
    gdat.priofactfdfnnorm = log(1. / (log(gdat.maxmfdfnnorm) - log(gdat.minmfdfnnorm)))

    # initialize the counter
    cntr = tdpy.util.cntr()
    
    # sample vector indices  
    gdat.indxsampnumbpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfdfnnorm = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    
    # dummy definitions
    gdat.indxsampfdfnslop = -1.
    gdat.indxsampfdfnbrek = -1
    gdat.indxsampfdfnsloplowr = -1
    gdat.indxsampfdfnslopuppr = -1

    if gdat.fdfntype == 'powr':
        gdat.indxsampfdfnslop = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    if gdat.fdfntype == 'brok':
        gdat.indxsampfdfnbrek = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
        gdat.indxsampfdfnsloplowr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
        gdat.indxsampfdfnslopuppr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsamppsfipara = arange(gdat.numbpsfipara) + cntr.incr(gdat.numbpsfipara)
    gdat.indxsampnormback = arange(gdat.numbback * gdat.numbener).reshape((gdat.numbback, gdat.numbener)) + cntr.incr(gdat.numbback * gdat.numbener)

    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.indxsampcompinit = amax(gdat.indxsampnormback) + 1
    
    # maximum number of parameters
    gdat.numbpara = int(gdat.indxsampcompinit + sum(gdat.maxmnumbcomp))
    gdat.indxpara = arange(gdat.numbpara)

    if gdat.numbburn == None:
        gdat.numbburn = min(1000000, gdat.numbswep - 1)
    if gdat.factthin == None:
        gdat.factthin = min(5 * gdat.numbpara, gdat.numbswep - gdat.numbburn)

    # run tag
    gdat.rtag = retr_rtag(gdat, None)
    
    # plots
    if gdat.makeplot:
        if (gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu') and getpass.getuser() == 'tansu':
            pathplotbase = '/n/pan/www/tansu/imag/pcat/'
        else:
            pathplotbase = gdat.pathdata + '/imag/'
        gdat.pathplot = pathplotbase + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
        cmnd = 'mkdir -p ' + gdat.pathplot
        os.system(cmnd)

    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc

    # rescale the positional update scale
    gdat.stdvlbhl /= 2. * gdat.maxmgang

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
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
    

    # input data
    if gdat.datatype == 'inpt':
        
        path = gdat.pathdata + '/' + gdat.strgexpr
        exprflux = pf.getdata(path)

        gdat.indxenerinclfull = arange(exprflux.shape[0])
        gdat.indxevttinclfull = arange(exprflux.shape[2])

        if gdat.pixltype == 'heal':
            gdat.numbpixlheal = exprflux.shape[1]
            gdat.numbsideheal = int(sqrt(gdat.numbpixlheal / 12))
        else:
            gdat.numbsidecart = exprflux.shape[1]
            exprflux = exprflux.reshape((exprflux.shape[0], gdat.numbsidecart**2, exprflux.shape[3]))
            
    else:
        
        if gdat.exprtype == 'ferm':
            gdat.indxenerinclfull = arange(5)
            gdat.indxevttinclfull = arange(4)
         
    if gdat.pixltype == 'heal':
        gdat.numbpixlheal = gdat.numbsideheal**2 * 12
        gdat.apix = 4. * pi / gdat.numbpixlheal
    if gdat.pixltype == 'cart':
        gdat.binslgalcart = linspace(gdat.minmlgal, gdat.maxmlgal, gdat.numbsidecart + 1)
        gdat.binsbgalcart = linspace(gdat.minmbgal, gdat.maxmbgal, gdat.numbsidecart + 1)
        gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
        gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
        gdat.apix = deg2rad(2. * gdat.maxmgang / gdat.numbsidecart)**2
        
    # temp
    gdat.tracsamp = False
    
    # center of the ROI
    if gdat.regitype == 'igal':
        gdat.cntrlghp, gdat.cntrbghp = 0., 0.
    else:
        gdat.cntrlghp, gdat.cntrbghp = 0., 90.
   
    # plot settings
    ## marker opacity
    gdat.mrkralph = 0.7
    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 250
    ## ROI
    gdat.exttrofi = array([gdat.minmlgal, gdat.maxmlgal, gdat.minmbgal, gdat.maxmbgal])
    if gdat.exprtype == 'sdss':
        gdat.exttrofi *= 3600.
        gdat.frambndr = gdat.maxmgang * 3600.
        gdat.frambndrmarg = gdat.maxmgangmarg * 3600.
    else:
        gdat.frambndr = gdat.maxmgang
        gdat.frambndrmarg = gdat.maxmgangmarg

    # FDM normalization prior limits
    gdat.factnormback = log(gdat.maxmnormback / gdat.minmnormback)

    # convenience variables
    gdat.fluxhistmodl = empty(gdat.numbflux)
    
    if gdat.exprtype == 'ferm':
        gdat.numbfluxprox = 3
    if gdat.exprtype == 'sdss':
        gdat.numbfluxprox = 1
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
    if gdat.exprtype == 'ferm':
        gdat.maxmangleval = empty(gdat.numbfluxprox)
        for h in gdat.indxfluxprox:
            if gdat.specfraceval == 0:
                gdat.maxmangleval[h] = 5 * gdat.maxmgang
            else:
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.fermpsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = rad2deg(psfnwdth[gdat.indxmaxmangl])
    if gdat.exprtype == 'sdss':
        gdat.maxmangleval = array([10. / 3600.])

    # pizelization
    if gdat.pixltype == 'heal':
        
        lgalheal, bgalheal, gdat.numbpixlheal, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)

        gdat.indxpixlrofi = where((abs(lgalheal) < gdat.maxmgang) & (abs(bgalheal) < gdat.maxmgang))[0]
        gdat.indxpixlrofimarg = where((abs(lgalheal) < gdat.maxmgangmarg + 300. / gdat.numbsideheal) & \
                (abs(bgalheal) < gdat.maxmgangmarg + 300. / gdat.numbsideheal))[0]
        
        gdat.lgalgrid = lgalheal[gdat.indxpixlrofi]
        gdat.bgalgrid = bgalheal[gdat.indxpixlrofi]
    else:
        isidecart = arange(gdat.numbsidecart)
        temp = meshgrid(isidecart, isidecart, indexing='ij')
        gdat.bgalgrid = gdat.bgalcart[temp[1].flatten()]
        gdat.lgalgrid = gdat.lgalcart[temp[0].flatten()]
        
    # store pixels as unit vectors
    gdat.xaxigrid, gdat.yaxigrid, gdat.zaxigrid = retr_unit(gdat.lgalgrid, gdat.bgalgrid)

    gdat.numbpixl = gdat.lgalgrid.size
    gdat.indxpixl = arange(gdat.numbpixl)
    
    if gdat.pixltype == 'heal':
        path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%03d.p' % (gdat.maxmgang)
        if os.path.isfile(path):
            fobj = open(path, 'rb')
            gdat.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            indxpixltemp = where((fabs(lgalheal) < gdat.maxmgangmarg + 2 * gdat.margsize) & (fabs(bgalheal) < gdat.maxmgangmarg + 2 * gdat.margsize))[0]
            gdat.pixlcnvt = zeros(gdat.numbpixlheal, dtype=int) - 1
            for k in range(indxpixltemp.size):
                dist = retr_angldistunit(gdat, lgalheal[indxpixltemp[k]], bgalheal[indxpixltemp[k]], gdat.indxpixl)
                gdat.pixlcnvt[indxpixltemp[k]] = argmin(dist)

            fobj = open(path, 'wb')
            cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
       
    gdat.indxcubefull = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
    
    gdat.indxcubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixl, gdat.indxevttincl, indexing='ij')
    
    gdat.numbpixlsave = min(1000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)


    if gdat.pixltype == 'heal':
        gdat.indxcubeheal = meshgrid(gdat.indxenerinclfull, gdat.indxpixlrofi, gdat.indxevttinclfull, indexing='ij')
        

    if gdat.datatype == 'inpt' and gdat.pixltype == 'heal':
        exprflux = exprflux[gdat.indxcubeheal]


    if gdat.datatype == 'inpt':
        exprflux = exprflux[gdat.indxcubefilt]
 
    if gdat.datatype == 'mock':
        if gdat.exprtype == 'ferm':
            gdat.mockpsfipara = gdat.fermpsfipara
        if gdat.exprtype == 'sdss':
            gdat.mockpsfipara = gdat.sdsspsfipara

    # exposure
    if gdat.strgexpo == 'unit':
        gdat.expo = ones((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    else:
        path = gdat.pathdata + '/' + gdat.strgexpo
        gdat.expo = pf.getdata(path)

        if gdat.pixltype == 'heal':
            gdat.expo = gdat.expo[gdat.indxcubeheal]
        else:
            gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))
            
        gdat.expo = gdat.expo[gdat.indxcubefilt]
   
    # backgrounds
    gdat.backflux = []
    gdat.backfluxmean = []
    for c in gdat.indxback:
        if gdat.strgback[c] == 'unit':
            backfluxtemp = ones((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        else:
            path = gdat.pathdata + '/' + gdat.strgback[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'heal':
                backfluxtemp = backfluxtemp[gdat.indxcubeheal]
            else:
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
            backfluxtemp = backfluxtemp[gdat.indxcubefilt]
        gdat.backflux.append(backfluxtemp)
        gdat.backfluxmean.append(mean(sum(gdat.backflux[c] * gdat.expo, 2) / sum(gdat.expo, 2), 1))

    gdat.minmcnts = gdat.minmflux * sum(mean(gdat.expo, 1), 1) * gdat.diffener
    gdat.maxmcnts = gdat.maxmflux * sum(mean(gdat.expo, 1), 1) * gdat.diffener
    gdat.binscnts = zeros((gdat.numbener, gdat.numbflux + 1))
    for i in gdat.indxener:
        gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbflux + 1) # [1]
        
    # get the experimental catalog
    if gdat.exprtype != None:
        retr_expr(gdat)

    # get count data
    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = exprflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]
    
    # temp
    if gdat.strgcnfg == 'pcat_ferm_expr_ngal':
        print 'CORRECTING THE EXPOSURE.'
        tempcorr = array([1., 1.2, 1.2])
        gdat.datacnts *= tempcorr[:, None, None]

    ## mock data
    if gdat.datatype == 'mock':

        if gdat.mocknumbpnts == None:
            gdat.mocknumbpnts = empty(gdat.numbpopl)
            for l in gdat.indxpopl:
                gdat.mocknumbpnts[l] = random_integers(gdat.minmnumbpnts, gdat.maxmnumbpnts[l])
            
        gdat.truefdfnnorm = gdat.mocknumbpnts
        
        # if mock FDF is not specified by the user, randomly seed it from the prior
        if gdat.mockfdfntype == 'powr':
            if gdat.mockfdfnslop == None:
                gdat.mockfdfnslop = empty(gdat.numbpopl)
                for l in gdat.indxpopl:
                    gdat.mockfdfnslop[l] = icdf_atan(rand(), gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
        if gdat.mockfdfntype == 'brok':
            if gdat.mockfdfnbrek == None:
                gdat.mockfdfnbrek = empty(gdat.numbpopl)
                for l in gdat.indxpopl:
                    gdat.mockfdfnbrek[l] = icdf_atan(rand(), gdat.minmfdfnbrek[l], gdat.factfdfnbrek[l])
            if gdat.mockfdfnsloplowr == None:
                gdat.mockfdfnsloplowr = empty(gdat.numbpopl)
                for l in gdat.indxpopl:
                    gdat.mockfdfnsloplowr[l] = icdf_atan(rand(), gdat.minmfdfnsloplowr[l], gdat.factfdfnsloplowr[l])
            if gdat.mockfdfnslopuppr == None:
                gdat.mockfdfnslopuppr = empty(gdat.numbpopl)
                for l in gdat.indxpopl:
                    gdat.mockfdfnslopuppr[l] = icdf_atan(rand(), gdat.minmfdfnslopuppr[l], gdat.factfdfnslopuppr[l])

        if gdat.mockfdfntype == 'powr':
            gdat.truefdfnslop = gdat.mockfdfnslop
        if gdat.mockfdfntype == 'brok':
            gdat.truefdfnbrek = gdat.mockfdfnbrek
            gdat.truefdfnsloplowr = gdat.mockfdfnsloplowr
            gdat.truefdfnslopuppr = gdat.mockfdfnslopuppr
    
        if gdat.mockpsfipara == None: 
            gdat.mockpsfntype = psfntpye
            numbmockpsfipara = gdat.numbpsfipara
            gdat.mockpsfipara = empty(numbmockpsfipara)
            for k in arange(numbmockpsfipara):
                gdat.mockpsfipara[k] = icdf_psfipara(gdat, rand(), k)
   
        if gdat.mocknormback == None:
            for c in gdat.indxback:
                gdat.mocknormback[c, :] = icdf_logt(rand(gdat.numbener), gdat.minmnormback[c], gdat.factnormback[c])

        mockcnts = [[] for l in gdat.indxpopl]
        mocklgal = [[] for l in gdat.indxpopl]
        mockbgal = [[] for l in gdat.indxpopl]
        mockspec = [[] for l in gdat.indxpopl]
        mocksind = [[] for l in gdat.indxpopl]
        for l in gdat.indxpopl:
            mocklgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            mockbgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
            mockspec[l] = empty((gdat.numbener, gdat.mocknumbpnts[l]))
            if gdat.mockfdfntype == 'powr':
                mockspec[l][gdat.indxenerfdfn, :] = icdf_flux_powr(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfdfnslop[l])
            if gdat.mockfdfntype == 'brok':
                mockspec[l][gdat.indxenerfdfn, :] = icdf_flux_brok(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfdfnbrek[l], gdat.mockfdfnsloplowr[l], gdat.mockfdfnslopuppr[l])
            mocksind[l] = icdf_eerr(rand(gdat.mocknumbpnts[l]), gdat.meansdfn[l], gdat.stdvsdfn[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
            mockspec[l] = retr_spec(gdat, mockspec[l][gdat.indxenerfdfn[0], :], mocksind[l])
            indxpixltemp = retr_indxpixl(gdat, mockbgal[l], mocklgal[l])
            mockcnts[l] = mockspec[l][:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]
        mockpntsflux = retr_pntsflux(gdat, concatenate(mocklgal), concatenate(mockbgal), concatenate(mockspec, axis=1), gdat.mockpsfipara, gdat.mockpsfntype)
        mocktotlflux = retr_rofi_flux(gdat, gdat.mocknormback, mockpntsflux, gdat.indxcubefull)
        mocktotlcnts = mocktotlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]

        gdat.datacnts = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for i in gdat.indxener:
            for k in range(gdat.numbpixl):
                for m in gdat.indxevtt:
                    gdat.datacnts[i, k, m] = poisson(mocktotlcnts[i, k, m])
              
        if gdat.trueinfo:
            
            gdat.truelgal = []
            gdat.truebgal = []
            gdat.truesind = []
            for l in gdat.indxpopl:
                gdat.truelgal.append(mocklgal[l])
                gdat.truebgal.append(mockbgal[l])
                gdat.truesind.append(mocksind[l])
                    
            gdat.truestrg = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgclss = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgassc = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.indxtruepntstimevari = [array([])] * gdat.numbpopl
                    
            gdat.truenumbpnts = gdat.mocknumbpnts
            if gdat.mockfdfntype == 'brok':
                gdat.mockfdfnsloplowr = gdat.mockfdfnsloplowr
                gdat.mockfdfnslopuppr = gdat.mockfdfnslopuppr
                gdat.mockfdfnbrek = gdat.mockfdfnbrek
            else:
                gdat.mockfdfnslop = gdat.mockfdfnslop
            gdat.truecnts = mockcnts
               
            gdat.truespec = []
            for l in gdat.indxpopl:
                gdat.truespectemp = empty((3, gdat.numbener, gdat.mocknumbpnts[l]))
                gdat.truespectemp[:] = mockspec[l][None, :, :]
                gdat.truespec.append(gdat.truespectemp)
            
            gdat.truepsfipara = gdat.mockpsfipara
            gdat.truepsfntype = gdat.mockpsfntype
            gdat.truenormback = gdat.mocknormback

    ## Real data
    # true data
    if gdat.trueinfo:
        if gdat.datatype == 'inpt':
            gdat.truenumbpnts = None
            gdat.truefdfnnorm = None
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
                
        if gdat.datatype == 'mock':
            gdat.truepsfn = retr_psfn(gdat, gdat.truepsfipara, gdat.indxener, gdat.binsangl, gdat.mockpsfntype)
        else:
            if gdat.exprtype == 'ferm':
                gdat.truepsfn = gdat.fermpsfn
            if gdat.exprtype == 'sdss':
                gdat.truepsfn = gdat.sdsspsfn
                
        truefwhm = 2. * retr_psfnwdth(gdat, gdat.truepsfn, 0.5)
        
        truebackcnts = []
        gdat.truesigm = []
        for l in gdat.indxpopl:
            indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
            truebackcntstemp = zeros((gdat.numbener, gdat.truenumbpnts[l], gdat.numbevtt))
            for c in gdat.indxback:
                truebackcntstemp += gdat.backflux[c][:, indxpixltemp, :] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None] * pi * truefwhm[:, None, :]**2 / 4.
            truebackcnts.append(truebackcntstemp)
            gdat.truesigm.append(gdat.truecnts[l] / sqrt(truebackcntstemp))
      
    # sanity checks
    if amax(abs(gdat.datacnts - gdat.datacnts.astype(int)) / gdat.datacnts) > 1e-3:
        print 'Fractional counts!'
        
    if amin(gdat.datacnts) < 0.:
        print 'Negative counts!'

    gdat.datafluxmean = mean(sum(gdat.datacnts, 2) / sum(gdat.expo, 2), 1) / gdat.apix / gdat.diffener
    gdat.datacntsmean = mean(sum(gdat.datacnts, 2), 1)
    gdat.datacntssatu = ceil((amax(sum(gdat.datacnts, 2), 1) - gdat.datacntsmean) * 0.05 + gdat.datacntsmean)
    gdat.resicntssatu = ceil(gdat.datacntssatu * 0.2)
    
    # auxiliary variables for plots
    if gdat.pixltype == 'heal':
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.datacntscarttemp = tdpy.util.retr_cart(gdat.datacnts[i, :, m], gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                    minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
                if i == 0 and m == 0:
                    gdat.datacntscart = zeros((gdat.datacntscarttemp.shape[0], gdat.datacntscarttemp.shape[1], gdat.numbener, gdat.numbevtt))
                gdat.datacntscart[:, :, i, m] = gdat.datacntscarttemp
        
    else:
        gdat.datacntscart = gdat.datacnts.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
        gdat.datacntscart = swapaxes(swapaxes(gdat.datacntscart, 0, 2), 0, 1)

    for i in gdat.indxener:
        indxdatacntscartsatu = where(gdat.datacntscart[:, :, i, :] > gdat.datacntssatu[i])
        gdat.datacntscart[indxdatacntscartsatu[0], indxdatacntscartsatu[1], i, indxdatacntscartsatu[2]] = gdat.datacntssatu[i]

    # make a look-up table of nearby pixels for each pixel
    if gdat.specfraceval == 0:
        path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%03d_%s.p' % (gdat.maxmgang, gdat.pixltype)
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%03d_%s_%.7g_%.7g_%02d_%04d_%04d.p' % (gdat.maxmgang, gdat.pixltype, \
                     gdat.minmflux, gdat.maxmflux, gdat.numbfluxprox, gdat.indxenerincl[gdat.indxmaxmangl[0]], gdat.indxevttincl[gdat.indxmaxmangl[1]])

    global indxpixlprox
    if os.path.isfile(path):
        print 'Retrieving previously computed pixel look-up table...'
        fobj = open(path, 'rb')
        indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        print 'Computing the look-up table...'
        indxpixlprox = [[] for h in range(gdat.numbfluxprox)]
        for j in gdat.indxpixl:
            dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
            dist[j] = 0.
            for h in range(gdat.numbfluxprox):
                indxpixlproxtemp = where(dist < deg2rad(gdat.maxmangleval[h]))[0]
                indxpixlproxtemp = indxpixlproxtemp[argsort(dist[indxpixlproxtemp])]
                indxpixlprox[h].append(indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()
        
    if gdat.verbtype > 1:
        print 'Memory budget: indxpixlprox'
        totl = 0.
        for h in gdat.indxfluxprox:
            for n in gdat.indxpixl:
                totl += sys.getsizeof(indxpixlprox[h][n]) / 2.**20
        print '%.4g MB' % totl


def init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, strgplot):

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

    axis.axvline(gdat.frambndr, ls='--', alpha=0.3, color='black')
    axis.axvline(-gdat.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(gdat.frambndr, ls='--', alpha=0.3, color='black')
    axis.axhline(-gdat.frambndr, ls='--', alpha=0.3, color='black')

    if indxevttplot == None:
        path = gdat.pathplot + strgplot + '%dA_' % gdat.indxenerincl[indxenerplot] + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    else:
        path = gdat.pathplot + strgplot + '%d%d_' % (gdat.indxenerincl[indxenerplot], gdat.indxevttincl[indxevttplot]) + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    
    return figr, axis, path


def supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot):

    # temp
    if True:
        if gdat.trueinfo:
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label='Sample', marker='+', linewidth=2, color='b')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablhits, marker='*', linewidth=2, color='g')
            # temp
            #if gdat.indxtruepntstimevari[indxpoplplot].size > 0:
            #    axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablvari, marker='*', linewidth=2, color='y')
        axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=2)
        
    # true catalog
    if gdat.trueinfo:
        ## get the true catalog
        mrkrsize = retr_mrkrsize(gdat, gdat.truespec[indxpoplplot][0, gdat.indxenerfdfn, :].flatten())
        lgal = copy(gdat.truelgal[indxpoplplot])
        bgal = copy(gdat.truebgal[indxpoplplot])
        if gdat.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
        numbpnts = int(gdat.truenumbpnts[indxpoplplot])
        
        ## associations
        ### missed
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].miss
        axis.scatter(lgal[indx], bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
        
        ### biased
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].bias[indxenerplot]
        axis.scatter(lgal[indx], bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
        
        ### hit
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].hits[indxenerplot]
        axis.scatter(lgal[indx], bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, label=gdat.truelablhits, marker='*', linewidth=2, color='g')
        
        ## time-variability
        # temp
        #if gdat.indxtruepntstimevari[indxpoplplot].size > 0:
        #    axis.scatter(lgal[gdat.indxtruepntstimevari[indxpoplplot]], bgal[gdat.indxtruepntstimevari[indxpoplplot]], \
        #                                s=mrkrsize[gdat.indxtruepntstimevari[indxpoplplot]], label=gdat.truelablvari, marker='*', linewidth=2, color='y')
                
        ## annotate
        # temp
        gdat.boolanot = False
        if gdat.boolanot:
            for a in range(numbpnts):
                strg = ''
                if gdat.truestrg[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrg[indxpoplplot][a]
                if gdat.truestrgassc[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrgassc[indxpoplplot][a]
                if gdat.truestrgclss[indxpoplplot][a] != None:
                    strg += '%s ' % gdat.truestrgclss[indxpoplplot][a]
                if strg != '':
                    axis.text(gdat.truelgal[indxpoplplot][a], gdat.truebgal[indxpoplplot][a] - gdat.offstext, strg, ha='center', va='center', color='g', fontsize=6)

    # model catalog
    mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[indxpoplplot][gdat.indxenerfdfn, :]])
    lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[indxpoplplot]]
    bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[indxpoplplot]]
    if gdat.exprtype == 'sdss':
        lgal *= 3600.
        bgal *= 3600.
    axis.scatter(lgal, bgal, s=mrkrsize, alpha=gdat.mrkralph, label='Sample', marker='+', linewidth=2, color='b')


def retr_levi(listllik):
    
    minmlistllik = amin(listllik)
    levi = log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    return levi


def retr_info(listllik, levi):
    
    info = mean(listllik) - levi

    return info


def retr_imag(gdat, axis, maps, thisindxener, thisindxevtt, logt=False, cmap='Reds', mean=False, satuuppr=None, satulowr=None, titl=''):

    if logt:
        maps = log10(abs(copy(maps)) + 1e-10)

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
        maps = tdpy.util.retr_cart(maps, \
            indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
            minmlgal=gdat.minmlgal, maxmlgal=gdat.maxmlgal, \
            minmbgal=gdat.minmbgal, maxmbgal=gdat.maxmbgal)
    else:
        maps = maps.reshape((gdat.numbsidecart, gdat.numbsidecart)).T
    
    # saturate the map
    if satulowr != None:
        maps[where(maps < satulowr[thisindxener])] = satulowr[thisindxener]
    if satuuppr != None:
        maps[where(maps > satuuppr[thisindxener])] = satuuppr[thisindxener]
   
    # plot
    imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='none')
    axis.set_title(titl)

    # make a color bar
    if thisindxevtt != None or thisindxener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    return axis, cbar


def retr_jcbn():
    
    lgl0, lgla, bgl0, bgla, flx0, flxa, snd0, snda = sympy.symbols('lgl0 lgla bgl0 bgla flx0 flxa snd0 snda')
    matr = sympy.Matrix([[1, 1 - 1 / (flx0 + flxa), 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [1,    -1 / (flx0 + flxa), 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [0,                     0, 1, 1 - aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,                     0, 1,    -aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,                     0, 0,               0,                     1,            0], \
                         [0,                     0, 0,               0,                    -1,            1]])

    jcbn = matr.det()
    return jcbn


def corr_catl(gdat, thisindxpopl, modllgal, modlbgal, modlspec):

    indxtruepntsassc = gdatstrt()
    indxtruepntsassc.miss = []
    indxtruepntsassc.bias = [[] for i in gdat.indxener]
    indxtruepntsassc.hits = [[] for i in gdat.indxener]
    indxtruepntsassc.mult = []
        
    # get the flux limit that delineates the biased associations and hits 
    fluxbias = empty((2, gdat.numbener, gdat.numbflux + 1))
    for i in gdat.indxener:
        fluxbias[:, i, :] = retr_fluxbias(gdat, gdat.binsspec, i)

    indxmodlpnts = zeros_like(gdat.truelgal[thisindxpopl], dtype=int) - 1
    specassc = zeros((gdat.numbener, gdat.truenumbpnts[thisindxpopl]), dtype=float)
    numbassc = zeros_like(gdat.truelgal[thisindxpopl], dtype=int)
    distassc = zeros_like(gdat.truelgal[thisindxpopl]) + deg2rad(3 * gdat.maxmgang)
    dir1 = array([gdat.truelgal[thisindxpopl], gdat.truebgal[thisindxpopl]])

    for k in range(modllgal.size):
        dir2 = array([modllgal[k], modlbgal[k]])
        dist = angdist(dir1, dir2, lonlat=True)
        thisindxtruepnts = where(dist < deg2rad(0.5))[0]
        
        if False:
            print 'k'
            print k
            print 'dist'
            print dist
            print 'thisindxtruepnts'
            print thisindxtruepnts
            print 
        if thisindxtruepnts.size > 0:
            
            # if there are multiple associated true PS, sort them
            indx = argsort(dist[thisindxtruepnts])
            dist = dist[thisindxtruepnts][indx]
            thisindxtruepnts = thisindxtruepnts[indx]
                
            if False:
                print 'thisindxtruepnts.size'
                print thisindxtruepnts.size
                print 'thisindxtruepnts'
                print thisindxtruepnts
                print 'dist'
                print dist

            # store the index of the model PS
            numbassc[thisindxtruepnts[0]] += 1
            if dist[0] < distassc[thisindxtruepnts[0]]:
                specassc[:, thisindxtruepnts[0]] = modlspec[:, k]
                distassc[thisindxtruepnts[0]] = dist[0]
                indxmodlpnts[thisindxtruepnts[0]] = k

            if False:
                print 'numbassc'
                print numbassc
                print 'specassc'
                print specassc
                print 'distassc'
                print distassc
                print 'indxmodlpnts'
                print indxmodlpnts
                print
    
    if False:
        print 'gdat.truespec[thisindxpopl][0, :, :]'
        print gdat.truespec[thisindxpopl][0, :, :]
        print 'specassc'
        print specassc
        print 'indxmodlpnts'
        print indxmodlpnts
        print 'numbassc'
        print numbassc
        print 'distassc'
        print distassc
        print 

    # check whether the flux of the associated model point source matches well with the flux of the deterministic point source
    for k in range(gdat.truenumbpnts[thisindxpopl]):
        if numbassc[k] == 0:
            indxtruepntsassc.miss.append(k)
        else:
            if numbassc[k] > 1:
                indxtruepntsassc.mult.append(k)
            for i in gdat.indxener:
                fluxbiasthis = interp(gdat.truespec[thisindxpopl][0, i, k], gdat.binsspec[i, :], fluxbias[0, i, :])
                boolbias = specassc[i, k] > fluxbiasthis or specassc[i, k] < gdat.truespec[thisindxpopl][0, i, k]**2 / fluxbiasthis 
                if boolbias:
                    indxtruepntsassc.bias[i].append(k)
                else:
                    indxtruepntsassc.hits[i].append(k)
   
    if False:
        print 'indxtruepntsassc.miss'
        print indxtruepntsassc.miss
        print 'indxtruepntsassc.mult'
        print indxtruepntsassc.mult
        print 'indxtruepntsassc.bias'
        print indxtruepntsassc.bias
        print 'indxtruepntsassc.hits'
        print indxtruepntsassc.hits
        print

    return indxmodlpnts, indxtruepntsassc


def retr_fluxbias(gdat, spec, indxenerthis):

    fluxbias = empty((2, spec.shape[1]))
    deno = log10(gdat.maxmspec[indxenerthis]) - log10(gdat.minmspec[indxenerthis])
    
    factlowr = 10.
    factuppr = 1.1
    
    offs = (log10(gdat.maxmspec[indxenerthis]) * log10(factlowr) - log10(gdat.minmspec[indxenerthis]) * log10(factuppr)) / deno
    
    slop = (log10(factuppr) + log10(gdat.maxmspec[indxenerthis]) - log10(factlowr) - log10(gdat.minmspec[indxenerthis])) / deno
    fluxbias[0, :] = 10**(slop * log10(spec[indxenerthis, :]) + offs)

    slop = (-log10(factuppr) + log10(gdat.maxmspec[indxenerthis]) + log10(factlowr) - log10(gdat.minmspec[indxenerthis])) / deno
    fluxbias[1, :] = 10**(slop * log10(spec[indxenerthis, :]) - offs)

    return fluxbias


