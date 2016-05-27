# common imports
from __init__ import *

class gdatstrt(object):
    
    def __init__(self):
        pass


def retr_psfnwdth(gdat, psfn, frac):
    
    fwhm = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            indxangldispgood = argsort(psfn[i, :, m])
            intpfwhm = max(frac * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, indxangldispgood, m]) and intpfwhm < amax(psfn[i, indxangldispgood, m]):
                fwhm[i, m] = interp1d(psfn[i, indxangldispgood, m], gdat.angldisp[indxangldispgood])(intpfwhm)
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
        indxsamplgaltemp = gdat.indxsampcompinit + gdat.maxmnumbcomp * l +             array(indxpntsfull[l], dtype=int) * gdat.numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], gdat.numbener, 0) +repeat(arange(gdat.numbener), len(indxpntsfull[l])).reshape(gdat.numbener, -1))
        if gdat.colrprio:
            indxsampsind.append(indxsamplgaltemp + 2 + gdat.numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, gdat.numbcomp) + tile(arange(gdat.numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(gdat, lgal, bgal, spec, psfipara, psfntype):
    
    numbpnts = lgal.size
    
    # calculate the distance to all pixels from each point source
    dist = empty((gdat.numbpixl, numbpnts))
    for k in range(numbpnts):
        dist[:, k] = retr_dist(gdat, lgal[k], bgal[k], gdat.lgalgrid, gdat.bgalgrid)

    # evaluate the PSF
    pntsflux = empty((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for k in range(numbpnts):
        psfn = retr_psfn(gdat, psfipara, gdat.indxener, dist[:, k], psfntype)
        pntsflux[k, :, :, :] = spec[:, k, None, None] * psfn

    # sum contributions from all PS
    pntsfluxtemp = sum(pntsflux, 0) 

    return pntsfluxtemp


def retr_rofi_flux(gdat, normback, pntsflux, tempindx):

    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += normback[c, :, None, None] * gdat.backflux[c][tempindx]        
    
    return modlflux


def cdfn_spec_brok(gdat, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, minmspec, maxmspec):

    norm = 1. / ((1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) + \
                    ((maxmspec / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))
    
    indxflux = where(flux >= fluxbrek)[0]
    fluxunit = norm / (1. - fdfnsloplowr) * ((flux / fluxbrek)**(1. - fdfnsloplowr) - (minmspec / fluxbrek)**(1. - fdfnsloplowr))
    fluxunit[indxflux] = norm / (1. - fdfnsloplowr) * (1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr)) - \
                                norm / (1. - fdfnslopuppr) * (1. - (flux[indxflux] / fluxbrek)**(1. - fdfnslopuppr))
      
    return fluxunit


def pdfn_spec_brok(gdat, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, minmspec, maxmspec):

    norm = 1. / ((1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) + \
                    ((maxmspec / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))
    
    indxflux = where(flux >= fluxbrek)[0]
    pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnsloplowr)
    pdfnflux = norm * (flux[indxflux] / fluxbrek)**(1. - fdfnslopuppr)
        
    return pdfnflux


def icdf_spec_brok(gdat, fluxunit, fdfnsloplowr, fdfnslopuppr, fluxbrek, minmspec, maxmspec):
    
    norm = 1. / ((1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) - \
                    (1. - (maxmspec / fluxbrek)**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    fluxunitbrek = norm / (1. - fdfnsloplowr) * (1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr))
    
    indxfluxunit = where(fluxunit >= fluxunitbrek)[0]
    flux = fluxbrek * (fluxunit * (1. - fdfnsloplowr) / norm + (minmspec / fluxbrek)**(1. - fdfnsloplowr))**(1. / (1. - fdfnsloplowr))
    flux[indxfluxunit] = fluxbrek * (1. - (norm / (1. - fdfnsloplowr) * (1. - (minmspec / fluxbrek)**(1. - fdfnsloplowr)) - \
                fluxunit[indxfluxunit]) * (1. - fdfnslopuppr) / norm)**(1. / (1. - fdfnslopuppr))

    return flux


def cdfn_spec(gdat, flux, fdfnslop, minmspec, maxmspec):
        
    fluxunit = (flux**(1. - fdfnslop) - minmspec**(1. - fdfnslop)) / (maxmspec**(1. - fdfnslop) - minmspec**(1. - fdfnslop))
        
    return fluxunit


def icdf_spec(gdat, fluxunit, fdfnslop, minmspec, maxmspec):
    
    flux = (fluxunit * (maxmspec**(1. - fdfnslop) - minmspec**(1. - fdfnslop)) + minmspec**(1. - fdfnslop))**(1. / (1. - fdfnslop))
    
    return flux


def pdfn_spec(gdat, flux, fdfnslop, minmspec, maxmspec):
  
    pdfnflux = (1. - fdfnslop) / (maxmspec**(1. - fdfnslop) - minmspec**(1. - fdfnslop)) * flux**(-fdfnslop)
          
    return pdfnflux


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


def retr_thisindxprop(gdat, samp):
    
    # choose the population to be modified
    gdat.indxpoplmodi = choice(gdat.indxpopl)
    
    numbpnts = gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]]
    if numbpnts == gdat.maxmnumbpnts[gdat.indxpoplmodi]:
        gdat.thisindxprop = choice(gdat.indxprop, p=gdat.probpropmaxm)
    elif numbpnts == gdat.minmnumbpnts:
        gdat.thisindxprop = choice(gdat.indxprop, p=gdat.probpropminm)
    else:
        gdat.thisindxprop = choice(gdat.indxprop, p=gdat.probprop)
        

def retr_indxpixl(gdat, bgal, lgal):

    if gdat.pixltype == 'heal':
        indxpixl = gdat.pixlcnvt[ang2pix(gdat.numbsideheal, deg2rad(90. - bgal), deg2rad(lgal))]
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


def retr_llik(gdat, listchrollik, init=False):

    if init:

        gdat.thisllik = gdat.datacnts * log(gdat.thismodlcnts) - gdat.thismodlcnts
        
    elif gdat.thisindxprop >= gdat.indxproppsfipara:

        # load convenience variables
        timebegn = time.time()
        
        ## save time by not reloading all PS into modi structure for PSF parameter updates
        if gdat.thisindxprop == gdat.indxproppsfipara:
            lgal = gdat.thissampvarb[concatenate(gdat.thisindxsamplgal)]
            bgal = gdat.thissampvarb[concatenate(gdat.thisindxsampbgal)]
            spec = gdat.thissampvarb[concatenate(gdat.thisindxsampspec)[gdat.indxenermodi, :]]
        if gdat.thisindxprop >= gdat.indxpropbrth:
            lgal = gdat.modilgal[0:gdat.numbmodipnts]
            bgal = gdat.modibgal[0:gdat.numbmodipnts]
            spec = gdat.modispec[meshgrid(gdat.indxenermodi, arange(gdat.numbmodipnts), indexing='ij')]
        
        timefinl = time.time()
        listchrollik[gdat.cntrswep, 0] = timefinl - timebegn
        
        # determine pixels over which to evaluate the log-likelihood
        timebegn = time.time()
        
        if gdat.thisindxprop == gdat.indxpropnormback:
            gdat.indxpixlmodi = gdat.indxpixl
        if gdat.thisindxprop >= gdat.indxpropbrth or gdat.thisindxprop == gdat.indxproppsfipara:
            thisindxpixlprox = []
            for k in range(gdat.numbmodipnts):
                
                # temp
                # this may not work for extreme color PS!
                # take the flux at the pivot energy
                if gdat.thisindxprop == gdat.indxproppsfipara:
                    spectemp = gdat.thissampvarb[concatenate(gdat.thisindxsampspec)[gdat.indxenerpivt, k]]
                else:
                    spectemp = gdat.modispec[gdat.indxenerpivt, k]

                # find the flux index
                indxspecproxtemp = amin(where(gdat.binsspecprox - abs(spectemp) > 0.)[0]) - 1
                indxpixltemp = retr_indxpixl(gdat, bgal[k], lgal[k])
                thisindxpixlprox.append(gdat.indxpixlprox[indxspecproxtemp][indxpixltemp])
            gdat.indxpixlmodi = unique(concatenate(thisindxpixlprox))
        
        timefinl = time.time()
        listchrollik[gdat.cntrswep, 1] = timefinl - timebegn

        # construct the mesh grid for likelihood evaluation
        timebegn = time.time()
        
        if gdat.thisindxprop >= gdat.indxproppsfipara:
            gdat.indxcubemodi = meshgrid(gdat.indxenermodi, gdat.indxpixlmodi, gdat.indxevtt, indexing='ij')

        timefinl = time.time()
        listchrollik[gdat.cntrswep, 2] = timefinl - timebegn

        # update the model point source flux map, if needed
        timebegn = time.time()
        
        if gdat.thisindxprop == gdat.indxproppsfipara or gdat.thisindxprop >= gdat.indxpropbrth:

            if gdat.thisindxprop == gdat.indxproppsfipara:
                gdat.nextpntsflux[gdat.indxcubemodi] = 0.
            else:
                gdat.nextpntsflux[gdat.indxcubemodi] = gdat.thispntsflux[gdat.indxcubemodi]
                
            for k in range(gdat.numbmodipnts):
                
                # calculate the distance to the pixels to be updated
                dist = retr_dist(gdat, lgal[k], bgal[k], gdat.lgalgrid[thisindxpixlprox[k]], gdat.bgalgrid[thisindxpixlprox[k]])

                # evaluate the PSF over the set of data cubes to be updated
                temppsfipara = copy(gdat.thissampvarb[gdat.indxsamppsfipara])
                if gdat.thisindxprop == gdat.indxproppsfipara:
                    temppsfipara[gdat.indxpsfiparamodi] = gdat.nextsampvarb[gdat.indxsampmodi]
                psfn = retr_psfn(gdat, temppsfipara, gdat.indxenermodi, dist, gdat.psfntype)
                
                # update the data cubes
                for i in range(gdat.indxenermodi.size):
                    gdat.nextpntsflux[gdat.indxenermodi[i], thisindxpixlprox[k], :] += spec[i, k] * psfn[i, :, :]
            
        timefinl = time.time()
        listchrollik[gdat.cntrswep, 3] = timefinl - timebegn

        # update the total model flux map
        timebegn = time.time()
       
        indxtemp = meshgrid(gdat.indxback, gdat.indxenermodi, indexing='ij')

        if gdat.thisindxprop == gdat.indxpropnormback:
            normback = gdat.nextsampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdat.thispntsflux
        if gdat.thisindxprop == gdat.indxproppsfipara or gdat.thisindxprop >= gdat.indxpropbrth:
            normback = gdat.thissampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdat.nextpntsflux
        gdat.nextmodlflux[gdat.indxcubemodi] = retr_rofi_flux(gdat, normback, pntsflux, gdat.indxcubemodi)

        timefinl = time.time()
        listchrollik[gdat.cntrswep, 4] = timefinl - timebegn

        # calculate the total model count map
        timebegn = time.time()
        
        gdat.nextmodlcnts[gdat.indxcubemodi] = gdat.nextmodlflux[gdat.indxcubemodi] * gdat.expo[gdat.indxcubemodi] * \
            gdat.apix * gdat.diffener[gdat.indxenermodi, None, None] # [1]
        
        timefinl = time.time()
        listchrollik[gdat.cntrswep, 5] = timefinl - timebegn

        # calculate the likelihood
        timebegn = time.time()
        
        gdat.nextllik[gdat.indxcubemodi] = gdat.datacnts[gdat.indxcubemodi] * log(gdat.nextmodlcnts[gdat.indxcubemodi]) - \
            gdat.nextmodlcnts[gdat.indxcubemodi]
            
        timefinl = time.time()
        listchrollik[gdat.cntrswep, 6] = timefinl - timebegn

        if not isfinite(gdat.nextllik[gdat.indxcubemodi]).any():
            warnings.warn('Log-likelihood went NAN!')
            
        gdat.deltllik = sum(gdat.nextllik[gdat.indxcubemodi] - gdat.thisllik[gdat.indxcubemodi])
    else:
        gdat.deltllik = 0.
        
    
def retr_fdfn(gdat, fdfnnorm, fdfnslop, i):
               
    fluxhistmodl = fdfnnorm / gdat.diffspec[i, gdat.indxspecpivt] *         gdat.diffspec[i, :] * (gdat.meanspec[i, :] / gdat.fluxpivt[i])**(-fdfnslop)
          
    return fluxhistmodl


def retr_lpri(gdat, init=False):
        
    if init:
        gdat.thislpri = zeros((gdat.numbpopl, gdat.numbener))
        
        for i in gdat.indxenerfdfn:
            for l in gdat.indxpopl:
                fdfnnorm = gdat.thissampvarb[gdat.indxsampfdfnnorm[l]]
                fdfnslop = gdat.thissampvarb[gdat.indxsampfdfnslop[l, i]]
                fluxhistmodl = retr_fdfn(gdat, fdfnnorm, fdfnslop, i)
                spec = gdat.thissampvarb[gdat.thisindxsampspec[l][i, :]]
                fluxhist = histogram(spec, gdat.binsspec[i, :])[0]
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                gdat.thislpri[l, i] = sum(lprbpois)            
      
        gdat.nextlpri = copy(gdat.thislpri)
                
    else:
        gdat.nextlpri = copy(gdat.thislpri)
        if gdat.thisindxprop == gdat.indxpropfdfnnorm or gdat.thisindxprop == gdat.indxpropfdfnslop \
            or gdat.thisindxprop >= gdat.indxpropbrth and gdat.thisindxprop <= gdat.indxpropmerg:
              
            if gdat.colrprio:
                indxenertemp = gdat.indxenerfdfn
            else:
                indxenertemp = gdat.indxenermodi
            for i in indxenertemp:
                if gdat.thisindxprop == gdat.indxpropfdfnnorm:
                    fdfnnorm = gdat.nextsampvarb[gdat.indxsampfdfnnorm[gdat.indxpoplmodi]]
                    fdfnslop = gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]]
                elif gdat.thisindxprop == gdat.indxpropfdfnslop:
                    fdfnnorm = gdat.thissampvarb[gdat.indxsampfdfnnorm[gdat.indxpoplmodi]]
                    fdfnslop = gdat.nextsampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]]
                else:
                    fdfnnorm = gdat.thissampvarb[gdat.indxsampfdfnnorm[gdat.indxpoplmodi]]
                    fdfnslop = gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]]
                fluxhistmodl = retr_fdfn(gdat, fdfnnorm, fdfnslop, i)
    
                fluxhist = histogram(gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][i, :]], gdat.binsspec[i, :])[0] 
                if gdat.thisindxprop == gdat.indxpropbrth:
                    fluxhist += histogram(gdat.modispec[i, 0], gdat.binsspec[i, :])[0]
                elif gdat.thisindxprop == gdat.indxpropdeth:
                    fluxhist -= histogram(-gdat.modispec[i, 0], gdat.binsspec[i, :])[0]
                elif gdat.thisindxprop == gdat.indxpropsplt:
                    fluxhist -= histogram(-gdat.modispec[i, 0], gdat.binsspec[i, :])[0]
                    fluxhist += histogram(gdat.modispec[i, 1:3], gdat.binsspec[i, :])[0]
                elif gdat.thisindxprop == gdat.indxpropmerg:
                    fluxhist -= histogram(-gdat.modispec[i, 0:2], gdat.binsspec[i, :])[0]
                    fluxhist += histogram(gdat.modispec[i, 2], gdat.binsspec[i, :])[0]
                
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                gdat.nextlpri[gdat.indxpoplmodi, i] = sum(lprbpois)

            gdat.deltlpri = sum(gdat.nextlpri[gdat.indxpoplmodi, gdat.indxenermodi] - gdat.thislpri[gdat.indxpoplmodi, gdat.indxenermodi])
        else:
            gdat.deltlpri = 0.
        
        
def retr_sampvarb(gdat, indxpntsfull, samp):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxsampnumbpnts] = samp[gdat.indxsampnumbpnts]
    sampvarb[gdat.indxsampfdfnnorm] = icdf_logt(samp[gdat.indxsampfdfnnorm], gdat.minmfdfnnorm, gdat.factfdfnnorm)
    for l in gdat.indxpopl:
        if gdat.fdfntype == 'brok':
            sampvarb[gdat.indxsampfdfnsloplowr[l]] = icdf_atan(samp[gdat.indxsampfdfnsloplowr[l]], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
            sampvarb[gdat.indxsampfdfnslopuppr[l]] = icdf_atan(samp[gdat.indxsampfdfnslopuppr[l]], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
            sampvarb[gdat.indxsampfdfnbrek[l]] = icdf_logt(samp[gdat.indxsampfdfnbrek[l]], gdat.minmfdfnbrek[l, i], gdat.factfdfnbrek[l, i])
        else:
            sampvarb[gdat.indxsampfdfnslop[l]] = icdf_atan(samp[gdat.indxsampfdfnslop[l]], gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
    for c in gdat.indxback:
        sampvarb[gdat.indxsampnormback[c, :]] = icdf_logt(samp[gdat.indxsampnormback[c, :]], gdat.minmnormback[c], gdat.factnormback[c])
    for k in gdat.indxpsfipara:
        sampvarb[gdat.indxsamppsfipara[k]] = icdf_psfipara(gdat, samp[gdat.indxsamppsfipara[k]], k)
    
    for l in gdat.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
        for i in gdat.indxenerfdfn:
            sampvarb[indxsampspec[l][i, :]] = icdf_spec(gdat, samp[indxsampspec[l][i, :]], sampvarb[gdat.indxsampfdfnslop[l, i]], gdat.minmspec[i], gdat.maxmspec[i])
        if gdat.colrprio:
            sampvarb[indxsampsind[l]] = icdf_eerr(samp[indxsampsind[l]], gdat.meansind[l], gdat.stdvsind[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
            sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampspec[l][gdat.indxenerfdfn, :]], sampvarb[indxsampsind[l]])
    
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


def retr_mrkrsize(gdat, spec, indxenerplot):
    if gdat.colrprio:
        indxenertemp = gdat.indxenerfdfn[0]
    else:
        indxenertemp = indxenerplot

    mrkrsize = (spec[indxenertemp, :] - gdat.minmspec[indxenertemp]) / (gdat.maxmspec[indxenertemp] - gdat.minmspec[indxenertemp]) * \
		(gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
        
    return mrkrsize


def retr_scalangl(gdat, angl):
    
    temp = sqrt(2. - 2. * cos(angl[None, :, None]))
    scalangl = 2. * arcsin(temp / 2.) / gdat.fermscalfact[:, None, :]
    
    return scalangl


def retr_fermpsfn(gdat):
   
    name = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    irfn = pf.getdata(name, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']

    gdat.numbfermscalpara = 3
    gdat.numbfermformpara = 5
    
    fermscal = zeros((gdat.numbevtt, gdat.numbfermscalpara))
    fermform = zeros((gdat.numbener, gdat.numbevtt, gdat.numbfermformpara))
    gdat.fermpsfipara = zeros((gdat.numbener * gdat.numbfermformpara * gdat.numbevtt))
    
    for m in gdat.indxevtt:
        fermscal[m, :] = pf.getdata(name, 2 + 3 * gdat.indxevttincl[m])['PSFSCALE']
        irfn = pf.getdata(name, 1 + 3 * gdat.indxevttincl[m])
        for k in range(gdat.numbfermformpara):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)

    factener = (10. * gdat.meanener[:, None])**fermscal[None, :, 2]
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * factener)**2 + fermscal[None, :, 1]**2)
    
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)

    #scalangl = retr_scalangl(gdat, gdat.angldisp)
    
    fermform[:, :, 1] = gdat.fermscalfact * fermform[:, :, 1]
    fermform[:, :, 3] = gdat.fermscalfact * fermform[:, :, 3]

    # store the fermi PSF parameters
    for m in gdat.indxevtt:
        for k in range(gdat.numbfermformpara):
            indxfermpsfiparatemp = m * gdat.numbfermformpara * gdat.numbener + gdat.indxener * gdat.numbfermformpara + k
            gdat.fermpsfipara[indxfermpsfiparatemp] = fermform[:, m, k]
        
    frac = fermform[:, :, 0]
    sigc = fermform[:, :, 1]
    gamc = fermform[:, :, 2]
    sigt = fermform[:, :, 3]
    gamt = fermform[:, :, 4]
   
    # temp
    angl = gdat.angldisp[None, :, None]
    #angl = scalangl
    gdat.fermpsfn = retr_doubking(angl, frac[:, None, :], sigc[:, None, :], gamc[:, None, :], sigt[:, None, :], gamt[:, None, :])


def retr_sdsspsfn(gdat):
   
    numbpsfiparaevtt = gdat.numbener * 3
    frac = array([[0.5, 0.5, 0.5]]).T
    sigc = array([[0.5, 1., 3.]]).T
    sigt = array([[0.5, 1., 3.]]).T
    gdat.sdsspsfipara = empty(numbpsfiparaevtt)
    gdat.sdsspsfn = retr_doubgaus(gdat.angldisp[None, :, None], frac[:, None, :], sigc[:, None, :], sigt[:, None, :])


def updt_samp(gdat):
    
    if gdat.thisindxprop == gdat.indxpropfdfnnorm:
        gdat.thissampvarb[gdat.indxsampfdfnnorm[gdat.indxpoplmodi]] = gdat.nextsampvarb[gdat.indxsampfdfnnorm[gdat.indxpoplmodi]]
        gdat.thislpri[gdat.indxpoplmodi, gdat.indxenermodi] = gdat.nextlpri[gdat.indxpoplmodi, gdat.indxenermodi]

    if gdat.thisindxprop == gdat.indxpropfdfnslop:
 
        # update the unit sample vector
        if gdat.colrprio:
            indxenertemp = gdat.indxenerfdfn
        else:
            indxenertemp = gdat.indxenermodi
        
        if False:
            print 'gdat.thissampvarb[gdat.thisindxsampspec]'
            for k in range(gdat.thisindxsampspec[0].shape[1]):
                print gdat.thissampvarb[gdat.thisindxsampspec[0][:, k]]
            print 'gdat.drmcsamp[gdat.thisindxsampspec]'
            for k in range(gdat.thisindxsampspec[0].shape[1]):
                print gdat.drmcsamp[gdat.thisindxsampspec[0][:, k]]
        
        gdat.drmcsamp[gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenermodi, :], -1] = \
            cdfn_spec(gdat, gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenermodi, :]], \
            gdat.nextsampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]], \
            gdat.minmspec[gdat.indxenermodi], gdat.maxmspec[gdat.indxenermodi])
       
        if False:
            print 'gdat.thissampvarb[gdat.thisindxsampspec]'
            for k in range(gdat.thisindxsampspec[0].shape[1]):
                print gdat.thissampvarb[gdat.thisindxsampspec[0][:, k]]
            print 'gdat.drmcsamp[gdat.thisindxsampspec]'
            for k in range(gdat.thisindxsampspec[0].shape[1]):
                print gdat.drmcsamp[gdat.thisindxsampspec[0][:, k]]
            print
            print

        # update the sample vector
        gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]] = \
            gdat.nextsampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]]
            
        # update the prior register
        gdat.thislpri[gdat.indxpoplmodi, gdat.indxenermodi] = gdat.nextlpri[gdat.indxpoplmodi, gdat.indxenermodi]

    # likelihood updates
    if gdat.thisindxprop >= gdat.indxproppsfipara:
        gdat.thisllik[gdat.indxcubemodi] = gdat.nextllik[gdat.indxcubemodi]
        gdat.thismodlcnts[gdat.indxcubemodi] = gdat.nextmodlcnts[gdat.indxcubemodi]
        
    if gdat.thisindxprop == gdat.indxproppsfipara:
        gdat.thissampvarb[gdat.indxsampmodi] = gdat.nextsampvarb[gdat.indxsampmodi]
        
    if gdat.thisindxprop == gdat.indxpropnormback:
        gdat.thissampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]] = \
            gdat.nextsampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]]
        
    if gdat.thisindxprop >= gdat.indxpropbrth or gdat.thisindxprop == gdat.indxproppsfipara:

        # temp
        figr, axis, path = init_fram(gdat, 0, 0, 'thispntsflux')
        axis, cbar = retr_imag(gdat, axis, gdat.thispntsflux, 0, 0, logt=False, titl='thispntsflux')
        plt.savefig(path)
        plt.close(figr)
        
        modipixl = zeros_like(gdat.thispntsflux)
        modipixl[0, gdat.indxpixlmodi, 0] = 1.
        figr, axis, path = init_fram(gdat, 0, 0, 'pixlmodi')
        axis, cbar = retr_imag(gdat, axis, modipixl, 0, 0, logt=False, titl='pixlmodi')
        plt.savefig(path)
        plt.close(figr)

        gdat.thispntsflux[gdat.indxcubemodi] = gdat.nextpntsflux[gdat.indxcubemodi]
        
        # temp
        figr, axis, path = init_fram(gdat, 0, 0, 'thispntsfluxprim')
        axis, cbar = retr_imag(gdat, axis, gdat.thispntsflux, 0, 0, logt=False, titl='thispntsfluxprim')
        plt.savefig(path)
        plt.close(figr)

    # transdimensinal updates
    if gdat.thisindxprop >= gdat.indxpropbrth and gdat.thisindxprop <= gdat.indxpropmerg:
        gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] = gdat.nextsampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]]
        gdat.thislpri[gdat.indxpoplmodi, gdat.indxenermodi] = gdat.nextlpri[gdat.indxpoplmodi, gdat.indxenermodi]
        
    # birth
    if gdat.thisindxprop == gdat.indxpropbrth:
        
        # update the PS index lists
        gdat.thisindxpntsfull[gdat.indxpoplmodi].append(gdat.thisindxpntsempt[gdat.indxpoplmodi][0])
        del gdat.thisindxpntsempt[gdat.indxpoplmodi][0]

        # update the components
        gdat.thissampvarb[gdat.indxsampmodi[0]] = gdat.modilgal[0]
        gdat.thissampvarb[gdat.indxsampmodi[1]] = gdat.modibgal[0]
        if gdat.colrprio:
            gdat.thissampvarb[gdat.indxsampmodi[2+gdat.indxener]] = gdat.modispec[:, 0]
            gdat.thissampvarb[gdat.indxsampmodi[2+gdat.numbener]] = gdat.modisind[0]
        else:
            gdat.thissampvarb[gdat.indxsampmodi[2:]] = gdat.modispec[gdat.indxenermodi, 0]
        
    # death
    if gdat.thisindxprop == gdat.indxpropdeth:
        
        # update the PS index lists
        gdat.thisindxpntsempt[gdat.indxpoplmodi].append(gdat.killindxpnts)
        gdat.thisindxpntsfull[gdat.indxpoplmodi].remove(gdat.killindxpnts)


    # split
    if gdat.thisindxprop == gdat.indxpropsplt:

        # update the PS index lists
        gdat.thisindxpntsfull[gdat.indxpoplmodi].append(gdat.thisindxpntsempt[gdat.indxpoplmodi][0])
        del gdat.thisindxpntsempt[gdat.indxpoplmodi][0]
        
        # update the components
        # first component
        gdat.thissampvarb[gdat.indxsampchd0] = gdat.modilgal[1]
        gdat.thissampvarb[gdat.indxsampchd0+1] = gdat.modibgal[1]
        gdat.thissampvarb[gdat.indxsampchd0+2:gdat.indxsampchd0+2+gdat.numbener] = gdat.modispec[:, 1]
  
        # second component
        gdat.thissampvarb[gdat.indxsampchd1] = gdat.modilgal[2]
        gdat.thissampvarb[gdat.indxsampchd1+1] = gdat.modibgal[2]
        gdat.thissampvarb[gdat.indxsampchd1+2:gdat.indxsampchd1+2+gdat.numbener] = gdat.modispec[:, 2]
        
    # merge
    if gdat.thisindxprop == gdat.indxpropmerg:
        
        # update the PS index lists
        gdat.thisindxpntsfull[gdat.indxpoplmodi].remove(gdat.mergindxchd1)
        gdat.thisindxpntsempt[gdat.indxpoplmodi].append(gdat.mergindxchd1)

        # update the component
        gdat.thissampvarb[gdat.indxsampmodi[0]] = gdat.modilgal[2]
        gdat.thissampvarb[gdat.indxsampmodi[1]] = gdat.modibgal[2]
        if gdat.colrprio:
            gdat.thissampvarb[gdat.indxsampmodi[2:2+gdat.indxener]] = gdat.modispec[:, 2]
            gdat.thissampvarb[gdat.indxsampmodi[2+gdat.numbener]] = gdat.modisind[2]
        else:
            gdat.thissampvarb[gdat.indxsampmodi[2:]] = gdat.modispec[:, 2]
        
    # component change
    if gdat.thisindxprop >= gdat.indxproplgal:  
        if gdat.thisindxprop == gdat.indxproplgal:
            gdat.thissampvarb[gdat.indxsampmodi] = gdat.modilgal[1]
        elif gdat.thisindxprop == gdat.indxpropbgal:
            gdat.thissampvarb[gdat.indxsampmodi] = gdat.modibgal[1]
        else:
            if gdat.colrprio:
                gdat.thissampvarb[gdat.indxsampmodispec] = gdat.modispec[:, 1]
                if gdat.thisindxprop == gdat.indxpropsind:
                    gdat.thissampvarb[gdat.indxsampmodi] = gdat.modisind[1]
            else:
                gdat.thissampvarb[gdat.indxsampmodi] = gdat.modispec[gdat.indxenermodi, 1]

    # update the unit sample vector
    if gdat.indxsampmodi.size > 0:
        gdat.drmcsamp[gdat.indxsampmodi, 0] = gdat.drmcsamp[gdat.indxsampmodi, 1]


def retr_postvarb(listvarb):

    shap = zeros(len(listvarb.shape), dtype=int)
    shap[0] = 3
    shap[1:] = listvarb.shape[1:]
    shap = list(shap)
    postvarb = zeros(shap)
    
    postvarb[0, :] = percentile(listvarb, 50., axis=0)
    postvarb[1, :] = percentile(listvarb, 16., axis=0)
    postvarb[2, :] = percentile(listvarb, 84., axis=0)

    return postvarb


def retr_errrvarb(postvarb):

    errr = abs(postvarb[0, :] - postvarb[1:3, :])

    return errr


def retr_pairlist(lgal, bgal):
    
    pairlist = []
    for k in range(lgal.size):
        indxpnts = k + 1 + where((lgal[k+1:] < lgal[k] + spmrlbhl) &     (lgal[k+1:] > lgal[k] - spmrlbhl) &     (bgal[k+1:] < bgal[k] + spmrlbhl) &     (bgal[k+1:] > bgal[k] - spmrlbhl))[0]
        for l in range(indxpnts.size):
            pairlist.append([k, indxpnts[l]])
            
    return pairlist


def retr_fgl3(gdat):
        
    path = os.environ["PCAT_DATA_PATH"] + '/gll_psc_v16.fit'

    fgl3 = pf.getdata(path)
    
    gdat.fgl3numbpnts = fgl3['glon'].size
    
    gdat.fgl3lgal = fgl3['glon']
    gdat.fgl3lgal = ((gdat.fgl3lgal - 180.) % 360.) - 180.

    gdat.fgl3bgal = fgl3['glat']

    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    gdat.fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    gdat.fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    gdat.fgl3sind = fgl3['Spectral_Index']
    
    gdat.fgl3spectype = fgl3['SpectrumType']
    gdat.fgl3scur = fgl3['beta']
    gdat.fgl3scut = fgl3['Cutoff'] * 1e-3
    
    gdat.fgl3timevari = fgl3['Variability_Index']
    gdat.indxfgl3timevari = where(gdat.fgl3timevari > 72.44)[0]
    
    fgl3spectemp = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], fgl3['Flux10000_100000']))
    fgl3spectemp = fgl3spectemp[gdat.indxenerincl, :] / gdat.diffener[:, None]
    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], fgl3['Unc_Flux10000_100000']))
    fgl3specstdvtemp = fgl3specstdvtemp[gdat.indxenerincl, :, :] / gdat.diffener[:, None, None]
    
    gdat.fgl3spec = zeros((3, gdat.numbener, gdat.fgl3numbpnts))
    gdat.fgl3spec[0, :, :] = fgl3spectemp
    gdat.fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdvtemp[:, :, 0]
    gdat.fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdvtemp[:, :, 1]
    
    indxpixlfgl3 = retr_indxpixl(gdat, gdat.fgl3bgal, gdat.fgl3lgal)
    gdat.fgl3cnts = gdat.fgl3spec[0, :, :, None] * gdat.expo[:, indxpixlfgl3, :] * gdat.diffener[:, None, None]
    gdat.fgl3gang = rad2deg(arccos(cos(deg2rad(gdat.fgl3lgal)) * cos(deg2rad(gdat.fgl3bgal))))
        
    # adjust 3FGL positions according to the ROI center
    if gdat.regitype == 'ngal':
        rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
        gdat.fgl3bgal, gdat.fgl3lgal = rad2deg(rttr(deg2rad(90. - gdat.fgl3bgal), deg2rad(gdat.fgl3lgal)))
        gdat.fgl3bgal = 90. - gdat.fgl3bgal

    # select the 3FGL point sources in the ROI
    gdat.indxfgl3rofi = arange(gdat.fgl3lgal.size, dtype=int)
    for i in gdat.indxener:
        gdat.indxfgl3rofi = intersect1d(where((gdat.fgl3spec[0, i, :] > gdat.minmspec[i]) & \
            (gdat.fgl3spec[0, i, :] < gdat.maxmspec[i]))[0], gdat.indxfgl3rofi)
    gdat.indxfgl3rofi = intersect1d(where((abs(gdat.fgl3lgal) < gdat.maxmgangmarg) & \
            (abs(gdat.fgl3bgal) < gdat.maxmgangmarg))[0], gdat.indxfgl3rofi)
    gdat.fgl3numbpntsrofi = gdat.indxfgl3rofi.size
    gdat.indxfgl3timevarirofi = where(gdat.fgl3timevari[gdat.indxfgl3rofi] > 72.44)[0]

    # sanity check
    for i in gdat.indxener:
        if (gdat.fgl3spec[0, i, gdat.indxfgl3rofi] > gdat.maxmspec[i]).any():
            print 'maxmspec %d is bad!' % i


def retr_rtag(gdat, indxprocwork):
    
    if indxprocwork == None:
        rtag = 'AA_%d_%d_%d_%d_%s_%s_%s' % (gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
    else:
        rtag = '%02d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
        
    return rtag


def retr_gaus(gdat, indxsamp, stdv):
    
    if gdat.fracrand > 0.:
        if rand() < gdat.fracrand:
            gdat.drmcsamp[indxsamp, 1] = rand()
        else:
            gdat.drmcsamp[indxsamp, 1] = gdat.drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        gdat.drmcsamp[indxsamp, 1] = gdat.drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_dist(gdat, lgal0, bgal0, lgal1, bgal1):
    
    if gdat.pixltype == 'heal':
        dir1 = array([lgal0, bgal0])
        dir2 = array([lgal1, bgal1])
        dist = angdist(dir1, dir2, lonlat=True)
    else:
        dist = deg2rad(sqrt((lgal0 - lgal1)**2 + (bgal0 - bgal1)**2))

    return dist
       
    
def retr_strgprop(gdat):
    
    gdat.strgprop = ['FDF Norm', 'FDF Slope', 'PSF', 'Back. Norm', 'Birth', 'Death', 'Split', 'Merge', \
  'Longitude Update', 'Latitude Update', 'Flux Update', 'Spectral Index Update']


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


def chsq_fdfnslop(gdat, para, i):

    fdfnslop = para[0]
    fdfnnormnorm = para[1]
    
    fluxhistmodl = fdfnnormnorm * gdat.diffspec[i, :] * pdfn_spec(gdat.meanspec[i, :], fdfnslop, gdat.minmspec[i], gdat.maxmspec[i])

    chsq = sum(((fluxhistmodl.flatten()[jspecfgl3] - gdat.fgl3spechist[i, jspecfgl3]) / gdat.fgl3spechist[i, jspecfgl3])**2)
    
    return chsq


def retr_strgfluxunit(gdat):
    
    if gdat.exprtype == 'ferm':
        strgfluxunit = r'[1/cm$^2$/s/GeV]'
    if gdat.exprtype == 'sdss':
        strgfluxunit = '[nMgy]'
        
    return strgfluxunit
     
    
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


def retr_prop(gdat):
  
    gdat.thisindxsamplgal, gdat.thisindxsampbgal,  gdat.thisindxsampspec, gdat.thisindxsampsind, \
            gdat.thisindxsampcomp = retr_indx(gdat, gdat.thisindxpntsfull)
    
    if gdat.verbtype > 2:
        print 'retr_prop(): '

        print 'drmcsamp'
        print gdat.drmcsamp
        
        print 'thissampvarb: '
        for k in range(gdat.thissampvarb.size):
            if k == gdat.indxsampcompinit:
                print
            if k > gdat.indxsampcompinit and (k - gdat.indxsampcompinit) % gdat.numbcomp == 0:
                print
            print gdat.thissampvarb[k]
        print
            
        print 'thisindxpntsfull: ', gdat.thisindxpntsfull
        print 'thisindxpntsempt: ', gdat.thisindxpntsempt  
        print 'thisindxsamplgal: ', gdat.thisindxsamplgal
        print 'thisindxsampbgal: ', gdat.thisindxsampbgal
        print 'thisindxsampspec: '
        print gdat.thisindxsampspec
        if gdat.colrprio:
            print 'thisindxsampsind: ', gdat.thisindxsampsind
        print 'thisindxsampcomp: ', gdat.thisindxsampcomp
        print
        
    # hyper-parameter changes
    # mean number of point sources
    if gdat.thisindxprop == gdat.indxpropfdfnnorm:
        gdat.indxsampmodi = gdat.indxsampfdfnnorm[gdat.indxpoplmodi]
        retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvfdfnnorm)
        gdat.nextsampvarb[gdat.indxsampfdfnnorm] = icdf_logt(gdat.drmcsamp[gdat.indxsampmodi, -1], gdat.minmfdfnnorm, gdat.factfdfnnorm)
        if gdat.colrprio:
            gdat.indxenermodi = gdat.indxenerfdfn
        else:
            gdat.indxenermodi = gdat.indxener
        
    # flux distribution function slope
    if gdat.thisindxprop == gdat.indxpropfdfnslop:
        if gdat.colrprio:
            gdat.indxenermodi = gdat.indxenerfdfn
        else:
            gdat.indxenermodi = choice(gdat.indxener)
        gdat.indxsampfdfnslopmodi = gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]
        retr_gaus(gdat, gdat.indxsampfdfnslopmodi, gdat.stdvfdfnslop)
        
        gdat.nextsampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]] = \
            icdf_atan(gdat.drmcsamp[gdat.indxsampfdfnslopmodi, -1], gdat.minmfdfnslop, gdat.factfdfnslop)
        if gdat.colrprio:
            gdat.indxsampmodi = concatenate((gdat.indxsampfdfnslopmodi, gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxener, :].flatten()))
        else:
            gdat.indxsampmodi = concatenate((array([gdat.indxsampfdfnslopmodi]), gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenermodi, :]))
           
        if gdat.verbtype > 2:
            print 'indxenermodi'
            print gdat.indxenermodi
            print 'indxpoplmodi'
            print gdat.indxpoplmodi
            print 'indxsampmodi'
            print gdat.indxsampmodi
            print 'drmcsamp[indxsampmodi, :]'
            print gdat.drmcsamp[gdat.indxsampmodi, :]
            print 'thissampvarb[gdat.indxsampfdfnslop]'
            print gdat.thissampvarb[gdat.indxsampfdfnslop]
            print 'nextsampvarb[gdat.indxsampfdfnslop]'
            print gdat.nextsampvarb[gdat.indxsampfdfnslop]
            print

    # PSF parameter change 
    if gdat.thisindxprop == gdat.indxproppsfipara:
        
        # index of the PSF parameter to change
        gdat.indxpsfiparamodi = choice(gdat.indxpsfipara)

        # the energy bin of the PS flux map to be modified
        gdat.indxenermodi = array([(gdat.indxpsfiparamodi % gdat.numbpsfiparaevtt) // gdat.numbformpara])
        
        # sample index to be modified
        gdat.indxsampmodi = gdat.indxsamppsfipara[gdat.indxpsfiparamodi]
        retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvpsfipara)
        gdat.nextsampvarb[gdat.indxsamppsfipara] = copy(gdat.thissampvarb[gdat.indxsamppsfipara])
        gdat.nextsampvarb[gdat.indxsamppsfipara[gdat.indxpsfiparamodi]] = \
            icdf_psfipara(gdat, gdat.drmcsamp[gdat.indxsampmodi, -1], gdat.indxpsfiparamodi)
            
        gdat.numbmodipnts = int(sum(gdat.thissampvarb[gdat.indxsampnumbpnts]))
        
        if gdat.verbtype > 1:
           
            print 'indxpsfiparamodi: ', gdat.indxpsfiparamodi
            print 'indxenermodi: ', gdat.indxenermodi
            print 'indxsampmodi: ', gdat.indxsampmodi
            print 'gdat.drmcsamp[gdat.indxsampmodi, -1]'
            print gdat.drmcsamp[gdat.indxsampmodi, -1]
            print 'icdf_psfipara(gdat, gdat.drmcsamp[gdat.indxsampmodi, -1], gdat.indxpsfiparamodi)'
            print icdf_psfipara(gdat, gdat.drmcsamp[gdat.indxsampmodi, -1], gdat.indxpsfiparamodi)
            print 'thispsfipara: ', gdat.thissampvarb[gdat.indxsamppsfipara]
            print 'thissampvarb[indxsampmodi]: ', gdat.thissampvarb[gdat.indxsampmodi]
            print 'nextpsfipara: ', gdat.nextsampvarb[gdat.indxsamppsfipara]
            print 'nextpsfipara[indxsampmodi]: ', gdat.nextsampvarb[gdat.indxsampmodi]
            print

        
    # background changes
    
    # diffuse model
    if gdat.thisindxprop == gdat.indxpropnormback:

        # determine the sample index to be changed
        gdat.indxenermodi = choice(gdat.indxener)
        gdat.indxbackmodi = choice(gdat.indxback)
        gdat.indxsampmodi = gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]
        
        # propose
        retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvback)

        # transform back from the unit space
        gdat.nextsampvarb[gdat.indxsampnormback] = copy(gdat.thissampvarb[gdat.indxsampnormback])
        gdat.nextsampvarb[gdat.indxsampmodi] = icdf_logt(gdat.drmcsamp[gdat.indxsampmodi, -1], \
            gdat.minmnormback[gdat.indxbackmodi], gdat.factnormback[gdat.indxbackmodi])

        if gdat.verbtype > 1:
            print 'indxenermodi: ', gdat.indxenermodi
            print 'indxbackmodi: ', gdat.indxbackmodi
            print 'indxsampmodi: ', gdat.indxsampmodi
            print 'thissampvarb[gdat.indxsampnormback]: ', gdat.thissampvarb[gdat.indxsampnormback]
            print 'nextsampvarb[gdat.indxsampnormback]: ', gdat.nextsampvarb[gdat.indxsampnormback]
            print 'gdat.drmcsamp[gdat.indxsampnormback, :]'
            print gdat.drmcsamp[gdat.indxsampnormback, :]
            print

    # birth
    if gdat.thisindxprop == gdat.indxpropbrth:

        # change the number of PS
        gdat.nextsampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] = gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxsampbrth = int(gdat.indxsampcompinit + gdat.maxmnumbcomp * gdat.indxpoplmodi + gdat.thisindxpntsempt[gdat.indxpoplmodi][0] * gdat.numbcomp)
        
        # sample auxiliary variables
        if gdat.colrprio:
            numbauxipara = gdat.numbcompcolr
        else:
            numbauxipara = gdat.numbcomp
        gdat.auxipara = rand(numbauxipara)

        if gdat.colrprio:
            gdat.drmcsamp[indxsampbrth:indxsampbrth+2, -1] = gdat.auxipara[0:2]
            gdat.drmcsamp[indxsampbrth+2+gdat.indxenerfdfn, -1] = gdat.auxipara[-2]
            gdat.drmcsamp[indxsampbrth+gdat.numbcomp-1, -1] = gdat.auxipara[-1]
        else:   
            gdat.drmcsamp[indxsampbrth:indxsampbrth+gdat.numbcomp, -1] = gdat.auxipara

        # sample indices to be modified
        gdat.indxsampmodi = arange(indxsampbrth, indxsampbrth + gdat.numbcomp, dtype=int)

        # modification catalog
        gdat.numbmodipnts = 1
        gdat.modilgal[0] = icdf_self(gdat.drmcsamp[indxsampbrth, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        gdat.modibgal[0] = icdf_self(gdat.drmcsamp[indxsampbrth+1, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        if gdat.colrprio:
            modiflux = icdf_spec(gdat, gdat.drmcsamp[indxsampbrth+2+gdat.indxenerfdfn, -1], \
                    gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenerfdfn]], \
                    gdat.minmspec[gdat.indxenerfdfn], gdat.maxmspec[gdat.indxenerfdfn])
            gdat.modisind[0] = icdf_eerr(gdat.drmcsamp[indxsampbrth+gdat.numbcomp-1, -1], gdat.meansind[gdat.indxpoplmodi], \
                        gdat.stdvsind[gdat.indxpoplmodi], gdat.sindcdfnnormminm[gdat.indxpoplmodi], gdat.sindcdfnnormdiff[gdat.indxpoplmodi])
            gdat.modispec[:, 0] = retr_spec(gdat, modiflux, gdat.modisind[0]).flatten()
        else:
            for i in gdat.indxener:
                gdat.modispec[i, 0] = icdf_spec(gdat, gdat.drmcsamp[indxsampbrth+2+i, -1], \
                    gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]], gdat.minmspec[i], gdat.maxmspec[i])
    
        if gdat.verbtype > 1:
            print 'auxipara: ', gdat.auxipara
            print 'modilgal: ', gdat.modilgal
            print 'modibgal: ', gdat.modibgal
            print 'modispec: '
            print gdat.modispec
            print
            
    # kill
    if gdat.thisindxprop == gdat.indxpropdeth:
        
        # change the number of PS
        gdat.nextsampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] = gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        gdat.killindxpnts = gdat.thisindxpntsfull[gdat.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        gdat.indxsampmodi = array([])
            
        # modification catalog
        gdat.numbmodipnts = 1
        gdat.modilgal[0] = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][killindxindxpnts]]
        gdat.modibgal[0] = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][killindxindxpnts]]
        gdat.modispec[:, 0] = -gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][:, killindxindxpnts]]

        if gdat.verbtype > 1:
            print 'killindxpnts: ', gdat.killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', gdat.modilgal
            print 'modibgal: ', gdat.modibgal
            print 'modispec: '
            print gdat.modispec
            print
            
  
    # split
    if gdat.thisindxprop == gdat.indxpropsplt:
        
        gdat.nextsampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] = gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] + 1
        
        # determine which point source to split
        spltindxindxpnts = choice(arange(gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]], dtype=int))
        spltindxpnts = gdat.thisindxpntsfull[gdat.indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        gdat.indxsampchd0 = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdat.indxpoplmodi + gdat.thisindxpntsfull[gdat.indxpoplmodi][spltindxindxpnts] * gdat.numbcomp
        indxfinl0 = gdat.indxsampchd0 + gdat.numbcomp
        gdat.indxsampchd1 = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdat.indxpoplmodi + gdat.thisindxpntsempt[gdat.indxpoplmodi][0] * gdat.numbcomp
        indxfinl1 = gdat.indxsampchd1 + gdat.numbcomp
        
        # determine the modified sample vector indices
        gdat.indxsampmodi = concatenate((arange(gdat.indxsampchd0, indxfinl0, dtype=int), arange(gdat.indxsampchd1, indxfinl1, dtype=int)))
        
        thislgal = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][spltindxindxpnts]]
        thisbgal = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][spltindxindxpnts]]
        thisspec = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][:, spltindxindxpnts]]
        
        if gdat.verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxsampchd0: ', gdat.indxsampchd0
            print 'indxfinl0: ', indxfinl0
            print 'indxsampchd1: ', gdat.indxsampchd1
            print 'indxfinl1: ', indxfinl1
            if pixltype == 'heal':
                print 'thislgal: ', thislgal
                print 'thisbgal: ', thisbgal
            else:
                print 'thislgal: ', 3600. * thislgal
                print 'thisbgal: ', 3600. * thisbgal
            print 'thisspec: ', thisspec
            
            
        # determine the new components
        gdat.auxipara = empty(gdat.numbcomp)
        gdat.auxipara[0:2] = rand(2) * spmrlbhl
        gdat.auxipara[2:] = (exp(rand(gdat.numbener)) - 1.) / (exp(1.) - 1.) * (gdat.maxmspec - gdat.minmspec) + gdat.minmspec
        
        if gdat.verbtype > 1:
            if pixltype == 'heal':
                print 'auxipara[0]: ', gdat.auxipara[0]
                print 'auxipara[1]: ', gdat.auxipara[1]
            else:
                print 'auxipara[0]: ', 3600. * gdat.auxipara[0]
                print 'auxipara[1]: ', 3600. * gdat.auxipara[1]
            print 'auxipara[2:]: ', gdat.auxipara[2:]
            print
            
        nextlgal0 = thislgal + gdat.auxipara[0]
        nextlgal1 = thislgal - gdat.auxipara[0]
        nextbgal0 = thisbgal + gdat.auxipara[1]
        nextbgal1 = thisbgal - gdat.auxipara[1]
        nextspec0 = (thisspec + gdat.auxipara[2:]) / 2.
        nextspec1 = (thisspec - gdat.auxipara[2:]) / 2.
        
        if gdat.verbtype > 1:
            if pixltype == 'heal':
                print 'nextlgal0: ', nextlgal0
                print 'nextlgal1: ', nextlgal1
                print 'nextbgal0: ', nextbgal0
                print 'nextbgal1: ', nextbgal1
            else:
                print 'nextlgal0: ', 3600. * nextlgal0
                print 'nextlgal1: ', 3600. * nextlgal1
                print 'nextbgal0: ', 3600. * nextbgal0
                print 'nextbgal1: ', 3600. * nextbgal1
            print 'nextspec0: ', nextspec0
            print 'nextspec1: ', nextspec1

        if abs(nextlgal0) > gdat.maxmgangmarg or abs(nextlgal1) > gdat.maxmgangmarg or \
           abs(nextbgal0) > gdat.maxmgangmarg or abs(nextbgal1) > gdat.maxmgangmarg or \
           where((nextspec0 > gdat.maxmspec) | (nextspec0 < gdat.minmspec))[0].size > 0 or \
           where((nextspec1 > gdat.maxmspec) | (nextspec1 < gdat.minmspec))[0].size > 0:
               gdat.reje = True
                
        if not gdat.reje:

            lgal = concatenate((array([nextlgal0, nextlgal1]), setdiff1d(gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgal0, nextbgal1]), setdiff1d(gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi]], thisbgal)))
            pairlist = retr_pairlist(lgal, bgal)

            ## first new component
            gdat.drmcsamp[gdat.indxsampchd0, -1] = cdfn_self(nextlgal0, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdat.drmcsamp[gdat.indxsampchd0+1, -1] = cdfn_self(nextbgal0, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            for i in gdat.indxener:
                gdat.drmcsamp[gdat.indxsampchd0+2+i, -1] = cdfn_spec(nextspec0[i], gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]], gdat.minmspec[i], gdat.maxmspec[i])

            ## second new component
            gdat.drmcsamp[gdat.indxsampchd1, -1] = cdfn_self(nextlgal1, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdat.drmcsamp[gdat.indxsampchd1+1, -1] = cdfn_self(nextbgal1, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            for i in gdat.indxener:
                gdat.drmcsamp[gdat.indxsampchd1+2+i, -1] = cdfn_spec(nextspec1[i], gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]], gdat.minmspec[i], gdat.maxmspec[i])

            gdat.numbmodipnts = 3
            ## component to be removed
            gdat.modilgal[0] = thislgal
            gdat.modibgal[0] = thisbgal
            gdat.modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            gdat.modilgal[1] = nextlgal0
            gdat.modibgal[1] = nextbgal0
            gdat.modispec[:, 1] = nextspec0.flatten()

            # second component to be added
            gdat.modilgal[2] = nextlgal1
            gdat.modibgal[2] = nextbgal1
            gdat.modispec[:, 2] = nextspec1.flatten()

    if gdat.thisindxprop == gdat.indxpropmerg:
        
        gdat.nextsampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] = gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi]], gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi]]])
            
        lgal = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi]]
        bgal = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi]]
        pairlist = retr_pairlist(lgal, bgal)
        
        if gdat.verbtype > 1:
            print 'lgal'
            print lgal
            print 'bgal'
            print bgal
            print 'pairlist'
            print pairlist
            
            
        if len(pairlist) == 0:
            gdat.reje = True
        else:
            gdat.reje = False
            jpair = choice(arange(len(pairlist)))
            mergindxindxpnts0 = pairlist[jpair][0]
            mergindxindxpnts1 = pairlist[jpair][1]
  
        if not gdat.reje:

            # fisrt PS index to be merged
            mergindxchd0 = gdat.thisindxpntsfull[gdat.indxpoplmodi][mergindxindxpnts0]
            mergindxsampinit0 = gdat.indxsampcompinit + mergindxchd0 * gdat.numbcomp

            # second PS index to be merged
            gdat.mergindxchd1 = gdat.thisindxpntsfull[gdat.indxpoplmodi][mergindxindxpnts1]
            mergindxsampinit1 = gdat.indxsampcompinit + gdat.mergindxchd1 * gdat.numbcomp

            # determine the modified sample vector indices
            gdat.indxsampchd0 = gdat.indxsampcompinit + gdat.numbcomp * mergindxchd0
            indxfinl0 = gdat.indxsampchd0 + gdat.numbcomp
            gdat.indxsampchd1 = gdat.indxsampcompinit + gdat.numbcomp * gdat.mergindxchd1
            indxfinl1 = gdat.indxsampchd1 + gdat.numbcomp

            gdat.indxsampmodi = arange(gdat.indxsampchd0, indxfinl0)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxchd0, gdat.mergindxchd1], dtype=int))

            thislgal0 = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][mergindxindxpnts0]]
            thisbgal0 = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][mergindxindxpnts0]]
            thisspec0 = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][:, mergindxindxpnts0]]

            thislgal1 = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][mergindxindxpnts1]]
            thisbgal1 = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][mergindxindxpnts1]]
            thisspec1 = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][:, mergindxindxpnts1]]

            # auxiliary component
            gdat.auxipara = zeros(gdat.numbcomp)
            gdat.auxipara[0] = (thislgal0 - thislgal1) / 2.
            gdat.auxipara[1] = (thisbgal0 - thisbgal1) / 2.
            gdat.auxipara[2:] = thisspec0 - thisspec1

            # merged PS
            nextlgal = (thislgal0 + thislgal1) / 2.
            nextbgal = (thisbgal0 + thisbgal1) / 2.
            nextspec = thisspec0 + thisspec1
            
            gdat.drmcsamp[gdat.indxsampchd0, -1] = cdfn_self(nextlgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdat.drmcsamp[gdat.indxsampchd0+1, -1] = cdfn_self(nextbgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            for i in gdat.indxener:
                gdat.drmcsamp[gdat.indxsampchd0+2+i, -1] = cdfn_spec(nextspec[i], gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, i]], gdat.minmspec[i], gdat.maxmspec[i])

            gdat.numbmodipnts = 3
            ## first component to be merged
            gdat.modilgal[0] = thislgal0
            gdat.modibgal[0] = thisbgal0
            gdat.modispec[:, 0] = -thisspec0.flatten()

            ## first component to be merged
            gdat.modilgal[1] = thislgal1
            gdat.modibgal[1] = thisbgal1
            gdat.modispec[:, 1] = -thisspec1.flatten()

            ## component to be added
            gdat.modilgal[2] = nextlgal
            gdat.modibgal[2] = nextbgal
            gdat.modispec[:, 2] = nextspec.flatten()

            if gdat.verbtype > 1:
                print 'mergindxchd0: ', mergindxchd0
                print 'mergindxindxpnts0: ', mergindxindxpnts0
                print 'mergindxchd1: ', gdat.mergindxchd1
                print 'mergindxindxpnts1: ', mergindxindxpnts1
                print 'indxsampchd0: ', gdat.indxsampchd0
                print 'indxfinl0: ', indxfinl0
                print 'indxsampchd1: ', gdat.indxsampchd1
                print 'indxfinl1: ', indxfinl1
                if pixltype == 'heal':
                    print 'thislgal0: ', thislgal0
                    print 'thisbgal0: ', thisbgal0
                    print 'thislgal1: ', thislgal1
                    print 'thisbgal1: ', thisbgal1
                else:
                    print 'thislgal0: ', 3600. * thislgal0
                    print 'thisbgal0: ', 3600. * thisbgal0
                    print 'thislgal1: ', 3600. * thislgal1
                    print 'thisbgal1: ', 3600. * thisbgal1 
                print 'thisspec0: ', thisspec0
                print 'thisspec1: ', thisspec1

                if pixltype == 'heal':
                    print 'nextlgal: ', nextlgal
                    print 'nextbgal: ', nextbgal
                    print 'auxipara[0]: ', gdat.auxipara[0]
                    print 'auxipara[1]: ', gdat.auxipara[1]
                else:
                    print 'nextlgal: ', 3600. * nextlgal
                    print 'nextbgal: ', 3600. * nextbgal
                    print 'auxipara[0]: ', 3600. * gdat.auxipara[0]
                    print 'auxipara[1]: ', 3600. * gdat.auxipara[1]
                print 'nextspec: ', nextspec
                print 'auxipara[2:]: ', gdat.auxipara[2:]
                print

    # component change
    if gdat.thisindxprop >= gdat.indxproplgal:     
        
        if gdat.thisindxprop == gdat.indxproplgal or gdat.thisindxprop == gdat.indxpropbgal:
            if gdat.thisindxprop == gdat.indxproplgal:
                gdat.indxcompmodi = 0
            else:
                gdat.indxcompmodi = 1
            gdat.indxenermodi = gdat.indxener
        else:
            if gdat.colrprio:
                gdat.indxenermodi = gdat.indxener
                if gdat.thisindxprop == gdat.indxpropspec:
                    gdat.indxcompmodi = 2 + gdat.indxenerfdfn
                elif gdat.thisindxprop == gdat.indxpropsind:
                    gdat.indxcompmodi = array([2 + gdat.numbener])
            else:
                gdat.indxenermodi = choice(gdat.indxener)
                gdat.indxcompmodi = gdat.indxenermodi + 2
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = gdat.thisindxpntsfull[gdat.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        gdat.indxsampmodiinit = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdat.indxpoplmodi + modiindxpnts * gdat.numbcomp
        
        # sample index to be modified
        gdat.indxsampmodi = gdat.indxsampmodiinit + gdat.indxcompmodi
        if gdat.colrprio:
            gdat.indxsampmodispec = gdat.indxsampmodiinit + 2 + gdat.indxener
        
        # propose
        if gdat.thisindxprop == gdat.indxproplgal or gdat.thisindxprop == gdat.indxpropbgal:
            retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvlbhl) 
        if gdat.thisindxprop == gdat.indxpropspec:
            retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvspec)
        if gdat.thisindxprop == gdat.indxpropsind:
            retr_gaus(gdat, gdat.indxsampmodi, gdat.stdvpropsind)

        if gdat.colrprio:
            thisflux = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenerfdfn, modiindxindxpnts]]
            thissind = gdat.thissampvarb[gdat.thisindxsampsind[gdat.indxpoplmodi][modiindxindxpnts]]
            thisspec = retr_spec(gdat, thisflux, thissind)
        else:
            thisspec = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenermodi, modiindxindxpnts]]
            
        gdat.numbmodipnts = 2
        gdat.modispec[:, 0] = -thisspec.flatten()
        if gdat.thisindxprop == gdat.indxproplgal or gdat.thisindxprop == gdat.indxpropbgal:
            if gdat.indxcompmodi == 0:
                gdat.modilgal[0] = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][modiindxindxpnts]]
                gdat.modilgal[1] = icdf_self(gdat.drmcsamp[gdat.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
                gdat.modibgal[:] = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][modiindxindxpnts]]
            else:
                gdat.modilgal[:] = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][modiindxindxpnts]]
                gdat.modibgal[0] = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][modiindxindxpnts]]
                gdat.modibgal[1] = icdf_self(gdat.drmcsamp[gdat.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdat.modispec[:, 1] = thisspec.flatten()
        else:
            gdat.modilgal[:] = gdat.thissampvarb[gdat.thisindxsamplgal[gdat.indxpoplmodi][modiindxindxpnts]]
            gdat.modibgal[:] = gdat.thissampvarb[gdat.thisindxsampbgal[gdat.indxpoplmodi][modiindxindxpnts]]
            if gdat.colrprio:
                if gdat.thisindxprop == gdat.indxpropspec:
                    modiflux = icdf_spec(gdat, gdat.drmcsamp[gdat.indxsampmodi, -1], \
                        gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenerfdfn]], 
                        gdat.minmspec[gdat.indxenerfdfn], gdat.maxmspec[gdat.indxenerfdfn])
                    gdat.modisind[1] = gdat.thissampvarb[gdat.thisindxsampsind[gdat.indxpoplmodi][modiindxindxpnts]]
                else:
                    modiflux = gdat.thissampvarb[gdat.thisindxsampspec[gdat.indxpoplmodi][gdat.indxenerfdfn, modiindxindxpnts]]
                    gdat.modisind[1] = icdf_eerr(gdat.drmcsamp[gdat.indxsampmodi, -1], gdat.meansind[gdat.indxpoplmodi], \
                            gdat.stdvsind[gdat.indxpoplmodi], gdat.sindcdfnnormminm[gdat.indxpoplmodi], gdat.sindcdfnnormdiff[gdat.indxpoplmodi])
                gdat.modispec[:, 1] = retr_spec(gdat, modiflux, gdat.modisind[1]).flatten()
            else:
                specunit = gdat.drmcsamp[gdat.indxsampmodi, -1]
                fdfnslop = gdat.thissampvarb[gdat.indxsampfdfnslop[gdat.indxpoplmodi, gdat.indxenermodi]]
                gdat.modispec[:, 1] = icdf_spec(gdat, specunit, fdfnslop, gdat.minmspec[gdat.indxenermodi], gdat.maxmspec[gdat.indxenermodi])

        # log
        if gdat.verbtype > 1:
            print 'modilgal: ', gdat.modilgal
            print 'modibgal: ', gdat.modibgal
            print 'modispec: '
            print gdat.modispec
            print 'indxcompmodi: ', gdat.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts


    # energy bin in which to evaluate the log-likelihood
    if gdat.indxpropbrth <= gdat.thisindxprop <= gdat.indxpropmerg:
        gdat.indxenermodi = arange(gdat.numbener)

    if type(gdat.indxenermodi) == int64:
        gdat.indxenermodi = array([gdat.indxenermodi])

    if gdat.verbtype > 1:
        print 'indxsampmodi: ', gdat.indxsampmodi
        print 'indxenermodi: ', gdat.indxenermodi

    # auxiliary variable density fraction and jacobian
    if (gdat.thisindxprop == gdat.indxpropsplt or gdat.thisindxprop == gdat.indxpropmerg) and not gdat.reje:

        spltgdat.combfact = log(gdat.thissampvarb[gdat.indxsampnumbpnts[gdat.indxpoplmodi]]**2 / len(pairlist))
        
        if gdat.thisindxprop == gdat.indxpropsplt:
            thisgdat.combfact = spltgdat.combfact 
            thisgdat.jcbnfact = spltgdat.jcbnfact
        else:
            thisgdat.combfact = -spltgdat.combfact 
            thisgdat.jcbnfact = -spltgdat.jcbnfact


        gdat.laccfrac = thisgdat.jcbnfact + thisgdat.combfact

        gdat.listgdat.numbpair[gdat.thisindxswep] = len(pairlist)
        gdat.listgdat.jcbnfact[gdat.thisindxswep] = thisgdat.jcbnfact
        gdat.listgdat.combfact[gdat.thisindxswep] = thisgdat.combfact
        gdat.listgdat.auxipara[gdat.thisindxswep, :] = gdat.auxipara
        gdat.listgdat.laccfrac[gdat.thisindxswep] = gdat.laccfrac

    else:
        gdat.laccfrac = 0.  
        
        
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
  
    thisangltemp = thisangl[None, :, None]

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
        psfn = retr_doubking(thisangltemp, frac, sigc, gamc, sigt, gamt)
            
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
    gdat.indxmodlpsfipara = arange(gdat.numbpsfipara)   

    minmformpara = zeros(gdat.numbformpara)
    maxmformpara = zeros(gdat.numbformpara)
    factformpara = zeros(gdat.numbformpara)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    if gdat.exprtype == 'ferm':
        minmanglpsfn = deg2rad(0.0001)
        maxmanglpsfn = deg2rad(5.)
        #minmanglpsfn = 0.01
        #maxmanglpsfn = 3.
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


def retr_propmodl(gdat):

    gdat.strgprop = ['fdfnnorm', 'fdfnslop', 'psfipara', 'normback', 'brth', 'deth', 'splt', 'merg', 'lgal', 'bgal', 'spec', 'sind']

    gdat.numbprop = len(gdat.strgprop)
    gdat.indxprop = arange(gdat.numbprop)
    gdat.indxpropfdfnnorm = 0
    gdat.indxpropfdfnslop = 1
    gdat.indxproppsfipara = 2
    gdat.indxpropnormback = 3
    gdat.indxpropbrth = 4
    gdat.indxpropdeth = 5
    gdat.indxpropsplt = 6
    gdat.indxpropmerg = 7
    gdat.indxproplgal = 8
    gdat.indxpropbgal = 9
    gdat.indxpropspec = 10
    gdat.indxpropsind = 11

    if gdat.probprop == None:

        probfdfnnorm = array([1.])
        probfdfnslop = array([1.])
            
        if gdat.proppsfn:
            probpsfipara = array([1.])
        else:
            probpsfipara = array([0.])
        probnormback = array([1.])
        
        probbrth = array([0.1 * sum(gdat.maxmnumbpnts)])
        probdeth = array([0.1 * sum(gdat.maxmnumbpnts)])
        probsplt = array([0. * sum(gdat.maxmnumbpnts)])
        probmerg = array([0. * sum(gdat.maxmnumbpnts)])
        
        problgal = array([sum(gdat.maxmnumbpnts) / 2.])
        probbgal = array([sum(gdat.maxmnumbpnts) / 2.])
        if gdat.colrprio:
            probspec = array([sum(gdat.maxmnumbpnts) / 2.])
            probsind = array([sum(gdat.maxmnumbpnts) / 2.])
        else:
            probspec = array([sum(gdat.maxmnumbpnts) / 2.])
            #probspec = array([sum(gdat.maxmnumbpnts) / 2.] * gdat.numbener)
            probsind = array([0.])
            
        gdat.probprop = concatenate((probfdfnnorm, probfdfnslop, probpsfipara, probnormback, probbrth, probdeth, \
                probsplt, probmerg, problgal, probbgal, probspec, probsind))
        gdat.probprop /= sum(gdat.probprop)
       

def retr_strgangl(gdat):

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

    
def retr_randunitpsfipara():

    while True:
        randunitpsfipara = rand(gdat.numbpsfipara)
        indxpar0 = 1
        if gdat.psfntype == 'doubgaus' or gdat.psfntype == 'gausking':
            indxpar1 = 2
        if gdat.psfntype == 'doubking':
            indxpar1 = 3
        thisbool = True
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                indx = m * numbpsfiparaevtt + i * numbformpara
                thisbool = thisbool and randunitpsfipara[indx+indxpar1] > randunitpsfipara[indx+indxpar0]
        if thisbool:
            break
        else:
            print 'Repeating the PSF parameter seed...'

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
        path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgexpo
        gdat.expo = pf.getdata(path)

        if gdat.pixltype == 'heal':
            gdat.expo = gdat.expo[gdat.indxcubeheal]
        else:
            gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))
            
        gdat.expo = gdat.expo[gdat.indxcubefilt]
    

def setp(gdat):
    
    # the number of processes (each with gdat.numbswep samples)
    if gdat.numbproc == None:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu':
            gdat.numbproc = 10
        else:
            gdat.numbproc = 1
        
    if gdat.exprtype == 'ferm':
        if gdat.strgback == None:
            gdat.strgback = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
        if gdat.lablback == None:
            gdat.lablback = [r'$\mathcal{I}$', r'$\mathcal{D}$']
        if gdat.nameback == None:
            gdat.nameback = ['normisot', 'normfdfm']

    gdat.strgfluxunit = retr_strgfluxunit(gdat)
    
    gdat.numbchrototl = 4
    gdat.numbchrollik = 7

    # number of bins
    gdat.numbspec = 10
    gdat.numbbins = 10

    gdat.minmnumbpnts = 1

    gdat.numbback = len(gdat.strgback)
    gdat.indxback = arange(gdat.numbback)
    
    gdat.numbevtt = gdat.indxevttincl.size
    gdat.indxevtt = arange(gdat.numbevtt)
    
    # energy axis
    gdat.numbener = gdat.indxenerincl.size
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.binsener = array([0.1, 0.3, 1., 3., 10., 100.])[gdat.indxenerinclbins]
    gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]

    gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
    gdat.indxener = arange(gdat.numbener, dtype=int)
    
    gdat.indxenerpivt = array([gdat.numbener / 2])
    if gdat.colrprio: 
        gdat.indxenerfdfn = gdat.indxenerpivt
    else:
        gdat.indxenerfdfn = gdat.indxener
    gdat.numbenerfdfn = gdat.indxenerfdfn.size
    gdat.enerfdfn = gdat.meanener[gdat.indxenerfdfn]
        
    gdat.minmspec = gdat.minmspec * (gdat.meanener[gdat.indxenerpivt] / gdat.meanener)**2
    gdat.maxmspec = gdat.maxmspec * (gdat.meanener[gdat.indxenerpivt] / gdat.meanener)**2
    if gdat.datatype == 'mock':
        if gdat.mockfdfntype == 'brok':
            gdat.mockfdfnsloplowr = tile(gdat.mockfdfnsloplowr, (1, gdat.numbener))
            gdat.mockfdfnslopuppr = tile(gdat.mockfdfnslopuppr, (1, gdat.numbener))
            gdat.mockfdfnbrek = tile(gdat.mockfdfnbrek, (1, gdat.numbener))
        else:
            gdat.mockfdfnslop = tile(gdat.mockfdfnslop, (1, gdat.numbener))

    gdat.indxspecpivt = gdat.numbspec / 2
    
    retr_strgprop(gdat)

    # temp
    if gdat.exprtype == 'sdss':
        gdat.diffener = ones(gdat.numbener)
        
    # angular gdat.deviation
    gdat.numbangl = 100
    if gdat.exprtype == 'sdss':
        gdat.maxmangldisp = deg2rad(10. / 3600.) # [rad]
    if gdat.exprtype == 'ferm':
        gdat.maxmangldisp = deg2rad(10.) # [rad]
    gdat.angldisp = linspace(0., gdat.maxmangldisp, gdat.numbangl) # [rad]
    maxmangl = deg2rad(3.5 * gdat.maxmgang) # [rad]
    angl = linspace(0., maxmangl, gdat.numbangl) # [rad]

    if gdat.exprtype == 'sdss':
        gdat.angldisptemp = rad2deg(gdat.angldisp) * 3600.
    if gdat.exprtype == 'ferm':
        gdat.angldisptemp = rad2deg(gdat.angldisp)

    if gdat.trueinfo:
        if gdat.datatype == 'mock':
            gdat.truelabl = 'Mock'
        if gdat.datatype == 'inpt':
            gdat.truelabl = '3FGL'
            
    if gdat.exprtype == 'ferm':
        gdat.strganglunit = '[deg]'
    if gdat.exprtype == 'sdss':
        gdat.strganglunit = '[arcsec]'
        
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
    gdat.numbcomp = 2 + gdat.numbener
    if gdat.colrprio:
        gdat.numbcomp += 1
    gdat.numbcompcolr = 4
    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
    
    # population index vector
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)

    # hyperpriors
    gdat.minmfdfnbrek = empty((gdat.numbpopl, gdat.numbener))
    gdat.maxmfdfnbrek = empty((gdat.numbpopl, gdat.numbener))
    for l in gdat.indxpopl:
        gdat.minmfdfnbrek[l, :] = gdat.minmspec
        gdat.maxmfdfnbrek[l, :] = gdat.maxmspec

    # convenience factors for CDF and ICDF transforms
    gdat.factfdfnbrek = log(gdat.maxmfdfnbrek / gdat.minmfdfnbrek)
    gdat.factfdfnnorm = log(gdat.maxmfdfnnorm / gdat.minmfdfnnorm)
    gdat.factfdfnslop = arctan(gdat.maxmfdfnslop) - arctan(gdat.minmfdfnslop)
    if gdat.colrprio:
        gdat.factsind = arctan(gdat.maxmsind) - arctan(gdat.minmsind)
        gdat.sindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.meansind) / gdat.stdvsind / sqrt(2.)) + 1.)
        gdat.sindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.meansind) / gdat.stdvsind / sqrt(2.)) + 1.)
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

    # modification catalog variables
    numbmodipnts = int(max(3, gdat.maxmnumbpnts))
    gdat.modilgal = empty(numbmodipnts)
    gdat.modibgal = empty(numbmodipnts)
    gdat.modisind = empty(numbmodipnts)
    gdat.modispec = empty((gdat.numbener, numbmodipnts))

    if gdat.regitype == 'igal':
        gdat.longlabl = '$l$'
        gdat.latilabl = '$b$'
    else:
        gdat.longlabl = r'$\nu$'
        gdat.latilabl = r'$\mu$'
        
    if gdat.exprtype == 'ferm':
        gdat.longlabl += r' [$^\circ$]'
        gdat.latilabl += r' [$^\circ$]'
    if gdat.exprtype == 'sdss':
        gdat.longlabl += ' [arcsec]'
        gdat.latilabl += ' [arcsec]'
    
    # construct the PSF model
    retr_psfimodl(gdat)

    # proposals
    retr_propmodl(gdat)
    
    # factors in the prior expression
    gdat.priofactlgalbgal = 2. * log(1. / 2. / gdat.maxmgang)
    gdat.priofactfdfnslop = gdat.numbener * log(1. / (arctan(gdat.maxmfdfnslop) - arctan(gdat.minmfdfnslop)))
    gdat.priofactfdfnnorm = log(1. / (log(gdat.maxmfdfnnorm) - log(gdat.minmfdfnnorm)))

    # sample vector indices  
    gdat.indxsampnumbpnts = arange(gdat.numbpopl)
    gdat.indxsampfdfnnorm = arange(gdat.numbpopl) + amax(gdat.indxsampnumbpnts) + 1
    if gdat.fdfntype == 'brok':
        gdat.indxsampfdfnsloplowr = arange(gdat.numbpopl * gdat.numbener).reshape((gdat.numbpopl, gdat.numbener)) + amax(gdat.indxsampfdfnnorm) + 1
        gdat.indxsampfdfnslopuppr = arange(gdat.numbpopl * gdat.numbener).reshape((gdat.numbpopl, gdat.numbener)) + amax(gdat.indxsampfdfnsloplowr) + 1
        gdat.indxsampfdfnbrek = arange(gdat.numbpopl * gdat.numbener).reshape((gdat.numbpopl, gdat.numbener)) + amax(gdat.indxsampfdfnslopuppr) + 1
        indxsamptemp = amax(gdat.indxsampfdfnbrek) + 1
    else:
        gdat.indxsampfdfnslop = arange(gdat.numbpopl * gdat.numbener).reshape((gdat.numbpopl, gdat.numbener)) + amax(gdat.indxsampfdfnnorm) + 1
        indxsamptemp = amax(gdat.indxsampfdfnslop) + 1
    gdat.indxsamppsfipara = arange(gdat.numbpsfipara) + indxsamptemp
    gdat.indxsampnormback = arange(gdat.numbback * gdat.numbener).reshape((gdat.numbback, gdat.numbener)) + amax(gdat.indxsamppsfipara) + 1

    gdat.fluxpivt = sqrt(gdat.minmspec * gdat.maxmspec)
    
    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.indxsampcompinit = amax(gdat.indxsampnormback) + 1
    
    # maximum number of parameters
    gdat.maxmsampsize = int(gdat.indxsampcompinit + gdat.maxmnumbcomp * gdat.numbpopl)
    
    if gdat.numbburn == None:
        gdat.numbburn = gdat.numbswep / 5
    if gdat.factthin == None:
        gdat.factthin = min(gdat.maxmsampsize * 5, gdat.numbswep / 2)


    # run tag
    gdat.rtag = retr_rtag(gdat, None)
    
    # plots
    if gdat.makeplot:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
            plotfold = '/n/pan/www/tansu/png/pcat/'
        else:
            plotfold = os.environ["PCAT_DATA_PATH"] + '/png/'
        gdat.plotpath = plotfold + gdat.strgtime + '_' + gdat.rtag + '/'
        cmnd = 'mkdir -p ' + gdat.plotpath
        os.system(cmnd)

    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin

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
    
    gdat.maxmgangmarg = gdat.maxmgang + gdat.margsize
    
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
        
        path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgexpr
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
    gdat.mrkralph = 0.8
    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 500
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

    # sky coordinates
    gdat.binslbhl = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbbins + 1)

    gdat.binsspec = zeros((gdat.numbener, gdat.numbspec + 1))
    for i in gdat.indxener:
        gdat.binsspec[i, :] = logspace(log10(gdat.minmspec[i]), log10(gdat.maxmspec[i]), gdat.numbspec + 1)
    gdat.meanspec = sqrt(gdat.binsspec[:, 1:] * gdat.binsspec[:, 0:-1])
    gdat.diffspec = gdat.binsspec[:, 1:] - gdat.binsspec[:, 0:-1]

    if gdat.colrprio:
        gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbbins + 1)
        gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    if gdat.exprtype == 'ferm':
        gdat.numbspecprox = 4
    if gdat.exprtype == 'sdss':
        gdat.numbspecprox = 1
    gdat.indxspecprox = arange(gdat.numbspecprox)
    gdat.binsspecprox = logspace(log10(gdat.minmspec[gdat.indxenerpivt]), log10(gdat.maxmspec[gdat.indxenerpivt]), gdat.numbspecprox + 1)
    gdat.meanspecprox = sqrt(gdat.binsspecprox[1:] * gdat.binsspecprox[:-1])
    gdat.specfraceval = 1e-2
    if gdat.exprtype == 'ferm':
        gdat.maxmangleval = empty(gdat.numbspecprox)
        for h in gdat.indxspecprox:
            frac = gdat.specfraceval * gdat.binsspecprox[0] / gdat.binsspecprox[h+1]
            gdat.maxmangleval[h] = rad2deg(amax(retr_psfnwdth(gdat, gdat.fermpsfn, frac)))
    if gdat.exprtype == 'sdss':
        gdat.maxmangleval = array([10. / 3600.])

    # pizelization
    if gdat.pixltype == 'heal':
        
        lgalheal, bgalheal, gdat.numbsideheal, gdat.numbpixlheal, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)

        gdat.indxpixlrofi = where((abs(lgalheal) < gdat.maxmgang) & (abs(bgalheal) < gdat.maxmgang))[0]
        gdat.indxpixlrofimarg = where((abs(lgalheal) < gdat.maxmgangmarg + 300. / gdat.numbsideheal) & \
                (abs(bgalheal) < gdat.maxmgangmarg + 300. / gdat.numbsideheal))[0]
        
        gdat.lgalgrid = lgalheal[gdat.indxpixlrofi]
        gdat.bgalgrid = bgalheal[gdat.indxpixlrofi]
        
        path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%03d.p' % (gdat.maxmgang)
        if os.path.isfile(path):
            fobj = open(path, 'rb')
            gdat.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            gdat.pixlcnvt = zeros(gdat.numbpixlheal, dtype=int)
            for k in range(gdat.indxpixlrofimarg.size):
                dist = retr_dist(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.lgalgrid, gdat.bgalgrid)
                gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)

            fobj = open(path, 'wb')
            cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
    else:
        isidecart = arange(gdat.numbsidecart)
        temp = meshgrid(isidecart, isidecart, indexing='ij')
        gdat.bgalgrid = gdat.bgalcart[temp[1].flatten()]
        gdat.lgalgrid = gdat.lgalcart[temp[0].flatten()]
        
    gdat.numbpixl = gdat.lgalgrid.size
    gdat.indxpixl = arange(gdat.numbpixl)
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
        path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgexpo
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
            path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgback[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'heal':
                backfluxtemp = backfluxtemp[gdat.indxcubeheal]
            else:
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
            backfluxtemp = backfluxtemp[gdat.indxcubefilt]
        gdat.backflux.append(backfluxtemp)
        gdat.backfluxmean.append(mean(sum(gdat.backflux[c] * gdat.expo, 2) / sum(gdat.expo, 2), 1))

    # count axis
    gdat.expotemp = mean(gdat.expo, 1)
    gdat.minmcnts = repeat(1e-1, gdat.numbener) # gdat.minmspec * amin(gdat.expotemp, 1) * gdat.diffener
    gdat.maxmcnts = repeat(1e4, gdat.numbener) # gdat.maxmspec * amax(gdat.expotemp, 1) * gdat.diffener
    gdat.binscnts = zeros((gdat.numbener, gdat.numbspec + 1))
    for i in gdat.indxener:
        gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbspec + 1) # [1]
        
    # get 3FGL catalog
    if gdat.exprtype == 'ferm':
        retr_fgl3(gdat)

    # get count data
    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = exprflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]

    ## mock data
    if gdat.datatype == 'mock':

        if gdat.mocknumbpnts == None:
            gdat.mocknumbpnts = empty(gdat.numbpopl)
            for l in gdat.indxpopl:
                gdat.mocknumbpnts[l] = random_integers(gdat.minmnumbpnts, gdat.maxmnumbpnts[l])
        
        if gdat.mockfdfntype == 'brok':
            pass
        else:
            if gdat.mockfdfnslop == None:
                gdat.mockfdfnslop = empty((gdat.numbpopl, gdat.numbenerfdfn))
                for l in gdat.indxpopl:
                    gdat.mockfdfnslop[l, :] = icdf_atan(rand(gdat.numbenerfdfn), gdat.minmfdfnslop[l], gdat.factfdfnslop[l])
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
        if gdat.colrprio:
            mocksind = [[] for l in gdat.indxpopl]
        for l in gdat.indxpopl:
            mocklgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            mockbgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
            mockspec[l] = empty((gdat.numbener, gdat.mocknumbpnts[l]))
            for i in gdat.indxenerfdfn:
                if gdat.mockfdfntype == 'brok':
                    mockspec[l][i, :] = icdf_spec_brok(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfdfnsloplowr[l, i], \
                        gdat.mockfdfnslopuppr[l, i], gdat.mockfdfnbrek[l, i], gdat.minmspec[i], gdat.maxmspec[i])
                else:
                    mockspec[l][i, :] = icdf_spec(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfdfnslop[l, i], gdat.minmspec[i], gdat.maxmspec[i])
            if gdat.colrprio:
                mocksind[l] = icdf_eerr(rand(gdat.mocknumbpnts[l]), gdat.meansind[l], gdat.stdvsind[l], \
                                                                                        gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
                mockspec[l] = retr_spec(gdat, mockspec[l][gdat.indxenerfdfn, :].flatten(), mocksind[l])
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
            gdat.truespec = []
            if gdat.colrprio:
                gdat.truesind = []
            for i in gdat.indxpopl:
                gdat.truelgal.append(mocklgal[l])
                gdat.truebgal.append(mockbgal[l])
                gdat.truespec.append(mockspec[l])
                if gdat.colrprio:
                    gdat.truesind.append(mocksind[l])
                    
            gdat.indxtruepntstimevari = [array([])] * gdat.numbpopl
                    
            gdat.truenumbpnts = gdat.mocknumbpnts
            if gdat.mockfdfntype == 'brok':
                gdat.truefdfnsloplowr = gdat.mockfdfnsloplowr
                gdat.truefdfnslopuppr = gdat.mockfdfnslopuppr
                gdat.truefdfnbrek = gdat.mockfdfnbrek
            else:
                gdat.truefdfnslop = gdat.mockfdfnslop
            gdat.truenormback = gdat.mocknormback
            gdat.truecnts = mockcnts
            
            gdat.truespec = []
            for l in gdat.indxpopl:
                gdat.truespectemp = empty((3, gdat.numbener, gdat.mocknumbpnts[l]))
                gdat.truespectemp[:] = mockspec[l][None, :, :]
                gdat.truespec.append(gdat.truespectemp)
                
            gdat.truepsfipara = gdat.mockpsfipara
            gdat.truepsfntype = gdat.mockpsfntype


    ## Real data
    # true data
    if gdat.trueinfo:
        if gdat.datatype == 'inpt':
            gdat.truenumbpnts = None
            gdat.truefdfnslop = None
            gdat.truenormback = None
    
            if gdat.exprtype == 'ferm':
                gdat.truenumbpnts = array([gdat.fgl3numbpntsrofi], dtype=int)
                gdat.truelgal = [gdat.fgl3lgal[gdat.indxfgl3rofi]]
                gdat.truebgal = [gdat.fgl3bgal[gdat.indxfgl3rofi]]
                gdat.truespec = [gdat.fgl3spec[:, :, gdat.indxfgl3rofi]]
                if gdat.colrprio:
                    gdat.truesind = [gdat.fgl3sind[gdat.indxfgl3rofi]]
                indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[0], gdat.truelgal[0])
                spec = gdat.fgl3spec[0, :, gdat.indxfgl3rofi]
                # temp
                spec = spec.T
                gdat.truecnts = [spec[:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]]
                gdat.indxtruepntstimevari = [gdat.indxfgl3timevarirofi]
                gdat.truepsfipara = gdat.fermpsfipara
                gdat.truepsfntype = 'doubking'
                
        if gdat.datatype == 'mock':
            gdat.truepsfn = retr_psfn(gdat, gdat.truepsfipara, gdat.indxener, gdat.angldisp, gdat.mockpsfntype)
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
                truebackcntstemp += gdat.backflux[c][:, indxpixltemp, :] * gdat.expo[:, indxpixltemp, :] * \
                    gdat.diffener[:, None, None] * pi * truefwhm[:, None, :]**2 / 4.
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
    gdat.resicntssatu = ceil(gdat.datacntssatu * 0.5)
    
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
    # temp
    # does not work for colrprio = False
    # needs to be generated for the lowest energy bin, i.e., indxenerincl = 0
    path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%03d_%s_%.5g_%.5g_%02d.p' % (gdat.maxmgang, gdat.pixltype, \
                                                               gdat.minmspec[gdat.indxenerfdfn[0]], gdat.maxmspec[gdat.indxenerfdfn[0]], gdat.numbspecprox)
    if os.path.isfile(path):
        print 'Retrieving previously computed pixel look-up table...'
        fobj = open(path, 'rb')
        gdat.indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        print 'Computing the look-up table...'
        gdat.indxpixlprox = [[] for h in range(gdat.numbspecprox)]
        for j in gdat.indxpixl:
            dist = retr_dist(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.lgalgrid, gdat.bgalgrid)
            for h in range(gdat.numbspecprox):
                gdat.indxpixlproxtemp = where(dist < deg2rad(gdat.maxmangleval[h]))[0]
                gdat.indxpixlprox[h].append(gdat.indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()


def init_fram(gdat, indxevttplot, indxenerplot, strgplot):

    figr, axis = plt.subplots(figsize=(12, 12))
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
        path = gdat.plotpath + strgplot + '%dA_' % gdat.indxenerincl[indxenerplot] + gdat.rtag + '_%09d.png' % gdat.cntrswep
    else:
        path = gdat.plotpath + strgplot + '%d%d_' % (gdat.indxenerincl[indxenerplot], gdat.indxevttincl[indxevttplot]) + gdat.rtag + '_%09d.png' % gdat.cntrswep
    
    return figr, axis, path


def supr_fram(gdat, axis, indxenerplot):

    for l in gdat.indxpopl:

        # true catalog
        if gdat.trueinfo:
            mrkrsize = retr_mrkrsize(gdat, gdat.truespec[l][0, :, :], indxenerplot)
            lgal = copy(gdat.truelgal[l])
            bgal = copy(gdat.truebgal[l])
            if gdat.exprtype == 'sdss':
                lgal *= 3600.
                bgal *= 3600.
            axis.scatter(lgal[gdat.trueindxpntsmiss], bgal[gdat.trueindxpntsmiss], s=mrkrsize[gdat.trueindxpntsmiss], \
                alpha=gdat.mrkralph, label=gdat.truelabl + ', missed', marker='x', linewidth=2, color='g')
            axis.scatter(lgal[gdat.trueindxpntsbias], bgal[gdat.trueindxpntsbias], s=mrkrsize[gdat.trueindxpntsbias], \
                alpha=gdat.mrkralph, label=gdat.truelabl + ', biased', marker='o', linewidth=2, color='g')
            indxpnts = setdiff1d(arange(gdat.truenumbpnts, dtype=int), concatenate((gdat.trueindxpntsbias, gdat.trueindxpntsmiss)))
            axis.scatter(lgal[indxpnts], bgal[indxpnts], s=mrkrsize[indxpnts], alpha=gdat.mrkralph, label=gdat.truelabl + ', hit', marker='D', linewidth=2, color='g')
            for l in gdat.indxpopl:
                if gdat.indxtruepntstimevari[l].size > 0:
                    axis.scatter(lgal[gdat.indxtruepntstimevari[l]], bgal[gdat.indxtruepntstimevari[l]], s=100, \
                        label=gdat.truelabl + ', variable', marker='*', linewidth=2, color='y')

        # model catalog
        mrkrsize = retr_mrkrsize(gdat, gdat.thissampvarb[gdat.thisindxsampspec[l]], indxenerplot)
        lgal = gdat.thissampvarb[gdat.thisindxsamplgal[l]]
        bgal = gdat.thissampvarb[gdat.thisindxsampbgal[l]]
        if gdat.exprtype == 'sdss':
            lgal *= 3600.
            bgal *= 3600.
            
        axis.scatter(lgal, bgal, s=mrkrsize, alpha=gdat.mrkralph, label='Sample', marker='+', linewidth=2, color='b')
    axis.legend(bbox_to_anchor=[0.12, 1.1], loc='center', ncol=2)


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


